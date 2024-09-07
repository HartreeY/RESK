using StatsBase, Distributions, Distributed, Random, SpecialFunctions, Serialization, Dates, SharedArrays
include("defaults.jl")


# Configuration
# ------------------------------------------------

# Migration probabilities for each dimensionality -> for each mode
MIGR_PROBS = [
    Dict(["ort" => (1, 0)]), # 1D
    Dict(["ort" => (1, 0), "hex" => (1, 0), "all" => (1 / 2, 1 / 2), "buffon1" => (2 / pi, 1 / pi), "buffon2" => (4 / 3 / pi, 1 / 3 / pi), "buffon3" => (0.4244132, 0.21221), "diag1/2" => (2 / 3, 1 / 3)]), # 2D. To add "hex"!
    Dict(["ort" => (1, 0), "all" => (1 / 2, 1 / 2), "buffon1" => (2 / pi, 1 / pi), "buffon2" => (4 / 3 / pi, 1 / 3 / pi), "buffon3" => (0.4244132, 0.21221), "diag1/2" => (2 / 3, 1 / 3)]) # 3D. To add "hex"! To confirm Buffon for 3d!
]

# Migration directions for each mode -> for each dimensionality
MIGR_DIRS_ORT = [
    [[1], [-1]], # 1D
    [[-1, 0], [0, -1], [0, 1], [1, 0]], # 2D
    [[-1, 0, 0], [1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, -1], [0, 0, 1]] # 3D
]
MIGR_DIRS_DIAG = [
    [[-1, -1], [-1, 1], [1, -1], [1, 1]], # 2D
    [[-1, -1, -1], [-1, -1, 0], [-1, -1, 1], [-1, 0, -1], [-1, 0, 1], [-1, 1, -1], [-1, 1, 0], [-1, 1, 1], [0, -1, -1], [0, -1, 0], [0, -1, 1], [0, 0, 0], [0, 1, -1], [0, 1, 1], [1, -1, -1], [1, -1, 0], [1, -1, 1], [1, 0, -1], [1, 0, 1], [1, 1, -1], [1, 1, 0], [1, 1, 1]] # 3D
]
MIGR_DIRS_HEX = [
    [], # 1D
    [[-1, 0], [0, -1], [-1, 1], [0, 1], [1, 0], [1, 1]], # 2D
    [] # 3D
]


# Common functions
# ------------------------------------------------

"""
Chooses all values except the specific value `n` at the last dimension.
"""
function repl(data,n)
    return data[repeat([:], length(size(data))-1)...,n]
end

"""
Chooses all values except the specific value `n` at the last dimension. Other name: `repl(data,n)`.
"""
function li(data,n)
    return repl(data,n)
end

"""
(Abbreviation of "inscribed square") Returns a range of coordinates of world centre (determined from `r_max_exp`) ± side length of the inscribed square of a circle with radius `r_max_burnin`.

Used in determining the starting fillup of demes in radial expansions.

---

`r_max_burnin`: radius that bounds the burn-in area

`r_max_exp`: radius that bounds the expansion area

---

Output: integer range of coordinates around the world centre
"""
function ins_sq(r_max_burnin, r_max_exp)
    return (trunc(Int, 1 + r_max_exp - r_max_burnin * 0.666)):(trunc(Int, 1 + r_max_exp + r_max_burnin * 0.666))
end

"""
(Abbreviation of "inscribed cube") Returns a range of coordinates of world centre (determined from `r_max_exp`) ± side length of the inscribed cube of a sphere with radius `r_max_burnin`.

Used in determining the starting fillup of demes in spherical expansions.

---

`r_max_burnin`: radius that bounds the burn-in area

`r_max_exp`: radius that bounds the expansion area

---

Output: integer range of coordinates around the world centre
"""
function ins_cb(r_max_burnin, r_max_exp)
    return (trunc(Int, 1 + r_max_exp - r_max_burnin * 0.577)):(trunc(Int, 1 + r_max_exp + r_max_burnin * 0.577))
end


"""
Calculates the number of offspring individuals in currently filled demes. Used when building the next generation.

---

`wld`: world array (a spatial array of demes that contain individuals' loci [Bool or Float] arrays)

`wld_stats`: world stats Dict

---

Output 1: array of deme coordinates to be filled in the next generation

Output 2: array of populations for the coordinates in Output 1
"""
function calc_offspring(wld, wld_stats; fixedrand=false)
    next_gen_posits = []
    next_gen_pops = fill(NaN, wld_stats["max"]...)
    
    for k in Iterators.product([1:n for n in wld_stats["max"]]...)
        if isassigned(wld, k...) && length(wld[k...]) > 0
            n_ppl_at_deme = length(wld[k...])
            expected_offspring = n_ppl_at_deme * (wld_stats["prolif_rate"] / (1 + (n_ppl_at_deme * (wld_stats["prolif_rate"] - 1)) / wld_stats["capacity"]))
            #println("A",n_ppl_at_deme)
            next_gen_pops[k...] = fixedrand ? expected_offspring : rand(Poisson(expected_offspring))
            if next_gen_pops[k...] > 0
                push!(next_gen_posits, [k...])
            end
            #println("B",next_gen_posits)
        end
    end
    return next_gen_posits, next_gen_pops
end

"""
Calculates an individual's migration distance. Used when building the next generation.

---

`deme`: individual's current deme coordinates

`wld_stats`: world stats Dict

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`max_migr`: a tuple of maximum migration area coordinates

`refl_walls`: if **true**, walls reflect migrants

`r_max_migr`: Int maximum migration radius. If *>0**, migration is kept within this radius. Can be used in addition to `max_migr`

`r_coords`: a tuple (array) of axes' ordinal numbers that the n-sphere with `r_max_migr` covers. For example:
- **(1,3)** - migration is bound within a disk at x and z axes
- **(1,2,3)** - migration is bound within a sphere at x, y and z axes

---

Output: array of the amount of demes moved per each coordinate
"""
function calc_migr_dist(deme, wld_stats, migr_mode, bottleneck, max_migr=wld_stats["max"], refl_walls=false, r_max_migr=0, r_coords=[1, 2])

    wlddim = wld_stats["wlddim"]
    max = wld_stats["max"]
    move = zeros(Int16, wlddim)
    if !(migr_mode in keys(MIGR_PROBS[wlddim]))
        migr_mode = "ort"
    end
    

    if rand() < wld_stats["migr_rate"]
        if migr_mode=="hex"
            dir = copy(sample(MIGR_DIRS_HEX[wlddim]))
        else
            p_lat, p_diag = MIGR_PROBS[wlddim][migr_mode]
            migr_res = p_lat==1 ? 0.5 : rand()
            if migr_res < p_lat
                dir = copy(sample(MIGR_DIRS_ORT[wlddim]))
            elseif migr_res < p_lat + p_diag
                dir = copy(sample(MIGR_DIRS_DIAG[wlddim]))
            end
        end

        # Raw migration results
        move = copy(dir)
        
        # Nullify migration on certain conditions
        #------------------------------------------
        # Inside certain radius check:
        if r_max_migr > 0
            r_arr = [(deme[i] - (max[i]-1)/2 + move[i] - 1)^2 for i in r_coords]
            r2 = sum(r_arr)
            
            if r2 > r_max_migr * r_max_migr
                #factor = r_max_migr*r_max_migr/r2
                move[r_coords] .= 0 # Do [trunc(Int16, factor * move[i]) for i in 1:wlddim] in the future (with multiple-deme jumps)
            end
        end

        # Inside certain square check:
        if isa(max_migr, Tuple)
            for i in 1:wlddim
                if !isnan(max_migr[i])
                    try_move = deme[i] + move[i]
                    if try_move > max_migr[i] || try_move < 1
                        move[i] = refl_walls ? -move[i] : 0
                    end
                end
            end
        end

        # Bottleneck barrier check:
        if isa(bottleneck, Tuple) && isa(bottleneck[2], Int) && bottleneck[2] > 0
            if bottleneck[1] == "midhole at x="
                common_cond = deme[1] + move[1] == bottleneck[2] && deme[2] + move[2] != ceil(max[2] / 2)
                if wlddim == 2 && common_cond
                    move .= 0
                elseif wlddim == 3 && common_cond && deme[3] + move[3] != ceil(max[3] / 2)
                    move .= 0
                end
            elseif bottleneck[1] == "midhole at y="
                common_cond = deme[2] + move[2] == bottleneck[2] && deme[1] + move[1] != ceil(max[1] / 2)
                if wlddim == 2 && common_cond
                    move .= 0
                elseif wlddim == 3 && common_cond && deme[3] + move[3] != ceil(max[3] / 2)
                    move .= 0
                end
            elseif bottleneck[1] == "midhole at z=" && wlddim == 3 && deme[3] + move[3] == bottleneck[2] && deme[1] + move[1] != ceil(max[1] / 2) && deme[2] + move[2] != ceil(max[2] / 2)
                move .= 0
            end
        end
        #------------------------------------------
    end
    #println(move,deme)
    return move
end


# Simulation functions (finite-sites)
# ------------------------------------------------

"""
Calculates the number of mutations and average fitness within a deme. Used when building the next generation in finite-sites expansions.

---

`deme_ms1`: array of individuals (within a deme), each a left monosome loci array

`deme_ms2`: array of individuals (within a deme), each a right monosome loci array

`domin_coef`: dominance coefficient (in heterozygous loci, new_fitness *= **1 -** `domin_coef` * `sel_coef`)

`loci`: array of selected coefficients for every locus

`sel_loci`: array of selected loci

---

Output 1: average number of selected AA mutations in this deme

Output 2: average number of selected Aa mutations in this deme

Output 3: average number of selected aa mutations in this deme

Output 4: average number of neutral AA mutations in this deme

Output 5: average number of neutral Aa mutations in this deme

Output 6: average number of neutral aa mutations in this deme

Output 7: average fitness in this deme
"""
function calc_muts_and_fitn_in_deme(deme_ms1, deme_ms2, domin_coef, loci, sel_loci=[])
    len = length(deme_ms1)
    muts_AAsel_total = 0
    muts_Aasel_total = 0
    muts_aasel_total = 0
    muts_AAneu_total = 0
    muts_Aaneu_total = 0
    muts_aaneu_total = 0
    fits = []

    for i in 1:len
        muts_AA_sel = 0
        muts_Aa_sel = 0
        muts_AA_neu = 0
        muts_Aa_neu = 0
        new_fitness = 1.0

        for j in 1:length(loci)
            if deme_ms1[i][j] == true && deme_ms2[i][j] == true
                if j in sel_loci
                    muts_AA_sel += 1
                    new_fitness *= 1 - loci[j]
                else
                    muts_AA_neu += 1
                end

            elseif deme_ms1[i][j] == true || deme_ms2[i][j] == true
                if j in sel_loci
                    muts_Aa_sel += 1
                    new_fitness *= 1 - domin_coef * loci[j]
                else
                    muts_Aa_neu += 1
                end
            end
        end

        push!(fits, new_fitness)
        muts_AAsel_total += muts_AA_sel
        muts_Aasel_total += muts_Aa_sel
        muts_AAneu_total += muts_AA_neu
        muts_Aaneu_total += muts_Aa_neu
        muts_aasel_total += length(sel_loci) - muts_AA_sel - muts_Aa_sel
        muts_aaneu_total += length(loci) - length(sel_loci) - muts_AA_neu - muts_Aa_neu
    end

    muts_AAsel_total /= len
    muts_Aasel_total /= len
    muts_aasel_total /= len
    muts_AAneu_total /= len
    muts_Aaneu_total /= len
    muts_aaneu_total /= len
    return muts_AAsel_total, muts_Aasel_total, muts_aasel_total, muts_AAneu_total, muts_Aaneu_total, muts_aaneu_total, fits
end

function calc_muts_and_fitn_in_deme(deme_ms1, deme_ms2, wld_stats)
    calc_muts_and_fitn_in_deme(deme_ms1, deme_ms2, wld_stats["domin_coef"], wld_stats["loci"], wld_stats["sel_loci"])
end

"""
Randomly mutates at selected loci. Used when building the next generation in finite-sites expansions.

---

`ms1`: an individual's left monosome array

`ms2`: an individual's right monosome array

`mut_rate`: genome-wide mutation rate

`n_loci`: number of loci
"""
function mutate(ms1, ms2, mut_rate, n_loci; mutratelocus=false)
    get_mutation_random = mutratelocus ? rand(Poisson(mut_rate*n_loci)) : rand(Poisson(mut_rate))
    @fastmath @inbounds for _ in 1:get_mutation_random
        pos_alter = sample(1:n_loci)

        if rand(1:2) == 1
            ms1[pos_alter] = true
        else
            ms2[pos_alter] = true
        end
    end
end

function mutate(ms1, ms2, wld_stats; mutratelocus=false)
    mutate(ms1, ms2, wld_stats["mut_rate"], wld_stats["n_loci"]; mutratelocus=mutratelocus)
end

"""
Recombines loci with a 1/2 chance. Used when building the next generation in finite-sites expansions.

---

`ms1`: an individual's left monosome array

`ms2`: an individual's right monosome array

`mut_rate`: genome-wide mutation rate

`n_loci`: number of loci
"""
function crossover(ms1, ms2, n_loci)
    for j in 1:n_loci
        lr = rand(1:2)
        ms1[j] = lr == 1 ? ms1[j] : ms2[j]
    end
end

"""
Creates a zygote from two individuals. Used when building the next generation in finite-sites expansions.

---

`ind1`: individual 1's left OR right monosome array

`ind2`: individual 2's left OR right monosome array

`mut_rate`: genome-wide mutation rate

`n_loci`: number of loci

---

Output: left OR right monosome array of a zygote
"""
function mate(ind1, ind2, n_loci)
    new_loci = vcat(ind1[1:n_loci], ind2[1:n_loci])
    return new_loci
end

"""
Builds the next generation in finite-sites expansions, i.e. advances two world arrays (left and right monosomes) by one generation and returns the new generation data for fitness, populations, mutation numbers.

---

`wld_ms1`: a spatial array of demes that contain individuals' left monosome [Bool] arrays

`wld_ms2`: a spatial array of demes that contain individuals' right monosome [Bool] arrays

`wld_stats`: world stats Dict

`fitn_out`: if **true**, the new generation data for fitness will be output

`pops_out`: if **true**, the new generation data for populations will be output

`sel_out`: if **true**, the new generation data for selected mutations will be output

`neu_out`: if **true**, the new generation data for neutral mutations will be output

`max_migr`: a tuple of maximum migration area coordinates

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`refl_walls`: if **true**, walls reflect migrants

`r_max_migr`: Int maximum migration radius. If *>0**, migration is kept within this radius. Can be used in addition to `max_migr`

`r_coords`: a tuple (array) of axes' ordinal numbers that the n-sphere with `r_max_migr` covers. For example:
- **(1,3)** - migration is bound within a disk at x and z axes
- **(1,2,3)** - migration is bound within a sphere at x, y and z axes

---

Output 1: a changed `wld_ms1` = a spatial array of demes that contain individuals' left monosome [Bool] arrays

Output 2: a changed `wld_ms2` = a spatial array of demes that contain individuals' right monosome [Bool] arrays

Output 3: a spatial array of demes with average fitness in the new generation

Output 4: a spatial array of demes with populations in the new generation

Output 5: a spatial array of demes with average selected AA mutation count in the new generation

Output 6: a spatial array of demes with average selected Aa mutation count in the new generation

Output 7: a spatial array of demes with average selected aa mutation count in the new generation

Output 8: a spatial array of demes with average neutral AA mutation count in the new generation

Output 9: a spatial array of demes with average neutral Aa mutation count in the new generation

Output 10: a spatial array of demes with average neutral aa mutation count in the new generation
"""
function build_next_gen(wld_ms1, wld_ms2, wld_stats, fitn_out=false, pops_out=false, sel_out=false, neu_out=false;
    max_migr=NaN, migr_mode=DEF_MIGR_MODE, bottleneck=NaN, refl_walls=false, r_max_migr=0, r_coords=[1, 2], mutratelocus=false)

    wlddim = wld_stats["wlddim"]

    # Determine the number of offspring for each deme
    next_gen_posits, next_gen_pops = calc_offspring(wld_ms1, wld_stats)

    # Define the world (as an array [=demes] of arrays [=individs] of two Bool arrays [=monosomes]) and the data arrays in the next generation
    wld_ms1_next = Array{Array{Array{Bool}},wlddim}(undef, wld_stats["max"]...)
    wld_ms2_next = Array{Array{Array{Bool}},wlddim}(undef, wld_stats["max"]...)
    mean_fitn_next = NaN
    pops_next = NaN
    muts_AAsel_next = NaN
    muts_Aasel_next = NaN
    muts_aasel_next = NaN
    muts_AAneu_next = NaN
    muts_Aaneu_next = NaN
    muts_aaneu_next = NaN
    all_birth_count = 0

    # Fill the next generation habitat
    if fitn_out
        mean_fitn_next = Array{Float32}(undef, wld_stats["max"]...)
        fill!(mean_fitn_next, NaN)
    end
    if pops_out
        pops_next = Array{Float32}(undef, wld_stats["max"]...)
        fill!(pops_next, NaN)
    end
    if sel_out
        muts_AAsel_next = Array{Float32}(undef, wld_stats["max"]...)
        muts_Aasel_next = Array{Float32}(undef, wld_stats["max"]...)
        muts_aasel_next = Array{Float32}(undef, wld_stats["max"]...)
        fill!(muts_AAsel_next, NaN)
        fill!(muts_Aasel_next, NaN)
        fill!(muts_aasel_next, NaN)
    end
    if neu_out
        muts_AAneu_next = Array{Float32}(undef, wld_stats["max"]...)
        muts_Aaneu_next = Array{Float32}(undef, wld_stats["max"]...)
        muts_aaneu_next = Array{Float32}(undef, wld_stats["max"]...)
        fill!(muts_AAneu_next, NaN)
        fill!(muts_Aaneu_next, NaN)
        fill!(muts_aaneu_next, NaN)
    end

    for deme in next_gen_posits
        ms1_at_pos = wld_ms1[deme...]
        ms2_at_pos = wld_ms2[deme...]

        fitns = []
        cnt_res_AAsel, cnt_res_Aasel, cnt_res_aasel, cnt_res_AAneu, cnt_res_Aaneu, cnt_res_aaneu, fitns =
            calc_muts_and_fitn_in_deme(ms1_at_pos, ms2_at_pos, wld_stats)

        if fitn_out
            mean_fitn_next[deme...] = mean(fitns)
        end
        if sel_out
            muts_AAsel_next[deme...] = cnt_res_AAsel
            muts_Aasel_next[deme...] = cnt_res_Aasel
            muts_aasel_next[deme...] = cnt_res_aasel
        end
        if neu_out
            muts_AAneu_next[deme...] = cnt_res_AAneu
            muts_Aaneu_next[deme...] = cnt_res_Aaneu
            muts_aaneu_next[deme...] = cnt_res_aaneu
        end

        next_generation_size = next_gen_pops[deme...]

        if next_generation_size > 0
            birth_count = 0
            for _ in 1:next_generation_size
                
                mom_ms1 = wsample(ms1_at_pos, Float32.(fitns))
                mom_ms2 = wsample(ms2_at_pos, Float32.(fitns))
                dad_ms1 = wsample(ms1_at_pos, Float32.(fitns))
                dad_ms2 = wsample(ms2_at_pos, Float32.(fitns))

                gamete_mom_ms1 = copy(mom_ms1)
                gamete_dad_ms1 = copy(dad_ms1)
                gamete_mom_ms2 = copy(mom_ms2)
                gamete_dad_ms2 = copy(dad_ms2)

                crossover(gamete_mom_ms1, gamete_mom_ms2, wld_stats["n_loci"])
                crossover(gamete_dad_ms1, gamete_dad_ms2, wld_stats["n_loci"])
                mutate(gamete_mom_ms1, gamete_mom_ms2, wld_stats; mutratelocus=mutratelocus)
                mutate(gamete_dad_ms1, gamete_dad_ms2, wld_stats; mutratelocus=mutratelocus)

                move = calc_migr_dist(deme, wld_stats, migr_mode, bottleneck, max_migr, refl_walls, r_max_migr, r_coords)

                indices = [deme[i] + move[i] for i in 1:wlddim]

                if !isassigned(wld_ms1_next, indices...)
                    wld_ms1_next[indices...] = []
                    wld_ms2_next[indices...] = []
                end
                push!(wld_ms1_next[indices...], gamete_mom_ms1)
                push!(wld_ms2_next[indices...], gamete_dad_ms2)

                birth_count += 1
                all_birth_count += 1
            end

            if pops_out
                pops_next[deme...] = birth_count
            end
        end
    end

    return wld_ms1_next, wld_ms2_next, mean_fitn_next, pops_next, muts_AAsel_next, muts_Aasel_next, muts_aasel_next, muts_AAneu_next, muts_Aaneu_next, muts_aaneu_next
end

"""
Creates an empty deme space, i.e. an N-dimensional lattice of demes that can house individuals with a finite-site genetic structure.
N is determined from the dimensionality of the `max` tuple (2-dimensional by default).

---

`max`: world extents

`name`: world name

`capacity`: capacity of each deme

`prolif_rate`: proliferation rate

`mut_rate`: indiv.genome-wide mutation rate per generation

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inversely proportional Buffon-Laplace

`sel_coef`: default selection coefficient

`domin_coef`: dominance coefficient (in heterozygous loci, new_fitness *= **1 -** `domin_coef` * `sel_coef`)

`prop_of_del_muts`: proportion of deleterious mutations in nature

`n_loci`: number of loci in each individual

`n_sel_loci`: number of selected loci in each individual

`loci`: array of selected coefficients for every locus

---

Output 1: a spatial array of demes that contain individuals' left monosome [Bool] arrays (all empty)

Output 2: a spatial array of demes that contain individuals' right monosome [Bool] arrays (all empty)

Output 3: world stats Dict

"""
function create_empty_world(max::Tuple=(DEF_X_MAX, DEF_Y_MAX); name::String=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), capacity=DEF_CAPACITY,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, migr_mode=DEF_MIGR_MODE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS,
    prolif_rate=DEF_PROLIF_RATE, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci))

    wld_ms1 = Array{Array{Array{Bool}}}(undef, max...) # array of left (in a pair) monosomes ("ms") of all individuals in space
    wld_ms2 = Array{Array{Array{Bool}}}(undef, max...) # array of right (in a pair) monosomes ("ms") of all individuals in space
    for k in Iterators.product([1:n for n in max]...)
        wld_ms1[k...] = Array{Bool,1}[]
        wld_ms2[k...] = Array{Bool,1}[]
    end

    wld_stats = Dict(
        "name" => name,
        #"min" => min,
        "max" => max,
        "capacity" => capacity,
        "prolif_rate" => prolif_rate,
        "loci" => loci,
        "n_loci" => n_loci,
        "n_sel_loci" => n_sel_loci,
        "mut_rate" => mut_rate,
        "migr_rate" => migr_rate,
        "migr_mode" => migr_mode,
        "sel_coef" => sel_coef,
        "domin_coef" => domin_coef,
        "prop_of_del_muts" => prop_of_del_muts,
        "wlddim" => length(max)
        #"rangeexps" => []
    )

    return wld_ms1, wld_ms2, wld_stats
end

"""
Fills random demes within given monosome arrays with finite-sites individuals. Usually used after an empty world is created.

---

`wld_ms1`: a spatial array of demes that contain individuals' left monosome [Bool] arrays

`wld_ms2`: a spatial array of demes that contain individuals' right monosome [Bool] arrays

`wld_stats`: world stats Dict

`fill`: an array of Int ranges of the coordinates that define the area within which to fill

`n_demes_to_fill`: number of demes to fill

"""
function fill_random_demes(wld_ms1, wld_ms2, wld_stats, fill::Vector{UnitRange{Int64}}, n_demes_to_fill=DEF_N_DEMES_STARTFILL) # ::Array{Array{Array{Bool}}}

    possible_init_coords = [collect(x) for x in Iterators.product(fill...)]
    init_coords = sample(possible_init_coords, n_demes_to_fill; replace=false)
    wld_stats["sel_loci"] = randperm(wld_stats["n_loci"])[1:wld_stats["n_sel_loci"]]

    for coord in init_coords
        if !isassigned(wld_ms1, coord...)
            wld_ms1[coord...] = []
            wld_ms2[coord...] = []
        end
        for _ in 1:wld_stats["capacity"]
            push!(wld_ms1[coord...], falses(wld_stats["n_loci"]))
            push!(wld_ms2[coord...], falses(wld_stats["n_loci"]))
        end
    end

    wld_stats["startfill"] = copy(fill)
    wld_stats["n_demes_startfill"] = n_demes_to_fill
end

"""
Simulates a range expansion `n_re` times.
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

---

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium

`n_gens_exp`: duration of the expansion

`n_re`: number of replicates

`max_burnin`: a tuple of maximum coordinates during burn-in

`max_exp`: a tuple of maximum coordinates during expansion

`max`: a tuple of maximum coordinates of space

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**fitn**)
- **P** - deme populations (**pops**)
- **S** - deme-average number of homo- and heterozygous selected loci (**AAsel**, **Aasel** and **aasel**)
- **N** - deme-average number of homo- and heterozygous neutral loci (**AAneu**, **Aaneu** and **aaneu**)

`name`: world name

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`r_max_burnin`: radius that bounds the burn-in area

`r_max_exp`: radius that bounds the expansion area

`r_coords`: a tuple (array) of axes' ordinal numbers that the n-sphere with `r_max_migr` covers. For example:
- **(1,3)** - migration is bound within a disk at x and z axes
- **(1,2,3)** - migration is bound within a sphere at x, y and z axes

`loci`: array of selected coefficients for every locus

`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start

`distributed`: if **true**, distribute to threads

If starting from existing world, also provide:

`wld_ms1`: a spatial array of demes that contain individuals, each of which is a Bool array representing left monosomes

`wld_ms2`: a spatial array of demes that contain individuals, each of which is a Bool array representing right monosomes

`wld_stats`: world stats Dict

---

Output: a Dict containing data after the expansion:
- **stats** - statistics array containing world and range expansion information
- **fitn**, **pops**, **AAsel**, **Aasel**, **aasel**, **AAneu**, **Aaneu**, **aaneu** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`
"""
function rangeexp(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; max_burnin=(DEF_X_MAX_BURNIN, DEF_Y_MAX), max_exp=(DEF_X_MAX_EXP, DEF_Y_MAX), max=(DEF_X_MAX, DEF_Y_MAX), migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, r_max_burnin=0, r_max_exp=0, r_coords=[1, 2], capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci), mutratelocus=false,
    startfill_range=NaN, distributed=true, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN)

    fitn_wld = Array{Float32}(undef,0)
    pops_wld = Array{Float32}(undef,0)
    muts_AAsel_wld = Array{Float32}(undef,0)
    muts_Aasel_wld = Array{Float32}(undef,0)
    muts_aasel_wld = Array{Float32}(undef,0)
    muts_AAneu_wld = Array{Float32}(undef,0)
    muts_Aaneu_wld = Array{Float32}(undef,0)
    muts_aaneu_wld = Array{Float32}(undef,0)

    #= if n_re>1
        procs = addprocs(n_re)
        for k in 1:length(procs)
            @spawnat procs[k] include("../resk.jl")
        end
    end =#
    
    if !(wld_ms1 isa Array{Array{Array{Bool}}})
        #println("No world provided. Creating a new world.")
        wld_ms1, wld_ms2, wld_stats = create_empty_world(max; name=name, capacity=capacity, prolif_rate=prolif_rate, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci,
            mut_rate=mut_rate, migr_rate=migr_rate, migr_mode=migr_mode, sel_coef=sel_coef, domin_coef=domin_coef)
        if !isa(startfill_range, Array) && !any(isnan, max_burnin)
            startfill_range = [1:upper for upper in max_burnin]
        end

        wld_ms1 = [copy(wld_ms1) for j in 1:n_re]
        wld_ms2 = [copy(wld_ms2) for j in 1:n_re]
        for j in 1:n_re
            fill_random_demes(wld_ms1[j], wld_ms2[j], wld_stats, startfill_range)
        end
    else
        wld_ms1 = [copy(wld_ms1) for j in 1:n_re]
        wld_ms2 = [copy(wld_ms2) for j in 1:n_re]
    end

    wlddim = wld_stats["wlddim"]
    n_gens = n_gens_burnin + n_gens_exp

    fitn_out = false
    pops_out = false
    sel_out = false
    neu_out = false
    if occursin("F", data_to_generate)
        fitn_out = true
        fitn_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
    end
    if occursin("P", data_to_generate)
        pops_out = true
        pops_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
    end
    if occursin("S", data_to_generate)
        sel_out = true
        muts_AAsel_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
        muts_Aasel_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
        muts_aasel_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
    end
    if occursin("N", data_to_generate)
        neu_out = true
        muts_AAneu_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
        muts_Aaneu_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
        muts_aaneu_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
    end


    if n_re>1 && distributed
        
        @sync begin
            @inbounds for j in 1:n_re
                @spawn begin
                    for g in 1:n_gens
                        if g <= n_gens_burnin
                            max_migr = max_burnin
                            r_max_migr = r_max_burnin
                        else
                            max_migr = max_exp
                            r_max_migr = r_max_exp
                        end
                        wld_ms1[j], wld_ms2[j], fitn_next, pops_next, muts_AAsel_next, muts_Aasel_next, muts_aasel_next, muts_AAneu_next,
                        muts_Aaneu_next, muts_aaneu_next = build_next_gen(wld_ms1[j], wld_ms2[j], wld_stats, fitn_out, pops_out, sel_out, neu_out;
                            max_migr=max_migr, migr_mode=migr_mode, bottleneck=bottleneck, r_max_migr=r_max_migr, r_coords=r_coords, mutratelocus=mutratelocus)
                        if fitn_out
                            fitn_wld[repeat([:],wlddim)...,g,j] = fitn_next
                        end
                        if pops_out
                            pops_wld[repeat([:],wlddim)...,g,j] = pops_next
                        end
                        if sel_out
                            muts_AAsel_wld[repeat([:],wlddim)...,g,j] = muts_AAsel_next
                            muts_Aasel_wld[repeat([:],wlddim)...,g,j] = muts_Aasel_next
                            muts_aasel_wld[repeat([:],wlddim)...,g,j] = muts_aasel_next
                        end
                        if neu_out
                            muts_AAneu_wld[repeat([:],wlddim)...,g,j] = muts_AAneu_next
                            muts_Aaneu_wld[repeat([:],wlddim)...,g,j] = muts_Aaneu_next
                            muts_aaneu_wld[repeat([:],wlddim)...,g,j] = muts_aaneu_next
                        end
                    end
                end
            end
        end
        #rmprocs(n_re)
    else
        @inbounds for j in 1:n_re, g in 1:n_gens
            if g <= n_gens_burnin
                max_migr = max_burnin
                r_max_migr = r_max_burnin
            else
                max_migr = max_exp
                r_max_migr = r_max_exp
            end
            wld_ms1[j], wld_ms2[j], fitn_next, pops_next, muts_AAsel_next, muts_Aasel_next, muts_aasel_next, muts_AAneu_next,
            muts_Aaneu_next, muts_aaneu_next = build_next_gen(wld_ms1[j], wld_ms2[j], wld_stats, fitn_out, pops_out, sel_out, neu_out;
                max_migr=max_migr, migr_mode=migr_mode, bottleneck=bottleneck, r_max_migr=r_max_migr, r_coords=r_coords, mutratelocus=mutratelocus)
            if fitn_out
                fitn_wld[repeat([:],wlddim)...,g,j] = fitn_next
            end
            if pops_out
                pops_wld[repeat([:],wlddim)...,g,j] = pops_next
            end
            if sel_out
                muts_AAsel_wld[repeat([:],wlddim)...,g,j] = muts_AAsel_next
                muts_Aasel_wld[repeat([:],wlddim)...,g,j] = muts_Aasel_next
                muts_aasel_wld[repeat([:],wlddim)...,g,j] = muts_aasel_next
            end
            if neu_out
                muts_AAneu_wld[repeat([:],wlddim)...,g,j] = muts_AAneu_next
                muts_Aaneu_wld[repeat([:],wlddim)...,g,j] = muts_Aaneu_next
                muts_aaneu_wld[repeat([:],wlddim)...,g,j] = muts_aaneu_next
            end
        end
    end

    #= append!(wld_stats["rangeexps"],Dict(
        "x_max_burnin" => x_max_burnin,
        "y_max_burnin" => DEF_Y_MAX,
        "n_gens_burnin" => n_gens_burnin,
        "n_gens_exp" => n_gens_exp,
        "n_gens" => n_gens)) =#
    wld_stats["max_burnin"] = max_burnin
    wld_stats["max_exp"] = max_exp
    wld_stats["n_gens_burnin"] = n_gens_burnin
    wld_stats["n_gens_exp"] = n_gens_exp
    wld_stats["n_gens"] = n_gens

    #= res = []
    function needed_data(symb,arr,res)
        if occursin("P", data_to_generate)
            push!(res,Ref(arr))
        end
    end
    needed_data("P",pops_wld) =#

    return Dict("stats" => wld_stats, "fitn" => Array(fitn_wld), "pops" => Array(pops_wld), "AAsel" => Array(muts_AAsel_wld), "Aasel" => Array(muts_Aasel_wld),
        "aasel" => Array(muts_aasel_wld), "AAneu" => Array(muts_AAneu_wld), "Aaneu" => Array(muts_Aaneu_wld), "aaneu" => Array(muts_aaneu_wld))
end

"""
Simulates a range expansion `n_re` times in 1D, starting from one side of a segment space.
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

---

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium

`n_gens_exp`: duration of the expansion

`n_re`: number of replicates

`x_max_burnin`: the outward x-coordinate bound for migration during burn-in

`x_max_exp`: the outward x-coordinate bound for migration during the expansion

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**fitn**)
- **P** - deme populations (**pops**)
- **S** - deme-average number of homo- and heterozygous selected loci (**AAsel**, **Aasel** and **aasel**)
- **N** - deme-average number of homo- and heterozygous neutral loci (**AAneu**, **Aaneu** and **aaneu**)

`name`: world name

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`distributed`: if **true**, distribute to threads

`capacity`: capacity of each deme

If starting from existing world, also provide:

`wld_ms1`: a spatial array of demes that contain individuals, each of which is a Bool array representing left monosomes

`wld_ms2`: a spatial array of demes that contain individuals, each of which is a Bool array representing right monosomes

`wld_stats`: world stats Dict

---

Output: a Dict containing data after the expansion:
- **stats** - statistics array containing world and range expansion information
- **fitn**, **pops**, **AAsel**, **Aasel**, **aasel**, **AAneu**, **Aaneu**, **aaneu** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`
"""
function rangeexp_ray(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; x_max_burnin=DEF_X_MAX_BURNIN, x_max_exp=DEF_X_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, distributed=true, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN, capacity=DEF_CAPACITY)

    rangeexp(n_gens_burnin, n_gens_exp, n_re; max_burnin=(x_max_burnin,), max_exp=(x_max_exp,), max=(x_max_exp,), startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_ms1=wld_ms1, wld_ms2=wld_ms2, wld_stats=wld_stats, name=name, bottleneck=bottleneck, distributed=distributed)
end

function rangeexp_linear(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, max=(r_max_exp * 2 + 1,), 
    migr_mode=DEF_MIGR_MODE, startfill_range=[(1-ceil(Int,r_max_burnin/2)+r_max_exp):(1+ceil(Int,r_max_burnin/2)+r_max_exp)], prolif_rate=DEF_PROLIF_RATE, max_exp=NaN, max_burnin=NaN,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, distributed=true, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN, capacity=DEF_CAPACITY)

    rangeexp(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, max=max, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci, r_coords=[1],
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_ms1=wld_ms1, wld_ms2=wld_ms2, wld_stats=wld_stats, name=name, bottleneck=bottleneck, distributed=distributed)
end


"""
Simulates a 2D strip range expansion, in which a population expands in the positive x direction (after an optional burn-in phase).
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

---

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium

`n_gens_exp`: duration of the expansion

`x_max_burnin`: the outward x-coordinate bound for migration during burn-in

`x_max_exp`: the outward x-coordinate bound for migration during the expansion

`y_max`: the upper y-coordinate bound (lower bound is always **1** currently)

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**fitn**)
- **P** - deme populations (**pops**)
- **S** - deme-average number of homo- and heterozygous selected loci (**AAsel**, **Aasel** and **aasel**)
- **N** - deme-average number of homo- and heterozygous neutral loci (**AAneu**, **Aaneu** and **aaneu**)

`name`: world name

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`capacity`: capacity of each deme

If starting from existing world, also provide:

`wld_ms1`: a spatial array of demes that contain individuals, each of which is a Bool array representing left monosomes

`wld_ms2`: a spatial array of demes that contain individuals, each of which is a Bool array representing right monosomes

`wld_stats`: world stats Dict

You can also further specify the space aside from `x_max_burnin`, `x_max_exp` and `y_max`:

`max_burnin`: a tuple of maximum coordinates during burn-in

`max_exp`: a tuple of maximum coordinates during expansion

`max`: a tuple of maximum coordinates of space

---

Output: a Dict containing data after the expansion:
- **stats** - statistics array containing world and range expansion information
- **fitn**, **pops**, **AAsel**, **Aasel**, **aasel**, **AAneu**, **Aaneu**, **aaneu** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`
"""
function rangeexp_strip(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; x_max_burnin=DEF_X_MAX_BURNIN, x_max_exp=DEF_X_MAX_EXP, y_max=DEF_Y_MAX, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=("midhole at x=", x_max_burnin * 2), prolif_rate=DEF_PROLIF_RATE,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, capacity=DEF_CAPACITY, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN, max_burnin=(x_max_burnin, y_max), max_exp=(x_max_exp, y_max), max=(x_max_exp, y_max))

    rangeexp(n_gens_burnin, n_gens_exp, n_re; max_burnin=max_burnin, max_exp=max_exp, max=max, startfill_range=startfill_range, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, capacity=capacity, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_ms1=wld_ms1, wld_ms2=wld_ms2, wld_stats=wld_stats, name=name, bottleneck=bottleneck)
end

"""
Simulates a range expansion, in which a population expands from the center of a 2D disk (after an optional burn-in phase).
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

---

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium

`n_gens_exp`: duration of the expansion

`r_max_burnin`: radius that bounds the burn-in area

`r_max_exp`: radius that bounds the expansion area

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**fitn**)
- **P** - deme populations (**pops**)
- **S** - deme-average number of homo- and heterozygous selected loci (**AAsel**, **Aasel** and **aasel**)
- **N** - deme-average number of homo- and heterozygous neutral loci (**AAneu**, **Aaneu** and **aaneu**)

`name`: world name

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`capacity`: capacity of each deme

If starting from existing world, also provide:

`wld_ms1`: a spatial array of demes that contain individuals, each of which is a Bool array representing left monosomes

`wld_ms2`: a spatial array of demes that contain individuals, each of which is a Bool array representing right monosomes

`wld_stats`: world stats Dict

You can also further specify the space aside from `x_max_burnin`, `x_max_exp` and `y_max`:

`max_burnin`: a tuple of maximum coordinates during burn-in

`max_exp`: a tuple of maximum coordinates during expansion

`max`: a tuple of maximum coordinates of space

---

Output: a Dict containing data after the expansion:
- **stats** - statistics array containing world and range expansion information
- **fitn**, **pops**, **AAsel**, **Aasel**, **aasel**, **AAneu**, **Aaneu**, **aaneu** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`
"""
function rangeexp_disk(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, max=(r_max_exp * 2 + 1, r_max_exp * 2 + 1), capacity=DEF_CAPACITY,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    max_exp=NaN, max_burnin=NaN, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN)

    if !isa(startfill_range, Array)
        ran = ins_sq(r_max_burnin, r_max_exp)
        startfill_range = [ran, ran]
    end

    rangeexp(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, max=max, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_ms1=wld_ms1, wld_ms2=wld_ms2, wld_stats=wld_stats, name=name, bottleneck=bottleneck)
end

function rangeexp_cylinder(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE,
    z_max_burnin=DEF_X_MAX_BURNIN, z_max_exp=DEF_X_MAX_EXP, max_burnin=(NaN, NaN, z_max_burnin), max_exp=(NaN, NaN, z_max_exp), max=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, z_max_exp), capacity=DEF_CAPACITY,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    data_to_generate=DEF_DATA_TO_GENERATE, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN)

    if !isa(startfill_range, Array)
        ran = ins_sq(r_max_burnin, r_max_exp)
        startfill_range = [ran, ran, 1:z_max_burnin]
    end

    rangeexp(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, max=max,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_ms1=wld_ms1, wld_ms2=wld_ms2, wld_stats=wld_stats, name=name, bottleneck=bottleneck, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci,
        startfill_range=startfill_range)
end

function rangeexp_sphere(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,
    max_burnin=NaN, max_exp=NaN, max=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, r_max_exp * 2 + 1), capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    data_to_generate=DEF_DATA_TO_GENERATE, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN)

    if !isa(startfill_range, Array)
        ran = ins_cb(r_max_burnin, r_max_exp)
        startfill_range = [ran, ran, ran]
    end

    rangeexp(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, max=max, r_coords=[1, 2, 3], prolif_rate=prolif_rate,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_ms1=wld_ms1, wld_ms2=wld_ms2, wld_stats=wld_stats, name=name, bottleneck=bottleneck, capacity=capacity,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci,
        startfill_range=startfill_range)
end



# Simulation functions (infinite-sites)
# ------------------------------------------------

"""
Randomly adds mutations and calculates the number of mutations in an individual. Used when building the next generation in infinite-sites expansions.

---

`person`: an individual's segr. regions (fitness) array

`mut_rate`: genome-wide mutation rate

`n_segr_regions`: number of segregating regions

`sel_coef`: selection coefficient

---

Output 1: number of deleterious mutations

Output 2: number of beneficial mutations

"""
function mutate_inf(person, mut_rate, n_segr_regions, sel_coef, prop_of_del_muts)
    muts_del = 0
    muts_ben = 0

    get_mutation_random = rand(Poisson(mut_rate))
    
    @fastmath @inbounds for _ in 1:get_mutation_random
        pos_alter = sample(1:n_segr_regions)
        if rand() < prop_of_del_muts
            person[pos_alter] *= 1 - sel_coef
            muts_del += 1
        else
            person[pos_alter] *= 1 + sel_coef
            muts_ben += 1
        end
    end

    return muts_del, muts_ben
end

function mutate_inf(person, wld_stats)
    mutate_inf(person, wld_stats["mut_rate"], wld_stats["n_segr_regions"], wld_stats["sel_coef"], wld_stats["prop_of_del_muts"])
end

"""
Recombines loci with a 1/2 chance. Used when building the next generation in infinite-sites expansions.

---

`person`: an individual's segr. regions (fitness) array

`n_segr_regions`: number of segregating regions
"""
function crossover_inf(person, n_segr_regions::Int)
    for i in 1:n_segr_regions
        lr = rand(1:2)
        person[i] = lr == 1 ? person[i] : person[i+n_segr_regions]
    end
end

function crossover_inf(person, wld_stats::Dict)
    crossover_inf(person, wld_stats["n_segr_regions"])
end

"""
Creates a zygote from two individuals. Used when building the next generation in infinite-sites expansions.

---

`ind1`: individual 1's segr. regions (fitness) array

`ind2`: individual 2's segr. regions (fitness) array

`n_segr_regions`: number of segregating regions

---

Output: segr. regions (fitness) array of a zygote

---
"""
function mate_inf(ind1, ind2, n_segr_regions; fixed_mate=false)
    lr1 = fixed_mate || (rand(1:2) == 1) ? (1:n_segr_regions) : ((n_segr_regions+1):(n_segr_regions*2))
    lr2 = fixed_mate || (rand(1:2) == 1) ? (1:n_segr_regions) : ((n_segr_regions+1):(n_segr_regions*2))
    return vcat(ind1[lr1], ind2[lr2])
end

function mate_cond(mom_fit,dad_fit,max_fitness)
    return (mom_fit > rand()*max_fitness) & (dad_fit > rand()*max_fitness)
end

"""
Builds the next generation in infinite-sites expansions, i.e. advances the world array (segr. regions' fitness array) by one generation and returns the new generation data for fitness, populations, mutation numbers.

---

`wld`: a spatial array of demes that contain individuals' segr. regions (fitness) [Float] arrays

`wld_stats`: world stats Dict

`fitn_out`: if **true**, the new generation data for fitness will be output

`pops_out`: if **true**, the new generation data for populations will be output

`mut_out`: if **true**, the new generation data for mutation counts (deleterious & beneficial) will be output

`max_migr`: a tuple of maximum migration area coordinates

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`refl_walls`: if **true**, walls reflect migrants

`r_max_migr`: Int maximum migration radius. If *>0**, migration is kept within this radius. Can be used in addition to `max_migr`

`r_coords`: a tuple (array) of axes' ordinal numbers that the n-sphere with `r_max_migr` covers. For example:
- **(1,3)** - migration is bound within a disk at x and z axes
- **(1,2,3)** - migration is bound within a sphere at x, y and z axes

---

Output 1: a changed `wld` = a spatial array of demes that contain individuals' segr. regions [Float] arrays

Output 2: a spatial array of demes with average fitness in the new generation

Output 3: a spatial array of demes with populations in the new generation

Output 4: a spatial array of demes with average deleterious mutation count in the new generation

Output 5: a spatial array of demes with average beneficial mutation count in the new generation
"""
function build_next_gen_inf(wld::Array{Array{Array{Float32}}}, wld_stats, fitn_out=false, pops_out=false, mut_out=false;
    max_migr=NaN, migr_mode=DEF_MIGR_MODE, bottleneck=NaN, refl_walls=false, r_max_migr=0, r_coords=[1, 2], weightfitn=true, condsel=false, fixed_mate=false, premutate=false)

    wlddim = wld_stats["wlddim"]

    # Determine the number of offspring for each deme
    next_gen_posits, next_gen_pops = calc_offspring(wld, wld_stats)
    
    # Define the habitat (world) and the data arrays in the next generation
    wld_next = Array{Array{Array{Float32}},wlddim}(undef, wld_stats["max"]...)
    mean_fitn_next = NaN
    pops_next = NaN
    muts_del_next = NaN
    muts_ben_next = NaN
    all_birth_count = 0

    # Fill the next generation habitat
    if fitn_out
        mean_fitn_next = Array{Float32}(undef, wld_stats["max"]...)
        fill!(mean_fitn_next, NaN)
    end
    if pops_out
        pops_next = Array{Float32}(undef, wld_stats["max"]...)
        fill!(pops_next, NaN)
    end
    if mut_out
        muts_del_next = Array{Float32}(undef, wld_stats["max"]...)
        muts_ben_next = Array{Float32}(undef, wld_stats["max"]...)
        fill!(muts_del_next, NaN)
        fill!(muts_ben_next, NaN)
    end


    for deme in next_gen_posits
        inds_at_pos = wld[deme...]
        fitns = prod.(inds_at_pos)

        if fitn_out
            mean_fitn_next[deme...] = mean(fitns)
        end

        next_generation_size = next_gen_pops[deme...]

        if next_generation_size > 0
            birth_count = 0
            
            while birth_count < next_generation_size
                mom = weightfitn ? wsample(inds_at_pos, fitns) : sample(inds_at_pos)
                dad = weightfitn ? wsample(inds_at_pos, fitns) : sample(inds_at_pos)

                if !condsel || mate_cond(prod(mom),prod(dad),maximum(fitns))
                    
                    gamete_mom = copy(mom)
                    gamete_dad = copy(dad)

                    crossover_inf(gamete_mom, wld_stats["n_segr_regions"])
                    crossover_inf(gamete_dad, wld_stats["n_segr_regions"])
                    
                    if premutate
                        muts_del, muts_ben = mutate_inf(gamete_mom, wld_stats)
                        muts_del, muts_ben = mutate_inf(gamete_dad, wld_stats)
                        mate_result = mate_inf(gamete_mom, gamete_dad, wld_stats["n_segr_regions"]; fixed_mate=fixed_mate)
                        
                    else
                        mate_result = mate_inf(gamete_mom, gamete_dad, wld_stats["n_segr_regions"]; fixed_mate=fixed_mate)
                        muts_del, muts_ben = mutate_inf(mate_result, wld_stats)
                    end

                    if mut_out
                        if isnan(muts_del_next[deme...])
                            muts_del_next[deme...] = 0
                        end
                        if isnan(muts_ben_next[deme...])
                            muts_ben_next[deme...] = 0
                        end
                        muts_del_next[deme...] += muts_del
                        muts_ben_next[deme...] += muts_ben
                    end

                    move = calc_migr_dist(deme, wld_stats, migr_mode, bottleneck, max_migr, refl_walls, r_max_migr, r_coords)
                    
                    indices = [deme[i] + move[i] for i in 1:wlddim]
                    if !isassigned(wld_next, indices...)
                        wld_next[indices...] = []
                    end
                    push!(wld_next[indices...], mate_result)

                    birth_count += 1
                    all_birth_count += 1
                end
            end

            if pops_out
                pops_next[deme...] = birth_count
            end
        end
    end

    return wld_next, mean_fitn_next, pops_next, muts_del_next, muts_ben_next
end

"""
Creates an empty world (deme space) with infinite-sites individual structure. 2-dimensional by default.

---

`max`: a tuple of maximal space bounds (coordinates)

`min`: a tuple of minimal space bounds (coordinates). Limited to (1,1) for now

`name`: world name

`capacity`: capacity of each deme

`prolif_rate`: proliferation rate

`n_segr_regions`: number of segregating regions

`mut_rate`: genome-wide mutation rate

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`sel_coef`: selection coefficient

`domin_coef`: dominance coefficient (in heterozygous loci, new_fitness *= **1 -** `domin_coef` * `sel_coef`)

---

Output 1: a spatial array of demes that contain individuals' segr. regions [Float] arrays (all empty)

Output 2: world stats Dict

"""
function create_empty_world_inf(max=(DEF_X_MAX, DEF_Y_MAX); min=(1, 1), name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), capacity=DEF_CAPACITY,
    prolif_rate=DEF_PROLIF_RATE, n_segr_regions=DEF_N_SEGR_REGIONS,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, migr_mode=DEF_MIGR_MODE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)

    wld = Array{Array{Array{Float32}}}(undef, max...) # array of fitness values of all individuals in space
    for k in Iterators.product([1:n for n in max]...)
        wld[k...] = Array{Float32,1}[]
    end

    wld_stats = Dict(
        "name" => name,
        "max" => max,
        "capacity" => capacity,
        "prolif_rate" => prolif_rate,
        "n_segr_regions" => n_segr_regions,
        "mut_rate" => mut_rate,
        "migr_rate" => migr_rate,
        "migr_mode" => migr_mode,
        "sel_coef" => sel_coef,
        "prop_of_del_muts" => prop_of_del_muts,
        "wlddim" => length(max)
        #"rangeexps" => []
    )

    return wld, wld_stats
end

"""
Fills random demes within given monosome arrays with infinite-sites individuals. Usually used after an empty world is created.

---

`wld`: a spatial array of demes that contain individuals' segr. regions [Float] arrays

`wld_stats`: world stats Dict

`fill`: an array of Int ranges of the coordinates that define the area within which to fill

`n_demes_to_fill`: number of demes to fill

"""
function fill_random_demes_inf(wld::Array{Array{Array{Float32}}}, wld_stats, fill::Vector{UnitRange{Int64}}, n_demes_to_fill=DEF_N_DEMES_STARTFILL)

    possible_init_coords = [collect(x) for x in Iterators.product(fill...)]
    init_coords = sample(possible_init_coords, n_demes_to_fill; replace=false)

    for coord in init_coords
        if !isassigned(wld, coord...)
            wld[coord...] = []
        end
        for _ in 1:wld_stats["capacity"]
            push!(wld[coord...], ones(wld_stats["n_segr_regions"] * 2))
        end
    end

    wld_stats["startfill"] = copy(fill)
    wld_stats["n_demes_startfill"] = n_demes_to_fill
end

"""

Simulates a range expansion with infinite-sites individuals `n_re` times.
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

---

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium

`n_gens_exp`: duration of the expansion

`n_re`: number of replicates

`max_burnin`: a tuple of maximum coordinates during burn-in

`max_exp`: a tuple of maximum coordinates during expansion

`max`: a tuple of maximum coordinates of space

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**fitn**)
- **P** - deme populations (**pops**)
- **M** - deme-average number of deleterious and beneficial mutations (**del** and **ben**)

`name`: world name

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`r_max_burnin`: radius that bounds the burn-in area

`r_max_exp`: radius that bounds the expansion area

`r_coords`: a tuple (array) of axes' ordinal numbers that the n-sphere with `r_max_migr` covers. For example:
- **(1,3)** - migration is bound within a disk at x and z axes
- **(1,2,3)** - migration is bound within a sphere at x, y and z axes

`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start

If starting from existing world, also provide:

`wld`: a spatial array of demes that contain individuals, each of which is a Float array of fitness values

`wld_stats`: world stats Dict

---

Output: a Dict containing data after the expansion:
- **stats** - statistics array containing world and range expansion information
- **fitn**, **pops**, **del**, **ben** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`
"""
function rangeexp_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; max_burnin=(DEF_X_MAX_BURNIN, DEF_Y_MAX), max_exp=(DEF_X_MAX_EXP, DEF_Y_MAX), max=(DEF_X_MAX, DEF_Y_MAX), migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, r_max_burnin=0, r_max_exp=0, capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, 
    n_segr_regions=DEF_N_SEGR_REGIONS,
    distributed=true, weightfitn=true, condsel=false, fixed_mate=false, premutate=false,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, r_coords=[1, 2], startfill_range=NaN, wld=NaN, wld_stats=NaN)

    fitn_wld = Array{Float32}(undef,0)
    pops_wld = Array{Float32}(undef,0)
    muts_del_wld = Array{Float32}(undef,0)
    muts_ben_wld = Array{Float32}(undef,0)

    if !(wld isa Array{Array{Array{Float32}}})
        #println("No world provided. Creating a new world.")
        wld, wld_stats = create_empty_world_inf(max; name=name, capacity=capacity, prolif_rate=prolif_rate, mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts,
            n_segr_regions=n_segr_regions)
        if !isa(startfill_range, Array) && !any(isnan, max_burnin)
            startfill_range = [1:upper for upper in max_burnin]
        end
        wld = [copy(wld) for j in 1:n_re]
        for j in 1:n_re
            fill_random_demes_inf(wld[j], wld_stats, startfill_range)
        end
        
    else
        wld = [copy(wld) for j in 1:n_re]
    end

    wlddim = wld_stats["wlddim"]
    n_gens = n_gens_burnin + n_gens_exp

    fitn_out = false
    pops_out = false
    mut_out = false
    if occursin("F", data_to_generate)
        fitn_out = true
        fitn_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
    end
    if occursin("P", data_to_generate)
        pops_out = true
        pops_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
    end
    if occursin("M", data_to_generate)
        mut_out = true
        muts_del_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
        muts_ben_wld = SharedArray{Float32}(wld_stats["max"]..., n_gens, n_re)
    end

    if n_re>1 && distributed
        @sync begin
            @inbounds for j in 1:n_re
                @spawn begin
                    for g in 1:n_gens

                        if g <= n_gens_burnin
                            max_migr = max_burnin
                            r_max_migr = r_max_burnin
                        else
                            max_migr = max_exp
                            r_max_migr = r_max_exp
                        end

                        wld[j], fitn_next, pops_next, muts_del_next, muts_ben_next = build_next_gen_inf(wld[j], wld_stats, fitn_out, pops_out, mut_out;
                            max_migr=max_migr, migr_mode=migr_mode, bottleneck=bottleneck, r_max_migr=r_max_migr, r_coords=r_coords, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate)
                        
                        if fitn_out
                            fitn_wld[repeat([:],wlddim)...,g,j] = fitn_next
                        end
                        if pops_out
                            pops_wld[repeat([:],wlddim)...,g,j] = pops_next
                        end
                        if mut_out
                            muts_del_wld[repeat([:],wlddim)...,g,j] = muts_del_next
                            muts_ben_wld[repeat([:],wlddim)...,g,j] = muts_ben_next
                        end
                    end
                end
            end
        end
        #rmprocs(n_re)
    else
        @inbounds for j in 1:n_re, g in 1:n_gens
            
            if g <= n_gens_burnin
                max_migr = max_burnin
                r_max_migr = r_max_burnin
            else
                max_migr = max_exp
                r_max_migr = r_max_exp
            end
            #println("gen",g)
            
            wld[j], fitn_next, pops_next, muts_del_next, muts_ben_next = build_next_gen_inf(wld[j], wld_stats, fitn_out, pops_out, mut_out;
                max_migr=max_migr, migr_mode=migr_mode, bottleneck=bottleneck, r_max_migr=r_max_migr, r_coords=r_coords, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate)
            #println("B"*string(rand()))

            if fitn_out
                fitn_wld[repeat([:],wlddim)...,g,j] = fitn_next
            end
            if pops_out
                pops_wld[repeat([:],wlddim)...,g,j] = pops_next
            end
            if mut_out
                muts_del_wld[repeat([:],wlddim)...,g,j] = muts_del_next
                muts_ben_wld[repeat([:],wlddim)...,g,j] = muts_ben_next
            end
        end
    end

    wld_stats["max_burnin"] = max_burnin
    wld_stats["max_exp"] = max_exp
    wld_stats["n_gens_burnin"] = n_gens_burnin
    wld_stats["n_gens_exp"] = n_gens_exp
    wld_stats["n_gens"] = n_gens

    return Dict("stats" => wld_stats, "fitn" => Array(fitn_wld), "pops" => Array(pops_wld), "del" => Array(muts_del_wld), "ben" => Array(muts_ben_wld))
end


"""
Simulates a range expansion with infinite-sites individuals `n_re` times in 1D, starting from one side of a segment space.
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

---

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium

`n_gens_exp`: duration of the expansion

`n_re`: number of replicates

`x_max_burnin`: the outward x-coordinate bound for migration during burn-in

`x_max_exp`: the outward x-coordinate bound for migration during the expansion

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**fitn**)
- **P** - deme populations (**pops**)
- **M** - deme-average number of deleterious and beneficial mutations (**del** and **ben**)

`name`: world name

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`distributed`: if **true**, distribute to threads

`capacity`: capacity of each deme

If starting from existing world, also provide:

`wld`: a spatial array of demes that contain individuals, each of which is a Float array of fitness values

`wld_stats`: world stats Dict

---

Output: a Dict containing data after the expansion:
- **stats** - statistics array containing world and range expansion information
- **fitn**, **pops**, **del**, **ben** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`
"""
function rangeexp_ray_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; x_max_burnin=DEF_X_MAX_BURNIN, x_max_exp=DEF_X_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, premutate=false, n_segr_regions=DEF_N_SEGR_REGIONS,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, distributed=true, wld=NaN, wld_stats=NaN, capacity=DEF_CAPACITY)

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; max_burnin=(x_max_burnin,), max_exp=(x_max_exp,), max=(x_max_exp,), startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate, n_segr_regions=n_segr_regions,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld=wld, wld_stats=wld_stats, name=name, bottleneck=bottleneck, distributed=distributed)
end

function rangeexp_linear_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, max=(r_max_exp * 2 + 1,),
    migr_mode=DEF_MIGR_MODE, startfill_range=[(1-ceil(Int,r_max_burnin/2)+r_max_exp):(1+ceil(Int,r_max_burnin/2)+r_max_exp)], prolif_rate=DEF_PROLIF_RATE, max_exp=NaN, max_burnin=NaN,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, premutate=false, n_segr_regions=DEF_N_SEGR_REGIONS,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, distributed=true, wld=NaN, wld_stats=NaN, capacity=DEF_CAPACITY)

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, max=max, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate, n_segr_regions=n_segr_regions, 
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld=wld, wld_stats=wld_stats, name=name, bottleneck=bottleneck, distributed=distributed, r_coords=[1])
end

"""
Simulates a 2D strip range expansion with infinite-sites individuals, in which a population expands in the positive x direction (after an optional burn-in phase).
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

---

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium

`n_gens_exp`: duration of the expansion

`x_max_burnin`: the outward x-coordinate bound for migration during burn-in

`x_max_exp`: the outward x-coordinate bound for migration during the expansion

`y_max`: the upper y-coordinate bound (lower bound is always **1** currently)

`migr_mode`: mode of migration. Possible values:
- **ort** - orthogonal directions only
- **all** - orthogonal and diagonal
- **hex** - hexagonal grid
- **diag1/2** - orthogonal and half-weighted diagonal
- **buffon1** - equidistant Buffon-Laplace (see documentation)
- **buffon2** - uniform Buffon-Laplace
- **buffon3** - inv.proportional Buffon-Laplace

`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start

`data_to_generate`: string of letters representing different data to output. Possible values:
- **F** - deme-average fitness (**fitn**)
- **P** - deme populations (**pops**)
- **M** - deme-average number of deleterious and beneficial mutations (**del** and **ben**)

`name`: world name

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`capacity`: capacity of each deme

If starting from existing world, also provide:

`wld`: a spatial array of demes that contain individuals, each of which is a Float array of fitness values

`wld_stats`: world stats Dict

You can also further specify the space aside from `x_max_burnin`, `x_max_exp` and `y_max`:

`max_burnin`: a tuple of maximum coordinates during burn-in

`max_exp`: a tuple of maximum coordinates during expansion

`max`: a tuple of maximum coordinates of space

---

Output: a Dict containing data after the expansion:
- **stats** - statistics array containing world and range expansion information
- **fitn**, **pops**, **del**, **ben** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`
"""
function rangeexp_strip_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; x_max_burnin=DEF_X_MAX_BURNIN, x_max_exp=DEF_X_MAX_EXP, y_max=DEF_Y_MAX, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,
    max_burnin=(x_max_burnin, y_max), max_exp=(x_max_exp, y_max), max=(x_max_exp, y_max), capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, n_segr_regions=DEF_N_SEGR_REGIONS,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, premutate=false,
    data_to_generate=DEF_DATA_TO_GENERATE, wld=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=("midhole at x=", x_max_burnin * 2))

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; max_burnin=max_burnin, max_exp=max_exp, max=max, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate, n_segr_regions=n_segr_regions,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld=wld, wld_stats=wld_stats, name=name, bottleneck=bottleneck)
end

function rangeexp_disk_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, max=(r_max_exp * 2 + 1, r_max_exp * 2 + 1),
    capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false,
    premutate=false, n_segr_regions=DEF_N_SEGR_REGIONS,
    data_to_generate=DEF_DATA_TO_GENERATE, wld=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, max_exp=NaN, max_burnin=NaN, startfill_range=NaN)

    if !isa(startfill_range, Array)
        ran = ins_sq(r_max_burnin, r_max_exp)
        startfill_range = [ran, ran]
    end

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, max=max, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate, n_segr_regions=n_segr_regions,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld=wld, wld_stats=wld_stats, name=name, bottleneck=bottleneck)
end

function rangeexp_cylinder_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE,
    z_max_burnin=DEF_X_MAX_BURNIN, z_max_exp=DEF_X_MAX_EXP, max_burnin=(NaN, NaN, z_max_burnin), max_exp=(NaN, NaN, z_max_exp), max=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, z_max_exp), capacity=DEF_CAPACITY,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, premutate=false, n_segr_regions=DEF_N_SEGR_REGIONS,
    data_to_generate=DEF_DATA_TO_GENERATE, wld=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN)

    if !isa(startfill_range, Array)
        ran = ins_sq(r_max_burnin, r_max_exp)
        startfill_range = [ran, ran, 1:z_max_burnin]
    end

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, max=max,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld=wld, wld_stats=wld_stats, name=name, bottleneck=bottleneck, capacity=capacity, prolif_rate=prolif_rate, n_segr_regions=n_segr_regions,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate,
        startfill_range=startfill_range)
end

function rangeexp_sphere_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,
    max_burnin=NaN, max_exp=NaN, max=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, r_max_exp * 2 + 1), capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, n_segr_regions=DEF_N_SEGR_REGIONS,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, premutate=false,
    data_to_generate=DEF_DATA_TO_GENERATE, wld=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN)

    if !isa(startfill_range, Array)
        ran = ins_cb(r_max_burnin, r_max_exp)
        startfill_range = [ran, ran, ran]
    end

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, max=max, r_coords=[1, 2, 3],
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld=wld, wld_stats=wld_stats, name=name, bottleneck=bottleneck, capacity=capacity, prolif_rate=prolif_rate, n_segr_regions=n_segr_regions,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate,
        startfill_range=startfill_range)
end

# Analysis functions
# ------------------------------------------------

"""
Finds the average values of `data` over the whole population for each generation.

---

`data`: array with dimensions (space + time)

`n_gens`: number of generations

`dims`: number of dimensions of `data`

---

Output: array of averages of `data` for every generation
"""
function average_all(data::Array, n_gens::Int)
    res = Array{Float32}(undef, 0)
    for j in 1:n_gens
        push!(res, mean(filter(!isnan, data[repeat([:],length(size(data))-2)...,j,:])))
    end
    return res
end

"""
Finds the average values of `dataname` in `re` over the whole population for each generation.

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

---

Output: array of averages of `re[dataname]` for every generation
"""
function average_all(re::Dict, dataname::String)
    average_all(re[dataname], re["stats"]["n_gens"])
end

"""
Finds the average value of `data` between all demes at the expansion front.

---

`data`: array with dimensions (space + time)

`n_gens`: number of generations

`greaterzero`: if **true**, **>0** values are considered when determining the front (**>=0** values if **false**)

`oneside`: if **true**, approach only from one side (i.e. from the positive direction in strip expansions)

`divide`: if **false**, find the sum instead of average

---

Output: array of averages of `data` for every generation
"""
function average_front(data, n_gens, x_max; greaterzero=false, oneside=false, divide=true)
    n_re = size(data, 3)
    av_arr = Array{Float32}(undef, n_gens, n_re)

    for i in 1:n_re, j in 1:n_gens
        a_sum = 0
        cnt = 0
        frontier = x_max
        while frontier != 1 && (isnan(data[frontier, j, i]) || (greaterzero && data[frontier, j, i] == 0))
            frontier -= 1
        end
        if data[frontier, j, i] >= 0 || (greaterzero && data[frontier, j, i] > 0)
            a_sum += data[frontier, j, i]
            cnt += 1
        end
        if !oneside
            frontier = 1
            while frontier != x_max && (isnan(data[frontier, j, i]) || (greaterzero && data[frontier, j, i] == 0))
                frontier += 1
            end
            if data[frontier, j, i] >= 0 || (greaterzero && data[frontier, j, i] > 0)
                a_sum += data[frontier, j, i]
                cnt += 1
            end
        end
        if divide
            a_sum /= cnt
        end
        av_arr[j, i] = a_sum
    end

    return av_arr
end

function average_front(data, n_gens, x_max, y_max; greaterzero=false, oneside=false, divide=true)
    n_re = size(data, 4)
    av_arr = Array{Float32}(undef, n_gens, n_re)

    for i in 1:n_re, j in 1:n_gens
        a_sum = 0
        cnt = 0
        # scanning every y: side 1
        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && (isnan(data[frontier_x, _y, j, i]) || (greaterzero && data[frontier_x, _y, j, i] == 0))
                frontier_x -= 1
            end
            if data[frontier_x, _y, j, i] >= 0 || (greaterzero && data[frontier_x, _y, j, i] > 0)
                a_sum += data[frontier_x, _y, j, i]
                cnt += 1
            end
        end
        # scanning every y: side 2
        if !oneside
            for _y in 1:y_max
                frontier_x = 1
                while frontier_x != x_max && (isnan(data[frontier_x, _y, j, i]) || (greaterzero && data[frontier_x, _y, j, i] == 0))
                    frontier_x += 1
                end
                if data[frontier_x, _y, j, i] >= 0 || (greaterzero && data[frontier_x, _y, j, i] > 0)
                    a_sum += data[frontier_x, _y, j, i]
                    cnt += 1
                end
            end
        end
        mean_both_sides_y = a_sum
        if divide
            mean_both_sides_y /= cnt
        end

        if !oneside
            a_sum = 0
            cnt = 0
            # scanning every x: side 1
            for _x in 1:x_max
                frontier_y = y_max
                while frontier_y != 1 && (isnan(data[_x, frontier_y, j, i]) || (greaterzero && data[_x, frontier_y, j, i] == 0))
                    frontier_y -= 1
                end
                if data[_x, frontier_y, j, i] >= 0 || (greaterzero && data[_x, frontier_y, j, i] > 0)
                    a_sum += data[_x, frontier_y, j, i]
                    cnt += 1
                end
            end
            # scanning every x: side 2
            for _x in 1:x_max
                frontier_y = 1
                while frontier_y != y_max && (isnan(data[_x, frontier_y, j, i]) || (greaterzero && data[_x, frontier_y, j, i] == 0))
                    frontier_y += 1
                end
                if data[_x, frontier_y, j, i] >= 0 || (greaterzero && data[_x, frontier_y, j] > 0)
                    a_sum += data[_x, frontier_y, j, i]
                    cnt += 1
                end
            end
            if divide
                mean_both_sides_x = a_sum / cnt
            end
            av_arr[j, i] = (mean_both_sides_x + mean_both_sides_y) / 2
        else
            av_arr[j, i] = mean_both_sides_y
        end
    end
    return av_arr
end

function average_front(data, n_gens, x_max, y_max, z_max; greaterzero=false, oneside=false, divide=true)
    n_re = size(data, 5)
    av_arr = Array{Float32}(undef, n_gens, n_re)

    for i in 1:n_re, j in 1:n_gens
        a_sum = 0
        cnt = 0
        # scanning every xy: side 1
        for _x in 1:x_max, _y in 1:y_max
            frontier_z = z_max
            while frontier_z != 1 && (isnan(data[_x, _y, frontier_z, j, i]) || (greaterzero && data[_x, _y, frontier_z, j, i] == 0))
                frontier_z -= 1
            end
            if data[_x, _y, frontier_z, j, i] >= 0 || (greaterzero && data[_x, _y, frontier_z, j, i] > 0)
                a_sum += data[_x, _y, frontier_z, j, i]
                cnt += 1
            end
        end
        # scanning every xy: side 2
        if !oneside
            for _x in 1:x_max, _y in 1:y_max
                frontier_z = 1
                while frontier_z != z_max && (isnan(data[_x, _y, frontier_z, j, i]) || (greaterzero && data[_x, _y, frontier_z, j, i] == 0))
                    frontier_z += 1
                end
                if data[_x, _y, frontier_z, j, i] >= 0 || (greaterzero && data[_x, _y, frontier_z, j, i] > 0)
                    a_sum += data[_x, _y, frontier_z, j, i]
                    cnt += 1
                end
            end
        end
        mean_both_sides_xy = a_sum
        if divide
            mean_both_sides_xy /= cnt
        end

        if !oneside
            a_sum = 0
            cnt = 0
            # scanning every yz: side 1
            for _y in 1:y_max, _z in 1:z_max
                frontier_x = x_max
                while frontier_x != 1 && (isnan(data[frontier_x, _y, _z, j, i]) || (greaterzero && data[frontier_x, _y, _z, j, i]== 0))
                    frontier_x -= 1
                end
                if data[frontier_x, _y, _z, j, i] >= 0 || (greaterzero && data[frontier_x, _y, _z, j, i] > 0)
                    a_sum += data[frontier_x, _y, _z, j, i]
                    cnt += 1
                end
            end
            # scanning every yz: side 2
            for _y in 1:y_max, _z in 1:z_max
                frontier_x = 1
                while frontier_x != x_max && (isnan(data[frontier_x, _y, _z, j, i]) || (greaterzero && data[frontier_x, _y, _z, j, i] == 0))
                    frontier_x += 1
                end
                if data[frontier_x, _y, _z, j, i] >= 0 || (greaterzero && data[frontier_x, _y, _z, j, i] > 0)
                    a_sum += data[frontier_x, _y, _z, j, i]
                    cnt += 1
                end
            end
            mean_both_sides_yz = a_sum
            if divide
                mean_both_sides_yz /= cnt
            end

            a_sum = 0
            cnt = 0
            # scanning every xz: side 1
            for _x in 1:x_max, _z in 1:z_max
                frontier_y = y_max
                while frontier_y != 1 && (isnan(data[_x, frontier_y, _z, j, i]) || (greaterzero && data[_x, frontier_y, _z, j, i] == 0))
                    frontier_y -= 1
                end
                if data[_x, frontier_y, _z, j, i] >= 0 || (greaterzero && data[_x, frontier_y, _z, j, i] > 0)
                    a_sum += data[_x, frontier_y, _z, j, i]
                    cnt += 1
                end
            end
            # scanning every yz: side 2
            for _x in 1:x_max, _z in 1:z_max
                frontier_y = 1
                while frontier_y != y_max && (isnan(data[_x, frontier_y, _z, j, i]) || (greaterzero && data[_x, frontier_y, _z, j, i] == 0))
                    frontier_y += 1
                end
                if data[_x, frontier_y, _z, j, i] >= 0 || (greaterzero && data[_x, frontier_y, _z, j, i] > 0)
                    a_sum += data[_x, frontier_y, _z, j, i]
                    cnt += 1
                end
            end
            mean_both_sides_xz = a_sum
            if divide
                mean_both_sides_xz /= cnt
            end

            av_arr[j, i] = (mean_both_sides_xy + mean_both_sides_yz + mean_both_sides_xz)/3
        else
            av_arr[j, i] = mean_both_sides_xy
        end
    end
    return av_arr
end

function average_front(re, dataname; greaterzero=false, oneside=false, divide=true)
    average_front(re[dataname], re["stats"]["n_gens"], re["stats"]["max"]...; greaterzero=greaterzero, oneside=oneside, divide=divide)
end

"""
Finds the front array of `dataname` in `re`.

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

`oneside`: if **true**, approach only from one side (i.e. from the positive direction in strip expansions)

---

Output: front array of the same dimensions as `re[dataname]` (space + time)
"""
function front_array(re, dataname; oneside=false)
    front_array(re[dataname], re["stats"]["n_gens"], re["stats"]["max"]...; oneside=oneside)
end

# 1D
function front_array(data, n_gens, x_max; oneside=false)
    front_arr = fill(NaN, x_max, n_gens)

    for j in 1:n_gens
        frontier = x_max
        while frontier != 1 && isnan(data[frontier, j])
            frontier -= 1
        end
        if !isnan(data[frontier, j])
            front_arr[frontier, j] = data[frontier, j]
        end
        if !oneside
            frontier = 1
            while frontier != x_max && isnan(data[frontier, j])
                frontier += 1
            end
            if !isnan(data[frontier, j])
                front_arr[frontier, j] = data[frontier, j]
            end
        end
    end
    return front_arr
end

# 2D
function front_array(data::Array, n_gens, x_max, y_max; oneside=false)
    front_arr = fill(NaN, x_max, y_max, n_gens)
    for j in 1:n_gens
        # scanning every y: side 1
        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && isnan(data[frontier_x, _y, j])
                frontier_x -= 1
            end
            if !isnan(data[frontier_x, _y, j])
                front_arr[frontier_x, _y, j] = data[frontier_x, _y, j]
            end
        end
        # scanning every y: side 2
        if !oneside
            for _y in 1:y_max
                frontier_x = 1
                while frontier_x != x_max && isnan(data[frontier_x, _y, j])
                    frontier_x += 1
                end
                if !isnan(data[frontier_x, _y, j])
                    front_arr[frontier_x, _y, j] = data[frontier_x, _y, j]
                end
            end
        end

        if !oneside
            # scanning every x: side 1
            for _x in 1:x_max
                frontier_y = y_max
                while frontier_y != 1 && isnan(data[_x, frontier_y, j])
                    frontier_y -= 1
                end
                if !isnan(data[_x, frontier_y, j])
                    front_arr[_x, frontier_y, j] = data[_x, frontier_y, j]
                end
            end
            # scanning every x: side 2
            for _x in 1:x_max
                frontier_y = 1
                while frontier_y != y_max && isnan(data[_x, frontier_y, j])
                    frontier_y += 1
                end
                if !isnan(data[_x, frontier_y, j] > 0)
                    front_arr[_x, frontier_y, j] = data[_x, frontier_y, j]
                end
            end
        end
    end
    return front_arr
end

# 3D
function front_array(data::Array, n_gens, x_max, y_max, z_max; oneside=false)
    front_arr = fill(NaN, x_max, y_max, z_max, n_gens)

    for j in 1:n_gens
        # scanning every xy: side 1
        for _x in 1:x_max, _y in 1:y_max
            frontier_z = z_max
            while frontier_z != 1 && isnan(data[_x, _y, frontier_z, j])
                frontier_z -= 1
            end
            if !isnan(data[_x, _y, frontier_z, j])
                front_arr[_x, _y, frontier_z, j] = data[_x, _y, frontier_z, j]
            end
        end
        # scanning every xy: side 2
        if !oneside
            for _x in 1:x_max, _y in 1:y_max
                frontier_z = 1
                while frontier_z != z_max && isnan(data[_x, _y, frontier_z, j])
                    frontier_z += 1
                end
                if !isnan(data[_x, _y, frontier_z, j])
                    front_arr[_x, _y, frontier_z, j] = data[_x, _y, frontier_z, j]
                end
            end
        end

        if !oneside
            # scanning every yz: side 1
            for _y in 1:y_max, _z in 1:z_max
                frontier_x = x_max
                while frontier_x != 1 && isnan(data[frontier_x, _y, _z, j])
                    frontier_x -= 1
                end
                if !isnan(data[frontier_x, _y, _z, j])
                    front_arr[frontier_x, _y, _z, j] = data[frontier_x, _y, _z, j]
                end
            end
            # scanning every yz: side 2
            for _y in 1:y_max, _z in 1:z_max
                frontier_x = 1
                while frontier_x != x_max && isnan(data[frontier_x, _y, _z, j])
                    frontier_x += 1
                end
                if !isnan(data[frontier_x, _y, _z, j])
                    front_arr[frontier_x, _y, _z, j] = data[frontier_x, _y, _z, j]
                end
            end

            # scanning every xz: side 1
            for _x in 1:x_max, _z in 1:z_max
                frontier_y = y_max
                while frontier_y != 1 && isnan(data[_x, frontier_y, _z, j])
                    frontier_y -= 1
                end
                if !isnan(data[_x, frontier_y, _z, j])
                    front_arr[_x, frontier_y, _z, j] = data[_x, frontier_y, _z, j]
                end
            end
            # scanning every yz: side 2
            for _x in 1:x_max, _z in 1:z_max
                frontier_y = 1
                while frontier_y != y_max && isnan(data[_x, frontier_y, _z, j, i])
                    frontier_y += 1
                end
                if !isnan(data[_x, frontier_y, _z, j])
                    front_arr[_x, frontier_y, _z, j] = data[_x, frontier_y, _z, j]
                end
            end
        end
    end
    return front_arr
end

# mean front fitness (or other data)
"""
Produces a normalised copy of a time series `ts` using the "maximum normalisation" method: starting from the onset generation (`n_gens_burnin`**+1**), divide `ts` by the maximum of all `av_data` in each generation.

---

`ts`: time series to normalise

`av_data`: deme data array of dimensions **spatial+1** (averaged over replicates)

`n_gens_burnin`: number of burn-in generations

---

Output: normalised array of the same dimensions as `re[dataname]` (space + time)

---
"""
function norm_maximum(ts, av_data, n_gens_burnin::Int)
    normal_array = copy(ts)
    start = n_gens_burnin+1
    for i in start:length(ts)
        normal_array[i] /= maximum(li(av_data,i))
    end
    return normal_array
end

"""
Produces a normalised copy of a time series `ts` using the "onset mean normalisation" method: starting from the onset generation (`n_gens_burnin`**+1**), divide `ts` by its value at `n_gens_burnin`**+1**.

---

`ts`: time series to normalise

`n_gens_burnin`: number of burn-in generations

---

Output: normalised array of the same dimensions as `re[dataname]` (space + time)
"""
function norm_onset_mean(ts, n_gens_burnin::Int)
    normal_array = copy(ts)
    start = n_gens_burnin+1
    normal_array[start:end] /= ts[start]
    return normal_array
end

"""
Finds an average (at each generation) over multiple similar time series (replicates).

---

`ts_arr`: time series array

`n_gens_burnin`: number of burn-in generations

---

Output: averaged time series
"""
function average_ts(ts_arr, n_gens)
    return [mean(ts_arr[i,:]) for i in 1:n_gens]
end


# Upcoming features
# ------------------------------------------------

#= """
Gives description for a method.
"""
macro d(x)
    quote
        display("text/markdown", @doc $x)
    end    
end =#

vc(x) = cat(eachslice(x, dims=4)..., dims=2)

function re_get_avrel(data::Array, x, gen, denom)
    nd = ndims(data)
    if nd == 4
        return mean(vc(data)[x, :, gen]) / denom
    elseif nd == 3
        return mean(data[x, :, gen]) / denom
    else
        println("Wrong data type.")
    end
end
function re_get_avrel(re::Dict, dataname::String, x, gen=Int(re["stats"]["n_gens"]); sel=true)
    denom = sel ? re["stats"]["n_sel_loci"] : re["stats"]["n_loci"] - re["stats"]["n_sel_loci"]
    return re_get_avrel(re[dataname], x, gen, denom)
end

function re_plot_avrelselneu(re::Dict, dataname::String, x_range=(1:Int(re["stats"]["x_max"])); x_scale_factor=1, sel=true, overlay=false)
    nd = ndims(re[dataname*"sel"])
    if nd == 4
        data1 = vc(re[dataname*"sel"])
        data2 = vc(re[dataname*"neu"])
    else
        data1 = re[dataname*"sel"]
        data2 = re[dataname*"neu"]
    end
    t = [re_get_avrel(data1, j, Int(re["stats"]["n_gens"]), re["stats"]["n_sel_loci"]) for j in x_range]
    t2 = [re_get_avrel(data2, j, Int(re["stats"]["n_gens"]), re["stats"]["n_loci"] - re["stats"]["n_sel_loci"]) for j in x_range]

    if haskey(re["stats"], "name")
        lbl1 = re["stats"]["name"] * "[selected $dataname]"
        lbl2 = re["stats"]["name"] * "[neutral $dataname]"
    else
        lbl1 = "selected $dataname"
        lbl2 = "neutral $dataname"
    end

    if overlay
        plot!(x_range * x_scale_factor, t, label=lbl1, xlabel="x")
    else
        plot(x_range * x_scale_factor, t, label=lbl1, xlabel="x")
    end
    #plot!(x_range*x_scale_factor,t2,label=lbl2)
end

function re_plot_avrelselneu!(re::Dict, dataname::String, x_range=(1:Int(re["stats"]["x_max"])); x_scale_factor=1, sel=true, overlay=false)
    re_plot_avrelselneu(re, dataname, x_range; x_scale_factor=x_scale_factor, sel=sel, overlay=true)
end
# ------------------------------------------------

println("RESK successfully loaded.")
using StatsBase, Distributions, Distributed, Random, Serialization, Dates, DataStructures, ThreadsX #SpecialFunctions
include("defaults.jl")

const RESK_VERSION = v"0.3.7"

# Constants
# ------------------------------------------------

# Migration probabilities for each dimensionality -> for each mode
MIGR_DIRS_ORT_1D = [[1], [-1]]
MIGR_DIRS_ORT_2D = [[-1, 0], [0, -1], [0, 1], [1, 0]]
MIGR_DIRS_ORT_3D = [[-1, 0, 0], [1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, -1], [0, 0, 1]]

MIGR_DIRS_DIAG_2D = [[-1, -1], [-1, 1], [1, -1], [1, 1]]
MIGR_DIRS_HEX_2D = [[-1, 0], [0, -1], [-1, 1], [0, 1], [1, 0], [1, 1]]

# Predefined probabilities instead of dictionary lookups
MIGR_PROB_ORT_1D = (1.0, 0.0)
MIGR_PROB_ORT_2D = (1.0, 0.0)
MIGR_PROB_ALL_2D = (0.5, 0.5)
MIGR_PROB_DIAG1_2D = (2.0/3.0, 1.0/3.0)

# Types
const typ_float = Float32
const typ_int = Int32
const typ_gt = Array{Array{Array{Bool}}} # Doesn't need to be constant though
const typ_gt_inf = Array{Array{Array{typ_float}}}
# Also, what about Array{Array{Array{Bool}},wlddim} !!!


# Common functions and structures
# ------------------------------------------------

function resk_version()
    return RESK_VERSION
end

"""
Returns the maximum in an array ignoring NaNs.
"""
function maxdef(a)
    return maximum(filter(!isnan,a))
end

"""
Returns the minimum in an array ignoring NaNs.
"""
function mindef(a)
    return minimum(filter(!isnan,a))
end

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
function ins_sq(r_max_burnin, x_max)
    return (trunc(Int, 1 + x_max/2 - r_max_burnin * 0.666)):(trunc(Int, 1 + x_max/2 + r_max_burnin * 0.666))
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
function ins_cb(r_max_burnin, x_max)
    return (trunc(Int, 1 + x_max/2 - r_max_burnin * 0.577)):(trunc(Int, 1 + x_max/2 + r_max_burnin * 0.577))
end

function sphere_lattice_iter(x::Real, y::Real, z::Real, r::Real)
    r ≥ 0 || throw(ArgumentError("r must be nonnegative"))
    xmin = max(1, ceil(Int, x - r)); xmax = floor(Int, x + r)
    ymin = max(1, ceil(Int, y - r)); ymax = floor(Int, y + r)
    zmin = max(1, ceil(Int, z - r)); zmax = floor(Int, z + r)
    r2 = r^2
    Iterators.filter(((i,j,k),)->(i-x)^2 + (j-y)^2 + (k-z)^2 <= r2,
                     Iterators.product(xmin:xmax, ymin:ymax, zmin:zmax))
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
function calc_offspring(wld, stats)
    next_gen_posits = Vector{Vector{Int}}()
    next_gen_pops = Array{Int16}(undef, stats["max"]...)
    fill!(next_gen_pops, -1)
    
    @inbounds for idx in CartesianIndices(stats["max"])
        if length(wld[idx]) > 0
            n_ppl_at_deme = length(wld[idx])
            expected_offspring = n_ppl_at_deme * (stats["prolif_rate"] / 
                               (1 + (n_ppl_at_deme * (stats["prolif_rate"] - 1)) / stats["capacity"]))

            next_gen_pops[idx] = rand(Poisson(expected_offspring))
            if next_gen_pops[idx] > 0
                push!(next_gen_posits, [Tuple(idx)...])
            end
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
function calc_migr_dist!(move::Vector{Int16}, deme::Vector{Int}, stats, 
                        max_migr::Tuple, r_max_migr::Int64=0, r_coords::Vector{Int}=[1, 2], bottleneck=NaN, verbose=false)
    wlddim = stats["wlddim"]
    fill!(move, 0)
    
    # Early return if no migration
    if rand() >= stats["migr_rate"]
        return
    end
    
    # Migration mode handling - simplified with fewer conditionals
    if wlddim == 1
        # 1D has only orthogonal directions
        dir = rand() < 0.5 ? MIGR_DIRS_ORT_1D[1] : MIGR_DIRS_ORT_1D[2]
    elseif wlddim == 2
        # 2D case
        if stats["migr_mode"] == "hex"
            dir_idx = rand(1:length(MIGR_DIRS_HEX_2D))
            dir = MIGR_DIRS_HEX_2D[dir_idx]
        elseif stats["migr_mode"] == "all"
            if rand() < 0.5
                dir_idx = rand(1:length(MIGR_DIRS_ORT_2D))
                dir = MIGR_DIRS_ORT_2D[dir_idx]
            else
                dir_idx = rand(1:length(MIGR_DIRS_DIAG_2D))
                dir = MIGR_DIRS_DIAG_2D[dir_idx]
            end
        else
            # Default to orthogonal
            dir_idx = rand(1:length(MIGR_DIRS_ORT_2D))
            dir = MIGR_DIRS_ORT_2D[dir_idx]
        end
    else
        # 3D case - just orthogonal for simplicity
        dir_idx = rand(1:length(MIGR_DIRS_ORT_3D))
        dir = MIGR_DIRS_ORT_3D[dir_idx]
    end
    
    # Copy direction to output
    for i in 1:length(dir)
        move[i] = dir[i]
    end
    if verbose
        println("deme ",deme)
        println(move)
    end
    
    # Apply radius constraint if needed
    if r_max_migr > 0
        r_arr = [(deme[i] - (stats["max"][i]-1)/2 + move[i] - 1)^2 for i in r_coords]
        r2 = sum(r_arr)
        
        if r2 > r_max_migr * r_max_migr
            #factor = r_max_migr*r_max_migr/r2
            move[r_coords] .= 0 # Do [trunc(Int16, factor * move[i]) for i in 1:wlddim] in the future (with multiple-deme jumps)
        end
    end
    if verbose
        println(move)
    end
    
    # Inside certain square check:
    for i in 1:wlddim
        try_move = deme[i] + move[i]
        if isa(max_migr, Tuple) && !isnan(max_migr[i])
            if try_move > max_migr[i]
                move[i] = 0
            end
        end
        if try_move < 1
            move[i] = 0
        end
    end
    if verbose
        println(move)
        println("--")
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
    #println(move,deme)
end


# Simulation functions (finite-sites)
# ------------------------------------------------

function calc_fitness(deme_ms1, deme_ms2, domin_coef, loci, sel_loci=[])
    fits = Vector{typ_float}(undef, length(deme_ms1))
    n_loci = length(loci)

    for i in 1:length(deme_ms1)
        new_fitness = 1.0

        for j in 1:n_loci
            if deme_ms1[i][j] == true && deme_ms2[i][j] == true
                if j in sel_loci
                    new_fitness *= 1 + loci[j]
                end

            elseif deme_ms1[i][j] == true || deme_ms2[i][j] == true
                if j in sel_loci
                    new_fitness *= 1 + domin_coef * loci[j]
                end
            end
        end
        fits[i]=new_fitness
    end
    return fits
end

function calc_fitness(deme_ms1, deme_ms2, wld_stats)
    calc_fitness(deme_ms1, deme_ms2, wld_stats["domin_coef"], wld_stats["loci"], wld_stats["sel_loci"])
end

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
function calc_fitness_and_states(deme_ms1, deme_ms2, domin_coef, loci, sel_loci=[])
    fits = Vector{typ_float}(undef, length(deme_ms1))
    n_loci = length(loci)

    cAA = zeros(n_loci)
    cAa = zeros(n_loci)
    caa = zeros(n_loci)

    for i in 1:length(deme_ms1)
        new_fitness = 1.0

        for j in 1:n_loci
            if deme_ms1[i][j] == true && deme_ms2[i][j] == true
                if j in sel_loci
                    new_fitness *= 1 + loci[j]
                end
                cAA[j] += 1

            elseif deme_ms1[i][j] == true || deme_ms2[i][j] == true
                if j in sel_loci
                    new_fitness *= 1 + domin_coef * loci[j]
                end
                cAa[j] += 1
            else
                caa[j] += 1
            end
        end
        fits[i]=new_fitness
    end

    return fits, cAA, cAa, caa
end

function calc_fitness_and_states(deme_ms1, deme_ms2, wld_stats)
    calc_fitness_and_states(deme_ms1, deme_ms2, wld_stats["domin_coef"], wld_stats["loci"], wld_stats["sel_loci"])
end

"""
Randomly mutates at selected loci. Used when building the next generation in finite-sites expansions.

---

`ms1`: an individual's left monosome array

`ms2`: an individual's right monosome array

`mut_rate`: genome-wide mutation rate

`n_loci`: number of loci
"""
function mutate!(ms1, ms2, mut_rate, n_loci; mutratelocus=false)
    muts = 0
    get_mutation_random = mutratelocus ? rand(Poisson(mut_rate*n_loci)) : rand(Poisson(mut_rate))
    @fastmath @inbounds for _ in 1:get_mutation_random
        pos_alter = sample(1:n_loci)

        if rand(1:2) == 1
            ms1[pos_alter] = true
        else
            ms2[pos_alter] = true
        end
        muts += 1
    end
    return muts
end

function mutate!(ms1, ms2, wld_stats; mutratelocus=false)
    return mutate!(ms1, ms2, wld_stats["mut_rate"], wld_stats["n_loci"]; mutratelocus=mutratelocus)
end

"""
Recombines loci with a 1/2 chance. Used when building the next generation in finite-sites expansions.

---

`ms1`: an individual's left monosome array

`ms2`: an individual's right monosome array

`mut_rate`: genome-wide mutation rate

`n_loci`: number of loci
"""
function crossover!(ms1, ms2, n_loci)
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

`wld_gt1`: a spatial array of demes that contain individuals' left monosome [Bool] arrays

`wld_gt2`: a spatial array of demes that contain individuals' right monosome [Bool] arrays

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

`premutate`: 

---

Output 1: a changed `wld_gt1` = a spatial array of demes that contain individuals' left monosome [Bool] arrays

Output 2: a changed `wld_gt2` = a spatial array of demes that contain individuals' right monosome [Bool] arrays

Output 3: a spatial array of demes with average fitness in the new generation

Output 4: a spatial array of demes with populations in the new generation

Output 5: a spatial array of demes with average selected AA mutation count in the new generation

Output 6: a spatial array of demes with average selected Aa mutation count in the new generation

Output 7: a spatial array of demes with average selected aa mutation count in the new generation

Output 8: a spatial array of demes with average neutral AA mutation count in the new generation

Output 9: a spatial array of demes with average neutral Aa mutation count in the new generation

Output 10: a spatial array of demes with average neutral aa mutation count in the new generation
"""
function build_next_gen!(g, gg, out_fields, wld_gt1::typ_gt, wld_gt2::typ_gt, stats, wld_fitn=NaN, wld_pops=NaN, wld_cAA=NaN, wld_cAa=NaN, wld_caa=NaN, wld_AA=NaN, wld_Aa=NaN, wld_aa=NaN;
    max_migr::Tuple=stats["max"], migr_mode=DEF_MIGR_MODE, r_max_migr::Int64=0, bottleneck=NaN, refl_walls=false, r_coords::Array{Int}=[1, 2], relcnt::Bool=false,
    weightfitn::Bool=true, condsel::Bool=false, mutratelocus::Bool=false,
    fixed_mate::Bool=false, SS::Bool=false, verbose=false)

    # Preallocate buffers for migration
    move_buffer = zeros(Int16, stats["wlddim"])
    wlddim = stats["wlddim"]
    #Random.seed!(1234)

    # Calculate offspring counts
    next_gen_posits, next_gen_pops = calc_offspring(wld_gt1[repeat([:],wlddim)...,gg], stats)
    if verbose
        println("g=",g, ", gg=",gg,", pops=", next_gen_pops)
    end

    # Process each deme
    for deme in next_gen_posits
        ms1_at_pos = wld_gt1[deme...,gg]
        ms2_at_pos = wld_gt2[deme...,gg]

        # Calculate fitnesses once
        fitnesses = calc_fitness(ms1_at_pos, ms2_at_pos, stats)
        
        next_generation_size = next_gen_pops[deme...]

        if next_generation_size > 0
            # Create births
            k = 0
            @inbounds while k < next_generation_size
                # Select parents
                mom_idx, dad_idx = 1, 1
                if weightfitn
                    mom_idx = wsample(1:length(ms1_at_pos), fitnesses)
                    dad_idx = wsample(1:length(ms1_at_pos), fitnesses)
                else
                    mom_idx = rand(1:length(ms1_at_pos))
                    dad_idx = rand(1:length(ms1_at_pos))
                end
                mom_ms1 = ms1_at_pos[mom_idx]
                mom_ms2 = ms2_at_pos[mom_idx]
                dad_ms1 = ms1_at_pos[dad_idx]
                dad_ms2 = ms2_at_pos[dad_idx]
                
                thing1 = rand()
                thing2 = rand()
                if verbose
                    print("Mate rand()=",thing1," rand()=", thing2," ")
                end
                if !condsel || (fitnesses[mom_idx] > thing1*maximum(fitnesses) && 
                    fitnesses[dad_idx] > thing2*maximum(fitnesses))
                    k+=1
                    # Create gametes with recombination
                    gamete_mom_ms1 = copy(mom_ms1)
                    gamete_dad_ms1 = copy(dad_ms1)
                    gamete_mom_ms2 = copy(mom_ms2)
                    gamete_dad_ms2 = copy(dad_ms2)
                    
                    # Apply recombination
                    crossover!(gamete_mom_ms1, gamete_mom_ms2, stats["n_loci"])
                    crossover!(gamete_dad_ms1, gamete_dad_ms2, stats["n_loci"])

                    mut_cnt = mutate!(gamete_mom_ms1, gamete_dad_ms2, stats; mutratelocus=mutratelocus)

#=                     if out_fields["mutcnt"][1]
                        if !out_fields["mutcnt"][2]
                            G = g
                        elseif g == stats["n_gens"]
                            G = 1
                        end
                        if !isnan(wld_mutcnt)
                            if isnan(wld_mutcnt[deme...,G])
                                wld_mutsben[deme...,G] = 0
                            end
                            wld_mutcnt[deme...,G] += mut_cnt
                        end
                    end =#

                    # Calculate migration
                    calc_migr_dist!(move_buffer, deme, stats, max_migr, r_max_migr, r_coords)

                    # Apply migration to calculate new position
                    indices = Vector{Int}(undef, wlddim)

                    for i in 1:wlddim
                        indices[i] = deme[i] + move_buffer[i]
                    end
                    
                    # Add offspring to next generation
                    push!(wld_gt1[indices...,gg+1], copy(gamete_mom_ms1))
                    push!(wld_gt2[indices...,gg+1], copy(gamete_dad_ms2))

                    if verbose
                        println("Completed birth #",k," rand()=",rand(1))
                    end
                elseif !SS
                    k+=1
                end
            end
            
        end
        if out_fields["pops"][1]
            G = out_fields["pops"][2] ? g+1 : 1
            wld_pops[deme...,G]=length(wld_gt1[deme...,gg+1])
        end
        if out_fields["fitn"][1] || out_fields["cAA"][1] || out_fields["AA"][1]
            fitn,cAA,cAa,caa = calc_fitness_and_states(wld_gt1[deme...,gg+1],wld_gt2[deme...,gg+1],stats)
            if out_fields["fitn"][1]
                G = out_fields["fitn"][2] ? g+1 : 1
                wld_fitn[deme...,G]=mean(fitn)
            end
            if out_fields["cAA"][1]
                G = out_fields["cAA"][2] ? g+1 : 1
                for ww in 1:stats["n_loci"]
                    wld_cAA[deme...,ww,G]=cAA[ww]
                    wld_cAa[deme...,ww,G]=cAa[ww]
                    wld_caa[deme...,ww,G]=caa[ww]
                end
            end
            if out_fields["AA"][1]
                G = out_fields["AA"][2] ? g+1 : 1
                wld_AA[deme...,G] = relcnt ? sum(cAA)/stats["n_loci"]/lenn : sum(cAA)/stats["n_loci"]
                wld_Aa[deme...,G] = relcnt ? sum(cAa)/stats["n_loci"]/lenn : sum(cAa)/stats["n_loci"]
                wld_aa[deme...,G] = relcnt ? sum(caa)/stats["n_loci"]/lenn : sum(caa)/stats["n_loci"]
            end
        end
    end

    if verbose
        println("--------------------")
    end
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
function create_empty_world(maxi=(DEF_X_MAX, DEF_Y_MAX); min=(1, 1), name="sim_$(Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"))", capacity=DEF_CAPACITY,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, migr_mode=DEF_MIGR_MODE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS,
    prolif_rate=DEF_PROLIF_RATE, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci))

    wld_gt1 = (typ_gt)(undef, maxi...) # array of left (in a pair) monosomes ("ms") of all individuals in space
    wld_gt2 = (typ_gt)(undef, maxi...) # array of right (in a pair) monosomes ("ms") of all individuals in space

    for k in CartesianIndices(maxi)
        wld_gt1[k] = Array{Bool,1}[]
        wld_gt2[k] = Array{Bool,1}[]
    end

    wld_stats = OrderedDict(
        "name" => name,
        #"min" => min,
        "max" => maxi,
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
        "wlddim" => length(maxi)
        #"rangeexps" => []
    )

    return wld_gt1, wld_gt2, wld_stats
end

"""
Fills random demes within given monosome arrays with finite-sites individuals. Usually used after an empty world is created.

---

`wld_gt1`: a spatial array of demes that contain individuals' left monosome [Bool] arrays

`wld_gt2`: a spatial array of demes that contain individuals' right monosome [Bool] arrays

`wld_stats`: world stats Dict

`fill`: an array of Int ranges of the coordinates that define the area within which to fill

`n_demes_to_fill`: number of demes to fill

"""
function fill_random_demes(wld_gt1, wld_gt2, wld_stats, fill_iter, n_demes_to_fill=DEF_N_DEMES_STARTFILL; fill_range=NaN, redims=NaN, verbose=false)
    if !isnothing(fill_iter)
        possible_init_coords = collect(fill_iter)
        #wld_stats["startfill"] = collect(fill_iter)
    else
        possible_init_coords = [collect(x) for x in Iterators.product(fill_range...)]
        wld_stats["startfill"] = copy(fill_range)
    end

    if length(possible_init_coords)<n_demes_to_fill
        println("n_demes_to_fill is more than existing demes!")
    end

    init_coords = sample(possible_init_coords, n_demes_to_fill; replace=false)
    wld_stats["sel_loci"] = randperm(wld_stats["n_loci"])[1:wld_stats["n_sel_loci"]]

    if haskey(wld_stats,"prop_of_del_muts")
        n_ids = Int(wld_stats["prop_of_del_muts"]*wld_stats["n_loci"])
        ids = randperm(wld_stats["n_loci"])[1:n_ids]
        wld_stats["loci"][ids] .*= -1
    end

    if verbose
        println("Filling ",first(wld_gt1,5)," & ",first(wld_gt2,5)," .")
    end

    for coord in init_coords
        for _ in 1:wld_stats["capacity"]
            push!(wld_gt1[coord...,redims...], falses(wld_stats["n_loci"]))
            push!(wld_gt2[coord...,redims...], falses(wld_stats["n_loci"]))
        end
    end
    if verbose
        println("Filled ",first(wld_gt1,5)," & ",first(wld_gt2,5)," .")
    end
    
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
- **C** - homo- and heterozygous counts in a deme per locus (**cAA**, **cAa** and **caa**)
- **A** - homo- and heterozygous counts averaged per deme (**AA**, **Aa** and **aa**). When used with relcount, these are observed homo- and heterozygosities

`name`: world name

`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates

`r_max_burnin`: radius that bounds the burn-in area

`r_max_exp`: radius that bounds the expansion area

`r_coords`: a tuple (array) of axes' ordinal numbers that the n-sphere with `r_max_migr` covers. For example:
- **(1,3)** - migration is bound within a disk at x and z axes
- **(1,2,3)** - migration is bound within a sphere at x, y and z axes

`loci`: array of selected coefficients for every locus

`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start

`multiproc`: if **true**, distribute to threads

If starting from existing world, also provide:

`wld_gt1`: a spatial array of demes that contain individuals, each of which is a Bool array representing left monosomes

`wld_gt2`: a spatial array of demes that contain individuals, each of which is a Bool array representing right monosomes

`wld_stats`: world stats Dict

---

Output: a Dict containing data after the expansion:
- **stats** - statistics array containing world and range expansion information
- **fitn**, **pops**, **AAsel**, **Aasel**, **aasel**, **AAneu**, **Aaneu**, **aaneu** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`
"""
function rangeexp(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; max_burnin=(DEF_X_MAX_BURNIN, DEF_Y_MAX), max_exp=(DEF_X_MAX_EXP, DEF_Y_MAX), maxi=(DEF_X_MAX, DEF_Y_MAX), migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, r_max_burnin=0, r_max_exp=0, r_coords=[1, 2], capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE,
    prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, verbose=false, condsel=false, fixed_mate=false, startfill_iter=nothing, relcnt=false,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci), mutratelocus=false, weightfitn=true, 
    startfill_range=NaN, multiproc=true, wld_gt1=NaN, wld_gt2=NaN, wld_stats=NaN, SS=true)
    
    if isnan(wld_gt1)
        is_fill_random_demes = true
        #println("No world provided. Creating a new world.")
        wld_gt1, wld_gt2, wld_stats = create_empty_world(maxi; name=name, capacity=capacity, prolif_rate=prolif_rate, 
        mut_rate=mut_rate, migr_rate=migr_rate, migr_mode=migr_mode, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, domin_coef=domin_coef,
            n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci)
            if isnothing(startfill_iter) && !isa(startfill_range, Array) && !any(isnan, max_burnin)
                startfill_range = [1:upper for upper in max_burnin]
            end
    end
    println("max_burnin ",max_burnin)
    println("startfill_range ",startfill_range)
    println("startfill_iter ",startfill_iter)
    wlddim = wld_stats["wlddim"]
    wld_stats["max_burnin"] = max_burnin
    wld_stats["max_exp"] = max_exp
    wld_stats["n_gens_burnin"] = n_gens_burnin
    wld_stats["n_gens_exp"] = n_gens_exp
    n_gens = n_gens_burnin + n_gens_exp
    wld_stats["n_gens"] = n_gens

    function check_what_output(letter)
        return (occursin(letter, data_to_generate),occursin(letter*"l", data_to_generate))
    end
    out_fields = OrderedDict{String, Any}("gt1" => check_what_output("G"), "gt2" => check_what_output("G"), "fitn" => check_what_output("F"), "pops" => check_what_output("P"),
        "cAA" => check_what_output("C"),"cAa" => check_what_output("C"),"caa" => check_what_output("C"), "AA" => check_what_output("A"), "Aa" => check_what_output("A"), "aa" => check_what_output("A"))

    function tsk_re(_proc_no=NaN)
        cond_n_gens = out_fields["gt1"][1] && out_fields["gt1"][2] ? n_gens+1 : 2
        rewld_gt1_local = (typ_gt)(undef, wld_stats["max"]..., cond_n_gens)
        rewld_gt2_local = (typ_gt)(undef, wld_stats["max"]..., cond_n_gens)

        rewld_gt1_local[repeat([:],wlddim)...,1] = deepcopy(wld_gt1)
        rewld_gt2_local[repeat([:],wlddim)...,1] = deepcopy(wld_gt2)
        for og in 2:cond_n_gens
            rewld_gt1_local[repeat([:],wlddim)...,og] = deepcopy(wld_gt1)
            rewld_gt2_local[repeat([:],wlddim)...,og] = deepcopy(wld_gt2)
        end
        if is_fill_random_demes
            fill_random_demes(rewld_gt1_local, rewld_gt2_local, wld_stats, startfill_iter; redims=(1), fill_range=startfill_range, verbose=verbose)
        end
        
        if out_fields["fitn"][1]
            cond_n_gens = out_fields["fitn"][2] ? n_gens+1 : 1
            rewld_fitn_local = Array{typ_float}(undef, wld_stats["max"]..., cond_n_gens)
            fill!(rewld_fitn_local,NaN)
            if out_fields["fitn"][2]
                for deme in CartesianIndices(maxi)
                    rewld_fitn_local[deme,1] = mean(calc_fitness(rewld_gt1_local[deme,1],rewld_gt2_local[deme,1],wld_stats))
                end
            end
        else
            rewld_fitn_local = NaN
        end
        if out_fields["pops"][1]
            cond_n_gens = out_fields["pops"][2] ? n_gens+1 : 1
            rewld_pops_local = Array{typ_float}(undef, wld_stats["max"]..., cond_n_gens)
            fill!(rewld_pops_local,NaN)
            if out_fields["pops"][2]
                for deme in CartesianIndices(maxi)
                    rewld_pops_local[deme,1] = length(rewld_gt1_local[deme,1])
                end
            end
        else
            rewld_pops_local = NaN
        end
        if out_fields["cAA"][1]
            cond_n_gens = out_fields["cAA"][2] ? n_gens+1 : 1
            rewld_cAA_local = Array{typ_float}(undef, wld_stats["max"]..., wld_stats["n_loci"], cond_n_gens)
            rewld_cAa_local = Array{typ_float}(undef, wld_stats["max"]..., wld_stats["n_loci"], cond_n_gens)
            rewld_caa_local = Array{typ_float}(undef, wld_stats["max"]..., wld_stats["n_loci"], cond_n_gens)
            fill!(rewld_cAA_local,NaN)
            fill!(rewld_cAa_local,NaN)
            fill!(rewld_caa_local,NaN)
        else
            rewld_cAA_local = NaN
            rewld_cAa_local = NaN
            rewld_caa_local = NaN
        end
        if out_fields["AA"][1]
            cond_n_gens = out_fields["AA"][2] ? n_gens+1 : 1
            rewld_AA_local = Array{typ_float}(undef, wld_stats["max"]..., cond_n_gens)
            rewld_Aa_local = Array{typ_float}(undef, wld_stats["max"]..., cond_n_gens)
            rewld_aa_local = Array{typ_float}(undef, wld_stats["max"]..., cond_n_gens)
            fill!(rewld_AA_local,NaN)
            fill!(rewld_Aa_local,NaN)
            fill!(rewld_aa_local,NaN)
        else
            rewld_AA_local = NaN
            rewld_Aa_local = NaN
            rewld_aa_local = NaN
        end

        @inbounds for g in 1:n_gens
            if g <= n_gens_burnin
                max_migr = max_burnin
                r_max_migr = r_max_burnin
            else
                max_migr = max_exp
                r_max_migr = r_max_exp
            end

            condi = out_fields["gt1"][1] && out_fields["gt1"][2]
            gg = condi ? g : 1

            if !condi
                for deme in CartesianIndices(maxi)
                    empty!(rewld_gt1_local[deme,2])
                    empty!(rewld_gt2_local[deme,2])
                end
            end

            build_next_gen!(g, gg, out_fields, rewld_gt1_local, rewld_gt2_local, wld_stats, rewld_fitn_local, rewld_pops_local, rewld_cAA_local, rewld_cAa_local, rewld_caa_local, rewld_AA_local, rewld_Aa_local, rewld_aa_local;
            max_migr=max_migr, migr_mode=migr_mode, r_max_migr=r_max_migr, bottleneck=bottleneck, r_coords=r_coords, relcnt=relcnt,
            mutratelocus=mutratelocus, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, 
            SS=SS, verbose=verbose)

            if !condi
                rewld_gt1_local[repeat([:],wlddim)...,1],rewld_gt1_local[repeat([:],wlddim)...,2] = rewld_gt1_local[repeat([:],wlddim)...,2],rewld_gt1_local[repeat([:],wlddim)...,1]
                rewld_gt2_local[repeat([:],wlddim)...,1],rewld_gt2_local[repeat([:],wlddim)...,2] = rewld_gt2_local[repeat([:],wlddim)...,2],rewld_gt2_local[repeat([:],wlddim)...,1]

            
            end
        end

        return Dict{String, Any}("gt1" => rewld_gt1_local, "gt2" => rewld_gt2_local, "fitn" => rewld_fitn_local, "pops" => rewld_pops_local, "AA" => rewld_AA_local, "Aa" => rewld_Aa_local, "aa" => rewld_aa_local,
        "cAA" => rewld_cAA_local, "cAa" => rewld_cAa_local, "caa" => rewld_caa_local)
    end

    extra_dim = 0

    if n_re > 1
        extra_dim += 1
        if multiproc
            npcs = workers()
            println("Running $n_re replicates on processes $npcs.")
            dicts_out,t = @timed pmap(tsk_re, 1:n_re)
            println("Simulated in ",t," s.")
            wld_stats["simtime"]=t
        else
            println("Using all threads.")
            dicts_out = ThreadsX.map(p->tsk_re(), 1:n_re)
        end
    else
        dicts_out,t = @timed [tsk_re()]
        println("Simulated in ",t," s.")
        wld_stats["simtime"]=t
    end

    res = OrderedDict{String, Any}("stats" => wld_stats)

    for key in keys(out_fields)
        if out_fields[key][1]
            if out_fields[key][2]
                if key=="cAA" || key=="cAa" || key=="caa" || key=="gt"
                    res[key] = cat([dictus[key] for dictus in dicts_out]...,dims=wlddim+extra_dim+1)
                else
                    res[key] = cat([dictus[key] for dictus in dicts_out]...,dims=wlddim+extra_dim+1)
                end
            else
                if key=="gt"
                    res[key] = cat([dictus[key][repeat([:],wlddim)...,1] for dictus in dicts_out]...,dims=wlddim+extra_dim)
                elseif key=="mutsdel" || key=="mutsben"
                    res[key] = cat([dictus[key][repeat([:],wlddim)...,1] for dictus in dicts_out]...,dims=wlddim+extra_dim)
                else
                    res[key] = cat([dictus[key][repeat([:],wlddim)...,1] for dictus in dicts_out]...,dims=wlddim+extra_dim)
                end
            end
        end
    end

    return res
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

`multiproc`: if **true**, distribute to threads

`capacity`: capacity of each deme

If starting from existing world, also provide:

`wld_gt1`: a spatial array of demes that contain individuals, each of which is a Bool array representing left monosomes

`wld_gt2`: a spatial array of demes that contain individuals, each of which is a Bool array representing right monosomes

`wld_stats`: world stats Dict

---

Output: a Dict containing data after the expansion:
- **stats** - statistics array containing world and range expansion information
- **fitn**, **pops**, **AAsel**, **Aasel**, **aasel**, **AAneu**, **Aaneu**, **aaneu** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`
"""
function rangeexp_ray(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; x_max_burnin=DEF_X_MAX_BURNIN, x_max_exp=DEF_X_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, weightfitn=true, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, multiproc=true, wld_gt1=NaN, wld_gt2=NaN, wld_stats=NaN, capacity=DEF_CAPACITY,
    prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, condsel=false, fixed_mate=false, SS=true, verbose=false)

    rangeexp(n_gens_burnin, n_gens_exp, n_re; max_burnin=(x_max_burnin,), max_exp=(x_max_exp,), maxi=(x_max_exp,), startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci, weightfitn=weightfitn, 
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt1=wld_gt1, wld_gt2=wld_gt2, wld_stats=wld_stats, name=name, bottleneck=bottleneck, multiproc=multiproc, condsel=condsel, fixed_mate=fixed_mate,
        SS=SS, verbose=verbose, prop_of_del_muts=prop_of_del_muts)
end

const rangeexp_1d = rangeexp_ray

function rangeexp_linear(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, maxi=(r_max_exp * 2 + 1,), 
    migr_mode=DEF_MIGR_MODE, startfill_range=[(1-ceil(Int,r_max_burnin/2)+r_max_exp):(1+ceil(Int,r_max_burnin/2)+r_max_exp)], prolif_rate=DEF_PROLIF_RATE, max_exp=(NaN,NaN), max_burnin=(NaN,NaN),
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, multiproc=true, wld_gt1=NaN, wld_gt2=NaN, wld_stats=NaN, capacity=DEF_CAPACITY)

    rangeexp(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, maxi=maxi, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci, r_coords=[1],
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt1=wld_gt1, wld_gt2=wld_gt2, wld_stats=wld_stats, name=name, bottleneck=bottleneck, multiproc=multiproc, condsel=condsel, fixed_mate=fixed_mate,
        SS=SS, verbose=verbose, prop_of_del_muts=prop_of_del_muts)
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

`wld_gt1`: a spatial array of demes that contain individuals, each of which is a Bool array representing left monosomes

`wld_gt2`: a spatial array of demes that contain individuals, each of which is a Bool array representing right monosomes

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
    prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, SS=true, verbose=false, relcnt=false,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=("midhole at x=", x_max_burnin * 2), prolif_rate=DEF_PROLIF_RATE,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, capacity=DEF_CAPACITY, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    wld_gt1=NaN, wld_gt2=NaN, wld_stats=NaN, max_burnin=(x_max_burnin, y_max), max_exp=(x_max_exp, y_max), maxi=(x_max_exp, y_max), multiproc=true)
    
    rangeexp(n_gens_burnin, n_gens_exp, n_re; max_burnin=max_burnin, max_exp=max_exp, maxi=maxi, startfill_range=startfill_range, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, capacity=capacity, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt1=wld_gt1, wld_gt2=wld_gt2, wld_stats=wld_stats, name=name, bottleneck=bottleneck, multiproc=multiproc, SS=SS, verbose=verbose,
        weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, prop_of_del_muts=prop_of_del_muts, relcnt=relcnt)
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

`wld_gt1`: a spatial array of demes that contain individuals, each of which is a Bool array representing left monosomes

`wld_gt2`: a spatial array of demes that contain individuals, each of which is a Bool array representing right monosomes

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
function rangeexp_disk(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, SS=true, verbose=false,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, maxi=(r_max_exp * 2 + 1, r_max_exp * 2 + 1), capacity=DEF_CAPACITY,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    max_exp=(NaN,NaN), max_burnin=(NaN,NaN), wld_gt1=NaN, wld_gt2=NaN, wld_stats=NaN)

    if !isa(startfill_range, Array)
        ran = ins_sq(r_max_burnin, maxi[1])
        startfill_range = [ran, ran]
    end

    rangeexp(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, maxi=maxi, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, SS=SS, verbose=verbose, 
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt1=wld_gt1, wld_gt2=wld_gt2, wld_stats=wld_stats, name=name, bottleneck=bottleneck)
end

function rangeexp_cylinder(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max=DEF_R_MAX_BURNIN, r_max_burnin=r_max, r_max_exp=r_max, migr_mode=DEF_MIGR_MODE, 
    startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE,
    z_max_burnin=DEF_X_MAX_BURNIN, z_max_exp=DEF_X_MAX_EXP, max_burnin=(NaN, NaN, z_max_burnin), max_exp=(NaN, NaN, z_max_exp), maxi=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, z_max_exp), capacity=DEF_CAPACITY,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    data_to_generate=DEF_DATA_TO_GENERATE, wld_gt1=NaN, wld_gt2=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, multiproc=true, 
    SS=true, verbose=false, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false)

    if !isa(startfill_range, Array)
        ran = ins_sq(r_max_burnin, maxi[1])
        startfill_range = [ran, ran, 1:z_max_burnin]
    end

    rangeexp(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, maxi=maxi,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt1=wld_gt1, wld_gt2=wld_gt2, wld_stats=wld_stats, name=name, bottleneck=bottleneck, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci,
        startfill_range=startfill_range, multiproc=multiproc, SS=SS, verbose=verbose, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate)
end

function rangeexp_sphere(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,startfill_iter=NaN,
    max_burnin=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, r_max_exp * 2 + 1), max_exp=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, r_max_exp * 2 + 1), maxi=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, r_max_exp * 2 + 1), capacity=DEF_CAPACITY,
     prolif_rate=DEF_PROLIF_RATE, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, SS=true, verbose=false,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, domin_coef=DEF_DOMIN_COEF, mutratelocus=false, n_loci=DEF_N_LOCI, n_sel_loci=ceil(Int,n_loci/2), loci=fill(sel_coef,n_loci),
    data_to_generate=DEF_DATA_TO_GENERATE, wld_gt1=NaN, wld_gt2=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN)

    if iterate(startfill_iter) !== nothing
        startfill_iter = sphere_lattice_iter(maxi[1]/2,maxi[2]/2,maxi[3]/2,r_max_burnin)
    end

    rangeexp(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, maxi=maxi, r_coords=[1, 2, 3], prolif_rate=prolif_rate,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt1=wld_gt1, wld_gt2=wld_gt2, wld_stats=wld_stats, name=name, bottleneck=bottleneck, capacity=capacity,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, domin_coef=domin_coef, mutratelocus=mutratelocus, n_loci=n_loci, n_sel_loci=n_sel_loci, loci=loci,
        startfill_range=startfill_range, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, SS=SS, verbose=verbose, startfill_iter=startfill_iter)
end



# Simulation functions (infinite-sites)
# ------------------------------------------------

@inline function calc_fitness_inf(person::Array{typ_float,1})
    return prod(person)
end

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
@inline function mutate_inf!(person::Array{typ_float,1}, mut_rate::Float64, 
                            n_segr_regions::Int, sel_coef::Float64, 
                            prop_of_del_muts::Float64; mutratelocus=false)
    muts_del = 0
    muts_ben = 0
    
    # Get number of mutations based on rate
    get_mutation_random = mutratelocus ? rand(Poisson(mut_rate*n_loci)) : rand(Poisson(mut_rate))
    
    @fastmath @inbounds for _ in 1:get_mutation_random
        pos_alter = rand(1:n_segr_regions)
        if rand() < prop_of_del_muts
            person[pos_alter] *= (1 - sel_coef)
            muts_del += 1
        else
            person[pos_alter] *= (1 + sel_coef)
            muts_ben += 1
        end
    end

    return muts_del, muts_ben
end

function mutate_inf!(person, wld_stats; mutratelocus=false)
    mutate_inf!(person, wld_stats["mut_rate"], wld_stats["n_segr_regions"], wld_stats["sel_coef"], wld_stats["prop_of_del_muts"];mutratelocus=mutratelocus)
end

"""
Recombines loci with a 1/2 chance. Used when building the next generation in infinite-sites expansions.

---

`person`: an individual's segr. regions (fitness) array

`n_segr_regions`: number of segregating regions
"""
@inline function crossover_inf!(person::Array{typ_float,1}, n_segr_regions::Int)
    @inbounds for i in 1:n_segr_regions
        if rand(Bool)  # More efficient than rand(1:2) == 1
            # Do nothing (keep first allele)
        else
            person[i] = person[i+n_segr_regions]
        end
    end
end

function crossover_inf!(person, wld_stats::Dict)
    crossover_inf!(person, wld_stats["n_segr_regions"])
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
@inline function mate_inf!(result::Array{typ_float,1}, ind1::Array{typ_float,1}, 
                          ind2::Array{typ_float,1}, n_segr_regions::Int, fixed_mate::Bool=false)
    @inbounds for i in 1:n_segr_regions
        # Choose alleles from parent 1
        from_first = fixed_mate || rand(Bool)
        result[i] = from_first ? ind1[i] : ind1[i+n_segr_regions]
        
        # Choose alleles from parent 2
        from_first = fixed_mate || rand(Bool)
        result[i+n_segr_regions] = from_first ? ind2[i] : ind2[i+n_segr_regions]
    end
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
function build_next_gen_inf!(g, gg, out_fields, wld_gt::typ_gt_inf, stats, wld_fitn=NaN, wld_pops=NaN, wld_mutsdel=NaN, wld_mutsben=NaN;
                        max_migr::Tuple=stats["max"], migr_mode=DEF_MIGR_MODE, r_max_migr::Int64=0, bottleneck=NaN, refl_walls=false, r_coords::Array{Int}=[1, 2], 
                        weightfitn::Bool=true, condsel::Bool=false, mutratelocus::Bool=false,
                        fixed_mate::Bool=false, premutate::Bool=false, SS::Bool=false, verbose=false)
    
    # Preallocate buffers for migration
    move_buffer = zeros(Int16, stats["wlddim"])
    wlddim = stats["wlddim"]
    #Random.seed!(1234)

    # Calculate offspring counts
    next_gen_posits, next_gen_pops = calc_offspring(wld_gt[repeat([:],wlddim)...,gg], stats)
    if verbose
        println("g=",g, ", gg=",gg,", pops=", next_gen_pops)
    end
    
    # Pre-allocate reused arrays for mating
    mate_result = ones(typ_float, stats["n_segr_regions"] * 2)
    
    # Process each deme
    for deme in next_gen_posits
        inds_at_pos = wld_gt[deme...,gg]
        
        # Calculate fitnesses once
        fitnesses = Vector{typ_float}(undef, length(inds_at_pos))
        @inbounds for i in 1:length(inds_at_pos)
            fitnesses[i] = calc_fitness_inf(inds_at_pos[i])
        end
        
        next_generation_size = next_gen_pops[deme...]
        
        if next_generation_size > 0
            # Create births
            k = 0
            @inbounds while k < next_generation_size
                # Select parents
                mom_idx, dad_idx = 1, 1
                if weightfitn
                    mom_idx = wsample(1:length(inds_at_pos), fitnesses)
                    dad_idx = wsample(1:length(inds_at_pos), fitnesses)
                else
                    mom_idx = rand(1:length(inds_at_pos))
                    dad_idx = rand(1:length(inds_at_pos))
                end
                
                mom = inds_at_pos[mom_idx]
                dad = inds_at_pos[dad_idx]
                
                # Fast path for non-conditional selection
                thing1 = rand()
                thing2 = rand()
                if verbose
                    print("Mate rand()=",thing1," rand()=", thing2," ")
                end
                if !condsel || (fitnesses[mom_idx] > thing1*maximum(fitnesses) && 
                               fitnesses[dad_idx] > thing2*maximum(fitnesses))
                    k+=1
                    # Create gametes with recombination
                    gamete_mom = copy(mom)
                    gamete_dad = copy(dad)
                    
                    # Apply recombination
                    crossover_inf!(gamete_mom, stats["n_segr_regions"])
                    crossover_inf!(gamete_dad, stats["n_segr_regions"])
                    
                    # Create zygote
                    if premutate
                        # Mutations in gametes
                        wld_mutsdel, wld_mutsben = mutate_inf!(gamete_mom, stats["mut_rate"], stats["n_segr_regions"], 
                                   stats["sel_coef"], stats["prop_of_del_muts"])
                        wld_mutsdel, wld_mutsben = mutate_inf!(gamete_dad, stats["mut_rate"], stats["n_segr_regions"], 
                                   stats["sel_coef"], stats["prop_of_del_muts"])
                        # Mating with fixed buffer
                        mate_inf!(mate_result, gamete_mom, gamete_dad, stats["n_segr_regions"], fixed_mate)
                    else
                        # Mating first
                        mate_inf!(mate_result, gamete_mom, gamete_dad, stats["n_segr_regions"], fixed_mate)
                        
                        # Then mutate
                        wld_mutsdel, wld_mutsben = mutate_inf!(mate_result, stats["mut_rate"], stats["n_segr_regions"], 
                                   stats["sel_coef"], stats["prop_of_del_muts"])
                    end

                    if out_fields["mutsdel"][1]
                        if !out_fields["mutsdel"][2]
                            G = g
                        elseif g == stats["n_gens"]
                            G = 1
                        end
                        if !isnan(wld_mutsdel) && !isnan(wld_mutsben)
                            if isnan(wld_mutsdel[deme...,G])
                                wld_mutsdel[deme...,G] = 0
                            end
                            if isnan(wld_mutsben[deme...,G])
                                wld_mutsben[deme...,G] = 0
                            end
                            wld_mutsdel[deme...,G] += wld_mutsdel
                            wld_mutsben[deme...,G] += wld_mutsben
                        end
                    end

                    # Calculate migration
                    calc_migr_dist!(move_buffer, deme, stats, max_migr, r_max_migr, r_coords)
                    
                    # Apply migration to calculate new position
                    indices = Vector{Int}(undef, wlddim)

                    for i in 1:wlddim
                        indices[i] = deme[i] + move_buffer[i]
                    end
                    
                    # Add offspring to next generation
                    push!(wld_gt[indices...,gg+1], copy(mate_result))

                    if verbose
                        println("Completed birth #",k," rand()=",rand(1))
                    end
                elseif !SS
                    k+=1
                end
            end
        end

        if out_fields["pops"][1]
            G = out_fields["pops"][2] ? g+1 : 1
            wld_pops[deme...,G]=length(wld_gt[deme...,gg+1])
        end
        if out_fields["fitn"][1]
            G = out_fields["fitn"][2] ? g+1 : 1
            wld_fitn[deme...,G]=mean(calc_fitness_inf.(wld_gt[deme...,gg+1]))
        end
    end

    if verbose
        println("--------------------")
    end

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
function create_empty_world_inf(maxi=(DEF_X_MAX, DEF_Y_MAX); min=(1, 1), name="sim_$(Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"))", capacity=DEF_CAPACITY,
    prolif_rate=DEF_PROLIF_RATE, n_segr_regions=DEF_N_SEGR_REGIONS, regions=fill(sel_coef,n_segr_regions),
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, migr_mode=DEF_MIGR_MODE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)

    wld_gt = (typ_gt_inf)(undef, maxi...)
    for k in CartesianIndices(maxi)
        wld_gt[k] = Array{typ_float,1}[]
    end
    #println("Created ",wld_gt)
    #println(size(wld_gt))

    wld_stats = OrderedDict(
        "name" => name,
        "max" => maxi,
        "capacity" => capacity,
        "prolif_rate" => prolif_rate,
        "regions" => regions,
        "n_segr_regions" => n_segr_regions,
        "mut_rate" => mut_rate,
        "migr_rate" => migr_rate,
        "migr_mode" => migr_mode,
        "sel_coef" => sel_coef,
        "prop_of_del_muts" => prop_of_del_muts,
        "wlddim" => length(maxi)
        #"rangeexps" => []
    )

    return wld_gt, wld_stats
end

"""
Fills random demes within given monosome arrays with infinite-sites individuals. Usually used after an empty world is created.

---

`wld`: a spatial array of demes that contain individuals' segr. regions [Float] arrays

`wld_stats`: world stats Dict

`fill`: an array of Int ranges of the coordinates that define the area within which to fill

`n_demes_to_fill`: number of demes to fill

"""
function fill_random_demes_inf(wld_gt, wld_stats, fill_iter, n_demes_to_fill=DEF_N_DEMES_STARTFILL; fill_range=NaN, redims=NaN, verbose=false)
    if !isnothing(fill_iter)
        possible_init_coords = collect(fill_iter)
        #wld_stats["startfill"] = collect(fill_iter)
    else
        possible_init_coords = [collect(x) for x in Iterators.product(fill_range...)]
        wld_stats["startfill"] = copy(fill_range)
    end
    if length(possible_init_coords)<n_demes_to_fill
        println("n_demes_to_fill is more than existing demes!")
    end
    init_coords = sample(possible_init_coords, n_demes_to_fill; replace=false)
    if verbose
        println("Filling ",wld_gt,".")
    end
    for coord in init_coords
        for _ in 1:wld_stats["capacity"]
            push!(wld_gt[coord...,redims...], ones(typ_float, wld_stats["n_segr_regions"] * 2))

            #println(size(wld_gt))
        end
    end
    if verbose
        println("Filled ",wld_gt,".")
    end
    
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
function rangeexp_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; max_burnin=(DEF_X_MAX_BURNIN, DEF_Y_MAX), max_exp=(DEF_X_MAX_EXP, DEF_Y_MAX), maxi=(DEF_X_MAX, DEF_Y_MAX), migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, r_max_burnin=0, r_max_exp=0, r_coords=[1, 2], capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, 
    multiproc=true, weightfitn=true, condsel=false, fixed_mate=false, premutate=false, mutratelocus=false,startfill_iter=nothing,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, verbose=false,
    n_segr_regions=DEF_N_SEGR_REGIONS, regions=fill(sel_coef,n_segr_regions), startfill_range=NaN, wld_gt=NaN, wld_stats=NaN, SS=true)

    if isnan(wld_gt)
        is_fill_random_demes = true
        #println("No world provided. Creating a new world.")
        wld_gt, wld_stats = create_empty_world_inf(maxi; name=name, capacity=capacity, prolif_rate=prolif_rate,
            mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts,
            n_segr_regions=n_segr_regions, regions=regions, migr_mode=migr_mode)
        if isnothing(startfill_iter) && !isa(startfill_range, Array) && !any(isnan, max_burnin)
            startfill_range = [1:upper for upper in max_burnin]
        end
    end
    println("max_burnin ",max_burnin)
    println("startfill_range ",startfill_range)
    println("startfill_iter ",startfill_iter)
    wlddim = wld_stats["wlddim"]
    wld_stats["max_burnin"] = max_burnin
    wld_stats["max_exp"] = max_exp
    wld_stats["n_gens_burnin"] = n_gens_burnin
    wld_stats["n_gens_exp"] = n_gens_exp
    n_gens = n_gens_burnin + n_gens_exp
    wld_stats["n_gens"] = n_gens

    function check_what_output(letter)
        return (occursin(letter, data_to_generate),occursin(letter*"l", data_to_generate))
    end
    out_fields = OrderedDict{String, Any}("gt" => check_what_output("G"), "fitn" => check_what_output("F"), "pops" => check_what_output("P"),
    "mutsdel" => check_what_output("M"), "mutsben" => check_what_output("M"))

    function tsk_re(_proc_no=NaN)
        cond_n_gens = out_fields["gt"][1] && out_fields["gt"][2] ? n_gens+1 : 2
        rewld_gt_local = (typ_gt_inf)(undef, wld_stats["max"]..., cond_n_gens)

        rewld_gt_local[repeat([:],wlddim)...,1] = deepcopy(wld_gt)
        for og in 2:cond_n_gens
            rewld_gt_local[repeat([:],wlddim)...,og] = deepcopy(wld_gt)
        end
        if is_fill_random_demes
            fill_random_demes_inf(rewld_gt_local, wld_stats, startfill_iter; redims=(1), fill_range=startfill_range, verbose=verbose)
        end
        
        if out_fields["fitn"][1]
            cond_n_gens = out_fields["fitn"][2] ? n_gens+1 : 1
            rewld_fitn_local = Array{typ_float}(undef, wld_stats["max"]..., cond_n_gens)
            fill!(rewld_fitn_local,NaN)
            if out_fields["fitn"][2]
                for deme in CartesianIndices(maxi)
                    rewld_fitn_local[deme,1] = mean(calc_fitness_inf.(rewld_gt_local[deme,1]))
                end
            end
        else
            rewld_fitn_local = NaN
        end
        if out_fields["pops"][1]
            cond_n_gens = out_fields["pops"][2] ? n_gens+1 : 1
            rewld_pops_local = Array{typ_float}(undef, wld_stats["max"]..., cond_n_gens)
            fill!(rewld_pops_local,NaN)
            if out_fields["pops"][2]
                for deme in CartesianIndices(maxi)
                    rewld_pops_local[deme,1] = length(rewld_gt_local[deme,1])
                end
            end
        else
            rewld_pops_local = NaN
        end
        if out_fields["mutsdel"][1]
            cond_n_gens = out_fields["mutsdel"][2] ? n_gens+1 : 1
            rewld_mutsdel_local = Array{typ_float}(undef, wld_stats["max"]..., wld_stats["n_loci"], cond_n_gens)
            rewld_mutsben_local = Array{typ_float}(undef, wld_stats["max"]..., wld_stats["n_loci"], cond_n_gens)
            fill!(rewld_mutsdel_local,NaN)
            fill!(rewld_mutsben_local,NaN)
        else
            rewld_mutsdel_local = NaN
            rewld_mutsben_local = NaN
        end

        @inbounds for g in 1:n_gens
            if g <= n_gens_burnin
                max_migr = max_burnin
                r_max_migr = r_max_burnin
            else
                max_migr = max_exp
                r_max_migr = r_max_exp
            end

            condi = out_fields["gt"][1] && out_fields["gt"][2]
            gg = condi ? g : 1

            if !condi
                for deme in CartesianIndices(maxi)
                    empty!(rewld_gt_local[deme,2])
                end
            end

            build_next_gen_inf!(g, gg, out_fields, rewld_gt_local, wld_stats, rewld_fitn_local, rewld_pops_local, rewld_mutsdel_local, rewld_mutsben_local;
                max_migr=max_migr, migr_mode=migr_mode, r_max_migr=r_max_migr, bottleneck=bottleneck, r_coords=r_coords,
                mutratelocus=mutratelocus, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate, SS=SS, verbose=verbose)

            if !condi
                rewld_gt_local[repeat([:],wlddim)...,1],rewld_gt_local[repeat([:],wlddim)...,2] = rewld_gt_local[repeat([:],wlddim)...,2],rewld_gt_local[repeat([:],wlddim)...,1]
            end
        end

        return Dict{String, Any}("gt" => rewld_gt_local, "fitn" => rewld_fitn_local, "pops" => rewld_pops_local,
        "mutsdel" => rewld_mutsdel_local, "mutsben" => rewld_mutsben_local)
    end

    extra_dim = 1

    if n_re > 1
        extra_dim += 1
        if multiproc
            npcs = workers()
            println("Running $n_re replicates on processes $npcs.")
            dicts_out,t = @timed pmap(tsk_re, 1:n_re)
            println("Simulated in ",t," s.")
            wld_stats["simtime"]=t
        else
            println("Using all threads.")
            dicts_out = ThreadsX.map(p->tsk_re(), 1:n_re)
        end
    else
        dicts_out,t = @timed [tsk_re()]
        println("Simulated in ",t," s.")
        wld_stats["simtime"]=t
    end

    res = OrderedDict{String, Any}("stats" => wld_stats)

    for key in keys(out_fields)
        if out_fields[key][1]
            if out_fields[key][2]
                if key=="mutsdel" || key=="mutsben" || key=="gt"
                    res[key] = cat([dictus[key] for dictus in dicts_out]...,dims=wlddim+extra_dim+1)
                else
                    res[key] = cat([dictus[key] for dictus in dicts_out]...,dims=wlddim+extra_dim+1)
                end
            else
                if key=="gt"
                    res[key] = cat([dictus[key][repeat([:],wlddim)...,1] for dictus in dicts_out]...,dims=wlddim+extra_dim)
                elseif key=="mutsdel" || key=="mutsben"
                    res[key] = cat([dictus[key][repeat([:],wlddim)...,1] for dictus in dicts_out]...,dims=wlddim+extra_dim)
                else
                    res[key] = cat([dictus[key][repeat([:],wlddim)...,1] for dictus in dicts_out]...,dims=wlddim+extra_dim)
                end
            end
        end
    end

    return res
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

`multiproc`: if **true**, distribute to threads

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
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, mutratelocus=false, weightfitn=true, condsel=false, fixed_mate=false, premutate=false, SS=true, verbose=false,
    n_segr_regions=DEF_N_SEGR_REGIONS, regions=fill(sel_coef,n_segr_regions),
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, multiproc=true, wld_gt=NaN, wld_stats=NaN, capacity=DEF_CAPACITY)

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; max_burnin=(x_max_burnin,), max_exp=(x_max_exp,), maxi=(x_max_exp,), migr_mode=migr_mode,
        data_to_generate=data_to_generate, name=name, bottleneck=bottleneck, capacity=capacity, prolif_rate=prolif_rate, multiproc=multiproc, weightfitn=weightfitn,
        condsel=condsel, fixed_mate=fixed_mate, premutate=premutate, mutratelocus=mutratelocus,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, SS=SS, verbose=verbose,
        n_segr_regions=n_segr_regions, regions=regions, startfill_range=startfill_range, 
        wld_gt=wld_gt, wld_stats=wld_stats)
end

const rangeexp_1d_inf = rangeexp_ray_inf

function rangeexp_linear_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, maxi=(r_max_exp * 2 + 1,),
    migr_mode=DEF_MIGR_MODE, startfill_range=[(1-ceil(Int,r_max_burnin/2)+r_max_exp):(1+ceil(Int,r_max_burnin/2)+r_max_exp)], prolif_rate=DEF_PROLIF_RATE, max_exp=(NaN,NaN), max_burnin=(NaN,NaN),
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, premutate=false, n_segr_regions=DEF_N_SEGR_REGIONS,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, multiproc=true, wld_gt=NaN, wld_stats=NaN, capacity=DEF_CAPACITY)

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, maxi=maxi, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate, n_segr_regions=n_segr_regions, 
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt=wld_gt, wld_stats=wld_stats, name=name, bottleneck=bottleneck, multiproc=multiproc, r_coords=[1])
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
    max_burnin=(x_max_burnin, y_max), max_exp=(x_max_exp, y_max), maxi=(x_max_exp, y_max), capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, n_segr_regions=DEF_N_SEGR_REGIONS,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, premutate=false,
    data_to_generate=DEF_DATA_TO_GENERATE, wld_gt=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=("midhole at x=", x_max_burnin * 2), SS=true, verbose=false)

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; max_burnin=max_burnin, max_exp=max_exp, maxi=maxi, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate, n_segr_regions=n_segr_regions,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt=wld_gt, wld_stats=wld_stats, name=name, bottleneck=bottleneck, SS=SS, verbose=verbose)
end

function rangeexp_disk_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, maxi=(r_max_exp * 2 + 1, r_max_exp * 2 + 1),
    capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false,
    premutate=false, n_segr_regions=DEF_N_SEGR_REGIONS,
    data_to_generate=DEF_DATA_TO_GENERATE, wld_gt=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, max_exp=(NaN,NaN), max_burnin=(NaN,NaN), startfill_range=NaN, SS=true, verbose=false)

    if !isa(startfill_range, Array)
        ran = ins_sq(r_max_burnin, maxi[1])
        startfill_range = [ran, ran]
#=         if isnan(max_burnin)
            max_burnin = ((1-ceil(Int,r_max_burnin/2)+r_max_exp),(1+ceil(Int,r_max_burnin/2)+r_max_exp))
        end =#
    end

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, maxi=maxi, startfill_range=startfill_range, capacity=capacity, prolif_rate=prolif_rate,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate, n_segr_regions=n_segr_regions,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt=wld_gt, wld_stats=wld_stats, name=name, bottleneck=bottleneck, SS=SS, verbose=verbose)
end

function rangeexp_cylinder_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max=DEF_R_MAX_BURNIN, r_max_burnin=r_max, r_max_exp=r_max, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE,
    z_max_burnin=DEF_X_MAX_BURNIN, z_max_exp=DEF_X_MAX_EXP, max_burnin=(NaN, NaN, z_max_burnin), max_exp=(NaN, NaN, z_max_exp), maxi=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, z_max_exp), capacity=DEF_CAPACITY,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, premutate=false, n_segr_regions=DEF_N_SEGR_REGIONS,
    data_to_generate=DEF_DATA_TO_GENERATE, wld_gt=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, SS=true, verbose=false)

    if !isa(startfill_range, Array)
        ran = ins_sq(r_max_burnin, maxi[1])
        startfill_range = [ran, ran, 1:z_max_burnin]
    end

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, maxi=maxi,
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt=wld_gt, wld_stats=wld_stats, name=name, bottleneck=bottleneck, capacity=capacity, prolif_rate=prolif_rate, n_segr_regions=n_segr_regions,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate,
        startfill_range=startfill_range, SS=SS, verbose=verbose)
end

function rangeexp_sphere_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, startfill_iter=NaN,
    max_burnin=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, r_max_exp * 2 + 1), max_exp=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, r_max_exp * 2 + 1), maxi=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, r_max_exp * 2 + 1), capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, n_segr_regions=DEF_N_SEGR_REGIONS,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, weightfitn=true, condsel=false, fixed_mate=false, premutate=false,
    data_to_generate=DEF_DATA_TO_GENERATE, wld_gt=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, SS=true, verbose=false)

    if iterate(startfill_iter) !== nothing
        startfill_iter = sphere_lattice_iter(maxi[1]/2,maxi[2]/2,maxi[3]/2,r_max_burnin)
    end

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; r_max_burnin=r_max_burnin, r_max_exp=r_max_exp, max_burnin=max_burnin, max_exp=max_exp, maxi=maxi, r_coords=[1, 2, 3],
        migr_mode=migr_mode, data_to_generate=data_to_generate, wld_gt=wld_gt, wld_stats=wld_stats, name=name, bottleneck=bottleneck, capacity=capacity, prolif_rate=prolif_rate, n_segr_regions=n_segr_regions,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate,
        startfill_range=startfill_range, SS=SS, verbose=verbose, startfill_iter=startfill_iter)
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
    res = Array{typ_float}(undef, 0)
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
Finds the average value of `data` between all demes at the expansion front, for all replicates.

---

`data`: array with dimensions (space + time + 1)

`n_gens`: number of generations

`greaterzero`: if **true**, **>0** values are considered when determining the front (**>=0** values if **false**)

`oneside`: if **true**, approach only from one side (i.e. from the positive direction in strip expansions)

`divide`: if **false**, find the sum instead of average

---

Output: array of averages of `data` for every generation and every replicate
"""
function average_front(data, n_gens, x_max; greaterzero=false, oneside=false, divide=true)
    n_re = size(data, 3)
    av_arr = Array{typ_float}(undef, n_gens, n_re)

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
    av_arr = Array{typ_float}(undef, n_gens, n_re)

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
    av_arr = Array{typ_float}(undef, n_gens, n_re)

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
Finds the front array of `dataname` in `re`, for every replicate.

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

`oneside`: if **true**, approach only from one side (i.e. from the positive direction in strip expansions)

---

Output: front array of the same dimensions as `re[dataname]` (space + time + 1)
"""
function front_array(re, dataname; oneside=false)
    front_array(re[dataname], re["stats"]["n_gens"], re["stats"]["max"]...; oneside=oneside)
end

# 1D
function front_array(data, n_gens, x_max; oneside=false)
    n_re = size(data, 3)
    front_arr = fill(NaN, x_max, n_gens, n_re)

    for i in 1:n_re, j in 1:n_gens
        frontier = x_max
        while frontier != 1 && isnan(data[frontier, j, i])
            frontier -= 1
        end
        if !isnan(data[frontier, j, i])
            front_arr[frontier, j, i] = data[frontier, j, i]
        end
        if !oneside
            frontier = 1
            while frontier != x_max && isnan(data[frontier, j, i])
                frontier += 1
            end
            if !isnan(data[frontier, j, i])
                front_arr[frontier, j, i] = data[frontier, j, i]
            end
        end
    end
    return front_arr
end

# 2D
function front_array(data::Array, n_gens, x_max, y_max; oneside=false)
    n_re = size(data, 4)
    front_arr = fill(NaN, x_max, y_max, n_gens)
    for i in 1:n_re, j in 1:n_gens
        # scanning every y: side 1
        for _y in 1:y_max
            frontier_x = x_max
            while frontier_x != 1 && isnan(data[frontier_x, _y, j, i])
                frontier_x -= 1
            end
            if !isnan(data[frontier_x, _y, j, i])
                front_arr[frontier_x, _y, j, i] = data[frontier_x, _y, j, i]
            end
        end
        # scanning every y: side 2
        if !oneside
            for _y in 1:y_max
                frontier_x = 1
                while frontier_x != x_max && isnan(data[frontier_x, _y, j, i])
                    frontier_x += 1
                end
                if !isnan(data[frontier_x, _y, j, i])
                    front_arr[frontier_x, _y, j, i] = data[frontier_x, _y, j, i]
                end
            end
        end

        if !oneside
            # scanning every x: side 1
            for _x in 1:x_max
                frontier_y = y_max
                while frontier_y != 1 && isnan(data[_x, frontier_y, j, i])
                    frontier_y -= 1
                end
                if !isnan(data[_x, frontier_y, j, i])
                    front_arr[_x, frontier_y, j, i] = data[_x, frontier_y, j, i]
                end
            end
            # scanning every x: side 2
            for _x in 1:x_max
                frontier_y = 1
                while frontier_y != y_max && isnan(data[_x, frontier_y, j, i])
                    frontier_y += 1
                end
                if !isnan(data[_x, frontier_y, j, i] > 0)
                    front_arr[_x, frontier_y, j, i] = data[_x, frontier_y, j, i]
                end
            end
        end
    end
    return front_arr
end

# 3D
function front_array(data::Array, n_gens, x_max, y_max, z_max; oneside=false)
    n_re = size(data, 5)
    front_arr = fill(NaN, x_max, y_max, z_max, n_gens, n_re)

    for i in 1:n_re, j in 1:n_gens
        # scanning every xy: side 1
        for _x in 1:x_max, _y in 1:y_max
            frontier_z = z_max
            while frontier_z != 1 && isnan(data[_x, _y, frontier_z, j, i])
                frontier_z -= 1
            end
            if !isnan(data[_x, _y, frontier_z, j, i])
                front_arr[_x, _y, frontier_z, j, i] = data[_x, _y, frontier_z, j, i]
            end
        end
        # scanning every xy: side 2
        if !oneside
            for _x in 1:x_max, _y in 1:y_max
                frontier_z = 1
                while frontier_z != z_max && isnan(data[_x, _y, frontier_z, j, i])
                    frontier_z += 1
                end
                if !isnan(data[_x, _y, frontier_z, j, i])
                    front_arr[_x, _y, frontier_z, j, i] = data[_x, _y, frontier_z, j, i]
                end
            end
        end

        if !oneside
            # scanning every yz: side 1
            for _y in 1:y_max, _z in 1:z_max
                frontier_x = x_max
                while frontier_x != 1 && isnan(data[frontier_x, _y, _z, j, i])
                    frontier_x -= 1
                end
                if !isnan(data[frontier_x, _y, _z, j, i])
                    front_arr[frontier_x, _y, _z, j, i] = data[frontier_x, _y, _z, j, i]
                end
            end
            # scanning every yz: side 2
            for _y in 1:y_max, _z in 1:z_max
                frontier_x = 1
                while frontier_x != x_max && isnan(data[frontier_x, _y, _z, j, i])
                    frontier_x += 1
                end
                if !isnan(data[frontier_x, _y, _z, j, i])
                    front_arr[frontier_x, _y, _z, j, i] = data[frontier_x, _y, _z, j, i]
                end
            end

            # scanning every xz: side 1
            for _x in 1:x_max, _z in 1:z_max
                frontier_y = y_max
                while frontier_y != 1 && isnan(data[_x, frontier_y, _z, j, i])
                    frontier_y -= 1
                end
                if !isnan(data[_x, frontier_y, _z, j, i])
                    front_arr[_x, frontier_y, _z, j, i] = data[_x, frontier_y, _z, j, i]
                end
            end
            # scanning every yz: side 2
            for _x in 1:x_max, _z in 1:z_max
                frontier_y = 1
                while frontier_y != y_max && isnan(data[_x, frontier_y, _z, j, i, i])
                    frontier_y += 1
                end
                if !isnan(data[_x, frontier_y, _z, j, i])
                    front_arr[_x, frontier_y, _z, j, i] = data[_x, frontier_y, _z, j, i]
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

`n_gens`: generation to stop at

---

Output: averaged time series
"""
function average_ts(ts_arr, n_gens=size(ts_arr))
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

function af_A(re,deme,locus,gen=re["stats"]["n_gens"],re_index=1)
    return (2*re["cAA"][deme...,locus,gen,re_index]+re["cAa"][deme...,locus,gen,re_index])/re["stats"]["n_loci"]
end

# TO-DO!!! Test if this == re["A"]

function af_a(re,deme,locus,gen=re["stats"]["n_gens"],re_index=1)
    return (2*re["caa"][deme...,locus,gen,re_index]+re["cAa"][deme...,locus,gen,re_index])/re["stats"]["n_loci"]
end

function twopq(re,deme,locus::Int,gen=re["stats"]["n_gens"],re_index=1)
    a = af_A(re,deme,locus,gen,re_index)
    return a*(1-a)*2
end

function twopq(re,deme,loci::Array,gen=re["stats"]["n_gens"],re_index=1)
    return mean([af_A(re,deme,l,gen,re_index)*(1-af_A(re,deme,l,gen,re_index))*2 for l in loci])
end

function twopq(re,deme,loci::UnitRange=1:test["stats"]["n_loci"],gen=re["stats"]["n_gens"],re_index=1)
    return mean([af_A(re,deme,l,gen,re_index)*(1-af_A(re,deme,l,gen,re_index))*2 for l in loci])
end

function pbar(re,deme_arr,loci=1:test["stats"]["n_loci"],gen=re["stats"]["n_gens"],re_index=1)
    return sum([af_A(re,deme,loci,gen,re_index)*re["pops"][deme...,gen,re_index] for deme in deme_arr])/sum([re["pops"][deme...,gen,re_index] for deme in deme_arr])
end

function H_S(re,deme_arr,loci=1:test["stats"]["n_loci"],gen=re["stats"]["n_gens"],re_index=1)
    return sum([twopq(re,deme,loci,gen,re_index)*re["pops"][deme...,gen,re_index] for deme in deme_arr])/sum([re["pops"][deme...,gen,re_index] for deme in deme_arr])
end

function H_T(re,deme_arr,loci=1:test["stats"]["n_loci"],gen=re["stats"]["n_gens"],re_index=1)
    return mean([2*pbar(re,deme_arr,l)*(1-pbar(re,deme_arr,l)) for l in loci])
end

function Ne(data,range=:)
    return [harmmean(k[range]) for k in data]
end

function F_ST(re,deme_arr,loci=1:re["stats"]["n_loci"],gen=re["stats"]["n_gens"],re_index=1;verbose=true)
    HT = H_T(re,deme_arr,loci,gen,re_index)
    println(HT)
    HS = H_S(re,deme_arr,loci,gen,re_index)
    println(HS)
    return (HT-HS)/HT
end

function Ne_avcuml(data,n_gens_burnin,range)
    N_e_av_cuml_set = []
    for i in range
        push!(N_e_av_cuml_set,mean([harmmean(k[n_gens_burnin:i]) for k in data]))
    end
    return N_e_av_cuml_set
end

function F_ST_old(data)
    F_ST_set = [1/(k+1) for k in data]
    return mean(F_ST_set)
end

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

function out_show(mtrx)
    open("out_show.txt", "w") do io
        show(io, "text/plain",  mtrx)
    end
end
# ------------------------------------------------

println("RESK successfully loaded.")
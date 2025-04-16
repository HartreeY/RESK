DEF_X_MAX = 500
DEF_Y_MAX = 10
DEF_Z_MAX = 10
DEF_N_DEMES_STARTFILL = 5
DEF_CAPACITY = 100
DEF_PROLIF_RATE = 2
DEF_MUT_RATE = 0.05 # Genome-wide mutation rate
DEF_MIGR_RATE = 0.05 # Migration rate
DEF_SEL_COEF = 0.005 # Selection coefficient
DEF_PROP_OF_DEL_MUTS = 0.9
DEF_N_SEGR_REGIONS = 20
DEF_X_MAX_BURNIN = 5
DEF_R_MAX_BURNIN = 3 # Radius that bounds the burn-in area
DEF_N_GENS_BURNIN = 10 # Number of burn-in generations
DEF_X_MAX_EXP = DEF_X_MAX
DEF_Y_MAX_EXP = DEF_Y_MAX
DEF_R_MAX_EXP = 20 # Radius that bounds the expansion area
DEF_N_GENS_EXP = 40 # Number of expansion generations
DEF_MIGR_MODE = "ort" # Migration mode
DEF_DATA_TO_GENERATE = "FP"

using StatsBase, Distributions, Random, SpecialFunctions, Serialization, Dates, DataStructures, Distributed

MIGR_PROBS = [
    Dict(["ort" => (1, 0)]), # 1D
    Dict(["ort" => (1, 0), "hex" => (1, 0), "all" => (1 / 2, 1 / 2), "buffon1" => (2 / pi, 1 / pi), "buffon2" => (4 / 3 / pi, 1 / 3 / pi), "buffon3" => (0.4244132, 0.21221), "diag1/2" => (2 / 3, 1 / 3)]), # 2D. To add "hex"!
    Dict(["ort" => (1, 0), "all" => (1 / 2, 1 / 2), "buffon1" => (2 / pi, 1 / pi), "buffon2" => (4 / 3 / pi, 1 / 3 / pi), "buffon3" => (0.4244132, 0.21221), "diag1/2" => (2 / 3, 1 / 3)]) # 3D. To add "hex"! To confirm Buffon for 3d!
]
MIGR_DIRS_ORT = [
    [[1], [-1]], # 1D
    [[-1, 0], [0, -1], [0, 1], [1, 0]], # 2D
    [[-1, 0, 0], [1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, -1], [0, 0, 1]] # 3D
]
MIGR_DIRS_DIAG = [
    [], # 1D
    [[-1, -1], [-1, 1], [1, -1], [1, 1]], # 2D
    [[-1, -1, -1], [-1, -1, 0], [-1, -1, 1], [-1, 0, -1], [-1, 0, 1], [-1, 1, -1], [-1, 1, 0], [-1, 1, 1], [0, -1, -1], [0, -1, 0], [0, -1, 1], [0, 0, 0], [0, 1, -1], [0, 1, 1], [1, -1, -1], [1, -1, 0], [1, -1, 1], [1, 0, -1], [1, 0, 1], [1, 1, -1], [1, 1, 0], [1, 1, 1]] # 3D
]
MIGR_DIRS_HEX = [
    [], # 1D
    [[-1, 0], [0, -1], [-1, 1], [0, 1], [1, 0], [1, 1]], # 2D
    [] # 3D
]

# Types
const typ_float = Float32
const typ_int = Int32
const typ_gt = Array{Array{Array{Bool}}}
const typ_gt_inf = Array{Array{Array{typ_float}}}

function maxdef(a)
    return maximum(filter(!isnan,a))
end

function mindef(a)
    return minimum(filter(!isnan,a))
end

function repl(data,n)
    return data[repeat([:], length(size(data))-1)...,n]
end

function li(data,n)
    return repl(data,n)
end

function ins_sq(r_max_burnin, r_max_exp)
    return (trunc(Int, 1 + r_max_exp - r_max_burnin * 0.666)):(trunc(Int, 1 + r_max_exp + r_max_burnin * 0.666))
end

function ins_cb(r_max_burnin, r_max_exp)
    return (trunc(Int, 1 + r_max_exp - r_max_burnin * 0.577)):(trunc(Int, 1 + r_max_exp + r_max_burnin * 0.577))
end

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

function mutate_inf(person, mut_rate, n_segr_regions, sel_coef, prop_of_del_muts; mutratelocus=false)
    muts_del = 0
    muts_ben = 0
    
    get_mutation_random = mutratelocus ? rand(Poisson(mut_rate*n_loci)) : rand(Poisson(mut_rate))
    
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

function mutate_inf(person, wld_stats; mutratelocus=false)
    mutate_inf(person, wld_stats["mut_rate"], wld_stats["n_segr_regions"], wld_stats["sel_coef"], wld_stats["prop_of_del_muts"])
end

function crossover_inf(person, n_segr_regions::Int)
    for i in 1:n_segr_regions
        lr = rand(1:2)
        person[i] = lr == 1 ? person[i] : person[i+n_segr_regions]
    end
end

function crossover_inf(person, wld_stats::Dict)
    crossover_inf(person, wld_stats["n_segr_regions"])
end

function mate_inf(ind1, ind2, n_segr_regions; fixed_mate=false)
    lr1 = fixed_mate || (rand(1:2) == 1) ? (1:n_segr_regions) : ((n_segr_regions+1):(n_segr_regions*2))
    lr2 = fixed_mate || (rand(1:2) == 1) ? (1:n_segr_regions) : ((n_segr_regions+1):(n_segr_regions*2))
    return vcat(ind1[lr1], ind2[lr2])
end

function mate_cond(mom_fit,dad_fit,max_fitness)
    return (mom_fit > rand()*max_fitness) & (dad_fit > rand()*max_fitness)
end

function build_next_gen_inf(wld_gt, wld_stats, fitn_out=false, pops_out=false, muts_out=false;
    max_migr=NaN, migr_mode=DEF_MIGR_MODE, bottleneck=NaN, refl_walls=false, r_max_migr=0, r_coords=[1, 2], weightfitn=true, condsel=false, fixed_mate=false, premutate=false,  mutratelocus=false)

    wlddim = wld_stats["wlddim"]

    # Determine the number of offspring for each deme
    next_gen_posits, next_gen_pops = calc_offspring(wld_gt, wld_stats)

    # Define the world (as an array [=demes] of arrays [=individs] of a Float array [=segr. regions]) and the data arrays in the next generation
    wld_gt_next = Array{Array{Array{typ_float}},wlddim}(undef, wld_stats["max"]...) #deepcopy(wld_gt)
    for k in Iterators.product([1:n for n in wld_stats["max"]]...)
        wld_gt_next[k...] = Array{typ_float,1}[]
    end
    
    wld_fitn_next = NaN
    wld_pops_next = NaN
    wld_mutsdel_next = NaN
    wld_mutsben_next = NaN
    all_birth_count = 0

    # Fill the next generation habitat
    if fitn_out
        wld_fitn_next = Array{typ_float}(undef, wld_stats["max"]...)
        fill!(wld_fitn_next, NaN)
    end
    if pops_out
        wld_pops_next = Array{typ_float}(undef, wld_stats["max"]...)
        fill!(wld_pops_next, NaN)
    end
    if muts_out # cnts?
        wld_mutsdel_next = Array{typ_float}(undef, wld_stats["max"]...)
        wld_mutsben_next = Array{typ_float}(undef, wld_stats["max"]...)
        fill!(wld_mutsdel_next, NaN)
        fill!(wld_mutsben_next, NaN)
    end


    for deme in next_gen_posits
        inds_at_pos = wld_gt[deme...]
        fitns = prod.(inds_at_pos)

        if fitn_out
            wld_fitn_next[deme...] = mean(fitns)
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
                        wld_mutsdel, wld_mutsben = mutate_inf(gamete_mom, wld_stats; mutratelocus=mutratelocus)
                        wld_mutsdel, wld_mutsben = mutate_inf(gamete_dad, wld_stats; mutratelocus=mutratelocus)
                        mate_result = mate_inf(gamete_mom, gamete_dad, wld_stats["n_segr_regions"]; fixed_mate=fixed_mate)
                        
                    else
                        mate_result = mate_inf(gamete_mom, gamete_dad, wld_stats["n_segr_regions"]; fixed_mate=fixed_mate)
                        wld_mutsdel, wld_mutsben = mutate_inf(mate_result, wld_stats; mutratelocus=mutratelocus)
                    end

                    if muts_out
                        if isnan(wld_mutsdel_next[deme...])
                            wld_mutsdel_next[deme...] = 0
                        end
                        if isnan(wld_mutsben_next[deme...])
                            wld_mutsben_next[deme...] = 0
                        end
                        wld_mutsdel_next[deme...] += wld_mutsdel
                        wld_mutsben_next[deme...] += wld_mutsben
                    end

                    move = calc_migr_dist(deme, wld_stats, migr_mode, bottleneck, max_migr, refl_walls, r_max_migr, r_coords)
                    
                    indices = [deme[i] + move[i] for i in 1:wlddim]
                    if !isassigned(wld_gt_next, indices...)
                        wld_gt_next[indices...] = []
                    end
                    push!(wld_gt_next[indices...], mate_result)

                    birth_count += 1
                    all_birth_count += 1
                end
            end

            if pops_out
                wld_pops_next[deme...] = birth_count
            end
        end
    end

    return wld_gt_next, wld_fitn_next, wld_pops_next, wld_mutsdel_next, wld_mutsben_next
end


function create_empty_world_inf(maxi=(DEF_X_MAX, DEF_Y_MAX); min=(1, 1), name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), capacity=DEF_CAPACITY,
    prolif_rate=DEF_PROLIF_RATE, n_segr_regions=DEF_N_SEGR_REGIONS, regions=fill(sel_coef,n_segr_regions),
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, migr_mode=DEF_MIGR_MODE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)

    wld_gt = (typ_gt_inf)(undef, maxi...)

    for k in Iterators.product([1:n for n in maxi]...)
        wld_gt[k...] = Array{typ_float,1}[]
    end

    wld_stats = Dict(
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
    )

    return wld_gt, wld_stats
end

function fill_random_demes_inf(wld_gt, wld_stats, fill::Vector{UnitRange{Int64}}, n_demes_to_fill=DEF_N_DEMES_STARTFILL; redims=NaN)

    possible_init_coords = [collect(x) for x in Iterators.product(fill...)]
    init_coords = sample(possible_init_coords, n_demes_to_fill; replace=false)

    for coord in init_coords
        if !isassigned(wld_gt, coord...)
            wld_gt[coord...,redims...] = []
        end
        for _ in 1:wld_stats["capacity"]
            push!(wld_gt[coord...,redims...], ones(wld_stats["n_segr_regions"] * 2))
        end
    end

    wld_stats["startfill"] = copy(fill)
    wld_stats["n_demes_startfill"] = n_demes_to_fill
end


function rangeexp_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; max_burnin=(DEF_X_MAX_BURNIN, DEF_Y_MAX), max_exp=(DEF_X_MAX_EXP, DEF_Y_MAX), maxi=(DEF_X_MAX, DEF_Y_MAX), migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, r_max_burnin=0, r_max_exp=0, r_coords=[1, 2], capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, 
    multiproc=true, weightfitn=false, condsel=true, fixed_mate=true, premutate=false, mutratelocus=false,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, n_segr_regions=DEF_N_SEGR_REGIONS, regions=fill(sel_coef,n_segr_regions), startfill_range=NaN, wld_gt=NaN, wld_stats=NaN)

    if isnan(wld_gt)
        is_fill_random_demes = true
        wld_gt, wld_stats = create_empty_world_inf(maxi; name=name, capacity=capacity, prolif_rate=prolif_rate,
            mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts,
            n_segr_regions=n_segr_regions, regions=regions, migr_mode=migr_mode)
        if !isa(startfill_range, Array) && !any(isnan, max_burnin)
            startfill_range = [1:upper for upper in max_burnin]
        end
    end

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
        cond_n_gens = out_fields["gt"][1] && !out_fields["gt"][2] ? 1 : n_gens
        rewld_gt_local = (typ_gt_inf)(undef, wld_stats["max"]..., cond_n_gens)

        rewld_gt_local[repeat([:],wlddim)...,1] = deepcopy(wld_gt)
        if is_fill_random_demes
            fill_random_demes_inf(rewld_gt_local, wld_stats, startfill_range; redims=(1))
        end
        
        if out_fields["fitn"][1]
            cond_n_gens = out_fields["fitn"][2] ? 1 : n_gens
            rewld_fitn_local = Array{typ_float}(undef, wld_stats["max"]..., cond_n_gens)
            fill!(rewld_fitn_local,NaN)
        else
            rewld_fitn_local = NaN
        end
        if out_fields["pops"][1]
            cond_n_gens = out_fields["pops"][2] ? 1 : n_gens
            rewld_pops_local = Array{typ_float}(undef, wld_stats["max"]..., cond_n_gens)
            fill!(rewld_pops_local,NaN)
        else
            rewld_pops_local = NaN
        end
        if out_fields["mutsdel"][1]
            cond_n_gens = out_fields["mutsdel"][2] ? 1 : n_gens
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

            condi = out_fields["gt"][1] && !out_fields["gt"][2]
            gg = condi ? g : 1

            wld_gt_next, wld_fitn_next, wld_pops_next, wld_mutsdel_next, wld_mutsben_next = build_next_gen_inf(
                rewld_gt_local[repeat([:],wlddim)...,max(gg-1,1)], wld_stats, 
                out_fields["fitn"][1], out_fields["pops"][1], out_fields["mutsdel"][1];
                max_migr=max_migr, migr_mode=migr_mode, bottleneck=bottleneck, r_max_migr=r_max_migr, r_coords=r_coords,
                mutratelocus=mutratelocus, weightfitn=weightfitn, condsel=condsel, fixed_mate=fixed_mate, premutate=premutate)

            rewld_gt_local[repeat([:],wlddim)...,gg] = wld_gt_next

            if out_fields["fitn"][1]
                if !out_fields["fitn"][2]
                    rewld_fitn_local[repeat([:],wlddim)...,g] = wld_fitn_next
                elseif g == n_gens
                    rewld_fitn_local[repeat([:],wlddim)...,1] = wld_fitn_next
                end
            end
            if out_fields["pops"][1]
                if !out_fields["pops"][2]
                    rewld_pops_local[repeat([:],wlddim)...,g] = wld_pops_next
                elseif g == n_gens
                    rewld_pops_local[repeat([:],wlddim)...,1] = wld_pops_next
                end
            end
            if out_fields["mutsdel"][1]
                if !out_fields["mutsdel"][2]
                    rewld_mutsdel_local[repeat([:],wlddim+1)...,g] = wld_mutsdel_next
                    rewld_mutsben_local[repeat([:],wlddim+1)...,g] = wld_mutsben_next
                elseif g == n_gens
                    rewld_mutsdel_local[repeat([:],wlddim+1)...,1] = wld_mutsdel_next
                    rewld_mutsben_local[repeat([:],wlddim+1)...,1] = wld_mutsben_next
                end
            end
        end

        return Dict{String, Any}("wld" => rewld_gt_local, "fitn" => rewld_fitn_local, "pops" => rewld_pops_local,
        "mutsdel" => rewld_mutsdel_local, "mutsben" => rewld_mutsben_local)
    end

    if n_re > 1
        if multiproc
            #addprocs(n_re)
            #@everywhere include(@__FILE__)
            dicts_out = pmap(tsk_re, 1:n_re) # uses workers, confirmed
            npcs = workers()
            println("Running $n_re replicates on $npcs processes")
        else
            dicts_out = ThreadsX.map(p->tsk_re(), 1:n_re)
            nt = Threads.nthreads()
            println("Using all $nt threads")
        end
    else
        dicts_out = [tsk_re()]
    end

    wld_stats["max_burnin"] = max_burnin
    wld_stats["max_exp"] = max_exp
    wld_stats["n_gens_burnin"] = n_gens_burnin
    wld_stats["n_gens_exp"] = n_gens_exp
    wld_stats["n_gens"] = n_gens

    res = OrderedDict{String, Any}("stats" => wld_stats)

    for key in keys(out_fields)
        if key=="mutsdel" || key=="mutsben" 
            res[key] = cat([dictus[key] for dictus in dicts_out]...,dims=wlddim+3)
        elseif out_fields[key][1]
            res[key] = cat([dictus[key] for dictus in dicts_out]...,dims=wlddim+2)
        end
    end

    return res
end

function rangeexp_ray_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; x_max_burnin=DEF_X_MAX_BURNIN, x_max_exp=DEF_X_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN, prolif_rate=DEF_PROLIF_RATE,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, mutratelocus=false, weightfitn=true, condsel=false, fixed_mate=false, premutate=false, 
    n_segr_regions=DEF_N_SEGR_REGIONS, regions=fill(sel_coef,n_segr_regions),
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, multiproc=true, wld_gt=NaN, wld_stats=NaN, capacity=DEF_CAPACITY)

    rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; max_burnin=(x_max_burnin,), max_exp=(x_max_exp,), maxi=(x_max_exp,), migr_mode=migr_mode,
        data_to_generate=data_to_generate, name=name, bottleneck=bottleneck, capacity=capacity, prolif_rate=prolif_rate, multiproc=multiproc, weightfitn=weightfitn,
        condsel=condsel, fixed_mate=fixed_mate, premutate=premutate, mutratelocus=mutratelocus,
        mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts,
        n_segr_regions=n_segr_regions, regions=regions, startfill_range=startfill_range, 
        wld_gt=wld_gt, wld_stats=wld_stats)
end

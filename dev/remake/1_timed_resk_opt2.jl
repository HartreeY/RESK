using StatsBase, Distributions, Random, SpecialFunctions, Serialization, Dates, DataStructures, Distributed
Random.seed!(1234)
# Define concrete types to improve type stability
struct WorldStats
    name::String
    max::Tuple
    capacity::Int
    prolif_rate::Float64
    n_segr_regions::Int
    mut_rate::Float64
    migr_rate::Float64
    migr_mode::String
    sel_coef::Float64
    prop_of_del_muts::Float64
    wlddim::Int
    
    # Optional fields with defaults
    max_burnin::Tuple
    max_exp::Tuple
    n_gens_burnin::Int
    n_gens_exp::Int
    n_gens::Int
    startfill::Vector{UnitRange{Int64}}
    n_demes_startfill::Int
end

# Pre-define migration direction constants
const MIGR_DIRS_ORT_1D = [[1], [-1]]
const MIGR_DIRS_ORT_2D = [[-1, 0], [0, -1], [0, 1], [1, 0]]
const MIGR_DIRS_ORT_3D = [[-1, 0, 0], [1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, -1], [0, 0, 1]]

const MIGR_DIRS_DIAG_2D = [[-1, -1], [-1, 1], [1, -1], [1, 1]]
const MIGR_DIRS_HEX_2D = [[-1, 0], [0, -1], [-1, 1], [0, 1], [1, 0], [1, 1]]

# Predefined probabilities instead of dictionary lookups
const MIGR_PROB_ORT_1D = (1.0, 0.0)
const MIGR_PROB_ORT_2D = (1.0, 0.0)
const MIGR_PROB_ALL_2D = (0.5, 0.5)
const MIGR_PROB_DIAG1_2D = (2.0/3.0, 1.0/3.0)

# More specific typing for arrays
const typ_float = Float32
const typ_int = Int32
const typ_gt_inf = Array{Array{Array{typ_float}}}

# Fast fitness calculation function - inlined for performance
@inline function calc_fitness(person::Array{typ_float,1})
    return prod(person)
end

# Optimized mutation function
@inline function mutate_inf!(person::Array{typ_float,1}, mut_rate::Float64, 
                            n_segr_regions::Int, sel_coef::Float64, 
                            prop_of_del_muts::Float64)
    muts_del = 0
    muts_ben = 0
    
    # Get number of mutations based on rate
    get_mutation_random = rand(Poisson(mut_rate))
    
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

# Optimized crossover function
@inline function crossover_inf!(person::Array{typ_float,1}, n_segr_regions::Int)
    @inbounds for i in 1:n_segr_regions
        if rand(Bool)  # More efficient than rand(1:2) == 1
            # Do nothing (keep first allele)
        else
            person[i] = person[i+n_segr_regions]
        end
    end
end

# Optimized mating function with fewer allocations
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

# Migration calculation with reduced branching
function calc_migr_dist!(move::Vector{Int16}, deme::Vector{Int}, stats::WorldStats, 
                        max_migr::Tuple, r_max_migr::Float64=0.0, r_coords::Vector{Int}=[1, 2])
    wlddim = stats.wlddim
    fill!(move, 0)
    
    # Early return if no migration
    if rand() >= stats.migr_rate
        return
    end
    
    # Migration mode handling - simplified with fewer conditionals
    if wlddim == 1
        # 1D has only orthogonal directions
        dir = rand() < 0.5 ? MIGR_DIRS_ORT_1D[1] : MIGR_DIRS_ORT_1D[2]
    elseif wlddim == 2
        # 2D case
        if stats.migr_mode == "hex"
            dir_idx = rand(1:length(MIGR_DIRS_HEX_2D))
            dir = MIGR_DIRS_HEX_2D[dir_idx]
        elseif stats.migr_mode == "all"
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
    
    # Apply radius constraint if needed
    if r_max_migr > 0
        r2 = 0.0
        for i in r_coords
            r2 += (deme[i] - (max_migr[i]-1)/2 + move[i] - 1)^2
        end
        
        if r2 > r_max_migr * r_max_migr
            for i in r_coords
                move[i] = 0
            end
        end
    end
    
    # Apply boundary constraints
    for i in 1:wlddim
        try_move = deme[i] + move[i]
        if try_move > max_migr[i] || try_move < 1
            move[i] = 0  # Just stop at boundary
        end
    end
end

# Optimized offspring calculator
function calc_offspring(wld::typ_gt_inf, stats::WorldStats)
    next_gen_posits = Vector{Vector{Int}}()
    next_gen_pops = Array{Int16}(undef, stats.max...)
    fill!(next_gen_pops, -1)
    
    @inbounds for idx in CartesianIndices(stats.max)
        if length(wld[idx]) > 0
            n_ppl_at_deme = length(wld[idx])
            expected_offspring = n_ppl_at_deme * (stats.prolif_rate / 
                               (1 + (n_ppl_at_deme * (stats.prolif_rate - 1)) / stats.capacity))
            
            next_gen_pops[idx] = rand(Poisson(expected_offspring))
            if next_gen_pops[idx] > 0
                push!(next_gen_posits, [Tuple(idx)...])
            end
        end
    end
    
    return next_gen_posits, next_gen_pops
end

# Main function for building next generation
function build_next_gen_inf!(genno, wld_gt_next::typ_gt_inf, wld_gt::typ_gt_inf, stats::WorldStats, pops_next, fitness_next;
                         max_migr::Tuple=stats.max, r_max_migr::Float64=0.0,
                         weightfitn::Bool=true, condsel::Bool=false, 
                         fixed_mate::Bool=false, premutate::Bool=false, SS::Bool=false,verbose=false)
    
    # Preallocate buffers for migration
    move_buffer = zeros(Int16, stats.wlddim)
    wlddim = stats.wlddim
    
    # Calculate offspring counts
    next_gen_posits, next_gen_pops = calc_offspring(wld_gt, stats)
    if verbose
        println(next_gen_pops)
    end
    
    # Pre-allocate reused arrays for mating
    mate_result = ones(typ_float, stats.n_segr_regions * 2)
    
    # Process each deme
    for deme in next_gen_posits
        inds_at_pos = wld_gt[deme...]
        
        # Calculate fitnesses once
        fitnesses = Vector{typ_float}(undef, length(inds_at_pos))
        @inbounds for i in 1:length(inds_at_pos)
            fitnesses[i] = calc_fitness(inds_at_pos[i])
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
                    print(" ",thing1," ", thing2," ")
                end
                if !condsel || (fitnesses[mom_idx] > thing1*maximum(fitnesses) && 
                               fitnesses[dad_idx] > thing2*maximum(fitnesses))
                    k+=1
                    # Create gametes with recombination
                    gamete_mom = copy(mom)
                    gamete_dad = copy(dad)
                    
                    # Apply recombination
                    crossover_inf!(gamete_mom, stats.n_segr_regions)
                    crossover_inf!(gamete_dad, stats.n_segr_regions)
                    
                    # Create zygote
                    if premutate
                        # Mutations in gametes
                        mutate_inf!(gamete_mom, stats.mut_rate, stats.n_segr_regions, 
                                   stats.sel_coef, stats.prop_of_del_muts)
                        mutate_inf!(gamete_dad, stats.mut_rate, stats.n_segr_regions, 
                                   stats.sel_coef, stats.prop_of_del_muts)
                        # Mating with fixed buffer
                        mate_inf!(mate_result, gamete_mom, gamete_dad, stats.n_segr_regions, fixed_mate)
                    else
                        # Mating first
                        mate_inf!(mate_result, gamete_mom, gamete_dad, stats.n_segr_regions, fixed_mate)
                        
                        # Then mutate
                        mutate_inf!(mate_result, stats.mut_rate, stats.n_segr_regions, 
                                   stats.sel_coef, stats.prop_of_del_muts)
                    end
                    if verbose
                        println(k," ",rand(1))
                    end
                    # Calculate migration
                    calc_migr_dist!(move_buffer, deme, stats, max_migr, r_max_migr)
                    
                    # Apply migration to calculate new position
                    indices = Vector{Int}(undef, wlddim)
                    for i in 1:wlddim
                        indices[i] = deme[i] + move_buffer[i]
                    end
                    
                    # Add offspring to next generation
                    push!(wld_gt_next[indices...], copy(mate_result))
                elseif !SS
                    k+=1
                end
            end
        end
        pops_next[deme...,genno]=length(wld_gt_next[deme...])
        fitness_next[deme...,genno]=mean(calc_fitness.(wld_gt_next[deme...]))
    end
end

# Main simulation function
function rangeexp_ray_inf(n_gens_burnin::Int=100, n_gens_exp::Int=400, n_re::Int=1;
                               x_max_burnin::Int=5, x_max_exp::Int=100, migr_mode::String="ort",
                               prolif_rate::Float64=2.0, capacity::Int=70,
                               mut_rate::Float64=0.08, migr_rate::Float64=0.2, 
                               sel_coef::Float64=0.01, prop_of_del_muts::Float64=0.9,
                               n_segr_regions::Int=20, weightfitn::Bool=false, 
                               condsel::Bool=true, fixed_mate::Bool=true, 
                               premutate::Bool=true, SS::Bool=true, verbose=false)
    
    # Create a typed stats structure
    stats = WorldStats(
        "sim_$(Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"))",
        (x_max_exp,),             # max dimensions
        capacity,
        prolif_rate,
        n_segr_regions,
        mut_rate,
        migr_rate,
        migr_mode,
        sel_coef,
        prop_of_del_muts,
        1,                       # 1D world
        (x_max_burnin,),         # max_burnin
        (x_max_exp,),            # max_exp
        n_gens_burnin,
        n_gens_exp,
        n_gens_burnin + n_gens_exp,
        Vector{UnitRange{Int64}}(undef,0),
        0
    )
    
    # Create the world
    wld_gt = (typ_gt_inf)(undef, stats.max...)
    for k in CartesianIndices(stats.max)
        wld_gt[k] = Array{typ_float,1}[]
    end

    # Fill initial demes
    possible_init_coords = collect(CartesianIndices((x_max_burnin,)))
    init_coords = sample(possible_init_coords, 5, replace=false)  # Start with 5 demes
    
    for coord in init_coords
        for _ in 1:stats.capacity
            push!(wld_gt[coord], ones(typ_float, stats.n_segr_regions * 2))
        end
    end
    
    # Create array for output
    wld_gt_next = (typ_gt_inf)(undef, stats.max...)
    for k in CartesianIndices(stats.max)
        wld_gt_next[k] = Array{typ_float,1}[]
    end
    
    
    final_pops = zeros(typ_float, stats.max...,n_gens_burnin + n_gens_exp)
    final_fitness = zeros(typ_float, stats.max...,n_gens_burnin + n_gens_exp)

    # Run burn-in phase
    for j in 1:n_gens_burnin
        # Clear next generation
        for k in CartesianIndices(stats.max)
            empty!(wld_gt_next[k])
        end
        
        # Build next generation
        build_next_gen_inf!(j, wld_gt_next, wld_gt, stats, final_pops, final_fitness;
                        max_migr=(x_max_burnin,), 
                        weightfitn=weightfitn,
                        condsel=condsel, 
                        fixed_mate=fixed_mate,
                        premutate=premutate,SS=SS,verbose=verbose)
        
        # Swap current and next generation
        wld_gt, wld_gt_next = wld_gt_next, wld_gt
    end
    
    # Run expansion phase
    for j in (n_gens_burnin+1):(n_gens_burnin+n_gens_exp)
        # Clear next generation
        for k in CartesianIndices(stats.max)
            empty!(wld_gt_next[k])
        end
        
        # Build next generation
        build_next_gen_inf!(j, wld_gt_next, wld_gt, stats, final_pops, final_fitness;
                        max_migr=(x_max_exp,), 
                        weightfitn=weightfitn,
                        condsel=condsel, 
                        fixed_mate=fixed_mate,
                        premutate=premutate,SS=SS,verbose=verbose)
        
        # Swap current and next generation
        wld_gt, wld_gt_next = wld_gt_next, wld_gt
    end
    
    return (fitness=final_fitness,pops=final_pops)
    #return (fitness=final_fitness)
end
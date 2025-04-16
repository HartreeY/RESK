# Performance Optimizations:
# 1. Use pre-allocated arrays and in-place operations
# 2. Replace Dict with NamedTuple for stats (immutable, faster)
# 3. Use StaticArrays for small fixed-size arrays
# 4. Minimize array allocations in hot loops
# 5. Use more efficient data structures

using StatsBase, Distributions, Distributed, Random, SpecialFunctions, Serialization, Dates
using StaticArrays, BenchmarkTools

const typ_float = Float32
const typ_int = Int32

# Replace Dict with NamedTuple for better performance
Base.@kwdef struct WorldStats
    name::String
    max::Tuple
    capacity::Int
    prolif_rate::Float64
    regions::Vector{Float64}
    n_segr_regions::Int
    mut_rate::Float64
    migr_rate::Float64
    migr_mode::String
    sel_coef::Float64
    prop_of_del_muts::Float64
    wlddim::Int
end

# Pre-compute migration directions using StaticArrays
const MIGR_DIRS_ORT = [
    [@SVector([1]), @SVector([-1])], # 1D
    [@SVector([-1, 0]), @SVector([0, -1]), @SVector([0, 1]), @SVector([1, 0])], # 2D
    [@SVector([-1, 0, 0]), @SVector([1, 0, 0]), @SVector([0, 1, 0]),
        @SVector([0, -1, 0]), @SVector([0, 0, -1]), @SVector([0, 0, 1])] # 3D
]

# Pre-allocate commonly used arrays
const THREAD_BUFFERS = Dict{Int,Vector{typ_float}}()

# Get or create thread-local buffer
function get_thread_buffer(size)
    tid = Threads.threadid()
    if !haskey(THREAD_BUFFERS, tid)
        THREAD_BUFFERS[tid] = Vector{typ_float}(undef, size)
    end
    return THREAD_BUFFERS[tid]
end

# Optimized mutation function using in-place operations
function mutate_inf!(person::Vector{typ_float}, stats::WorldStats)
    muts_del = muts_ben = 0
    n_muts = rand(Poisson(stats.mut_rate))

    @inbounds for _ in 1:n_muts
        pos = rand(1:stats.n_segr_regions)
        if rand() < stats.prop_of_del_muts
            person[pos] *= (1 - stats.sel_coef)
            muts_del += 1
        else
            person[pos] *= (1 + stats.sel_coef)
            muts_ben += 1
        end
    end
    return muts_del, muts_ben
end

# Optimized crossover using pre-allocated buffer
function crossover_inf!(person::Vector{typ_float}, stats::WorldStats)
    buffer = get_thread_buffer(stats.n_segr_regions)
    @inbounds for i in 1:stats.n_segr_regions
        buffer[i] = rand(Bool) ? person[i] : person[i+stats.n_segr_regions]
    end
    copyto!(person, 1, buffer, 1, stats.n_segr_regions)
end

# Optimized offspring calculation using views
function calc_offspring!(next_gen_posits::Vector{Vector{Int}}, next_gen_pops::Array{Float64},
    wld::Array, stats::WorldStats)
    @inbounds for idx in CartesianIndices(wld)
        if !isempty(wld[idx])
            n_ppl = length(wld[idx])
            expected = n_ppl * (stats.prolif_rate / (1 + (n_ppl * (stats.prolif_rate - 1)) / stats.capacity))
            next_gen_pops[idx] = rand(Poisson(expected))
            if next_gen_pops[idx] > 0
                push!(next_gen_posits, [Tuple(idx)...])
            end
        end
    end
end

# Main simulation function with optimizations
function build_next_gen_inf!(wld_gt::Array, stats::WorldStats,
    next_gen::Array, fitn_out::Array, pops_out::Array)
    next_gen_posits = Vector{Vector{Int}}()
    next_gen_pops = similar(wld_gt, Float64)
    fill!(next_gen_pops, NaN)

    # Calculate offspring counts
    calc_offspring!(next_gen_posits, next_gen_pops, wld_gt, stats)

    # Process each populated location
    Threads.@threads for pos in next_gen_posits
        process_deme!(pos, wld_gt, next_gen, fitn_out, pops_out,
            next_gen_pops, stats)
    end

    return next_gen, fitn_out, pops_out
end

# Helper function to process a single deme
function process_deme!(pos, wld_gt, next_gen, fitn_out, pops_out,
    next_gen_pops, stats)
    inds = wld_gt[pos...]
    if isempty(inds)
        return
    end

    # Pre-calculate fitness once
    fitnesses = [prod(ind) for ind in inds]
    n_offspring = round(Int, next_gen_pops[pos...])

    if !isnothing(fitn_out)
        fitn_out[pos...] = mean(fitnesses)
    end

    # Generate offspring
    offspring = Vector{typ_float}[]
    sizehint!(offspring, n_offspring)

    for _ in 1:n_offspring
        # Select parents using pre-calculated fitness
        parent1 = sample(inds, Weights(fitnesses))
        parent2 = sample(inds, Weights(fitnesses))

        # Create offspring using thread-local buffers
        child = similar(parent1)
        mate_and_mutate!(child, parent1, parent2, stats)
        push!(offspring, child)
    end

    if !isnothing(pops_out)
        pops_out[pos...] = length(offspring)
    end

    next_gen[pos...] = offspring
end

# Optimized mating and mutation
function mate_and_mutate!(child::Vector{typ_float},
    parent1::Vector{typ_float},
    parent2::Vector{typ_float},
    stats::WorldStats)
    # Perform crossing over
    @inbounds for i in 1:stats.n_segr_regions
        child[i] = rand(Bool) ? parent1[i] : parent2[i+stats.n_segr_regions]
    end

    # Mutate in-place
    mutate_inf!(child, stats)
end

# Continued optimizations for remaining functions

# Pre-allocate commonly used arrays and matrices
function create_empty_world_inf(maxi=(DEF_X_MAX, DEF_Y_MAX); min=(1, 1), name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), capacity=DEF_CAPACITY,
    prolif_rate=DEF_PROLIF_RATE, n_segr_regions=DEF_N_SEGR_REGIONS, regions=fill(sel_coef,n_segr_regions),
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, migr_mode=DEF_MIGR_MODE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)
    #stats_params...)
    # Preallocate world array with empty vectors
    wld_gt = Array{Vector{typ_float}}(undef, maxi...)
    @inbounds for i in CartesianIndices(wld_gt)
        wld_gt[i] = Vector{typ_float}()
    end

    # Create WorldStats struct (faster than Dict)
    wld_stats = WorldStats(
        name=name,
        max=maxi,
        wlddim=length(maxi),
        capacity=capacity,
        prolif_rate=prolif_rate,
        regions=regions,
        n_segr_regions=n_segr_regions,
        mut_rate=mut_rate,
        migr_rate=migr_rate,
        migr_mode=migr_mode,
        sel_coef=sel_coef,
        prop_of_del_muts=prop_of_del_muts
        #stats_params...
    )

    return wld_gt, wld_stats
end

# Optimized deme filling with pre-allocation
function fill_random_demes_inf!(wld_gt::Array{Vector{typ_float}},
    stats::WorldStats,
    fill::Vector{UnitRange{Int64}},
    n_demes_to_fill::Int=DEF_N_DEMES_STARTFILL;
    redims::Union{Tuple,Number}=NaN)

    # Create cartesian product once
    possible_coords = vec([CartesianIndex(x...) for x in Iterators.product(fill...)])

    # Sample coordinates efficiently
    init_coords = sample(possible_coords, n_demes_to_fill; replace=false)

    # Pre-allocate individual array for reuse
    individual = ones(typ_float, stats.n_segr_regions * 2)

    # Fill selected demes
    @inbounds for coord in init_coords
        if isnan(redims)
            deme = wld_gt[coord]
        else
            deme = wld_gt[coord, redims...]
        end

        empty!(deme)
        sizehint!(deme, stats.capacity)

        for _ in 1:stats.capacity
            push!(deme, copy(individual))
        end
    end
end

# Optimized range expansion with better memory management
function rangeexp_inf(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; max_burnin=(DEF_X_MAX_BURNIN, DEF_Y_MAX), max_exp=(DEF_X_MAX_EXP, DEF_Y_MAX), maxi=(DEF_X_MAX, DEF_Y_MAX), migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, r_max_burnin=0, r_max_exp=0, r_coords=[1, 2], capacity=DEF_CAPACITY, prolif_rate=DEF_PROLIF_RATE, 
    multiproc=true, weightfitn=true, condsel=false, fixed_mate=false, premutate=false, mutratelocus=false,
    mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, sel_coef=DEF_SEL_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS, n_segr_regions=DEF_N_SEGR_REGIONS, regions=fill(sel_coef,n_segr_regions), startfill_range=NaN, wld_gt=NaN, wld_stats=NaN)

    if isnan(wld_gt)
        is_fill_random_demes = true
        #println("No world provided. Creating a new world.")
        wld_gt, wld_stats = create_empty_world_inf(maxi; name=name, capacity=capacity, prolif_rate=prolif_rate,
            mut_rate=mut_rate, migr_rate=migr_rate, sel_coef=sel_coef, prop_of_del_muts=prop_of_del_muts,
            n_segr_regions=n_segr_regions, regions=regions, migr_mode=migr_mode)
        if !isa(startfill_range, Array) && !any(isnan, max_burnin)
            startfill_range = [1:upper for upper in max_burnin]
        end
    end

    wlddim = wld_stats.wlddim
    wld_stats.max_burnin = max_burnin
    wld_stats.max_exp = max_exp
    wld_stats.n_gens_burnin = n_gens_burnin
    wld_stats.n_gens_exp = n_gens_exp
    n_gens = n_gens_burnin + n_gens_exp
    wld_stats.n_gens = n_gens

    function check_what_output(letter)
        return (occursin(letter, data_to_generate),occursin(letter*"l", data_to_generate))
    end
    out_fields = OrderedDict{String, Any}("gt" => check_what_output("G"), "fitn" => check_what_output("F"), "pops" => check_what_output("P"),
        "mutsdel" => check_what_output("M"), "mutsben" => check_what_output("M"))

    # Create thread-local storage for temporary arrays
    temp_arrays = [create_temp_arrays(params) for _ in 1:Threads.nthreads()]

    # Process each replicate
    Threads.@threads for re in 1:n_re
        process_replicate!(outputs, re, n_gens_burnin, n_gens_exp, params, temp_arrays[Threads.threadid()])
    end

    return finalize_outputs(outputs)
end

# Helper function to process a single replicate
function process_replicate!(outputs, re, n_gens_burnin, n_gens_exp, params, temp_arrays)
    n_gens = n_gens_burnin + n_gens_exp

    # Initialize world for this replicate
    wld_gt, wld_stats = initialize_world(params)

    # Process generations
    @inbounds for g in 1:n_gens
        # Update migration parameters based on burn-in/expansion phase
        if g <= n_gens_burnin
            wld_stats.max_migr = max_burnin
            wld_stats.r_max_migr = r_max_burnin
        else
            wld_stats.max_migr = max_exp
            wld_stats.r_max_migr = r_max_exp
        end

        # Process one generation
        process_generation!(wld_gt, wld_stats, outputs, re, g, migr_params, temp_arrays)

        # Periodic garbage collection to manage memory
        if g % 10 == 0
            GC.gc()
        end
    end
end

# Optimized 1D range expansion
function rangeexp_ray_inf(n_gens_burnin::Int=DEF_N_GENS_BURNIN,
    n_gens_exp::Int=DEF_N_GENS_EXP,
    n_re::Int=1;
    x_max_burnin::Int=DEF_X_MAX_BURNIN,
    x_max_exp::Int=DEF_X_MAX_EXP,
    params...)

    # Convert 1D parameters to general format
    modified_params = Dict{Symbol,Any}(params)
    modified_params[:max_burnin] = (x_max_burnin,)
    modified_params[:max_exp] = (x_max_exp,)
    modified_params[:maxi] = (x_max_exp,)

    return rangeexp_inf(n_gens_burnin, n_gens_exp, n_re; modified_params...)
end

# Helper structs and functions for temporary storage
struct TempArrays
    offspring_counts::Vector{Int}
    fitness_values::Vector{typ_float}
    migration_buffer::Vector{Int}
    individual_buffer::Vector{typ_float}
end

function create_temp_arrays(params)
    max_deme_size = params[:capacity]
    n_segr_regions = params[:n_segr_regions]

    return TempArrays(
        Vector{Int}(undef, max_deme_size),
        Vector{typ_float}(undef, max_deme_size),
        Vector{Int}(undef, params[:wlddim]),
        Vector{typ_float}(undef, n_segr_regions * 2)
    )
end

# Optimized output initialization
function initialize_outputs(params, n_gens, n_re)
    data_types = params[:data_to_generate]
    max_dims = params[:maxi]

    outputs = Dict{String,Any}()

    if contains(data_types, "G")
        outputs["gt"] = Array{Vector{typ_float}}(undef, max_dims..., n_gens, n_re)
    end

    if contains(data_types, "F")
        outputs["fitn"] = Array{typ_float}(undef, max_dims..., n_gens, n_re)
    end

    if contains(data_types, "P")
        outputs["pops"] = Array{typ_float}(undef, max_dims..., n_gens, n_re)
    end

    return outputs
end

# Memory-efficient processing of a single generation
function process_generation!(wld_gt, wld_stats, outputs, re, gen, temp_arrays)
    # Calculate next generation
    next_gen, fitn_out, pops_out = build_next_gen_inf!(wld_gt, wld_stats, temp_arrays)

    # Update outputs
    update_outputs!(outputs, next_gen, re, gen, temp_arrays)

    # Swap generations
    copyto!(wld_gt, next_gen)
end

# Efficient output update
function update_outputs!(outputs, current_gen, re, gen, temp_arrays)
    @inbounds for (key, data) in outputs
        if key == "gt"
            data[.., gen, re] = current_gen
        elseif key == "fitn"
            calculate_fitness!(data[.., gen, re], current_gen, temp_arrays.fitness_values)
        elseif key == "pops"
            calculate_populations!(data[.., gen, re], current_gen)
        end
    end
end
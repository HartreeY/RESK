# ------------------------------------------------
# parameters changed from indiv_based/1d_r_gridrefl2_BH.jl, which runs identically to 1d_r_gridrefl2 (confirmed), but with id_counter and age for comparison with 1d_r_gridrefl2_IB
# Full name: 1d_revised_grid_reflecting_version2_Beverton-Holt_model
# ------------------------------------------------

using StatsBase, Distributions, Random

# Input parameters
# ------------------------------------------------
const BURN_IN_GEN_N = 200
const TOTAL_GEN_N = 500

# Max coordinates of the population bounding space
const X_MAX_BURN_IN = 5
const X_MAX = 500

const X_START = X_MAX_BURN_IN

# Dimensions of the whole space
const X_DIM = X_MAX

# Population parameters
#const INIT_PERSON_N = 30
const DEMES_FULL_AT_START = 5
const K_CAPACITY = 100
const R_PROLIF_RATE = 2
const r_LOG_PROLIF_RATE = log(2)

# Gene parameters
const LOCI_N = 20
const MUT_RATE = 0.05
const M_MIG_RATE = 0.05
const MUT_DELETER_RATE = 0.9
const S_SELECT_COEF = 0.005

# Main program
# ------------------------------------------------
id_counter = 1
x_range = 1:Int(X_START)
init_coords = sample(x_range,DEMES_FULL_AT_START;replace=false)

world = Array{Array{Array{Float32}}}(undef,X_DIM)

for k in 1:X_DIM
    world[k] = Array{Float32,1}[]
end
for coord in init_coords
#=     if !isassigned(world,coord)
        world[coord] = []
    end =#
    for i in 1:K_CAPACITY
        push!(world[coord],vcat(ones(LOCI_N*2),id_counter,0))
        #push!(world[coord],id_counter)
        global id_counter+=1
    end
end

@inbounds function multi_fitn_in_person(person)
    return prod(person[1:(end-2)])
end

@inbounds function max_fitn(persons_at_pos)
    return maximum(multi_fitn_in_person.(persons_at_pos))
end

@inbounds function mean_fitn(persons_at_pos)
    return mean(multi_fitn_in_person.(persons_at_pos))
end

@inbounds function mutate(person)
    get_mutation_random = rand(Poisson(MUT_RATE))
    @fastmath @inbounds for _ in 1:get_mutation_random
        pos_alter = sample(1:LOCI_N)
        if rand() < MUT_DELETER_RATE
            person[pos_alter] *= 1 - S_SELECT_COEF
        else
            person[pos_alter] *= 1 + S_SELECT_COEF
        end
    end
end

@inbounds function recombine(person)
    for i in 1:LOCI_N
        lr = rand(1:2)
        person[i] = lr==1 ? person[i] : person[i+LOCI_N]
    end
end

@inbounds function mate_cond(mom_fit,dad_fit,max_fitness)
    return (mom_fit > rand()*max_fitness) & (dad_fit > rand()*max_fitness)
end

@inbounds function mate(person1,person2)
    global id_counter
    id_counter+=1
    new_loci = vcat(person1[1:LOCI_N],person2[1:LOCI_N], id_counter, 0)
    return new_loci
end

@inbounds @inbounds function build_next_gen(wld,x_max_migrate)
    # Determine the number of offspring for each deme
    next_gen_pops = zeros(Int16,X_DIM)
    birth_chances = zeros(Float32,X_DIM)
    next_gen_posits = []
    fill!(next_gen_pops,-1)
    for x in 1:X_DIM
        if isassigned(world,x) && length(world[x])>0
            n_ppl_at_deme = length(world[x])
            expected_offspring = n_ppl_at_deme * (R_PROLIF_RATE/(1 + (n_ppl_at_deme*(R_PROLIF_RATE-1))/K_CAPACITY))
            next_gen_pops[x] = rand(Poisson(expected_offspring))
            #birth_chances[x] = 1 - expected_offspring/K_CAPACITY/R_PROLIF_RATE
            #println("x: $x, ",birth_chances[x])
            if next_gen_pops[x]>0
                push!(next_gen_posits,x)
            end
        end
    end
    

    # Define the world (habitat)
    wld_next = Array{Array{Array{Float32}}}(undef,X_DIM) #deepcopy(wld)
    for k in 1:X_DIM
        wld_next[k] = Array{Float32,1}[]
    end
    
    all_birth_count = 0

    # Main generation cycle (algorithm)
    mean_fitn_wld = Array{Float32}(undef,X_DIM)
    fill!(mean_fitn_wld,-1)
    pops_wld = zeros(Int32,X_DIM)

    for deme in next_gen_posits
        curr_persons_at_pos = wld[deme]
        mean_fitn_wld[deme] = mean_fitn(curr_persons_at_pos)
        max_fitness =  max_fitn(curr_persons_at_pos)

        next_generation_size = next_gen_pops[deme]
        
        if next_generation_size > 0
            birth_count = 0
            while birth_count < next_generation_size
                mom = curr_persons_at_pos[rand(1:end)]
                dad = curr_persons_at_pos[rand(1:end)]
                mom_fit = multi_fitn_in_person(mom)
                dad_fit = multi_fitn_in_person(dad)

                if mate_cond(mom_fit,dad_fit,max_fitness)
                    gamete_mom = copy(mom) # technically a person, but we'll only use the first half of loci in the mate function
                    gamete_dad = copy(dad) # technically a person, but we'll only use the first half of loci in the mate function
                    recombine(gamete_mom)
                    recombine(gamete_dad)
                    mutate(gamete_mom)
                    mutate(gamete_dad)
                    mate_result = mate(gamete_mom,gamete_dad)

                    wv = [M_MIG_RATE/2,1-M_MIG_RATE,M_MIG_RATE/2]
                    move_x = sample(-1:1,Weights(wv))
                    if deme[1]+move_x > x_max_migrate || deme[1]+move_x < 1
                        #move_x = 0
                        move_x = -move_x
                    end
                    if !isassigned(wld_next,deme[1]+move_x)
                        wld_next[deme[1]+move_x] = []
                    end
                    push!(wld_next[deme[1]+move_x],mate_result)

                    birth_count += 1
                    all_birth_count += 1
                end
            end
            pops_wld[deme] = birth_count
        end
    end
    return wld_next,mean_fitn_wld,pops_wld
end

# Iterate the main cycle and save the output
# ------------------------------------------------

meanf_world_all = Array{Float32}(undef,X_DIM,TOTAL_GEN_N,0)

function rangeexp(n_re=1)
    for _ in 1:n_re
        meanf_world = Array{Float32}(undef,X_DIM,0)

        @inbounds for _ in 1:BURN_IN_GEN_N
            global world,meanf,pops = build_next_gen(world,X_MAX_BURN_IN)
            meanf_world = cat(meanf_world,meanf, dims=2)
        end

        @inbounds for _ in (BURN_IN_GEN_N+1):TOTAL_GEN_N
            global world,meanf,pops = build_next_gen(world,X_MAX)
            meanf_world = cat(meanf_world,meanf, dims=2)
        end

        global meanf_world_all = cat(meanf_world_all, meanf_world, dims=3)
    end
end

# ------------------- Tests ----------------------
# 1. Time for 1 re

#a, b = @timed rangeexp()
#b

# 2. Time for 10 re

a, b = @timed rangeexp(3)
b

# 3. Plots

#= function average_front(data, n_gens, x_max; greaterzero=false, oneside=false, divide=true)
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

function average_ts(ts_arr, n_gens=size(ts_arr))
    return [mean(ts_arr[i,:]) for i in 1:n_gens]
end

function norm_onset_mean(ts, n_gens_burnin::Int)
    normal_array = copy(ts)
    start = n_gens_burnin+1
    normal_array[start:end] /= ts[start]
    return normal_array
end

using Plots
array_nan = replace(meanf_world, -1.0 => NaN)
heatmap(array_nan[:,:,1])

A_fitn_frontav = average_front(reshape(meanf_world, size(meanf_world)..., 1),500,500;oneside=true)
A_fitn_frontav_mean = average_ts(A_fitn_frontav)
A_fitn_frontav_meanN = norm_onset_mean(A_fitn_frontav_mean, BURN_IN_GEN_N+1)
Plots.plot(A_fitn_frontav_meanN[(BURN_IN_GEN_N+2):end],label="Onset mean normalisation",xlabel="Generation") =#
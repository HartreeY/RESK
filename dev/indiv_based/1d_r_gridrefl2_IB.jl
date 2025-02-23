# ------------------------------------------------

# Full name: 1d_revised_grid_reflecting_version2_individual-based_model
# ------------------------------------------------

using StatsBase, Distributions, Random

# Input parameters
# ------------------------------------------------
Random.seed!(1234)
const BURN_IN_GEN_N = 300
const TOTAL_GEN_N = 1000

# Max coordinates of the population bounding space
const X_MAX_BURN_IN = 5
const X_MAX = 100

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
        push!(world[coord],vcat(ones(LOCI_N*2),id_counter,rand(0:2))) # option to randomise
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
    new_loci = vcat(person1[1:LOCI_N],person2[1:LOCI_N], id_counter, 0)
    id_counter+=1
    return new_loci
end

@inbounds function build_next_gen(wld,x_max_migrate)
    # Define the world (habitat)
    wld_next = deepcopy(wld)
    all_birth_count = 0

    # Main generation cycle (algorithm)
    mean_fitn_wld = Array{Float32}(undef,X_DIM)
    fill!(mean_fitn_wld,-1)
    pops_wld = zeros(Int32,X_DIM)
    
    for deme in 1:X_DIM
        if length(wld[deme])>0
            n_ppl_at_deme = length(wld[deme])
            expected_offspring = n_ppl_at_deme * (R_PROLIF_RATE/(1 + (n_ppl_at_deme*(R_PROLIF_RATE-1))/K_CAPACITY))
            birth_chance = 0.5 - expected_offspring/K_CAPACITY/R_PROLIF_RATE
            curr_persons_at_pos = wld[deme]
            mean_fitn_wld[deme] = mean_fitn(curr_persons_at_pos)
            max_fitness =  max_fitn(curr_persons_at_pos)
            
            birth_count = 0
            del_guys = []
            for (ind_i, ind) in pairs(curr_persons_at_pos)


                # Calculate future migration beforehand
                wv = [M_MIG_RATE/2,1-M_MIG_RATE,M_MIG_RATE/2]
                move_x = sample(-1:1,Weights(wv))
                if deme[1]+move_x > x_max_migrate || deme[1]+move_x < 1
                    #move_x = 0
                    move_x = 0
                end
                if move_x!=0
                    push!(wld_next[deme[1]+move_x],ind)
                    push!(del_guys, ind_i)
                end

                # Giving birth
                if ind[end] > 0 && rand() < birth_chance
                    mom = curr_persons_at_pos[rand(1:end)]
                    dad = curr_persons_at_pos[rand(1:end)]
                    mom_fit = multi_fitn_in_person(mom)
                    dad_fit = multi_fitn_in_person(dad)
                    mate_cond_res = mate_cond(mom_fit,dad_fit,max_fitness)
                    while !mate_cond_res
                        mom = curr_persons_at_pos[rand(1:end)]
                        dad = curr_persons_at_pos[rand(1:end)]
                        mom_fit = multi_fitn_in_person(mom)
                        dad_fit = multi_fitn_in_person(dad)
                        mate_cond_res = mate_cond(mom_fit,dad_fit,max_fitness)
                    end
                    
                    gamete_mom = copy(mom) # technically a person, but we'll only use the first half of loci in the mate function
                    gamete_dad = copy(dad) # technically a person, but we'll only use the first half of loci in the mate function
                    recombine(gamete_mom)
                    recombine(gamete_dad)
                    mutate(gamete_mom)
                    mutate(gamete_dad)
                    mate_result = mate(gamete_mom,gamete_dad)

                    push!(wld_next[deme[1]],mate_result)

                    birth_count += 1
                    all_birth_count += 1
                end

                # Death
                ind[end] += 1
                if ind[end]>=3
                    if !(ind_i in del_guys)
                        push!(del_guys, ind_i)
                    end
                end
            end
            #births_wld[deme] = birth_count
            deleteat!(wld_next[deme],del_guys)
            pops_wld[deme] = length(wld_next[deme])
        end
    end
    return wld_next,mean_fitn_wld,pops_wld
end

# Iterate the main cycle and save the output
# ------------------------------------------------

meanf_world = Array{Float32}(undef,X_DIM,0)
pops_world = Array{Float32}(undef,X_DIM,0)

@inbounds for _ in 1:BURN_IN_GEN_N
    global world,meanf,pops = build_next_gen(world,X_MAX_BURN_IN)
    global meanf_world = cat(meanf_world,meanf, dims=2)
    global pops_world = cat(pops_world,pops, dims=2)
end

@inbounds for _ in (BURN_IN_GEN_N+1):TOTAL_GEN_N
    global world,meanf,pops = build_next_gen(world,X_MAX)
    global meanf_world = cat(meanf_world,meanf, dims=2)
    global pops_world = cat(pops_world,pops, dims=2)
end

println("Alive indivs: ",length(collect(Iterators.flatten(world))))
println("Max indiv ID: ",maximum([maximum([k[end-1] for k in i]) for i in filter(!isempty,world)]))

# For the use on HPC
# ---------------------------------
#using Serialization
#procid = myid()-1
#serialize("output/1d/r_gridrefl2_$procid-world.dat",world)
#serialize("output/1d/r_gridrefl2_$procid-pops.dat",pops_world)
#serialize("output/1d/r_gridrefl2_$procid-meanf.dat",meanf_world)
using Plots
heatmap(meanf_world',clim=(0.9,1.1))
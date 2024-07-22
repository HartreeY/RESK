# ------------------------------------------------
# Script for one process of the HPC
# Full name: 2d_revised_axial_grid_reflecting
# ("2d strip" in the paper)
# ------------------------------------------------

using StatsBase
using Distributions

# Input parameters
# ------------------------------------------------

const BURN_IN_GEN_N = 1000
const TOTAL_GEN_N = 2000

# Max coordinates of the population bounding space
# (population = disk)
const X_MAX_BURN_IN = 5
const X_MAX = 250
const Y_MAX = 10

 # (population reaches space bounds)
const X_START = X_MAX_BURN_IN
const Y_START = Y_MAX
# (population starts out twice smaller)
#const X_START = trunc(X_MAX/2)
#const Y_START = trunc(Y_MAX/2)

# Dimensions of the whole space
const X_DIM = X_MAX
const Y_DIM = Y_MAX

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
const possible_dirs = [[-1,0],[0,-1],[0,1],[1,0]] #[[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]]

# Main program
# ------------------------------------------------

x_range = 1:Int(X_START)
y_range = 1:Int(Y_START)
possible_init_coords = [collect(x) for x in Iterators.product(x_range, y_range)]
init_coords = sample(possible_init_coords,DEMES_FULL_AT_START;replace=false)

world = Array{Array{Array{Float32}}}(undef,X_DIM,Y_DIM)

for coord in init_coords
    if !isassigned(world,coord...)
        world[coord...] = []
    end
    for i in 1:K_CAPACITY
        push!(world[coord...],ones(LOCI_N*2))
    end
end

@inbounds function multi_fitn_in_person(person)
    return prod(person)
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
    new_loci = vcat(person1[1:LOCI_N],person2[1:LOCI_N])
    return new_loci
end

@inbounds function build_next_gen(wld,x_max_migrate)
    # Determine the number of offspring for each deme
    next_gen_pops = zeros(Int16,X_DIM,Y_DIM)
    next_gen_posits = []
    fill!(next_gen_pops,-1)
    for x in 1:X_DIM,y in 1:Y_DIM
        if isassigned(world,x,y) && length(world[x,y])>0
            n_ppl_at_deme = length(world[x,y])
            expected_offspring = n_ppl_at_deme * (R_PROLIF_RATE/(1 + (n_ppl_at_deme*(R_PROLIF_RATE-1))/K_CAPACITY))
            next_gen_pops[x,y] =  rand(Poisson(expected_offspring))
            if next_gen_pops[x,y]>0
                push!(next_gen_posits,[x,y])
            end
        end
    end
    

    # Define the world (habitat)
    #wld_next = Array{Person}(undef,sum(next_gen_pops))
    wld_next = Array{Array{Array{Float32}}}(undef,X_DIM,Y_DIM)
    
    all_birth_count = 0

    # Main generation cycle (algorithm)
    mean_fitn_wld = Array{Float32}(undef,X_DIM,Y_DIM)
    fill!(mean_fitn_wld,-1)

    for deme in next_gen_posits
        curr_persons_at_pos = wld[deme...]
        mean_fitn_wld[deme...] = mean_fitn(curr_persons_at_pos)
        max_fitness =  max_fitn(curr_persons_at_pos)

        next_generation_size = next_gen_pops[deme...]
        
        if next_generation_size > 0
            birth_count = 0
            while birth_count < next_generation_size
                mom = curr_persons_at_pos[rand(1:end)]
                dad = curr_persons_at_pos[rand(1:end)]
                mom_fit = multi_fitn_in_person(mom)
                dad_fit = multi_fitn_in_person(dad)




                #print(mate_cond(relative_extract_xx_ind,relative_extract_xy_ind,max_fitness))
                if mate_cond(mom_fit,dad_fit,max_fitness)
                    gamete_mom = copy(mom) # technically a person, but we'll only use the first half of loci in the mate function
                    gamete_dad = copy(dad) # technically a person, but we'll only use the first half of loci in the mate function
                    recombine(gamete_mom)
                    recombine(gamete_dad)
                    mutate(gamete_mom)
                    mutate(gamete_dad)
                    mate_result = mate(gamete_mom,gamete_dad)

                    move_x = 0
                    move_y = 0
                    if rand()<M_MIG_RATE
                        dir = sample(possible_dirs)
                        move_x = dir[1]
                        move_y = dir[2]
                        if deme[1]+move_x > x_max_migrate || deme[1]+move_x < 1
                            #move_x = 0
                            move_x = -move_x
                        end
                        if deme[2]+move_y > Y_MAX || deme[2]+move_y < 1
                            #move_y = 0
                            move_y = -move_y
                        end
                    end

                    if !isassigned(wld_next,deme[1]+move_x,deme[2]+move_y)
                        wld_next[deme[1]+move_x,deme[2]+move_y] = []
                    end
                    push!(wld_next[deme[1]+move_x,deme[2]+move_y],mate_result)

                    birth_count += 1
                    all_birth_count += 1
                end
            end
        end
    end
    return wld_next,mean_fitn_wld
end

# Iterate the main cycle and save the output
# ------------------------------------------------

meanf_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
#pops_world = Array{Float32}(undef,X_DIM,Y_DIM,0)

@inbounds for _ in 1:BURN_IN_GEN_N
    global world,meanf = build_next_gen(world,X_MAX_BURN_IN)
    global meanf_world = cat(meanf_world,meanf, dims=3)
    #global pops_world = cat(pops_world,pops, dims=3)
end

@inbounds for _ in (BURN_IN_GEN_N+1):TOTAL_GEN_N
    global world,meanf = build_next_gen(world,X_MAX)
    global meanf_world = cat(meanf_world,meanf, dims=3)
    #global pops_world = cat(pops_world,pops, dims=3)
end

# For the use on HPC
# ---------------------------------
using Serialization
procid = myid()-1
#serialize("data/2d_axial/lrs_gridrefl_$procid-world.dat",world)
#serialize("data/2d_axial/lrs_gridrefl_$procid-pop.dat",pops_world)
serialize("data/2d_axial/lrs_gridrefl_$procid-meanf.dat",meanf_world)
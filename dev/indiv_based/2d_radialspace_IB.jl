# ------------------------------------------------

# Full name: 2d_space_individual-based_model
# ------------------------------------------------

using StatsBase, Distributions, Random

# Input parameters
# ------------------------------------------------
Random.seed!(1234)
const BURN_IN_GEN_N = 1
const TOTAL_GEN_N = 500

# Max coordinates of the population bounding space
# (population = disk)
const R_MAX_BURN_IN = 0
const R_MAX = 10

const X_MAX = R_MAX
const Y_MAX = R_MAX
const X_START = R_MAX_BURN_IN
const Y_START = R_MAX_BURN_IN

# Dimensions of the whole space
const X_DIM = 2*X_MAX+1
const Y_DIM = 2*Y_MAX+1

# Population parameters
const INDIVS_AT_START = 3 # monoecious
const DEMES_FULL_AT_START = 1
const K_CAPACITY = 20
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
id_counter = 1
x_range = (X_MAX+1-X_START):(X_MAX+1+X_START)
y_range = (Y_MAX+1-Y_START):(Y_MAX+1+Y_START)
possible_init_coords = [collect(x) for x in Iterators.product(x_range, y_range)]
init_coords = sample(possible_init_coords,DEMES_FULL_AT_START;replace=false)

world = Array{Array{Array{Float32}}}(undef,X_DIM,Y_DIM)

for k in 1:X_DIM, j in 1:Y_DIM
    world[k,j] = Array{Float32,1}[]
end
for coord in init_coords
#=     if !isassigned(world,coord...)
        world[coord...] = []
    end =#
    for i in 1:INDIVS_AT_START
        push!(world[coord...],vcat(ones(LOCI_N*2),id_counter,0))
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

@inbounds function build_next_gen(wld,r_max_migrate)
    # Determine the number of offspring for each deme
    next_gen_pops = zeros(Int16,X_DIM,Y_DIM)
    birth_chances = zeros(Float32,X_DIM,Y_DIM)
    next_gen_posits = []
    fill!(next_gen_pops,-1)
    for x in 1:X_DIM,y in 1:Y_DIM
        if isassigned(world,x,y) && length(world[x,y])>0
            n_ppl_at_deme = length(world[x,y])
            expected_offspring = n_ppl_at_deme * (R_PROLIF_RATE/(1 + (n_ppl_at_deme*(R_PROLIF_RATE-1))/K_CAPACITY))
            birth_chances[x,y] = maximum([0,0.5-expected_offspring/K_CAPACITY/R_PROLIF_RATE])
            #println("x: $x, ",birth_chances[x])
            if next_gen_pops[x,y]>0
                push!(next_gen_posits,[x,y])
            end
        end
    end
    

    # Define the world (habitat)
    #wld_next = Array{Person}(undef,sum(next_gen_pops))
    wld_next = deepcopy(wld)
    
    all_birth_count = 0

    # Main generation cycle (algorithm)
    mean_fitn_wld = Array{Float32}(undef,X_DIM,Y_DIM)
    fill!(mean_fitn_wld,-1)
    pops_wld = zeros(Int32,X_DIM,Y_DIM)

    for deme in next_gen_posits
        curr_persons_at_pos = wld[deme...]
        mean_fitn_wld[deme...] = mean_fitn(curr_persons_at_pos)
        max_fitness =  max_fitn(curr_persons_at_pos)

        birth_count = 0
        del_guys = []
        for (ind_i, ind) in pairs(curr_persons_at_pos)

            # Calculate future migration beforehand
            res_x = deme[1]
            res_y = deme[2]
            if rand()<M_MIG_RATE
                dir = sample(possible_dirs)
                move_x = dir[1]
                move_y = dir[2]
                x2 = (deme[1]-X_MAX-1+move_x)*(deme[1]-X_MAX-1+move_x)
                y2 = (deme[2]-Y_MAX-1+move_y)*(deme[2]-Y_MAX-1+move_y)
                r2 = x2+y2
                if r2 > r_max_migrate*r_max_migrate
                    res_x = X_MAX+1-trunc(Int16,(r_max_migrate*r_max_migrate/r2) * move_x)
                    res_y = Y_MAX+1-trunc(Int16,(r_max_migrate*r_max_migrate/r2) * move_y)
                else
                    res_x += move_x
                    res_y += move_y
                end

                push!(wld_next[res_x,res_y],ind)
                push!(del_guys, ind_i)
            end

            # Giving birth
            if ind[end] > 0 && rand() < birth_chances[deme...]
                dad = curr_persons_at_pos[rand(1:end)]
                mom_fit = multi_fitn_in_person(ind)
                dad_fit = multi_fitn_in_person(dad)
                mate_cond_res = mate_cond(mom_fit,dad_fit,max_fitness)
                while !mate_cond_res
                    mom = curr_persons_at_pos[rand(1:end)]
                    dad = curr_persons_at_pos[rand(1:end)]
                    mom_fit = multi_fitn_in_person(mom)
                    dad_fit = multi_fitn_in_person(dad)
                    mate_cond_res = mate_cond(mom_fit,dad_fit,max_fitness)
                end

                gamete_mom = copy(ind) # copy and only use one half in recombine() (haplotype)
                gamete_dad = copy(dad) # copy and only use one half in recombine() (haplotype)
                recombine(gamete_mom)
                recombine(gamete_dad)
                mutate(gamete_mom)
                mutate(gamete_dad)
                mate_result = mate(gamete_mom,gamete_dad)

                push!(wld_next[res_x,res_y],mate_result)

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
        pops_wld[deme...] = birth_count
        deleteat!(wld_next[deme...],del_guys)

    end
    return wld_next,mean_fitn_wld
end

# Iterate the main cycle and save the output
# ------------------------------------------------

meanf_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
#pops_world = Array{Float32}(undef,X_DIM,Y_DIM,0)

@inbounds for _ in 1:BURN_IN_GEN_N
    global world,meanf = build_next_gen(world,R_MAX_BURN_IN)
    global meanf_world = cat(meanf_world,meanf, dims=3)
    #global pops_world = cat(pops_world,pops, dims=3)
end

@inbounds for _ in (BURN_IN_GEN_N+1):TOTAL_GEN_N
    global world,meanf = build_next_gen(world,R_MAX)
    global meanf_world = cat(meanf_world,meanf, dims=3)
    #global pops_world = cat(pops_world,pops, dims=3)
end

println("Alive indivs: ",length(collect(Iterators.flatten(world))))
println("Max indiv ID: ",maximum([maximum([k[end-1] for k in i]) for i in filter(!isempty,world)]))

# For the use on HPC
# ---------------------------------
#using Serialization
#procid = myid()-1
#serialize("output/2d_radial/rrs_gridrefl_$procid-world.dat",world)
#serialize("output/2d_radial/rrs_gridrefl_$procid-pop.dat",pops_world)
#serialize("output/2d_radial/rrs_gridrefl_lat_$procid-meanf.dat",meanf_world)
using Plots
#heatmap(meanf_world[:,:,end],clim=(0.9,1.0))

slow_down = 1
gen_start = 1
gen_end = TOTAL_GEN_N

@gif for i=gen_start:(gen_end*slow_down-1)
    gen_no = trunc(Int,i/slow_down)+1
    heatmap(meanf_world[:,:,gen_no],aspect_ratio=1,yticks=false,clims=(0.9,1.0),xlabel="gen=$gen_no")
end
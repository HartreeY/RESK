# ------------------------------------------------

# Full name: 2d_space_individual-based_model
# ------------------------------------------------

using StatsBase, Distributions, Random

# Input parameters
# ------------------------------------------------
Random.seed!(1234)
const BURN_IN_GEN_N = 1
const TOTAL_GEN_N = 20

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
const LOCI_N = 1000
const MUT_DELETER_RATE = 0.85
const SEL_LOCI_N = trunc(Int,LOCI_N * MUT_DELETER_RATE)
const MUT_RATE = 0.05
const M_MIG_RATE = 0.05
const S_SELECT_COEF = 0.005
const H_DOMIN_COEF = 0.5
const possible_dirs = [[-1,0],[0,-1],[0,1],[1,0]] #[[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]]

# Main program
# ------------------------------------------------
id_counter = 1
x_range = (X_MAX+1-X_START):(X_MAX+1+X_START)
y_range = (Y_MAX+1-Y_START):(Y_MAX+1+Y_START)
possible_init_coords = [collect(x) for x in Iterators.product(x_range, y_range)]
init_coords = sample(possible_init_coords,DEMES_FULL_AT_START;replace=false)

world_ms1 = Array{Array{Array{Int}}}(undef,X_DIM,Y_DIM) # array of all first homologous chromosomes (monosomes="ms") in space
world_ms2 = Array{Array{Array{Int}}}(undef,X_DIM,Y_DIM)

for coord in init_coords
    if !isassigned(world_ms1,coord...)
        world_ms1[coord...] = []
        world_ms2[coord...] = []
    end
    for i in 1:INDIVS_AT_START
        push!(world_ms1[coord...],vcat(rand(0:1,LOCI_N),id_counter,0))
        push!(world_ms2[coord...],vcat(rand(0:1,LOCI_N),id_counter,0))
    end
end
selected_loci = randperm(LOCI_N)[1:SEL_LOCI_N]

@inbounds function muts_by_sel_nonsel(monosomes1,monosomes2)
    len = length(monosomes1)
    muts1s = 0
    muts2s = 0
    muts3s = 0
    muts1ns = 0
    muts2ns = 0
    muts3ns = 0
    fits = []
    
    for i in 1:len
        muts_AA_sel = 0
        muts_Aa_sel = 0
        muts_AA_nonsel = 0
        muts_Aa_nonsel = 0
        new_fitness = 1.0

        for j in 1:LOCI_N
            if monosomes1[i][j]==1 && monosomes2[i][j]==1
                if j in selected_loci
                    muts_AA_sel += 1
                    new_fitness *= 1 - S_SELECT_COEF
                else
                    muts_AA_nonsel += 1
                    new_fitness *= 1 + S_SELECT_COEF
                end

            elseif monosomes1[i][j]==1 || monosomes2[i][j]==1
                if j in selected_loci
                    muts_Aa_sel += 1
                    new_fitness *= 1 - H_DOMIN_COEF * S_SELECT_COEF
                else
                    muts_Aa_nonsel += 1
                    new_fitness *= 1 + H_DOMIN_COEF * S_SELECT_COEF
                end
            end
        end

        push!(fits,new_fitness)
        muts1s += muts_AA_sel
        muts2s += muts_Aa_sel
        muts1ns += muts_AA_nonsel
        muts2ns += muts_Aa_nonsel
        muts3s = SEL_LOCI_N - muts_AA_sel - muts_Aa_sel
        muts3ns = LOCI_N-SEL_LOCI_N - muts_AA_nonsel - muts_Aa_nonsel
    end

    muts1s /= len
    muts2s /= len
    muts3s /= len
    muts1ns /= len
    muts2ns /= len
    muts3ns /= len
    return muts1s,muts2s,muts3s,muts1ns,muts2ns,muts3ns,fits
end

@inbounds function mutate(monosome1,monosome2)
    @fastmath @inbounds for h in 1:LOCI_N

        if rand()<MUT_RATE
            monosome1[h] = true
        end
        if rand()<MUT_RATE
            monosome2[h] = true
        end
    end
end

@inbounds function crossover(monosome1,monosome2)
    for j in 1:LOCI_N
        lr = rand(1:2)
        monosome1[j] = lr==1 ? monosome1[j] : monosome2[j]
    end
end

@inbounds function mate_cond(mom_fit,dad_fit,max_fitness)
    return (mom_fit > rand()*max_fitness) & (dad_fit > rand()*max_fitness)
end


@inbounds function build_next_gen(wld_ms1,wld_ms2,r_max_migrate)
    # Determine the number of offspring for each deme
    next_gen_pops = zeros(Int16,X_DIM,Y_DIM)
    next_gen_posits = []
    fill!(next_gen_pops,-1)
    for x in 1:X_DIM,y in 1:Y_DIM
        if isassigned(wld_ms1,x,y) && length(wld_ms1[x,y])>0
            n_ppl_at_deme = length(wld_ms1[x,y])
            expected_offspring = n_ppl_at_deme * (R_PROLIF_RATE/(1 + (n_ppl_at_deme*(R_PROLIF_RATE-1))/K_CAPACITY))
            next_gen_pops[x,y] =  rand(Poisson(expected_offspring))
            #birth_chances[x] = 1 - expected_offspring/K_CAPACITY/R_PROLIF_RATE
            #println("x: $x, ",birth_chances[x])
            if next_gen_pops[x,y]>0
                push!(next_gen_posits,[x,y])
            end
        end
    end
    

    # Define the world (habitat)
    #wld_next = Array{Person}(undef,sum(next_gen_pops))
    wld_ms1_next = Array{Array{Array{Int}}}(undef,X_DIM,Y_DIM)
    wld_ms2_next = Array{Array{Array{Int}}}(undef,X_DIM,Y_DIM)
    for k in 1:X_DIM, j in 1:Y_DIM
        wld_ms1_next[k,j] = Array{Int,1}[]
        wld_ms2_next[k,j] = Array{Int,1}[]
    end
    
    all_birth_count = 0

    # Main generation cycle (algorithm)
    mean_fitn_wld = Array{Float32}(undef,X_DIM,Y_DIM)
    fill!(mean_fitn_wld,-1)
    pops_wld = zeros(Int32,X_DIM,Y_DIM)
    muts_AAsel_wld = zeros(Float32,X_DIM,Y_DIM)
    muts_Aasel_wld = zeros(Float32,X_DIM,Y_DIM)
    muts_aasel_wld = zeros(Float32,X_DIM,Y_DIM)
    muts_AAnonsel_wld = zeros(Float32,X_DIM,Y_DIM)
    muts_Aanonsel_wld = zeros(Float32,X_DIM,Y_DIM)
    muts_aanonsel_wld = zeros(Float32,X_DIM,Y_DIM)

    for deme in next_gen_posits
        monosomes1_at_pos = wld_ms1[deme...]
        monosomes2_at_pos = wld_ms2[deme...]
        fitns = []
        muts_AAsel_wld[deme...],muts_Aasel_wld[deme...],muts_aasel_wld[deme...],muts_AAnonsel_wld[deme...],muts_Aanonsel_wld[deme...],muts_aanonsel_wld[deme...],fitns = muts_by_sel_nonsel(monosomes1_at_pos,monosomes2_at_pos)
        mean_fitn_wld[deme...] = mean(fitns)
        max_fitness = maximum(fitns)
        sum_fitn = sum(fitns)
        fitns /= sum_fitn

        next_generation_size = next_gen_pops[deme...]
        
        if next_generation_size > 0
            birth_count = 0
            while birth_count < next_generation_size
                guyslen = length(monosomes1_at_pos)
                mom_id = rand(1:guyslen)
                dad_id = rand(1:guyslen)
                mom_fit = fitns[mom_id]
                dad_fit = fitns[dad_id]
                #println(mom_fit,"  ",dad_fit)

                #print(mate_cond(relative_extract_xx_ind,relative_extract_xy_ind,max_fitness))
                if mate_cond(mom_fit,dad_fit,max_fitness)
                    mom1 = copy(monosomes1_at_pos[mom_id])
                    mom2 = copy(monosomes2_at_pos[mom_id])
                    dad1 = copy(monosomes1_at_pos[dad_id])
                    dad2 = copy(monosomes2_at_pos[dad_id])
                    crossover(mom1,mom2)
                    crossover(dad1,dad2)
                    mutate(mom1,mom2)
                    mutate(dad1,dad2)

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
                    end

                    
                    global id_counter
                    mom1[end-1] = id_counter
                    dad2[end-1] = id_counter
                    id_counter += 1

                    push!(wld_ms1_next[res_x,res_y],mom1)
                    push!(wld_ms2_next[res_x,res_y],dad2)

                    birth_count += 1
                    all_birth_count += 1
                end
            end
            pops_wld[deme...] = birth_count
        end
    end
    return wld_ms1_next,wld_ms2_next,mean_fitn_wld,muts_Aasel_wld
end

# Iterate the main cycle and save the output
# ------------------------------------------------

meanf_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_AAsel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_Aasel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_aasel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_AAnonsel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_Aanonsel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_aanonsel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
pops_world = Array{Int32}(undef,X_DIM,Y_DIM,0)

@inbounds for _ in 1:BURN_IN_GEN_N
    global world_ms1, world_ms2,meanf,muts2 = build_next_gen(world_ms1,world_ms2,R_MAX_BURN_IN)
    global meanf_world = cat(meanf_world,meanf, dims=3)
    global muts_Aasel_world = cat(muts_Aasel_world, muts2, dims=3)
end

@inbounds for _ in (BURN_IN_GEN_N+1):TOTAL_GEN_N
    global world_ms1, world_ms2,meanf,muts2  = build_next_gen(world_ms1,world_ms2,R_MAX)
    global meanf_world = cat(meanf_world,meanf, dims=3)
    global muts_Aasel_world = cat(muts_Aasel_world, muts2, dims=3)
end

println("Alive indivs: ",length(collect(Iterators.flatten(world_ms1))))
println("Max indiv ID: ",maximum([maximum([k[end-1] for k in i]) for i in filter(!isempty,world_ms1)]))

# For the use on HPC
# ---------------------------------
#using Serialization
#procid = myid()-1
#serialize("output/2d_radial/rrs_gridrefl_$procid-world.dat",world)
#serialize("output/2d_radial/rrs_gridrefl_$procid-pop.dat",pops_world)
#serialize("output/2d_radial/rrs_gridrefl_lat_$procid-meanf.dat",meanf_world)
using Plots, StatsBase
#heatmap(meanf_world[:,:,end],clim=(0.9,1.0))
maxfitn = maximum(meanf_world)
minfitn = minimum(meanf_world)
#meanf_world_norm = rescale(meanf_world, (0,1))

slow_down = 1
gen_start = 1
gen_end = TOTAL_GEN_N

@gif for i=gen_start:(gen_end*slow_down-1)
    gen_no = trunc(Int,i/slow_down)+1
    heatmap(muts_Aasel_world[:,:,gen_no],aspect_ratio=1,yticks=false,xlabel="gen=$gen_no")
end

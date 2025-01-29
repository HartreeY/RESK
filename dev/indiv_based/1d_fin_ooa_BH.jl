# ------------------------------------------------

# Full name: 1d_finite_out-of-Africa_Beverton-Holt_model
# ------------------------------------------------

using StatsBase, Distributions, Random

# Input parameters
# ------------------------------------------------
Random.seed!(1234)
const BURN_IN_GEN_N = 124
const TOTAL_GEN_N = 244

# Max coordinates of the population bounding space
const X_MAX_BURN_IN = 5
const X_MAX = 60

const X_START = X_MAX_BURN_IN

# Dimensions of the whole space
const X_DIM = X_MAX

# Population parameters
#const INIT_PERSON_N = 30
const DEMES_FULL_AT_START = 5
const K_CAPACITY = 35
const R_PROLIF_RATE = 1.8

# Gene parameters
const LOCI_N = 1000
const SEL_LOCI_N = 312
const MUT_RATE = 1.0 # genome-wide #0.05 * LOCI_N
const M_MIG_RATE = 0.145
const MUT_DELETER_RATE = 0.9
const S_SELECT_COEF = 0.01
const H_DOMIN_COEF = 0.5

# Main program
# ------------------------------------------------
id_counter = 1
x_range = 1:Int(X_START)
init_coords = sample(x_range,DEMES_FULL_AT_START;replace=false)

world_ms1 = Array{Array{Array{Int32}}}(undef,X_DIM)
world_ms2 = Array{Array{Array{Int32}}}(undef,X_DIM)

for k in 1:X_DIM
    world_ms1[k] = Array{Int32,1}[]
    world_ms2[k] = Array{Int32,1}[]
end
for coord in init_coords
#=     if !isassigned(world,coord)
        world[coord] = []
    end =#
    for i in 1:K_CAPACITY
        push!(world_ms1[coord],vcat(falses(LOCI_N),id_counter,0))
        push!(world_ms2[coord],vcat(falses(LOCI_N),id_counter,0))
        #push!(world[coord],id_counter)
        global id_counter+=1
    end
end

selected_loci = randperm(LOCI_N)[1:SEL_LOCI_N]

@inbounds function calc_muts_and_fitn_in_deme(deme_ms1, deme_ms2)
    fits = []
    muts1s = 0
    muts2s = 0
    muts3s = 0
    muts1ns = 0
    muts2ns = 0
    muts3ns = 0
    len = length(deme_ms1)

    for i in 1:len
        muts_AA_sel = 0
        muts_Aa_sel = 0
        muts_AA_nonsel = 0
        muts_Aa_nonsel = 0
        new_fitness = 1.0

        for j in 1:LOCI_N
            if deme_ms1[i][j] == true && deme_ms2[i][j] == true
                if j in selected_loci
                    muts_AA_sel += 1
                    new_fitness *= 1 - S_SELECT_COEF
                else
                    muts_AA_nonsel += 1
                end

            elseif deme_ms1[i][j] == true || deme_ms2[i][j] == true
                if j in selected_loci
                    muts_Aa_sel += 1
                    new_fitness *= 1 - H_DOMIN_COEF * S_SELECT_COEF
                else
                    muts_Aa_nonsel += 1
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

@inbounds function multi_fitn_in_person(person)
    return prod(person[1:(end-2)])
end

@inbounds function max_fitn(persons_at_pos)
    return maximum(multi_fitn_in_person.(persons_at_pos))
end

@inbounds function mean_fitn(persons_at_pos)
    return mean(multi_fitn_in_person.(persons_at_pos))
end

@inbounds function mutate(monosome1,monosome2)
    get_mutation_random = rand(Poisson(MUT_RATE))
    @fastmath @inbounds for _ in 1:get_mutation_random
        pos_alter = sample(1:LOCI_N)

        if rand(1:2)==1
            monosome1[pos_alter] = true
        else
            monosome2[pos_alter] = true
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

@inbounds @inbounds function build_next_gen(wld_ms1,wld_ms2,x_max_migrate)
    # Determine the number of offspring for each deme
    next_gen_pops = zeros(Int16,X_DIM)
    birth_chances = zeros(Float32,X_DIM)
    next_gen_posits = []
    fill!(next_gen_pops,-1)
    for x in 1:X_DIM
        if isassigned(world_ms1,x) && length(world_ms1[x])>0
            n_ppl_at_deme = length(world_ms1[x])
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
    wld_ms1_next = Array{Array{Array{Int32}}}(undef,X_DIM) #deepcopy(wld)
    wld_ms2_next = Array{Array{Array{Int32}}}(undef,X_DIM) #deepcopy(wld)
    for k in 1:X_DIM
        wld_ms1_next[k] = Array{Int32,1}[]
        wld_ms2_next[k] = Array{Int32,1}[]
    end
    
    all_birth_count = 0

    # Main generation cycle (algorithm)
    mean_fitn_wld = Array{Float32}(undef,X_DIM)
    fill!(mean_fitn_wld,-1)
    pops_wld = zeros(Int32,X_DIM)
    muts_AAsel_wld = zeros(Float32,X_DIM)
    muts_Aasel_wld = zeros(Float32,X_DIM)
    muts_aasel_wld = zeros(Float32,X_DIM)
    muts_AAnonsel_wld = zeros(Float32,X_DIM)
    muts_Aanonsel_wld = zeros(Float32,X_DIM)
    muts_aanonsel_wld = zeros(Float32,X_DIM)

    for deme in next_gen_posits
        monosomes1_at_pos = wld_ms1[deme]
        monosomes2_at_pos = wld_ms2[deme]
        fitns = []
        muts_AAsel_wld[deme],muts_Aasel_wld[deme],muts_aasel_wld[deme],muts_AAnonsel_wld[deme],muts_Aanonsel_wld[deme],muts_aanonsel_wld[deme],fitns = calc_muts_and_fitn_in_deme(monosomes1_at_pos,monosomes2_at_pos)
        mean_fitn_wld[deme] = mean(fitns)
        max_fitness =  maximum(fitns)
        sum_fitn = sum(fitns)
        fitns /= sum_fitn

        next_generation_size = next_gen_pops[deme]
        
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
                    mom1 = deepcopy(monosomes1_at_pos[mom_id])
                    mom2 = deepcopy(monosomes2_at_pos[mom_id])
                    dad1 = deepcopy(monosomes1_at_pos[dad_id])
                    dad2 = deepcopy(monosomes2_at_pos[dad_id])
                    crossover(mom1,mom2)
                    crossover(dad1,dad2)
                    mutate(mom1,mom2)
                    mutate(dad1,dad2)

                    wv = [M_MIG_RATE/2,1-M_MIG_RATE,M_MIG_RATE/2]
                    move_x = sample(-1:1,Weights(wv))
                    if deme[1]+move_x > x_max_migrate || deme[1]+move_x < 1
                        move_x = 0
                    end
                    if !isassigned(wld_ms1_next,deme[1]+move_x)
                        wld_ms1_next[deme[1]+move_x] = []
                        wld_ms2_next[deme[1]+move_x] = []
                    end

                    global id_counter
                    id_counter+=1
                    mom1[end-1] = id_counter
                    dad2[end-1] = id_counter
                    push!(wld_ms1_next[deme[1]+move_x],mom1)
                    push!(wld_ms2_next[deme[1]+move_x],dad2)


                    birth_count += 1
                    all_birth_count += 1
                end
            end
            pops_wld[deme] = birth_count
        end
    end
    return wld_ms1_next,wld_ms2_next,mean_fitn_wld,muts_Aasel_wld,muts_Aanonsel_wld
end

# Iterate the main cycle and save the output
# ------------------------------------------------

meanf_world = Array{Float32}(undef,X_DIM,0)
pops_world = Array{Float32}(undef,X_DIM,0)
#muts_AAsel_world = Array{Float32}(undef,X_DIM,0)
muts_Aasel_world = Array{Float32}(undef,X_DIM,0)
#muts_aasel_world = Array{Float32}(undef,X_DIM,0)
#muts_AAnonsel_world = Array{Float32}(undef,X_DIM,0)
muts_Aanonsel_world = Array{Float32}(undef,X_DIM,0)
#muts_aanonsel_world = Array{Float32}(undef,X_DIM,0)

@inbounds for _ in 1:BURN_IN_GEN_N
    global world_ms1,world_ms2,meanf,muts2,muts5 = build_next_gen(world_ms1,world_ms2,X_MAX_BURN_IN)
    global meanf_world = cat(meanf_world,meanf, dims=2)
    #global pops_world = cat(pops_world,pops, dims=2)
    global muts_Aasel_world = cat(muts_Aasel_world, muts2, dims=2)
    global muts_Aanonsel_world = cat(muts_Aanonsel_world, muts5, dims=2)
end

@inbounds for _ in (BURN_IN_GEN_N+1):TOTAL_GEN_N
    global world_ms1,world_ms2,meanf,muts2,muts5 = build_next_gen(world_ms1,world_ms2,X_MAX)
    global meanf_world = cat(meanf_world,meanf, dims=2)
    #global pops_world = cat(pops_world,pops, dims=2)
    global muts_Aasel_world = cat(muts_Aasel_world, muts2, dims=2)
    global muts_Aanonsel_world = cat(muts_Aanonsel_world, muts5, dims=2)
end

println("Alive indivs: ",length(collect(Iterators.flatten(world_ms1))))
println("Max indiv ID: ",maximum([maximum([k[end-1] for k in i]) for i in filter(!isempty,world_ms1)]))

# For the use on HPC
# ---------------------------------
#using Serialization
#procid = myid()-1
#serialize("output/1d/r_gridrefl2_$procid-world.dat",world)
#serialize("output/1d/r_gridrefl2_$procid-pops.dat",pops_world)
#serialize("output/1d/r_gridrefl2_$procid-meanf.dat",meanf_world)
using Plots
#= slow_down = 1
gen_start = 1
gen_end = TOTAL_GEN_N

@gif for i=gen_start:(gen_end*slow_down-1)
    gen_no = trunc(Int,i/slow_down)+1
    heatmap(meanf_world[:,gen_no],aspect_ratio=1,yticks=false,clims=(0.9,1.0),xlabel="gen=$gen_no")
end

@gif for i=gen_start:(gen_end*slow_down-1)
    gen_no = trunc(Int,i/slow_down)+1
    heatmap(muts_Aasel_world[:,gen_no],aspect_ratio=1,clims=(0.5,3.0),yticks=false,xlabel="gen=$gen_no")
end =#
heatmap((muts_Aasel_world+muts_Aanonsel_world)')
function av_by_dist(obj,gen)
    return [mean(obj[x,gen]) for x in 1:size(obj)[1]]
end
plot(av_by_dist((muts_Aanonsel_world+muts_Aasel_world),1:TOTAL_GEN_N))
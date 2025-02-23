# ------------------------------------------------
# Full name: 2d_revised_axial_Out_of_Africa_migrationpattern_0
# (same parameters as 2d_ra_gridE_ramabuff1roul)
# ------------------------------------------------

using StatsBase
using Distributions
using Random

# Input parameters

const FOLDERNAME = "2d_axial_real"
const TESTNAME = "ooa_m0"
const BURN_IN_GEN_N = 124
const TOTAL_GEN_N = 244

# Max coordinates of the population bounding space
# (population = disk)
const X_MAX_BURN_IN = 5
const X_MAX = 60
const Y_MAX = 5

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
const R_PROLIF_RATE = 1.8
const r_LOG_PROLIF_RATE = log(R_PROLIF_RATE)

# Gene parameters
const LOCI_N = 1000
const SEL_LOCI_N = 500
const MUT_RATE = 0.7567 # genome-wide
const M_MIG_RATE = 0.27
const S_SELECT_COEF = 0.002
const H_DOMIN_COEF = 0
const MUT_DELETER_RATE = 0.9 # only for theory
const lat_dirs = [[-1,0],[0,-1],[0,1],[1,0]]
const diag_dirs = [[-1,-1],[-1,1],[1,-1],[1,1]]
x_range = 1:Int(X_START)
y_range = 1:Int(Y_START)
possible_init_coords = [collect(x) for x in Iterators.product(x_range, y_range)]
init_coords = sample(possible_init_coords,DEMES_FULL_AT_START;replace=false)

world_ms1 = Array{Array{Array{Bool}}}(undef,X_DIM,Y_DIM) # array of all first homologous chromosomes (monosomes="ms") in space
world_ms2 = Array{Array{Array{Bool}}}(undef,X_DIM,Y_DIM) # array of all second homologous chromosomes (monosomes="ms") in space

for coord in init_coords
    if !isassigned(world_ms1,coord...)
        world_ms1[coord...] = []
        world_ms2[coord...] = []
    end
    for i in 1:K_CAPACITY
        push!(world_ms1[coord...],falses(LOCI_N))
        push!(world_ms2[coord...],falses(LOCI_N))
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
            if monosomes1[i][j]==true && monosomes2[i][j]==true
                if j in selected_loci
                    muts_AA_sel += 1
                    new_fitness *= 1 - S_SELECT_COEF
                else
                    muts_AA_nonsel += 1
                end

            elseif monosomes1[i][j]==true || monosomes2[i][j]==true
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

@inbounds function mate(person1,person2)
    new_loci = vcat(person1[1:LOCI_N],person2[1:LOCI_N])
    return new_loci
end

@inbounds function build_next_gen(wld_ms1,wld_ms2,x_max_migrate;migr_method=1)
    # Function 1
    next_gen_pops = zeros(Int16,X_DIM,Y_DIM)
    next_gen_posits = []
    fill!(next_gen_pops,-1)
    for x in 1:X_DIM,y in 1:Y_DIM
        if isassigned(wld_ms1,x,y) && length(wld_ms1[x,y])>0
            n_ppl_at_deme = length(wld_ms1[x,y])
            expected_offspring = n_ppl_at_deme * (R_PROLIF_RATE/(1 + (n_ppl_at_deme*(R_PROLIF_RATE-1))/K_CAPACITY))
            next_gen_pops[x,y] =  rand(Poisson(expected_offspring))
            if next_gen_pops[x,y]>0
                push!(next_gen_posits,[x,y])
            end
        end
    end
    

    # Function 2
    wld_ms1_next = Array{Array{Array{Bool}}}(undef,X_DIM,Y_DIM)
    wld_ms2_next = Array{Array{Array{Bool}}}(undef,X_DIM,Y_DIM)
    
    all_birth_count = 0

    # Function 3
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
        sum_fitn = sum(fitns)
        fitns /= sum_fitn

        next_generation_size = next_gen_pops[deme...]
        
        if next_generation_size > 0
            birth_count = 0
            for _ in 1:next_generation_size
                mom_ms1 = wsample(monosomes1_at_pos,fitns)
                mom_ms2 = wsample(monosomes2_at_pos,fitns)
                dad_ms1 = wsample(monosomes1_at_pos,fitns)
                dad_ms2 = wsample(monosomes2_at_pos,fitns)
                gamete_mom_ms1 = copy(mom_ms1)
                gamete_dad_ms1 = copy(dad_ms1)
                gamete_mom_ms2 = copy(mom_ms2)
                gamete_dad_ms2 = copy(dad_ms2)

                crossover(gamete_mom_ms1,gamete_mom_ms2)
                crossover(gamete_dad_ms1,gamete_dad_ms2)
                mutate(gamete_mom_ms1,gamete_mom_ms2)
                mutate(gamete_dad_ms1,gamete_dad_ms2)
                migr_res = rand()

                if migr_method == 0
                    p_lat = 1
                    p_diag = 0
                elseif migr_method == 1
                    p_lat = 2/pi
                    p_diag = 1/pi
                elseif migr_method == 2
                    p_lat = 4/3/pi
                    p_diag = 1/3/pi
                elseif migr_method == 3
                    p_lat = 0.4244132
                    p_diag = 0.21221
                elseif migr_method == 4
                    p_lat = 1/2
                    p_diag = 1/2
                elseif migr_method == 5
                    p_lat = 2/3
                    p_diag = 1/3
                end

                move_x = 0
                move_y = 0
                if migr_res < p_lat+p_diag && rand()<M_MIG_RATE
                    if migr_res < p_lat
                        dir = sample(lat_dirs)
                    elseif migr_res < p_lat+p_diag
                        dir = sample(diag_dirs)
                    end
                    move_x = dir[1]
                    move_y = dir[2]
                    if (deme[1]+move_x==X_MAX_BURN_IN+1 && deme[2]+move_y!=ceil(Y_MAX/2)) # bottleneck barrier
                        move_x = 0
                        move_y = 0
                    else
                        if deme[1]+move_x > x_max_migrate || deme[1]+move_x < 1
                            move_x = 0
                            #move_x = -move_x
                        end
                        if deme[2]+move_y > Y_MAX || deme[2]+move_y < 1
                            move_y = 0
                            #move_y = -move_y
                        end
                    end
                end

                if !isassigned(wld_ms1_next,deme[1]+move_x,deme[2]+move_y)
                    wld_ms1_next[deme[1]+move_x,deme[2]+move_y] = []
                    wld_ms2_next[deme[1]+move_x,deme[2]+move_y] = []
                end
                push!(wld_ms1_next[deme[1]+move_x,deme[2]+move_y],gamete_mom_ms1)
                push!(wld_ms2_next[deme[1]+move_x,deme[2]+move_y],gamete_dad_ms2)

                birth_count += 1
                all_birth_count += 1
            end
            pops_wld[deme...] = birth_count
        end
    end
    return pops_wld,wld_ms1_next,wld_ms2_next,mean_fitn_wld,muts_AAsel_wld,muts_Aasel_wld,muts_aasel_wld,muts_AAnonsel_wld,muts_Aanonsel_wld,muts_aanonsel_wld
end

meanf_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_AAsel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_Aasel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_aasel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_AAnonsel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_Aanonsel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
muts_aanonsel_world = Array{Float32}(undef,X_DIM,Y_DIM,0)
pops_world = Array{Int32}(undef,X_DIM,Y_DIM,0)

@inbounds function run_sim(xmax)
    global pops,world_ms1,world_ms2,meanf,muts1,muts2,muts3,muts4,muts5,muts6 = build_next_gen(world_ms1,world_ms2,xmax;migr_method=0)
    global meanf_world = cat(meanf_world, meanf, dims=3)
    global muts_AAsel_world = cat(muts_AAsel_world, muts1, dims=3)
    global muts_Aasel_world = cat(muts_Aasel_world, muts2, dims=3)
    global muts_aasel_world = cat(muts_aasel_world, muts3, dims=3)
    global muts_AAnonsel_world = cat(muts_AAnonsel_world, muts4, dims=3)
    global muts_Aanonsel_world = cat(muts_Aanonsel_world, muts5, dims=3)
    global muts_aanonsel_world = cat(muts_aanonsel_world, muts6, dims=3)
    global pops_world = cat(pops_world, pops, dims=3)
end

@inbounds for heh in 1:BURN_IN_GEN_N
    run_sim(X_MAX_BURN_IN)
end
@inbounds for heh in (BURN_IN_GEN_N+1):TOTAL_GEN_N
    run_sim(X_MAX)
end

res = [meanf_world,
    muts_AAsel_world, 
    muts_Aasel_world, 
    muts_aasel_world,
    muts_AAnonsel_world,
    muts_Aanonsel_world,
    muts_aanonsel_world,
    pops_world]

info = [BURN_IN_GEN_N, TOTAL_GEN_N, X_MAX_BURN_IN, X_MAX, Y_MAX, X_START, Y_START, X_DIM, Y_DIM, DEMES_FULL_AT_START, K_CAPACITY,
    R_PROLIF_RATE, r_LOG_PROLIF_RATE, LOCI_N, SEL_LOCI_N, MUT_RATE, M_MIG_RATE, S_SELECT_COEF, H_DOMIN_COEF, MUT_DELETER_RATE, 0]

function av_by_dist(obj,gen)
    return [mean(obj[x,:,gen]) for x in 1:size(obj)[1]]
end
using Plots
plot(av_by_dist((res[3]+res[6])/LOCI_N,1:TOTAL_GEN_N))
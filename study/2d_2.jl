@everywhere include("../resk.jl")

test, time = @timed rangeexp_disk_inf(1000,500,20; r_max_burnin=3, r_max_exp=100, data_to_generate="F",capacity=100,prolif_rate=2,
    weightfitn=false,condsel=true,premutate=true,fixed_mate=true,
    bottleneck=NaN,mut_rate=0.05,migr_rate=0.05,sel_coef=0.005,startfill_range=[(100+1-3):(100+1+3),(100+1-3):(100+1+3)])
serialize("data/2d_2.re", test)
println(time)
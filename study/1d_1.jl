@everywhere include("../resk.jl")

test = rangeexp_1d_inf(10000,2500,70;x_max_burnin=5,x_max_exp=500,data_to_generate="F",capacity=100,prolif_rate=2,mut_rate=0.05,migr_rate=0.05,sel_coef=0.005,weightfitn=false,condsel=true,fixed_mate=true)
serialize("data/1d/1.re",test)
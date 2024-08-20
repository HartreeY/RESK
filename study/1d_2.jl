@everywhere include("../resk.jl")

test = rangeexp_1d_inf(5000,1250,10;x_max_burnin=5,x_max_exp=500,data_to_generate="F",k_capacity=100,r_prolif_rate=2,mut_rate=0.05,migr_rate=0.05,s_sel_coef=0.005,prop_of_del_muts=0.)
serialize("data/1d/1_2.re",test)
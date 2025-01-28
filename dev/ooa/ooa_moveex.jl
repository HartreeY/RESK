@everywhere include("resk.jl")

test, time = @timed rangeexp_strip(500,2500,50; x_max_burnin=5, x_max_exp=500, y_max=5, data_to_generate="SN",capacity=35, prolif_rate=1.8, n_loci=1000,
    n_sel_loci=312, mut_rate=1,migr_rate=0.145,sel_coef=0.01,bottleneck=NaN)
serialize("data/ooa_moveex.re", test)
println(time)
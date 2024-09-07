@everywhere include("resk.jl")

test, time = @timed rangeexp_strip(120,200,40; x_max_burnin=5, x_max_exp=100, y_max=5, data_to_generate="FPSN",capacity=35, prolif_rate=1.8,
    n_loci=1000, migr_mode="ort",
    n_sel_loci=312,mut_rate=1,migr_rate=0.15,sel_coef=0.01,bottleneck=NaN)
serialize("data/hextest_ort.re", test)
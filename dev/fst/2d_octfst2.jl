@everywhere include("resk.jl")

test, time = @timed rangeexp_strip(50,400,10; x_max_burnin=3, x_max_exp=130, y_max=5,
        data_to_generate="FPCM",capacity=35, prolif_rate=1.8, n_loci=1000,
        n_sel_loci=312, mut_rate=1,migr_rate=0.145,sel_coef=0.01,bottleneck=NaN,migr_mode="diag1/2")
serialize("data/2d_octfst.re", test)
println(time)
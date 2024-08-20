@everywhere include("resk.jl")

test, time = @timed rangeexp_strip_inf(1000, 1000, 20; x_max_burnin=5, x_max_exp=250, y_max=10, 
    data_to_generate="F", capacity=100, prolif_rate=2, mut_rate=0.05, migr_rate=0.05, sel_coef=0.005, weightfitn=false, condsel=true, fixed_mate=true, premutate=true, bottleneck=NaN)
serialize("data/2d_1.re", test)
println(time)
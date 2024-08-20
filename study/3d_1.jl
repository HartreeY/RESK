@everywhere include("resk.jl")

test, time = @timed rangeexp_cylinder_inf(0,300,10; r_max_burnin=10, r_max_exp=10, z_max_burnin=5, z_max_exp=250, data_to_generate="F",capacity=100,
    prolif_rate=2,weightfitn=false,condsel=true,premutate=true,
    fixed_mate=true, bottleneck=NaN, mut_rate=0.05, migr_rate=0.05, sel_coef=0.005, startfill_range=[1:21,1:21,1:5])
serialize("data/3d_1nb.re", test)
println(time)
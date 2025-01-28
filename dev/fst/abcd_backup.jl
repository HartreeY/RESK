include("resk.jl")
re_res,re_time = @timed rangeexp_cylinder(30,60,1;data_to_generate="FPC",z_max_burnin=4, z_max_exp=30,r_max_burnin=3,r_max_exp=3,bottleneck=NaN,cap$n_sel_loci=312, mut_rate=1,migr_rate=0.145,sel_coef=0.01,multiproc=false)
serialize("data/test_fst.re", re_res)
re_time

include("resk.jl")
re_res,re_time = @timed rangeexp_cylinder(30,60,1;data_to_generate="FPC",z_max_burnin=4, z_max_exp=30,r_max_burnin=3,r_max_exp=3,bottleneck=NaN,cap$n_sel_loci=312, mut_rate=1,migr_rate=0.145,sel_coef=0.01,multiproc=false)
serialize("data/test_fst2.re", re_res)
re_time

include("resk.jl")
re_res,re_time = @timed rangeexp_cylinder(30,60,1;data_to_generate="FPC",z_max_burnin=4, z_max_exp=30,r_max_burnin=3,r_max_exp=3,bottleneck=NaN,cap$n_sel_loci=312, mut_rate=1,migr_rate=0.145,sel_coef=0.01,multiproc=false)
serialize("data/test_fst3.re", re_res)
re_time

include("resk.jl")
re_res,re_time = @timed rangeexp_cylinder(30,60,1;data_to_generate="FPC",z_max_burnin=4, z_max_exp=30,r_max_burnin=3,r_max_exp=3,bottleneck=NaN,cap$n_sel_loci=312, mut_rate=1,migr_rate=0.145,sel_coef=0.01,multiproc=false)
serialize("data/test_fst4.re", re_res)
re_time
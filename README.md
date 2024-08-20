# Range Expansions Simulation Kit (RESK)
A set of programs in Julia designed to efficiently simulate range expansions and study their population dynamics.
This set of programs has been used in the study "The evolution of fitness during range expansions in multiple dimensions". You can find the preprint at https://www.biorxiv.org/content/10.1101/2023.12.29.573608v2.

## How to use
### Prerequisites
To begin using this set of tools, you need to have Julia 1.9+ installed, along with the required packages. To install the packages, please run the initialisation script: access the *programs* folder and run the script *init.jl*. This should only take a couple of minutes.

Once you have the required packages, you can use the documented methods of this package on your own, or you can follow one of several examples in the *programs/examples* folder.

### Main files
RESK consists of just a couple of files:
- *resk.jl*: main methods
- *reskplots.jl*: visualisation methods (comes as a separate file due to some machines not supporting graphical output)
- *init.jl*: initialisation script
- *defaults.jl*: easy-to-change list of default constants

### Main use cases
Include the *resk.jl* script and use its methods.

To simulate a range expansion once, use `rangeexp` or methods that start with `rangeexp`:

- 1D: `rangeexp_1d`
- 2D: `rangeexp_disk`,`rangeexp_strip`
- 3D: `rangeexp_cylinder`,`rangeexp_sphere`.

Running this with default options, a world (= habitat) will be created, seeded with individuals, and an expansion will be run on it.

The `rangeexp` functions output a fixed dictionary that includes metadata (**stats**) and expansion data. The types of expansion data within it are determined by the *data_to_generate* argument. It can take on the following values:
- **F** - **fitn** (deme-average fitness)
- **P** - **pops** (deme populations)
- **S** - **AAsel**, **Aasel** and **aasel** (deme-average number of homo- [**AA**, **aa**] and heterozygous [**Aa**] selected loci)
- **N** - **AAneu**, **Aaneu** and **aaneu** (deme-average number of homo- and heterozygous neutral loci)
The above can be combined and should be passed in a string. For example,
```
test = rangeexp_strip(15, 30; data_to_generate="SF", y_max=5)
```
will output
```
Dict{String, Any} with 9 entries:
  "AAsel" => Float32[0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.…
  "Aasel" => Float32[0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.…
  "aasel" => Float32[0.0 0.0 … 0.0 0.0; 25.0 25.0 … 0.0 25.0; … ; 0.0 0.0 … 0.0…
  "pops"  => Float32[]
  "fitn" => Float32[-1.0 -1.0 … -1.0 -1.0; 0.05 0.05 … -1.0 0.05; … ; -1.0 -1.…
  "aaneu" => Float32[]
  "Aaneu" => Float32[]
  "AAneu" => Float32[]
  "stats" => Dict{String, Any}("y_max_burnin"=>10, "x_max"=>100, "migr_mode"=>[…
```

These expansion data can be plotted and worked with. To plot expansion data, include the *reskplots.jl* file and use the unique plotting functions that start with *re_*. For example,
```
re_heatmap_AAsel(test; log_factor=1.02)
```
will output the average number of selected homozygous mutant loci in a deme:
![alt text](https://github.com/HartreeY/RESK/blob/master/img/readme0.gif?raw=true)

Here's an example of a deme-average fitness heatmap of a longer axial simulation in 2D:
```
test = rangeexp_strip(100,1000;data_to_generate="FPSN",y_max=8,migr_mode="diag1/2")
```
![alt text](https://github.com/HartreeY/RESK/blob/master/img/readme1.gif?raw=true)

The above examples use the finite-sites model for individual genotypes. For every `rangeexp` function, there is also an infinite-sites equivalent (e.g. `rangeexp_inf`). The `_inf` functions are computationally faster, but due to the nature of the model, they cannot output hetero- and homozygosities. Instead, they can output **del** and **ben** (deme-average number of deleterious and beneficial mutations). The finite-sites function currently only simulate deleterious mutations. Here are some visualised examples of running `_inf` functions:

`rangeexp_cylinder_inf`:

![alt text](https://github.com/HartreeY/RESK/blob/master/img/readme2.gif?raw=true)

`rangeexp_sphere_inf`:

![alt text](https://github.com/HartreeY/RESK/blob/master/img/readme3.gif?raw=true)

### Distributed and batch simulation

RESK is also to generate multiple replicates in one go, and employs Julia's standard library's powerful parallel processing in those cases by default. The `rangeexp` methods feature the `distributed` option to toggle distributed processing, and the `n_re` to set the number of replicates. The processes are automatically added and removed via the `addprocs` and `rmprocs` methods.

Here is a benchmark, which you can find in *programs/examples*, showing the benefit of distributed simulations:

<img src="https://github.com/HartreeY/RESK/blob/master/img/readme4.png" width="650"/>

## Main methods

### create_empty_world
`(max::Tuple=(DEF_X_MAX, DEF_Y_MAX); name::String=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), k_capacity=DEF_K_CAPACITY, r_prolif_rate=DEF_R_PROLIF_RATE, n_loci=DEF_N_LOCI, n_sel_loci=DEF_N_SEL_LOCI,  mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, migr_mode=DEF_MIGR_MODE, s_sel_coef=DEF_S_SEL_COEF, h_domin_coef=DEF_H_DOMIN_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)`

Creates an empty deme space, i.e. an N-dimensional lattice of demes that can house individuals with a finite-site genetic structure.
N is determined from the dimensions of the `max` tuple (2-dimensional by default).

`max`: world extents \
`name`: world name \
`k_capacity`: capacity of each deme \
`r_prolif_rate`: proliferation rate \
`n_loci`: number of loci \
`n_sel_loci`: number of selected loci \
`mut_rate`: genome-wide mutation rate \
<a name="migr"></a>`migr_mode`: mode of migration. Possible values: \
&nbsp;&nbsp;&nbsp;**ort** - orthogonal directions only \
&nbsp;&nbsp;&nbsp;**all** - orthogonal and diagonal \
&nbsp;&nbsp;&nbsp;**hex** - hexagonal grid \
&nbsp;&nbsp;&nbsp;**diag1/2** - orthogonal and half-weighted diagonal \
&nbsp;&nbsp;&nbsp;**buffon1** - equidistant Buffon-Laplace (see documentation) \
&nbsp;&nbsp;&nbsp;**buffon2** - uniform Buffon-Laplace \
&nbsp;&nbsp;&nbsp;**buffon3** - inv.proportional Buffon-Laplace \
`s_sel_coef`: selection coefficient \
`h_domin_coef`: dominance coefficient (in heterozygous loci, new_fitness *= **1 -** `h_domin_coef` * `s_sel_coef`) \
`prop_of_del_muts`: proportion of deleterious mutations in nature

Output 1: a spatial array of demes that contain individuals' left monosome [Bool] arrays (all empty) \
Output 2: a spatial array of demes that contain individuals' right monosome [Bool] arrays (all empty) \
Output 3: world stats Dict

---

### fill_random_demes
`(pnt_wld_ms1, pnt_wld_ms2, pnt_wld_stats, fill::Vector{UnitRange{Int64}}, n_demes_to_fill=DEF_N_DEMES_STARTFILL)`

Fills random demes within given monosome arrays with finite-sites individuals. Usually used after an empty world is created.

`pnt_wld_ms1`: a spatial array of demes that contain individuals' left monosome [Bool] arrays \
`pnt_wld_ms2`: a spatial array of demes that contain individuals' right monosome [Bool] arrays \
`pnt_wld_stats`: world stats Dict \
`fill`: an array of Int ranges of the coordinates that define the area within which to fill \
`n_demes_to_fill`: number of demes to fill \

---

### rangeexp
`(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; max_burnin=(DEF_X_MAX_BURNIN, DEF_Y_MAX), max_exp=(DEF_X_MAX_EXP, DEF_Y_MAX), max=(DEF_X_MAX, DEF_Y_MAX), migr_mode=DEF_MIGR_MODE,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, r_max_burnin=0, r_max_exp=0, r_coords=[1, 2],
    startfill_range=NaN, distributed=true, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN)`
    
Simulates a range expansion `n_re` times.
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium \
`n_gens_exp`: duration of the expansion \
`n_re`: number of replicates \
`max_burnin`: a tuple of maximum coordinates during burn-in \
`max_exp`: a tuple of maximum coordinates during expansion \
`max`: a tuple of maximum coordinates of space \
`migr_mode`: mode of migration ([possible values](#migr)) \
<a name="dtg"></a>`data_to_generate`: string of letters representing different data to output. Possible values: \
&nbsp;&nbsp;&nbsp;**F** - deme-average fitness (**fitn**) \
&nbsp;&nbsp;&nbsp;**P** - deme populations (**pops**) \
&nbsp;&nbsp;&nbsp;**S** - deme-average number of homo- and heterozygous selected loci (**AAsel**, **Aasel** and **aasel**) \
&nbsp;&nbsp;&nbsp;**N** - deme-average number of homo- and heterozygous neutral loci (**AAneu**, **Aaneu** and **aaneu**) \
`name`: world name \
`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates \
`r_max_burnin`: radius that bounds the burn-in area \
`r_max_exp`: radius that bounds the expansion area \
`r_coords`: a tuple (array) of axes' ordinal numbers that the n-sphere with `r_max_migr` covers. For example: \
&nbsp;&nbsp;&nbsp;**(1,3)** - migration is bound within a disk at x and z axes \
&nbsp;&nbsp;&nbsp;**(1,2,3)** - migration is bound within a sphere at x, y and z axes \
`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start \
`distributed`: if **true**, distribute to threads \
If starting from existing world, also provide: \
`wld_ms1`: a spatial array of demes that contain individuals' left monosome [Bool] arrays \
`wld_ms2`: a spatial array of demes that contain individuals' right monosome [Bool] arrays \
`wld_stats`: world stats Dict

Output: a Dict containing data after the expansion: \
&nbsp;&nbsp;&nbsp;**stats** - statistics array containing world and range expansion information \
&nbsp;&nbsp;&nbsp;**fitn**, **pops**, **AAsel**, **Aasel**, **aasel**, **AAneu**, **Aaneu**, **aaneu** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`

---

### rangeexp_1d
`(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP, n_re=1; x_max_burnin=DEF_X_MAX_BURNIN, x_max_exp=DEF_X_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, distributed=true, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN)`

Simulates a range expansion `n_re` times in 1D, starting from one side of a segment space.

Changes from `rangeexp`:

`x_max_burnin`: the outward x-coordinate bound for migration during burn-in \
`x_max_exp`: the outward x-coordinate bound for migration during the expansion

---

### rangeexp_strip
`(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP; x_max_burnin=DEF_X_MAX_BURNIN, x_max_exp=DEF_X_MAX_EXP, y_max=DEF_Y_MAX, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=("midhole at x=", x_max_burnin * 2),
    wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN, max_burnin=(x_max_burnin, y_max), max_exp=(x_max_exp, y_max), max=(x_max_exp, y_max))`

Simulates a 2D strip range expansion, in which a population expands in the positive x direction.

Changes from `rangeexp`:

`x_max_burnin`: the outward x-coordinate bound for migration during burn-in \
`x_max_exp`: the outward x-coordinate bound for migration during the expansion \
`y_max`: the upper y-coordinate bound (lower bound is always **1** currently)

---

### rangeexp_disk
`(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,
    data_to_generate=DEF_DATA_TO_GENERATE, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN, max=(r_max_exp * 2 + 1, r_max_exp * 2 + 1), 
    max_exp=NaN, max_burnin=NaN, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN)`

Simulates a range expansion, in which a population expands from the radially from the center by default, and that is bound by a disk.

Changes from `rangeexp`:

`r_max_burnin`: radius that bounds the burn-in area \
`r_max_exp`: radius that bounds the expansion area

---

### rangeexp_cylinder
`(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,
    z_max_burnin=DEF_X_MAX_BURNIN, z_max_exp=DEF_X_MAX_EXP, max_burnin=(NaN, NaN, z_max_burnin), max_exp=(NaN, NaN, z_max_exp), max=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, z_max_exp),
    data_to_generate=DEF_DATA_TO_GENERATE, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN)`
    
Simulates a 3D range expansion, in which a population expands in the positive z direction, as well as radially from the z axis, and is bound by a cylinder.

Changes from `rangeexp`:

`r_max_burnin`: XY radius of a cylinder that bounds migration during burn-in \
`r_max_exp`: XY radius of a cylinder that bounds migration during expansion \
`z_max_burnin`: the outward z-coordinate bound for migration during burn-in \
`z_max_exp`: the outward z-coordinate bound for migration during the expansion

---

### rangeexp_sphere
`(n_gens_burnin=DEF_N_GENS_BURNIN, n_gens_exp=DEF_N_GENS_EXP; r_max_burnin=DEF_R_MAX_BURNIN, r_max_exp=DEF_R_MAX_EXP, migr_mode=DEF_MIGR_MODE, startfill_range=NaN,
    max_burnin=NaN, max_exp=NaN, max=(r_max_exp * 2 + 1, r_max_exp * 2 + 1, r_max_exp * 2 + 1),
    data_to_generate=DEF_DATA_TO_GENERATE, wld_ms1=NaN, wld_ms2=NaN, wld_stats=NaN, name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), bottleneck=NaN)`

Simulates a 3D range expansion, in which a population expands radially from the center by default, and that is bound by a sphere.

Changes from `rangeexp`:

`r_max_burnin`: radius of a sphere that bounds migration during burn-in \
`r_max_exp`: radius of a sphere that bounds migration during expansion

---

### <a name="reh1"></a>re_heatmap
`(data::Array, dims::Int, gen_start=1, gen_end=DEF_N_GENS_BURNIN + DEF_N_GENS_EXP, re_index::Int = 1; n_gens_sub=0, slow_factor=1, log_base=-1, clim=:default, kwargs...)`

Shows an animated heatmap of `data` from `gen_start` to `gen_end` in 1D or 2D.

`data`: array with dimensions (space + time) \
`gen_start`: start generation \
`gen_end`: end generation \
`re_index`: which replicate to plot \
`n_gens_sub`: number of generations to subtract; e.g. set this as the number of burn-in gen-s if you wish to display the burn-in gen-s in negative numbers \
`slow_factor`: number of animation frames per generation \
`log_base`: if not **-1**, color shows log values with this as base \
`clim`: color bounds (Plots.jl's `clim` parameter) \
`kwargs...`: any Plots.jl parameters

### re_heatmap
`(re::Dict, dataname::String, gen_start=1, gen_end=re["stats"]["n_gens"]; n_gens_sub=re["stats"]["n_gens_burnin"], slow_factor=1, log_base=-1, clim=:default, kwargs...)`

Shows an animated heatmap of `dataname` in `re` from `gen_start` to `gen_end` in 1D or 2D.

Changes from `re_heatmap`([1](#reh1)):

`re`: range expansion results dictionary \
`dataname`: name of data in `re`

---

### re_heatmap_[dataname]

Shows *[dataname]* data of `re` from `gen_start` to `gen_end`. For example, `re_heatmap_pops`. This is a useful function since it uses graph options that are optimised for each data type.

See `re_heatmap` for more.

## To do
- implement beneficial mutations for finite-sites
- include the possibility of multiple and partial range expansions
- include total mutation count as output too
- add proper error catching
- add elitism, proper epistasis, fitness landscapes (rugged etc.), assortative mating etc.
- implement migration with multiple-deme leaps 

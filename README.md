# Range Expansions Simulation Kit (RESK)
A set of programs in Julia designed to efficiently simulate range expansions and study their population dynamics.
This set of programs has been used in the study "The evolution of fitness during range expansions in multiple dimensions". You can find the preprint at https://www.biorxiv.org/content/10.1101/2023.12.29.573608v2.

## How to use
### Prerequisites
To begin using this set of tools, you need to have Julia 1.9+ installed, along with the required packages. To install the packages, please run the initialisation script: access the *programs* folder and run the script *init.jl*. This should only take a couple of minutes.

Once you have the required packages, you can use the documented methods of this package on your own, or you can follow one of several examples in the *programs/examples* folder.

### Main use cases
Include the *resk.jl* script and use its methods.

To simulate a range expansion once, use `rangeexp` or the methods that start with `rangeexp`:

- 1D: `rangeexp_1d`
- 2D: `rangeexp_disk`,`rangeexp_strip`
- 3D: `rangeexp_cylinder`,`rangeexp_sphere`.

Running this with default options, a world (= habitat) will be created, seeded with individuals, and the expansion will be run on it.

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
  "pops"  => NaN
  "fitn" => Float32[-1.0 -1.0 … -1.0 -1.0; 0.05 0.05 … -1.0 0.05; … ; -1.0 -1.…
  "aaneu" => NaN
  "Aaneu" => NaN
  "AAneu" => NaN
  "stats" => Dict{String, Any}("y_max_burnin"=>10, "x_max"=>100, "migr_mode"=>[…
```

These expansion data can be plotted and worked with. To plot expansion data, use unique plotting functions that start with *re_*. For example,
```
re_heatmap_AAsel(test; log_factor=1.02)
```
will output the average number of selected homozygous mutant loci in a deme:
![alt text](https://github.com/HartreeY/RESK/blob/main/animations/readme0.gif?raw=true)

Here's an example of a deme-average fitness heatmap of a longer axial simulation in 2D:
```
test = rangeexp_strip(100,1000;data_to_generate="FPSN",y_max=8,migr_mode="diag1/2")
```
![alt text](https://github.com/HartreeY/RESK/blob/main/animations/readme1.gif?raw=true)

The above examples use the finite-sites model for individual genomes. For every `rangeexp` function, there is also an infinite-sites equivalent (e.g. `rangeexp_inf`). Practically, `_inf` functions are computationally faster. Some examples of `_inf` functions:

`rangeexp_cylinder_inf`:

![alt text](https://github.com/HartreeY/RESK/blob/main/animations/readme2.gif?raw=true)

`rangeexp_sphere_inf`:

![alt text](https://github.com/HartreeY/RESK/blob/main/animations/readme3.gif?raw=true)

## Main methods

### create_new_world
`(max=(DEF_X_MAX, DEF_Y_MAX); min=(1, 1), name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), k_capacity=DEF_K_CAPACITY, r_prolif_rate=DEF_R_PROLIF_RATE, n_loci=DEF_N_LOCI, n_sel_loci=DEF_N_SEL_LOCI,  mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, migr_mode=DEF_MIGR_MODE, s_sel_coef=DEF_S_SEL_COEF, h_domin_coef=DEF_H_DOMIN_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)`

Creates an empty world (deme space) with finite-sites individual structure. 2-dimensional by default.

`max`: a tuple of maximal space bounds (coordinates) \
`min`: a tuple of minimal space bounds (coordinates). Limited to (1,1) for now \
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
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium \
`n_gens_exp`: duration of the expansion \
`n_re`: number of replicates \
`x_max_burnin`: the outward x-coordinate bound for migration during burn-in \
`x_max_exp`: the outward x-coordinate bound for migration during the expansion \
`migr_mode`: mode of migration ([possible values](#migr)) \
`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start \
`data_to_generate`: string of letters representing different data to output ([possible values](#dtg)) \
`name`: world name \
`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates \
`distributed`: if **true**, distribute to threads \
If starting from existing world, also provide: \
`wld_ms1`: a spatial array of demes that contain individuals' left monosome [Bool] arrays \
`wld_ms2`: a spatial array of demes that contain individuals' right monosome [Bool] arrays \
`wld_stats`: world stats Dict

Output: a Dict containing data after the expansion: \
&nbsp;&nbsp;&nbsp;**stats** - statistics array containing world and range expansion information \
&nbsp;&nbsp;&nbsp;**fitn**, **pops**, **AAsel**, **Aasel**, **aasel**, **AAneu**, **Aaneu**, **aaneu** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`

---

### rangeexp_strip

Simulates a strip range expansion, in which a population expands in the positive x direction (after an optional burn-in phase).
If no world is provided, generates a world and seeds it with `DEF_N_DEMES_STARTFILL` demes filled with individuals.\

`n_gens_burnin`: duration of the burn-in phase, used to reach mutation-selection equilibrium \
`n_gens_exp`: duration of the expansion \
`n_re`: number of replicates \
`x_max_burnin`: the outward x-coordinate bound for migration during burn-in \
`x_max_exp`: the outward x-coordinate bound for migration during the expansion \
`y_max`: the upper y-coordinate bound (lower bound is always **1** currently)
`migr_mode`: mode of migration ([possible values](#migr)) \
`startfill_range`: an array of Int ranges of the coordinates that define the area to fill with individuals at start \
`data_to_generate`: string of letters representing different data to output ([possible values](#dtg)) \
`name`: world name \
`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates \
`distributed`: if **true**, distribute to threads \
If starting from existing world, also provide: \
`wld_ms1`: a spatial array of demes that contain individuals' left monosome [Bool] arrays \
`wld_ms2`: a spatial array of demes that contain individuals' right monosome [Bool] arrays \
`wld_stats`: world stats Dict

Output: a Dict containing data after the expansion: \
&nbsp;&nbsp;&nbsp;**stats** - statistics array containing world and range expansion information \
&nbsp;&nbsp;&nbsp;**fitn**, **pops**, **AAsel**, **Aasel**, **aasel**, **AAneu**, **Aaneu**, **aaneu** - data array with dimensions (space+time) that are generated if they were selected in `data_to_generate`

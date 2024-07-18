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
test = rangeexp_strip(100,1000;data_to_generate="FPSN",prop_of_sel_loci=0.8,y_max=8,migr_mode="diag1/2")
```
![alt text](https://github.com/HartreeY/RESK/blob/main/animations/readme1.gif?raw=true)

The above examples use the finite-sites model for individual genomes. For every `rangeexp` function, there is also an infinite-sites equivalent (e.g. `rangeexp_inf`). Practically, `_inf` functions are computationally faster. Some examples of `_inf` functions:

`rangeexp_cylinder_inf`:

![alt text](https://github.com/HartreeY/RESK/blob/main/animations/readme2.gif?raw=true)

`rangeexp_sphere_inf`:

![alt text](https://github.com/HartreeY/RESK/blob/main/animations/readme3.gif?raw=true)

## Main methods

**create_empty_world** \
(max=(DEF_X_MAX, DEF_Y_MAX); min=(1, 1), name=Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS"), k_capacity=DEF_K_CAPACITY, r_prolif_rate=DEF_R_PROLIF_RATE, n_loci=DEF_N_LOCI, n_sel_loci=DEF_N_SEL_LOCI,  mut_rate=DEF_MUT_RATE, migr_rate=DEF_MIGR_RATE, migr_mode=DEF_MIGR_MODE, s_sel_coef=DEF_S_SEL_COEF, h_domin_coef=DEF_H_DOMIN_COEF, prop_of_del_muts=DEF_PROP_OF_DEL_MUTS)

Builds the next generation in finite-sites expansions, i.e. advances two world arrays (left and right monosomes) by one generation and returns the new generation data for fitness, populations, mutation numbers.

---

`pnt_wld_ms1`: a spatial array of demes that contain individuals' left monosome [Bool] arrays \
`pnt_wld_ms2`: a spatial array of demes that contain individuals' right monosome [Bool] arrays \
`pnt_wld_stats`: world stats Dict \
`fitn_out`: if **true**, the new generation data for fitness will be output \
`pops_out`: if **true**, the new generation data for populations will be output \
`sel_out`: if **true**, the new generation data for selected mutations will be output \
`neu_out`: if **true**, the new generation data for neutral mutations will be output \
`max_migr`: a tuple of maximum migration area coordinates \
`migr_mode`: mode of migration. Possible values: \
&nbsp;&nbsp;&nbsp;**ort** - orthogonal directions only \
&nbsp;&nbsp;&nbsp;**all** - orthogonal and diagonal \
&nbsp;&nbsp;&nbsp;**hex** - hexagonal grid \
**diag1/2** - orthogonal and half-weighted diagonal \
**buffon1** - equidistant Buffon-Laplace (see documentation) \
**buffon2** - uniform Buffon-Laplace \
**buffon3** - inv.proportional Buffon-Laplace \
`bottleneck`: if not **NaN**, a tuple of bottleneck coordinates \
`refl_walls`: if **true**, walls reflect migrants \
`r_max_migr`: Int maximum migration radius. If **>0**, migration is kept within this radius. Can be used in addition to `max_migr` \
`r_coords`: a tuple (array) of axes' ordinal numbers that the n-sphere with `r_max_migr` covers. For example: \
**(1,3)** - migration is bound within a disk at x and z axes \
**(1,2,3)** - migration is bound within a sphere at x, y and z axes\

---
Output 1: a changed `pnt_wld_ms1` = a spatial array of demes that contain individuals' left monosome [Bool] arrays \
Output 2: a changed `pnt_wld_ms2` = a spatial array of demes that contain individuals' right monosome [Bool] arrays \
Output 3: a spatial array of demes with average fitness in the new generation \
Output 4: a spatial array of demes with populations in the new generation \
Output 5: a spatial array of demes with average selected AA mutation count in the new generation \
Output 6: a spatial array of demes with average selected Aa mutation count in the new generation \
Output 7: a spatial array of demes with average selected aa mutation count in the new generation \
Output 8: a spatial array of demes with average neutral AA mutation count in the new generation \
Output 9: a spatial array of demes with average neutral Aa mutation count in the new generation \
Output 10: a spatial array of demes with average neutral aa mutation count in the new generation

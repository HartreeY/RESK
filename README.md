# Range Expansions Simulation Kit (RESK)
A Julia library designed to efficiently simulate range expansions and analyse, as well as visualise their population dynamics.
This library has been presented in the study "RESK: An easy-to-use library for performing population genetic simulations". You can find the preprint at https://www.biorxiv.org/zzzz.

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

This repository also includes simple examples (*examples/*) and specific programs used in the aforementioned study (*study/*).

### Main use cases
Include the *resk.jl* script and use its methods.

To simulate a range expansion once, use `rangeexp` or methods that start with `rangeexp`:

- 1D: `rangeexp_ray` (=`rangeexp_1d`),`rangeexp_linear`
- 2D: `rangeexp_strip`,`rangeexp_disk`
- 3D: `rangeexp_cylinder`,`rangeexp_sphere`.

Running this with default options, a world (= habitat) will be created, seeded with individuals, and an expansion will be run on it.

The `rangeexp` functions output a fixed dictionary that includes metadata (**stats**) and expansion data. The types of expansion data within it are determined by the *data_to_generate* argument. It can take on the following values:
- **G** - **gt1**, **gt2** (genotype arrays)
- **F** - **fitn** (deme-average fitness)
- **P** - **pops** (deme populations)
- **C** - **cAA**, **cAa** and **caa** (counts of homo- [**cAA**, **caa**] and heterozygous [**cAa**] loci in a deme)
- **A** - **AA**, **Aa** and **aa** (counts of homo- [**AA**, **aa**] and heterozygous [**Aa**] loci in a deme, averaged over all loci)
The above can be combined and should be passed in a string. For example,
```
test = rangeexp_strip(15, 30; data_to_generate="ClFl", y_max=5)
```
will output:
```
OrderedDict{String, Any} with 5 entries:
  "stats" => OrderedDict{String, Any}("name"=>"2025-09-16_21-51-15", "max"=>(10…
  "fitn"  => Float32[NaN NaN … 1.0 NaN; NaN NaN … NaN NaN; … ; NaN NaN … NaN Na…
  "cAA"   => Float32[NaN NaN … NaN NaN; NaN NaN … NaN NaN; … ; NaN NaN … NaN Na…
  "cAa"   => Float32[NaN NaN … NaN NaN; NaN NaN … NaN NaN; … ; NaN NaN … NaN Na…
  "caa"   => Float32[NaN NaN … NaN NaN; NaN NaN … NaN NaN; … ; NaN NaN … NaN Na…
```
Here, **l** means **long**, and needs to be written after every output data type in *data_to_generate* if the user wants to receive the corresponding data for all generations. (Without **l**, the output will be limited to the last generation.)
These expansion data sets can be plotted and worked with. To plot expansion data, include the *reskplots.jl* file and use the unique plotting functions that start with *re_*. For example,
```
re_heatmap(test, "AA"; log_factor=1.02)
```
will output the average number of selected mutant homozygous loci in a deme:
![alt text](https://github.com/HartreeY/RESK/blob/master/img/readme0.gif?raw=true)

Here's an example of a deme-average fitness heatmap of a longer axial simulation in 2D:
```
test = rangeexp_strip(100,1000;data_to_generate="Fl",y_max=8,migr_mode="diag1/2")
re_heatmap(test, "fitn")
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

### More
For more information, please reference the [RESK wiki](https://github.com/HartreeY/RESK/wiki).


## To do
- include the possibility of multiple and partial range expansions
- add elitism, proper epistasis, fitness landscapes (rugged etc.), assortative mating
- implement migration with multiple-deme leaps
- decide on the graphical interface

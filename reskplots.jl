using Plots, GLMakie

# Plot configuration
default(margin=6Plots.mm) # Add sufficient margins to prevent cutoff

# Plotting functions
# ------------------------------------------------

"""
Shows an animated heatmap of `data` from `gen_start` to `gen_end`.

---

`data`: array with dimensions (space + time)

`gen_start`: start generation

`gen_end`: end generation

`n_gens_sub`: number of generations to subtract; e.g. set this as the number of burn-in gen-s if you wish to display the burn-in gen-s in negative numbers

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any other Plots.jl parameters
"""
function re_heatmap(data::Array, gen_start=1, gen_end=DEF_N_GENS_BURNIN + DEF_N_GENS_EXP; n_gens_sub=0, slow_factor=1, log_base=-1, clim=:default, kwargs...)
    dims = length(size(data))

    # Override default clim, since it's too hectic in animations
    if clim == :default
        no_nans = filter(!isnan, data)
        clim = minimum(no_nans)!=maximum(no_nans) ? (minimum(no_nans), maximum(no_nans)) : (0, maximum(no_nans))
    end

    @gif for i in gen_start:(gen_end*slow_factor-1)
        gen_no = trunc(Int, i / slow_factor) + 1

        if all(isnan, li(data, gen_no))
            println("No values found in any deme.")
        end

        if log_base > 0 && log_base == 1
            obj = log.(log_base, li(data, gen_no)')
        else
            obj = li(data, gen_no)'
        end

        if dims == 2 # Including time
            Plots.heatmap(obj, ylabel="Generation $(gen_no-n_gens_sub)", size=(1200, 200), yshowaxis=false, clim=clim; kwargs...)
        else
            Plots.heatmap(obj, ylabel="Generation $(gen_no-n_gens_sub)", clim=clim; kwargs...)
        end

    end
end

"""
Shows an animated heatmap of `dataname` in `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

`gen_start`: start generation

`gen_end`: end generation

`re_index`: which replicate to plot

`n_gens_sub`: number of generations to subtract; e.g. set this as the number of burn-in gen-s if you wish to display the burn-in gen-s in negative numbers

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap(re::Dict, dataname::String, gen_start=1, gen_end=re["stats"]["n_gens"]; re_index::Int = 1, n_gens_sub=re["stats"]["n_gens_burnin"], slow_factor=1, log_base=-1, clim=:default, kwargs...)
    if !isa(re[dataname], Array)
        println("This data was not generated.")
    else
        wlddim = re["stats"]["wlddim"]
        re_heatmap(re[dataname][repeat([:],wlddim)...,:,re_index], gen_start, gen_end; n_gens_sub=n_gens_sub, slow_factor=slow_factor, log_base=log_base, clim=clim, title=dataname*" #$re_index", kwargs...)
    end
end


"""
Shows population data of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`n_gens_sub`: number of generations to subtract; e.g. set this as the number of burn-in gen-s if you wish to display the burn-in gen-s in negative numbers

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_pops(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; n_gens_sub=re["stats"]["n_gens_burnin"], slow_factor=1, clim=(0, re["stats"]["capacity"]), log_base=-1, kwargs...)
    re_heatmap(re, "pops", gen_start, gen_end; n_gens_sub=n_gens_sub, slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows deme-average fitness data of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`n_gens_sub`: number of generations to subtract; e.g. set this as the number of burn-in gen-s if you wish to display the burn-in gen-s in negative numbers

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_fitn(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; slow_factor=1, clim=(0, 1), log_base=-1, kwargs...)
    re_heatmap(re, "fitn", gen_start, gen_end; slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows deme-average selected AA (homozygosity) mutations of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`n_gens_sub`: number of generations to subtract; e.g. set this as the number of burn-in gen-s if you wish to display the burn-in gen-s in negative numbers

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_AAsel(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; slow_factor=1, clim=(0, length(re["stats"]["sel_loci"])), log_base=-1, kwargs...)
    re_heatmap(re, "AAsel", gen_start, gen_end; slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows deme-average selected Aa (heterozygosity) mutations of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_Aasel(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; slow_factor=1, clim=(0, length(re["stats"]["sel_loci"])), log_base=-1, kwargs...)
    re_heatmap(re, "Aasel", gen_start, gen_end; slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows deme-average selected aa (homozygosity) mutations of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_aasel(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; slow_factor=1, clim=(0, length(re["stats"]["sel_loci"])), log_base=-1, kwargs...)
    re_heatmap(re, "aasel", gen_start, gen_end; slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows deme-average neutral AA (homozygosity) mutations of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_AAneu(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; slow_factor=1, clim=(0, re["stats"]["n_loci"] - length(re["stats"]["sel_loci"])), log_base=-1, kwargs...)
    re_heatmap(re, "AAneu", gen_start, gen_end; slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows deme-average neutral Aa (heterozygosity) mutations of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_Aaneu(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; slow_factor=1, clim=(0, re["stats"]["n_loci"] - length(re["stats"]["sel_loci"])), log_base=-1, kwargs...)
    re_heatmap(re, "Aaneu", gen_start, gen_end; slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows deme-average neutral aa (homozygosity) mutations of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_aaneu(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; slow_factor=1, clim=(0, re["stats"]["n_loci"] - length(re["stats"]["sel_loci"])), log_base=-1, kwargs...)
    re_heatmap(re, "aaneu", gen_start, gen_end; slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows deme-average deleterious mutations count of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_del(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; slow_factor=1, clim=(0, 50), log_base=-1, kwargs...)
    re_heatmap(re, "del", gen_start, gen_end; slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows deme-average beneficial mutations count of `re` from `gen_start` to `gen_end`.

---

`re`: range expansion results dictionary

`gen_start`: start generation

`gen_end`: end generation

`slow_factor`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap_bensel(re::Dict, gen_start=1, gen_end=re["stats"]["n_gens"]; slow_factor=1, clim=(0, 50), log_base=-1, kwargs...)
    re_heatmap(re, "ben", gen_start, gen_end; slow_factor=slow_factor, clim=clim, log_base=log_base, kwargs...)
end

"""
Shows a heatstack (3d heatmap) of `data`. 

---

`data`: a 3d array of data

`x_range`: Int range of x coordinates to show

`z_range`: Int range of z coordinates to show

`clim`: color bounds (Plots.jl's `clim` parameter)

`title`: if not empty string, define a custom title

`scene`: if specified, use a custom GLMakie scene
"""
function re_heatstack_frame(data::Array, x_range, z_range, clim; title="", scene=Figure())
    toNaN(x) = x < 0 ? NaN : x
    data = toNaN.(data)
    ax = Axis3(scene[1, 1], aspect=(1, 1, 1), elevation=π / 6)
    if title != ""
        ax.title = title
    end

    for i in z_range
        hm = GLMakie.heatmap!(ax, x_range, x_range, li(data, i), colorrange=clim, colormap=(:thermal, 0.25))
        GLMakie.translate!(hm, 0, 0, z_range[i])

        i == 1 && Colorbar(scene[1, 2], hm) # Add the colorbar once
    end

    GLMakie.zlims!(ax, minimum(z_range), maximum(z_range))
    scene
end

"""
Shows a heatstack (3d heatmap) of `dataname` in `re`. 

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

`x_range`: Int range of x coordinates to show

`z_range`: Int range of z coordinates to show

`clim`: color bounds (Plots.jl's `clim` parameter)

`title`: if not empty string, define a custom title

`scene`: if specified, use a custom GLMakie scene
"""
function re_heatstack_frame(re::Dict, dataname::String, x_range=1:re["stats"]["max"][1], z_range=1:re["stats"]["max"][3], clim=(minimum(filter(!isnan, re[dataname])), maximum(filter(!isnan, re[dataname]))); title="", scene=Figure())
    re_heatstack_frame(re[dataname], x_range, z_range, clim; title=title, scene=scene)
end

"""
Shows an animated heatstack (3d heatmap) of `data`. 

---

`data`: array with dimensions (3 + 1)

`clim`: color bounds (Plots.jl's `clim` parameter)

`x_range`: Int range of x coordinates to show

`z_range`: Int range of z coordinates to show

`title`: if not empty string, define a custom title

`n_gens_burnin`: number of burn-in next_generation_size
"""
function re_heatstack(data::Array, gen_start, gen_end; clim=NaN, x_range=1:1, z_range=1:1, title="", n_gens_burnin=0)
    scene = Figure()

    if isnan(clim)
        no_nans = filter(!isnan, data)
        clim = (minimum(no_nans), maximum(no_nans))
        println(clim)
    end

    record(scene,  Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH-MM-SS") * ".mp4") do io
        for i in gen_start:gen_end
            if title == ""
                ti = "Generation " * string(i - n_gens_burnin)
            end
            re_heatstack_frame(li(data, i), x_range, z_range, clim; title=ti, scene=scene)
            recordframe!(io)
            empty!(scene)
        end
    end
end

"""
Shows an animated heatstack (3d heatmap) of `dataname` in `re`. 

---

`re`: range expansion results dictionary

`dataname`: name of data in `re`

`clim`: color bounds (Plots.jl's `clim` parameter)

`x_range`: Int range of x coordinates to show

`z_range`: Int range of z coordinates to show

`title`: if not empty string, define a custom title

`n_gens_burnin`: number of burn-in next_generation_size

`kwargs...`: any Plots.jl parameters
"""
function re_heatstack(re::Dict, dataname::String, gen_start=1, gen_end=re["stats"]["n_gens"]; re_index::Int = 1, clim=NaN, x_range=1:re["stats"]["max"][1], z_range=1:re["stats"]["max"][3], title="", n_gens_burnin=re["stats"]["n_gens_burnin"], kwargs...)
    if !isa(re[dataname], Array)
        println("This data was not generated.")
    else
        wlddim = re["stats"]["wlddim"]
        re_heatstack(re[dataname][repeat([:],wlddim)...,:,re_index], gen_start, gen_end; clim=clim, x_range=x_range, z_range=z_range, title=title, n_gens_burnin=n_gens_burnin, kwargs...)
    end
end

# ------------------------------------------------

println("RESKPlots successfully loaded.")
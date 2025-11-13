module RESK2D

using RESK
using Plots
using Dates



# Plotting functions
# ------------------------------------------------


"""
Shows an animated heatmap of `data` from `gen_start` to `gen_end`.

---

`data`: array with dimensions (space + time)

`gen_start`: start generation

`gen_end`: end generation

`n_gens_sub`: number of generations to subtract; e.g. set this as the number of burn-in gen-s if you wish to display the burn-in gen-s in negative numbers

`animspeed`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`clim`: color bounds (Plots.jl's `clim` parameter)

`kwargs...`: any other Plots.jl parameters
"""
function RESK.re_heatmap(data::Array, gen_start=1, gen_end=last(size(data)); fileout::Union{String,Nothing}=nothing, n_gens_sub=0, animspeed=1, hex=false, log_base=-1, clim=:default, kwargs...)
    dims = length(size(data))
    
    # Override default clim, since it's not fixed in animations
    if clim == :default
        no_nans = filter(!isnan, data)
        clim = minimum(no_nans)!=maximum(no_nans) ? (minimum(no_nans), maximum(no_nans)) : (0, maximum(no_nans))
        println(clim)
    end

    if hex
        RESK.re_heatmap_hex(data, gen_start, gen_end; fileout=fileout, n_gens_sub=n_gens_sub, animspeed=animspeed, log_base=log_base, clim=clim, kwargs...)
    else
        an = @animate for gen_no in gen_start:round(Int,animspeed):gen_end

            if all(isnan, li(data, gen_no))
                println("No values found in any deme.")
                continue
            end

            if log_base > 0 && log_base == 1
                obj = log.(log_base, li(data, gen_no)')
            else
                obj = li(data, gen_no)'
            end
            
            if dims == 2 # Including time
                Plots.heatmap(obj, ylabel="Generation $(gen_no-n_gens_sub)", size=(1200, 200), yshowaxis=false, clim=clim, margin=6Plots.mm; kwargs...)
            else
                Plots.heatmap(obj, ylabel="Generation $(gen_no-n_gens_sub)", clim=clim, margin=6Plots.mm; kwargs...)
            end
            
        end
        if isnothing(fileout)
            gif(an, fps=20)
        else
            # write the gif to the given file
            gif(an, fileout, fps=20) |> _ -> nothing
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

`animspeed`: number of animation frames per generation

`log_base`: if not **-1**, color shows log values with this as base

`defc`: if **true**, overrides automatically defined `clim` (colour bounds) with RESK's predefined ones, depending on the type of data

`kwargs...`: any Plots.jl parameters
"""
function RESK.re_heatmap(re, dataname::String, gen_start=1, gen_end=last(size(re[dataname])); re_index::Int = 1, n_gens_sub=re["stats"]["n_gens_burnin"], hex=re["stats"]["migr_mode"]=="hex" ? true : false, animspeed=1, log_base=-1, defc=false, clim=:default, kwargs...)
    if !isa(re[dataname], Array)
        println("This data does not exist / was not selected for generation.")
    else
        wlddim = re["stats"]["wlddim"]

        if defc
            if dataname=="pops"
                clim=(0, re["stats"]["capacity"])
            elseif dataname=="fitn"
                clim=(0, 1)
            elseif dataname=="AAsel" || dataname=="Aasel" || dataname=="aasel" || dataname=="AAneu" || dataname=="Aaneu" || dataname=="aaneu"
                clim=(0, length(re["stats"]["sel_loci"]))
            elseif dataname=="del"
                clim=(0, re["stats"]["n_loci"]*re["stats"]["prop_of_del_muts"])
            elseif dataname=="ben"
                clim=(0, re["stats"]["n_loci"]*(1-re["stats"]["prop_of_del_muts"]))
            end
        end

        re_heatmap(re[dataname][repeat([:],wlddim)...,:,re_index], gen_start, gen_end; n_gens_sub=n_gens_sub, animspeed=animspeed, log_base=log_base, clim=clim, title=dataname*" #$re_index", 
        hex = hex,
        kwargs...)
    end
end


end
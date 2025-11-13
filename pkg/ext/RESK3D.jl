module RESK3D

using RESK
using GLMakie
using Dates

"""
Shows a heatstack (3d heatmap) of `data`. 

---

`data`: a 3d array of data

`x_range`: Int range of x coordinates to show

`z_range`: Int range of z coordinates to show

`clim`: color bounds (Plots.jl's `clim` parameter) NOT ACTUALLY PLOTS.JL

`title`: if not empty string, define a custom title

`scene`: if specified, use a custom GLMakie scene
"""
function RESK.re_heatstack_frame(data::Array, x_range, z_range, clim; title="", scene=Figure())
    
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
function RESK.re_heatstack_frame(re, dataname::String, x_range=1:re["stats"]["max"][1], z_range=1:re["stats"]["max"][3], defc=false, 
    clim=(minimum(filter(!isnan, re[dataname])), maximum(filter(!isnan, re[dataname]))); title="", scene=Figure())

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
function RESK.re_heatstack(data::Array, gen_start=1, gen_end=last(size(data)); clim::Union{Tuple,Nothing}=nothing, x_range=1:size(data,1), z_range=1:size(data,3), title="", n_gens_burnin=0)
    scene = Figure()

    if isnothing(clim)
        no_nans = filter(!isnan, data)
        clim = minimum(no_nans)==maximum(no_nans) ? (0,1) : (minimum(no_nans), maximum(no_nans))
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
function RESK.re_heatstack(re, dataname::String, gen_start=1, gen_end=re["stats"]["n_gens"]; re_index::Int = 1, defc=false, clim=nothing, x_range=1:re["stats"]["max"][1], z_range=1:re["stats"]["max"][3], title="", n_gens_burnin=re["stats"]["n_gens_burnin"], kwargs...)
    if !isa(re[dataname], Array)
        println("This data was not generated.")
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

        re_heatstack(re[dataname][repeat([:],wlddim)...,:,re_index], gen_start, gen_end; clim=clim, x_range=x_range, z_range=z_range, title=title, n_gens_burnin=n_gens_burnin, kwargs...)
    end
end

end
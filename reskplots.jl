using Plots, GLMakie, Luxor

# Plot configuration
default(margin=6Plots.mm) # Add sufficient margins to prevent cutoff

toNaN(x) = x < 0 ? NaN : x
toNaNzero(x) = x <= 0 ? NaN : x

# Plotting functions
# ------------------------------------------------

"""
Custom function for ticks in hexagonal plots in `re_heatmap`.
    
---

`xx`: if **true**, the ticks are for the x-axis, else - for the y-axis
"""
function hextickfun(n, pos, xx=true; startnumber=30, finishnumber=40, nticks=10)
    @layer begin
        Luxor.translate(pos)
        fontsize(12)
        ticklength = get_fontsize()
        if xx
            line(O, O + polar(ticklength, 3π/2), :stroke)
            offs = (0, -ticklength*2)
        else
            line(O + polar(-ticklength, 3π/2), O, :stroke)
            offs = (0, ticklength*2)
        end
        k = rescale(n, 0, nticks - 1, startnumber, finishnumber)
        Luxor.text("$(trunc(Int,k))",
            O + offs,
            halign=:center,
            valign=:middle)
    end
end
hextickfun_y(n, pos; startnumber=30, finishnumber=40, nticks=10) = hextickfun(n, pos, false; startnumber=startnumber, finishnumber=finishnumber, nticks=nticks)



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
function re_heatmap(data::Array, gen_start=1, gen_end=DEF_N_GENS_BURNIN + DEF_N_GENS_EXP; n_gens_sub=0, slow_factor=1, log_base=-1, hex=false, clim=:default, kwargs...)
    dims = length(size(data))

    # Override default clim, since it's not fixed in animations
    if clim == :default
        no_nans = filter(!isnan, data)
        clim = minimum(no_nans)!=maximum(no_nans) ? (minimum(no_nans), maximum(no_nans)) : (0, maximum(no_nans))
        println(clim)
    end

    if hex
        side = 30
        shortd = 1.73205*side
        if dims == 3 # Including time
            xmax = size(data)[1]
            ymax = size(data)[2]

            function frame(scene, framenumber)
                gen_no = trunc(Int, framenumber / slow_factor) + 1
                background("white")
                origin(0,0)
                tickline(Luxor.Point(104, 61), Luxor.Point(104+shortd*(xmax-1), 61),major=xmax-2,startnumber=1, finishnumber=xmax,  major_tick_function = hextickfun) # x-axis
                tickline(Luxor.Point(51, 92), Luxor.Point(51, 92+1.5*side*(ymax-1)),major=ymax-2,startnumber=1, finishnumber=ymax, major_tick_function = hextickfun_y) # y-axis
                for q in 1:xmax # vertical
                    for r in 1:ymax # horizontal
                        pgon = hextile(Luxor.HexagonOffsetEvenR(q+1, r+1, 30))
                        sethue(Luxor.HSB((toNaNzero(data[q,r,framenumber])-clim[1])*100/(clim[2]-clim[1]),1,1))
                        Luxor.poly(pgon, :fill)
                    end
                end
                
                cb = blend(Luxor.Point(135,145+1.5*side*(ymax-1)), Luxor.Point(185,145+1.5*side*(ymax-1)), Luxor.HSB(0,1,1), Luxor.HSB(100,1,1))
                setblend(cb)
                polysmooth(box(Luxor.Point(160,145+1.5*side*(ymax-1)), 50, 20, :fill), 4, action = :stroke)

                setcolor("black")
                fontsize(18)
                Luxor.text("$(round(clim[1];digits=2))", Luxor.Point(120,150+1.5*side*(ymax-1)),halign=:right)
                Luxor.text("$(round(clim[2];digits=2))", Luxor.Point(195,150+1.5*side*(ymax-1)))

                Luxor.text("Generation $(gen_no-n_gens_sub)", Luxor.Point(295,150+1.5*side*(ymax-1)))
            end
            
            demo = Movie(trunc(Int,122+shortd*xmax),trunc(Int,90+shortd*ymax), "test.gif", 1:gen_end)            
            Luxor.animate(demo, Luxor.Scene(demo, frame, gen_start:(gen_end*slow_factor-1)); creategif=true)
        end
    else

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

`defc`: if **true**, overrides automatically defined `clim` (colour bounds) with RESK's predefined ones, depending on the type of data

`kwargs...`: any Plots.jl parameters
"""
function re_heatmap(re::Dict, dataname::String, gen_start=1, gen_end=re["stats"]["n_gens"]; re_index::Int = 1, n_gens_sub=re["stats"]["n_gens_burnin"], slow_factor=1, log_base=-1, defc=false, clim=:default, kwargs...)
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

        re_heatmap(re[dataname][repeat([:],wlddim)...,:,re_index], gen_start, gen_end; n_gens_sub=n_gens_sub, slow_factor=slow_factor, log_base=log_base, clim=clim, title=dataname*" #$re_index", 
        hex = re["stats"]["migr_mode"]=="hex" ? true : false,
        kwargs...)
    end
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
function re_heatstack_frame(re::Dict, dataname::String, x_range=1:re["stats"]["max"][1], z_range=1:re["stats"]["max"][3], defc=false, 
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
function re_heatstack(data::Array, gen_start=1, gen_end=size(data,4); clim=NaN, x_range=1:size(data,1), z_range=1:size(data,3), title="", n_gens_burnin=0)
    scene = Figure()

    if isnan(clim)
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
function re_heatstack(re::Dict, dataname::String, gen_start=1, gen_end=re["stats"]["n_gens"]; re_index::Int = 1, defc=false, clim=NaN, x_range=1:re["stats"]["max"][1], z_range=1:re["stats"]["max"][3], title="", n_gens_burnin=re["stats"]["n_gens_burnin"], kwargs...)
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

# ------------------------------------------------

println("RESKPlots successfully loaded.")
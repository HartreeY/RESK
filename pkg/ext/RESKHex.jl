module RESKHex

using RESK
using Luxor
using Dates

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

function RESK.re_heatmap_hex(data::Array, gen_start=1, gen_end=last(size(data)); fileout::Union{String,Nothing}=nothing, n_gens_sub=0, animspeed=1, log_base=-1, clim=:default, kwargs...)
    
    side = 30
    shortd = 1.73205*side

    xmax = size(data)[1]
    ymax = size(data)[2]

    function frame(scene, framenumber)
        gen_no = trunc(Int, framenumber / animspeed) + 1
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
    Luxor.animate(demo, Luxor.Scene(demo, frame, gen_start:round(Int,animspeed):gen_end); creategif=true)


end

end
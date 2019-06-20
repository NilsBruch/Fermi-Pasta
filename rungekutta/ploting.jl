using Plots
using DelimitedFiles

theme(:juno)

xreso = 1280
yreso = 720

filename="rungekutta/data/malsehen"

sol=readdlm(filename*"_y", '\t', Float64, '\n'; skipstart=1)
energydata = readdlm(filename*"_E", '\t', Float64, '\n'; skipstart=1)
zeit = energydata[:,1]
E_full= energydata[:,2]
T= energydata[:,3]
Pot= energydata[:,4]

N, α, β, t, Δt=readdlm(filename*"_p", '\t', Float64, '\n'; skipstart=1)
N=Int(N)

Δt_per_f=1000
fps=30
50/Δt/Δt_per_f

anim = @animate for t in 1:Δt_per_f:length(zeit)
    p1 = plot(sol[t, 1:N], ylims = (-1,1),
    label = "Atoms",
    xlabel = "x",
    ylabel = "Displacement",
    linewidth = 0,
    markershape = :hexagon,
    markersize = 15,
    line = nothing,
    size=(xreso,yreso))

    plotlength = Δt_per_f*fps*Δt

    #if t<plotlength+1
    #    tp = 1
    #else
    #    tp = t-plotlength
    #end

    p2 = plot(zeit[1:t],E_full[1:t], label = "Total Energy",
    ylims = (0,E_full[1]+0.2),
    xlims = (-plotlength*3/4+Δt*t, plotlength*1/4+Δt*t),
    xlabel = "Time",
    ylabel = "Energy",
    legend = :bottom,
    size=(xreso,yreso))
    plot!(zeit[1:t],T[1:t], label = "Kintetic Energy")
    plot!(zeit[1:t], Pot[1:t], label = "Potential Energy")

    plot(p1,p2,layout = (2,1))
end

gif(anim, "first_gif.gif", fps = 30)


plot(zeit, E_full, label = "Total Energy", ylims = (0,1.5),xlabel = "Time",ylabel = "Energy", legend = :bottom, size=(xreso,yreso))
plot!(zeit, T, label = "Kintetic Energy")
plot!(zeit, Pot, label = "Potential Energy")

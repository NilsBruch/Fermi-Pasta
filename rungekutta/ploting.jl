using Plots
using DelimitedFiles

theme(:juno)

xreso = 1280
yreso = 720

filename="test" #"./rungekutta/data/test"

sol=readdlm(filename*"_y", '\t', Float64, '\n'; skipstart=1)
energydata = readdlm(filename*"_E", '\t', Float64, '\n'; skipstart=1)
zeit = energydata[:,1]
E_full= energydata[:,2]
T= energydata[:,3]
Pot= energydata[:,4]

N, α, β, t, Δt=readdlm(filename*"_p", '\t', Float64, '\n'; skipstart=1)
N=Int(N)

@gif for t in 1:length(zeit)
    p1 = plot(sol[t, 1:N], ylims = (-1,1),
    label = "Atoms",
    xlabel = "x",
    ylabel = "Displacement",
    linewidth = 0,
    markershape = :hexagon,
    markersize = 15,
    line = nothing,
    size=(xreso,yreso))

    plotlength = 2000
    if t<plotlength+1
        tp = 1
    else
        tp = t-plotlength
    end
    xaxis = range(tp,stop = t,step = 1)
    p2 = plot(xaxis,E_full[tp:t], label = "Total Energy",
    ylims = (0,1.5),
    xlabel = "Time",
    ylabel = "Energy",
    legend = :bottom,
    size=(xreso,yreso))
    plot!(xaxis,T[tp:t], label = "Kintetic Energy")
    plot!(xaxis, Pot[tp:t], label = "Potential Energy")



    plot(p1,p2,layout = (2,1))
end

plot(zeit, E_full, label = "Total Energy", ylims = (0,1.5),xlabel = "Time",ylabel = "Energy", legend = :bottom, size=(xreso,yreso))
plot!(zeit, T, label = "Kintetic Energy")
plot!(zeit, Pot, label = "Potential Energy")

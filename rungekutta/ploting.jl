using Plots
using DelimitedFiles

theme(:juno)

xreso = 1280
yreso = 1000

filename="rungekutta/data/spectrum_harmonic"

sol=readdlm(filename*"_y", '\t', Float64, '\n'; skipstart=1)
energydata = readdlm(filename*"_E", '\t', Float64, '\n'; skipstart=1)
zeit = energydata[:,1]
E_full= energydata[:,2]
T= energydata[:,3]
Pot= energydata[:,4]
spec = readdlm(filename*"_spec", '\t', Float64, '\n'; skipstart=1)

#for t in 1:length(zeit)
#    for k in 1:length(spec[1,:])
#        if


N, α, β, t, Δt=readdlm(filename*"_p", '\t', Float64, '\n'; skipstart=1)
N=Int(N)

anim = @animate for t in 1:length(zeit)
    #Plot 
    p1 = plot(sol[t, :], ylims = (-1,1),
    label = "Atoms",
    xlabel = "x",
    ylabel = "Displacement",
    linewidth = 0,
    markershape = :hexagon,
    markersize = 1,
    line = nothing,
    size=(xreso,yreso))

    p2 = plot(spec[t, :],
    xlabel="Wave number k",
    ylabel="Occupation",
    ylim=(0, 20),
    line=nothing,
    legend=nothing,
    marker=:circle,
    markersize=2)

    fps=30
    Δt_per_f=zeit[2]-zeit[1]
    plotlength = Δt_per_f*fps*Δt

    p3 = plot(zeit[1:t],E_full[1:t], label = "Total Energy",
    ylims = (0,E_full[1]+0.2),
    xlims = (-plotlength*9/10+zeit[t], plotlength*1/10+zeit[t]),
    xlabel = "Time",
    ylabel = "Energy",
    legend = :bottom,
    size=(xreso,yreso))
    plot!(zeit[1:t],T[1:t], label = "Kintetic Energy")
    plot!(zeit[1:t], Pot[1:t], label = "Potential Energy")

    plot(p1,p2,p3,layout = (3,1))
    println(t)
end

gif(anim, "rungekutta/figures/thermalization_harmonic.gif", fps = 30)


plot(zeit, E_full, label = "Total Energy", ylims = (0,1.5),xlabel = "Time",ylabel = "Energy", legend = :bottom, size=(xreso,yreso))
plot!(zeit, T, label = "Kintetic Energy")
plot!(zeit, Pot, label = "Potential Energy")

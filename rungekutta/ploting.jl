xreso = 1280
yreso = 720


sol=readdlm("./rungekutta/data/test_y", '\t', Float64, '\n'; skipstart=2)
energydata = readdlm("./rungekutta/data/test_E", '\t', Float64, '\n'; skipstart=1)
time=
E_full=
T=
Pot=

@gif for t in 1:M
    t_sim = t*T/M

    p1 = plot(sol(t_sim)[1:N], ylims = (-1,1),
    label = "Atoms",
    xlabel = "x",
    ylabel = "Displacement",
    linewidth = 0,
    markershape = :hexagon,
    markersize = 15,
    line = nothing,
    size=(xreso,yreso))

    plotlength = 20
    if t<plotlength+1
        tp = 1
    else
        tp = t-plotlength
    end
    xaxis = range(tp,stop = t,step = 1)
    p2 = plot(xaxis,energyValues[tp:t], label = "Total Energy",
    ylims = (0,0.3),
    xlabel = "Time",
    ylabel = "Energy",
    legend = :bottom,
    size=(xreso,yreso))
    plot!(xaxis,kinteticValues[tp:t], label = "Kintetic Energy")
    plot!(xaxis,potentialValues[tp:t], label = "Potential Energy")



    plot(p1,p2,layout = (2,1))
end


@sync begin
    @async TimeEvo(cfg_harm, t)
    @async TimeEvo(cfg_anharm, t)
end

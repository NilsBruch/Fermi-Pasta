using DifferentialEquations
using Plots
using Distributed
include("structs.jl")
include("energy.jl")
theme(:juno)

addprocs(2)
#Hello Nils
function main()

    N = 10
    t = 10
    Δt = 0.0001

    cfg_harm = ChainConstructor(NumberOfAtoms = N, T = t).getHarmonic;
    cfg_anharm = ChainConstructor(NumberOfAtoms = N, T = t).getAnHarmonic;


    function OneDChain(dy,y,p,t)
        α,β = p[1],p[2]

        for i in 1:N
            if i == 1
                dy[i]       = dx₁ = 0
                dy[i + N]   = dv₁ = 0
            elseif i == N
                dy[i]       = dxₙ = 0
                dy[i + N]   = dvₙ = 0
            else
                dy[i]       = dxᵢ = y[i + N]
                dy[i + N]   = dvᵢ = α*(y[i+1]+y[i-1]-2*y[i]) + β*((y[i+1]-y[i])^2-(y[i]-y[i-1])^2)
            end
        end
    end


    function TimeEvo(cfg::config, T)

        N = cfg.NumberOfAtoms
        M = Int(T/Δt)
        tspan = cfg.time
        u0 = cfg.u0
        p = [cfg.α,cfg.β]

        energyValues = zeros(M)
        kinteticValues = zeros(M)
        potentialValues = zeros(M)

        prob = ODEProblem(OneDChain,u0,tspan,p)
        sol = solve(prob)

        for t in 1:M
            t_sim = t*T/M
            displacement = sol(t_sim)[1:N]
            velocity     = sol(t_sim)[N+1:2*N]
            E = Energy(displacement,velocity,cfg)

            kinteticValues[t]  = E[1]
            potentialValues[t] = E[2]
            energyValues[t]    = E[3]
        end

        xreso = 1280
        yreso = 720
        @gif for t in 1:2000:M

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
    end

    @sync begin
        @async TimeEvo(cfg_harm, t)
        @async TimeEvo(cfg_anharm, t)
    end
end




main()

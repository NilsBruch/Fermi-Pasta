using DifferentialEquations
using Plots
include("structs.jl")
include("energy.jl")
theme(:juno)

#Hello Nils
function main()

    N = 10
    t = 100
    Δt = 0.1

    cfg = ChainConstructor(NumberOfAtoms = N, T = t).getHarmonic;


    function OneDChain(dy,y,p,t)
        for i in 1:N
            if i == 1
                dy[i] = 0
                dy[i + N] = 0
            elseif i == N
                dy[i] = 0
                dy[i + N] = 0
            else
                dy[i] = y[i + N]
                dy[i + N] = (y[i+1]+y[i-1]-2*y[i])
            end
        end
    end


    function TimeEvo(cfg::config, T)

        N = cfg.NumberOfAtoms
        M = Int(T/Δt)
        tspan = cfg.time
        u0 = cfg.u0

        energyValues = zeros(M)
        kinteticValues = zeros(M)
        potentialValues = zeros(M)

        prob = ODEProblem(OneDChain,u0,tspan)
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

        @gif for t in 1:M
            t_sim = t*T/M

            p1 = plot(sol(t_sim)[1:N], ylims = (-1,1),
            label = "Atoms",
            xlabel = "x",
            ylabel = "Displacement",
            linewidth = 0,
            markershape = :hexagon,
            markersize = 10,
            line = nothing)

            p2 = plot(energyValues[1:t], label = "Total Energy",
            ylims = (0,0.3),
            xlabel = "Time",
            ylabel = "Energy")
            plot!(kinteticValues[1:t], label = "Kintetic Energy")
            plot!(potentialValues[1:t], label = "Potential Energy")



            display(plot(p1,p2,layout = (2,1)))
        end
    end


    TimeEvo(cfg, t)
end

main()

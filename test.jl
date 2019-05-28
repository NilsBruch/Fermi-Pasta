using DifferentialEquations
using ParameterizedFunctions

include("structs.jl")
include("energy.jl")


N = 10
t = 100
Δt = 0.1

cfg = ChainConstructor(NumberOfAtoms = N, T = t).getHarmonic;


tspan = cfg.time
u0 = cfg.u0



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
            dy[i + N]   = dvᵢ = α*(y[i+1]+y[i-1]-2*y[i]) + β*0
        end
    end
end


p = [cfg.α,cfg.β]


prob = ODEProblem(OneDChain,u0,tspan,p)
sol = solve(prob)

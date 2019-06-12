include("structs_runge.jl")
using DelimitedFiles

function main()
    #Initialize parameters
    N = 10      #Number of atoms
    α = 1.0     #harmonic coupling
    β = 0.0       #anharmonic coupling

    t = 2       #time to solve the equation for
    Δt = 0.1    #timestep width
    nstep=Int(t/Δt)

    filename="./rungekutta/data/test" #Output filename
    #Starting conditions
    u0=zeros(2*N) #u0[1:N]: x0[1:N],   u0[N+1:2N]=v0
    u0[2]=1

    #Configuration
    cfg = config(N, α, β, u0)
    #timestep(Δt, cfg)
    #Solve equation and output to file
    io_y = open(filename*"_y", "w") #Output cfg.y
    write(io_y, "Parameters in first line (N, alpha, beta, t_max, delta t). Solution after second line( u_1(t1), ..., u_N(t1)) \n") #Header
    writedlm(io, [N, α, β, t, Δt]') #Parameters

    io_E = open(filename*"_E", "w") #Output Energies
    write(io_E, "t | E_total | E_kin |E_pot \n")

    for i in 1:nstep
        t_sim = (i-1)*Δt
        T, Pot, E = Energy(cfg) #kinetic, potential and total energy
        writedlm(io_y, cfg.y') #Write out y
        writedlm(io_E, [t_sim, E, T, Pot]')
        timestep(Δt, cfg) #Calculate cfg.y(t+Δt)
    end
    close(io_y)
end

main()

include("structs_runge.jl")
using DelimitedFiles
using FFTW

function main()
    #Initialize parameters
    N = 50      #Number of atoms
    α = 1.0     #harmonic coupling
    β = 1.0       #anharmonic coupling

    t = 200      #time to solve the equation for
    Δt = 0.0001    #timestep width
    Δt_per_f=1000  #Timesteps per output/frame
    filename="rungekutta/data/spectrum" #Output filename




    #Starting conditions
    u0=zeros(2*N) #u0[1:N]: x0[1:N],   u0[N+1:2N]=v0
    x=range(0, N-1; step=1)
    u0[1:N]=sin.(16*π/(N-1)*x)

    #Initialize configuration
    cfg = config(N, α, β, u0)

    #Prepare output files:
    io_par = open(filename*"_p", "w")  #Parameters
    write(io_par, "N, alpha, beta, t, deltat \n")
    writedlm(io_par, [N, α, β, t, Δt]')
    close(io_par)

    io_y = open(filename*"_y", "w") #Displacements
    write(io_y, "(u_1(t1), ..., u_N(t1)) \n") #Header

    io_E = open(filename*"_E", "w") #Energies
    write(io_E, "t | E_total | E_kin |E_pot \n") #Header

    io_spec = open(filename*"_spec", "w") #Spectrum
    write(io_spec, "Absoloute value of DFT of y[1:N-1] \n")

    nstep=Int(t/Δt)  #Numer of timesteps
    for i in 1:nstep
        if rem(i-1, Δt_per_f) == 0
            t_sim = (i-1)*Δt  #Time in simulation
            T, Pot, E = Energy(cfg) #kinetic, potential and total energy
            spec=abs.(fft(cfg.y[1:N-1]))

            writedlm(io_y, cfg.y[1:N]') #Write out y
            writedlm(io_E, [t_sim, E, T, Pot]')  #write out energy
            writedlm(io_spec, spec[1:Int(N/2)+1]') #Write out spectrum

        end
        timestep(Δt, cfg) #Calculate cfg.y(t+Δt)
    end

    #Close output streams
    close(io_y)
    close(io_E)
    close(io_spec)
end

main()

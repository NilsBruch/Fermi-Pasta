include("structs_runge.jl")
using DelimitedFiles

#Initialize parameters
N = 10      #Number of atoms
α = 1.0     #harmonic coupling
β = 0.0       #anharmonic coupling

t = 2       #time to solve the equation for
Δt = 0.1    #timestep width
nstep=Int(t/Δt)

filename="./rungekutta/data/filename.txt" #Output filename

#Starting conditions
u0=zeros(2*N) #u0[1:N]: x0[1:N],   u0[N+1:2N]=v0
u0[2]=1

#Configuration
cfg = config(N, α, β, u0)

#Solve equation and output to file
io = open(filename, "w") #Output stream
write(io, "Parameters in first line (N, alpha, beta, t_max, delta t). Solution in second line( u_1(t1), ..., u_N(t1)) \n") #Header
writedlm(io, [N, α, β, t, Δt]') #Parameters
for i in 1:nstep
    writedlm(io, cfg.y[1:N]')
    timestep(Δt, cfg)
end
close(io)

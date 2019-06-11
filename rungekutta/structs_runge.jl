mutable struct config #type holding configuration and parameters
    NumberOfAtoms::Int64
    α ::Float64
    β ::Float64
    y ::Array #y[1:N]=x, y[N+1:N]=v
end

#Write out right hand side of y'=f(y), dy will contain the derivative
function f(dy, y, N, α, β)

    dy[1]=dy[N]=dy[N+1]=dy[2*N]=0 #Fixed boundary conditions (at 0)

    for i in 2:N-1
        dy[i]=y[i+N] #x'=v
        dy[i+N]=α*(y[i+1]+y[i-1]-2*y[i])+β*((y[i+1]-y[i])^3-(y[i]-y[i-1])^3) #x''= e.o.m...
    end
end

#Calculate timestep using the Runge-Kutta method
function timestep(h, cfg::config)
    #Calculate coefficients
    dy=zeros(2*cfg.NumberOfAtoms)
    f(dy, cfg.y, cfg.NumberOfAtoms, cfg.α, cfg.β)
    k1=dy

    f(dy, cfg.y+h/2*k1, cfg.NumberOfAtoms, cfg.α, cfg.β)
    k2=dy

    f(dy, cfg.y+h/2*k2, cfg.NumberOfAtoms, cfg.α, cfg.β)
    k3=dy

    f(dy, cfg.y+h*k3, cfg.NumberOfAtoms, cfg.α, cfg.β)
    k4=dy

    #Reasign value
    cfg.y+=h/6*(k1+2*k2+2*k3+k4)
end

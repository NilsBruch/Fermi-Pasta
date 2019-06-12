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

function Energy(cfg::config)
    N=cfg.NumberOfAtoms

    function ScalarProduct(v,w)
        return sum((v.*w))
    end

    function KineticEnergy(Geschwindigkeiten::Array)
        return 1/2*ScalarProduct(Geschwindigkeiten,Geschwindigkeiten)
    end



    function Potential(Auslenkungen)
        function V(i)
            return cfg.α/2*(Auslenkungen[i]-Auslenkungen[i+1])^2+cfg.β/4*(Auslenkungen[i]-Auslenkungen[i+1])^4
        end

        s = 0
        for i in 1:N-1
            s = s + V(i)
        end
        return s
    end

    T = KineticEnergy(cfg.y[N+1:2*N])
    Pot = Potential(cfg.y[1:N])
    E = T+Pot
    return [T,Pot,E]
end

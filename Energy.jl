function Energy(Auslenkungen::Array,Geschwindigkeiten::Array,cfg::config)
    N = length(Auslenkungen)

    α = cfg.α
    β = cfg.β
    function ScalarProduct(v,w)
        return sum((v.*w))
    end

    function KineticEnergy(Geschwindigkeiten::Array)
        return 1/2*ScalarProduct(Geschwindigkeiten,Geschwindigkeiten)
    end

    function V(i)
        return α/2*(Auslenkungen[i]-Auslenkungen[i+1])^2
    end

    function Potential(Auslenkungen)
        s = 0
        for i in 1:N-1
            s = s + V(i)
        end
        return s
    end
    T = KineticEnergy(Geschwindigkeiten)
    Pot = Potential(Auslenkungen)
    E = KineticEnergy(Geschwindigkeiten) + Potential(Auslenkungen)
    return [T,Pot,E]
end

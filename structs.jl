
struct config
    NumberOfAtoms::Int
    α ::Float64
    β ::Float64

    u0 ::Array
    time ::Tuple{Float64,Float64}

    f(x) = sin(π*x/(N-1))

    function config(; NumberOfAtoms = 10, α = 1, β = 1,t = 10)
        new(NumberOfAtoms,α,β,append!(f.(range(0,stop=N-1)),zeros(N)),(0.0,t))
    end

end


struct ChainConstructor
    getHarmonic ::config
    getAnHarmonic ::config

    function getHarmonic(N,T)
        cfg = config(NumberOfAtoms =  N, α = 1, β = 0, t = T)
        #cfg.u0[2] = 0.5
        return cfg
    end

    function getAnHarmonic(N,T)
        cfg = config(NumberOfAtoms =  N, α = 1, β = 0.3, t = T)
        #cfg.u0[2] = 0.5
        return cfg
    end

    function ChainConstructor(; NumberOfAtoms = 10, T = 10)
        new(getHarmonic(NumberOfAtoms, T), getAnHarmonic(NumberOfAtoms, T))
    end
end

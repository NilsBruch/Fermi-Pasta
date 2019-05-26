
using Plots
theme(:juno)
N = 10
Δt = 0.1
α = 0.1
β = 0

Auslenkungen = zeros(N)
Auslenkungen[2] = 0.8
Geschwindigkeiten = zeros(N)
function Euler(Auslenkungen::Array, Geschwindigkeiten::Array)


    function CreateTimeStep(Auslenkungen::Array)
        f = zeros(N)
        f[1] = 0 #α*(Auslenkungen[2] - Auslenkungen[1])
        f[N] = 0 #α*(-Auslenkungen[N]+Auslenkungen[N-1])
        for i in 2:N-1
            f[i] = α*(Auslenkungen[i+1] + Auslenkungen[i-1] - 2*Auslenkungen[i])
        end
        return f
    end

    Auslenkungen = copy(Auslenkungen)
    Geschwindigkeiten = copy(Geschwindigkeiten)

    y_i = [Auslenkungen, Geschwindigkeiten]
    f = [Geschwindigkeiten,CreateTimeStep(Auslenkungen)]
    y_p = y_i + Δt*f

    return y_p
end

function TimeEvoEuler(Time::Int, Ausl_Initial::Array, Geschw_Initial::Array)

    x = Euler(Ausl_Initial,Geschw_Initial)

    for t in 1:Time-1
        newAuslenkung = x[1:N]
        newGeschwindigkeit = x[N+1:2*N]
        x = Euler(newAuslenkung,newGeschwindigkeit)
        display(plot(x[1:N],ylims= (-1,1)))
    end
end
#TimeEvoEuler(300,Auslenkungen,Geschwindigkeiten)



function RungeKutta(Auslenkungen::Array, Geschwindigkeiten::Array)


    function CreateTimeStep(Auslenkungen::Array)

        f = zeros(N)
        f[1] = 0 #α*(Auslenkungen[2] - Auslenkungen[1])
        f[N] = 0 #α*(-Auslenkungen[N]+Auslenkungen[N-1])
        for i in 2:N-1
            f[i] = α*(Auslenkungen[i+1] + Auslenkungen[i-1] - 2*Auslenkungen[i])
        end
        return f
    end

    Auslenkungen = copy(Auslenkungen)
    Geschwindigkeiten = copy(Geschwindigkeiten)


    y_i = [Auslenkungen, Geschwindigkeiten]
    f1 = CreateTimeStep(Auslenkungen)
    f2 = CreateTimeStep(Auslenkungen+Δt*f1)
    f = [Geschwindigkeiten,f1+f2]
    y_p = y_i + Δt/2*f

    return y_p
end



function TimeEvo(Time::Int, Ausl_Initial::Array, Geschw_Initial::Array)


    function Energy(Auslenkungen::Array,Geschwindigkeiten::Array)

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
        K = KineticEnergy(Geschwindigkeiten)
        T = Potential(Auslenkungen)
        E = KineticEnergy(Geschwindigkeiten) + Potential(Auslenkungen)
        return [K,T,E]
    end

    x = Euler(Ausl_Initial,Geschw_Initial)
    energyValues = zeros(0)
    kintetic = zeros(0)
    potential = zeros(0)
    V0 = Energy(Ausl_Initial,Geschw_Initial)[2]

    @gif for t in 1:Time-1
        newAuslenkung = x[1]
        newGeschwindigkeit = x[2]
        x = Euler(newAuslenkung,newGeschwindigkeit)
        push!(kintetic, Energy(newAuslenkung,newGeschwindigkeit)[1])
        push!(potential, Energy(newAuslenkung,newGeschwindigkeit)[2])
        push!(energyValues, Energy(newAuslenkung,newGeschwindigkeit)[3])

        p1 = plot(newAuslenkung,
                     ylims=(-1,1),
                     label = "Atoms",
                     linewidth = 0,
                     markershape = :hexagon,
                     line = nothing)
        p2 = plot(energyValues,
        ylims = (0,2*V0),
        xtick = false,
        label = "Total Energy")
        plot!(kintetic, label = "Kinetic Energy")
        plot!(potential, label = "Potential Energy")

        display(plot(p1,p2,layout=2))


    end
end

TimeEvo(100,Auslenkungen,Geschwindigkeiten)

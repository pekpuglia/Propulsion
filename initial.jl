using Unitful, NonlinearSolve, Parameters
##
@with_kw struct ThermodynamicProperties
    P
    z
    T
end

const Rmolar = 8.3144598u"J/mol/K"

function residues(;P, z, T)
    P - z*T*ustrip(u"J/mol/K", Rmolar)
end

function solveTP(;kwargs...)
    allvars = fieldnames(ThermodynamicProperties)
    if length(kwargs) != 2
        error("2 thermodynamic properties needed, $(length(kwargs)) given: $(keys(kwargs))")
    end
    if any([!(k in allvars) for k in keys(kwargs)])
        error("expected keys from $allvars, got: $(keys(kwargs)) ")
    end
    #https://www.geeksforgeeks.org/sets-in-julia/
    missingvar = (setdiff(Set(allvars), Set(keys(kwargs))) |> collect)[1]

    prob = NonlinearProblem((value, params) -> residues(;Dict(missingvar => value, kwargs...)...), 1.0)
    sol = solve(prob, NewtonRaphson())

    ThermodynamicProperties(;
    Dict(missingvar => sol.u, 
    Dict(k => Float64(v) for (k, v) in kwargs)...)...)
end
##
struct MassProperties
    tp::ThermodynamicProperties
    MM
    rho
    R
end

struct CalorificProperties
    mp::MassProperties
    cv
    cp
    gamma
    a
end

struct FlowProperties
    cp::CalorificProperties
    M
    v
    T0
    rho0
    P0
    a0
end

struct Quasi1dimflowProperties
    fp::FlowProperties
    mdot
    A
    Astar
end
##
#assume que kwargs tá completo
function residues(;kwargs...)
    p = kwargs[:p]
    rho = kwargs[:rho]
    T = kwargs[:T]
    a = kwargs[:a]
    gamma = kwargs[:gamma]
    
    p0 = kwargs[:p0]
    rho0 = kwargs[:rho0]
    T0 = kwargs[:T0]
    a0 = kwargs[:a0]

    mdot = kwargs[:mdot]
    A = kwargs[:A]
    v = kwargs[:v]
    M = kwargs[:M]
    
    

    [
        1 + (gamma - 1) / 2 * M^2 - T0/T
        (1 + (gamma -  1) / 2 * M^2)^(gamma/(gamma-1)) - p0/p
        (1 + (gamma -  1) / 2 * M^2)^(  1  /(gamma-1)) - rho0/rho
        M*a - v
        a^2 - gamma*p/rho
        a^2/(gamma-1) + v^2/2 - a0^2/(gamma-1)
        mdot - rho*v*A
    ]
end
##
params = Dict(
    :p => 1u"Pa",
    :rho => 1u"kg/m^3",
    :T => 1u"K",
    :a => 1u"m/s",
    :gamma => 1.1,
    :p0 => 2u"Pa",
    :rho0 => 2u"kg/m^3",
    :T0 => 2u"K",
    :a0 => 2u"m/s",
    :mdot => 1u"kg/s",
    :A => 1u"m^2",
    :v => 1u"m/s",
    :M => 1
)
##
residues(;params...)
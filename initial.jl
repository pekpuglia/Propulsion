#https://docs.juliahub.com/Flatten/hpRkL/0.4.0/ - útil?
using Unitful, NonlinearSolve
##
abstract type PhysicalProperties end

dof(::Type{<:PhysicalProperties}) = error("PhysicalProperties types must implement dof")

function variables(T::Type{<:PhysicalProperties})
    vars = fieldnames(T)
    types = fieldtypes(T)

    indexes_to_recurse = findall(types .<: PhysicalProperties)

    cat(((i in indexes_to_recurse) ? variables(types[i]) : var for (i, var) in enumerate(vars))..., dims=1)
end

residues(::T) where T <: PhysicalProperties = error("PhysicalProperties types must implement residues")
##
struct ThermodynamicProperties <: PhysicalProperties
    P
    z
    T
end

function ThermodynamicProperties(; P, z, T)
    ThermodynamicProperties(P, z, T)
end

#garantir length(residues) + dof = length(fieldnames)?
dof(::Type{ThermodynamicProperties}) = 2

const Rmolar = 8.3144598u"J/mol/K"

function residues(tp::ThermodynamicProperties)
    [tp.P - tp.z*tp.T*ustrip(u"J/mol/K", Rmolar)]
end
##
struct MassProperties <: PhysicalProperties
    tp::ThermodynamicProperties
    MM
    rho
    R
end

function MassProperties(; MM, rho, R, kwargs...)
    MassProperties(ThermodynamicProperties(;kwargs...), MM, rho, R)
end

dof(::Type{MassProperties}) = dof(ThermodynamicProperties) + 1

function residues(mp::MassProperties)
    [
        residues(mp.tp)
        mp.R * mp.MM - ustrip(u"J/mol/K", Rmolar)
        mp.rho - mp.MM * mp.tp.z
    ]
end
##
struct CalorificProperties <: PhysicalProperties
    mp::MassProperties
    cv
    cp
    gamma
    a
end

##
struct FlowProperties <: PhysicalProperties
    cp::CalorificProperties
    M
    v
    T0
    rho0
    P0
    a0
end
##
struct Quasi1dimflowProperties <: PhysicalProperties
    fp::FlowProperties
    mdot
    A
    Astar
end

##
function solve_params(T::Type; kwargs...)
    allvars = variables(T)
    if length(kwargs) != dof(T)
        error("$(dof(T)) thermodynamic properties needed, $(length(kwargs)) given: $(keys(kwargs))")
    end
    if any([!(k in allvars) for k in keys(kwargs)])
        error("expected keys from $allvars, got: $(keys(kwargs)) ")
    end
    #https://www.geeksforgeeks.org/sets-in-julia/
    missingvars = (setdiff(Set(allvars), Set(keys(kwargs))) |> collect)

    prob = NonlinearProblem(
        (values, p) -> residues(T(;
            Dict(
                Dict(missingvar => value for (missingvar, value) in zip(missingvars, values))..., 
                kwargs...
            )...)), ones(size(missingvars)), p=())
    
            sol = solve(prob, NewtonRaphson())

    T(;
        Dict(
            Dict(missingvar => u for (missingvar, u) in zip(missingvars, sol.u))...,
            Dict(k => Float64(v) for (k, v) in kwargs)...
        )...)
end
##
solve_params(ThermodynamicProperties, P= 1.0, T = 10.0)
##
MassProperties(; MM = 1, rho = 2, R = 3, P = 1.0, T = 10.0, z= 3.0)
##
solve_params(MassProperties, P=1.0, T=10.0, rho = 2.0)
##
#assume que kwargs tá completo
# function residues(;kwargs...)
#     p = kwargs[:p]
#     rho = kwargs[:rho]
#     T = kwargs[:T]
#     a = kwargs[:a]
#     gamma = kwargs[:gamma]
    
#     p0 = kwargs[:p0]
#     rho0 = kwargs[:rho0]
#     T0 = kwargs[:T0]
#     a0 = kwargs[:a0]

#     mdot = kwargs[:mdot]
#     A = kwargs[:A]
#     v = kwargs[:v]
#     M = kwargs[:M]
    
    

#     [
#         1 + (gamma - 1) / 2 * M^2 - T0/T
#         (1 + (gamma -  1) / 2 * M^2)^(gamma/(gamma-1)) - p0/p
#         (1 + (gamma -  1) / 2 * M^2)^(  1  /(gamma-1)) - rho0/rho
#         M*a - v
#         a^2 - gamma*p/rho
#         a^2/(gamma-1) + v^2/2 - a0^2/(gamma-1)
#         mdot - rho*v*A
#     ]
# end
# ##
# params = Dict(
#     :p => 1u"Pa",
#     :rho => 1u"kg/m^3",
#     :T => 1u"K",
#     :a => 1u"m/s",
#     :gamma => 1.1,
#     :p0 => 2u"Pa",
#     :rho0 => 2u"kg/m^3",
#     :T0 => 2u"K",
#     :a0 => 2u"m/s",
#     :mdot => 1u"kg/s",
#     :A => 1u"m^2",
#     :v => 1u"m/s",
#     :M => 1
# )
# ##
# residues(;params...)
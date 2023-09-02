#https://docs.juliahub.com/Flatten/hpRkL/0.4.0/ - útil?
using Unitful, NonlinearSolve
##
abstract type PhysicalProperties end

function (T::Type{<:PhysicalProperties})(;kwargs...)
    vars = fieldnames(T)
    types = fieldtypes(T)

    indexes_to_recurse = findall(types .<: PhysicalProperties)

    parameters = cat(((i in indexes_to_recurse) ? types[i](;kwargs...) : kwargs[var] for (i, var) in enumerate(vars))..., dims=1)
    
    T(parameters...)
end

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

#garantir length(residues) + dof = length(fieldnames)?
dof(::Type{ThermodynamicProperties}) = 2

const Rmolar = 8.3144598u"J/mol/K"

unit_adapter(::Type, val) = ustrip(val)
unit_adapter(T::Type{<:Quantity}, val) = uconvert(unit(T), val)

function residues(tp::ThermodynamicProperties)
    [tp.P - unit_adapter(typeof(tp.P), tp.z*tp.T*Rmolar)]
end
##
struct MassProperties <: PhysicalProperties
    tp::ThermodynamicProperties
    MM
    rho
    R
end

dof(::Type{MassProperties}) = dof(ThermodynamicProperties) + 1

function residues(mp::MassProperties)
    [
        residues(mp.tp)
        mp.R * mp.MM - unit_adapter(typeof(mp.R*mp.MM), Rmolar)
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

dof(::Type{CalorificProperties}) = dof(MassProperties) + 1

function residues(cp::CalorificProperties)
    [
        residues(cp.mp)
        #cp = cv + R
        cp.cv + cp.mp.R - cp.cp
        #gamma = cp/cv
        cp.gamma * cp.cv - cp.cp
        #a = sqrt(gamma R T)
        cp.a^2 - cp.gamma * cp.mp.R * cp.mp.tp.T
    ]
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

dof(::Type{FlowProperties}) = dof(CalorificProperties) + 1

function residues(fp::FlowProperties)
    [
        residues(fp.cp)
        #M*a - v
        fp.M * fp.cp.a - fp.v
        #1 + (gamma - 1) / 2 * M^2 - T0/T
        1 + (fp.cp.gamma - 1) / 2 * fp.M^2 - fp.T0 / fp.cp.mp.tp.T
        #(1 + (gamma -  1) / 2 * M^2)^(gamma/(gamma-1)) - p0/p
        (1 + (fp.cp.gamma -1) / 2 * fp.M^2) ^ (fp.cp.gamma / (fp.cp.gamma - 1)) - fp.P0 / fp.cp.mp.tp.P
        #(1 + (gamma -  1) / 2 * M^2)^(  1  /(gamma-1)) - rho0/rho
        (1 + (fp.cp.gamma - 1) / 2 * fp.M^2) ^ (1 / (fp.cp.gamma - 1)) - fp.rho0 / fp.cp.mp.rho
        #a^2/(gamma-1) + v^2/2 - a0^2/(gamma-1)
        fp.cp.a^2 / (fp.cp.gamma - 1) + fp.v^2 / 2 - fp.a0^2 / (fp.cp.gamma - 1)
    ]
end
##
struct Quasi1dimflowProperties <: PhysicalProperties
    fp::FlowProperties
    mdot
    A
    Astar
end

dof(::Type{Quasi1dimflowProperties}) = dof(FlowProperties) + 1

function residues(qp::Quasi1dimflowProperties)
    [
        residues(qp.fp)
        qp.mdot - qp.A * qp.fp.v * qp.fp.cp.mp.rho
        #p 682 anderson
        (qp.A / qp.Astar) ^ 2 - 
            1 / qp.fp.M^2 * (
                2 / (qp.fp.cp.gamma + 1) * (
                    1 + (qp.fp.cp.gamma - 1) / 2 * qp.fp.M ^ 2
            )) ^ ((qp.fp.cp.gamma + 1) / (qp.fp.cp.gamma - 1))
    ]
end

##
#adicionar testes, gerenciar unidades, refatorar essa função
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
            )...)), 
        ones(size(missingvars)), p=()
    )
    
    sol = solve(prob, NewtonRaphson())

    T(;
        Dict(
            Dict(missingvar => u for (missingvar, u) in zip(missingvars, sol.u))...,
            Dict(k => Float64(v) for (k, v) in kwargs)...
        )...)
end
##
ThermodynamicProperties(P = 1, T = 10.0, z= 3.0)
##
solve_params(ThermodynamicProperties, P= 1.0, T = 10.0)
##
MassProperties(; MM = 1, rho = 2, R = 3, P = 1.0, T = 10.0, z= 3.0)
##
#testar que sistema não está super/sub restringido?
solve_params(MassProperties, P=1.0, MM=10.0, rho = 2.0)
##
CalorificProperties(; MM = 1, rho = 2, R = 3, P = 1.0, T = 10.0, z= 3.0, cv= 1.0, cp=1.0, gamma = 1.0, a = 3.0)
##
solve_params(CalorificProperties, P=1e5, T=10.0, rho = 2.0, gamma = 1.4)
##
solve_params(FlowProperties, P=1e5, T=10.0, rho = 2.0, gamma = 1.4, M = 1.5)
##
q1dparams = solve_params(Quasi1dimflowProperties, P=1e5, T=10.0, rho = 2.0, gamma = 1.4, M = 1.5, A = 1.0)
##
residues(Quasi1dimflowProperties(
    P = (q1dparams.fp.cp.mp.tp.P)u"Pa",
    T = (q1dparams.fp.cp.mp.tp.T)u"K",
    z = (q1dparams.fp.cp.mp.tp.z)u"mol/m^3",
    rho = (q1dparams.fp.cp.mp.rho)u"kg/m^3",
    MM = (q1dparams.fp.cp.mp.MM)u"kg/mol",
    R = (q1dparams.fp.cp.mp.R)u"J/kg/K",
    a = (q1dparams.fp.cp.a)u"m/s",
    cp = (q1dparams.fp.cp.cp)u"J/kg/K",
    cv = (q1dparams.fp.cp.cv)u"J/kg/K",
    gamma = (q1dparams.fp.cp.gamma),
    M = (q1dparams.fp.M),
    a0 = (q1dparams.fp.a0)u"m/s",
    P0 = (q1dparams.fp.P0)u"Pa",
    rho0 = (q1dparams.fp.rho0)u"kg/m^3",
    T0 = (q1dparams.fp.T0)u"K",
    v = (q1dparams.fp.v)u"m/s",
    A = (q1dparams.A)u"m^2",
    Astar = (q1dparams.Astar)u"m^2",
    mdot = (q1dparams.mdot)u"kg/s"
))
##
unit_thermo_props = ThermodynamicProperties(1u"Pa", 1u"mol/m^3", 1u"K")
plain_thermo_props = ThermodynamicProperties(1.0, 1, 1.0)
##
should_fail = ThermodynamicProperties(1, 1.0, 1.0u"K")
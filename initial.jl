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

# units(::Type) = error("PhysicalProperties types must implement units")

#melhorar
units(T::Type{<:PhysicalProperties}) = Dict(var => NoUnits for var in variables(T))

residues(::T) where T <: PhysicalProperties = error("PhysicalProperties types must implement residues")
##
abstract type AbstractUnitMarker end

struct WithUnits <: AbstractUnitMarker end
struct WithOutUnits <: AbstractUnitMarker end

struct ThermodynamicProperties{UnitMarker <: AbstractUnitMarker} <: PhysicalProperties
    P
    z
    T
    function ThermodynamicProperties(P::Real, z::Real, T::Real)
        new{WithOutUnits}(P, z, T)
    end
    function ThermodynamicProperties(P::Unitful.Pressure, z::Unitful.Molarity, T::Unitful.Temperature)
        new{WithUnits}(P, z, T)
    end
end


#garantir length(residues) + dof = length(fieldnames)?
dof(::Type{<:ThermodynamicProperties}) = 2

units(::Type{ThermodynamicProperties{WithUnits}}) = Dict(
    :P => u"Pa",
    :z => u"mol/m^3",
    :T => u"K"
)

# units(T::Type{ThermodynamicProperties{WithOutUnits}}) = Dict(var => NoUnits for var in variables(T))

const Rmolar = 8.3144598u"J/mol/K"

r_molar(::Type{WithUnits}) = Rmolar
r_molar(::Type{WithOutUnits}) = ustrip(u"J/mol/K", Rmolar)

function residues(tp::ThermodynamicProperties{T}) where T
    [tp.P - tp.z*tp.T*r_molar(T)]
end
##
@derived_dimension MolarMass Unitful.𝐌/Unitful.𝐍 true
@derived_dimension SpecificHeatCapacity Unitful.𝐋^2 * Unitful.𝐓^-2 /Unitful.𝚯 true
struct MassProperties{UnitMarker <: AbstractUnitMarker} <: PhysicalProperties
    tp::ThermodynamicProperties{UnitMarker}
    MM
    rho
    R
    function MassProperties(tp::ThermodynamicProperties{WithOutUnits}, MM::Real, rho::Real, R::Real)
        new{WithOutUnits}(tp, MM, rho, R)
    end
    function MassProperties(tp::ThermodynamicProperties{WithUnits}, MM::MolarMass, rho::Unitful.Density, R::SpecificHeatCapacity)
        new{WithUnits}(tp, MM, rho, R)
    end
end

dof(::Type{MassProperties}) = dof(ThermodynamicProperties) + 1

units(::Type{MassProperties{WithUnits}}) = Dict(
    :MM => u"kg/mol",
    :rho => u"kg/m^3",
    :R => u"J/kg/K",
    units(ThermodynamicProperties{WithUnits})...
)

# units(T::Type{MassProperties{WithOutUnits}}) = Dict(
#     var => NoUnits for var in variables(T)
# )

function residues(mp::MassProperties{T}) where T
    [
        residues(mp.tp)
        mp.R * mp.MM - r_molar(T)
        mp.rho - mp.MM * mp.tp.z
    ]
end
##
struct CalorificProperties{UnitMarker <: AbstractUnitMarker} <: PhysicalProperties
    mp::MassProperties{UnitMarker}
    cv
    cp
    gamma
    a
    function CalorificProperties(mp::MassProperties{WithOutUnits}, cv::Real, cp::Real, gamma::Real, a::Real)
        new{WithOutUnits}(mp, cv, cp, gamma, a)
    end
    function CalorificProperties(mp::MassProperties{WithUnits}, cv::SpecificHeatCapacity, cp::SpecificHeatCapacity, gamma::Real, a::Unitful.Velocity)
        new{WithUnits}(mp, cv, cp, gamma, a)
    end
end

dof(::Type{CalorificProperties}) = dof(MassProperties) + 1

units(::Type{CalorificProperties{WithUnits}}) = Dict(
    :cv => u"J/kg/K",
    :cp => u"J/kg/K",
    :gamma => NoUnits,
    :a => u"m/s",
    units(MassProperties{WithUnits})...
)

# units(T::Type{CalorificProperties{WithOutUnits}}) = Dict(
#     var => NoUnits for var in variables(T)
# )

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
struct FlowProperties{UnitMarker <: AbstractUnitMarker} <: PhysicalProperties
    cp::CalorificProperties{UnitMarker}
    M
    v
    T0
    rho0
    P0
    a0
    function FlowProperties(cp::CalorificProperties{WithOutUnits}, 
        M::Real, v::Real, T0::Real, rho0::Real, P0::Real, a0::Real)
        new{WithOutUnits}(cp, M, v, T0, rho0, P0, a0)
    end
    function FlowProperties(cp::CalorificProperties{WithUnits}, 
        M::Real, v::Unitful.Velocity, T0::Unitful.Temperature, 
        rho0::Unitful.Density, P0::Unitful.Pressure, a0::Unitful.Velocity)
        new{WithUnits}(cp, M, v, T0, rho0, P0, a0)
    end
end

dof(::Type{FlowProperties}) = dof(CalorificProperties) + 1

units(::Type{FlowProperties{WithUnits}}) = Dict(
    :M => NoUnits,
    :v => u"m/s",
    :T0 => u"K",
    :rho0 => u"kg/m^3",
    :P0 => u"Pa",
    :a0 => u"m/s",
    units(CalorificProperties{WithUnits})...
)

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
struct Quasi1dimflowProperties{UnitMarker <: AbstractUnitMarker} <: PhysicalProperties
    fp::FlowProperties{UnitMarker}
    mdot
    A
    Astar
    function Quasi1dimflowProperties(fp::FlowProperties{WithOutUnits}, mdot::Real, A::Real, Astar::Real)
        new{WithOutUnits}(fp, mdot, A, Astar)
    end
    function Quasi1dimflowProperties(fp::FlowProperties{WithUnits}, mdot::Unitful.MassFlow, 
        A::Unitful.Area, Astar::Unitful.Area)
        new{WithUnits}(fp, mdot, A, Astar)
    end
end

dof(::Type{Quasi1dimflowProperties}) = dof(FlowProperties) + 1

units(::Type{Quasi1dimflowProperties{WithUnits}}) = Dict(
    :mdot => u"kg/s",
    :A => u"m^2",
    :Astar => u"m^2",
    units(FlowProperties{WithUnits})...
)

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
#todo
#abstrair construtores
#adicionar testes, refatorar essa função
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

    #decide if solver should use units or not
    if all(typeof.(values(kwargs) |> collect) .<: Real)
        allunits = units(T{WithOutUnits})
    else
        allunits = units(T{WithUnits})
    end

    unitless_kwargs = Dict(key => ustrip(allunits[key], val) for (key, val) in kwargs)

    prob = NonlinearProblem(
        (values, p) -> residues(T(;
            Dict(
                Dict(missingvar => value for (missingvar, value) in zip(missingvars, values))..., 
                unitless_kwargs...
            )...)), 
        ones(size(missingvars)), p=()
    )
    
    sol = solve(prob, NewtonRaphson())

    #bind units back to variables
    T(;
        Dict(
            Dict(missingvar => u*allunits[missingvar] for (missingvar, u) in zip(missingvars, sol.u))...,
            Dict(k => Float64(v)*allunits[k] for (k, v) in unitless_kwargs)...
        )...)
end
##
# ThermodynamicProperties(P = 1, T = 10.0, z= 3.0)
# ##
# solve_params(ThermodynamicProperties, P= 1.0, T = 10.0)
# ##
# solve_params(ThermodynamicProperties, P= 1.0u"Pa", T = 10.0u"K")
# ##
# MassProperties(; MM = 1, rho = 2, R = 3, P = 1.0, T = 10.0, z= 3.0)
# ##
# MassProperties(; MM = 1u"g/mol", rho = 2u"kg/cm^3", R = 3u"kJ/kg/K", P = 1.0u"Pa", T = 10.0u"K", z= 3.0u"mol/mm^3")
# ##
# #testar que sistema não está super/sub restringido?
# solve_params(MassProperties, P=1.0, MM=10.0, rho = 2.0)
# ##
# solve_params(MassProperties, P=1.0u"Pa", MM=10.0u"g/mol", rho = 2.0u"kg/m^3")
# ##
# CalorificProperties(; MM = 1, rho = 2, R = 3, P = 1.0, T = 10.0, z= 3.0, cv= 1.0, cp=1.0, gamma = 1.0, a = 3.0)
# ##
# solve_params(CalorificProperties, P=1e5, T=10.0, rho = 2.0, gamma = 1.4)
# ##
# solve_params(FlowProperties, P=1e5, T=10.0, rho = 2.0, gamma = 1.4, M = 1.5)
# ##
# q1dparams = solve_params(Quasi1dimflowProperties, P=1e5, T=10.0, rho = 2.0, gamma = 1.4, M = 1.5, A = 1.0)
# ##
# residues(Quasi1dimflowProperties(
#     P = (q1dparams.fp.cp.mp.tp.P)u"Pa",
#     T = (q1dparams.fp.cp.mp.tp.T)u"K",
#     z = (q1dparams.fp.cp.mp.tp.z)u"mol/m^3",
#     rho = (q1dparams.fp.cp.mp.rho)u"kg/m^3",
#     MM = (q1dparams.fp.cp.mp.MM)u"kg/mol",
#     R = (q1dparams.fp.cp.mp.R)u"J/kg/K",
#     a = (q1dparams.fp.cp.a)u"m/s",
#     cp = (q1dparams.fp.cp.cp)u"J/kg/K",
#     cv = (q1dparams.fp.cp.cv)u"J/kg/K",
#     gamma = (q1dparams.fp.cp.gamma),
#     M = (q1dparams.fp.M),
#     a0 = (q1dparams.fp.a0)u"m/s",
#     P0 = (q1dparams.fp.P0)u"Pa",
#     rho0 = (q1dparams.fp.rho0)u"kg/m^3",
#     T0 = (q1dparams.fp.T0)u"K",
#     v = (q1dparams.fp.v)u"m/s",
#     A = (q1dparams.A)u"m^2",
#     Astar = (q1dparams.Astar)u"m^2",
#     mdot = (q1dparams.mdot)u"kg/s"
# ))
# ##
# unit_thermo_props = ThermodynamicProperties(1u"Pa", 1u"mol/m^3", 1u"K")
# plain_thermo_props = ThermodynamicProperties(1.0, 1, 1.0)
# ##
# should_fail = ThermodynamicProperties(1u"Pa", 1.0, 1.0)
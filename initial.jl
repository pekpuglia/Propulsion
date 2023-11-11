#todo
#pipe flow
#participation matrix
#sub/supersonic choice - initial guesses depend on the variables
#more tests
#normal shocks?
#actually use this - do As/Ac x Pt/Ps plot

#https://docs.juliahub.com/Flatten/hpRkL/0.4.0/ - útil?
using Revise, Unitful, NonlinearSolve, Symbolics
##
abstract type PhysicalProperties end

function phys_prop_from_kwargs(T::Type{<:PhysicalProperties};kwargs...)
    vars = fieldnames(T)
    types = fieldtypes(T)

    indexes_to_recurse = findall(types .<: PhysicalProperties)

    parameters = cat(((i in indexes_to_recurse) ? phys_prop_from_kwargs(types[i];kwargs...) : kwargs[var] for (i, var) in enumerate(vars))..., dims=1)
    T(parameters...)
end

function add_units(pp::T, unit_dict) where T <: PhysicalProperties
    vars = fieldnames(T)
    types = fieldtypes(T)
    indexes_to_recurse = findall(types .<: PhysicalProperties)
    parameters = cat(((i in indexes_to_recurse) ? add_units(getfield(pp, var), unit_dict) : getfield(pp, var) * unit_dict[var] for (i, var) in enumerate(vars))..., dims=1)

    T(parameters...)
end


dof(::Type{<:PhysicalProperties}) = error("PhysicalProperties types must implement dof")

function variables(T::Type{<:PhysicalProperties})
    vars = fieldnames(T)
    types = fieldtypes(T)

    indexes_to_recurse = findall(types .<: PhysicalProperties)

    cat(((i in indexes_to_recurse) ? variables(types[i]) : var for (i, var) in enumerate(vars))..., dims=1)
end

Base.propertynames(::T) where T <: PhysicalProperties = variables(T)

##ONLY WORKS FOR SINGLE RECURSION CASE
function Base.getproperty(pp::T, s::Symbol) where T <: PhysicalProperties
    vars = fieldnames(T)
    types = fieldtypes(T)

    index_to_recurse = findfirst(types .<: PhysicalProperties)

    (s in vars) ? getfield(pp, s) : getproperty(getfield(pp, vars[index_to_recurse]), s)
end

#melhorar
units(T::Type{<:PhysicalProperties}) = Dict(var => NoUnits for var in variables(T))

residues(::T) where T <: PhysicalProperties = error("PhysicalProperties types must implement residues")

## symbolic analysis tools
function generate_sym_var_dict(T::Type{<:PhysicalProperties})
    Dict(var => (@variables $var)[1] for var in variables(T))
end

function participation_vector(T::Type{<:PhysicalProperties})
    sym_var_dict = generate_sym_var_dict(T)
    part_vector_symbolic = Symbolics.get_variables.(residues(phys_prop_from_kwargs(T; sym_var_dict...)))
    
    rev_svd = Dict(values(sym_var_dict) .=> keys(sym_var_dict))

    map(v -> getindex.([rev_svd], v), part_vector_symbolic)
end
##

struct ThermodynamicProperties <: PhysicalProperties
    P
    z
    T
    function ThermodynamicProperties(P::Real, z::Real, T::Real)
        new(P, z, T)
    end
    function ThermodynamicProperties(P::Unitful.Pressure, z::Unitful.Molarity, T::Unitful.Temperature)
        new(P, z, T)
    end
end

dof(::Type{<:ThermodynamicProperties}) = 2

units(::Type{ThermodynamicProperties}) = Dict(
    :P => u"Pa",
    :z => u"mol/m^3",
    :T => u"K"
)

const Rmolar = 8.3144598u"J/mol/K"

r_molar(::Number) = Rmolar
r_molar(::Real) = ustrip(u"J/mol/K", Rmolar)

function residues(tp::ThermodynamicProperties)
    [tp.P - tp.z*tp.T*r_molar(tp.P)]
end
##
@derived_dimension MolarMass Unitful.𝐌/Unitful.𝐍 true
@derived_dimension SpecificHeatCapacity Unitful.𝐋^2 * Unitful.𝐓^-2 /Unitful.𝚯 true
struct MassProperties <: PhysicalProperties
    tp::ThermodynamicProperties
    MM
    rho
    R
    function MassProperties(tp::ThermodynamicProperties, MM::Real, rho::Real, R::Real)
        new(tp, MM, rho, R)
    end
    function MassProperties(tp::ThermodynamicProperties, MM::MolarMass, rho::Unitful.Density, R::SpecificHeatCapacity)
        new(tp, MM, rho, R)
    end
end

dof(::Type{MassProperties}) = dof(ThermodynamicProperties) + 1

units(::Type{MassProperties}) = Dict(
    :MM => u"kg/mol",
    :rho => u"kg/m^3",
    :R => u"J/kg/K",
    units(ThermodynamicProperties)...
)

function residues(mp::MassProperties)
    [
        residues(mp.tp)
        mp.R * mp.MM - r_molar(mp.MM)
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
    function CalorificProperties(mp::MassProperties, cv::Real, cp::Real, gamma::Real, a::Real)
        new(mp, cv, cp, gamma, a)
    end
    function CalorificProperties(mp::MassProperties, cv::SpecificHeatCapacity, cp::SpecificHeatCapacity, gamma::Real, a::Unitful.Velocity)
        new(mp, cv, cp, gamma, a)
    end
end

dof(::Type{CalorificProperties}) = dof(MassProperties) + 1

units(::Type{CalorificProperties}) = Dict(
    :cv => u"J/kg/K",
    :cp => u"J/kg/K",
    :gamma => NoUnits,
    :a => u"m/s",
    units(MassProperties)...
)

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
    function FlowProperties(cp::CalorificProperties, 
        M::Real, v::Real, T0::Real, rho0::Real, P0::Real, a0::Real)
        new(cp, M, v, T0, rho0, P0, a0)
    end
    function FlowProperties(cp::CalorificProperties, 
        M::Real, v::Unitful.Velocity, T0::Unitful.Temperature, 
        rho0::Unitful.Density, P0::Unitful.Pressure, a0::Unitful.Velocity)
        new(cp, M, v, T0, rho0, P0, a0)
    end
end

dof(::Type{FlowProperties}) = dof(CalorificProperties) + 1

units(::Type{FlowProperties}) = Dict(
    :M => NoUnits,
    :v => u"m/s",
    :T0 => u"K",
    :rho0 => u"kg/m^3",
    :P0 => u"Pa",
    :a0 => u"m/s",
    units(CalorificProperties)...
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
struct Quasi1dimflowProperties <: PhysicalProperties
    fp::FlowProperties
    mdot
    A
    Astar
    function Quasi1dimflowProperties(fp::FlowProperties, mdot::Real, A::Real, Astar::Real)
        new(fp, mdot, A, Astar)
    end
    function Quasi1dimflowProperties(fp::FlowProperties, mdot::Unitful.MassFlow, 
        A::Unitful.Area, Astar::Unitful.Area)
        new(fp, mdot, A, Astar)
    end
end

dof(::Type{Quasi1dimflowProperties}) = dof(FlowProperties) + 1

units(::Type{Quasi1dimflowProperties}) = Dict(
    :mdot => u"kg/s",
    :A => u"m^2",
    :Astar => u"m^2",
    units(FlowProperties)...
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
function internal_solver(T::Type, input_data::Dict{Symbol, <:Real})
    allvars = variables(T)

    missingvars = (setdiff(Set(allvars), Set(keys(input_data))) |> collect)
    
    prob = NonlinearProblem(
        (values, p) -> residues(phys_prop_from_kwargs(T;
            Dict(
                Dict(missingvar => value for (missingvar, value) in zip(missingvars, values))..., 
                input_data...
            )...)), 
        ones(size(missingvars)), p=()
    )

    sol = solve(prob, NewtonRaphson())

    phys_prop_from_kwargs(T;
        Dict(
            Dict(missingvar => u for (missingvar, u) in zip(missingvars, sol.u))...,
            Dict(k => Float64(v) for (k, v) in input_data)...
        )...)
end

function internal_solver(T::Type, input_data::Dict{Symbol, <:Number})
    internal_units = units(T)

    unitless_kwargs = Dict(key => ustrip(internal_units[key], val) for (key, val) in input_data)

    unitless_solution = internal_solver(T, unitless_kwargs)

    add_units(unitless_solution, internal_units)
end
##
function (T::Type)(; kwargs...)
    allvars = variables(T)

    if any([!(k in allvars) for k in keys(kwargs)])
        error("expected keys from $allvars, got: $(keys(kwargs)) ")
    end

    if length(kwargs) != dof(T)
        error("$(dof(T)) thermodynamic properties needed, $(length(kwargs)) given: $(keys(kwargs))")
    end
    
    

    internal_solver(T, Dict(kwargs...))
    
end
# q1dparams = Quasi1dimflowProperties(P=1e5, T=10.0, rho = 2.0, gamma = 1.4, Astar = 0.85, A = 1.0)
##
#don't work???
#Quasi1dimflowProperties(P = 1u"bar", R = 287u"J/kg/K", gamma = 1.4, mdot = 1u"kg/s", T = 300u"K", M=1.5)
#Quasi1dimflowProperties(P = 1u"bar", R = 287u"J/kg/K", gamma = 1.4, mdot = 1u"kg/s", T = 300u"K", A=0.3u"m^2")
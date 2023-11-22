export ThermodynamicProperties
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
export MassProperties

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
export CalorificProperties
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
        cp.a - √(cp.gamma * cp.mp.R * cp.mp.tp.T)
    ]
end
##
export FlowProperties
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
        (1 + (fp.cp.gamma - 1) / 2 * fp.M^2) * fp.cp.mp.tp.T - fp.T0
        #(1 + (gamma -  1) / 2 * M^2)^(gamma/(gamma-1)) - p0/p
        ((1 + (fp.cp.gamma - 1) / 2 * fp.M^2) ^ (fp.cp.gamma / (fp.cp.gamma - 1))) * fp.cp.mp.tp.P - fp.P0
        #(1 + (gamma -  1) / 2 * M^2)^(  1  /(gamma-1)) - rho0/rho
        (1 + (fp.cp.gamma - 1) / 2 * fp.M^2) ^ (1 / (fp.cp.gamma - 1)) * fp.cp.mp.rho - fp.rho0
        #a^2/(gamma-1) + v^2/2 - a0^2/(gamma-1)
        fp.cp.a^2 / (fp.cp.gamma - 1) + fp.v^2 / 2 - fp.a0^2 / (fp.cp.gamma - 1)
    ]
end
##
export Quasi1dimflowProperties
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
        qp.A * qp.fp.M -
            qp.Astar * (
                2 / (qp.fp.cp.gamma + 1) * (
                    1 + (qp.fp.cp.gamma - 1) / 2 * qp.fp.M ^ 2
            )) ^ ((qp.fp.cp.gamma + 1) / 2(qp.fp.cp.gamma - 1))
    ]
end

#expand for N sections
#refactor
#test solver !!!!! analyze coherence of input, residues, etc
export NozzleFlowProperties
struct NozzleFlowProperties <: PhysicalProperties
    sec1::Quasi1dimflowProperties
    sec2::Quasi1dimflowProperties
    F
    function NozzleFlowProperties(sec1::Quasi1dimflowProperties, sec2::Quasi1dimflowProperties, F::Real)
        new(sec1, sec2, F)
    end
    function NozzleFlowProperties(sec1::Quasi1dimflowProperties, sec2::Quasi1dimflowProperties, F::Unitful.Force)
        new(sec1, sec2, F)
    end
end

dof(::Type{NozzleFlowProperties}) = dof(Quasi1dimflowProperties) + 1

function variables(T::Type{NozzleFlowProperties})
    [
        variables(Quasi1dimflowProperties) .|> string .|> (str -> str*"_1") .|> Symbol
        variables(Quasi1dimflowProperties) .|> string .|> (str -> str*"_2") .|> Symbol
        :F
    ]
end

units(::Type{NozzleFlowProperties}) = Dict(
    Dict(
        Symbol(string(pair.first)*suff) => pair.second 
        for pair in Propulsion.units(Quasi1dimflowProperties) 
        for suff in ["_1", "_2"]
    )...,
    :F => u"N"
)

function residues(nfp::NozzleFlowProperties)
    [
        residues(nfp.sec1)
        residues(nfp.sec2)
        nfp.sec1.fp.cp.gamma - nfp.sec2.fp.cp.gamma
        nfp.sec1.R - nfp.sec2.R
        nfp.sec1.Astar - nfp.sec2.Astar
        nfp.sec1.mdot - nfp.sec2.mdot
        nfp.sec1.T0 - nfp.sec2.T0
        nfp.F - (
            nfp.sec2.mdot * nfp.sec2.v - nfp.sec2.P * nfp.sec2.A
        ) + (
            nfp.sec1.mdot * nfp.sec1.v - nfp.sec1.P * nfp.sec1.A
        )
    ]
end

Base.getproperty(nfp::NozzleFlowProperties, s::Symbol) = getfield(nfp, s)

function Base.getindex(nfp::NozzleFlowProperties, i)
    if !(i == 1 || i == 2)
        throw(BoundsError(nfp, i))
    end

    (i == 1) ? nfp.sec1 : nfp.sec2
end

function select_and_remove_kwarg_suffix(suff::String, kwargs)
    Dict(
        Symbol(string(pair.first)[1:(end-2)]) => pair.second
        for pair in filter(x -> endswith(string(x.first), suff), kwargs)
    )
end

function phys_prop_from_kwargs(T::Type{NozzleFlowProperties}; kwargs...)
    
    kw1 = select_and_remove_kwarg_suffix("_1", kwargs)
    kw2 = select_and_remove_kwarg_suffix("_2", kwargs)

    NozzleFlowProperties(
        phys_prop_from_kwargs(Quasi1dimflowProperties; kw1...),
        phys_prop_from_kwargs(Quasi1dimflowProperties; kw2...),
        kwargs[:F]
    )
end

default_initial_guesses(::Type{NozzleFlowProperties}) = Dict(
    :gamma_2 => 1.4
)
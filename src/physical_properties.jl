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
        cp.a^2 - (cp.gamma * cp.mp.R * cp.mp.tp.T)
    ]
end
##
export FlowProperties
struct FlowProperties <: PhysicalProperties
    cal_prop::CalorificProperties
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
        residues(fp.cal_prop)
        #M*a - v
        fp.M * fp.cal_prop.a - fp.v
        #1 + (gamma - 1) / 2 * M^2 - T0/T
        (1 + (fp.cal_prop.gamma - 1) / 2 * fp.M^2) * fp.cal_prop.mp.tp.T - fp.T0
        #(1 + (gamma -  1) / 2 * M^2)^(gamma/(gamma-1)) - p0/p
        ((1 + (fp.cal_prop.gamma - 1) / 2 * fp.M^2) ^ (fp.cal_prop.gamma / (fp.cal_prop.gamma - 1))) * fp.cal_prop.mp.tp.P - fp.P0
        #(1 + (gamma -  1) / 2 * M^2)^(  1  /(gamma-1)) - rho0/rho
        (1 + (fp.cal_prop.gamma - 1) / 2 * fp.M^2) ^ (1 / (fp.cal_prop.gamma - 1)) * fp.cal_prop.mp.rho - fp.rho0
        #a^2/(gamma-1) + v^2/2 - a0^2/(gamma-1)
        fp.cal_prop.a^2 / (fp.cal_prop.gamma - 1) + fp.v^2 / 2 - fp.a0^2 / (fp.cal_prop.gamma - 1)
    ]
end
##
#gives negative answers sometimes
export NormalShockProperties
struct NormalShockProperties <: PhysicalProperties
    fp1::FlowProperties
    fp2::FlowProperties
end

function variables(::Type{NormalShockProperties})
    [
        variables(FlowProperties) .|> string .|> (str -> str*"_1") .|> Symbol
        variables(FlowProperties) .|> string .|> (str -> str*"_2") .|> Symbol
    ]
end

units(T::Type{NormalShockProperties}) = Dict(
    Symbol(string(v)*suff) => u for (v, u) in units(FlowProperties) for suff in ["_1", "_2"]
)

function residues(nsp::NormalShockProperties)
    [
        residues(nsp.fp1)
        residues(nsp.fp2)
        nsp.fp1.gamma - nsp.fp2.gamma
        nsp.fp1.R - nsp.fp2.R
        (nsp.fp1.gamma * nsp.fp1.M^2 - (nsp.fp1.gamma - 1) / 2) * (nsp.fp2.gamma * nsp.fp2.M^2 - (nsp.fp2.gamma - 1) / 2) - ((nsp.fp1.gamma + 1)/2)^2
        # nsp.fp2.M - sqrt((1 + (nsp.fp1.gamma - 1)/2 * nsp.fp1.M^2) / (nsp.fp1.gamma * nsp.fp1.M^2 - (nsp.fp1.gamma - 1)/2))
        nsp.fp1.T0 - nsp.fp2.T0
        nsp.fp2.P - nsp.fp1.P * (1 + 2nsp.fp1.gamma*(nsp.fp1.M^2 - 1) / (nsp.fp1.gamma + 1))
        # nsp.fp1.M - sqrt((nsp.fp2.P / nsp.fp1.P - 1) * (nsp.fp1.gamma + 1) / (2*nsp.fp1.gamma) + 1)
    ]
end

function NormalShockProperties(data_dict::Dict)
    
    dict1 = select_and_remove_dict_key_suffix("_1", data_dict)
    dict2 = select_and_remove_dict_key_suffix("_2", data_dict)

    NormalShockProperties(
        FlowProperties(dict1),
        FlowProperties(dict2),
    )
end

function Base.getproperty(nsp::NormalShockProperties, s::Symbol)
    if s in fieldnames(NormalShockProperties)
        getfield(nsp, s)
    else
        if string(s)[(end-1):end] == "_1"
            getproperty(nsp.fp1, Symbol(string(s)[1:(end-2)]))
        elseif string(s)[(end-1):end] == "_2"
            getproperty(nsp.fp2, Symbol(string(s)[1:(end-2)]))
        else
            error("property $s cannot be found among fields or variables")
        end
    end
end

default_initial_guesses(::Type{NormalShockProperties}) = Dict(
    :gamma_2 => 1.4
)
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

units(::Type{Quasi1dimflowProperties}) = Dict(
    :mdot => u"kg/s",
    :A => u"m^2",
    :Astar => u"m^2",
    units(FlowProperties)...
)

function residues(qp::Quasi1dimflowProperties)
    [
        residues(qp.fp)
        qp.mdot - qp.A * qp.fp.v * qp.fp.cal_prop.mp.rho
        #p 682 anderson
        qp.A * qp.fp.M -
            qp.Astar * (
                2 / (qp.fp.cal_prop.gamma + 1) * (
                    1 + (qp.fp.cal_prop.gamma - 1) / 2 * qp.fp.M ^ 2
            )) ^ ((qp.fp.cal_prop.gamma + 1) / 2(qp.fp.cal_prop.gamma - 1))
    ]
end

##
#expand for N sections
#refactor
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
        nfp.sec1.fp.cal_prop.gamma - nfp.sec2.fp.cal_prop.gamma
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

function select_and_remove_dict_key_suffix(suff::String, dict)
    Dict(
        Symbol(string(pair.first)[1:(end-2)]) => pair.second
        for pair in filter(x -> endswith(string(x.first), suff), dict)
    )
end

function NozzleFlowProperties(data_dict::Dict)
    
    dict1 = select_and_remove_dict_key_suffix("_1", data_dict)
    dict2 = select_and_remove_dict_key_suffix("_2", data_dict)

    NozzleFlowProperties(
        Quasi1dimflowProperties(dict1),
        Quasi1dimflowProperties(dict2),
        data_dict[:F]
    )
end

default_initial_guesses(::Type{NozzleFlowProperties}) = Dict(
    :gamma_2 => 1.4
)

function add_units(nfp::NozzleFlowProperties, unit_dict)
    NozzleFlowProperties(
        add_units(nfp.sec1, select_and_remove_dict_key_suffix("_1", unit_dict)),
        add_units(nfp.sec2, select_and_remove_dict_key_suffix("_2", unit_dict)),
        nfp.F * unit_dict[:F]
    )
end
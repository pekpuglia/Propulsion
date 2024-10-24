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

const DEF_PRESSURE_UNIT = u"atm"
const DEF_MOLARCONC_UNIT = u"mol/m^3"
const DEF_TEMPERATURE_UNIT = u"K"

units(::Type{ThermodynamicProperties}) = Dict(
    :P => DEF_PRESSURE_UNIT,
    :z => DEF_MOLARCONC_UNIT,
    :T => DEF_TEMPERATURE_UNIT
)

const Rmolar = 1u"R"
const DEF_MOLARCALORCAP_UNIT = DEF_PRESSURE_UNIT / (DEF_MOLARCONC_UNIT * DEF_TEMPERATURE_UNIT) 

r_molar(::Number) = Rmolar
r_molar(::Real) = ustrip(DEF_MOLARCALORCAP_UNIT, Rmolar)

function residues(tp::ThermodynamicProperties)
    [tp.z*tp.T*r_molar(tp.P) - tp.P]
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

const DEF_MOLARMASS_UNIT = u"g/mol"
const DEF_DENSITY_UNIT = DEF_MOLARMASS_UNIT * DEF_MOLARCONC_UNIT
const DEF_MASSCALORCAP_UNIT = DEF_MOLARCALORCAP_UNIT / DEF_MOLARMASS_UNIT

units(::Type{MassProperties}) = Dict(
    :MM => DEF_MOLARMASS_UNIT,
    :rho => DEF_DENSITY_UNIT,
    :R => DEF_MASSCALORCAP_UNIT,
    units(ThermodynamicProperties)...
)

function residues(mp::MassProperties)
    [
        residues(mp.tp)
        mp.R * mp.MM / r_molar(mp.MM) - 1
        mp.MM * mp.tp.z - mp.rho
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

const DEF_SPEED_UNIT = √(DEF_MASSCALORCAP_UNIT * DEF_TEMPERATURE_UNIT)

units(::Type{CalorificProperties}) = Dict(
    :cv => DEF_MASSCALORCAP_UNIT,
    :cp => DEF_MASSCALORCAP_UNIT,
    :gamma => NoUnits,
    :a => DEF_SPEED_UNIT,
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
    :v => DEF_SPEED_UNIT,
    :T0 => DEF_TEMPERATURE_UNIT,
    :rho0 => DEF_DENSITY_UNIT,
    :P0 => DEF_PRESSURE_UNIT,
    :a0 => DEF_SPEED_UNIT,
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
        nsp.fp1.cal_prop.gamma - nsp.fp2.cal_prop.gamma
        nsp.fp1.cal_prop.mp.R - nsp.fp2.cal_prop.mp.R
        (nsp.fp1.cal_prop.gamma * nsp.fp1.M^2 - (nsp.fp1.cal_prop.gamma - 1) / 2) * (nsp.fp2.cal_prop.gamma * nsp.fp2.M^2 - (nsp.fp2.cal_prop.gamma - 1) / 2) - ((nsp.fp1.cal_prop.gamma + 1)/2)^2
        # nsp.fp2.M - sqrt((1 + (nsp.fp1.gamma - 1)/2 * nsp.fp1.M^2) / (nsp.fp1.gamma * nsp.fp1.M^2 - (nsp.fp1.gamma - 1)/2))
        nsp.fp1.T0 - nsp.fp2.T0
        nsp.fp2.cal_prop.mp.tp.P - nsp.fp1.cal_prop.mp.tp.P * (1 + 2nsp.fp1.cal_prop.gamma*(nsp.fp1.M^2 - 1) / (nsp.fp1.cal_prop.gamma + 1))
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

# function Base.getproperty(nsp::NormalShockProperties, s::Symbol)
#     if s in fieldnames(NormalShockProperties)
#         getfield(nsp, s)
#     else
#         if string(s)[(end-1):end] == "_1"
#             getproperty(nsp.fp1, Symbol(string(s)[1:(end-2)]))
#         elseif string(s)[(end-1):end] == "_2"
#             getproperty(nsp.fp2, Symbol(string(s)[1:(end-2)]))
#         else
#             error("property $s cannot be found among fields or variables")
#         end
#     end
# end

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

const DEF_AREA_UNIT = u"m^2"
const DEF_MASS_FLOW_UNIT = DEF_AREA_UNIT * DEF_DENSITY_UNIT * DEF_SPEED_UNIT

units(::Type{Quasi1dimflowProperties}) = Dict(
    :mdot => DEF_MASS_FLOW_UNIT,
    :A => DEF_AREA_UNIT,
    :Astar => DEF_AREA_UNIT,
    units(FlowProperties)...
)

function residues(qp::Quasi1dimflowProperties)
    [
        residues(qp.fp)
        qp.A * qp.fp.v * qp.fp.cal_prop.mp.rho / qp.mdot - 1
        #p 682 anderson
        qp.A * qp.fp.M -
            qp.Astar * (
                2 / (qp.fp.cal_prop.gamma + 1) * (
                    1 + (qp.fp.cal_prop.gamma - 1) / 2 * qp.fp.M ^ 2
            )) ^ ((qp.fp.cal_prop.gamma + 1) / 2(qp.fp.cal_prop.gamma - 1))
    ]
end

##
using StaticArrays
#make wrapper for kwarg constructor that computes the number of sections
export NozzleFlowProperties
struct NozzleFlowProperties4{N} <: PhysicalProperties
    secs::SVector{N, Quasi1dimflowProperties}
    F
    function NozzleFlowProperties4(secs::Vector{Quasi1dimflowProperties}, F::Real)
        new{length(secs)}(secs, F)
    end
    function NozzleFlowProperties4(secs::Vector{Quasi1dimflowProperties}, F::Unitful.Force)
        new{length(secs)}(secs, F)
    end
end

NozzleFlowProperties = NozzleFlowProperties4

function variables(T::Type{NozzleFlowProperties{N}}) where N
    vcat(
        [variables(Quasi1dimflowProperties) .|> string .|> (str -> str*"_$(i)") .|> Symbol for i in 1:N]...,
        :F
    )
end

#precisa ser igual a DEF_PRESSURE_UNIT*DEF_AREA_UNIT!
const DEF_FORCE_UNIT = DEF_MASS_FLOW_UNIT * DEF_SPEED_UNIT

function units(T::Type{NozzleFlowProperties{N}}) where N
    qfp_unit_dict = units(Quasi1dimflowProperties)
    Dict(
        Dict(
            var => qfp_unit_dict[q1dfp_var] for (var, q1dfp_var) in zip(variables(T), repeat(variables(Quasi1dimflowProperties), N))
        )...,
        :F => DEF_FORCE_UNIT
    )
end

function residues(nfp::NozzleFlowProperties{N}) where N
    vcat(
        (residues(sec) for sec in nfp.secs)...,
        ([
            nfp.secs[1].fp.cal_prop.gamma - seci.fp.cal_prop.gamma
            nfp.secs[1].fp.cal_prop.mp.R - seci.fp.cal_prop.mp.R
            nfp.secs[1].Astar - seci.Astar
            nfp.secs[1].mdot - seci.mdot
            nfp.secs[1].fp.T0 - seci.fp.T0
        ] for seci in nfp.secs[2:end])...,
        nfp.F - (
            nfp.secs[end].mdot * nfp.secs[end].fp.v - nfp.secs[end].fp.cal_prop.mp.tp.P * nfp.secs[end].A
        ) + (
            nfp.secs[1].mdot * nfp.secs[1].fp.v - nfp.secs[1].fp.cal_prop.mp.tp.P * nfp.secs[1].A
        )
    )
end

function Base.getindex(nfp::NozzleFlowProperties{N}, i::Int) where N
    nfp.secs[i]
end

function select_and_remove_dict_key_suffix(suff::String, dict)
    Dict(
        Symbol(string(pair.first)[1:(end-2)]) => pair.second
        for pair in filter(x -> endswith(string(x.first), suff), dict)
    )
end

function NozzleFlowProperties4{N}(data_dict::Dict) where N
    
    # max_ind = maximum(parse(Int, last(split_res)) for split_res in split.(String.(keys(data_dict)), "_") if length(split_res) > 1)

    dicts = (select_and_remove_dict_key_suffix("_$i", data_dict) for i in 1:N)

    NozzleFlowProperties(
        Quasi1dimflowProperties.(dicts),
        data_dict[:F]
    )
end

default_initial_guesses(::Type{NozzleFlowProperties}) = Dict(
    :gamma_2 => 1.4
)

function add_units(nfp::NozzleFlowProperties{N}, unit_dict) where N
    NozzleFlowProperties(
        [add_units(seci, select_and_remove_dict_key_suffix("_$i", unit_dict)) for (i, seci) in enumerate(nfp.secs)],
        nfp.F * unit_dict[:F]
    )
end
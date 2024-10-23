using Unitful


abstract type PhysicalProperties end

function (T::Type{<:PhysicalProperties})(data_dict::AbstractDict)
    vars = fieldnames(T)
    types = fieldtypes(T)

    indexes_to_recurse = findall(types .<: PhysicalProperties)

    parameters = [
        (i in indexes_to_recurse) ? types[i](data_dict) : data_dict[var] 
        for (i, var) in enumerate(vars)
    ]

    T(parameters...)
end

function Base.Dict(pp::T) where T <: PhysicalProperties
    Dict(var => getindex(pp, var) for var in variables(T))
end

export dof
function dof(T::Type{<:PhysicalProperties})
    (T |> variables |> length)- (T |> sym_vars |> T |> residues |> length)
end

export variables
function variables(T::Type{<:PhysicalProperties})
    vars = fieldnames(T)
    types = fieldtypes(T)

    indexes_to_recurse = findall(types .<: PhysicalProperties)

    cat(((i in indexes_to_recurse) ? variables(types[i]) : var for (i, var) in enumerate(vars))..., dims=1)
end

function add_units(pp::T, unit_dict) where T <: PhysicalProperties
    # vars = fieldnames(T)
    # types = fieldtypes(T)
    # indexes_to_recurse = findall(types .<: PhysicalProperties)
    # parameters = cat(((i in indexes_to_recurse) ? add_units(getfield(pp, var), unit_dict) : getfield(pp, var) * unit_dict[var] for (i, var) in enumerate(vars))..., dims=1)

    parameters = Dict(
        var => getindex(pp, var) * unit_dict[var]
        for var in variables(T)
    )

    T(parameters)
end

function Base.getindex(pp::T, s::Symbol) where T <: PhysicalProperties
    fields = fieldnames(T)
    types = fieldtypes(T)

    subtype_indices = findall(types .<: PhysicalProperties)

    if s in fields
        return getfield(pp, s)
    else
        #single recursion case
        if length(subtype_indices) == 1
            return getindex(getfield(pp, fields[subtype_indices[1]]), s)
        #strip index and recurse
        else
            subsym, ind = split(String(s), "_")
            subsym = Symbol(subsym)
            ind = parse(Int, ind)
            subtype_index = subtype_indices[ind]
            return getindex(getfield(pp, fields[subtype_index]), subsym)
        end
    end
end

units(T::Type{<:PhysicalProperties}) = Dict(var => NoUnits for var in variables(T))

export residues
residues(::T) where T <: PhysicalProperties = error("PhysicalProperties types must implement residues")

default_initial_guesses(::Type{<:PhysicalProperties}) = Dict()

export SymbolicParticipationVariable

struct SymbolicParticipationVariable <: Real
    variables::Vector{Symbol}
end

function SymbolicParticipationVariable(x::Number)
    SymbolicParticipationVariable([])
end

for op in [(:+), (:-), (:*), (:/), (:^)]
    @eval import Base.$op
    @eval $op(spvs::SymbolicParticipationVariable...) = SymbolicParticipationVariable(vcat(getfield.(spvs, :variables)...))
end

import Base.promote_rule
promote_rule(::Type{SymbolicParticipationVariable}, ::Type{<:Number}) = SymbolicParticipationVariable

## symbolic analysis tools
function sym_vars(T::Type{<:PhysicalProperties})
    Dict(var => SymbolicParticipationVariable([var]) for var in variables(T))
end

function participation_vector(T::Type{<:PhysicalProperties})
    T |> sym_vars |> T |> residues .|> (x->x.variables)
end

#incluir unidade?
#default VD for specific variables
#constructors: [lb, ub], IG
struct VariableData
    initial_guess
    lower_bound
    upper_bound
end

export VariableData
function Unitful.uconvert(unit, vd::VariableData)
    VariableData(uconvert(unit, vd.initial_guess), uconvert(unit, vd.lower_bound), uconvert(unit, vd.upper_bound))
end

function Unitful.ustrip(vd::VariableData)
    VariableData(ustrip(vd.initial_guess), ustrip(vd.lower_bound), ustrip(vd.upper_bound))
end

function VariableData()
    VariableData(1, eps(), Inf)
end

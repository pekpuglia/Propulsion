using Unitful, NonlinearSolve, Symbolics


abstract type PhysicalProperties end

#make T -> var dict for convinience
#make into Dict and PhysicalProperties constructors
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

function Base.Dict(pp::PhysicalProperties)
    Dict(var => getproperty(pp, var) for var in variables(pp))
end

export dof
dof(::Type{<:PhysicalProperties}) = error("PhysicalProperties types must implement dof")

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
        var => getproperty(pp, var) * unit_dict[var]
        for var in variables(T)
    )

    T(parameters)
end

Base.propertynames(::T) where T <: PhysicalProperties = variables(T)

##ONLY WORKS FOR SINGLE RECURSION CASE
#does not allow for accessing inner PhysicalProperties!!!
function Base.getproperty(pp::T, s::Symbol) where T <: PhysicalProperties
    vars = fieldnames(T)
    types = fieldtypes(T)

    index_to_recurse = findall(types .<: PhysicalProperties)

    @assert length(index_to_recurse) <= 1 "Cannot delegate property access due to ambiguity in field name $s"

    (s in vars) ? getfield(pp, s) : getproperty(getfield(pp, vars[index_to_recurse[1]]), s)
end

#melhorar
units(T::Type{<:PhysicalProperties}) = Dict(var => NoUnits for var in variables(T))

export residues
residues(::T) where T <: PhysicalProperties = error("PhysicalProperties types must implement residues")

default_initial_guesses(::Type{<:PhysicalProperties}) = Dict()

## symbolic analysis tools
function generate_sym_var_dict(T::Type{<:PhysicalProperties})
    Dict(var => (@variables $var)[1] for var in variables(T))
end

function participation_vector(T::Type{<:PhysicalProperties})
    sym_var_dict = generate_sym_var_dict(T)
    part_vector_symbolic = Symbolics.get_variables.(residues(T(sym_var_dict)))
    
    rev_svd = Dict(values(sym_var_dict) .=> keys(sym_var_dict))

    map(v -> getindex.([rev_svd], v), part_vector_symbolic)
end


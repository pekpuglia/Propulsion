#debug help
export sym_substitution_dict
function sym_substitution_dict(T::Type{<:PhysicalProperties}, input_variables::Dict{Symbol, <:Real}, default_value=1.0)
    Dict(sym_vars(T)[var] => (var in keys(input_variables)) ? input_variables[var] : (
        (var in keys(default_initial_guesses(T))) ? default_initial_guesses(T)[var] : default_value
    ) for var in variables(T))
end
export adjacency_list
function adjacency_list(T::Type{<:PhysicalProperties})
    pv = participation_vector(T)
    vars = variables(T)
    Dict(
        v => unique(
            filter(
                !=(v),
                vcat(
                    pv[
                        findall(eq_members -> v ∈ eq_members, pv)
                    ]...
                )
            )
        )
        for (i, v) in enumerate(vars)
    )
end

function neighbors(vertex::Symbol, adj_list::Vector{Tuple{Symbol, Symbol}})
    adj_list |>
        filter(pair -> vertex ∈ pair) .|>
        filter(!=(vertex)) .|>
        first
end

export connected_subgraphs
function connected_subgraphs(T::Type{<:PhysicalProperties}, size::Int, remaining_variables_per_equation=nothing)
    #previous attempts
    # #[given_vars; filter(∉(given_vars), variables(T))]
    remaining_variables_per_equation = something(remaining_variables_per_equation, participation_vector(T))
    reordered_variables = vcat(sort(remaining_variables_per_equation, by = length)...) |> unique

    connected_subgraphs(adjacency_list(T), size, [], reordered_variables, true)
end

#connected_subgraphs(MassProperties, 3) |> collect should be 
# [:P, :T, :z]
# [:P, :z, :rho]
# [:P, :z, :MM]
# [:z, :rho, :MM]
# [:z, :MM, :R]
# T, z, MM
#adapted from https://stackoverflow.com/questions/75727217/fast-way-to-find-all-connected-subgraphs-of-given-size-in-python
function connected_subgraphs(
    adj_list::Dict{Symbol, Vector{Symbol}}, 
    size::Int,
    excluded, #should be unique
    possible, #should be unique
    first_call=false)
    #starting the iteration
    possible = filter(p -> p ∉ excluded, possible)
    
    if size == 0
        return [[]]
    end
    (
        [vertex, other_subgraphs...]
        for (i, vertex) in enumerate(possible)
            for other_subgraphs in connected_subgraphs(
                adj_list, size-1, 
                unique([excluded..., possible[1:i]...]), 
                (first_call) ? adj_list[vertex] : unique([adj_list[vertex]..., possible...])
                )
    )
end

#change to single clique repr
@enum CliqueDiagnostic TooManyEquations CliqueFound TooFewEquations NoCliqueFound
struct CliqueResult
    expected_clique_order::Int
    clique_equations::Vector{Int}
    clique_vars::Vector{Symbol}
    diagnostic::CliqueDiagnostic
    #assumes all elements have the same length
    function CliqueResult(expected_clique_order, clique_equations, clique_vars)
        diagnostic = if length(clique_vars) > length(clique_equations)
                TooFewEquations
            elseif length(clique_vars) == length(clique_equations) == expected_clique_order
                CliqueFound
            elseif length(clique_vars) < length(clique_equations)
                TooManyEquations
            else
                NoCliqueFound
            end
        new(
            expected_clique_order,
            clique_equations,
            clique_vars,
            diagnostic
        )
    end
end

#https://github.com/JuliaLang/julia/issues/43737
findeach(testf::Function, A) = (first(p) for p in enumerate(A) if testf(last(p)))
selectindices(v, indices) = Iterators.map(pair -> last(pair),
    Iterators.filter(
        i_vars -> first(i_vars) in indices, 
        enumerate(v))
)

#inverter participation vector?
#buscar subgrafos com equações como nós?

#return first found, reject wrong sizes
#N cliques de ordem M são um clique de ordem N*M
#nem todas as variáveis precisam aparecer em todas as equações
#a + b = 0
#b + c = 0
#a + c = 0
#é um clique de ordem 3
using Base.Iterators
using IterTools
export find_clique
function find_clique(
    T::Type{<:PhysicalProperties}, 
    given_vars::AbstractVector{Symbol}, #remover! - redundante com remaining_variables_per_equation
    clique_order::Int,
    remaining_variables_per_equation=nothing
) #list of sets of equations with the same remaining variables
    pv = participation_vector(T)
    allvars = variables(T)
    missingvars = setdiff(allvars, given_vars)

    if isnothing(remaining_variables_per_equation)
        remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
    end

    variable_participation_dict = Dict(
        var => findall(vars_in_eq -> var ∈ vars_in_eq, pv)
        for var in allvars
    )
    #remove input variables
    variable_subgraphs = connected_subgraphs(T, clique_order, remaining_variables_per_equation)

    ret = CliqueResult(clique_order, [], [])

    nsubgraphs = 0
    #var_subgraph guaranteed to have clique_order variables
    for var_subgraph in variable_subgraphs
        possible_vars = filter(var -> var ∈ missingvars, var_subgraph)
        engaged_equations = unique(
            #filter equations that only have vars ∈ possible_vars
            filter(
                eq -> all(var_eq ∈ var_subgraph for var_eq in remaining_variables_per_equation[eq]),
                vcat(
                    getindex.([variable_participation_dict], possible_vars)...
                )
            )
        )
        ret = CliqueResult(clique_order, engaged_equations, possible_vars)
        nsubgraphs += 1
        if ret.diagnostic == CliqueFound || ret.diagnostic == TooManyEquations
            break
        end
    end
    # display(nsubgraphs)
    return ret
end
#find_clique(MassProperties, [:P, :MM, :T, :rho], 1) - should return :z is overconstrained
#find_clique(FlowProperties, [:P, :MM, :rho, :M, :gamma, :R, :z, :P0, :rho0, :T, :a, :T0, :v, :a0], 2) - find cp, cv

struct ResidueFunctionParameters
    other_vars_dict::Dict
    vars_to_solve_for::Vector{Symbol}
    residue_indices::Vector{Int}
end

using Optimization, OptimizationOptimJL
function solve_step(T, cr::CliqueResult, known_data, missingvars, full__missing_variable_data::Dict{Symbol, VariableData})
    residue_index = cr.clique_equations

    #3 groups of variables: known, to_solve, other. rewrite residues so others is not needed?
    vars_to_solve_for = cr.clique_vars

    #guess for missing vars which won't be solved this iteration - needed?
    other_missing_vars = setdiff(missingvars, vars_to_solve_for, keys(known_data))
    other_variable_data = getindex.(Ref(full__missing_variable_data), other_missing_vars)
    other_missing_values = getfield.(other_variable_data, :initial_guess)
    other_vars_dict = Dict(known_data..., Dict(other_var => other_value 
        for (other_var, other_value) in zip(
                                    other_missing_vars, 
                                    other_missing_values
    ))...)


    residue_function = (values, p::ResidueFunctionParameters) -> residues(
        T(Dict(
            Dict(
                var => value for (var, value) in zip(p.vars_to_solve_for, values)
            )..., p.other_vars_dict...
        )))[p.residue_indices]
    
    opt_function = OptimizationFunction(
        (values, p) -> residue_function(values, p) .^ 2 |> sum,
        AutoForwardDiff()
    )
    
    variable_data_to_solve = getindex.(Ref(full__missing_variable_data), vars_to_solve_for)
    initial_guess = getfield.(variable_data_to_solve, :initial_guess) .|> Float64
    lower_bound = getfield.(variable_data_to_solve, :lower_bound) .|> Float64
    upper_bound = getfield.(variable_data_to_solve, :upper_bound) .|> Float64

    opt_problem = OptimizationProblem(
        opt_function, initial_guess,
        ResidueFunctionParameters(other_vars_dict, vars_to_solve_for, residue_index),
        lb = lower_bound,
        ub = upper_bound)
    sol = solve(opt_problem, Optim.IPNewton())
    @assert Bool(sol.retcode)

    #return new known_data and setdiff
    Dict(
        known_data...,
        [var => val for (var, val) in zip(vars_to_solve_for, sol.u)]...
    ), setdiff(missingvars, vars_to_solve_for)
end
#initial_guesses must be missingvars
#melhorar
function internal_solver(T::Type, input_variables::Union{Dict{Symbol, <:Real}, Vector{Symbol}}, variable_data::Dict=Dict())
    pv = participation_vector(T)
    allvars = variables(T)
    sym_residues = residues(T(sym_vars(T)))

    numerically_solve = (input_variables isa Dict)

    given_vars = (numerically_solve) ? keys(input_variables) |> collect : input_variables
    missingvars = setdiff(allvars, given_vars) |> collect
    
    remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)

    if any(isempty.(remaining_variables_per_equation))
        over_constrained_index = isempty.(remaining_variables_per_equation)
        over_constrained_equation_variables = participation_vector(T)[over_constrained_index][1]
        error("Over-constrained equation: $(sym_residues[over_constrained_index]). Must not specify $over_constrained_equation_variables all at once")
    end

    known_data = deepcopy(input_variables)

    if numerically_solve
        #pass var to default VariableData constructor
        full__missing_variable_data = Dict(
            var => (var in keys(variable_data)) ? variable_data[var] : VariableData() for var in missingvars
        )
    end


    max_clique_order = length(missingvars)

    clique_order = 1
    for attempts_at_finding_solvable_var in 1:max_clique_order
        #overconstraint_validation used given_vars for the 2nd arg - is this wrong?
        clique_result = find_clique(T, 
            (numerically_solve) ? keys(known_data) |> collect : known_data, 
            clique_order, remaining_variables_per_equation)
        
        if clique_result.diagnostic == TooManyEquations
            over_constrained_variables = clique_result.clique_vars
            over_constrained_equations = sym_residues[clique_result.clique_equations]

            error("variables $over_constrained_variables are overconstrained\nby the following equations: $over_constrained_equations")
        end
        # display(clique_result)
        if clique_result.diagnostic == CliqueFound
            #solve for var list
            if numerically_solve
                #union(keys(known_data), missing_vars) = variables(T)
                known_data, missingvars = solve_step(T, clique_result, known_data, missingvars, full__missing_variable_data)
            else
                push!(known_data, clique_result.clique_vars...)
                filter!(∈(clique_result.clique_vars), missingvars)
            end
            remaining_variables_per_equation = [
                filter(∉(clique_result.clique_vars), rem_vars)
                for rem_vars in remaining_variables_per_equation
            ]
        end
        
        #try higher order cliques
        if clique_order < max_clique_order && clique_result.diagnostic != CliqueFound
            clique_order += 1
        end
        #try again from the simple equations
        if clique_order > 1 && clique_result.diagnostic == CliqueFound
            clique_order = 1
        end
    end
    (numerically_solve) ? T(known_data) : all(remaining_variables_per_equation .|> isempty)
end

function internal_solver(T::Type, input_variables::Dict{Symbol, <:Number}, variable_data::Dict)
    internal_units = units(T)

    unitless_kwargs = Dict(key => ustrip(internal_units[key], val) for (key, val) in input_variables)

    unitless_guesses = Dict(key => ustrip(internal_units[key], val) for (key, val) in variable_data)

    unitless_solution = internal_solver(T, unitless_kwargs, unitless_guesses)
     
    add_units(unitless_solution, internal_units)
end
##
function (T::Type{<:PhysicalProperties})(; kwargs...)
    allvars = variables(T)
    kwargs = Dict(kwargs)

    data_keys = filter(s -> endswith(string(s), "_data"), keys(kwargs))
    # data_kwargs = filter(x -> x.first ∈ data_keys, kwargs)
    data_kwargs = Dict{Symbol, VariableData}(Symbol(string(key)[1:(end-length("_data"))]) => kwargs[key] for key in data_keys)


    input_kwargs = filter(p -> p.first ∉ data_keys, kwargs)
    if all(isa.(values(input_kwargs), Real))
        input_kwargs = Dict{Symbol, Real}(input_kwargs)
    else
        input_kwargs = Dict{Symbol, Number}(input_kwargs)
    end

    #correct parameters validation
    if any([!(k in allvars) for k in keys(input_kwargs)])
        error("expected keys from $allvars, got: $(keys(input_kwargs)) ")
    end

    #dof validation
    if length(input_kwargs) != dof(T)
        error("$(dof(T)) thermodynamic properties needed, $(length(input_kwargs)) given: $(keys(input_kwargs))")
    end
    
    # overconstraint_validation(T, keys(kwargs) |> collect)

    internal_solver(T, input_kwargs, data_kwargs)
end
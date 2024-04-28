#debug help
export sym_substitution_dict
function sym_substitution_dict(T::Type{<:PhysicalProperties}, input_data::Dict{Symbol, <:Real}, default_value=1.0)
    Dict(sym_vars(T)[var] => (var in keys(input_data)) ? input_data[var] : (
        (var in keys(default_initial_guesses(T))) ? default_initial_guesses(T)[var] : default_value
    ) for var in variables(T))
end
export adjacency_list
function adjacency_list(T::Type{<:PhysicalProperties})
    pv = participation_vector(T)
    vars = variables(T)
    [
        (v, other)
        for (i, v) in enumerate(vars[1:(end-1)]) for other in unique(
            filter(
                cand -> cand ∈ vars[(i+1):end],
                vcat(
                    pv[
                        findall(eq_members -> v ∈ eq_members, pv)
                    ]...
                )
            )
        )
    ]
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
            elseif length(clique_vars) == length(clique_equations) == 0
                NoCliqueFound
            elseif length(clique_vars) < length(clique_equations)
                TooManyEquations
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
    given_vars::AbstractVector{Symbol}, 
    clique_order::Int,
    remaining_variables_per_equation=nothing
) #list of sets of equations with the same remaining variables

    if isnothing(remaining_variables_per_equation)
        allvars = variables(T)
        missingvars = setdiff(allvars, given_vars)
        pv = participation_vector(T)
        remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
    end

    clique_candidate_subsets = subsets(
        findall(!isempty, remaining_variables_per_equation), clique_order)

    unique_vars = Iterators.map(eq_subset -> unique(
        cat(remaining_variables_per_equation[eq_subset]..., dims=1)), clique_candidate_subsets)
    
    clique_equations = []
    clique_vars = []
    
    nonempty_indices = findeach(x -> !isempty(x), unique_vars)
    unique_vars = selectindices(unique_vars, nonempty_indices)

    clique_candidate_subsets = selectindices(clique_candidate_subsets, nonempty_indices)

    for (equation_subset, unique_var) in zip(clique_candidate_subsets, unique_vars)
        
        clique_equations = equation_subset
        clique_vars = unique_var

        #usar equation subset?
        if length(unique_var) == clique_order
            #check that there is no other subset with the same unique variables
            #in this case, include the other equations as well
            subsets_with_only_these_variables = selectindices(clique_candidate_subsets, findeach(==(unique_var), unique_vars))
            clique_equations = cat(subsets_with_only_these_variables..., dims=1)
            break
        end
    end
    
    CliqueResult(clique_order, clique_equations, clique_vars)
end
#find_clique(MassProperties, [:P, :MM, :T, :rho], 1) - should return :z is overconstrained
#find_clique(FlowProperties, [:P, :MM, :rho, :M, :gamma, :R, :z, :P0, :rho0, :T, :a, :T0, :v, :a0], 2) - find cp, cv

using Optimization, OptimizationOptimJL
function solve_step(T, cr::CliqueResult, known_data, missingvars, initial_guesses_dict)
    residue_index = cr.clique_equations
    vars_to_solve_for = cr.clique_vars
    #i have known_data
    #guess for missing vars which won't be solved this iteration
    other_missing_vars = setdiff(missingvars, vars_to_solve_for, keys(known_data))
    other_missing_values = getindex.(Ref(initial_guesses_dict), other_missing_vars)
    
    residue_function = (values, p) -> residues(T(Dict(
        Dict(var => value for (var, value) in zip(vars_to_solve_for, values))..., 
        Dict(other_var => other_value 
            for (other_var, other_value) in zip(
                                        other_missing_vars, 
                                        other_missing_values
        ))...,
        known_data...
    )))[residue_index]
    
    opt_function = OptimizationFunction(
        (values, p) -> residue_function(values, p) .^ 2 |> sum,
        AutoForwardDiff()
    )
    #missing: initial guess for variables to solve now
    initial_guess = getindex.(Ref(initial_guesses_dict), vars_to_solve_for)
    opt_problem = OptimizationProblem(
        opt_function, initial_guess,
        lb = √eps() * ones(size(initial_guess)),
        ub = fill(Inf, size(initial_guess)))
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
function internal_solver(T::Type, input_data::Union{Dict{Symbol, <:Real}, Vector{Symbol}}, input_initial_guesses::Dict=Dict())
    pv = participation_vector(T)
    allvars = variables(T)
    sym_residues = residues(T(sym_vars(T)))

    numerically_solve = (input_data isa Dict)

    given_vars = (numerically_solve) ? keys(input_data) |> collect : input_data
    missingvars = setdiff(allvars, given_vars) |> collect
    
    remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)

    if any(isempty.(remaining_variables_per_equation))
        over_constrained_index = isempty.(remaining_variables_per_equation)
        over_constrained_equation_variables = participation_vector(T)[over_constrained_index][1]
        error("Over-constrained equation: $(sym_residues[over_constrained_index]). Must not specify $over_constrained_equation_variables all at once")
    end

    known_data = deepcopy(input_data)

    initial_guesses_dict = Dict(mv => 
        Float64((mv ∈ keys(input_initial_guesses)) ? input_initial_guesses[mv] : 1.0)
        for mv in missingvars
    )


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

        if clique_result.diagnostic == CliqueFound
            #solve for var list
            if numerically_solve
                known_data, missingvars = solve_step(T, clique_result, known_data, missingvars, initial_guesses_dict)
            end
            remaining_variables_per_equation = remaining_variables_per_equation = [
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

function internal_solver(T::Type, input_data::Dict{Symbol, <:Number}, input_initial_guesses::Dict)
    internal_units = units(T)

    unitless_kwargs = Dict(key => ustrip(internal_units[key], val) for (key, val) in input_data)

    unitless_guesses = Dict(key => ustrip(internal_units[key], val) for (key, val) in input_initial_guesses)

    unitless_solution = internal_solver(T, unitless_kwargs, unitless_guesses)
     
    add_units(unitless_solution, internal_units)
end
##
function (T::Type{<:PhysicalProperties})(; kwargs...)
    allvars = variables(T)
    kwargs = Dict(kwargs)

    initial_keys = filter(s -> startswith(string(s), "initial_"), keys(kwargs))

    data_kwargs = filter(p -> p.first ∉ initial_keys, kwargs)
    if all(isa.(values(data_kwargs), Real))
        data_kwargs = Dict{Symbol, Real}(data_kwargs)
    else
        data_kwargs = Dict{Symbol, Number}(data_kwargs)
    end
    initial_kwargs = Dict(
        default_initial_guesses(T)...,
        Dict(Symbol(string(init_key)[9:end]) => kwargs[init_key]
        for init_key in initial_keys)...
    )

    #correct parameters validation
    if any([!(k in allvars) for k in keys(data_kwargs)])
        error("expected keys from $allvars, got: $(keys(data_kwargs)) ")
    end

    #dof validation
    if length(data_kwargs) != dof(T)
        error("$(dof(T)) thermodynamic properties needed, $(length(data_kwargs)) given: $(keys(data_kwargs))")
    end
    
    # overconstraint_validation(T, keys(kwargs) |> collect)


    internal_solver(T, data_kwargs, initial_kwargs)
end
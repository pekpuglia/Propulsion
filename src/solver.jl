#debug help
export sym_substitution_dict
function sym_substitution_dict(T::Type{<:PhysicalProperties}, input_data::Dict{Symbol, <:Real}, default_value=1.0)
    Dict(sym_vars(T)[var] => (var in keys(input_data)) ? input_data[var] : (
        (var in keys(default_initial_guesses(T))) ? default_initial_guesses(T)[var] : default_value
    ) for var in variables(T))
end

export opt_func_residues
function opt_func_residues(T, missingvars, input_data)
    (values, p) -> residues(T(Dict(
            Dict(missingvar => value for (missingvar, value) in zip(missingvars, values))..., 
            input_data...
        )))
end

using Optimization, OptimizationOptimJL, NonlinearSolve

export DEFAULT_OPT_PROB_GENERATOR
DEFAULT_OPT_PROB_GENERATOR = (f_value_p, u0, lb, ub) -> begin
    opt_func = OptimizationFunction(
        (values, p) -> f_value_p(values, p) .^2 |> sum,
        AutoForwardDiff()
    )

    OptimizationProblem(
        opt_func, u0, 
        iterations = 10000,
        lb = lb, ub = ub)
end

#initial_guesses must be missingvars
"""
opt_prob_generator(residue_value_p, initial_guesses_vec) -> prob
"""
function internal_solver(T::Type, input_data::Dict{Symbol, <:Real}, input_initial_guesses::Dict, 
        opt_prob_generator::Union{Function, Nothing} = nothing, 
        solver = nothing, 
        return_sol=false
    )
    pv = participation_vector(T)
    allvars = variables(T)
    given_vars = keys(input_data) |> collect
    missingvars = setdiff(allvars, given_vars) |> collect
    remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)

    known_data = deepcopy(input_data)

    initial_guesses_dict = Dict(mv => 
        Float64((mv ∈ keys(input_initial_guesses)) ? input_initial_guesses[mv] : 1.0)
        for mv in missingvars
    )


    max_clique_order = length(missingvars)

    clique_order = 1
    for attempts_at_finding_solvable_var in 1:max_clique_order
        clique_result = find_clique(T, keys(known_data) |> collect, clique_order, remaining_variables_per_equation)
        
        if clique_result.diagnostic == CliqueFound
            #solve for var list
            residue_index = clique_result.clique_equations
            vars_to_solve_for = clique_result.clique_vars
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

            #updating loop variables
            known_data = Dict(
                known_data...,
                [var => val for (var, val) in zip(vars_to_solve_for, sol.u)]...
            )
            
            missingvars = setdiff(missingvars, vars_to_solve_for)

            remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
        end
        
        #try higher order cliques
        if clique_order < max_clique_order && clique_result.diagnostic == NoCliqueFound
            clique_order += 1
        end
        #try again from the simple equations
        if clique_order > 1 && clique_result.diagnostic == CliqueFound
            clique_order = 1
        end
    end
    T(known_data)
end

function internal_solver(T::Type, input_data::Dict{Symbol, <:Number}, input_initial_guesses::Dict,
        opt_prob_generator = nothing, 
        solver = nothing,
        return_sol=false
    )
    internal_units = units(T)

    unitless_kwargs = Dict(key => ustrip(internal_units[key], val) for (key, val) in input_data)

    unitless_guesses = Dict(key => ustrip(internal_units[key], val) for (key, val) in input_initial_guesses)

    unitless_solution = internal_solver(T, unitless_kwargs, unitless_guesses, opt_prob_generator, solver, return_sol)
    
    return_sol = something(return_sol, false)
    
    if return_sol
        add_units(unitless_solution[1], internal_units), unitless_solution[2]
    else
        add_units(unitless_solution, internal_units)
    end
end
##
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
#update logic for detecting overconstraint - think about this
#N cliques de ordem M são um clique de ordem N*M
#nem todas as variáveis precisam aparecer em todas as equações
#a + b = 0
#b + c = 0
#a + c = 0
#é um clique de ordem 3
using Base.Iterators
using IterTools
function find_clique(
    T::Type{<:PhysicalProperties}, 
    given_vars::AbstractVector{Symbol}, 
    clique_order::Int,
    remaining_variables_per_equation=nothing
) #list of sets of equations with the same remaining variables

    allvars = variables(T)
    missingvars = setdiff(allvars, given_vars)
    if isnothing(remaining_variables_per_equation)
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




#clean
#make iterator on cliques
#TooFewEquations = continuar
#TooManyEquations = parar
export overconstraint_validation
function overconstraint_validation(T::Type{<:PhysicalProperties}, given_vars::AbstractVector{Symbol})
    pv = participation_vector(T)
    allvars = variables(T)
    sym_residues = residues(T(sym_vars(T)))
    missingvars = setdiff(allvars, given_vars)

    remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
    #try and remove this?
    if any(isempty.(remaining_variables_per_equation))
        over_constrained_index = isempty.(remaining_variables_per_equation)
        over_constrained_equation_variables = participation_vector(T)[over_constrained_index][1]
        error("Over-constrained equation: $(sym_residues[over_constrained_index]). Must not specify $over_constrained_equation_variables all at once")
    end


    max_clique_order = length(missingvars)

    #signal over/under-constrain!
    clique_order = 1
    for attempts_at_finding_solvable_var in 1:max_clique_order
        clique_result = find_clique(T, given_vars, clique_order, remaining_variables_per_equation)
        if clique_result.diagnostic == TooManyEquations
            over_constrained_variables = clique_result.clique_vars
            over_constrained_equations = sym_residues[clique_result.clique_equations]

            error("variables $over_constrained_variables are overconstrained\nby the following equations: $over_constrained_equations")
        end
        if clique_result.diagnostic == CliqueFound
            remaining_variables_per_equation = [
                filter(∉(clique_result.clique_vars), rem_vars)
                for rem_vars in remaining_variables_per_equation
            ]
        end
        #try higher order cliques
        if clique_order < max_clique_order && clique_result.diagnostic == NoCliqueFound
            clique_order += 1
        end
        #try again from the simple equations
        if clique_order > 1 && clique_result.diagnostic == CliqueFound
            clique_order = 1
        end
    end
    #should always be true, if it returns false 
    #and no error pops up there's a bug 🤷
    all(remaining_variables_per_equation .|> isempty)
end

##
function (T::Type{<:PhysicalProperties})(; kwargs...)
    allvars = variables(T)
    kwargs = Dict(kwargs)

    opt_prob_gen = pop!(kwargs, :opt_prob_gen, nothing)
    solver = pop!(kwargs, :solver, nothing)
    return_sol = pop!(kwargs, :return_sol, nothing)

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
    
    overconstraint_validation(T, keys(kwargs) |> collect)


    internal_solver(T, data_kwargs, initial_kwargs, opt_prob_gen, solver, return_sol)
end
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
        
        found_clique_indices = findall(==(CliqueFound), clique_result.diagnostic)
        found_clique_equations = clique_result.clique_equations[found_clique_indices]
        found_clique_variables = clique_result.clique_vars[found_clique_indices]
        
        for (var_list, eq_list) in zip(found_clique_variables, found_clique_equations)
            #solve for var list
            residue_index = eq_list
            vars_to_solve_for = var_list
            #i have known_data
            #guess for missing vars which won't be solved this iteration
            other_missing_vars = setdiff(missingvars, vars_to_solve_for, keys(known_data))
            other_missing_values = getindex.(Ref(initial_guesses_dict), other_missing_vars)
            
            residue_function = (values, p) -> residues(T(Dict(
                Dict(var => value for (var, value) in zip(var_list, values))..., 
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
            initial_guess = getindex.(Ref(initial_guesses_dict), var_list)
            opt_problem = OptimizationProblem(
                opt_function, initial_guess,
                lb = √eps() * ones(size(initial_guess)),
                ub = fill(Inf, size(initial_guess)))
            sol = solve(opt_problem, Optim.IPNewton())
            @assert Bool(sol.retcode)

            #updating loop variables
            known_data = Dict(
                known_data...,
                [var => val for (var, val) in zip(var_list, sol.u)]...
            )
            
            missingvars = setdiff(missingvars, var_list)

            remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
        end
        
        #try higher order cliques
        if clique_order < max_clique_order && isempty(found_clique_equations)
            clique_order += 1
        end
        #try again from the simple equations
        if clique_order > 1 && !isempty(found_clique_equations)
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
@enum CliqueDiagnostic TooManyEquations CliqueFound TooFewEquations
struct CliqueResult
    expected_clique_order::Int
    clique_equations::Vector{Vector{Int}}
    clique_vars::Vector{Vector{Symbol}}
    diagnostic::Vector{CliqueDiagnostic}
    #assumes all elements have the same length
    function CliqueResult(expected_clique_order, clique_equations, clique_vars)
        diagnostic = [
            if length(vars) > length(eqs)
                TooFewEquations
            elseif length(vars) == length(eqs) == expected_clique_order
                CliqueFound
            elseif length(vars) < length(eqs)
                TooManyEquations
            end
            for (vars, eqs) in zip(clique_vars, clique_equations)
        ]

        #if found cliques, filter for them
        #otherwise return the full list of not found cliques
        if any(diagnostic .== CliqueFound)
            indices = findall(diagnostic .== CliqueFound)
            clique_equations = clique_equations[indices]
            clique_vars = clique_vars[indices]
            diagnostic = diagnostic[indices]
        end

        new(
            expected_clique_order,
            clique_equations,
            clique_vars,
            diagnostic
        )
    end
end

#N cliques de ordem M são um clique de ordem N*M

using IterTools
#nem todas as variáveis precisam aparecer em todas as equações
#a + b = 0
#b + c = 0
#a + c = 0
#é um clique de ordem 3
export find_clique
function find_clique(
        T::Type{<:PhysicalProperties}, 
        given_vars::AbstractVector{Symbol}, 
        clique_order::Int,
        remaining_variables_per_equation=nothing
    ) :: CliqueResult #list of sets of equations with the same remaining variables
    
    allvars = variables(T)
    missingvars = setdiff(allvars, given_vars)
    if isnothing(remaining_variables_per_equation)
        pv = participation_vector(T)
        remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
    end

    #wrong
    indices_equations_clique_order_remaining_variables = findall(rem_vars -> length(rem_vars) == clique_order, remaining_variables_per_equation)

    clique_equations = Vector{Vector{Int}}()
    clique_vars = Vector{Vector{Symbol}}()

    for ind in indices_equations_clique_order_remaining_variables
        rem_vars_i = remaining_variables_per_equation[ind]
        equations_with_those_variables = findall(rv -> Set(rem_vars_i) == Set(rv), remaining_variables_per_equation)
        if rem_vars_i ∉ clique_vars
            push!(clique_equations, equations_with_those_variables) 
            push!(clique_vars, rem_vars_i)
        end
    end
    CliqueResult(clique_order, clique_equations, clique_vars)
end

#https://github.com/JuliaLang/julia/issues/43737
findeach(testf::Function, A) = (first(p) for p in enumerate(A) if testf(last(p)))
selectindices(v, indices) = Iterators.map(pair -> last(first(pair)),
    Iterators.filter(
        vars_ind_pair -> first(first(vars_ind_pair)) in last(vars_ind_pair), 
        zip(enumerate(v), indices))
)

#return first found, reject wrong sizes
#update logic for detecting overconstraint - think about this
using Base.Iterators
function find_clique2(
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

    clique_candidate_subsets = subsets(1:lastindex(remaining_variables_per_equation), clique_order)

    unique_vars = Iterators.map(eq_subset -> unique(
        cat(remaining_variables_per_equation[eq_subset]..., dims=1)), clique_candidate_subsets)
    
    display("unique vars found")
    clique_equations = Vector{Vector{Int}}()
    clique_vars = Vector{Vector{Symbol}}()
    
    nonempty_indices = findeach(x -> !isempty(x), unique_vars)
    unique_vars = selectindices(unique_vars, nonempty_indices)

    clique_candidate_subsets = selectindices(clique_candidate_subsets, nonempty_indices)

    unique_vars = collect(unique_vars)
    clique_candidate_subsets = collect(clique_candidate_subsets)
    #check that no other subset has the same vars - remove this
    i = 0
    for (equation_subset, unique_var) in zip(clique_candidate_subsets, unique_vars)
        if i % 50 == 0
            display(i)
        end
        i += 1
        #usar equation subset?
        if unique_var ∉ clique_vars
            subset_indices_with_these_variables = findall(other_unique_variable_list -> 
                Set(other_unique_variable_list) == Set(unique_var), 
                unique_vars
            )

            new_found_clique_equations = cat(
                (clique_candidate_subsets)[subset_indices_with_these_variables]...,
                dims=1) |> unique

            filter!(eq -> !isempty(remaining_variables_per_equation[eq]), new_found_clique_equations)

            if length(new_found_clique_equations) == clique_order || length(unique_var) == clique_order
                display("found clique: $unique_var")
                push!(clique_equations, new_found_clique_equations)
                push!(clique_vars, unique_var)
            end
        end
    end
    
    CliqueResult(clique_order, clique_equations, clique_vars)
end
#find_clique(MassProperties, [:P, :MM, :T, :rho, :R], 1) - should return :z is overconstrained
#find_clique(FlowProperties, [:P, :MM, :rho, :M, :gamma, :R, :z, :P0, :rho0, :T, :a, :T0, :v, :a0], 2)
function test_find_clique_1_var(find_clique_function)
    clique_res = find_clique_function(MassProperties, [:P, :z, :MM], 1)
    clique_res.clique_equations == [[1], [2], [3]] &&
    clique_res.clique_vars == [[:T], [:R], [:rho]] &&
    all(clique_res.diagnostic .== CliqueFound)
end

function test_find_clique_2_var(find_clique_function)
    clique_res = find_clique_function(CalorificProperties, [:P, :R, :gamma], 2)
    i = findfirst(==(CliqueFound), clique_res.diagnostic)
    clique_res.clique_equations[i] == [4, 5] &&
    Set(clique_res.clique_vars[i]) == Set([:cp, :cv])
end

#clean
export overconstraint_validation
function overconstraint_validation(T::Type{<:PhysicalProperties}, given_vars::AbstractVector{Symbol})
    pv = participation_vector(T)
    allvars = variables(T)
    sym_residues = residues(T(sym_vars(T)))
    missingvars = setdiff(allvars, given_vars)

    remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
    if any(isempty.(remaining_variables_per_equation))
        over_constrained_index = isempty.(remaining_variables_per_equation)
        over_constrained_equation_variables = participation_vector(T)[over_constrained_index][1]
        error("Over-constrained equation: $(sym_residues[over_constrained_index]). Must not specify $over_constrained_equation_variables all at once")
    end


    max_clique_order = length(missingvars)

    #signal over/under-constrain!
    clique_order = 1
    for attempts_at_finding_solvable_var in 1:max_clique_order
        display(attempts_at_finding_solvable_var)
        clique_result = find_clique(T, given_vars, clique_order, remaining_variables_per_equation)
        if !isempty(clique_result.diagnostic) && any(clique_result.diagnostic .== TooManyEquations)
            over_constrained_indices = findall(==(TooManyEquations), clique_result.diagnostic)
            
            over_constrained_equation_indices = cat(clique_result.clique_equations[over_constrained_indices]..., dims=1)
            over_constrained_variables = cat(clique_result.clique_vars[over_constrained_indices]..., dims=1)
            
            over_constrained_equations = sym_residues[over_constrained_equation_indices]

            error("variables $over_constrained_variables are overconstrained\nby the following equations: $over_constrained_equations")
        end
        found_clique_indices = findall(==(CliqueFound), clique_result.diagnostic)
        found_clique_equations = clique_result.clique_equations[found_clique_indices]
        found_clique_variables = clique_result.clique_vars[found_clique_indices]
        for var_list in found_clique_variables
            remaining_variables_per_equation = [
                filter(∉(var_list), rem_vars)
                for rem_vars in remaining_variables_per_equation
            ]
        end
        #try higher order cliques
        if clique_order < max_clique_order && isempty(found_clique_equations)
            clique_order += 1
        end
        #try again from the simple equations
        if clique_order > 1 && !isempty(found_clique_equations)
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
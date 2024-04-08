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
    allvars = variables(T)

    missingvars = (setdiff(Set(allvars), Set(keys(input_data))) |> collect)
    
    initial_guesses_vec = [
        (mv ∈ keys(input_initial_guesses)) ? input_initial_guesses[mv] : 1.0
        for mv in missingvars
    ] .|> Float64

    #sol.original.minimum
    opt_prob_generator = something(opt_prob_generator, DEFAULT_OPT_PROB_GENERATOR)
    solver = something(solver, Optim.IPNewton())
    sol = solve(
        opt_prob_generator(
            opt_func_residues(T, missingvars, input_data), 
            initial_guesses_vec,
            √eps() * ones(size(initial_guesses_vec)),
            fill(Inf, size(initial_guesses_vec))
            ), solver)

    ret = T(Dict(
        Dict(missingvar => u for (missingvar, u) in zip(missingvars, sol.u))...,
        Dict(k => Float64(v) for (k, v) in input_data)...
    ))
    return_sol = something(return_sol, false)
    (return_sol) ? (ret, sol) : ret
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
#find all sets of 1 equation with 1 remaining variable
#later: find all sets of 2 equations with the same 2 remaining variables
#later: find all sets of N equations with the same N remaining variables
export find_clique
function find_clique(
        T::Type{<:PhysicalProperties}, 
        given_vars::AbstractVector{Symbol}, 
        clique_order::Int,
        remaining_variables_per_equation=nothing
    )
    if isnothing(remaining_variables_per_equation)
        pv = participation_vector(T)
        allvars = variables(T)
        missingvars = setdiff(allvars, given_vars)
        remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
    end
    
    indices_equations_clique_order_remaining_variables = findall(rem_vars -> length(rem_vars) == clique_order, remaining_variables_per_equation)

    #check if count is equal to clique order
    filter(
        i -> count(
            rem_vars -> Set(rem_vars) == Set(remaining_variables_per_equation[i]), 
            remaining_variables_per_equation) == clique_order,
        indices_equations_clique_order_remaining_variables
    )
end

function test_find_clique_0_var()
    find_clique(ThermodynamicProperties, [:P, :z, :T], 0) == [1]
end

function test_find_clique_1_var()
    find_clique(MassProperties, [:P, :z, :MM], 1) == [1, 2, 3]
end

function test_find_clique_2_var()
    find_clique(CalorificProperties, [:P, :R, :gamma], 2) == [4, 5]
end

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
    for clique_order = 1:max_clique_order
        for nth_order_clique_attempts = 1:length(missingvars)
            #index of equations with clique_order remaining vars
            n_remaining_var_eq_index = findfirst(==(clique_order), length.(remaining_variables_per_equation))
            if n_remaining_var_eq_index |> isnothing
                continue
            end
    
            newly_found_vars_vec = remaining_variables_per_equation[n_remaining_var_eq_index]
    
            #check if variable is overconstrained
            var_particip_indices = findall(==(newly_found_vars_vec), remaining_variables_per_equation)
            if length(var_particip_indices) > clique_order
                overconstrained_equations = sym_residues[var_particip_indices]
                overconstraint_order = length(var_particip_indices) - clique_order
                error("variables $newly_found_vars_vec are overconstrained, \nwith $overconstraint_order extra variable supplied. Equations: $overconstrained_equations")
            end
    
            filter!(∉(newly_found_vars_vec), missingvars)

            remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
        end
    end

    true
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
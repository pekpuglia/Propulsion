
#initial_guesses must be missingvars
function internal_solver(T::Type, input_data::Dict{Symbol, <:Real}, initial_guesses::Dict)
    allvars = variables(T)

    missingvars = (setdiff(Set(allvars), Set(keys(input_data))) |> collect)
    
    initial_guesses_vec = [
        (mv ∈ keys(initial_guesses)) ? initial_guesses[mv] : 1.0
        for mv in missingvars
    ]

    prob = NonlinearProblem(
        (values, p) -> residues(T(Dict(
                Dict(missingvar => value for (missingvar, value) in zip(missingvars, values))..., 
                input_data...
            ))), 
        initial_guesses_vec, p=()
    )

    sol = solve(prob, NewtonRaphson())

    T(Dict(
            Dict(missingvar => u for (missingvar, u) in zip(missingvars, sol.u))...,
            Dict(k => Float64(v) for (k, v) in input_data)...
    ))
end

function internal_solver(T::Type, input_data::Dict{Symbol, <:Number}, initial_guesses::Dict)
    internal_units = units(T)

    unitless_kwargs = Dict(key => ustrip(internal_units[key], val) for (key, val) in input_data)

    unitless_guesses = Dict(key => ustrip(internal_units[key], val) for (key, val) in initial_guesses)

    unitless_solution = internal_solver(T, unitless_kwargs, unitless_guesses)

    add_units(unitless_solution, internal_units)
end
##
export overconstraint_validation
function overconstraint_validation(T::Type{<:PhysicalProperties}, given_vars::AbstractVector{Symbol})
    pv = participation_vector(T)
    allvars = variables(T)
    missingvars = setdiff(allvars, given_vars)

    remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)

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
                overconstrained_equations = residues(T(generate_sym_var_dict(T)))[var_particip_indices]
                overconstraint_order = length(var_particip_indices) - clique_order
                return "variables $newly_found_vars_vec are overconstrained, with $overconstraint_order extra variable supplied. Equations: $overconstrained_equations"
            end
    
            filter!(∉(newly_found_vars_vec), missingvars)

            remaining_variables_per_equation = map(v -> v[v .∈ [missingvars]], pv)
        end
    end

    remaining_variables_per_equation
end

##
function (T::Type{<:PhysicalProperties})(; kwargs...)
    allvars = variables(T)

    initial_keys = filter(s -> startswith(string(s), "initial_"), keys(kwargs))

    data_kwargs = filter(p -> p.first ∉ initial_keys, kwargs)
    initial_kwargs = Dict(
        default_initial_guesses(T)...,
        Dict(Symbol(string(init_key)[9:end]) => kwargs[init_key]
        for init_key in initial_keys)...
    )

    #correct parameters validation
    if any([!(k in allvars) for k in keys(data_kwargs)])
        error("expected keys from $allvars, got: $(keys(data_kwargs)) ")
    end

    #dof validation - review dof calculation?
    if length(data_kwargs) != dof(T)
        error("$(dof(T)) thermodynamic properties needed, $(length(data_kwargs)) given: $(keys(data_kwargs))")
    end
    
    remaining_vars = map(v -> v[v .∉ [keys(data_kwargs)]], participation_vector(T))
    if any(isempty.(remaining_vars))
        over_constrained_equation_variables = participation_vector(T)[isempty.(remaining_vars)][1]
        error("Over-constrained equation. Must not specify $over_constrained_equation_variables all at once")
    end


    internal_solver(T, data_kwargs, initial_kwargs)
end
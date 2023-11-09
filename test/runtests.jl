include("../initial.jl")

using Test

##############################################################################
#usage tests
function test_access_property()
    mp = MassProperties(P = 1, T = 2, MM = 20)
    mp.P == 1
end

@test test_access_property()
############################################################################
#internal coherence tests/unit tests
function test_dof_variable_count(type)
    vars = variables(type)
    dummy_values = ones(size(vars))
    dummy_instance = phys_prop_from_kwargs(type; Dict(var => dummy_val for (var, dummy_val) in zip(vars, dummy_values))...)
    res = residues(dummy_instance)
    length(res) + dof(type) == length(vars)
end

for t in [
    ThermodynamicProperties,
    MassProperties,
    CalorificProperties,
    FlowProperties,
    Quasi1dimflowProperties]

    @test test_dof_variable_count(t)
end

function test_internal_solver()
    sol = internal_solver(ThermodynamicProperties, Dict(:P => 1.0, :T => 10.0))
    true
end

@test test_internal_solver()

function test_residue_unit_coherence()
    residues(phys_prop_from_kwargs(Quasi1dimflowProperties,
        P = 1u"Pa",
        T = 1u"K",
        z = 1u"mol/m^3",
        rho = 1u"kg/m^3",
        MM = 1u"kg/mol",
        R = 1u"J/kg/K",
        a = 1u"m/s",
        cp = 1u"J/kg/K",
        cv = 1u"J/kg/K",
        gamma = 1.4,
        M = 1.1,
        a0 = 1u"m/s",
        P0 = 1u"Pa",
        rho0 = 1u"kg/m^3",
        T0 = 1u"K",
        v = 1u"m/s",
        A = 1u"m^2",
        Astar = 1u"m^2",
        mdot = 1u"kg/s"
    ))
    true
end

@test test_residue_unit_coherence()

function test_reject_invalid_units()
    ThermodynamicProperties(P = 1u"Pa", T = 1.0)
end

@test_throws Unitful.DimensionError test_reject_invalid_units()
#############################################################################
# correctness

#anderson
function test_flow_properties()
    solution = FlowProperties(P = 1u"atm", T = 320u"K", v = 1000u"m/s", gamma = 1.4, R = 287u"J/kg/K")

    solution2 = FlowProperties(P0 = 2220u"lbf/ft^2", P=1455.6u"lbf/ft^2", T = 483.04u"Ra", gamma = 1.4, R = 287u"J/kg/K")

    isapprox(solution.T0, 817.8u"K", atol=1e-1u"K") &&
    isapprox(solution.P0, 26.7u"atm", atol=1e-1u"atm") &&
    isapprox(solution2.T0, 544.9u"Ra", atol=0.1u"Ra") &&
    isapprox(solution2.v, 862u"ft/s", atol=1u"ft/s")
end

@test test_flow_properties()
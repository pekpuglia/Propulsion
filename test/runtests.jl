include("../initial.jl")

using Test

##############################################################################
#usage tests



############################################################################
#internal coherence tests
function test_dof_variable_count()
    
end

@test test_dof_variable_count()

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
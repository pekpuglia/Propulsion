include("../initial.jl")

using Test
#anderson
function test_example_7_6()
    solution = solve_params(FlowProperties, P = 1u"atm", T = 320u"K", v = 1000u"m/s", gamma = 1.4, R = 287u"J/kg/K")

    isapprox(solution.T0, 817.8u"K", atol=1e-1u"K") &&
    isapprox(solution.P0, 26.7u"atm", atol=1e-1u"atm")
end

@test test_example_7_6()

function test_example_7_7()
    solution = solve_params(FlowProperties, P0 = 2220u"lb/ft^2", rho=1455.6u"lb/ft^2", T = 483.04u"Ra", gamma = 1.4)

    isapprox(solution.T0, 544.9u"Ra", atol=0.1u"Ra") &&
    isapprox(solution.v, 862u"ft/s", atol=1u"ft/s")
end

@test test_example_7_7()
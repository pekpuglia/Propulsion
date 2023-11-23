include("../src/Propulsion.jl")

using .Propulsion

using Test, Unitful

##############################################################################
#usage tests
function test_access_property()
    MassProperties(P = 1, T = 2, MM = 20)
end

@test test_access_property().P == 1

function test_reject_invalid_units()
    ThermodynamicProperties(P = 1u"Pa", T = 1.0)
end

@test_throws Unitful.DimensionError test_reject_invalid_units()

function test_reject_overconstrained_system()
    MassProperties(P = 1, z = 1, T = 1)
end

@test_throws ErrorException test_reject_overconstrained_system()

function test_use_initial_value()
    Quasi1dimflowProperties(
        P=1e5, 
        T=10.0, 
        rho = 2.0, 
        gamma = 1.4, 
        Astar = 0.85, 
        A = 1.0, 
        initial_M=5
    )
end

@test test_use_initial_value().M > 1
############################################################################
#internal coherence tests/unit tests
function test_dof_variable_count(type)
    vars = variables(type)
    dummy_values = ones(size(vars))
    dummy_instance = Propulsion.phys_prop_from_kwargs(type; Dict(var => dummy_val for (var, dummy_val) in zip(vars, dummy_values))...)
    res = residues(dummy_instance)
    length(res) + dof(type) == length(vars)
end

@test test_dof_variable_count(ThermodynamicProperties)
@test test_dof_variable_count(MassProperties)
@test test_dof_variable_count(CalorificProperties)
@test test_dof_variable_count(FlowProperties)
@test test_dof_variable_count(Quasi1dimflowProperties)
@test test_dof_variable_count(NozzleFlowProperties)


function test_internal_solver()
    Propulsion.internal_solver(ThermodynamicProperties, Dict(:P => 1.0, :T => 10.0), Dict(:z => 2.0))
end

#just test it built something
@test test_internal_solver() isa ThermodynamicProperties

function test_nozzle_flow_from_kwargs()
    Propulsion.phys_prop_from_kwargs(NozzleFlowProperties,
    F = 2u"N",
    P_1 = 1u"Pa",
    T_1 = 1u"K",
    z_1 = 1u"mol/m^3",
    rho_1 = 1u"kg/m^3",
    MM_1 = 1u"kg/mol",
    R_1 = 1u"J/kg/K",
    a_1 = 1u"m/s",
    cp_1 = 1u"J/kg/K",
    cv_1 = 1u"J/kg/K",
    gamma_1 = 1.4,
    M_1 = 1.1,
    a0_1 = 1u"m/s",
    P0_1 = 1u"Pa",
    rho0_1 = 1u"kg/m^3",
    T0_1 = 1u"K",
    v_1 = 1u"m/s",
    A_1 = 1u"m^2",
    Astar_1 = 1u"m^2",
    mdot_1 = 1u"kg/s",
    P_2 = 1u"Pa",
    T_2 = 1u"K",
    z_2 = 1u"mol/m^3",
    rho_2 = 1u"kg/m^3",
    MM_2 = 1u"kg/mol",
    R_2 = 1u"J/kg/K",
    a_2 = 1u"m/s",
    cp_2 = 1u"J/kg/K",
    cv_2 = 1u"J/kg/K",
    gamma_2 = 1.4,
    M_2 = 1.1,
    a0_2 = 1u"m/s",
    P0_2 = 1u"Pa",
    rho0_2 = 1u"kg/m^3",
    T0_2 = 1u"K",
    v_2 = 1u"m/s",
    A_2 = 1u"m^2",
    Astar_2 = 1u"m^2",
    mdot_2 = 1u"kg/s"
    )
end

#testing it didn't error
@test test_nozzle_flow_from_kwargs() isa NozzleFlowProperties

function test_residue_unit_coherence()
    residues(Propulsion.phys_prop_from_kwargs(NozzleFlowProperties,
        F = 2u"N",
        P_1 = 1u"Pa",
        T_1 = 1u"K",
        z_1 = 1u"mol/m^3",
        rho_1 = 1u"kg/m^3",
        MM_1 = 1u"kg/mol",
        R_1 = 1u"J/kg/K",
        a_1 = 1u"m/s",
        cp_1 = 1u"J/kg/K",
        cv_1 = 1u"J/kg/K",
        gamma_1 = 1.4,
        M_1 = 1.1,
        a0_1 = 1u"m/s",
        P0_1 = 1u"Pa",
        rho0_1 = 1u"kg/m^3",
        T0_1 = 1u"K",
        v_1 = 1u"m/s",
        A_1 = 1u"m^2",
        Astar_1 = 1u"m^2",
        mdot_1 = 1u"kg/s",
        P_2 = 1u"Pa",
        T_2 = 1u"K",
        z_2 = 1u"mol/m^3",
        rho_2 = 1u"kg/m^3",
        MM_2 = 1u"kg/mol",
        R_2 = 1u"J/kg/K",
        a_2 = 1u"m/s",
        cp_2 = 1u"J/kg/K",
        cv_2 = 1u"J/kg/K",
        gamma_2 = 1.4,
        M_2 = 1.1,
        a0_2 = 1u"m/s",
        P0_2 = 1u"Pa",
        rho0_2 = 1u"kg/m^3",
        T0_2 = 1u"K",
        v_2 = 1u"m/s",
        A_2 = 1u"m^2",
        Astar_2 = 1u"m^2",
        mdot_2 = 1u"kg/s"
    ))
end

#testing it didn't error
@test test_residue_unit_coherence() isa Vector

#############################################################################
# correctness

#anderson
function test_flow_properties()
    #example 7.6
    solution = FlowProperties(P = 1u"atm", T = 320u"K", v = 1000u"m/s", gamma = 1.4, R = 287u"J/kg/K")
    #example 7.7
    solution2 = FlowProperties(P0 = 2220u"lbf/ft^2", P=1455.6u"lbf/ft^2", T = 483.04u"Ra", gamma = 1.4, R = 287u"J/kg/K")
    solution, solution2
end

let
    solution, solution2 = test_flow_properties()
    @test isapprox(solution.T0, 817.8u"K", atol=1e-1u"K")
    @test isapprox(solution.P0, 26.7u"atm", atol=1e-1u"atm")
    @test isapprox(solution2.T0, 544.9u"Ra", atol=0.1u"Ra")
    @test isapprox(solution2.v, 862u"ft/s", atol=1u"ft/s")
end

#Quasi1dimflowProperties tests
function test_example_10_1()
    Quasi1dimflowProperties(
        A = 10.25u"m^2", 
        Astar=1u"m^2",
        P0 = 5u"atm",
        T0 = 600u"Ra",
        gamma = 1.4,
        R = 287u"J/kg/K",
        initial_M = 5
    )
end

let sol = test_example_10_1()
    @test isapprox(sol.M, 3.95, atol=0.01)
    @test isapprox(sol.P, 0.035u"atm", atol=1e-3u"atm")
    @test isapprox(sol.T, 145.6u"Ra", atol=0.1u"Ra")
end

#NozzleFlowProperties test
function test_example_10_2()
    nfp_supersonic = NozzleFlowProperties(
        P0_1 = 1u"atm",
        T0_1 = 288u"K",
        M_1 = 1,
        gamma_1 = 1.4,
        R_1 = 287u"J/kg/K",
        Astar_1 = 1u"m^2",
        A_2 = 2u"m^2", 
        initial_M_2 = 5)

    nfp_subsonic = NozzleFlowProperties(
        P0_1 = 1u"atm",
        T0_1 = 288u"K",
        M_1 = 1,
        gamma_1 = 1.4,
        R_1 = 287u"J/kg/K",
        Astar_1 = 1u"m^2",
        A_2 = 2u"m^2", 
        initial_M_2 = 0.1)
        
    nfp_supersonic, nfp_subsonic
end

let
    nfp_supersonic, nfp_subsonic = test_example_10_2()
    
    @test isapprox(nfp_supersonic[1].P, 0.528u"atm", atol=1e-3u"atm")
    @test isapprox(nfp_supersonic[1].T, 240u"K", atol=1u"K")
    @test isapprox(nfp_supersonic[2].M, 2.2, atol=0.1)
    #small numerical difference
    @test isapprox(nfp_supersonic[2].P, 0.0935u"atm", atol=2e-3u"atm")
    @test isapprox(nfp_supersonic[2].T, 146u"K", atol=1u"K")
    @test isapprox(nfp_subsonic[2].M, 0.3, atol=0.1)
    @test isapprox(nfp_subsonic[2].P, 0.94u"atm", atol=1e-2u"atm")
    @test isapprox(nfp_subsonic[2].T, 282.9u"K", atol=0.2u"K")
end

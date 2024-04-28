using Propulsion

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
        P=1, 
        T=10.0, 
        MM=29, 
        gamma = 1.4, 
        Astar = 0.85, 
        A = 1.0, 
        initial_M=5
    )
end

@test test_use_initial_value().M > 1

#constraint tests
function test_mp_correct_input()
    right_given_vars_mp = [:P, :MM, :T]
    Propulsion.internal_solver(MassProperties, right_given_vars_mp)
end


function test_mp_wrong_input()
    wrong_given_vars_mp = [:P, :MM, :T, :rho]
    Propulsion.internal_solver(MassProperties, wrong_given_vars_mp)
end


function test_fp_right_input()
    right_given_vars_fp = [:P, :MM, :rho, :M, :gamma]
    Propulsion.internal_solver(FlowProperties, right_given_vars_fp)
end


function test_fp_wrong_input()
    wrong_given_vars_fp = [:P, :MM, :rho, :M, :T]
    Propulsion.internal_solver(FlowProperties, wrong_given_vars_fp)
end

@testset "overconstraint validation" begin
    @test test_mp_correct_input()
    @test_throws ErrorException test_mp_wrong_input()
    @test test_fp_right_input()
    @test_throws ErrorException test_fp_wrong_input()
end

############################################################################
#numerical tests found through bugs
#overconstraint_validation takes forever & fails
function test_normal_shock_inf_or_nan()
    nsp = NormalShockProperties(
           v_1 = 680u"m/s",
           T_1 = 288u"K",
           P_1 = 1u"atm",
           gamma_1 = 1.4,
           a_2 = 441.79u"m/s"
       )
end
#broken bc overconstraint_validation takes forever
@test_skip test_normal_shock_inf_or_nan() isa NormalShockProperties
############################################################################
#internal coherence tests/unit tests

function test_internal_solver()
    Propulsion.internal_solver(ThermodynamicProperties, Dict(:P => 1.0, :T => 10.0), Dict(:z => 2.0))
end

#just test it built something
@test test_internal_solver() isa ThermodynamicProperties

function test_nozzle_flow_from_dict()
    NozzleFlowProperties(Dict(
        :F => 2u"N",
        :P_1 => 1u"Pa",
        :T_1 => 1u"K",
        :z_1 => 1u"mol/m^3",
        :rho_1 => 1u"kg/m^3",
        :MM_1 => 1u"kg/mol",
        :R_1 => 1u"J/kg/K",
        :a_1 => 1u"m/s",
        :cp_1 => 1u"J/kg/K",
        :cv_1 => 1u"J/kg/K",
        :gamma_1 => 1.4,
        :M_1 => 1.1,
        :a0_1 => 1u"m/s",
        :P0_1 => 1u"Pa",
        :rho0_1 => 1u"kg/m^3",
        :T0_1 => 1u"K",
        :v_1 => 1u"m/s",
        :A_1 => 1u"m^2",
        :Astar_1 => 1u"m^2",
        :mdot_1 => 1u"kg/s",
        :P_2 => 1u"Pa",
        :T_2 => 1u"K",
        :z_2 => 1u"mol/m^3",
        :rho_2 => 1u"kg/m^3",
        :MM_2 => 1u"kg/mol",
        :R_2 => 1u"J/kg/K",
        :a_2 => 1u"m/s",
        :cp_2 => 1u"J/kg/K",
        :cv_2 => 1u"J/kg/K",
        :gamma_2 => 1.4,
        :M_2 => 1.1,
        :a0_2 => 1u"m/s",
        :P0_2 => 1u"Pa",
        :rho0_2 => 1u"kg/m^3",
        :T0_2 => 1u"K",
        :v_2 => 1u"m/s",
        :A_2 => 1u"m^2",
        :Astar_2 => 1u"m^2",
        :mdot_2 => 1u"kg/s"
    ))
end

#testing it didn't error
@test test_nozzle_flow_from_dict() isa NozzleFlowProperties

function test_residue_unit_coherence()
    residues(NozzleFlowProperties(Dict(
        :F => 2u"N",
        :P_1 => 1u"Pa",
        :T_1 => 1u"K",
        :z_1 => 1u"mol/m^3",
        :rho_1 => 1u"kg/m^3",
        :MM_1 => 1u"kg/mol",
        :R_1 => 1u"J/kg/K",
        :a_1 => 1u"m/s",
        :cp_1 => 1u"J/kg/K",
        :cv_1 => 1u"J/kg/K",
        :gamma_1 => 1.4,
        :M_1 => 1.1,
        :a0_1 => 1u"m/s",
        :P0_1 => 1u"Pa",
        :rho0_1 => 1u"kg/m^3",
        :T0_1 => 1u"K",
        :v_1 => 1u"m/s",
        :A_1 => 1u"m^2",
        :Astar_1 => 1u"m^2",
        :mdot_1 => 1u"kg/s",
        :P_2 => 1u"Pa",
        :T_2 => 1u"K",
        :z_2 => 1u"mol/m^3",
        :rho_2 => 1u"kg/m^3",
        :MM_2 => 1u"kg/mol",
        :R_2 => 1u"J/kg/K",
        :a_2 => 1u"m/s",
        :cp_2 => 1u"J/kg/K",
        :cv_2 => 1u"J/kg/K",
        :gamma_2 => 1.4,
        :M_2 => 1.1,
        :a0_2 => 1u"m/s",
        :P0_2 => 1u"Pa",
        :rho0_2 => 1u"kg/m^3",
        :T0_2 => 1u"K",
        :v_2 => 1u"m/s",
        :A_2 => 1u"m^2",
        :Astar_2 => 1u"m^2",
        :mdot_2 => 1u"kg/s"
    )))
end

#testing it didn't error
@test test_residue_unit_coherence() isa Vector

function test_adjacency_list()
    adjacency_list(MassProperties)
end

@test test_adjacency_list() == [
    (:P, :T), 
    (:P, :z), 
    (:z, :T), 
    (:z, :rho), 
    (:z, :MM), 
    (:MM, :R), 
    (:MM, :rho)
]

function test_find_clique_1_var()
    clique_res = find_clique(MassProperties, [:P, :z, :MM], 1)
    any(clique_res.clique_equations .∈ [[1], [2], [3]]) &&
    any(clique_res.clique_vars .∈ [[:T], [:R], [:rho]]) &&
    clique_res.diagnostic == Propulsion.CliqueFound
end

@test test_find_clique_1_var()

function test_find_clique_2_var()
    clique_res = find_clique(CalorificProperties, [:P, :R, :gamma], 2)
    clique_res.clique_equations == [4, 5] &&
    Set(clique_res.clique_vars) == Set([:cp, :cv])
end

@test test_find_clique_2_var()

#############################################################################
# correctness
function test_air_molar_mass()
    MassProperties(P = 1u"bar", T = 300u"K", R = 287u"J/kg/K")
end

@test isapprox(test_air_molar_mass().MM, 28.9u"g/mol", atol=0.1u"g/mol")

#anderson
function test_flow_properties()
    #example 7.6
    solution = FlowProperties(P = 1u"atm", T = 320u"K", v = 1000u"m/s", gamma = 1.4, R = 287u"J/kg/K")
    #example 7.7
    solution2 = FlowProperties(P0 = 2220u"lbf/ft^2", P=1455.6u"lbf/ft^2", T = 483.04u"Ra", gamma = 1.4, R = 287u"J/kg/K")
    solution, solution2
end

@testset "examples 7.6 and 7.7" begin
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

@testset "example 10.1" begin 
    sol = test_example_10_1()
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

@testset "example 10.2" begin
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

function test_example_10_4()
    Quasi1dimflowProperties(
        P0=30u"atm",
        T0=3500u"K",
        Astar=0.4u"m^2",
        A=0.4u"m^2",
        gamma=1.22,
        R = 520u"J/kg/K"
    )
end

@testset "example 10.4" begin
    sol = test_example_10_4()
    @test isapprox(sol.rho0, 1.665u"kg/m^3", atol=10u"g/m^3")
    @test isapprox(sol.rho, 1.036u"kg/m^3", atol=10u"g/m^3")
    @test isapprox(sol.T, 3154u"K", atol=1u"K")
    @test isapprox(sol.a, 1415u"m/s", atol=1u"m/s")
    @test isapprox(sol.mdot, 586.4u"kg/s", atol=2u"kg/s")
end

function test_example_8_11()
    nsp = NormalShockProperties(
        v_1 = 680u"m/s",
        T_1 = 288u"K",
        P_1 = 1u"atm",
        gamma_1 = 1.4,
        R_1 = 287u"J/kg/K"
    )
end

@testset "example 8.11" begin 
    nsp = test_example_8_11()

    @test isapprox(nsp.P_2, 4.5u"atm", atol=0.1u"atm")
    @test isapprox(nsp.T_2, 486u"K", atol=1u"K")
    @test isapprox(nsp.v_2, 255u"m/s", atol=1u"m/s")
    @test ustrip(nsp.a_2 ) > 0
    @test isapprox(nsp.M_2, 0.577, atol = 1e-3)
end
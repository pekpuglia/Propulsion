module Propulsion
#todo
#pipe flow
#prioritize using intensive properties?
#sub/supersonic choice - initial guesses depend on the variables
#more tests
#normal shocks?
#actually use this - do As/Ac x Pt/Ps plot
#discober when initial guess for M is required
#wrap and use supersonic::Bool?
#add de Laval Nozzle - multiple sections
# q1dparams = Quasi1dimflowProperties(P=1e5, T=10.0, rho = 2.0, gamma = 1.4, Astar = 0.85, A = 1.0, initial_M=5)
# q1dparams.M
# res = Quasi1dimflowProperties(P = 1, T = 300, R = 287, gamma = 1.4, M=1.5, mdot=1)
# res = Quasi1dimflowProperties(P = 100, T = 300, R = 287, gamma = 1.4, M=1.5, mdot=1)
##
# Quasi1dimflowProperties(P = 1u"bar", R = 287u"J/kg/K", gamma = 1.4, mdot = 1u"kg/s", T = 300u"K", A=0.3u"m^2", initial_M = 10).M

include("physical_properties_base.jl")
include("physical_properties.jl")
include("solver.jl")

end # module propulsion

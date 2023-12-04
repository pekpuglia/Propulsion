module Propulsion
#todo
#pipe flow
#prioritize using intensive properties?

#sub/supersonic choice - initial guesses depend on the variables
#discober when initial guess for M is required
#wrap and use supersonic::Bool?

#more tests

#normal shocks?

#make "default values" for initial_guesses

#actually use this - do As/Ac x Pt/Ps plot

#add de Laval Nozzle - multiple sections

#FlowProperties(P=1u"atm", MM=29u"g/mol", rho=1.225u"kg/m^3", M = 1.8, T=300u"K", initial_gamma = 1.4) - singular exception
include("physical_properties_base.jl")
include("physical_properties.jl")
include("solver.jl")

end # module propulsion

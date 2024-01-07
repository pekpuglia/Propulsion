module Propulsion
#todo
#pipe flow
#prioritize using intensive properties?

#sub/supersonic choice - initial guesses depend on the variables
#discober when initial guess for M is required
#wrap and use supersonic::Bool?
#add de Laval Nozzle - multiple sections
#make "default values" for initial_guesses

#more tests

#make "default values" for initial_guesses

#actually use this - do As/Ac x Pt/Ps plot

#add de Laval Nozzle - multiple sections

#nsp = NormalShockProperties(
#     v_1 = 680u"m/s",
#     T_1 = 288u"K",
#     P_1 = 1u"atm",
#     gamma_1 = 1.4,
#     a_2 = 441.79u"m/s"
# ) - matrix contains Infs or NaNs

include("physical_properties_base.jl")
include("physical_properties.jl")
include("solver.jl")

end # module propulsion

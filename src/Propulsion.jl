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
#make "default values" for initial_guesses

include("physical_properties_base.jl")
include("physical_properties.jl")
include("solver.jl")

end # module propulsion

module Propulsion
#todo
#prioritize using intensive properties?

#sub/supersonic choice - initial guesses depend on the variables
#discover when initial guess for M is required
#wrap and use supersonic::Bool?
#add de Laval Nozzle - multiple sections
#make "default values" for initial_guesses

#more tests

#make "default values" for initial_guesses

#actually use this - do As/Ac x Pt/Ps plot

#pipe flow
#add de Laval Nozzle - multiple sections

include("physical_properties_base.jl")
include("physical_properties.jl")
include("solver.jl")

end # module propulsion

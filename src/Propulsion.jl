module Propulsion
#todo

#https://github.com/ksil/LFPSQP.jl
#use kg intead of g (speed scaling)

#prioritize using intensive properties?
#make all equations unitless

#sub/supersonic choice - initial guesses depend on the variables
#discover when initial guess for M is required
#wrap and use supersonic::Bool?


#make "default values" for initial_guesses

#actually use this - do As/Ac x Pt/Ps plot

#pipe flow
#add de Laval Nozzle - multiple sections

#unit and interval refactor!
#refactor find_clique
#add ReTest
#add more find_clique tests

include("physical_properties_base.jl")
include("physical_properties.jl")
include("solver.jl")

end # module propulsion

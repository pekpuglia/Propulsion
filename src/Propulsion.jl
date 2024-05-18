module Propulsion
#todo

#make "default values" for initial_guesses?
#unit (use better units)

#pipe flow
#add de Laval Nozzle - multiple sections
#actually use this - do As/Ac x Pt/Ps plot
#sub/supersonic choice - initial guesses depend on the variables
#discover when initial guess for M is required
#wrap and use supersonic::Bool?

#interval refactor! - allow negative force
#maybe this solves NormalShock?

#add ReTest

#refactor find_clique - search only subgraphs
#add more find_clique tests
#inject clique guess

include("physical_properties_base.jl")
include("physical_properties.jl")
include("solver.jl")

end # module propulsion

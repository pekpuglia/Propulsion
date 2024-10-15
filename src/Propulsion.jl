module Propulsion
#todo

#unit (use better units)
#test simple P0_1 T_1 A_1 R_1 gamma_1 A_2 P_2
#remove double constructors? - better decouple units

#actually use this - do As/Ac x Pt/Ps plot
#discover when initial guess for M is required
#wrap and use supersonic::Bool?

#variables have:
#initial_guess
#interval
#units

#maybe this solves NormalShock?
#add entropy

#add de Laval Nozzle - multiple sections

#add ReTest


include("physical_properties_base.jl")
include("physical_properties.jl")
include("solver.jl")

end # module propulsion

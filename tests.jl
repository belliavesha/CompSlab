using DelimitedFiles

include("./cs.jl")
using .IsothermalComptonAtmosphere: init_atmosphere, compute_slab
# println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.01,0.02,0.5,0.5,Theta))
# println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.1,0.01,0.5,0.05,Theta))
# println("Angular symmetry: ",CheckAngularSymmetry(Compton_redistribution_aa,0.02,0.02,-0.4,0.5))


IsothermalComptonAtmosphere.init_τ(50, 0.6)
IsothermalComptonAtmosphere.init_x(150)
IsothermalComptonAtmosphere.init_μ(9)
IsothermalComptonAtmosphere.set_ScatterNum(12)
init_atmosphere()
println(IsothermalComptonAtmosphere.mu)
I = compute_slab()

open("I_150_9_50.txt", "a") do io
            writedlm(io, I[:,:,1])
        end



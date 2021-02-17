
# print(a)

include("./cs.jl")
using .IsothermalComptonAtmosphere: init_atmosphere, compute_slab
# println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.01,0.02,0.5,0.5,Theta))
# println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.1,0.01,0.5,0.05,Theta))
# println("Angular symmetry: ",CheckAngularSymmetry(Compton_redistribution_aa,0.02,0.02,-0.4,0.5))


IsothermalComptonAtmosphere.init_τ(20)
IsothermalComptonAtmosphere.init_x(50)
# IsothermalComptonAtmosphere.init_τ(20,0.6 )
IsothermalComptonAtmosphere.init_μ(6)
IsothermalComptonAtmosphere.set_ScatterNum(18)
init_atmosphere()
println(IsothermalComptonAtmosphere.mu)
I = compute_slab()

println(I[:,10,1])


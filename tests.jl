
# print(a)

@time include("./cs.jl")
using .IsothermalComptonAtmosphere: init_atmosphere, compute_slab,collect_args
# println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.01,0.02,0.5,0.5,Theta))
# println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.1,0.01,0.5,0.05,Theta))
# println("Angular symmetry: ",CheckAngularSymmetry(Compton_redistribution_aa,0.02,0.02,-0.4,0.5))

# @time IsothermalComptonAtmosphere.s_functions(1)
# @time IsothermalComptonAtmosphere.s_functions(1)




begin
    IsothermalComptonAtmosphere.init_τ(2) # 20
    IsothermalComptonAtmosphere.init_x(2) #50
    # IsothermalComptonAtmosphere.init_τ(20,0.6 )
    IsothermalComptonAtmosphere.init_μ(2)
    IsothermalComptonAtmosphere.set_ScatterNum(1)
    # @time init_atmosphere()
    @time init_atmosphere()
    # println(IsothermalComptonAtmosphere.mu)
    println("the first call, the function compiles:")
    @time I = compute_slab()
    println("the next call, the function has been compiled, works a thousand times faster:")
    @time I = compute_slab()
    # @time I = compute_slab_impure()
    # @time I = compute_slab_impure()

    "done"
end
begin
    IsothermalComptonAtmosphere.init_τ(20) # 20
    IsothermalComptonAtmosphere.init_x(50) #50
    # IsothermalComptonAtmosphere.init_τ(20,0.6 )
    IsothermalComptonAtmosphere.init_μ(10)
    IsothermalComptonAtmosphere.set_ScatterNum(2)
    # @time init_atmosphere()

    @time init_atmosphere()
    # println(IsothermalComptonAtmosphere.mu)
    @time args = collect_args()
    @time I = compute_slab(args...)
    # @time I = compute_slab_impure()

    "done"

end

# # println(I[:,10,1])
# # [1 2; 3 4]

show(IsothermalComptonAtmosphere.to)


# @time IsothermalComptonAtmosphere.s_functions(1)
# @time IsothermalComptonAtmosphere.s_functions(1)

# show(IsothermalComptonAtmosphere.to)


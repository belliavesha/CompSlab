# print(a)

@time include("./csol.jl")
using .IsothermalComptonAtmosphere: init_atmosphere, compute_slab,CheckAngularSymmetry,CheckFrequencySymmetry,Compton_redistribution_aa
    
# IsothermalComptonAtmosphere.init_Θ()
# IsothermalComptonAtmosphere.Theta
# begin
#     println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.01,0.02,0.5,0.5))
#     println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.1,0.01,0.5,0.05))
#     # println("Angular symmetry: ",CheckAngularSymmetry(Compton_redistribution_aa,0.02,0.02,-0.4,0.5))
# end
# @time IsothermalComptonAtmosphere.s_functions(1)
# @time IsothermalComptonAtmosphere.s_functions(1)




begin
    IsothermalComptonAtmosphere.init_τ(2,1.0) # 20
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
    IsothermalComptonAtmosphere.init_Θ(0.1,0.002)
    IsothermalComptonAtmosphere.init_τ(50,1.0) # 20
    IsothermalComptonAtmosphere.init_x(100) #50
    # IsothermalComptonAtmosphere.init_τ(20,0.6 )
    IsothermalComptonAtmosphere.init_μ(9)
    IsothermalComptonAtmosphere.set_ScatterNum(12)
    # @time init_atmosphere()

    @time init_atmosphere()
    # println(IsothermalComptonAtmosphere.mu)
    # @time args = IsothermalComptonAtmosphere.collect_args()
    # @time I = compute_slab(args...)
    @time I = compute_slab()
    # @time I = compute_slab_impure()

    "done"

end



using PyPlot

# Intensity = zeros((NEnergy,NZenith,2)) #
begin 
    cols = ["tab:red","tab:orange","tab:green","tab:blue","tab:cyan"]

    fig, ax= plt.subplots(2,2, sharex = "col",sharey = "row",figsize=(7,7) ) 
    μ_grid=IsothermalComptonAtmosphere.μ_grid
    x_grid=IsothermalComptonAtmosphere.x_grid
    
    NEnergy, x, x_weight = x_grid
    NMu, NZenith, mu, mu_weight = μ_grid 
    
    # NEnergy, x, x_weight = x_grid
    # NMu, NZenith, mu, mu_weight = μ_grid 
    # NDepth, tau_T, tau, tau_weight = τ_grid
    # RedistributionMatrix, σ = R_grid

    J = I .*  1.235593147556e+20 # mc^2/hz
    for i in 1:5
        ax[1,1].plot(mu[NMu+1:end], x[27+6*i].*J[27+6*i,NMu+1:end,1],color=cols[i])
        ax[2,1].plot(mu[NMu+1:end],100 .* (I[27+6*i,NMu+1:end,2]./I[27+6*i,NMu+1:end,1]),color=cols[i])
        println(x[27+6*i]*500)
    end

    st = 24
    en = 64
    for i in 1:4
        ax[1,2].plot(x[st:en].*511,x[st:en] .* J[st:en,NMu+2*i,1],color=cols[i])
        ax[2,2].plot(x[st:en].*511,100 .* (I[st:en,NMu+2*i,2]./I[st:en,NMu+2*i,1]),color=cols[i])
    end 

    ax[1,2].set_yscale("log")
    ax[1,1].set_yscale("log")
    ax[1,2].set_xscale("log")
    ax[2,2].set_xscale("log")
    
    fig.tight_layout(h_pad=0.,w_pad=0.)
           
    fig.savefig("figure.pdf")
        
end 


# # println(I[:,10,1])
# # [1 2; 3 4]

show(IsothermalComptonAtmosphere.to)

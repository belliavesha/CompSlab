using Pkg
using Plots
gr()  # Initialize the GR backend for plotting

# List of packages used
const packages = ["SpecialFunctions", "FastGaussQuadrature", "LinearAlgebra", "TimerOutputs", "StaticArrays", "DelimitedFiles", "Plots", "PyPlot"]

# Function to check and add missing packages
function ensure_packages_installed()
    installed_packages = keys(Pkg.dependencies())
    for pkg in packages
        if !(pkg in installed_packages)
            println("Installing ", pkg)
            Pkg.add(pkg)
        end
    end
end

# Check and install packages
ensure_packages_installed()

println("Hello, world!")

using DelimitedFiles

include("./cs.jl")
using .IsothermalComptonAtmosphere: init_atmosphere, compute_slab
using .IsothermalComptonAtmosphere: set_absorption
set_absorption(0.1)

IsothermalComptonAtmosphere.init_Θ(0.2, 0.00002)
IsothermalComptonAtmosphere.init_τ(50, 3.0)
IsothermalComptonAtmosphere.init_x(100, -6.3, 0.3)
IsothermalComptonAtmosphere.init_μ(9)
IsothermalComptonAtmosphere.set_ScatterNum(30)
init_atmosphere()
println(IsothermalComptonAtmosphere.μ_grid)



using .IsothermalComptonAtmosphere: init_atmosphere, compute_slab,CheckFrequencySymmetry,CheckAngularSymmetry,Compton_redistribution_aa,Theta
println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.01,0.02,0.5,0.5))
println("Energy symmetry: ",CheckFrequencySymmetry(Compton_redistribution_aa,0.1,0.01,0.5,0.05))
println("Angular symmetry: ",CheckAngularSymmetry(Compton_redistribution_aa,0.02,0.02,-0.4,0.5))


I = compute_slab()

open("I_corona_100_9_50_02-5_3_fill_abs01.txt", "a") do io
    writedlm(io, I[:,:,1])
end
# Plot the spectrum I(x) and polarization degree for different mu
function plot_spectrum(I, mu, x)
    # Create a 2x1 subplot layout
    p = plot(layout=(2,1), size=(800,1000))

    # Plot intensity spectrum
    x_min, x_max = minimum(x), maximum(x)
    I_min, I_max = 0.0001, maximum(I[:, :, 1])
    
    subplot = plot!(p[1], x, I[:, 1, 1], xscale=:log10, yscale=:log10, 
                   legend=false, xlims=(x_min, x_max), ylims=(I_min, I_max))
    for i in 2:length(mu)
        plot!(p[1], x, I[:, i, 1])
    end
    xlabel!(p[1], "x")
    ylabel!(p[1], "I_x") 
    title!(p[1], "Spectrum")

    # Plot polarization degree
    pol = @. abs(I[:,:,2]/I[:,:,1])
    plot!(p[2], x, pol[:,1], xscale=:log10, legend=false, xlims=(x_min, x_max))
    for i in 2:length(mu)
        plot!(p[2], x, pol[:,i])
    end
    xlabel!(p[2], "x")
    ylabel!(p[2], "Polarization degree")
    title!(p[2], "Polarization")

    display(p)
end

# Assuming x and mu are defined and I is computed
x = IsothermalComptonAtmosphere.x_grid[2]  # Use x from x_grid
mu = IsothermalComptonAtmosphere.μ_grid[3]  # Example mu values
plot_spectrum(I, mu, x)

display(plot)

savefig("spectrum_plot_3_fill.png")

using ArgParse
using FITSIO

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--te"
            help = "give T_e"
            arg_type = Float64
            default = 0.08
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
#set values
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
    eval(Meta.parse("$(arg)=$(val)"))
end 

T__e = te/100.0
# T__e = 0.05

include("./cs.jl")

using .IsothermalComptonAtmosphere: init_atmosphere, compute_slab


t__bb = range(0.001, step=0.0001, stop=0.002) 
tau__T = range(0.5, step=0.1, stop=4.0) 

IsothermalComptonAtmosphere.init_x(150) 
IsothermalComptonAtmosphere.init_μ(9) 
IsothermalComptonAtmosphere.set_ScatterNum(12) #20? orders of scattering icluded, 
init_atmosphere()

f = FITS("CompSlab_$T__e.fits", "r+")

for ii in tau__T
    for iii in t__bb
        IsothermalComptonAtmosphere.init_τ(50, ii) #here come tau(nubmer of optical depth levels) and then tau_t, Thomson optical depth of thermalization
        IsothermalComptonAtmosphere.init_Θ(T__e, iii) # theta=kTe/mec2 and t=kTbb/mec2
        I = compute_slab()
        data = Dict("I"=>I[:,10:end,1], "Q"=>I[:,10:end,2])
        write(f,data)
    end
end

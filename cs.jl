
module IsothermalComptonAtmosphere

    include("./sigma.jl")

    using SpecialFunctions: besselk
    using FastGaussQuadrature
    using LinearAlgebra
    using .S_Functions: σ_Maxwell

    ScatterNum = 10 # total number of scatterings
    global NGamma= 10# number of Lorenz factor points (\gamma)
    global NAzimuth= 10 # 12 # numbers of azimuth angles (\phi) [0,pi]
    NEnergy = 20# 200 # 50# 101 # number of energy points (x)
    NDepth = 10 # 101  # number of optical depth levels (\tau)
    NMu = 5 # 20# 15 # number of propagation zenith angle cosines (\mu) [0,1]
    NZenith = 2*NMu # number of propagation zenith angles (z) [0,pi]
    
    # Atmosphere parameters: 
    tau_T= .6 # Thomson optical depth of thermalization 
    x_l, x_u = -3.7 , .3 # lower and upper bounds of the log_10 energy span
    Theta = 0.08  # dimensionless electron gas temperature (Theta = k T_e / m_e c^2) # it's about 0.1 
    T = 0.002 # 10/evere #  dimensionless photon black body temperature T = k T_bb / m_e c^2
    
    x = 10 .^ ( range(x_l,stop=x_u,length=NEnergy))
    x_weight = log(x[2]/x[1])
    
    mu = zeros(NZenith)
    mu_weight = zeros(NZenith)
    m2,mw = gausslegendre(NMu)
    mu[1:NMu] = (m2 .- 1)./2
    mu[NMu+1:2NMu] = (m2 .+ 1)./2
    mu_weight[1:NMu] = (mw)./2
    mu_weight[NMu+1:2NMu] = (mw)./2
    println(mu)
    println(mu_weight)
    tau = range(0, stop=tau_T, length=NDepth)
    tau_weight = tau_T/NDepth

    K2Y = besselk(2,1/Theta) # second modified Bessel function of reversed dimensionless temperature       
    

    σ = zeros(NEnergy)
    for e in 1:NEnergy
        σ[e] = σ_Maxwell(x[e],Theta)[1]
        println(σ[e],"   ",x[e])
    end
    # """   Planck function for Intensity of black body radiation 
    # The only argument x is the energy of a photon in units of electron rest energy ( hν/ m_e c^2 ) 
    # The photon temperature is given by T also in units of electron rest mass
    # Planck returns the intensity of  BB radiation
    # """
    function Planck(x)
        e = x/T
        C = 2*6.6261e-27/2.9979245e-10^2/4.135666e-18^3#/1.6021773e-9
        #C = 2h/c^2/[kev]^3
        C*e^3/expm1(e) 
    end

    # function Delta(x)
    #       C=2e4
    #       I=C*exp(-1e2*(x-T)^2/T/T)
    #       return I   
    # end   

    function sigma_cs(x) # not averaged on electron distribution 
        #   """ This function compute the Compton scattering cross-section in electron rest frame 
        #   x is the energy of the scattering photon in units of electron rest energy
        #   this function approaches the mean compton cross-section when electron gas temperature is small
        #   """ 
        if x<.1
            a = 3/8
            n = 0
            s = 0
            while (abs(a)*(n+2)^2 > 1e-11) # Taylor series sum of the formula below
                s=s+a*(n+2+2/(n+1)+8/(n+2)-16/(n+3))
                n=n+1
                a=-2*x*a
            end
            s
        else 
            3(2 - (1/x + 1 - x/2)log(1 + 2x))/4x^2 + 3(1+x) / 4(1+2*x)^2
        end
    end



    function Compton_redistribution_m(x1,x2,mu,gamma)
    #   """   Compton redistribution matrix for monoenergetic electron gas
    #   The arguements are:
    #   x1 and x2 - photon energies in units of electron rest energy ( h \\nu / m_e c^2 ) 
    #   mu - cosine of a scattering angle 
    #   gamma - energy of each electron in the gas in units of the electron rest mass
        
    #   This fuctions returns (R,RI,RQ,RU)
    #   which are scalar redistribution functions for isotropic monoenergetic gas
    #   and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectively 
    #   R44 or RV is also not equal to zero but we never need it 
    #   """

        #the next variables' names are adopted from J. Poutanen & O. Vilhu 1993
        r = ( 1 + mu )/( 1 - mu )
        a1 = √( (gamma-x1)^2 + r )
        a2 = √( (gamma+x2)^2 + r )
        v = a1*a2
        u = a2-a1  #(x1+x2)*(2.*gamma+x2-x1)/(a1+a2) # the formulas probably give always the same numbers 
        
        q = x1*x2*( 1 - mu )
        Q = √( x1*x1 + x2*x2 - 2*x1*x2*mu ) # √( (x1-x2)^2 +2.*q ) # the formula probably gives the same also
        gammaStar = ( x1 - x2 + Q*√( 1 + 2/q ) )/2

        # print(  gamma-gammaStar, gammaStar,gamma,Q,u,(u-Q)/(Q+u))

        if gamma < gammaStar 
            println(" Compton_redistribution_m got gamma < gammaStar " )
            # I belive, in the case fucntion just won't ever be called
            return  (0.,0.,0.,0.)
        else

            Ra = u*( u*u - Q*Q )*( u*u + 5*v )/2/q/q/v/v/v + u*Q*Q/q/q/v/v
            Rb = 2/Q + u/v*( 1 - 2/q )
            Rc = u/q/v*( ( u*u - Q*Q )/r/q - 2)
            Rd = 2/Q + 2*(u-Q)/r/q*((u-Q)/r/q*(2*Q+u) - 4) + 2*u/v/q 

            R = Ra + Rb
            RI = Ra + Rc
            RU = Rd + 2Rc
            RQ = RU + Ra
            
            #print(r,v,u,q,Q,gammaStar,Ra,Rb,Rc,R)
            #print(x1,x2,mu,gamma,R,RI,RQ,RU)
            
            return (R,RI,RQ,RU)
        end
    end


    function Maxwell_r(gamma)
        # """The normalized relativistic Maxwellian distribution
        # the density of particles in the dimensionless momentum volume (4 \pi z^2 dz) is nomalized to unity
        # Theta is the dimensionless electron gas temperature (Theta = k * T_e / m_e c^2)
        # gamma is electron energy in units of the electron rest mass
        # The fuction returns the momentum dencity value ( f(\gamma) )
        # """
        r = 0.25/pi/Theta*exp(-gamma/Theta)/K2Y
        return r
    end

    function Compton_redistribution(x1,x2,mu) # if distribution is not Maxwellian the function must be modified.
    #   """    Thermal Compton redistribution matrix (integrated with electron distribution function)
    #   And the distribution is maxwellian (if it's not the function must be modified)
    #   The arguments are:
    #   x1 and x2 - photon energies in units of electron rest energy ( h \\nu / m_e c^2 ) 
    #   mu - cosine of a scattering angle 
      
    #   This fuctions returns (R,RI,RQ,RU)
    #   which are scalar redistribution functions for Maxwellian relativistic gas
    #   and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectively 
    #   R44 or RV is also not equal to zero but we never need it  
    #   """
        q = x1*x2*(1 - mu)
        Q = √( x1*x1 + x2*x2 - 2*x1*x2*mu )
        gammaStar = (x1-x2+Q*√( 1 + 2/q ) )/2 # lower bound of integration 
        C=3/8*Theta*Maxwell_r(gammaStar)
        gamma, gamma_weight = gausslaguerre(NGamma) 
        
        R=zeros(4)
        for i in 1:NGamma
            T=Compton_redistribution_m(x1,x2,mu,Theta*gamma[i]+gammaStar)
            for j in 1:4
                R[j] += C*gamma_weight[i]*T[j]
            end
        end
        R
    end
          
    function Compton_redistribution_aa(x1,x2,mu1,mu2)
        # """   Azimuth-avereged Compton redistribution matrix 
        # for computing of electron scattering source function 
        # The arguements are:
        # x1 and x2 are photon energies in units of electron rest energy ( h \\nu / m_e c^2 ) 
        # mu1 and mu2 are cosines of angles between photon propagation directions and fixed direction
        
        # This function returns R11 R12 R21 R22 matrix elements
        # We need only 2x2 matrix in the upper left corner of the general matrix,
        # becouse U and V components on the Stokes vector are zero in this case.
        # """
        
        eta1 = 1 - mu1*mu1  # squared sinuses of the angles 
        eta2 = 1 - mu2*mu2  
        
        # phi, phi_weight = IntAzimuth
        phi, phi_weight = gausslegendre( NAzimuth*2 );

        phi = pi .* phi
        az_c = cos.(phi)    # array of azimuth cosines
        az_s = sin.(phi) .^ 2  # array of azimuth square sinuses
        sc_c = (mu1*mu2) .- sqrt(eta1*eta2) .* az_c # array of scattering angles' cosines
        sc_s = 1 .- sc_c .^ 2 # array of scattering angles' squared sinuses
        cos2chi1 = 2 .* (mu1 .* sc_c .- mu2) .* (mu1 .* sc_c .- mu2) ./ eta1 ./ sc_s .- 1  # array[ cos( 2 \chi_1 ) ]
        cos2chi2 = 2 .* (mu1 .- mu2 .* sc_c) .* (mu1 .- mu2 .* sc_c) ./ eta2 ./ sc_s .- 1  # array[ cos( 2 \chi_2 ) ]
        sin2chiP = 4 .* (mu1 .- mu2 .* sc_c) .* (mu1 .* sc_c .- mu2) .* az_s ./ sc_s .^ 2  # array[ sin( 2 \chi_1 )*sin( 2 \chi_2 ) ]

        R=zeros( (2,2) )
        for i in 1:NAzimuth*2
            (C,I,Q,U)=Compton_redistribution(x1,x2,sc_c[i])
            R[1,1] += C*pi*phi_weight[i]
            R[1,2] += I*pi*cos2chi2[i]*phi_weight[i]
            R[2,1] += I*pi*cos2chi1[i]*phi_weight[i]
            R[2,2] += pi*(Q*cos2chi1[i]*cos2chi2[i]+U*sin2chiP[i])*phi_weight[i]         
        end 
        # print(x1,x2,mu1,mu2,R)
        return R .* (x1*x1/x2)
    end

    function CheckAngularSymmetry(r,x1,x2,mu1,mu2)
        eps=1e-10
        one = r(x1,x2,mu1,mu2)
        two = r(x1,x2,mu2,mu1)
        three = r(x1,x2,-mu1,-mu2)
        four = r(x1,x2,-mu2,-mu1) 
        v=true
        a,b,c = two[1,1]/one[1,1]-1,three[1,1]/one[1,1]-1,four[1,1]/two[1,1]-1 # 0 0 0
        println(a," ",b," ",c)
        if abs(a)>eps || abs(b)>eps || abs(c)>eps 
            v=false
        end 
        a,b,c = two[2,1]/one[1,2]-1,three[2,1]/one[2,1]-1,four[1,2]/two[1,2]-1 # 0 0 0
        println(a," ",b," ",c)
        if abs(a)>eps || abs(b)>eps || abs(c)>eps 
            v=false
        end
        a,b,c = two[1,2]/one[2,1]-1,three[1,2]/one[1,2]-1,four[2,1]/two[2,1]-1 # 0 0 0
        println(a," ",b," ",c)
        if abs(a)>eps || abs(b)>eps || abs(c)>eps 
            v=false
        end
        a,b,c = two[2,2]/one[2,2]-1,three[2,2]/one[2,2]-1,four[2,2]/two[2,2]-1 # 0 0 0
        println(a," ",b," ",c)
        if abs(a)>eps || abs(b)>eps || abs(c)>eps 
            v=false
        end
        return v

    end 
# # frequency symmetry: CHECK [v]
    function CheckFrequencySymmetry(r,x1,x2,mu1,mu2,Theta)
        eps =1e-10
        one = r(x1,x2,mu1,mu2)
        two = r(x2,x1,mu1,mu2)
        v=true
        ratio  = x1^3 / x2^3 * exp((x2-x1)/Theta)
        a = abs( ratio*two[1,1]/one[1,1]-1)# 1e-15
        println(a)
        v=v && a < eps
        a = abs( ratio*two[2,1]/one[2,1]-1)# 1e-12
        println(a)
        v=v && a < eps
        a = abs( ratio*two[1,2]/one[1,2]-1)# 1e-12
        println(a)
        v=v && a < eps
        a = abs( ratio*two[2,2]/one[2,2]-1)# 1e-14
        println(a)
        v=v && a < eps
        v
    end


    function CRM()
        
        sigma=zeros(NEnergy)
        RedistributionMatrix = ones( (NEnergy,NEnergy,NZenith,NZenith,2,2) )
        percent=0.0
        for e in 1:NEnergy # x [-\infty,\infty]
            percent+=100/NEnergy
            println((percent))
            
            for e1 in e:NEnergy # x1 [x,\infty]
                # percent+=200/NEnergy/(NEnergy+1)
                # println((percent))
                for d in 1:NMu # mu [-1,0]
                    for d1 in d:NMu # mu1 [-1,mu]
                        md= NZenith-d+1 # -mu
                        md1=NZenith-d1+1 # -mu1
                        w=mu_weight[d1]*x_weight*mu_weight[d]
                        t=d1>d
                        f=e1>e

                        r=Compton_redistribution_aa(x[e],x[e1],mu[d],mu[d1])
                        rm=Compton_redistribution_aa(x[e],x[e1],mu[d],mu[md1])
                        sigma[e1]+=(r[1,1]+rm[1,1])*w
                        RedistributionMatrix[e,e1,d,d1,:,:]=r
                        RedistributionMatrix[e,e1,md,md1,:,:]=r
                        RedistributionMatrix[e,e1,d,md1,:,:]=rm
                        RedistributionMatrix[e,e1,md,d1,:,:]=rm
                        if f # frequency symmetry
                            m=exp((x[e]-x[e1])/Theta)*x[e1]^3/x[e]^3
                            rf=r .* m  # when Maxwellian or Wein distributions
                            rmf=rm .* m  
                            sigma[e]+=(rf[1,1]+rmf[1,1])*w
                            RedistributionMatrix[e1,e,d,d1,:,:]=rf
                            RedistributionMatrix[e1,e,md,md1,:,:]=rf
                            RedistributionMatrix[e1,e,d,md1,:,:]=rmf
                            RedistributionMatrix[e1,e,md,d1,:,:]=rmf
                        end
                        if t # angular symmetry
                            r[1,2],r[2,1] = r[2,1],r[1,2]
                            rm[1,2],rm[2,1] = rm[2,1],rm[2,1]
                            sigma[e1]+=(r[1,1]+rm[1,1])*w
                            RedistributionMatrix[e,e1,d1,d,:,:]=r
                            RedistributionMatrix[e,e1,md1,md,:,:]=r
                            RedistributionMatrix[e,e1,md1,d,:,:]=rm
                            RedistributionMatrix[e,e1,d1,md,:,:]=rm
                            if f # both symmeties 
                                rf=r .* m
                                rmf=rm .* m
                                sigma[e]+=(rf[1,1]+rmf[1,1])*w
                                RedistributionMatrix[e1,e,d1,d,:,:]=rf
                                RedistributionMatrix[e1,e,md1,md,:,:]=rf
                                RedistributionMatrix[e1,e,md1,d,:,:]=rmf
                                RedistributionMatrix[e1,e,d1,md,:,:]=rmf
                            end
                        end
                    end
                end
            end
        end
        println(RedistributionMatrix[2,2,2,2,:,:])
        RedistributionMatrix, sigma
    end




    function Slab() 

        RedistributionMatrix = CRM()[1]
        # Initializing Stokes vectors arrays, computiong scatterings 
        # Iin=Planck # Delta # initial photon distribution 
        Source = zeros((ScatterNum,NDepth,NEnergy,NZenith,2)) # source function                 
        Stokes = zeros((ScatterNum,NDepth,NEnergy,NZenith,2)) # intensity Stokes vector
        Stokes_out = zeros((ScatterNum+1,NEnergy,NZenith,2)) # outgoing Stokes vector of each scattering
        Stokes_in = zeros((NDepth,NEnergy,NZenith,2)) # Stokes vector of the initial raiation (0th scattering) 
        Intensity = zeros((NEnergy,NZenith,2)) # total intensity of all scattering orders from the slab suface 
        for e in 1:NEnergy
            for d in NMu+1:NZenith
                for t in 1:NDepth
                    Stokes_in[t,e,d,1]=Planck(x[e])*exp(-tau[t]*σ[e]/mu[d]) 
                end
                Stokes_out[1,e,d,1]=Planck(x[e])*exp(-tau_T*σ[e]/mu[d]) 
            end
        end

        for k in 1:ScatterNum # do ScatterNum scattering iterations
            for t in 1:NDepth # S_k= R I_{k-1}
                for e in 1:NEnergy
                    for d in 1:NZenith
                        S=zeros(2)  
                        for e1 in 1:NEnergy
                            for d1 in 1:NZenith
                                w = mu_weight[d1]*x_weight # total weight
                                r = RedistributionMatrix[e,e1,d,d1,:,:]  # 
                                if k>1 
                                    I = Stokes[k-1,t,e1,d1,:] 
                                else
                                    I = Stokes_in[t,e1,d1,:]
                                end
                                S[1] += w*( I[1]*r[1,1] + I[2]*r[1,2] ) # 
                                S[2] += w*( I[1]*r[2,1] + I[2]*r[2,2] ) #
                            end
                        end
                        Source[k,t,e,d,:] += S #     
                    end
                end
            end
            for t in 1:NDepth# I_k= integral S_k
                for e in 1:NEnergy 
                    for d in 1:NZenith
                        I = Source[k,t,e,d,:] .* (tau_weight/2)
                        if mu[d]>0
                            for t1 in 1:(t-1) #
                                S = Source[k,t1,e,d,:] #
                                I += tau_weight .* S .* exp(σ[e]*(tau[t1]-tau[t])/mu[d])
                            end
                            S = Source[k,1,e,d,:] #
                            I -= (tau_weight)/2 .*S .* exp(σ[e]*(-tau[t])/mu[d])
                        else
                            for t1 in (t+1):NDepth
                                S = Source[k,t1,e,d,:] #
                                I += tau_weight .* S .* exp(σ[e]*(tau[t1]-tau[t])/mu[d])
                            end
                            S = Source[k,NDepth,e,d,:] #
                            I -= (tau_weight/2) .* S .* exp(σ[e]*(tau_T-tau[t])/mu[d])
                        end
                        Stokes[k,t,e,d,:] += I ./ abs(mu[d]) #abs
                    end
                end
            end
            println("order ",k)
        end
                    
        Intensity += Stokes_out[1,:,:,:]
        for k in 1:ScatterNum
            Stokes_out[k+1,:,:,:] += Stokes[k,end,:,:,:]
            Intensity += Stokes[k,end,:,:,:]
        end
        Intensity
    end
end

module IsothermalComptonAtmosphere

    using SpecialFunctions: besselk
    using FastGaussQuadrature
    using LinearAlgebra

    function init_Θ(θ = 0.08, t = 0.002) # dimensionless photon black body temperature T = k T_bb / m_e c^2
        global Θ = θ
        global Theta = θ
        global T = t
        global K2Y = besselk(2,1/Theta) # second modified Bessel function of reversed dimensionless temperature       
    end; init_Θ()

    function init_τ(n=10, t = 0.6 )
        global NDepth = n # 101  # number of optical depth levels (\tau)
        global tau_T= t # Thomson optical depth of thermalization 
        global τ_T = t
        global tau = range(0, stop=tau_T, length=NDepth)
        global tau_weight = tau[2]-tau[1]
    end; init_τ();

    function init_μ(n = 5)
        global NMu = n # 20# 15 # number of propagation zenith angle cosines (\mu) [0,1]
        global NZenith = 2*NMu # number of propagation zenith angles (z) [0,pi]
        global mu = Array{Float64}(undef,NZenith)
        global mu_weight = Array{Float64}(undef,NZenith)
        m2,mw = gausslegendre(NMu)
        mu[1:NMu] = (m2 .- 1)./2
        mu[NMu+1:2NMu] = (m2 .+ 1)./2
        mu_weight[1:NMu] = (mw)./2
        mu_weight[NMu+1:2NMu] = (mw)./2
    end; init_μ();

    function init_x(n = 20, x_l = -3.7 , x_u = .3 )# lower and upper bounds of the log_10 energy span
        global NEnergy = n# 200 # 50# 101 # number of energy points (x)
        global x = 10 .^ ( range(x_l,stop=x_u,length=NEnergy))
        global x_weight = log(x[2]/x[1])
    end; init_x();

    function set_ScatterNum(n = 10)
        global ScatterNum = n # total number of scatterings
    end; set_ScatterNum()

    function init_g(ng = 30, nγ = 10, na = 10 )
        # global NGauss = ng
        # global NGamma= nγ # number of Lorenz factor points (\gamma)
        # global NAzimuth= na # 12 # numbers of azimuth angles (\phi) [0,pi]
        global glγ = gausslaguerre(nγ) 
        global glϕ = gausslegendre(na*2)
        global glξ = gausslaguerre(ng) 
    end; init_g()

    # for formulas see Nagirner and Poutanen 1994 (https://users.utu.fi/jurpou/papers/1994_single_cs.pdf)
    # C  S- functions. if K=0, then only  SS -- mean cross-section; 
    function s_functions(x)
        if x<0.25 # Asymptotic expansion
            s_0 = 0
            s_1 = 0
            s_2 = 0
            S_1 = 0
            S_2 = 0
            S_3 = 0
            S_4 = 0
            a = 1
            b = -2x
            n = 0
            while abs(4a)*(n+3)^3>eps(1.0)
                n1 = n + 1
                n2 = n + 2
                n3 = n + 3
                n4 = n + 4
                n5 = n + 5
                d1 = 1/n1
                d2 = 1/n2
                d3 = 1/n3
                d4 = 1/n4
                d5 = 1/n5
                s_0 += (n2 + 2d1 + 8d2 - 16d3)*a # (4.1.1)
                s_1 += ((n+5)n + 24d3)*a # (4.1.2)
                S_1 += (n1*n3 - 6 - 6d2 - 24d3 + 72d4 )*a # (4.1.3)
                S_2 += (n1*n4 - 6 - 12d3 - 72d4 + 144d5)*a
                s_2 += (0.25*n2*n3*n4 + 2 - n)*a # (4.1.3) but not literally
                S_3 += (0.25*n1*n3*n4 + 7 - n - 24d4)*a
                S_4 += (0.25*n1*n2*n4 + 12 -n + 6d3 + 24d4 - 96d5)*a
                a *= b
                n += 1
            end
            # println(x," n =  ",n)
            s_0 *= 0.375
            s_1 *= 0.125
            s_2 *= 0.125 
            S_1 *= 0.25
            S_2 *= 0.25
            S_3 *= 0.25
            S_4 *= 0.5
        else 
            R = 1/(1+2x) # (3.2.13)
            l = - log(R) # (3.2.13) 
            s_0 = 0.375/x^2 *(4 + (x - 2 - 2/x)*l + 2*(1+x)*(x*R)^2)  # (3.2.12)
            s_1 = 0.125/x^3 *(3l + 4x^2 - 4.5x - x*R*(1.5 + x*R^2)) # (3.2.14)
            s_2 = (1+(4+(7+4.5x)x)x)R^4 # (3.2.15)
            S_1 = (s_0 - s_1)/x  # (3.2.18) 
            S_2 = (s_1 - S_1)/x
            S_3 = (s_1 - s_2)/x
            S_4 = (S_1 - S_3)/x           
        end
        S_5 = 3S_4 - 4S_3     # (3.2.18) 
        S_7 = S_3 - 0.5S_4                                 
        S_6 = s_2 - 3S_7     
        return (s_0,s_1,s_2,S_1,S_2,S_3,S_4,S_5,S_6,S_7)    
    end

    function σ_Maxwell(x) # Averaging over Maxwellian distribution
        Y = 1/Θ # ΘY = 1 
        D = Θ*exp(-Y)/(2*K2Y) # basically, the factor befor the integrals
        σ = 0 # cross section, in units of Thomson cross-section;  
        X = 0 # mean energy of scattered photon;
        Q = 0 # dispersion energy of scattered photon;
        P = 0 # radiative pressure, 
        ξ, weight = glξ
        for l = 1:length(ξ)
            for j in [1,2] # substitutions for γ < 1 + Θ and γ > 1 + Θ
                if j == 1
                    g = 1+ξ[l] # substitution γ = 1 + Θ*(1 + ξ)
                    m = weight[l]*D/ℯ
                else
                    g = exp(-ξ[l]) # substitution γ = 1 + Θ*exp(-ξ)
                    m = weight[l]*D*exp(-g)
                end
                γ = 1+Θ*g
                z = √(Θ*g*(γ+1))
                for sig in [-1,1]  # u = γ + z and u = γ - z
                    u = γ + sig*z
                    s = s_functions(x*u) 
                    σ += u^2*s[1]*m/z # (3.4.2)
                    X += u^3*( (γ + Θ)*s[4] + x*s[5])*m/z # (3.4.2)
                    Q += u^2*(s[9] - u*( (γ + Θ)*s[8] + u*(s[10] - (γ^2 + 2γ*Θ+ 2Θ^2)*s[7]) ) ) *m/z # (3.4.2)
                    P += u^4*s[4]*m/z # (5.3.8)
                end
            end
        end
        σ, X, Q, P
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



    function Compton_redistribution_m(x1,x2,μ,γ)
    #   """   Compton redistribution matrix for monoenergetic electron gas
    #   The arguements are:
    #   x1 and x2 - photon energies in units of electron rest energy ( h \\nu / m_e c^2 ) 
    #   μ - cosine of a scattering angle 
    #   γ - energy of each electron in the gas in units of the electron rest mass
        
    #   This fuctions returns (R,RI,RQ,RU)
    #   which are scalar redistribution functions for isotropic monoenergetic gas
    #   and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectively 
    #   R44 or RV is also not equal to zero but we never need it 
    #   """

        #the next variables' names are adopted from J. Poutanen & O. Vilhu 1993
        r = ( 1 + μ )/( 1 - μ )
        a1 = √( (γ-x1)^2 + r )
        a2 = √( (γ+x2)^2 + r )
        v = a1*a2
        u = a2-a1  #(x1+x2)*(2.*γ+x2-x1)/(a1+a2) # the formulas probably give always the same numbers 
        
        q = x1*x2*( 1 - μ )
        Q = √( x1*x1 + x2*x2 - 2*x1*x2*μ ) # √( (x1-x2)^2 +2.*q ) # the formula probably gives the same also
        γStar = ( x1 - x2 + Q*√( 1 + 2/q ) )/2

        # print(  γ-γStar, γStar,γ,Q,u,(u-Q)/(Q+u))

        if γ < γStar 
            println(" Compton_redistribution_m got γ < γStar " )
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
            
            #print(r,v,u,q,Q,γStar,Ra,Rb,Rc,R)
            #print(x1,x2,μ,γ,R,RI,RQ,RU)
            
            return (R,RI,RQ,RU)
        end
    end


    function Maxwell_r(γ)
        # """The normalized relativistic Maxwellian distribution
        # the density of particles in the dimensionless momentum volume (4 \pi z^2 dz) is nomalized to unity
        # Theta is the dimensionless electron gas temperature (Theta = k * T_e / m_e c^2)
        # γ is electron energy in units of the electron rest mass
        # The fuction returns the momentum dencity value ( f(\γ) )
        # """
        r = 0.25/pi/Theta*exp(-γ/Theta)/K2Y
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
        γStar = (x1-x2+Q*√( 1 + 2/q ) )/2 # lower bound of integration 
        C=3/8*Theta*Maxwell_r(γStar)
        γ, weight = glγ
        
        R=zeros(4)
        for i in 1:length(γ)
            T=Compton_redistribution_m(x1,x2,mu,Theta*γ[i]+γStar)
            for j in 1:4
                R[j] += C*weight[i]*T[j]
            end
        end
        R
    end
        
    function Compton_redistribution_aa(x1,x2,μ1,μ2)
        # """   Azimuth-avereged Compton redistribution matrix 
        # for computing of electron scattering source function 
        # The arguements are:
        # x1 and x2 are photon energies in units of electron rest energy ( h \\nu / m_e c^2 ) 
        # μ1 and μ2 are cosines of angles between photon propagation directions and fixed direction
        
        # This function returns R11 R12 R21 R22 matrix elements
        # We need only 2x2 matrix in the upper left corner of the general matrix,
        # becouse U and V components on the Stokes vector are zero in this case.
        # """
        
        η1 = 1 - μ1*μ1  # squared sinuses of the angles 
        η2 = 1 - μ2*μ2  
        
        # phi, phi_weight = IntAziμth
        ϕ, weight = glϕ

        ϕ = pi .* ϕ
        az_c = cos.(ϕ)    # array of azimuth cosines
        az_s = sin.(ϕ) .^ 2  # array of azimuth square sinuses
        sc_c = (μ1*μ2) .- sqrt(η1*η2) .* az_c # array of scattering angles' cosines
        sc_s = 1 .- sc_c .^ 2 # array of scattering angles' squared sinuses
        cos2χ1 = 2 .* (μ1 .* sc_c .- μ2) .* (μ1 .* sc_c .- μ2) ./ η1 ./ sc_s .- 1  # array[ cos( 2 χ_1 ) ]
        cos2χ2 = 2 .* (μ1 .- μ2 .* sc_c) .* (μ1 .- μ2 .* sc_c) ./ η2 ./ sc_s .- 1  # array[ cos( 2 χ_2 ) ]
        sin2χP = 4 .* (μ1 .- μ2 .* sc_c) .* (μ1 .* sc_c .- μ2) .* az_s ./ sc_s .^ 2  # array[ sin( 2 χ_1 )*sin( 2 χ_2 ) ]

        R=zeros( (2,2) )
        for i in 1:length(ϕ)
            (C,I,Q,U)=Compton_redistribution(x1,x2,sc_c[i])
            R[1,1] += C*pi*weight[i]
            R[1,2] += I*pi*cos2χ2[i]*weight[i]
            R[2,1] += I*pi*cos2χ1[i]*weight[i]
            R[2,2] += pi*(Q*cos2χ1[i]*cos2χ2[i]+U*sin2χP[i])*weight[i]         
        end 
        # print(x1,x2,μ1,μ2,R)
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
        # sigma=zeros(NEnergy)
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
                        # sigma[e1]+=(r[1,1]+rm[1,1])*w
                        RedistributionMatrix[e,e1,d,d1,:,:]=r
                        RedistributionMatrix[e,e1,md,md1,:,:]=r
                        RedistributionMatrix[e,e1,d,md1,:,:]=rm
                        RedistributionMatrix[e,e1,md,d1,:,:]=rm
                        if f # frequency symmetry
                            m=exp((x[e]-x[e1])/Theta)*x[e1]^3/x[e]^3
                            rf=r .* m  # when Maxwellian or Wein distributions
                            rmf=rm .* m  
                            # sigma[e]+=(rf[1,1]+rmf[1,1])*w
                            RedistributionMatrix[e1,e,d,d1,:,:]=rf
                            RedistributionMatrix[e1,e,md,md1,:,:]=rf
                            RedistributionMatrix[e1,e,d,md1,:,:]=rmf
                            RedistributionMatrix[e1,e,md,d1,:,:]=rmf
                        end
                        if t # angular symmetry
                            r[1,2],r[2,1] = r[2,1],r[1,2]
                            rm[1,2],rm[2,1] = rm[2,1],rm[1,2]
                            # sigma[e1]+=(r[1,1]+rm[1,1])*w
                            RedistributionMatrix[e,e1,d1,d,:,:]=r
                            RedistributionMatrix[e,e1,md1,md,:,:]=r
                            RedistributionMatrix[e,e1,md1,d,:,:]=rm
                            RedistributionMatrix[e,e1,d1,md,:,:]=rm
                            if f # both symmeties 
                                rf=r .* m
                                rmf=rm .* m
                                # sigma[e]+=(rf[1,1]+rmf[1,1])*w
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
        RedistributionMatrix
    end

    function init_atmosphere()
        global RedistributionMatrix = CRM()
        global σ = Array{Float64}(undef,NEnergy)
        for e in 1:NEnergy
            σ[e] = σ_Maxwell(x[e])[1]
        end
    end


    function compute_slab() 
    
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
        # println(Stokes_out[1,:,:,1])
        Intensity += Stokes_out[1,:,:,:]
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

            Stokes_out[k+1,:,:,:] += Stokes[k,end,:,:,:]
            Intensity += Stokes[k,end,:,:,:]
            contribution, e = findmax(Stokes[k,end,:,end,1])
            println("order ",k, " contribution : ", contribution/Intensity[e,end,1])

        end
                    
        Intensity
    end

end

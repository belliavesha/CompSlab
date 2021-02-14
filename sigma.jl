

module S_Functions
    using FastGaussQuadrature: gausslaguerre
    using SpecialFunctions: besselk

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



    function σ_Maxwell(x,Θ) # Averaging over Maxwellian distribution
        Y = 1/Θ # ΘY = 1 
        K2Y =  besselk(2,Y) 
        D = Θ*exp(-Y)/(2*K2Y) # basically, the factor befor the integrals
        NG = 30
        ξ, w = gausslaguerre(NG) 
        σ = 0 # cross section, in units of Thomson cross-section;  
        X = 0 # mean energy of scattered photon;
        Q = 0 # dispersion energy of scattered photon;
        P = 0 # radiative pressure, 
        for l = 1:NG
            for j in [1,2] # substitutions for γ < 1 + Θ and γ > 1 + Θ
                if j == 1
                    g = 1+ξ[l] # substitution γ = 1 + Θ*(1 + ξ)
                    m = w[l]*D/ℯ
                else
                    g = exp(-ξ[l]) # substitution γ = 1 + Θ*exp(-ξ)
                    m = w[l]*D*exp(-g)
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
                                            
    # σ_Maxwell(0.25,0.08) 
    #  0.666953277      0.675436571       0.753782020       0.595729683
end                        

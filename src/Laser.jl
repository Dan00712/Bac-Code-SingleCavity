module Laser
using LinearAlgebra

using ForwardDiff

∂ = ForwardDiff.derivative

using ..Constants

export Ec, Et, αeq_c, δc, γ, Hs, ∂zHs, HHs

function E0c(Δω)
    # Δω = ωc - ω0
    ωc = Δω + ω0

    sqrt(ħ*ωc/2/ϵ0/Vc)
end


function Ec(x, y, z, Δω)
    kc = let
        ωc = Δω + ω0
        ωc/c
    end
    E0c(Δω) * exp(-(x^2+z^2)/Wc^2) * ζ * cos(kc*y)
end

function ϕt(x, y, z)
    -atan(z/zR) + k0*z/2 * (x^2+y^2)/(z^2+zR^2)
end

function W(z)
    Wt * sqrt(1+(z/zR)^2)
end

expi(z) = exp(im*z)

function Et(x, y, z)
    E0/2 * expi(k0*z+ϕt(x, y, z)) * Wt/W(z) * exp(-(y/Ay/W(z))^2) * exp(-(x/Ax/W(z))^2)
end

function δc(x, y, z, Δω)
    Δω - α/ħ * (Ec(x, y, z, Δω))^2
end

function αeq_c(x, y, z, Δω, κ)
    α/ħ * Et(x, y, z) * Ec(x, y, z, Δω)/(δc(x, y, z, Δω)-im*κ/2)
end



function Hs(x, y, z, Δω, κ)
    -abs2(
         Et(x,y,z) + αeq_c(x,y,z, Δω, κ)*Ec(x,y,z, Δω)
        )
end

function ∂zHs(z, Δω, κ)
    a = αeq_c(0,0,z, Δω, κ)
    ωc = ω0 + Δω
    ∂(r->
      ħ*ωc*abs2(a) - abs2(Et(0,0,r) + a*Ec(0, 0, r, Δω)),
      z)
end

function map_real2imag(H, n)
    N = size(H,1)
    Hc = ComplexF64.(H)

    T = Matrix{ComplexF64}(I, N, N)

    T[n,   n]   = 1/2
    T[n,   n+1] = im/2
    T[n+1,n]   = 1/2
    T[n+1,n+1] =  -im/2

    Hnew = T * Hc * T'

    keep = [1:n; n+2:N]

    return Hnew[keep, keep]
end

# Hessian of the Hamiltonian
function HHs(x, y, z, Δω, κ)
    αeq = αeq_c(x,y,z, Δω, κ)
    function f(r, p)
        ωc = ω0 + Δω
        a = r[4] + im * r[5]
        as = p[4] + im * p[5]
        Etot2 = real((Et(r[1], r[2], r[3]) + a*Ec(r[1], r[2], r[3], Δω)) *
						(conj(Et(r[1], r[2], r[3])) + as*Ec(r[1], r[2], r[3], Δω))
					)
                #Etot * conj(Etot)
        
        (p[1]^2+p[2]^2+p[3]^2)/2m + ħ*ωc *real(a*as) - α*Etot2
    end
    #@assert ∂zHs(z, Δω, κ) .< 1e10
	cr = 1e9
	ca = 1e16
	dr = 1/cr
	da = 1e16
	
    HH_ = ForwardDiff.hessian(η->f(
			[cr*η[1], cr*η[2], cr*η[3], ca*η[4], ca*η[5]], 
			[dr*η[6], dr*η[7], dr*η[8], da*η[9], da*η[10]]),
                        [x, y, z, real(αeq), imag(αeq), 
						 0, 0, 0, imag(αeq), -imag(αeq)]
                       )
	
	HH_ = map_real2imag(HH_, 9)	
	HH_ = map_real2imag(HH_, 4)

	HH_
end

end # module

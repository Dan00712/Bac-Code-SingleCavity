module Laser

using LinearAlgebra

using ForwardDiff

using ..Constants

export Ec, Et, αeq_c, δc, γ, Hs, ∂zHs, L, ηdot, isstable, isstable2

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


# function to minimize
function L(z, Δω, κ)
    f0 = α/ħ * Ec(0, 0, z, Δω)^2
    δz = Δω - f0
    fz = im*k0 - im*zR/(z^2+zR^2) - z/(z^2+zR^2)
    gz = -2*z/Wc^2

    real((1+f0/(δz+im*κ/2))*(fz+gz*f0/(δz-im*κ/2)))
end

function map_real2imag(H, n)
    N = size(H,1)
    Hc = ComplexF64.(H)

    T = Matrix{ComplexF64}(I, N, N)

    T[n,   n]   = 1/2
    T[n,   n+1] = -im/2
    T[n+1,n]   = 1/2
    T[n+1,n+1] =  im/2

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

        (p[1]^2+p[2]^2+p[3]^2)/2m + ħ*Δω *real(a*as) - α*Etot2
    end
    #@assert ∂zHs(z, Δω, κ) .< 1e10

    HH_ = ForwardDiff.hessian(η->f(
			[0, 0, η[3], η[4], η[5]],
			[0, 0, η[8], η[9], η[10]]),
                        [0, 0, z, real(αeq), imag(αeq),
						 0, 0, 0, imag(αeq), -imag(αeq)]
                       )

	HH_ = map_real2imag(HH_, 9)
	HH_ = map_real2imag(HH_, 4)

	HH_
end

isstable(z, Δω, κ) = begin #isposdef(HHs(0,0,z, Δω, κ))
	S = HHs(0,0, z, Δω, κ)
	N = Int(size(S, 1) // 2)
	J = let
		A = I(N) .|> ComplexF64
		A[N, N] = 1/im /ħ
		[zeros(N,N)  A ; -A zeros(N,N)]
	end
	K = let
		B = zeros(N, N) .|> ComplexF64
		B[N,N] = κ/2
		[
			B 			zeros(N,N);
			zeros(N,N) 	B
		]
	end

	hx = let
		δ = δc(0, 0, z, Δω)
		γ = α/ħ * Ec(0, 0, z, Δω)^2/(δ^2 + κ^2/4)

		fx = (im * k0/(z^2 + zR^2) - 2/Ax^2 /W(z)^2)
		gx = -2/Wc^2

		(1+γ*δ) * (real(fx) + gx *γ*δ) + γ*κ/2 * (imag(fx) + gx*γ*κ/2)
	end
	hy = let
		δ = δc(0, 0, z, Δω)
		γ = α/ħ * Ec(0, 0, z, Δω)^2/(δ^2 + κ^2/4)

		fy = (im * k0/(z^2 + zR^2) - 2/Ax^2 /W(z)^2)
		kc = let
			ωc = Δω + ω0

			ωc/c
		end
		gy = -kc

		(1+γ*δ) * (real(fy) + gy *γ*δ) + γ*κ/2 * (imag(fy) + gy*γ*κ/2)
	end

	A = J*S - K
	V = eigen(A).values
	ReV = real.(V)

	#@show ma, mi
	foo(x) = x < 0
	(foo.(ReV) |> all) && hx >= 0 && hy >= 0
end

function ηdot(q; Δ, κ)
    z = q[1]
    pz = 0
    a = q[2]
    Ec_ = Ec(0, 0, z, Δ)
    Et_ = Et(0, 0, z)

    adot = -im * (Δ - Ec_^2/ħ) * a + im*α/ħ * Ec_ * Et_ - κ/2*a
    asdot = conj(adot)

    zdot = pz/m
    pdot = α/2 * abs2(Et(0, 0, z)) * L(z, Δ, κ)

    # [̇q, ̇p]
    [zdot, adot, pdot, asdot]
end

function isstable2(z; Δ, κ)
    a = αeq_c(0, 0, z, Δ, κ)
    q = [z, a]
    function stable_pertubation(δq)
        dη = ηdot(q+δq; Δ = Δ, κ = κ)
        dq = dη[1:2]
        ddq = dη[3:end]/m

        # theta:=angle(Δq, δq) >= 90°) <=> cos(theta) <= 0 
        #cos(theta) :=
        (real(δq⋅dq)/norm(δq)/norm(dq) <= 0 && real(δq⋅ddq)/norm(δq)/norm(ddq) <= 0)
    end

    map([
        [BigFloat(z/1e3), 0],
        [BigFloat(-z/1e3), 0],
        #[0, real(a)/1e3],
        #[0, imag(a)/1e3]
    ]) do δq
        stable_pertubation(δq)
    end |> all
end

end # module

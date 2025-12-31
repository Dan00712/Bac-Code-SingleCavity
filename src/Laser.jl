module Laser

using LinearAlgebra

using ..Constants

export Ec, Et, αeq_c, δc, γ, Hs, ∂zHs, L, ηdot, isstable

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

function ηdot(q; Δ, κ)
    z = q[1]
    pz = 0
    a = q[2]
    Ec_ = Ec(0, 0, z, Δ)
    Et_ = Et(0, 0, z)

    adot = -im * (Δ - Ec_^2/ħ) * a + im*α/ħ * Ec_ * Et_
    asdot = conj(adot)

    zdot = α/2 * abs2(Et(0, 0, z)) * L(z, Δ, κ)
    pdot = pz/m

    # [̇q, ̇p]
    [zdot, adot, pdot, asdot]
end

function isstable(z; Δ, κ)
    a = αeq_c(0, 0, z, Δ, κ)
    q = [z, a]
    function stable_pertubation(δq)
        Δη = ηdot(q+δq; Δ = Δ, κ = κ)

        # theta:=angle(Δq, δq) >= 90°) <=> cos(theta) <= 0
        #cos(theta) :=
        real(δq⋅Δη[1:2])/norm(δq)/norm(Δη[1:2]) <= 0
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

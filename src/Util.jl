module Util

using Dates

using Roots

using ..Laser

export now_nodots, get_zeroes

function now_nodots()
    now() |> string |> x->split(x, ".")[1] |> x -> replace(x, ":" => "-")
end

function get_zeroes(f, Z)
    guesses = let
        guesses = []
        zprev = Z[1]
        fprev = f(Z[1])

        for z in Z[2:end]
            fz = f(z)

            if fz * fprev < 0
                push!(guesses, (zprev, z))
            end
            zprev = z
            fprev = fz
        end
        guesses
    end

    zmins = []
    for guess in guesses
        push!(zmins, find_zero(f, guess))
    end

    zmins
end

end #module

using DrWatson
@quickactivate "SingleCavity"

include(scriptsdir("shared_code.jl"))

global_logger(ConsoleLogger(Info))

Ω = vcat(
    range(0, 400, length = 100) .* 1e3*2π,  #.|> x-> 10^x,
) |> sort
Z = range(start = -2-6, stop = log10(50)-6, length = 200) |> x -> 10 .^ x #|> x-> BigFloat.(x)
const κ = 2π * 18e4

p = plot(;
    title = "κ=$((κ/2π/1000)) 2πkHz",
    xlabel = "Δ/(2π kHz)",
    ylabel = "z/μm",
    yaxis = :log,
    grid = true,
    minorgrid = true,
    formatter = :plain,
)

@info "checking equilibrium positions ∀ ω ∈ Ω"
zmins = []
ωs = []
for z in ProgressBar(Z)
    ωs_ = get_zeroes(ω->L(z, ω, κ), Ω)
    append!(zmins, [z for _ in ωs_])
    append!(ωs, ωs_)
end


@info "check stability of equilibrium positions"
stabilities = [
    if isstable(zmin, ω, κ)
        "stable"
    else
        "unstable"
    end for (ω, zmin) in zip(ωs, zmins)
]

@info "save generated data"
date = now_nodots()
M = (Δ = ωs, zmin = zmins, stability = stabilities)
mkpath(datadir("sims", "stability"))
@save datadir("sims", "stability", "s-$date.jld2") M

@info "create plot"
for (g, marker, color) in zip(["stable", "unstable"], markers, [:blue, :orange])
    mask = stabilities .== g
    if !any(mask)
        continue
    end
    scatter!(
        p,
        ωs[mask]/2π/1e3,
        zmins[mask] .* 1e6,
        marker = marker,
        color = color,
        label = string(g),
    )
end

@info "save plot"
mkpath(plotsdir("stabilities"))
savefig(p, plotsdir("stabilities", "stab2-$date.png"))
savefig(p, plotsdir("stabilities.pdf"))
if isinteractive()
    gui(p)
end

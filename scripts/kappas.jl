using DrWatson
@quickactivate "SingleCavity"

include(scriptsdir("shared_code.jl"))
global_logger(ConsoleLogger(Info))

Ω = vcat(
    range(0, 800, length = 100) .* 1e3*2π,  #.|> x-> 10^x,
) |> sort
Z = range(start = -2-6, stop = log10(50)-6, length = 100) |> x -> 10 .^ x #|> x-> BigFloat.(x)
Κ = 2π * [1.5e6, 1e6, 18e4]

p = plot(;
    #         title="κ=$((κ/2π/1000)) 2πkHz",
    xlabel = "Δ/(2π kHz)",
    ylabel = "z/μm",
    yaxis = :log,
    grid = true,
    minorgrid = true,
    formatter = :plain,
)

@info "generating data ∀ κ ∈ Κ"
for (κ, color, marker) in zip(Κ, colors, markers)
    zmins = []
    ωs = []
    for z in ProgressBar(Z)
        newmins = get_zeroes(ω->L(z, ω, κ), Ω)
        append!(zmins, [z for _ in newmins])
        append!(ωs, newmins)
    end

    scatter!(
        p,
        ωs ./(2π*1e3),
        zmins .* 1e6,
        label = "κ=$(κ/2π /1e3) 2πkHz",
        marker = marker,
        color = color,
    )
end

@info "save plot"

mkpath(plotsdir("kappas"))
savefig(p, plotsdir("kappas", "kap-$(now_nodots()).png"))
savefig(p, plotsdir("kappas.pdf"))
if isinteractive()
    gui(p)
end

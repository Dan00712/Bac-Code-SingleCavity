using LinearAlgebra
using Logging

using ProgressBars

using Plots
if isinteractive()
    plotlyjs()
else
    gr()
end

using DelimitedFiles
using JLD2

using SingleCavity.Util
using SingleCavity.Constants
using SingleCavity.Laser


const colors = ["orange", "green", "blue"]
const markers = [:cross, :xcross, :dtriangle]

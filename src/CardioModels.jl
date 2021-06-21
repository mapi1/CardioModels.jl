module CardioModels

# deps
using ControlSystems
using Random, Distributions, LinearAlgebra
using DSP
using HypothesisTests
using Parameters
using ControlSystemIdentification
using ControlSystems
using MonteCarloMeasurements
using RecipesBase
# using Plots
# using KardioUtils

include("Baselli.jl")
include("Identification.jl")
include("Karemaker.jl")
include("DeBoer.jl")
include("utils.jl")

export BaselliModel
export predict, predict!
export DeBoerModel
export KaremakerModel, vagal_balance
export getModel
export fitBaselli, postprocess
export phenylephrine
export psdplot, psdplot!


end # module


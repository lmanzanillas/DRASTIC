__precompile__(true)

module DRASTIC

using Images
using ImageFeatures
using FileIO
using Plots
using ImageComponentAnalysis
using ImageFiltering
using Statistics
using StatsBase
using HDF5

export color_dist 
export moving_window_filter
export find_peaks
export get_horizontal_pitch
export get_maximums

include("DRASTIC_functions.jl")
include("DRASTIC_segmentation.jl")

end # module DRASTIC

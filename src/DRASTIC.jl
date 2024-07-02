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
export collect_groups
export divide_img_sq
export get_selection
export get_max_h
export get_average_color
export merge_divided_binary_img

include("DRASTIC_functions.jl")
include("DRASTIC_segmentation.jl")

end # module DRASTIC

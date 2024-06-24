"""
function color_dist(a::Color, b::Color)
function to get the distance between two colors in HSV color space
Each color is converted to HSV color space using cylindrical coordinates and then the cartesian distance is computed
"""
function color_dist(a::Color, b::Color)
    @fastmath begin
        ca = convert(HSV, a)
        cb = convert(HSV, b)
        x1 = ca.s * cos(deg2rad(ca.h))
        y1 = ca.s * sin(deg2rad(ca.h))
        z1 = ca.v
        x2 = cb.s * cos(deg2rad(cb.h))
        y2 = cb.s * sin(deg2rad(cb.h))
        z2 = cb.v
        dist = sqrt( (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2 )
    end
end

"""
	moving_window_filter(w_in::Vector, window::Int)
Function to apply a moving average window filter
It accepts a waveform and the size of the window to perform the average
It returns a new waveform
"""
function moving_window_filter(w_in::Vector{<:Real}, window::Int)
    w_out = zeros(length(w_in))
    half_window = window ÷ 2
    for i = 1 : 1 : half_window
        w_out[i] = sum(w_in[i:i+half_window]) / (half_window + 1)
    end
    for i = half_window + 1 : 1 : length(w_in)-half_window
        w_out[i] = sum(w_in[i-half_window:i+half_window]) / (window + 1)
    end
    for i = length(w_in)- half_window + 1 : 1 : length(w_in)
        w_out[i] = sum(w_in[i-half_window:i]) / (half_window + 1)
    end
    return w_out
end

"""
	function find_peaks(x::AbstractVector{T};
    nups::Int=1,
    ndowns::Int=nups,
    zerostr::Char='0',
    peakpat=nothing,
    minpeakheight=typemin(T),
    minpeakdistance::Int=1,
    threshold=zero(T),
    npeaks::Int=0,
    sortstr=false) where T
Function to find peaks in a given array using severla paramters as input as: minimum height of peaks (minpeakheight), threshold, distance bettwen peaks
the minimal working example for a trapezoidal filter is:
find_peaks(x;minpeakheight=350,minpeakdistance=10,threshold=0)	which will return all peaks with a threshold of 350 units (~keV) and separated of at least 10 samples
"""
function find_peaks(x::AbstractVector{T};
    nups::Int=1,
    ndowns::Int=nups,
    zerostr::Char='0',
    peakpat=nothing,
    minpeakheight=typemin(T),
    minpeakdistance::Int=1,
    threshold=zero(T),
    npeaks::Int=0,
    sortstr=false) where T

    zerostr ∉ ('0', '+', '-') && error("zero must be one of `0`, `-` or `+`")

    # generate the peak pattern with no of ups and downs or use provided one
    peakpat = Regex(peakpat === nothing ? "[+]{$nups,}[-]{$ndowns,}" : peakpat)

    # transform x into a "+-+...-+-" character string
    xs = String(map(diff(x)) do e
        e < 0 && return '-'
        e > 0 && return '+'
        return zerostr
    end)

    # find index positions and maximum values
    peaks = map(findall(peakpat, xs)) do m
        v, i = findmax(@view x[m])
        (;value=v, idx=first(m) + i - 1, start=first(m), stop=last(m) + 1)
    end

    # eliminate peaks that are too low
    filter!(peaks) do p
        p.value >= minpeakheight && p.value - max(x[p.start], x[p.stop]) >= threshold
    end

    # sort according to peak height
    if sortstr || minpeakdistance > 1
        sort!(peaks, by=x -> x.value; rev=true)
    end

    # find peaks sufficiently distant
    if minpeakdistance > 1
        removal = falses(length(peaks))
        for i in 1:length(peaks)
            removal[i] && continue
            for j in 1:length(peaks)
                removal[j] && continue
                dist = abs(peaks[i].idx - peaks[j].idx)
                removal[j] = 0 < dist < minpeakdistance
            end
        end
        deleteat!(peaks, removal)
    end

    # Return only the first 'npeaks' peaks
    npeaks > 0 && resize!(peaks, min(length(peaks), npeaks))

    return peaks
end

"""
function get_horizontal_pitch(Distance::Matrix{<:Real})
function to get the pattern of given structure in a matrix of colors
The function has been optimized to get the column positions of holes. It assumes that the holes are aligned vertically
"""
function get_horizontal_pitch(Distance::Matrix{<:Real})
    Band = Float64[]
    for i = 1 : 1 :length(Distance[1,:])
        mu = mean(Distance[:,i])
        push!(Band,mu)
    end
    Band_filter = moving_window_filter(Band,30)
    ph = mean(Band_filter)/2
    vpitch = find_peaks(Band_filter;minpeakheight=ph,minpeakdistance=50)

    pks = []
    for i = 1:length(vpitch)
        if length(vpitch) == 0
            continue
        else
            amplitude = getproperty(vpitch[i], :value)
            position = getproperty(vpitch[i], :idx)
            push!(pks,[amplitude position])
        end
    end
    peaks = vcat(pks...)
    return Int.(sort(peaks[:,2]))
end

"""
function get_maximums(Distance::Matrix{<:Real},peak::Int)
Function used to get the center of holes, It will takes a matrix of distances using the center of holes as reference color
and an integer (you first need to run get_horizontal_pitch function ) to scan a band corresponding to the row +/- range 
"""
function get_maximums(Distance::Matrix{<:Real},peak::Int,range::Int=7)
    Band = Float64[]
    for i = 1 : 1 :length(Distance[:,1])
        mu = mean(Distance[i,peak-range:peak+range])
        push!(Band,-mu)
    end
    Band_filter = moving_window_filter(Band,40)
    pm = mean(Band_filter) + std(Band_filter)

    vpitch = find_peaks(Band_filter;minpeakheight=pm,minpeakdistance=150)

    pks = []
    for i = 1:length(vpitch)
        if length(vpitch) == 0
            continue
        else
            amplitude = getproperty(vpitch[i], :value)
            position = getproperty(vpitch[i], :idx)
            push!(pks,[amplitude position])
        end
    end
    peaks = vcat(pks...)
    return Int.(sort(peaks[:,2]))
end

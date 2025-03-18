"""
    color_dist(a::Color, b::Color) -> Float64

Computes the Euclidean distance between two colors `a` and `b` 
in the HSV (Hue, Saturation, Value) color space.

# Arguments:
- `a::Color`: The first color.
- `b::Color`: The second color.

# Returns:
- The Euclidean distance between the two colors in HSV space.
"""
function color_dist(a::Color, b::Color)::Float64
    @fastmath begin
        ca = convert(HSV, a)
        cb = convert(HSV, b)
        
        x1, y1, z1 = ca.s * cos(deg2rad(ca.h)), ca.s * sin(deg2rad(ca.h)), ca.v
        x2, y2, z2 = cb.s * cos(deg2rad(cb.h)), cb.s * sin(deg2rad(cb.h)), cb.v

        return sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
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
    maxs = Int[]
    for i = 1 : 1 :length(Distance[:,1])
        mu = mean(Distance[i,peak-range:peak+range])
        push!(Band,-mu)
    end
    Band_filter = moving_window_filter(Band,40)
    pm = mean(Band_filter) + std(Band_filter)
    #check if pm is larger than maximum and adjust if needed
    if pm > maximum(Band_filter)
        pm = mean(Band_filter) + std(Band_filter)/2
    end
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
    if length(peaks) > 0
        maxs = Int.(sort(peaks[:,2]))
    end
    return maxs
end

"""
function correct_slope_x(datos::Matrix{<:Real}, slope_x::Float64)
It assumes your data is according to: x y p
Function to correct the slope when the camera is not completely horizontal
You first need to find the slope and then give it as input to the function
"""
function correct_slope_x(datos::Matrix{<:Real}, slope_x::Real)
    new_signal = Vector{AbstractFloat}(undef, length(datos[:,3]))
    new_signal[1] = datos[1,3]
    for i =2:1:length(datos[:,3])
        new_signal[i] =  datos[i,3] - datos[i,1]* slope_x
    end
    return [datos[:,1] datos[:,2] new_signal]
end

"""
function correct_slope_y(datos::Matrix{<:Real}, slope_y::Float64)
It assumes your data is according to: x y p
Function to correct the slope when the camera is not completely horizontal
You first need to find the slope and then give it as input to the function
"""
function correct_slope_y(datos::Matrix{<:Real}, slope_y::Real)
    new_signal = Vector{AbstractFloat}(undef, length(datos[:,3]))
    new_signal[1] = datos[1,3]
    for i =2:1:length(datos[:,3])
        new_signal[i] =  datos[i,3] - datos[i,2]* slope_y
    end
    return [datos[:,1] datos[:,2] new_signal]
end


"""
function get_shadow_correction(cx::Float64,cy::Float64,x1::Float64,y1::Float64,calibration::Float64,d_camera_anode::Float64)
Apply a shadow correction to determine the transparency of anodes, all values are in mm
"""
function get_shadow_correction(cx::Float64,cy::Float64,x1::Float64,y1::Float64,calibration::Float64,d_camera_anode::Float64)
    r = 1.2 #expected hole diameter
    Area_perfect = pi * r^2
    Anode_tickness = 3.2 #in mm
    camera_to_anode = d_camera_anode # in mm to be measured in your setup
    d_to_center = sqrt( (x1-cx)^2 + (y1-cy)^2 )/calibration ## in mm
    Δ = (Anode_tickness/camera_to_anode)*d_to_center 
    Area = 2*r^2*acos(Δ/2*r) - Δ*sqrt(r^2 - (Δ/2)^2) #Expected correction from trigonometric considerations
    return sqrt(Area/Area_perfect)
end

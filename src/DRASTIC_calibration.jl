"""
function get_calibration_factor(calib_img::AbstractArray{RGB{N0f8}},red_calib_color::HSV{Float32},threshold::Float64=0.15)
function to obtain the calibration factor to convert pixels to mm
It accepts an image containing the red circles for calibration. It assumes the red circles are 10 mm diameter, but can be adjuste
It also need the red color of that circles, a threshold which if not given is taken as 0.15 which is a good aproximation
It will return the mean value of the red circles used for the calibration, in general two
"""
function get_calibration_factor(calib_img::AbstractArray{RGB{N0f8}},red_calib_color::HSV{Float32},threshold::Float64=0.15,d_red_circles::Float64= 10.0)
    distance = color_dist.(calib_img,red_calib_color)
    binary_img =  distance .< threshold
    components = Images.label_components(binary_img)
    measurements = analyze_components(components, BasicMeasurement(area = true, perimeter = true))
    min_area_red = π*(100)^2 #Need to estimate min radius for security check
    max_area_red = π*(100*15.0/2)^2 #Need to estimate maximum radius for security check
    filter_result = measurements[ min_area_red .< measurements[!,:area] .< max_area_red, :]
    diam = []
    if length(filter_result[!,:l]) > 2
        @warn "more than 2 circles found, check calibration manually"
    end
    for comp = 1 : 1 : length(filter_result[!,:l])
        r = sqrt.(filter_result[comp,:area]/pi)#pixels
        d = 2*r/d_red_circles;
        push!(diam,d)
    end
    return mean(diam)
end

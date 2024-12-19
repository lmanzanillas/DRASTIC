"""
   get_calibration_factor(calib_img::AbstractArray{RGB{N0f8}},red_calib_color::HSV{Float32},threshold::Float64=0.15)
function to obtain the calibration factor to convert pixels to mm
It accepts an image containing the red circles for calibration. It assumes the red circles are 10 mm diameter, but can be adjusted.
It also needs the red color of those circles, a threshold which if not given is taken as 0.15 which is a good approximation.
It will return the mean value of the red circles used for the calibration, in general two.
"""
function get_calibration_factor(calib_img::AbstractArray{RGB{N0f8}}, red_calib_color::HSV{<:AbstractFloat}, threshold::T=0.15, d_red_circles::T=10.0) where {T<: AbstractFloat}
    distance = color_dist.(calib_img, red_calib_color)
    binary_img = distance .< threshold
    components = Images.label_components(binary_img)
    measurements = analyze_components(components, BasicMeasurement(area=true, perimeter=true))
    
    min_area_red = π * 100^2  # Precompute π * 100^2
    max_area_red = π * (100 * 15.0 / 2)^2  # Precompute π * (100 * 15.0 / 2)^2
    
    filter_result = measurements[(min_area_red .< measurements[!, :area]) .& (measurements[!, :area] .< max_area_red), :]
    
    if length(filter_result[!, :l]) > 2
        @warn "more than 2 circles found, check calibration manually"
    end

    diam = (filter_result[!, :area] ./ π) .|> sqrt .|> (r -> 2r / d_red_circles)
    #The previous line is equivalent to:
    #areas = filter_result[!, :area] ./ π  # Step 1: Divide area by π
    #radii = sqrt.(areas)                 # Step 2: Take element-wise square root
    #diam = radii .|> (r -> 2r / d_red_circles)  # Step 3: Apply the function element-wise
    
    return mean(diam)
end

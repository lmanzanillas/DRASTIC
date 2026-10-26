"""
    indiv_photo_analysis(photo_path, pixel_per_mm=29.93, pitch_dir_angle=150., method="Otsu";
                         tile=448, hole_diameter_mm=2.4, area_range_mm=(1.5, 2.9),
                         pitch_window=0.6, outlier_band=0.4)
 
Detect holes in `photo_path` and compute the pitch along `pitch_dir_angle`.
 
Returns `(p, μ)`:
  * `p`: matrix with columns x [px], y [px], pitch [mm], outliers removed
  * `μ`: mean pitch [mm] of the remaining points (`NaN` if nothing was found)
"""
function indiv_photo_analysis(photo_path::AbstractString,
                              pixel_per_mm = 29.93,
                              pitch_dir_angle = 150.,
                              method = "Otsu";
                              tile::Int = 448,
                              hole_diameter_mm = 2.4,
                              area_range_mm = (1.5, 2.9),
                              pitch_window = 0.6,
                              outlier_band = 0.4)
    @info "Starting individual photo analysis" photo_path
    calib = pixel_per_mm
    min_area_hole = π*(calib*area_range_mm[1]/2)^2
    max_area_hole = π*(calib*area_range_mm[2]/2)^2
 
    anode_image = load(photo_path)
 
    # overlap of about one hole diameter
    overlap = round(Int, hole_diameter_mm*calib)
    full_bin = binarize_tiled(anode_image; tile = tile, overlap = overlap, method = method)
 
    # flip vertically so that y points up in the plot
    full_bin = full_bin[end:-1:1, :]
 
    results = get_holes_basic_info(full_bin, min_area_hole, max_area_hole)
    data = rotate(results, pitch_dir_angle)
    p = get_pitch(data, calib, pitch_window)
    p = rotate(p, -pitch_dir_angle)   # rotate back by the same angle
 
    if isempty(p)
        @warn "No pitch values found" photo_path
        return p, NaN
    end
 
    # robust outlier filter around the median, then recompute the mean
    m = median(p[:, 3])
    p = p[abs.(p[:, 3] .- m) .< outlier_band, :]
    μ = isempty(p) ? NaN : mean(p[:, 3])
    return p, μ
end



"""
    binarize_tiled(img; tile=448, overlap=75, method="Otsu", sigma=2,
                   window_size=109, percentage=15.0) -> BitMatrix
 
Binarize `img` tile by tile (local thresholds, to cope with uneven lighting),
but compute each tile's blur and threshold on a padded window that extends
`overlap` pixels beyond the tile on every side (clamped at the image border).
Only the central `tile x tile` core is written to the output, so:
 
  * the Gaussian blur never sees an artificial border inside the image,
  * neighbouring tiles get nearly the same threshold, because their padded
    windows share a wide band of pixels,
  * no separate divide/merge step or index bookkeeping is needed.
 
Returns `true` where a hole (dark region) is detected, in the same orientation
as `img`.
 
Choose `overlap` of about one hole diameter (2.4 mm * pixel_per_mm, about 72 px
here), and at least 3*sigma.
"""
function binarize_tiled(img::AbstractMatrix; tile::Int = 448, overlap::Int = 75,
                        method::AbstractString = "Otsu", sigma = 2,
                        window_size::Int = 109, percentage::Real = 15.0)
    method in ("Otsu", "AdaptiveThreshold") || error("Unknown method: $method")
    H, W = size(img)
    full_bin = falses(H, W)
    kern = Kernel.gaussian(sigma)
 
    for r0 in 1:tile:H, c0 in 1:tile:W
        # core region (what we keep)
        r1 = min(r0 + tile - 1, H)
        c1 = min(c0 + tile - 1, W)
        # padded region (what we process)
        pr0 = max(r0 - overlap, 1);  pr1 = min(r1 + overlap, H)
        pc0 = max(c0 - overlap, 1);  pc1 = min(c1 + overlap, W)
 
        section = img[pr0:pr1, pc0:pc1]
        blurred = imfilter(section, kern)
 
        bw = if method == "Otsu"
            .!Bool.(binarize(blurred, Otsu()))
        else
            .!Bool.(binarize(blurred, AdaptiveThreshold(percentage = percentage,
                                                        window_size = window_size)))
        end
 
        # crop the padded result back to the core
        rr = (r0 - pr0 + 1):(r1 - pr0 + 1)
        cc = (c0 - pc0 + 1):(c1 - pc0 + 1)
        full_bin[r0:r1, c0:c1] = bw[rr, cc]
    end
    return full_bin
end



"""
function collect_groups(labels)
function to collect the labels and improve in efficiency
"""
function collect_groups(labels)
    groups = [Int[] for i = 1:maximum(labels)]
    for (i,l) in enumerate(labels)
        if l != 0
            push!(groups[l], i)
        end
    end
    groups
end

"""
function divide_img_sq(full_img::AbstractArray{RGB{N0f8}},nPixels::Int = 500)
function to divide an image in small sections for processing
The function receives an image of a given size and an integer indicated the size in pixels of the small sections
"""
function divide_img_sq(full_img::AbstractArray{RGB{N0f8}},nPixels::Int = 500)
    vertical_pixels = size(full_img)[1]
    horizontal_pixels = size(full_img)[2]
    vertical_chunks = vertical_pixels ÷ nPixels
    horizontal_chunks = horizontal_pixels ÷ nPixels
    v_starts = zeros(Int,vertical_chunks)
    v_ends = zeros(Int,vertical_chunks)
    Sections = fill(Matrix{RGB{N0f8}},vertical_chunks,horizontal_chunks)
    ##Find start and end of sections in vertical axis
    for i = 1 : 1 : vertical_chunks
        v_s = (i- 1)*nPixels + 1
        v_e = nPixels * i
        if i == vertical_chunks
            v_e = vertical_pixels
        end
        v_starts[i] = v_s
        v_ends[i] = v_e
        #println("vertical start: ",v_starts," ends: ",v_ends)
    end
    ##Find start and end of sections in horizontal axis
    h_starts = zeros(Int,horizontal_chunks)
    h_ends = zeros(Int,horizontal_chunks)
    for i = 1 : 1 : horizontal_chunks
        h_s = (i- 1)*nPixels + 1
        h_e = nPixels * i
        if i == horizontal_chunks
            h_e = horizontal_pixels
        end
        h_starts[i] = h_s
        h_ends[i] = h_e
    end
    ##Do the division
    s_img = []
    index_v_h = zeros(Int,vertical_chunks,horizontal_chunks)
    for i = 1 : 1 : vertical_chunks
        for j = 1 : 1 : horizontal_chunks
            segment_img = full_img[v_starts[i]:v_ends[i],h_starts[j]:h_ends[j]]
            push!(s_img,segment_img)
            index_v_h[i,j] = (i - 1)*horizontal_chunks + j
        end
    end
    return s_img, index_v_h
end

"""
function get_selection(dist::Matrix{<:Real},max_h::Float64 = 0.2)
function to get the cut used for the selection of the holes
It get the dist matrix and the value of the peak distance to copper, which is obtained
using an histogram
"""
function get_selection(dist::Matrix{<:Real},max_h::T = 0.2, nSigmas::T = 3.0) where {T<: AbstractFloat}
    v_dist = vcat(dist...)
    cut_l = max_h - 0.05
    cut_h = max_h + 0.05
    f_v_dist = filter(x-> cut_l < x < cut_h,v_dist)
    std_dev_dist = std(f_v_dist)
    return max_h - nSigmas*std_dev_dist
end


"""
function get_max_h(h::Histogram, min_value::Float64)
function to determine the coordinate of the bin containing the maximum of counts in a histogram
It has been optimized assuming the hole color is used as reference
The distance will maximize when copper regions are found and minimize when hole regions are found
In general min_value of 0.3 is a good starting point to avoit the peaks close to zero
"""
function get_max_h(h::Histogram, min_value::T)where {T<: AbstractFloat}
    observed_counts = h.weights
    bin_edges = h.edges[1]
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2
    start_h = Int(min_value ÷ bin_widths[1]) 
    _max = findmax(observed_counts[start_h:end])[2] + start_h - 1
    return bin_centers[_max]
end

"""
function get_average_color(SectionImg::AbstractArray{RGB{N0f8}},Dist::Matrix{<:Real},hPeaks::Vector{Int},range::Int=5)
Function to get the average color of the holes in a given image
It take the image, the matrix distance, and the horizontal pattern hPeaks, and a range that will be used for computing the mean of each hole
It returns the average hole color of the image
"""
function get_average_color(SectionImg::AbstractArray{RGB{N0f8}},Dist::Matrix{<:Real},hPeaks::Vector{Int},range::Int=5)
    average_color_vector = []
    for h = 1 : 1 : length(hPeaks)
        #if too close to border then do not consider the point
        if abs(hPeaks[h] - size(Dist)[2]) < 2*range || abs(hPeaks[h] - range) < 2*range
            continue
        end
        local_max = get_maximums(Dist,hPeaks[h])
        for lm in local_max
            if abs(lm - size(Dist)[1]) > 2*range && abs(lm -range) > 2*range
                push!(average_color_vector,mean(SectionImg[lm - range:lm + range,hPeaks[h] - range:hPeaks[h] + range]))
            end
        end
    end
    return HSV{Float32}(mean(average_color_vector))
end

"""
function get_average_color_cu(SectionImg::AbstractArray{RGB{N0f8}},Dist::Matrix{<:Real},hPeaks::Vector{Int},range::Int=4)
same as get_average_color but use instead copper color
"""
function get_average_color_cu(SectionImg::AbstractArray{RGB{N0f8}},Dist::Matrix{<:Real},hPeaks::Vector{Int},range::Int=4)
    average_color_vector = []
    for h = 1 : 1 : length(hPeaks)
        #if too close to border then do not consider the point
        if abs(hPeaks[h] - size(Dist)[2]) < 2*range || abs(hPeaks[h] - range) < 2*range
            continue
        end
        local_max = get_maximums(-Dist,hPeaks[h])
        for lm in local_max
            if abs(lm - size(Dist)[1]) > 2*range && abs(lm -range) > 2*range
                push!(average_color_vector,mean(SectionImg[lm - range:lm + range,hPeaks[h] - range:hPeaks[h] + range]))
            end
        end
    end
    return HSV{Float32}(mean(average_color_vector))
end

"""
function merge_divided_binary_img(v_sections::Vector,index_sections::Matrix{Int})
function to merge all sections of a segmented image
it accepts a vector containing the sections and a matrix containing the indices of the sections that 
are required to the merge
"""
function merge_divided_binary_img(v_sections::Vector,index_sections::Matrix{Int})
    merg_h = []
    for i = 1 : 1 : length(index_sections[:,1])
        push!(merg_h,hcat(v_sections[index_sections[i,:]]...))
    end
    full_img = vcat(merg_h...)
    return full_img
end

"""
function get_holes_info(comp::BitMatrix, min_area = 5000, max_area = 15000)
function to get the hole measurements, it take a binary image and obtain the required info
It will return the coordinates of the center of each hole ant its diameter
"""
function get_holes_info(comp::BitMatrix, min_area = 5000, max_area = 15000)
    components = Images.label_components(comp)
    measurements = analyze_components(components, BasicMeasurement(area = true, perimeter = true))
    centroids = component_centroids(components)
    filter_result = measurements[ min_area .< measurements[!,:area] .< max_area , :]
    centroids = component_centroids(components)
    filter_centroides = Tuple{Float64, Float64}[]
    diam = Float64[]
    for comp = 1 : 1 : length(filter_result[!,:l])
        r = sqrt.(filter_result[comp,:area]/pi)#pixels
        d = 2*r;
        push!(diam,d)
        push!(filter_centroides,centroids[filter_result[comp,:l]])
    end
    Centers = hcat(first.(filter_centroides), last.(filter_centroides))
    Diameters = reshape(diam, length(diam), 1)
    [Centers[:,2] Centers[:,1] Diameters]
end

"""
function get_holes_basic_info(comp::BitMatrix, min_area = 5000, max_area = 15000)
Similar to get_holes_info but much faster since it uses a much simpler algorith
results should be comparable
"""
function get_holes_basic_info(comp::BitMatrix, min_area = 5000, max_area = 15000)
    components = Images.label_components(comp)
    centroids = component_centroids(components)
    filter_centroides = Tuple{Float64, Float64}[]
    groups = collect_groups(components)
    diameters = Float64[]
    for g = 1 : 1 : length(groups)
        area_group = length(groups[g])
        if min_area < area_group && area_group < max_area
            push!(diameters,2*sqrt.(area_group/pi))
            push!(filter_centroides,centroids[g])
        end
    end
    Centers = hcat(first.(filter_centroides), last.(filter_centroides))
    Diameters = reshape(diameters, length(diameters), 1)
    return [Centers[:,2] Centers[:,1] Diameters]
end

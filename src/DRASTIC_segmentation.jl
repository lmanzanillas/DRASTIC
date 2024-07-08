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
function get_selection(dist::Matrix{<:Real},max_h::Float64 = 0.2,nSigmas::Int=3)
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
function get_max_h(h::Histogram, min_value::Float64)
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

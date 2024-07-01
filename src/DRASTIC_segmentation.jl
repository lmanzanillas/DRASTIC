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
function get_selection(h::Histogram, min_value::Float64)
function to determine the cut selection according to the distance to a given color
It has been optimized assuming the hole color is used as reference
The distance will maximize when copper regions are found and minimize when hole regions are found
In general min_value of 0.3 is a good starting point
"""
function get_selection(h::Histogram, min_value::Float64)
    observed_counts = h.weights
    bin_edges = h.edges[1]
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2
    min_check = 0.2
    start_check = Int(min_check ÷ bin_widths[1])
    max_check = findmax(observed_counts[start_check:end])[2] + start_check
    start_h = Int(min_value ÷ bin_widths[1])
    #security check
    if abs(max_check - start_h) < 10
        start_h = start_check
    end
    println(max_check," ",start_h)
    _max = findmax(observed_counts[start_h:end])
    vmax = _max[1]
    bmax = _max[2]
    _min = findmin(observed_counts[start_h:start_h+bmax])
    vmin = _min[1]
    bmin = _min[2]
    cut = 0.
    for i = start_h + bmin : 1 : start_h + bmax
        if observed_counts[i] > 1.3 * vmin
            cut = i * bin_widths[1]
            break
        end
    end
    return cut
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
        if abs(hPeaks[h] - size(Dh)[2]) < 2*range || abs(hPeaks[h] - range) < 2*range
            continue
        end
        local_max = get_maximums(Dist,hPeaks[h])
        for lm in local_max
            if abs(lm - size(Dh)[1]) > 2*range || abs(lm -range) > 2*range
                push!(average_color_vector,mean(SectionImg[lm - range:lm + range,hPeaks[h] - range:hPeaks[h] + range]))
            end
        end
    end
    return HSV{Float32}(mean(average_color_vector))
end

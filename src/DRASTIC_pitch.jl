"""
function rotate(Data::Matrix{<:Real},Angle::Real)
function to rotate in a given drirection data in an X-Y plane
It accepts a matrix and it will consider that the first column corresponds to the x coordinate
while the second column corresponds to the Y c orrdinate
"""
function rotate(Data::Matrix{<:Real},Angle::Real)
    θ = deg2rad(Angle)
    Xnew = Data[:,1] * cos(θ) .- Data[:,2] * sin(θ)
    Ynew = Data[:,1] * sin(θ) .+ Data[:,2] * cos(θ)
    return [Xnew Ynew Data[:,3]]
end

"""
function get_index_rows(Data::Matrix{<:Real},Calib=33.9)
function to get the index of rows in a given 2d distribution of holes
it assumes that the holes have been sorted in the y axis
if the distance in y between 2 holes a new row will be identified
it will return an array with the indexes of the first and last hole of each row
"""
function get_index_rows(Data::Matrix{<:Real},Calib=33.9)
    Index_transition = []
    row_ends = -1
    row_starts = 1
    for h = 2 : 1 : length(Data[:,2])
        if h == length(Data[:,2])
            row_ends = h
        end
        if Data[h,2] - Data[h-1,2] >  Calib* 1.0/2
            row_ends = h - 1
            push!(Index_transition,[row_starts row_ends])
            row_starts = h
        end
    end
    return vcat(Index_transition...)
end

"""
function get_pitch(Data::Matrix{<:Real},calib=45.0,tolerance = 0.2)
function to get the pitch ina 2d distribution of holes
only holes with a diameter +/- tolerance will be used, this avoid using holes with wrong reconstruction
It will return an array with the coordiantes and pitch of each hole, except the first hole since the pitch is measured between two holes, then
in each row of n holes we will have n-1 measurements
"""
function get_pitch(Data::Matrix{<:Real},calib=45.0,tolerance = 0.2)
    #sort by y coordinate
    sorted_data = Data[sortperm(Data[:, 2]), :]
    #identify rows
    Index_rows = get_index_rows(sorted_data,calib)
    PITCH = []
    hole_diameter = 2.4 #expected value 
    for i = 1 : 1 : length(Index_rows[:,1])
        ROW = sorted_data[Index_rows[i,1] : Index_rows[i,2],:]
        #sort row by x axis
        sorted_ROW = ROW[sortperm(ROW[:, 1]), :]
        for h = 2 : 1 : length(sorted_ROW[:,1])
            if abs(sorted_ROW[h,3]/calib - hole_diameter) > tolerance
                continue
            end
            x1 = sorted_ROW[h-1,1]
            x2 = sorted_ROW[h,1]
            y1 = sorted_ROW[h-1,2]
            y2 = sorted_ROW[h,2]
            pitch = sqrt( (x2-x1)^2 + (y2 -y1)^2)/calib
            push!(PITCH,[x2 y2 pitch])
        end
    end
    return vcat(PITCH...)
end

"""
function get_section_merged(base_dir::String,photo_list::Vector{Int64}, dx::Real,x_size::Real,y_size::Real,my_cte::Real=31.0,rot_angle::Real=0., expected_pitch::Real = 2.94, tolerance::Real=0.3)
Function to merge photos of aa given photo shooting campaing
The list of photos is needed to do the merge. This list is provided by the person doing the photo shooting
"""
function get_section_merged(base_dir::String,photo_list::Vector{Int64}, dx::Real,x_size::Real,y_size::Real,my_cte::Real=31.0,rot_angle::Real=0., expected_pitch::Real = 2.94, tolerance::Real=0.3)
    h5_files = base_dir*"/".*readdir(base_dir);
    counter = 0
    Column = []
    d_cut_low = expected_pitch - tolerance
    d_cut_high = expected_pitch + tolerance
    x_cut_low = 150.
    x_cut_high = x_size - 150.
    for photo in photo_list
        counter += 1
        my_file = filter(x->occursin("0"*string(photo)*".h5",x),  h5_files)
        if length(my_file) == 0
            continue
        end
        dy = y_size*counter
	datos = h5read(my_file[1], "Diameter")
	#rotate the data to find the pitch in a given direction
        data = DRASTIC.rotate(datos,rot_angle)
        pitch = get_pitch(data,my_cte,tolerance)
        pitch = pitch[ d_cut_low .< pitch[:,3] .<  d_cut_high, :]
        #rotate the data to return to the original position
        pitch = DRASTIC.rotate(pitch,-rot_angle)
        pitch = pitch[ x_cut_low .< pitch[:,1] .<  x_cut_high, :]

        new_y = dy .+ pitch[:,2]
        new_x = pitch[:,1] .+ dx
        new_z = pitch[:,3]

        col = [new_x new_y new_z]
        push!(Column,col)
    end
    return vcat(Column...)
end

"""
function get_section_merged_diameter(base_dir::String,photo_list::Vector{Int64}, dx::Real,x_size::Real,y_size::Real,my_cte::Real=31.0,rot_angle::Real=0., expected_diam::Real = 2.94, tolerance::Real=0.3, cut_x_h::Real=150., cut_x_l::Real=150., cut_y_h::Real=150., cut_y_l::Real=150.)
function to merge photos using diameters of holes
"""
function get_section_merged_diameter(base_dir::String,photo_list::Vector{Int64}, dx::Real,x_size::Real,y_size::Real,my_cte::Real=31.0,rot_angle::Real=0., expected_diam::Real = 2.94, tolerance::Real=0.3, cut_x_h::Real=150., cut_x_l::Real=150., cut_y_h::Real=150., cut_y_l::Real=150.)
    h5_files = base_dir*"/".*readdir(base_dir);
    counter = 0
    Column = []
    d_cut_low = expected_diam - tolerance
    d_cut_high = expected_diam + tolerance
    x_cut_low = cut_x_l
    x_cut_high = x_size - cut_x_h
    y_cut_low = cut_y_h
    y_cut_high = y_size - cut_y_l
    for photo in photo_list
        counter += 1
        my_file = filter(x->occursin(string(photo)*".h5",x),  h5_files)
        if length(my_file) == 0
            continue
        end
        dy = y_size*counter
        datos = h5read(my_file[1], "Diameter")

        datos = datos[ d_cut_low .< datos[:,3]/my_cte .<  d_cut_high, :]
        datos = datos[ x_cut_low .< datos[:,1] .<  x_cut_high, :]
        datos = datos[ y_cut_low .< datos[:,2] .<  y_cut_high, :]

        new_y = dy .+ datos[:,2]
        new_x = datos[:,1] .+ dx
        new_z = datos[:,3]/my_cte

        col = [new_x new_y new_z]
        push!(Column,col)
    end
    return vcat(Column...)
end

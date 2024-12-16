using DelimitedFiles  # For reading text files
using Interpolations  # For interpolation
using LinearAlgebra   # For matrix operations if needed
using Dierckx
using MRIReco, FileIO, MRIFiles, MRIBase, MRICoilSensitivities
using Plots

include("plot_gradients.jl")

function interpolate_nom_kspace(k_nom)

    # Parameters for measured k-space
    nS_k_meas = 16384  # Number of samples in measured k-space
    dt_k_meas = 2.5e-6  # Dwell time of measured k-space trajectories

    # Time axes
    t_nom = (1:size(G_nom_data, 1)) .* 1e-5 .- 5e-6  # Time axis for G_nom
    t_meas = (1:nS_k_meas) .* dt_k_meas .- dt_k_meas / 2  # Time axis for measured signals
    t_out = t_meas[1:end-1] .+ dt_k_meas / 2  # Time axis for calculated input signals

end


function interpolate_nom_gradient(filename::String)
    # Load nominal gradient file
    G_nom_data = readdlm(filename)

    # Restrict size of G_nom
    size_res = 4096
    idx_crop = 3775 # 140mm1p71mmR1 , Matlab create_all_trajectoryfiles.m
    # Find index of last data point with smaller time step (assuming time step of 2.5 us)
    t_end = idx_crop *10/2.5 

    # Parameters for measured k-space
    nS_k_meas = 16384  # Number of samples in measured k-space
    dt_k_meas = 2.5e-6  # Dwell time of measured k-space trajectories

    # Time axes
    t_nom = (1:size(G_nom_data, 1)) .* 1e-5 .- 5e-6  # Time axis for G_nom
    t_meas = (1:nS_k_meas) .* dt_k_meas .- dt_k_meas / 2  # Time axis for measured signals
    t_out = t_meas[1:end-1] .+ dt_k_meas / 2  # Time axis for calculated input signals

    # Extend G_nom with zeros to avoid errors during extrapolation
    G_nom_length = size(G_nom_data, 1)
    ext = size_res - G_nom_length  # Calculate extension size
    G_nom = vcat(G_nom_data, zeros(ext, 3))  # Append zeros to extend G_nom
    t_nom = vcat(t_nom, t_nom[end] .+ (1:ext) .* 10e-6)  # Extend t_nom accordingly

    # Interpolate G_nom to t_out
    G_nom_interp = hcat(
        [LinearInterpolation(t_nom, G_nom[:, col], extrapolation_bc = Line())(t_out) for col in 1:3]...
    )

   # Crop G_nom_interp to the length tEnd
   G_nom_interp_cropped = G_nom_interp[1:Int(t_end)-1, :]

   # Rearrange the columns as per your requirement
   # - First column -> G_nom_interp_cropped[:, 2] (x-direction)
   # - Second column -> G_nom_interp_cropped[:, 1] (y-direction)
   # - Third column stays the same -> G_nom_interp_cropped[:, 3] (z-direction)
   G_nom = hcat(
       G_nom_interp_cropped[:, 1],  # First column: G_nom_interp_cropped[:, 2] (x-direction)
       G_nom_interp_cropped[:, 2],  # Second column: G_nom_interp_cropped[:, 1] (y-direction)
       G_nom_interp_cropped[:, 3]   # Third column stays the same: G_nom_interp_cropped[:, 3] (z-direction)
   )
   
    # Return interpolated gradient data
    return G_nom
end



## ---- IDEA: modified function for readding nominal gradient file ---

"""
# Arguments
* `filename` - filename (with full path) of text file with gradient waveform information
* `fov:: Tuple{Int64,Int64,Int64}` - size of reconstructed image (trailing dimension 1 for 2D acquisitions)
        from params[:fov] = [192,192,1] [mm]
Keyword arguments:
* `delay` - delay in seconds from the nominal first sampling point to the actual first sampling point
*  'plot' -Bool ; plotting graddients and k-space trajectory or not

function read_gradient_text_file(filename, reconsize, delay)
"""

function read_nominal_gradient_file(filename,fov; delay = 0.00, doplot = false)

    ## Parameter given as variables to function
    fov_m = fov.* 1e-3 #[m]
    # Keep the same as fov for now
    reconsize  = fov
    voxelsize = fov_m ./reconsize

    G_nom_data = readdlm(filename)
    #G_nom_data = interpolate_nom_gradient(filename)


    ####
    n_samples_tot = length(G_nom_data)
    n_samples = length(G_nom_data[:,1,1])
    dwell_time = 1.0e-5 # s
    acq_duration = n_samples * dwell_time # s

    # Parameters, from matlab file


    gradient_dict = Dict{Symbol,Any}()
    gradient_dict[:version_number] = "#4"
    gradient_dict[:num_samples] = n_samples_tot
    gradient_dict[:dwell_time] = dwell_time # [seconds]
    gradient_dict[:samples_per_interleave] = n_samples
    gradient_dict[:num_interleaves] = 1
    gradient_dict[:num_dims] = 2 # G_x and G_y
    gradient_dict[:time_to_center_kspace] = 0.0 # [seconds]
    gradient_dict[:acq_duration] = acq_duration #[seconds]
    # USIKKER PÅ DISSE
    gradient_dict[:samples_per_acq] = 0 
    gradient_dict[:num_acq] = 1
    gradient_dict[:acq_TR] = 0.0
    gradient_dict[:gradient_acq_start_delay] = 0.0
    ####
    gradient_dict[:echo_time_shift_samples] = 0.0
    gradient_dict[:fov] = fov_m # [m]
    gradient_dict[:voxel_dims] = voxelsize # [m] USIKKER
    gradient_dict[:gradient_strength_factor] = 1 # [mT/m] ENDRE
    gradient_dict[:is_binary] = 0
    gradient_dict[:gamma] = 42577.478 # [Hz/mT] CAN CHANGE
    gradient_dict[:field_strength] = 7 # [T] CAN CHANGE

    #print(gradient_dict)
    # Extract G_x and G_y gradients
    G_x = G_nom_data[:, 1]
    G_y = G_nom_data[:, 2]

    # Initialize the gradient array
    gradient_array = Array{Float64, 3}(undef, gradient_dict[:samples_per_interleave], gradient_dict[:num_interleaves], gradient_dict[:num_dims]) # [mT/m]
    # Populate the gradient array
    gradient_array[:, 1, 1] = G_x  # Populate G_x in the first dimension
    gradient_array[:, 1, 2] = G_y  # Populate G_y in the second dimension


    planned_times = gradient_dict[:dwell_time] .* (0:(gradient_dict[:samples_per_interleave]-1))
    delayed_times = planned_times .- delay .- gradient_dict[:dwell_time] ./ 2 # seconds (dwell_time/2 compensates for integration)


    gradient_array_new = Array{Float64,3}(undef, size(gradient_array))

    ## Loop over all of the unique excitation trajectories and create an interpolant of the gradient
    for dim = 1:gradient_dict[:num_dims]

        for l = 1:gradient_dict[:num_interleaves]

            #print((dim,l),"\n")

            sp = Spline1D(planned_times, gradient_array[:, l, dim], w = ones(length(planned_times)), k = 1, bc = "zero", s = 0.0)

            # evaluate the interpolant at the sampling times of the kspace data
            gradient_array_new[:, l, dim] = sp(delayed_times)

            #print(interleave_gradient_array_new[:,l,dim][end],"\n")

        end

    end


    ## cumulative summation and numerical integration of the gradient data, resulting in the kspace trajectory
    kspace_trajectory_array_new = gradient_dict[:gamma] * gradient_dict[:dwell_time] * cumsum(gradient_array_new, dims = 1) # [rad/m]

    ## Conversion to the trajectory scaling convention in MRIReco.jl
    #  Currently only 2d Trajectories
    converted_kspace_trajectory_array_new = kspace_trajectory_array_new
    converted_kspace_trajectory_array_new[:, :, 1] *= gradient_dict[:fov][1] ./ reconsize[1]
    converted_kspace_trajectory_array_new[:, :, 2] *= gradient_dict[:fov][2] ./ reconsize[2]



    ## Construction of the trajectory object ##

    ## Reshaping of the array to the format expected by the Trajectory constructor in MRIReco.jl
    # - dim 1 = kspace dimension
    # - dim 2 = kspace position (with interleaves/profiles arranged consecutively)
    permuted_trajectory =
        permutedims(reshape(converted_kspace_trajectory_array_new, 
        gradient_dict[:samples_per_interleave] * gradient_dict[:num_interleaves], gradient_dict[:num_dims]), [2, 1])

    ## Construction of the trajectory
    # - Note: timing vectors are automatically generated - seems to be consistent with the dwell time
    trajectory_object = Trajectory(
        permuted_trajectory,
        gradient_dict[:num_interleaves],
        gradient_dict[:samples_per_interleave],
        TE = gradient_dict[:echo_time_shift_samples],
        AQ = gradient_dict[:acq_duration],
        numSlices = 1,
        cartesian = false,
        circular = false,
    )

    if doplot
        @info "Plotting gradients"
        # Extract the gradients for the first interleave
        G_x_nom = gradient_array[:, 1, 1]  # Gradient in x-direction
        G_y_nom = gradient_array[:, 1, 2]  # Gradient in y-direction

        G_x_new = gradient_array_new[:,1,1]
        G_y_new = gradient_array_new[:,1,2]
        # Define the time axis
        time_points = delayed_times  # This assumes delayed_times is precomputed

        # Example Usage
        # Full plot
        plot_gradients(time_points, G_x_nom, G_y_nom, G_x_new, G_y_new, range = 1:50)
        @info "Plotting k-space trajectory"
        ### Plot trajectory

        # Extract k-space trajectory for the first interleave and dimensions
        k_x = kspace_trajectory_array_new[:, 1, 1]  # k-space trajectory in x-direction
        k_y = kspace_trajectory_array_new[:, 1, 2]  # k-space trajectory in y-direction


        # Plot the 2D spiral trajectory
        p2 = Plots.plot(
            k_x, k_y,
            label = "Spiral Trajectory",
            xlabel = "k_x (rad/m)",
            ylabel = "k_y (rad/m)",
            linewidth = 2,
            title = "2D Spiral k-Space Trajectory ",
            aspect_ratio = 1,  # Ensures equal scaling of axes
            legend = :topright

        )
        display(p2)
    end
    return trajectory_object
end



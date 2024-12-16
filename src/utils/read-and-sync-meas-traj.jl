# Test utils script for modifying gradient_reader


"""
    read_gradient_text_file(filename, recon_size, delay)
Reads in text file containing gradient waveform information

# Arguments
* `filename` - filename (with full path) of trajectory file : trajectory calculated by method of Marens thesis
* `recon_size::Tuple{Int64,Int64,Int64}` - size of reconstructed image (trailing dimension 1 for 2D acquisitions)

"""
function read_gradient_text_file_mod(filename, recon_size)

    # gradient_data acquired by meethod described in section 2.4 Measuring k-space trajectories in Marens master
    gradient_data = readdlm(filename,',','\n')

    ## Read in the header data of the gradient text file (lines 1 to 21)
    ## 
    gradient_dict = Dict{Symbol,Any}()
    gradient_dict[:B0_phase] = gradient_data[:,1]; # in [rad], this is the B0- eddy-current phase
    gradient_dict[:k_traj] = gradient_data[:,2:end]./(2*pi); # in [1/m]

    gradient_dict[:gamma] = 42577.478 # [Hz/mT] CAN CHANGE
    gradient_dict[:field_strength] = 7 # [T] CAN CHANGE


    ## Conversion to the trajectory scaling convention in MRIReco.jl (dimensionless form used by MRIReco (k = [-0.5, 0.5]))
    #  Currently only 2d Trajectories
    gradient_dict[:k_traj][:,1] = gradient_dict[:k_traj][:,1].* recon_size[1] / recon_size[1] / 1000; 
    gradient_dict[:k_traj][:, 2] = gradient_dict[:k_traj][:, 2] .* recon_size[2] / recon_size[2] / 1000;

    ## Construction of the trajectory object ##

    ## Reshaping of the array to the format expected by the Trajectory constructor in MRIReco.jl
    # - dim 1 = kspace dimension
    # - dim 2 = kspace position (with interleaves/profiles arranged consecutively)
    permuted_trajectory = permutedims(gradient_dict[:k_traj] , [2, 1])

    ## Construction of the trajectory
    # - Note: timing vectors are automatically generated - seems to be consistent with the dwell time
    trajectory_object = Trajectory(
        permuted_trajectory,
        1,
        size(gradient_dict[:k_traj],2),
        TE = 0.0 ,
        AQ = size(gradient_dict[:k_traj],2)*1e-6 ,
        numSlices = 1,
        cartesian = false,
        circular = true,
    )

    return trajectory_object, gradient_dict

end

function sync_traj_and_data_mod!(raw_data, traj, idx_crop, interleave, B0_phase)

    N_signal = size(raw_data.profiles[1].data, 1) #size of signal vector

    N_traj = size(traj.nodes, 2) #size of trajectory vector = number of sampling points in k-space

    
    # define dwell times for trajectory (dt_k) and signal sampling (dt_s)
    dt_s = 2.5*10^(-6)  # [s]
    dt_k = 2.5*10^(-6)  # [s]

    t_s = (1:N_signal)*dt_s .- dt_s/2

    t_k = (1:N_traj)*dt_k .- dt_k/2 #traj type = "meas"

    ### Interpolate trajectory onto the same sample times as the sampled signal
    trajNodes_interpolated_X = Spline1D(t_k, traj.nodes[1,:], w=ones(length(traj.nodes[1,:])), k=3)
    trajNodes_interpolated_Y = Spline1D(t_k, traj.nodes[2,:], w=ones(length(traj.nodes[1,:])), k=3)

    B0_phase_spline = Spline1D(t_k, B0_phase, w=ones(length(B0_phase)), k=3)
    B0_ph_eddy = B0_phase_spline(t_s)
    
    

    ### Concatenate the trajectory node kx and ky positions
    adjustedTraj = vcat(trajNodes_interpolated_X(t_s)', trajNodes_interpolated_Y(t_s)')
    # synchronize trajectory data and the kspace data

    

    ## Add trajectory to RawAcquisitionData.profiles[i]

    if idx_crop != 0
        for i = 1:length(raw_data.profiles)
        raw_data.profiles[i].traj = adjustedTraj[:, 1:idx_crop]
        raw_data.profiles[i].data = raw_data.profiles[i].data[1:idx_crop, :]
        end
        t_s = t_s[1:idx_crop]
        B0_ph_eddy = B0_ph_eddy[1:idx_crop]
    else
        for i = 1:length(raw_data.profiles)
        raw_data.profiles[i].traj = adjustedTraj
        end
    end
    

    # Return the vector of sampling times
    length_B0 = size(B0_ph_eddy)
    @info "Length B0_ph_eddy : $length_B0"
    return dt_s * (0:idx_crop-1), B0_ph_eddy


end

function adjust_header_mod!(raw_data, num_samples)
    # Adjust header (see adjustHeader! in Utils.jl:223)
    for i = 1:length(raw_data.profiles)

        # Set the discard post to 0 (don't discard any samples from the end of the acquisition)
        raw_data.profiles[i].head.discard_post = 0
    
        # Set the discard pre to 0 (don't discard any samples from the beginning of the acqusition)
        raw_data.profiles[i].head.discard_pre = 0
    
        ## Set the number of samples properly 
        raw_data.profiles[i].head.number_of_samples = num_samples 
    
        ## Set center sample to 0 (only for spiral scans)
        raw_data.profiles[i].head.center_sample = 0
    
    end

end

"""
    merge_raw_interleaves_mod(params, output_raw)
Does what the original merge_raw_interleaves 

# Arguments
* `params`          - Dictionary
* `resolution`      - Constant
* 'applyk0'         - Bool
"""


function merge_raw_interleaves_mod(params, resolution, mm_shift, applyk0corr::Bool = false)


    # @info "indices = $interleave_complement" #DEBUG

    # read in the data file from the ISMRMRD format
    data_file = ISMRMRDFile(params[:interleave_data_filenames][params[:interleave]])

    # Get the trajectory object
    traj, gradient_dictionary  = read_gradient_text_file_mod(params[:traj_filename], params[:recon_size])

    # Read in raw data from the data_file
    raw_data = RawAcquisitionData(data_file)

    # @info "indices = $ic" #DEBUG
    # The B0-eddycurrent-phase
    B0_phase_eddy = gradient_dictionary[:B0_phase]

    # set up time vector for tracking all of the interleaves
    
    time_track_vector = []

    # synchronize trajectory data and the kspace data onto the same time 
    times, B0_eddy_ph_correct =  sync_traj_and_data_mod!(raw_data, traj, params[:num_samples], 1, B0_phase_eddy)

    adjust_header_mod!(raw_data, params[:num_samples])

    # add the times to the time tracking vector
    append!(time_track_vector, times)
    
    # Convert RawAcquisitionData to AcquisitionData
    @info "Converting RawAcquisitionData to AcquisitionData"
    acq_data = AcquisitionData(raw_data, estimateProfileCenter=true)

    ## Update acq_data.traj, from Utils-Maren.jl
    for i = 1:length(acq_data.traj)
    acq_data.traj[i].times = time_track_vector # set times to the total time vector
    acq_data.traj[i].TE = 0.00 # set the TE to 0
    acq_data.traj[i].AQ = time_track_vector[end] # set the acquisition time to the last element of the time vector (should be the latest time)
    acq_data.traj[i].circular = false # set whether to use a circular filter on the kspace data
    end

    # Shift image by mm_shift
 
    k = acq_data.traj[1].nodes[2, :]
    for s = 1:size(acq_data.kdata, 3)
        for c = 1:size(acq_data.kdata[s], 2)
            acq_data.kdata[s][:, c] = acq_data.kdata[s][:, c] .* exp.((im*2*pi*mm_shift/resolution) .* k)
        end
    end

    if applyk0corr
        @info "Correct for B0-phase"
        for numRep = 1:20
            acq_data.kdata[numRep] = acq_data.kdata[numRep] .* exp.(-1im .* B0_eddy_ph_correct)
        end
        
        #figure()
        #plot(B0_phase)
        #title("B0-phase")
        #ylabel("[rad]")
    end

    #Check acquisition nodes
    acq_data.traj[1].nodes[acq_data.traj[1].nodes[:] .> 0.5] .= 0.5
    acq_data.traj[1].nodes[acq_data.traj[1].nodes[:] .< -0.5] .= -0.5

    return acq_data
end


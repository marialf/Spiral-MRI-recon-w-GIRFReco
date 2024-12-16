#-----------------------------------------------------------------------------------
# # Script for running both the cartesian reconstruction with the field-map scan and calculate SENSE maps and use this in the spiral reconstruction 
# # 

#-----------------------------------------------------------------------------------

#=
# Uses:
- the configuration file 'recon_config_mod.jl' for setting general parameters
- cartesian scan (dual echo):
    - In vivo: inVivo_fieldMap220mm.h5
    - Phantom: 32ch_phantom_fm_converted.h5
- spiral scan files: 
    - In vivo scan: inVivo_140mm1p71mmR1.h5
    - Phantom scan: 32ch_140mm1p71mmR1.h5
- trajectory file (name for the spiral measured gradient file) : k_meas_140mm1p71mmR1.txt
=#

#=
## 1. Setup

The necessary Julia packages needed for spiral reconstruction.
=#


#Our developed packages
using GIRFReco, MRIGradients

#MRIReco and its sub-packages
using MRIReco, FileIO, MRIFiles, MRIBase, MRICoilSensitivities

using RegularizedLeastSquares

using ImageTransformations

using Dierckx

using MosaicViews

using DelimitedFiles, FourierTools, ROMEO, Unitful, ImageView

#import GIRFReco: load_map, save_map
#=
## 2. Configurations for reconstruction

The following file, [`recon_config_joss_demo.jl`](@__REPO_ROOT_URL__/docs/lit/examples/recon_config_joss_demo.jl),
includes general configuration for spiral reconstruction.
It is necessary to execute this file to make sure all parameters are loaded.
Sample Data that works with this script can be found [here](https://doi.org/10.5281/zenodo.7779044).
Please download, extract and set the `root_project_path` as the top level folder (should be something like `/your/path/joss_data_zenodo/`)
=#

# Set the root path to the data (also where the results will be stored)
data_root_project_path = "/Volumes/MasterB/MariaThesis/DATA_test_recon/IN-VIVO" 
# Include configuration file
include("in_vivo_config.jl")
# Include utils file. 
include("/Users/mariafoyen/Documents/GitHub/Specialization-project/simple_recon_30-09/GIRFReco_dependencies/src/utils/my_utils.jl")
# Include file for 
include("/Users/mariafoyen/Documents/GitHub/Specialization-project/simple_recon_30-09/GIRFReco_dependencies/src/utils/read-and-sync-meas-traj.jl")
### Configurations

num_total_diffusion_directions = params_general[:num_total_diffusion_directions]
## Determine to reconstruct single-interleave data, or one interleave out of multi-interleave data.
is_single_interleave = ~(length(params_general[:scan_fullpath]) > 1)

start_idx_interleave = 1;
slice_selection = [1]

#############################################

### Cartesian reconstruction :
##  calculation of sensitivity maps and B0 map from dual echo cartesian can

#############################################

# Only run when coil/B0 maps have not been calculated
if ~(params_general[:do_load_maps] && isfile(params_general[:b0_map_save_fullpath]))
    @info "Running cartesian_recon to retrieve maps (cartesian_sensitivity and b0_maps)"
    run_cartesian_recon(params_general) 
end

# Load the SENSE maps from the previously calculated NIfTI files.
@info "Loading SENSE and B0 maps from $(params_general[:sensitivity_save_fullpath])"
cartesian_sensitivity = load_map(params_general[:sensitivity_save_fullpath]; do_split_phase = true)
b0_maps = load_map(params_general[:b0_map_save_fullpath])
num_slices = size(b0_maps, 3)

#############################################

### Parameters for the spiral data 

#############################################

reload_spiral_data = true; # Set true if we need to reload raw data compulsively.


#=
The first step is to select the part of spiral k-space data that we 
would like to reconstruct. This include selecting slices, diffusion directions, 
and averages that we want.
=#

params_spiral = Dict{Symbol,Any}() 
params_spiral[:recon_size] = Tuple(params_general[:recon_size])
params_spiral[:interleave] = start_idx_interleave
params_spiral[:num_samples] = params_general[:num_adc_samples]
params_spiral[:delay] = 0.00000
# Full paths of raw k-space data files of spiral acquisition
params_spiral[:interleave_data_filenames] = params_general[:scan_fullpath] 
# Full paths of k-space trajectory txt file
params_spiral[:traj_filename] = params_general[:gradient_fullpath]
params_spiral[:do_multi_interleave] = !is_single_interleave
params_spiral[:do_odd_interleave] = false
params_spiral[:num_interleaves] = is_single_interleave ? 1 : length(params_spiral[:interleave_data_filenames]) # one interleaf per file, count files, if filenames are array of strings (not only one string)
params_spiral[:single_slice] = true

resolution = 1.71
num_slices = 1
mm_shift = 7

# Calculate slice index with spatially ascending order in the raw kspace data file.
raw_temp  = RawAcquisitionData(ISMRMRDFile(params_general[:scan_fullpath][1]))
slice_idx_array_spiral = get_slice_order(raw_temp, num_slices, (num_slices+1)*2, 2)

### Load the spiral data and adjust the sensititvy array (the adjustment is not really necessary considering we only are dealing with one slice here)

if reload_spiral_data || !(@isdefined imaging_acq_data)
    @info "Reading spiral data and merging interleaves"
    imaging_acq_data = merge_raw_interleaves_mod(params_spiral, resolution , mm_shift, true)
    b0_maps = b0_maps[:, :, invperm(slice_idx_array_spiral)]
    cartesian_sensitivity = cartesian_sensitivity[:, :, invperm(slice_idx_array_spiral), :]
    # Load the b0 maps only to be able to call the plotting function
end

#shift_kspace!(imaging_acq_data, params_general[:fov_shift])

#=
#### 3.2.5 Processing Coil Sensitivity Maps

We need to preprocess the coil sensitivity maps before reconstruction. 
This includes resizing the coil maps to the size of output encoding matrix size; 
compress the channels according to user's setting to achieve a faster reconstruction.
=#
sensitivity = mapslices(x -> imresize(x, params_spiral[:recon_size][1], params_spiral[:recon_size][2]), cartesian_sensitivity, dims = [1, 2])
resized_b0_maps = mapslices(x -> imresize(x, params_spiral[:recon_size][1], params_spiral[:recon_size][2]), b0_maps, dims = [1, 2])

# Optional: Plot the sensitivity maps of each coil on a given slice.
if params_general[:do_plot_recon]
    plotlyjs(size=(1000, 800))
    plot_sense_maps(sensitivity, 32)
end

#############################################

### Parameters for the spiral reconstruction

#############################################
selected_slice = 1

@info "Setting parameters for simple reconstruction"
params_recon_simple = Dict{Symbol,Any}()
params_recon_simple[:reco] = "direct"
params_recon_simple[:reconSize] = params_spiral[:recon_size][1:2] # cannot avoid camel-case here since it is defined by MRIReco.jl and RegularizedLeastSquares.jl
params_recon_simple[:regularization] = "L2"
params_recon_simple[:λ] = 1e-3
params_recon_simple[:iterations] = params_general[:num_recon_iterations]
params_recon_simple[:solver] = "cgnr"
params_recon_simple[:solverInfo] = SolverInfo(ComplexF32, store_solutions = false)


@info "Performing spiral recon without SENSE"
@time reco_simple = reconstruction(imaging_acq_data, params_recon_simple)

@info "Setting Reconstruction Parameters for the reconstruction with SENSE"
params_recon_sense = deepcopy(params_recon_simple)
params_recon_sense[:reco] = "multiCoil"
params_recon_sense[:senseMaps] = ComplexF32.(sensitivity[:, :, :, :])

@info "Performing Spiral Reconstruction with SENSE"
@time reco_sense = reconstruction(imaging_acq_data, params_recon_sense)

############################################## PLOTTING FOR SENSE AND SIMPLE ##################

@info "Plotting sense reconstruction"
plotlyjs(size=(1000, 800))
plot_reconstruction(
    reco_sense[:, :, 1, 1, 1, 1],
    1:length(selected_slice),
    resized_b0_maps,
    is_slice_interleaved = false,
    rotation = 0,
)
 

# Plot the difference between simple and sense reconstruction
plot_difference_recon(reco_simple[:, :, 1, 1, 1, 1],reco_sense[:, :, 1, 1, 1, 1], "Simple ", "SENSE", 1 ,270,  is_slice_interleaved = false, plot_phase_difference = false)

############################################## RECON WITH B0 ##################

@info "Setting reconstruction parameters for reconstruction with SENSE and B0 correction"
params_recon_B0 = deepcopy(params_recon_sense)
params_recon_B0[:correctionMap] = ComplexF32.(-1im .* resized_b0_maps[:, :, selected_slice])


@info "Performing Spiral Reconstruction with SENSE and B0 map correction"
@time reco_sense_b0_meas = reconstruction(imaging_acq_data, params_recon_B0)


@info "Plotting sense reconstruction"
plotlyjs(size=(800, 800))
plot_reconstruction(
    reco_sense_b0[:, :, 1, 1, 1, 1],
    1:length(selected_slice),
    resized_b0_maps,
    is_slice_interleaved = false,
    rotation = 0,
)

# Plot the difference between sense and sense+b0 reconstruction
plot_difference_recon(reco_sense[:, :, 1, 1, 1, 1],reco_sense_b0[:, :, 1, 1, 1, 1], "Simple ", "SENSE", 1 ,270,  is_slice_interleaved = false, plot_phase_difference = false)

############################################## RECON WITH NOMINAL TRAJECTORY ##################
params_spiral_nom = deepcopy(params_spiral)
params_spiral_nom[:traj_filename] = params_general[:gradient_nom_fullpath]

if reload_spiral_data || !(@isdefined imaging_acq_data)
    @info "Reading spiral data and merging interleaves"
    imaging_acq_data_nom = merge_raw_interleaves_mod(params_spiral_nom, resolution , mm_shift, true)
end

@info "Performing Spiral Reconstruction with SENSE"
@time reco_sense_nom = reconstruction(imaging_acq_data_nom, params_recon_sense)

@info "Performing Spiral Reconstruction with SENSE and B0 map correction for nominal trajectory"
@time reco_sense_b0 = reconstruction(imaging_acq_data_nom, params_recon_B0)

# Plot the difference between sense and sense+b0 reconstruction
plot_difference_recon(reco_sense_b0[:, :, 1, 1, 1, 1],reco_sense_b0_meas[:, :, 1, 1, 1, 1], "Simple ", "SENSE", 1 ,270,  is_slice_interleaved = false, plot_phase_difference = false)

@info "Plotting sense reconstruction"
plotlyjs(size=(800, 800))
plot_reconstruction(
    reco_sense_b0[:, :, 1, 1, 1, 1],
    1:length(selected_slice),
    resized_b0_maps,
    is_slice_interleaved = false,
    rotation = 0,
)
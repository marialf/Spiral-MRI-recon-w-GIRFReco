## This recon_config.jl file describes all reconstruction parameters, as well as data locations and selections for an iterative non-Cartesian reconstruction that relies 
#  on an external reference scan (Cartesian) to estimate calibration maps (coil sensitivities, B0 maps)

using Dates

"""
For this modified version: simplify, without girf correction, coil compression etc....
Setting majority of parameters to 0 (default."""

params_general = Dict{Symbol,Any}()
# Gyromagnetic ratio, in unit of Hz
params_general[:gamma] = 42577478;
## General options for recon script
params_general[:do_load_maps] = true             # if true, reloads B0/SENSE maps instead of recalculating
params_general[:do_save_recon] = true            # if true, saves reconstruction and all auxiliary image data (maps) as NIfTI files
params_general[:do_plot_recon] = true             # if true, plots intermediate debugging and output recon figures (needs graphics, not recommended in multi-thread mode due to PyPlot)
params_general[:do_process_map_scan] = true         # if true, compute sensitivity and B0 maps from reconstructed Cartesian scan   
params_general[:do_save_processed_map_scan] = false; # save ISMRMD file of preprocessed Cartesian data (before recon)

## Reconstruction Parameters
# update time stamp for new recon, otherwise keep fixed, will create a new recon/<recon_id> directory
#params_general[:recon_id] = Dates.format(Dates.now(), "yyyy-mm-dd_HH_MM_SS") # recon ID is recon_id
# params_general[:recon_id] = "2022-10-20_09_07_07"
params_general[:recon_id] = "2024-11-27_15_00_00";
params_general[:do_correct_with_b0_map] = false
params_general[:do_correct_with_girf_k1] = false
params_general[:do_correct_with_girf_k0] = false

params_general[:num_virtual_coils] = 0;
params_general[:do_coil_compression] = false;
params_general[:fov_shift] = [0, 0]; # Unit: number of voxels (mmshift /resolution = 10 mm/2 mm = 5)

## Scan parameters, Additional acquisition information, e.g., slice distance etc.
params_general[:slice_distance_factor_percent] = 000 # 400

#Total number of ADC points BEFORE the rewinder at the end of the spiral readout. For gradient 508, use 15655 (out of 16084); for gradient 511, use 15445 (out of 15624).
params_general[:num_adc_samples] = 15100 
# Matrix size of the reconstructed image. For gradient 508 with all 4 interleaves, use 200 for high resolution image; otherwise consider using 112 or 84 for a lower resolution. The FOV is 220 mm for both gradients 508 and 511.
params_general[:recon_size] = [220, 220, 1] # 192 for all phantom images
params_general[:num_recon_iterations] = 10; # number of recon iterations (for both Cartesian and Spiral recon)
params_general[:b0_map_beta] = 0.1 # for estimate_b0_maps, * `β` - Regularization parameter controlling roughness penalty (larger = smoother, default 5e-4)
params_general[:do_normalize_recon] = false # set max abs to 1
params_general[:saving_scalefactor] = 1.0e9 # 1 # typical range of recon intensities is 1e-7, rescale when saving, e.g., to 0...1000 roughly for fMRI analysis
params_general[:num_total_diffusion_directions] = 0;          # Need to specify total diffusion directions included in the raw data
params_general[:slice_distance_factor_percent] = 0 # Scan parameters, Additional acquisition information, e.g., slice distance etc.
# TEST: istedenfor å lagre sense_mapet, hvis det blir noe feil i load_map()
params_general[:sense_maps] = 0
# Data selector
#  Choose diffusion direction; starting from 0 (b=0) to the total number in MDDW protocol, e.g. for 6 diffusion directions, 1-6 stands for 6 DWIs)
# boolean is_called_from_global_recon is true, if this RunReconLoop is active
# If is_called_from_global_recon is false or not defined, the data selector needs to be defined here.
if !(@isdefined is_called_from_global_recon) || !is_called_from_global_recon
    global selector = Dict{Symbol,Any}()
    selector[:avg] = 1
    selector[:seg] = 1
    selector[:dif] = 0
end

#=
### Specifying Directories
=#

params_general[:project_path] = data_root_project_path # Root path for the project

#Path to ISMRMRD files (raw k-space data) [Input]
params_general[:data_path] = params_general[:project_path]
#Path to spiral readout gradient files [Input]
params_general[:gradients_path] = joinpath(params_general[:project_path], "Trajectories")
#Path to middle results (coil and B₀ maps) files [Output]
params_general[:results_path] = joinpath(params_general[:project_path], "results", "in-vivo")
#Path to final reconstructed spiral images [Output]
params_general[:recon_save_path] = joinpath(params_general[:results_path], "recon", params_general[:recon_id]);


#=
### Specifying File Names
=#
#=
- cartesian scan (dual echo):
    - IN VIVO: inVivo_fieldMap220mm.h5
    - PHANTOM: 32ch_phantom_fm_converted.h5
- spiral scan files: 
    - IN VIVO: inVivo_140mm1p71mmR1.h5
    - PHANTOM: 32ch_140mm1p71mmR1.h5
- trajectory file (name for the spiral measured gradient file) : k_meas_140mm1p71mmR1.txt
=#

## Input files: 

# Map scan file (Cartesian multi-echo file)
params_general[:map_scan_filename] = "Cartesian-fieldmap-scan/inVivo_fieldMap_220mm.h5" # Cartesian dual-echo file, for coil and B₀ maps calculation [Input]
params_general[:map_scan_filename_stem] = "inVivo_fieldMap_220mm.h5"
# The different TEs from 220523_Maren.pdf / Contrast - Common
params_general[:mapTEs_ms] = [4.08, 5.1] 

 # File name for the spiral gradient [Input]
params_general[:gradient_filename] = joinpath("k_meas_140mm1p71mmR1.txt");
 # File name for the NOMINAL spiral gradient [Input]
 params_general[:gradient_nom_filename] = joinpath("k_nom_140mm1p71mmR1.txt");
#  non-Cartesian (Spiral) scan file: MDDW30
params_general[:scan_filename] = ["Spirals/inVivo_140mm1p71mmR1.h5"] # ISMRMRD Raw k-space data for spiral acquisition [Input]

#Output files

params_general[:scan_filename_stem] = "in-vivo_140mm1p71mmR1.h5" # Main file name when saving the result
params_general[:processed_map_scan_filename] = "in-vivo_processed_cartesian_data.h5" # file name for preprocessed data (remove oversampling, permute dimensions wrt MRIReco) [Output]
params_general[:map_save_filename] = splitext(params_general[:map_scan_filename_stem])[1] * "_reconmap.nii" # File name for reconstructed dual-echo Cartesian images [Output]
params_general[:sensitivity_save_filename] = splitext(params_general[:map_scan_filename_stem])[1] * "_sensemap.nii" # File name for calculated coil sensitivity maps [Output]
params_general[:b0_map_save_filename] = splitext(params_general[:map_scan_filename_stem])[1] * "_b0map.nii"; # File name for calculated off-resonance (B₀) maps [Output]

#=
File name for the final reconstructed spiral image.
If we reconstructing multiple spiral data files (e.g. multiple interleaves) through `RunReconLoop.jl`, 
the file name for the final reconstructed image is concatenated from multiple scan file names. 
Otherwise, just append `_recon.nii` as suffix to file name.
=#

if isa(params_general[:scan_filename], AbstractVector)
    # for multiple files, concatenate recon name from scan file names, e.g., 508_124_2_508_126_2_508_128_2_508_130_2_recon.nii
    params_general[:recon_save_filename] =
        join([(x[1] * "_") for x in splitext.(params_general[:scan_filename])]) *
        "dif$(selector[:dif])_" *
        "itl$(selector[:seg])_" *
        "avg$(selector[:avg])_" *
        "recon.nii"
else
    # otherwise, just concat _recon.nii to file name
    params_general[:recon_save_filename] = splitext(params_general[:scan_filename])[1] * "_recon.nii"
end

#=
### Assembling Full Paths

Assembling directories and file names for final full pathes. 
These are automated operations.
=#
params_general[:gradient_fullpath] = joinpath(params_general[:gradients_path], params_general[:gradient_filename]) # Full paths of spiral readout gradients
params_general[:gradient_nom_fullpath] = joinpath(params_general[:gradients_path], params_general[:gradient_nom_filename]) # Full paths of spiral readout gradients
#params_general[:girf_fullpath] = joinpath.(params_general[:girf_path], params_general[:girf_filename]) # Full paths of GIRF files
params_general[:map_scan_fullpath] = joinpath(params_general[:data_path], params_general[:map_scan_filename]) # Full path of dual-echo Cartesian data
params_general[:scan_fullpath] = joinpath.(params_general[:data_path], params_general[:scan_filename]) # Full paths of raw k-space data files of spiral acquisition
params_general[:processed_map_scan_fullpath] = joinpath(params_general[:recon_save_path], params_general[:processed_map_scan_filename]) # Full paths of pre-processed Cartesian dual-echo data [Output]
params_general[:recon_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:recon_save_filename]) # Full paths of the reconstructed spiral image [Output]
params_general[:map_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:map_save_filename]) # Full paths of reconstructed dual-echo Cartesian images [Output]
params_general[:sensitivity_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:sensitivity_save_filename]) # Full paths of calculated coil sensitivity maps [Output]
params_general[:b0_map_save_fullpath] = joinpath(params_general[:recon_save_path], params_general[:b0_map_save_filename]); # Full paths of calculated off-resonance (B₀) maps [Output]

#=
## Final Steps

Optional: If the path for results writing is not existing, create it.

As the last step of configuration, copy this config file 
to the recon path for further checking and debugging purposes.
=#


if ~ispath(params_general[:recon_save_path])
    mkpath(params_general[:recon_save_path])
end

# copies this config file to the recon path for later checks of parameter functions
cp(@__FILE__, joinpath(params_general[:recon_save_path], "recon_config.jl"); force = true)


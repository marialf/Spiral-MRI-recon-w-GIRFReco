using PyPlot
using MosaicViews
using ImageUtils

function load_map(filename; do_split_phase::Bool = false)

    # if separate mag and phase are saved, load and combine them
    if do_split_phase

        magnitude_filename = splitext(filename)[1] * "_magn.nii"
        magnitude_image = loadImage(magnitude_filename) # map is needed, because abs.(im) would convert AxisArray back into basic array

        phase_filename = splitext(filename)[1] * "_phase.nii"
        phase_image = loadImage(phase_filename)

        calib_map = (magnitude_image.data) .* exp.(1im .* (phase_image.data))

    else

        I = loadImage(filename)
        calib_map = I.data

    end

    # squeeze singleton dimensions of 6-dim array
    ## DO THIS DIFFEREENTLY THAN ORIGINAL; considering we have single slice:
    # Check if the 5th and 6th dimensions are singleton (i.e., size 1)
    dims_to_drop = []
    if size(calib_map, 5) == 1
        push!(dims_to_drop, 5)
    end
    if size(calib_map, 6) == 1
        push!(dims_to_drop, 6)
    end

    # Drop the 5th and 6th dimensions only if they are singleton
    if !isempty(dims_to_drop)
        calib_map = dropdims(calib_map, dims = tuple(dims_to_drop...))
    end
    # calib_map = dropdims(calib_map, dims = tuple(findall(size(calib_map) .== 1)...))

    return calib_map

end

data_root_project_path = "/Volumes/MasterB/MariaThesis/DATA_test_recon/PHANTOM"
include("phantom_config.jl")

# Load the SENSE maps from the previously calculated NIfTI files.
@info "Loading SENSE and B0 maps from $(params_general[:sensitivity_save_fullpath])"
cartesian_sensitivity = load_map(params_general[:sensitivity_save_fullpath]; do_split_phase = true)
b0_maps = load_map(params_general[:b0_map_save_fullpath])
# Load the cartesian recon also, for plotting
cart_recon = load_map("/Volumes/MasterB/MariaThesis/DATA_test_recon/PHANTOM/results/phantom/recon/2024-10-17_12_00_00/32ch_phantom_fm_converted_reconmap.nii"; do_split_phase = true)
num_slices = size(b0_maps, 3)

# Resize B0-maps and signal mask to match encoding size of data matrix
plotB0 = mapslices(x->imresize(x, (96, 96)), b0_maps, dims=[1,2])
resized_b0_maps = mapslices(x -> imresize(x, Tuple(params_general[:recon_size])[1], Tuple(params_general[:recon_size])[2]), b0_maps, dims = [1, 2])


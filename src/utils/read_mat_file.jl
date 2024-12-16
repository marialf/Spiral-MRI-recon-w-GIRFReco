using MAT  # For reading MATLAB .mat files

# Specify the path to the .mat file
matfile_path = "/Volumes/MasterB/MarenThesis/master_thesis/results/SpiralMeasurements/spiral_results_32ch_34mm_ny.mat"

# Load the .mat file
mat_contents = matread(matfile_path)

# Access the fields of the struct in the .mat file
gradient_nominal = mat_contents["gradient"]["nominal"]
gradient_nominal_interp = mat_contents["gradient"]["nominal_interp"]
gradient_nominal_shifted = mat_contents["gradient"]["nominal_shifted"]
gradient_measured = mat_contents["gradient"]["measured"]
gradient_measured_filtered = mat_contents["gradient"]["measured_filtered"]
gradient_B0 = mat_contents["gradient"]["B0"]
gradient_B0_filtered = mat_contents["gradient"]["B0_filtered"]
gradient_predicted = mat_contents["gradient"]["predicted"]
gradient_t = mat_contents["gradient"]["t"]
gradient_tNom = mat_contents["gradient"]["tNom"]

k_nominal = mat_contents["k"]["nominal"]
k_nominal_interp = mat_contents["k"]["nominal_interp"]
k_nominal_shifted = mat_contents["k"]["nominal_shifted"]
k_measured = mat_contents["k"]["measured"]
k_measured_filtered = mat_contents["k"]["measured_filtered"]
k_B0phase = mat_contents["k"]["B0phase"]
k_B0phase_filtered = mat_contents["k"]["B0phase_filtered"]
k_predicted = mat_contents["k"]["predicted"]
k_t = mat_contents["k"]["t"]
k_tNom = mat_contents["k"]["tNom"]

tshiftx = mat_contents["tshiftx"]
tshifty = mat_contents["tshifty"]
specs = mat_contents["specs"]
nom_length = mat_contents["nom_length"]

println(gradient_nominal_interp == G_nom_interp./1e-3)

println(gradient_nominal_interp[:3,1,8])
# File to test the read_nominal_gradient.jl functions and compare with files from matlab output

using MAT

include("read_nominal_gradient.jl")

fov = Tuple([192,192,1])

# Gradient filen
fname = "/Volumes/MasterB/MariaThesis/DATA_test_recon/PHANTOM/Trajectories/ArbGradientRO_140mm1p71mmR1.txt"

# Egen funksjon, sjekke outputen

G_nom_interp = interpolate_nom_gradient(fname)

G_nom = read_nominal_gradient_file(fname,fov; delay = 0.00, doplot = true)

# Specify the path to the .mat file
matfile_path = "/Volumes/MasterB/MarenThesis/master_thesis/results/SpiralMeasurements/spiral_results_32ch_34mm_ny.mat"

# Load the .mat file og sammen likne medd G-nom_interp
mat_contents = matread(matfile_path)


gradient_nominal = mat_contents["gradient"]["nominal"]

# Gradienten soom korresponderer til vår spiral sekvens (140mm1p71mmR1)
current_nom_grad = gradient_nominal[:,:,8]

gradient_nominal_interp = mat_contents["gradient"]["nominal_interp"]
current_nominal_grad = gradient_nominal_interp[:, :, 8]
t_end = 3775 *10/2.5
current_nominal_grad_interp = current_nominal_grad[1:Int(t_end)-1, :].*1e3
G_nom_interp

# SJEKKE OM DE ER LIKE
println(isapprox(current_nominal_grad_interp,G_nom_interp)) #TRUE 26.11.24


### Teste funksjon
trajectory_test = read_nominal_gradient_file(fname,fov; delay = 0.00, doplot = false)

num_sampl = trajectory_test.numSamplingPerProfile

interleave_extr_vec =  num_sampl * (1 - 1) .+ (1:num_sampl)

traj = trajectory_test.nodes[:,interleave_extr_vec]

println(traj[:,1:10].*10e-3)
#### Reprod noe av Marens matlab kode

G_nom_data = readdlm(G_nom_fname)
n_samples_tot = size(G_nom_data,1)
tIn = (1:size(G_nom_data, 1)) * 1e-5 .- 5e-6

size_res = 4096
ext = size_res - size(G_nom_data, 1)

# Extend G_nom with zeros
G_nom_data = vcat(G_nom_data, zeros(ext, size(G_nom_data, 2)))

# Extend tIn with the additional time points
tIn = vcat(tIn, tIn[end] .+ (1:ext) .* 10e-6)



######
fov = Tuple([192,192,1])
nom_traj = read_nominal_gradient_file(G_nom_fname,fov; delay = 0.00, doplot = true)


#### Unddersøke k_nom og k_meas Matlab filene
k_nom_name = "/Volumes/MasterB/MarenThesis/master_thesis/code/Julia/data/Trajectories/k_nom_140mm1p71mmR1.txt"
k_meas_name = "/Volumes/MasterB/MarenThesis/master_thesis/code/Julia/data/Trajectories/k_meas_140mm1p71mmR1.txt"

k_nom_mat_data =  readdlm(k_nom_name,',','\n')
k_traj_nom_mat = k_nom_mat_data[:,2:end][:,1]

#####
include("/Users/mariafoyen/Documents/GitHub/Specialization-project/test_joss_recon_owndata_08_30-09/sync-traj-utils-test.jl")

traj_nom, dict_nom = read_gradient_text_file_mod(k_nom_name, fov)
traj_meas, dict_meas = read_gradient_text_file_mod(k_meas_name, fov)

k_x_mat_nom = traj_nom.nodes[1,:]
k_y_mat_nom = traj_nom.nodes[2,:]

k_x_mat_meas = traj_nom.nodes[1,:]
k_y_mat_meas = traj_nom.nodes[2,:]

k_x_new = G_nom.nodes[1,:]
k_y_new = G_nom.nodes[2,:]

# Plot the 2D spiral trajectory
p2 = Plots.plot(
    k_x_mat, k_y_mat,
    label = "Nominal Matlab constructed trajectory",
    xlabel = "k_x (rad/m)",
    ylabel = "k_y (rad/m)",
    linewidth = 2,
    title = "2D Spiral k-Space Trajectory ",
    aspect_ratio = 1,  # Ensures equal scaling of axes
    legend = :topright

)

plot!(k_x_new, k_y_new,
label = "New trajectory",
xlabel = "k_x (rad/m)",
ylabel = "k_y (rad/m)",
linewidth = 2,
aspect_ratio = 1,  # Ensures equal scaling of axes
legend = :topright)

display(p2)
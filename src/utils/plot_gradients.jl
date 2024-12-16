using Plots

function plot_gradients(time_points, G_x_old, G_y_old, G_x_new, G_y_new; range = :all)
    """
    Plots gradients for old and new gradient arrays.

    Arguments:
    - `time_points`: Vector of time points.
    - `G_x_old`, `G_y_old`: Old gradients (x and y components).
    - `G_x_new`, `G_y_new`: New gradients (x and y components).
    - `range`: Specifies the range of time points to plot (:all or a range, e.g., 1:1000).
    """
    # Select the range of points to plot
    if range == :all
        indices = 1:length(time_points)
    else
        indices = range
    end

    # Subset the data
    times_subset = time_points[indices]
    G_x_old_subset = G_x_old[indices]
    G_y_old_subset = G_y_old[indices]
    G_x_new_subset = G_x_new[indices]
    G_y_new_subset = G_y_new[indices]

    # Plot
    p = Plots.plot(
        times_subset, G_x_old_subset,
        label = "G_x Old",
        xlabel = "Time (s)",
        ylabel = "Gradient Amplitude (mT/m)",
        linewidth = 2,
        title = "Gradients vs Time",
    )
    plot!(times_subset, G_y_old_subset, label = "G_y Old", linewidth = 2)
    plot!(times_subset, G_x_new_subset, label = "G_x New", linewidth = 2)
    plot!(times_subset, G_y_new_subset, label = "G_y New", linewidth = 2)

    display(p)
end
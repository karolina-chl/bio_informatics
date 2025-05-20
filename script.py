from parametter_fitting_final import *


if __name__ == "__main__":

    visualise_model((10**(-6), 2, 0.1, 9, 100, 2.6), "figures/solved_model.png")
    #solve_differential_equation(gc_pathway_exit, (0.1,5,0.1),(10**(-6), 2, 0.1, 9, 100, 2.6))
    #load_the_data("C:/for_python/bio_inf/bio_informatics/data.xlsx")

    # # #visualize the data
    # estimated_parameters_paper = (10**(-6), 2, 0.1, 9, 100, 2.6)
    # visualize_fitting(estimated_parameters_paper, 4, "all data")
    # visualize_data("C:/for_python/bio_inf/bio_informatics/data.xlsx", ["CB", "CC"])
    # visualize_data("C:/for_python/bio_inf/bio_informatics/data.xlsx", ["PC"])

    # # #parameter_fitting
    # bounds = [(0, 100),(0, 100), (0,100), (0,100), (0,100), (0,100)]  # bounds for: μ_p, μ_b, μ_r, σ_p, σ_b, σ_r (in this order)
    # parameter_fitting(bounds, bistability_check = False)

    # # #visualise estimated parameters, bounds equal to [0, 100] for all parameters
    # estimated_parameters0_100 = np.load("estimated_parameters.npy")
    # visualize_fitting(estimated_parameters0_100, 30, "mean")
    # visualize_fitting(estimated_parameters0_100, 4, "all data")

    # # #visualise estimated parameters, smaller bounds 
    # bounds_smaller = [(0, 100),(0, 100), (0,1), (0,100), (0, 100), (0,3)]  # bounds for: μ_p, μ_b, μ_r, σ_p, σ_b, σ_r (in this order)
    # parameter_fitting(bounds_smaller)
    # estimated_parameters = np.load("estimated_parameters_smaller_bounds.npy")
    # visualize_fitting(estimated_parameters, 30, "mean")
    # visualize_fitting(estimated_parameters, 4, "all data")

    #parameter_fitting_bistability
    # bounds = [(0, 100),(0, 100), (0,0.101), (0,100), (0, 100), (1.79,2.62)]  # bounds for: μ_p, μ_b, μ_r, σ_p, σ_b, σ_r (in this order)
    # parameter_fitting(bounds, bistability_check=True)

    # # #visualise estimated parameters, bounds equal to [0, 100] for all parameters, bistability
    # estimated_parameters0_100 = np.load("estimated_parameters_bistability[0,100].npy")
    # visualize_fitting(estimated_parameters0_100, 30, "mean")
    # visualize_fitting(estimated_parameters0_100, 4, "all data")

    # # #visualise estimated parameters, bistability
    # estimated_parameters0_100 = np.load("estimated_parameters_bistability_1.npy")
    # visualize_fitting(estimated_parameters0_100, 30, "mean")
    # visualize_fitting(estimated_parameters0_100, 4, "all data")
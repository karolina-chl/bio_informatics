from parametter_fitting_final import *
from Experiment_bounds import *

def run_experiments():

    #integrate the system of ODE's 
    visualise_model((10**(-6), 2, 0.1, 9, 100, 2.6), "figures/solved_model.png")

    # visualize the data
    estimated_parameters_paper = (10**(-6), 2, 0.1, 9, 100, 2.6)
    visualize_fitting(estimated_parameters_paper, 4, "all data")
    visualize_data("C:/for_python/bio_inf/bio_informatics/data.xlsx", ["CB", "CC"])
    visualize_data("C:/for_python/bio_inf/bio_informatics/data.xlsx", ["PC"])

    # visualise estimated parameters, without bistability
    estimated_parameters_nb = np.load("data/estimated_parameters.npy")
    visualize_fitting(estimated_parameters_nb, 30, "mean")
    #visualize_fitting(estimated_parameters_nb, 4, "all data")

    #visualise the data with bistability 
    estimated_parameters_b = np.load("data/estimated_parameters_bistability.npy")
    visualize_fitting(estimated_parameters_b, 30, "mean")
    #visualize_fitting(estimated_parameters_b, 4, "all data")

    #experiment with the bound
    visualise_the_experiment("data/bounds_experiment_results.xlsx")

if __name__ == "__main__":

    run_experiments()
from parametter_fitting_final import *

np.random.seed(42)

def parametter_fitting_data():

    # # parameter_fitting without the bistability conditions - ok
    # bounds = [(0, 10),(0,10), (0,10), (0,10), (90, 101), (0,10)]  # bounds for: μ_p, μ_b, μ_r, σ_p, σ_b, σ_r (in this order)
    # parameter_fitting(bounds,"data/estimated_parameters.npy", bistability_check = False)

    # parameter_fitting without the bistability conditions - ok
    bounds = [(0, 100),(0,100), (0,100), (0,100), (0, 100), (0,100)]  # bounds for: μ_p, μ_b, μ_r, σ_p, σ_b, σ_r (in this order)
    parameter_fitting(bounds,"data/estimated_parameters_new.npy", bistability_check = False)

    # #parameter fitting with bistability conditions 
    # bounds = [(0, 100),(0, 100), (0,0.101), (0,100), (0, 100), (1.79,2.62)]  # bounds for: μ_p, μ_b, μ_r, σ_p, σ_b, σ_r (in this order)
    # parameter_fitting(bounds, "data/estimated_parameters_bistability.npy", bistability_check=True)

if __name__ == "__main__":
    parametter_fitting_data()
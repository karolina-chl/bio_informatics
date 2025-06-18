from parametter_fitting_final import *

def generate_bounds(true_values, percent):
    bounds = []
    for value in true_values:
        change = value * percent
        lower = max(0, value - change)
        upper = value + change
        bounds.append((lower, upper))
    return bounds

def final_experiment(true_parameter_value, save_location_path, increase_factors=(0.01,0.05, 0.1, 0.2, 0.3,0.4,0.5)):
    results = []

    for percent in increase_factors:
        print(f"\nRunning experiment with ±{percent}% perturbation")

        bounds = generate_bounds(true_parameter_value, percent)

        result = parameter_fitting(bounds, bistability_check=True)
        estimated_params = result.x
        final_error = result.fun

        results.append({
            "Percent": percent,
            "Error": final_error,
            "mu_p": estimated_params[0],
            "mu_b": estimated_params[1],
            "mu_r": estimated_params[2],
            "sigma_p": estimated_params[3],
            "sigma_b": estimated_params[4],
            "sigma_r": estimated_params[5]
        })

    # Save all results to CSV
    data = pd.DataFrame(results)
    data.to_excel(save_location_path, index = False)

    print(f"\nResults saved to: {save_location_path}")

def visualise_the_experiment(location_path, true_params = (1e-6, 2, 0.1, 9, 100, 2.6)):
    mu_p, mu_b,mu_r, sigma_p, sigma_b, sigma_r = true_params
    data = pd.read_excel(location_path)
    results = pd.DataFrame(data)
    error = results["Error"].tolist()
    percent = results["Percent"].tolist()
    mu_p_est = results["mu_p"].tolist()
    #difference_mu_p = [np.abs(value - mu_p) for value in mu_p_est]
    mu_b_est = results["mu_b"].tolist()
    mu_r_est = results["mu_r"].tolist()
    sigma_p_est = results["sigma_p"].tolist()
    sigma_b_est = results["sigma_b"].tolist()
    sigma_r_est = results["sigma_r"].tolist()

    #visualise all parameters 
    fig, ax = plt.subplots(2, 3, figsize=(12, 8))  # Create a 2x3 grid of subplots

    # Plot each set of data on its own subplot
    ax[0, 0].plot(percent, mu_p_est, marker="o")
    ax[0, 0].plot(percent, [mu_p]*len(percent), color = "gray")
    ax[0, 0].set_title("mu_p")

    ax[0, 1].plot(percent, mu_b_est, marker="o")
    ax[0, 1].plot(percent, [mu_b]*len(percent), color = "gray")
    ax[0, 1].set_title("mu_b")

    ax[0, 2].plot(percent, mu_r_est, marker="o")
    ax[0, 2].plot(percent, [mu_r]*len(percent), color = "gray")
    ax[0, 2].set_title("mu_r")

    ax[1, 0].plot(percent, sigma_p_est, marker="o")
    ax[1, 0].plot(percent, [sigma_p]*len(percent), color = "gray")
    ax[1, 0].set_title("sigma_p")

    ax[1, 1].plot(percent, sigma_b_est, marker="o")
    ax[1, 1].plot(percent, [sigma_b]*len(percent), color = "gray")
    ax[1, 1].set_title("sigma_b")

    ax[1, 2].plot(percent, sigma_r_est, marker="o")
    ax[1, 2].plot(percent, [sigma_r]*len(percent), color = "gray")
    ax[1, 2].set_title("sigma_r")

    fig.suptitle("Estimates vs. Percent")
    plt.tight_layout()
    plt.savefig("figures/experiment.png")
    plt.show()
    return 


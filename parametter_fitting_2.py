import numpy as np
import pandas as pd
from scipy.optimize import fsolve, dual_annealing
import matplotlib.pyplot as plt

######################################### Load data and visualize it #################################

def load_expression_data(file_path):
    # Load the Excel file
    df = pd.read_excel(file_path)

    gc_data = df[df['Sample'].isin(['CB', 'CC'])]
    pc_data = df[df['Sample'] == 'PC']

    # Compute mean expression values for each gene
    gc_mean = gc_data[['BLIMP1', 'BCL6', 'IRF4']].mean().values
    pc_mean = pc_data[['BLIMP1', 'BCL6', 'IRF4']].mean().values
 
    return gc_mean, pc_mean

def load_pc_expression(file_path):
    df = pd.read_excel(file_path)
    df.columns = df.columns.str.strip()
    pc_data = df[df['Sample'] == 'PC']
    pc_mean = pc_data[['BLIMP1', 'BCL6', 'IRF4']].mean().values
    return pc_mean

def visualize_data(file_path, Sample): 
    data = pd.read_excel(file_path)
    PC_data = data[data["Sample"].isin(Sample)]
    BLIMP1 = PC_data['BLIMP1']
    BCL6 = PC_data['BCL6']
    IRF4 = PC_data['IRF4']

    #plotting 
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))  

    # Plot each boxplot on a separate subplot
    axes[0].boxplot(BLIMP1)
    axes[0].set_title('BLIMP1')

    axes[1].boxplot(BCL6)
    axes[1].set_title('BCL6')

    axes[2].boxplot(IRF4)
    axes[2].set_title('IRF4')
    fig.suptitle(Sample)
    plt.savefig(f"C:/for_python/bio_inf/bio_informatics/figures/{Sample}")


######################################### simulate the system in steady state #################################


def steady_state_system(params, t_0, u_p, u_b, u_r, sigma_p, sigma_b, sigma_r):
    P, B, R = params
    k_p = k_b = k_r = 1
    lambda_p = lambda_b = lambda_r = 1
    a_bcr = 15
    a_cd = 5
    mean_cd = 60
    mean_bcr = 40
    s = 1.1
    bcr_0 = a_bcr*np.exp(-(t_0 - mean_bcr)**2/s**2)
    cd_0 = a_cd*np.exp(-(t_0 - mean_cd)**2/s**2)
    
    BCR = bcr_0*(k_b**2/(k_b**2+B**2)) 
    CD40 = cd_0*(k_b**2/(k_b**2+B**2))

    dpdt = u_p + sigma_p*(k_b**2/(k_b**2+B**2)) + sigma_p*(R**2/(k_r**2 + R**2)) - lambda_p*P
    dbdt = u_b +sigma_b*(k_p**2/(k_p**2+P**2))*(k_b**2/(k_b**2 +B**2))*(k_r**2/(k_r**2 +R**2)) -(lambda_b + BCR)*B
    drdt = u_r +sigma_r*(R**2/(k_r**2+R**2)) + CD40 - lambda_r*R

    return [dpdt, dbdt, drdt]

def solve_steady_state(mu_p, mu_b, mu_r, sigma_p, sigma_b, sigma_r):
    initial_guess = [9, 2, 2]

    solution, _, ier, _ = fsolve(
        steady_state_system, initial_guess,
        args=(150, mu_p, mu_b, mu_r, sigma_p, sigma_b, sigma_r),
        full_output=True
    )

    if ier != 1:
        raise RuntimeError("Failed to converge to steady state")
    return np.array(solution)

def objective(params, pc_data):
    mu_p, mu_b, mu_r, sigma_p, sigma_b, sigma_r = params

    try:
        #pred_gc = solve_steady_state(mu_p, mu_b, mu_r, sigma_p, sigma_b, sigma_r)
        pred_pc = solve_steady_state(mu_p, mu_b, mu_r, sigma_p, sigma_b, sigma_r)

        #gc_error = np.sum((pred_gc - gc_data) ** 2)
        pc_error = np.sum((pred_pc - pc_data) ** 2)

        return pc_error

    except Exception as e:
        return 1e6

    
######################################### parameter fitting #################################

def run_parameter_fitting(file_path):

    pc_expr = load_pc_expression(file_path)

    # Parameter bounds
    bounds = [(0, 0.1),(1, 2.5), (0,0.101), (0, 10), (90, 100), (1.79, 2.62)]  # bounds for: μ_p, μ_b, μ_r, σ_p, σ_b, σ_r (in this order)

    #global optimization by simulated annealing
    result = dual_annealing(objective, bounds, args=(pc_expr,))

    print("Optimal Parameters:")
    labels = ["mu_p", "mu_b", "mu_r", "sigma_p", "sigma_b", "sigma_r"]
    for name, val in zip(labels, result.x):
        print(f"{name}: {val:.5f}")

    print(f"\nFinal Error: {result.fun:.5f}")
    return result


if __name__ == "__main__":
    visualize_data("C:/for_python/bio_inf/bio_informatics/data.xlsx", ["CB", "CC"])
    fitted = run_parameter_fitting("C:/for_python/bio_inf/bio_informatics/data.xlsx")

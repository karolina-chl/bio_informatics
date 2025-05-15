import numpy as np
import pandas as pd
from scipy.optimize import fsolve, dual_annealing
from scipy.integrate import odeint
from numba import jit
import matplotlib.pyplot as plt

###### load the data from the excel file ############

def load_the_data(file_path, mode = "mean"):

    data = pd.read_excel(file_path)
    if mode == "mean":
        #geminal center relevant data
        gc_data = data[data["Sample"].isin(["CB", "CC"])]
        blimp1_gc = np.mean(gc_data["BLIMP1"])
        bcl6_gc = np.mean(gc_data["BCL6"])
        irf4_gc = np.mean(gc_data["IRF4"])

        #plasma cell relevant data 
        pc_data = data[data["Sample"].isin(["PC"])]
        blimp1_pc = np.mean(pc_data["BLIMP1"])
        bcl6_pc = np.mean(pc_data["BCL6"])
        irf4_pc = np.mean(pc_data["IRF4"])
    
    elif mode == "full":
        #geminal center relevant data
        gc_data = data[data["Sample"].isin(["CB", "CC"])]
        blimp1_gc = gc_data["BLIMP1"]
        bcl6_gc = gc_data["BCL6"]
        irf4_gc = gc_data["IRF4"]

        #plasma cell relevant data 
        pc_data = data[data["Sample"].isin(["PC"])]
        blimp1_pc = pc_data["BLIMP1"]
        bcl6_pc = pc_data["BCL6"]
        irf4_pc = pc_data["IRF4"]
    
    else: 
        print("Data mode not supported. Choose mean or full")

    return blimp1_gc, bcl6_gc, irf4_gc, blimp1_pc, bcl6_pc, irf4_pc


############### solve the differential equations ##################

def gc_pathway_exit(y, t, arguments): 
    u_p, u_b,u_r, sigma_p, sigma_b, sigma_r = arguments
    P,B,R = y
    #constant model parameters (from literature)
    k_p = k_b = k_r = 1
    lambda_p = lambda_b = lambda_r = 1

    #values for BRC and CD40 calculation
    a_bcr = 15
    a_cd = 5
    mean_cd = 60
    mean_bcr = 40
    s = 1.1
    bcr_0 = a_bcr*np.exp(-(t - mean_bcr)**2/s**2)
    cd_0 = a_cd*np.exp(-(t - mean_cd)**2/s**2)
    
    #BCR and CD40 calculations 
    BCR = bcr_0*(k_b**2/(k_b**2+B**2)) 
    CD40 = cd_0*(k_b**2/(k_b**2+B**2))

    #model equations 
    dpdt = u_p + sigma_p*(k_b**2/(k_b**2+B**2)) + sigma_p*(R**2/(k_r**2 + R**2)) - lambda_p*P
    dbdt = u_b +sigma_b*(k_p**2/(k_p**2+P**2))*(k_b**2/(k_b**2 +B**2))*(k_r**2/(k_r**2 +R**2)) -(lambda_b + BCR)*B
    drdt = u_r +sigma_r*(R**2/(k_r**2+R**2)) + CD40 - lambda_r*R

    return dpdt, dbdt, drdt

def solve_differential_equation(gc_pathway_exit, arguments):
    initial_conditions = (0.1,5,0.1)
    #solve 
    time = np.linspace(0,200,200)
    solution = odeint(gc_pathway_exit, initial_conditions, time, args = (arguments,))
    P, B, R = solution.T
    gc_equilibrium = {"BLIMP1": P[25], "BCL6": B[25], "IRF4":R[25]}
    pc_equilibrium = {"BLIMP1": P[175], "BCL6": B[175], "IRF4":R[175]}

    return gc_equilibrium, pc_equilibrium, P, B, R, time

def error_function(arguments):
    gc_equilibrium, pc_equilibrium, _, _, _, _= solve_differential_equation(gc_pathway_exit, arguments)
    blimp1_gc, bcl6_gc, irf4_gc, blimp1_pc, bcl6_pc, irf4_pc = load_the_data("C:/for_python/bio_inf/bio_informatics/data.xlsx")

    error_gc = np.sqrt((gc_equilibrium["BLIMP1"] - blimp1_gc)**2) + np.sqrt((gc_equilibrium["BCL6"] - bcl6_gc)**2) + np.sqrt((gc_equilibrium["IRF4"] - irf4_gc)**2)
    error_pc = np.sqrt((pc_equilibrium["BLIMP1"] - blimp1_pc)**2) + np.sqrt((pc_equilibrium["BCL6"] - bcl6_pc)**2) + np.sqrt((pc_equilibrium["IRF4"] - irf4_pc)**2)

    overall_error = error_gc + error_pc

    return overall_error
    

################# fitting function ###########################

def parameter_fitting():    
    bounds = [(0, 1),(1, 3), (0,0.101), (8,10), (90, 101), (1.79, 2.62)]  # bounds for: μ_p, μ_b, μ_r, σ_p, σ_b, σ_r (in this order)
    result = dual_annealing(error_function, bounds)

    print("Optimal Parameters:")
    labels = ["mu_p", "mu_b", "mu_r", "sigma_p", "sigma_b", "sigma_r"]
    np.save("estimated_parameters.npy", result.x)
    
    for name, val in zip(labels, result.x):
        print(f"{name}: {val}")
    print(f"\nFinal Error: {result.fun:.5f}")

    return result

#################### visualize the model with the estimated parameters #########################

def visualize(arguments, dotsize, mode = "all data"):
    _,_, P, B, R, time = solve_differential_equation(gc_pathway_exit, arguments)

    #plotting the model with estimated parameters
    plt.figure()
    plt.plot(time, P, "blue", label="BLIMP1")
    plt.plot(time, B, "green", label = "BCL6")
    plt.plot(time, R, "yellow", label = "IRF4")
    plt.xlabel("Time")
    plt.legend()

    if mode == "all data":
        #plotting the data on the figure 
        blimp1_gc, bcl6_gc, irf4_gc, blimp1_pc, bcl6_pc, irf4_pc = load_the_data("C:/for_python/bio_inf/bio_informatics/data.xlsx", "full")
        plt.scatter([25]*len(blimp1_gc), blimp1_gc, marker ="o", color = "blue", s=dotsize)
        plt.scatter([25]*len(bcl6_gc), bcl6_gc, marker ="o", color = "green", s=dotsize)
        plt.scatter([25]*len(irf4_gc), irf4_gc, marker ="o", color = "yellow", s=dotsize)
        plt.scatter([175]*len(blimp1_pc), blimp1_pc, marker = "o", color = "blue", s=dotsize)
        plt.scatter([175]*len(bcl6_pc), bcl6_pc, marker ="o", color = "green", s=dotsize)
        plt.scatter([175]*len(irf4_pc), irf4_pc, marker ="o", color = "yellow", s=dotsize)
        plt.savefig(f"C:/for_python/bio_inf/bio_informatics/figures/{mode}.jpg")
        plt.show()
    elif mode == "mean":
        #adding mean on the figure 
        blimp1_gc, bcl6_gc, irf4_gc, blimp1_pc, bcl6_pc, irf4_pc = load_the_data("C:/for_python/bio_inf/bio_informatics/data.xlsx")
        plt.scatter(25, blimp1_gc, marker ="*", color = "blue", s=dotsize)
        plt.scatter(25, bcl6_gc, marker ="*", color = "green", s=dotsize)
        plt.scatter(25, irf4_gc, marker ="*", color = "yellow", s=dotsize)
        plt.scatter(175, blimp1_pc, marker = "*", color = "blue", s=dotsize)
        plt.scatter(175, bcl6_pc, marker ="*", color = "green", s=dotsize)
        plt.scatter(175, irf4_pc, marker ="*", color = "yellow", s=dotsize)
        plt.savefig(f"C:/for_python/bio_inf/bio_informatics/figures/{mode}, {blimp1_gc}.jpg")
        plt.show()


if __name__ == "__main__":
    #solve_differential_equation(gc_pathway_exit, (0.1,5,0.1),(10**(-6), 2, 0.1, 9, 100, 2.6))
    #load_the_data("C:/for_python/bio_inf/bio_informatics/data.xlsx")
    parameter_fitting()
    estimated_parameters = np.load("estimated_parameters.npy")
    visualize(estimated_parameters, 30, "mean")
    visualize(estimated_parameters, 4, "all data")


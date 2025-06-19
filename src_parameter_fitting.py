import numpy as np
import pandas as pd
from scipy.optimize import dual_annealing
from scipy.integrate import odeint
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

def solve_differential_equation(gc_pathway_exit, arguments, initial_conditions = (0.1,5,0.1)):
    #solve 
    time = np.linspace(0,200,200)
    solution = odeint(gc_pathway_exit, initial_conditions, time, args = (arguments,))
    P, B, R = solution.T
    gc_equilibrium = {"BLIMP1": P[25], "BCL6": B[25], "IRF4":R[25]}
    pc_equilibrium = {"BLIMP1": P[175], "BCL6": B[175], "IRF4":R[175]}

    return gc_equilibrium, pc_equilibrium, P, B, R, time

def visualise_model(arguments, save_path):
    _, _,  P, B, R, _ = solve_differential_equation(gc_pathway_exit, arguments)
    time = np.linspace(0,200,200)
    plt.figure(figsize=(8,4))
    plt.plot(time, P, "blue", label="BLIMP1")
    plt.plot(time, B, "green", label = "BCL6")
    plt.plot(time, R, "yellow", label = "IRF4")
    plt.xlabel("Time")
    plt.ylabel("Level")
    plt.legend()
    plt.savefig(save_path)
    plt.show()


def error_function(arguments):
    
    gc_equilibrium, pc_equilibrium, _, _, _, _= solve_differential_equation(gc_pathway_exit, arguments)
    blimp1_gc, bcl6_gc, irf4_gc, blimp1_pc, bcl6_pc, irf4_pc = load_the_data("C:/for_python/bio_inf/bio_informatics/data.xlsx")

    error_gc = np.sqrt((gc_equilibrium["BLIMP1"] - blimp1_gc)**2) + np.sqrt((gc_equilibrium["BCL6"] - bcl6_gc)**2) + np.sqrt((gc_equilibrium["IRF4"] - irf4_gc)**2)
    error_pc = np.sqrt((pc_equilibrium["BLIMP1"] - blimp1_pc)**2) + np.sqrt((pc_equilibrium["BCL6"] - bcl6_pc)**2) + np.sqrt((pc_equilibrium["IRF4"] - irf4_pc)**2)

    overall_error = error_gc + error_pc

    return overall_error

def error_function_bistability(arguments):

    _, _,u_r, _, _, sigma_r = arguments

    if bistability_conditions(u_r,sigma_r) == False:
        return 1e6
    
    gc_equilibrium, pc_equilibrium, _, _, _, _= solve_differential_equation(gc_pathway_exit, arguments)
    blimp1_gc, bcl6_gc, irf4_gc, blimp1_pc, bcl6_pc, irf4_pc = load_the_data("C:/for_python/bio_inf/bio_informatics/data.xlsx")

    error_gc = np.sqrt((gc_equilibrium["BLIMP1"] - blimp1_gc)**2) + np.sqrt((gc_equilibrium["BCL6"] - bcl6_gc)**2) + np.sqrt((gc_equilibrium["IRF4"] - irf4_gc)**2)
    error_pc = np.sqrt((pc_equilibrium["BLIMP1"] - blimp1_pc)**2) + np.sqrt((pc_equilibrium["BCL6"] - bcl6_pc)**2) + np.sqrt((pc_equilibrium["IRF4"] - irf4_pc)**2)

    overall_error = error_gc + error_pc

    return overall_error

################# fitting function ###########################

def bistability_conditions(u_r, sigma_r):
    CD40 = 0 
    lambda_r = k_r = 1
    beta = (u_r + CD40 + sigma_r)/lambda_r*k_r

    if beta < np.sqrt(3):
        return False
    
    condition_S10 = beta**3 - np.sqrt((beta**2 -3)**3) + 9*beta
    condition_S11 = beta**3 + np.sqrt((beta**2 -3)**3) + 9*beta
    treshold = (27*sigma_r)/(2*lambda_r*k_r)

    if  condition_S10 < treshold and condition_S11 > treshold: 
        return True
    else:
        return False 
    
def parameter_fitting(bounds, file_name, bistability_check = True):

    if bistability_check:
        result = dual_annealing(error_function_bistability, bounds)
    else: 
        result = dual_annealing(error_function, bounds)

    # Save both parameters and error
    results_to_save = {
        "parameters": result.x,
        "error": result.fun
    }
    np.save(file_name, results_to_save)

    labels = ["mu_p", "mu_b", "mu_r", "sigma_p", "sigma_b", "sigma_r"]
    print("Optimal Parameters:")
    for name, val in zip(labels, result.x):
        print(f"{name}: {val}")
    print(f"\nFinal Error: {result.fun:.5f}")

    return result

#################### visualize the model with the estimated parameters #########################

def visualize_fitting(arguments, dotsize, mode = "all data"):
    _,_, P, B, R, time = solve_differential_equation(gc_pathway_exit, arguments)

    #plotting the model with estimated parameters
    plt.figure(figsize=(8,4))
    plt.plot(time, P, "blue", label="BLIMP1")
    plt.plot(time, B, "green", label = "BCL6")
    plt.plot(time, R, "yellow", label = "IRF4")
    plt.xlabel("Time")
    plt.ylabel("Level")
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
        plt.savefig(f"C:/for_python/bio_inf/bio_informatics/figures/{mode}, {arguments[1]}.jpg")
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
        plt.savefig(f"C:/for_python/bio_inf/bio_informatics/figures/{mode}, {arguments[1]}.jpg")
        plt.show()

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

    for ax in axes:
        ax.set_xticklabels([])
    plt.savefig(f"C:/for_python/bio_inf/bio_informatics/figures/{Sample}")
    plt.show()
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

data = pd.read_excel("C:/for_python/bio_inf/bio_informatics/data.xlsx")


#parameters 
C_0 = 10**(-8) #unit of concentration
t_0 = 4 #the unit of time set to a typical protein half-life
u_p = 10**(-6) # u stand for basal production rate of each protein
u_b = 2
u_r = 0.1
sigma_p = 9
sigma_b = 100
sigma_r = 2.6
k_p = 1
k_b = 1
k_r = 1
lambda_p = 1 #lambda degradation rate
lambda_b = 1
lambda_r = 1
bcr_0 = 5 #not constant
cd_0 = 0.09

#Initial conditions 
B_0 = 5 #protein level of BCL6
R_0 = 0.1 #protein level of IRF4
P_0 = 0.1  #protein level of BLIMP1

a_bcr = 15
a_cd = 5
mean_cd = 60
mean_bcr = 40
s = 1.1

time = np.linspace(0,200,100)

def gc_steady_state(y, u_p, u_b, u_r, sigma_p, sigma_b, sigma_r, BCR, CD40):
    P, B, R = y
    
    # Constants
    k_b = k_p = k_r = 1
    lambda_p = lambda_b = lambda_r = 1

    # BCL6-dependent modulation
    dpdt = u_p + sigma_p * (k_b**2 / (k_b**2 + B**2)) + sigma_p * (R**2 / (k_r**2 + R**2)) - lambda_p * P
    dbdt = u_b + sigma_b * (k_p**2 / (k_p**2 + P**2)) * (k_b**2 / (k_b**2 + B**2)) * (k_r**2 / (k_r**2 + R**2)) - (lambda_b + BCR) * B
    drdt = u_r + sigma_r * (R**2 / (k_r**2 + R**2)) + CD40 - lambda_r * R
    return [dpdt, dbdt, drdt]


#initial conditions 
y_0 = (P_0, B_0, R_0)

params = {
    'u_p': 1e-6,
    'u_b': 2.0,
    'u_r': 0.1,
    'sigma_p': 9.0,
    'sigma_b': 100.0,
    'sigma_r': 2.6,
    'BCR': 5.0,
    'CD40': 0.5
}

# Solve
steady_state = fsolve(gc_steady_state, y_0,
                      args=(params['u_p'], params['u_b'], params['u_r'],
                            params['sigma_p'], params['sigma_b'], params['sigma_r'],
                            params['BCR'], params['CD40']))

print("Steady-state values [BLIMP1, BCL6, IRF4]:", steady_state)
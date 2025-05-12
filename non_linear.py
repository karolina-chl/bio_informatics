#Parameter fitting using non-linear least squares (or another method) and using the bistability conditions provided by equations S9, S10 and S11.
import numpy as np
from scipy.optimize import least_squares

def bistability_conditions(CD40, u_r, sigma_r):
    lambda_r = k_r = 1
    beta = (u_r + CD40 + sigma_r)/lambda_r*k_r
    condition_S10 = beta**3 - np.sqrt((beta**2 -3)**3) + 9*beta
    condition_S11 = beta**3 + np.sqrt((beta**2 -3)**3) + 9*beta
    treshold = (27*sigma_r)/(2*lambda_r*k_r)

    if beta > np.sqrt(3): 
        print("Condition S9 met")
    else:
        print("Condition S9 NOT met")
    if condition_S10 < treshold: 
        print("Condition S10 met") 
    else:
        print("Condition S10 NOT met")
    if condition_S11 > treshold:
        print("Condition S11 met") 
    else:
        print("Condition S11 NOT met")

bistability_conditions(1,0.1,2.62)

def IRF4(u_r, sigma_r,R, CD40):
    lambda_r = k_r = 1
    return u_r +sigma_r*(R**2/(k_r**2+R**2)) + CD40 - lambda_r*R
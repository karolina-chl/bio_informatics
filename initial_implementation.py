import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

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
bcr_0 = 5*(C_0/t_0)
cd_0 = 0.09

#check
critical_value = sigma_b/lambda_b
if critical_value >20: print("True")

#Initial conditions 
B_0 = 5 #protein level of BCL6
R_0 = 0.1 #protein level of IRF4
P_0 = 0.1  #protein level of BLIMP1

#beta
CD40 = cd_0*(k_b**2/(k_b**2+B_0**2))
beta = (u_r + CD40 + sigma_r)/lambda_r*k_r
condition = beta**3 - np.sqrt((beta**2 -3)**3) + 9*beta
treshold = (27*sigma_r)/(2*lambda_r*k_r)
if beta > np.sqrt(3): print("True")
if condition < (27*sigma_r)/(2*lambda_r*k_r): 
    print("True") 
else:
    print("False")
print(beta)
print(f"Conditions is {condition} and theshold is {treshold}")




#time
time = np.linspace(0,200,100)

BCR = bcr_0*(k_b**2/(k_b**2+B_0**2)) #constant or changing?
CD40 = cd_0*(k_b**2/(k_b**2+B_0**2))

def gc_pathway_exit(y,t, BCR, CD40):
    P,B,R = y
    dpdt = u_p + sigma_p*(k_b**2/(k_b**2+B**2)) + sigma_p*(R**2/(k_r**2 + R**2)) - lambda_p*P
    dbdt = u_b +sigma_b*(k_p**2/(k_p**2+P**2))*(k_b**2/(k_b**2 +B**2))*(k_r**2/(k_r**2 +R**2)) -(lambda_b + BCR)*B
    drdt = u_r +sigma_r*(R**2/(k_r**2+R**2)) + CD40 - lambda_r*R

    return dpdt, dbdt, drdt

#initial conditions 
y_0 = (P_0, B_0, R_0)

#solving 
solution = odeint(gc_pathway_exit, y_0, time, (BCR, CD40))
P,B,R = solution.T

#plotting 

plt.figure()
plt.plot(time, P, "blue", label="BLIMP1")
plt.plot(time, B, "green", label = "BCL6")
plt.plot(time, R, "yellow", label = "IRF4")
plt.legend()
plt.show()

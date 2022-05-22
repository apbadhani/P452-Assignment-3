import math
import numpy as np
import matplotlib.pyplot as plt
from MNL import *

#This is a 2-point boundary value problem
#where there is a infinite well in schrodinger's equation
#is to be solved using shooting method with RK4 integrator to
#the solution in the tolerance limit

#---Potential function---#
def V(x):
    return 0.0
#Constants
m = 9.10e-31  #mass of electron
h_b = 1.05e-34 # h/2*pi
e = 1.6e-19
n = 1000
l = 5.29e-11 #Bohr radius
h = l/n

#Defining Shooter function with RK4 integrator
#to solve the boundary value problems
#As an input it takes energy intervals and nodes
#and returns the normalized psi, it's energy.
def ode_shoot_rk4(e_lim, nodes):
    psi_0 = 0.0
    phi_0 = 1.0
    psi_in = np.array([psi_0, phi_0])
    h_mesh = 1.0/100.0
    x_r = np.arange(0.0, 1.0+h_mesh, h_mesh)
    v_ip = np.zeros(len(x_r))
    e1tu, e2tu = tune_energ(e_lim[0], e_lim[1], nodes, psi_in, x_r, v_ip)
    psi = rk4(Rhs_schr, psi_in, x_r, v_ip, e2tu)[:, 0]
    return e1tu, norm_ed(psi), x_r

#Function that defines Schrodinger equation to solve
#It returns an array containing psi and RHS to Schrodinger eqn.
def Rhs_schr(y, r, V, E):
    psi, phi = y
    dphidx = [phi, (V-E)*psi]
    return np.array(dphidx)

#This function is used to to refine or tune energy values after every step
#To initilize it takes input of refining variable's upper and lower limit
#of energy, nodes, initial psi0 and returns the tuned energy limits.
def tune_energ(Nodes,E_low,E_up, psi0, x, V):
#Tolerance
    tol = 1e-12
#Upper limit
    E1 = E_up
#Lower limit
    E2 = E_low
#Initialization
    psi = [1]
#This loop runs till the refined energy limit difference is less than
#the tolerance value or till the absolute value is greater than 1e-3.
    while (abs(E2-E1) > tol or abs(psi[-1]) > 1e-3):
        iE = (E1+E2)/2.0
#Calling the RK4 integrator function
        psi = rk4(Rhs_schr, psi0, x, V, iE)[:, 0]
        node_s = len(findZeros(psi))-1
        if node_s > Nodes+1:
            E1 = iE
            continue
        if node_s < Nodes-1:
            E2 = iE
            continue
        if (node_s % 2 == 0):
            if((psi[len(psi)-1] <= 0.0)):
                E1 = iE
            else:
                E2 = iE
        elif node_s > 0:
            if((psi[len(psi)-1] <= 0.0)):
                E2 = iE
            else:
                E1 = iE
        elif node_s < 0:
            E2 = iE
    return E2, E1


e_initial = [1.0, 500.0]
nodes_arr = [0, 1]
L = 0.0
N = 1.0
figipw = plt.figure()
e1, psi_ipw1, x_ipw1 = ode_shoot_rk4(e_initial, 0)
print('Energy of ground state:',e1)
e2, psi_ipw, x_ipw = ode_shoot_rk4(e_initial, 2)
print('Energy of first excited state:',e2)

#-----Plotting----#
plt.plot(x_ipw1, psi_ipw1,color="red")
plt.plot(x_ipw, psi_ipw,color="green")
plt.title('$\psi(x)$ for ground state and first excited state')
plt.xlabel('x')
plt.ylabel(r'$\psi$')
plt.show()






import matplotlib.pyplot as plt
import numpy as np

#Boundary Conditions for Laplace Equation solver
p_t = 0 # y = 1
p_b = 1 #y = 0
p_l = 0 #x=1
p_r = 0 #x=0
#Guess phi
p_g = 0.5
#-----Solves Laplace eqn----#
def lapla_sol(phi, tol=1e-6):
    dif = 1.0
    while (dif > tol):
        phi1 = phi.copy()
        phi[1:-1, 1:-1] = (phi1[1:-1, :-2] + phi1[1:-1, 2:] + phi1[:-2, 1:-1] + phi1[2:, 1:-1]) / 4
        dif = np.sqrt(np.sum((phi - phi1) ** 2) / np.sum(phi1 ** 2))
    return phi
#Initializing x,y vectors
x = np.linspace(0, 1, 10)
y = np.linspace(0, 1, 10)
#Initializing potential as 10x10 matrix
pot = np.zeros((10, 10))
#Boundary condition for potential
pot[:, 0] = 1
#Solving the Laplace equation
pot = lapla_sol(pot)
#---Plotting the solution of Laplace equation as a contour---#
X, Y = np.meshgrid(x, y)
plt.title("Voltage contour ")
plt.contourf(X/10, Y/10, pot, 50, cmap="jet")

# Set Colorbar
plt.colorbar()
plt.show()
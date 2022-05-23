import math
import matplotlib.pyplot as plt

#Defining Schrodinger equation
#takes V, E as an input and returns an array
def sEqn(y, r, V, E):
    psi, phi = y
    dphidx = [phi, (V-E)*psi]
    return np.array(dphidx)
#---Normalization of psi---#
def norm(psi):
    return psi / np.max(psi)
#Defining Shooter function RK4 integrator
def shoot_with_rk4(f, a, b, y0, y1):
    tol = 1e-5
    h = 0.01
    z_i = 2
    X, Y = rk42o(f, a, y0, z_i, b, h)
    if abs(Y[-1]-y1) <= tol:
        return X, Y
    else:
        z_n = z_i
        y_n = Y[-1]
        z0 = float(input())
        X, Y = rk42o(f, a, y0, z0, b, h)
        if (y_n > y1 and Y[-1] < y1) or (y_n < y1 and Y[-1] > y1):
            z_i = z_i + (z_n-z_i)*(y1 - Y[-1])/(y_n - Y[-1])
            X,Y = rk42o(f, a, y0, z0, b, h)
        return X, Y



#----Plotting----#





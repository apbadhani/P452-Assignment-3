import matplotlib as plt
import numpy as np

#-----Solves Laplace eqn----#
def pde_lapla_solve(u_bou,L,N):
    h = L/(N+1)
    u = [[0 for i in range(N+2)] for j in range(N+2)]
    for i in range(N+2):
        u[0][i] = u_bou(0, i*h)[0][0]
        u[N+1][i] = u_bou(L, i*h)[0][1]
        u[i][0] = u_bou(i*h, 0)[1][0]
        u[i][N+1] = u_bou(i*h, L)[1][1]
    for iter in range(1000):
        for i in range(1,N+1):
            for j in range(1, N+1):
                u[i][j] = (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1])/4
    plt.imshow(u)
    plt.show()

    return u

x = np.linspace(0, 1, 10)
y = np.linspace(0, 1, 10)
#potential
pot = np.zeros((10, 10))
#Boundary condition for potential
pot[:, 0] = 1

pot = pde_lapla_solve(pot,tol=1e-6)


#---Plotting---#
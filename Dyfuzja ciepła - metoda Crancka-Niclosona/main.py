import matplotlib.pyplot as plt
from numba import njit
import numpy as np
import scipy.linalg
import scipy



def zad():
    delta = 1
    dt = 1
    TA = 40.0
    TB = 0.0
    TD = 0.0
    TC = 30.0
    KB = 0.1
    KD = 0.6
    nx = 40
    ny = 40
    N = nx*ny + nx + ny + 1
    IT_MAX = 2000
    A = np.zeros((N, N))
    B = np.zeros((N, N))
    C = np.zeros(N)
    for i in range(1, nx):
        for j in range(1, ny):
            l = i + j * (nx+1)
            A[l][l-nx-1] = dt/(2*delta**2)
            A[l][l-1] = dt/(2*delta**2)
            A[l][l+1] = dt/(2*delta**2)
            A[l][l+nx+1] = dt/(2*delta**2)
            A[l][l] = -2.0*dt/(delta**2) - 1
            B[l][l - nx - 1] = -dt / (2 * delta ** 2)
            B[l][l - 1] = -dt / (2 * delta ** 2)
            B[l][l + 1] = -dt / (2 * delta ** 2)
            B[l][l + nx + 1] = -dt / (2 * delta ** 2)
            B[l][l] = 2.0 * dt / (delta ** 2) - 1


    print("\n\n\n\n")
    for j in range(ny+1):
        l = j*(nx+1)
        A[l][l] = 1.0
        B[l][l] = 1.0
        C[l] = 0.0


    for j in range(ny+1):
        l = nx + j*(nx+1)
        A[l][l] = 1.0
        B[l][l] = 1.0
        C[l] = 0.0


    for i in range(1, nx):
        l = i + ny*(nx+1)
        A[l][l-nx-1] = -1.0/(KB*delta)
        A[l][l] = 1 + 1.0/(KB*delta)
        C[l] = TB
        for k in range(N):
            B[l][k] = 0.0

    for i in range(1, nx):
        l = i
        A[l][l+nx+1] = -1.0/(KD*delta)
        A[l][l] = 1 + 1/(KD*delta)
        C[l] = TD
        for k in range(N):
            B[l][k] = 0.0

    T = np.zeros(N)
    for j in range(ny+1):
        T[j*(nx+1)] = TA
        T[nx+j*(nx+1)] = TC

    P = np.zeros(N)
    LU = np.zeros((N,N))
    LU, P = scipy.linalg.lu_factor(A)
    print(LU.shape)
    for it in range(IT_MAX+1):
        d = B@T + C
        T1 = scipy.linalg.lu_solve((LU, P), d)
        if it == 100 or it == 200 or it == 500 or it == 1000 or it == 2000:
            plt.pcolor(T.reshape(nx+1, ny+1))
            plt.xlabel("x")
            plt.ylabel("y")
            plt.title("Temperatura it = " + str(it))
            plt.colorbar()
            plt.show()
            plt.pcolor(((T1 - T)/dt).reshape(nx+1, ny+1))
            plt.xlabel("x")
            plt.ylabel("y")
            plt.title("Dyfuzja it = " + str(it))
            plt.colorbar()
            plt.show()
        T = T1


zad()



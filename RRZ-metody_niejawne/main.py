import matplotlib.pyplot as plt
import numpy as np
import math


beta = 0.001
N = 500
gamma = 0.1
t_max = 100
dt = 0.1
maxiter = 20
TOL = 10**-6
u0 = 1
alfa = beta*N - gamma

def f(u):
    return (beta*N - gamma)*u - beta*u**2


def fpoch(u):
    return beta*N - gamma - 2*beta*u


def piccard_method():
    t = 0.1
    u1 = u0
    u_arr = []
    z_arr = []
    u_arr.append(u1)
    z_arr.append(N - u1)
    while t < t_max:
        mi_it = 0
        u = u1
        umax = 0.0
        while math.fabs(u1-umax) > TOL and mi_it < maxiter:
            mi_it += 1
            umax = u1
            u1 = u + (dt/2)*(f(u1)+f(umax))
        t += dt
        u_arr.append(u1)
        z_arr.append(N - u1)
    u_arr.pop()
    z_arr.pop()
    xaxis = np.linspace(0, t_max, int(t_max / dt))
    plt.plot(xaxis, u_arr, label='u(t)')
    plt.plot(xaxis, z_arr, label='z(t)')
    plt.legend()
    plt.show()

def newton_iter():
    t = 0.1
    u1 = u0
    u_arr = []
    z_arr = []
    u_arr.append(u1)
    z_arr.append(N - u1)
    while t < t_max:
        mi_it = 0
        mi_it += 1
        u = u1
        umax = 0.0
        while math.fabs(u1-umax) > TOL and mi_it < maxiter:
            umax = u1
            u1 = umax - (umax - u - dt/2*(f(u1)+f(umax)))/(1 - (dt/2)*fpoch(umax))
        t += dt
        u_arr.append(u1)
        z_arr.append(N - u1)
    u_arr.pop()
    z_arr.pop()
    xaxis = np.linspace(0, t_max, int(t_max / dt))
    plt.plot(xaxis, u_arr, label='u(t)', color = 'red')
    plt.plot(xaxis, z_arr, label='z(t)', color = 'green')
    plt.legend()
    plt.show()


def mxy(u1, a11):
    return (-1)*dt*a11*(alfa - 2*beta*u1)


def mxx(u1, a11):
    return 1 - dt * a11 * (alfa - 2 * beta * u1)


def F1(U1, U2, u,a11, a12):
    return U1 - u - dt * (a11 * (alfa * U1 - beta * U1 * U1) + a12 * (alfa * U2 - beta * U2 * U2))


def F2(U1, U2, u,a21, a22):
    return U2 - u - dt * (a21 * (alfa * U1 - beta * U1 * U1) + a22 * (alfa * U2 - beta * U2 * U2))



def methodRKK():
    a11 = 0.25
    a12 = 0.25 - math.sqrt(3)/6
    a21 = 0.25 + math.sqrt(3)/6
    a22 = a11
    b = 0.5
    u = u0
    u1 = u
    u_arr = []
    z_arr = []
    u_arr.append(u1)
    z_arr.append(N - u1)
    t = 0.1
    while t < t_max:
        mi_it = 0
        mi_it += 1
        u = u1
        U1tmp = 0.0
        U2tmp = 0.0
        U1 = u
        U2 = u
        while (math.fabs(U1 - U1tmp) > TOL or math.fabs(U2 - U2tmp) > TOL) and mi_it < maxiter:
            U1tmp = U1
            U2tmp = U2
            dtU1 = (F2(U1, U2, u, a21, a22)*mxy(U2, a12) - F1(U1, U2, u, a11, a12)*mxx(U2, a22))/(mxx(U1, a11)*mxx(U2, a22) - mxy(U2, a12)*mxy(U1, a21))
            dtU2 = (F1(U1, U2, u, a11, a12)*mxy(U1, a21)-F2(U1, U2, u, a21, a22)*mxx(U1, a11))/(mxx(U1, a11)*mxx(U2, a22)-mxy(U2, a12)*mxy(U1, a21))
            U1 = U1tmp + dtU1
            U2 = U2tmp + dtU2
            u1 = u + dt * (b * f(U1) + b * f(U2))
        u_arr.append(u1)
        z_arr.append(N - u1)
        t += dt
    u_arr.pop()
    z_arr.pop()
    xaxis = np.linspace(0, t_max, int(t_max / dt))
    plt.plot(xaxis, u_arr, label='u(t)', color='violet')
    plt.plot(xaxis, z_arr, label='z(t)', color='blue')
    plt.legend()
    plt.show()


piccard_method()
newton_iter()
methodRKK()

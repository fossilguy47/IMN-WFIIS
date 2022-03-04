import matplotlib.pyplot as plt
import numpy as np
import math


x0 = 0.01
v0 = 0.0
dt0 = 1.0
S = 0.75
p = 2.0


def g(x, v, alfa):
    return alfa * (1 - x**2) * v - x


def RK2_method(x, v, dt, alfa):
    k1x = v
    k1v = g(x, v, alfa)
    k2x = v + dt*k1v
    k2v = alfa * (1 - (x+dt*k1x)**2) * (v+dt*k1v) - (x+dt*k1x)
    x_re = x + dt/2*(k1x+k2x)
    v_re = v + dt/2*(k1v+k2v)
    return x_re, v_re


def F(x1, x, v1, v, dt):
    return x1 - x - dt/2*(v + v1)


def G(x1, x, v1, v, dt, alfa):
    return v1 - v - dt / 2 * (g(x, v, alfa) + g(x1, v1, alfa))


a11 = 1.0


def a12(dt):
    return (-1)*dt/2


def a21(dt, x, v, alfa):
    return (-1)*dt/2*((-2)*alfa*x*v-1)


def a22(dt, x, alfa):
    return 1 - dt/2*alfa*(1-x**2)


def trapezMethod(x, v, dt, alfa):
    delta = 10**-10
    xn1 = x
    vn1 = v
    dx = 1
    dv = 1
    while math.fabs(dx) > delta or math.fabs(dv) > delta:
        vn1k = vn1
        xn1k = xn1
        dx = ((-1) * F(xn1, x, vn1, v, dt) * a22(dt, xn1k, alfa) - (-1) * G(xn1, x, vn1, v, dt, alfa) * a12(dt))/(a11 * a22(dt, xn1k, alfa) - a12(dt) * a21(dt, xn1k, vn1k, alfa))
        dv = (a11 * (-1) * G(xn1, x, vn1, v, dt, alfa) - a21(dt, xn1k, vn1k, alfa) * (-1) * F(xn1, x, vn1, v, dt)) /(a11 * a22(dt, xn1k, alfa) - a12(dt) * a21(dt, xn1k, vn1k, alfa))
        xn1 = xn1k + dx
        vn1 = vn1k + dv
    x_res = x + (dt / 2) * (v + vn1)
    v_res = v + (dt / 2) * (g(x, v, alfa) + g(xn1, vn1, alfa))

    return x_res, v_res


def krok_count(way_of_solution, TOL, alfa):
    dt = dt0
    t = dt
    x = x0
    v = v0
    tmax = 40.0
    x_res_arr = []
    v_res_arr = []
    dt_res_arr = []
    t_res_arr = []
    dt_res_arr.append(dt)
    t_res_arr.append(t)
    x_res_arr.append(x)
    v_res_arr.append(v)
    while t < tmax:
        tmpx1, tmpv1 = way_of_solution(x, v, dt, alfa)
        tmpx1, tmpv1 = way_of_solution(tmpx1, tmpv1, dt, alfa)
        tmpx2, tmpv2 = way_of_solution(x, v, 2*dt, alfa)
        Ex = (tmpx1 - tmpx2)/(2**p - 1)
        Ev = (tmpv1 - tmpv2)/(2**p - 1)
        if max(math.fabs(Ex), math.fabs(Ev)) < TOL:
            t = t + 2*dt
            x = tmpx1
            v = tmpv1
            x_res_arr.append(x)
            v_res_arr.append(v)
            dt_res_arr.append(dt)
            t_res_arr.append(t)
        dt = ((S*TOL/max(math.fabs(Ex), math.fabs(Ev)))**(1.0/(p+1.0)))*dt

    return dt_res_arr, t_res_arr, x_res_arr, v_res_arr

#zad1
tmax = 40.0
alfa = 5
dt_res_arr1, t_res_arr1, x_res_arr1, v_res_arr1 = krok_count(RK2_method, 10**-2, alfa)
dt_res_arr2, t_res_arr2, x_res_arr2, v_res_arr2 = krok_count(RK2_method, 10**-5, alfa)
plt.plot(t_res_arr1, x_res_arr1, label='TOL=10^-2', color='red')
plt.plot(t_res_arr2, x_res_arr2, label='TOL=10^-5', color='green')
plt.legend()
plt.title("RK2 x(t)")
plt.show()

plt.plot(t_res_arr1, v_res_arr1, label='TOL=10^-2', color='red')
plt.plot(t_res_arr2, v_res_arr2, label='TOL=10^-5', color='green')
plt.legend()
plt.title("RK2 v(t)")
plt.show()

plt.plot(x_res_arr1, v_res_arr1, label='TOL=10^-2', color='red')
plt.plot(x_res_arr2, v_res_arr2, label='TOL=10^-5', color='green')
plt.legend()
plt.title("RK2 V(x)")
plt.show()

plt.plot(t_res_arr1, dt_res_arr1, label='TOL=10^-2', color='red')
plt.plot(t_res_arr2, dt_res_arr2, label='TOL=10^-5', color='green')
plt.legend()
plt.title("RK2 dt(t)")
plt.show()


dt_res_arr1, t_res_arr1, x_res_arr1, v_res_arr1 = krok_count(trapezMethod, 10**-2, alfa)
dt_res_arr2, t_res_arr2, x_res_arr2, v_res_arr2 = krok_count(trapezMethod, 10**-5, alfa)
plt.plot(t_res_arr1, x_res_arr1, label='TOL=10^-2', color='red')
plt.plot(t_res_arr2, x_res_arr2, label='TOL=10^-5', color='green')
plt.legend()
plt.title("Trapez x(t)")
plt.show()


plt.plot(t_res_arr1, v_res_arr1, label='TOL=10^-2', color='red')
plt.plot(t_res_arr2, v_res_arr2, label='TOL=10^-5', color='green')
plt.legend()
plt.title("Trapez v(t)")
plt.show()

plt.plot(x_res_arr1, v_res_arr1, label='TOL=10^-2', color='red')
plt.plot(x_res_arr2, v_res_arr2, label='TOL=10^-5', color='green')
plt.legend()
plt.title("Trapez V(x)")
plt.show()

plt.plot(t_res_arr1, dt_res_arr1, label='TOL=10^-2', color='red')
plt.plot(t_res_arr2, dt_res_arr2, label='TOL=10^-5', color='green')
plt.legend()
plt.title("Trapez dt(t)")
plt.show()

import matplotlib.pyplot as plt
import numpy as np
import math


def metodaJawnaEulera():

    y0 = 1
    lam = -1
    t = [0.01, 0.1, 1.0]
    estimatedy1 = []
    estimatedy2 = []
    estimatedy0 = []
    analiticvalue0 = []
    analiticvalue1 = []
    analiticvalue2 = []
    estimatedy1.append(y0)
    estimatedy2.append(y0)
    estimatedy0.append(y0)
    _i = 0.0

    while _i <= 5.0:
        analiticvalue0.append(estimatedy0[-1] - math.exp(lam * _i))
        estimatedy0.append(estimatedy0[-1] + t[0] * lam * estimatedy0[-1])
        _i = _i + t[0]
    _i = 0.0

    while _i <= 5.0:
        analiticvalue1.append(estimatedy1[-1] - math.exp(lam * _i))
        estimatedy1.append(estimatedy1[-1] + t[1] * lam * estimatedy1[-1])
        _i = _i + t[1]
    _i = 0.0

    while _i <= 5.0:
        analiticvalue2.append(estimatedy2[-1] - math.exp(lam * _i))
        estimatedy2.append(estimatedy2[-1] + t[2] * lam * estimatedy2[-1])
        _i = _i + t[2]
    _i = 0.0
# wypisywanie
    xaxis0 = np.linspace(0, 5, len(estimatedy0))
    xaxis1 = np.linspace(0, 5, len(estimatedy1))
    xaxis2 = np.linspace(0, 5, len(estimatedy2))
    analiticanalitic = []
    while _i <= 5.0:
        analiticanalitic.append(math.exp(lam * _i))
        _i += t[0]
    xaxisanalitic = np.linspace(0,5, len(analiticanalitic))
    plt.plot(xaxis0, estimatedy0, label='t=0.01')
    plt.plot(xaxis1, estimatedy1, label='t=0.1')
    plt.plot(xaxis2, estimatedy2, label='t=1.0')
    plt.plot(xaxisanalitic, analiticanalitic, label='analitic')
    plt.legend()
    plt.show()

    xaxisanal0 = np.linspace(0, 5, len(analiticvalue0))
    xaxisanal1 = np.linspace(0, 5, len(analiticvalue1))
    xaxisanal2 = np.linspace(0, 5, len(analiticvalue2))
    plt.plot(xaxisanal0, analiticvalue0, label='t=0.01')
    plt.plot(xaxisanal1, analiticvalue1, label='t=0.1')
    plt.plot(xaxisanal2, analiticvalue2, label='t=1.0')
    plt.legend()
    plt.show()


def RK2():
    y0 = 1
    lam = -1
    t = [0.01, 0.1, 1.0]
    estimatedy1 = []
    estimatedy2 = []
    estimatedy0 = []
    analiticvalue0 = []
    analiticvalue1 = []
    analiticvalue2 = []
    estimatedy1.append(y0)
    estimatedy2.append(y0)
    estimatedy0.append(y0)
    _i = 0.0

    while _i <= 5.0:
        k1 = lam*estimatedy0[-1]
        k2 = lam*(estimatedy0[-1] + k1*t[0])
        analiticvalue0.append(estimatedy0[-1] - math.exp(lam * _i))
        estimatedy0.append(estimatedy0[-1] + t[0]/2 * (k1+k2))
        _i = _i + t[0]
    _i = 0.0

    while _i <= 5.0:
        k1 = lam * estimatedy1[-1]
        k2 = lam * (estimatedy1[-1] + k1 * t[1])
        analiticvalue1.append(estimatedy1[-1] - math.exp(lam * _i))
        estimatedy1.append(estimatedy1[-1] + t[1] / 2 * (k1 + k2))
        _i = _i + t[1]
    _i = 0.0

    while _i <= 5.0:
        k1 = lam * estimatedy2[-1]
        k2 = lam * (estimatedy2[-1] + k1 * t[2])
        analiticvalue2.append(estimatedy2[-1] - math.exp(lam * _i))
        estimatedy2.append(estimatedy2[-1] + t[2]/ 2 * (k1 + k2))
        _i = _i + t[2]

    # wypisywanie
    xaxis0 = np.linspace(0, 5, len(estimatedy0))
    xaxis1 = np.linspace(0, 5, len(estimatedy1))
    xaxis2 = np.linspace(0, 5, len(estimatedy2))
    _i = 0.0
    analiticanalitic = []
    while _i <= 5.0:
        analiticanalitic.append(math.exp(lam * _i))
        _i += t[0]
    xaxisanalitic = np.linspace(0,5, len(analiticanalitic))
    plt.plot(xaxis0, estimatedy0, label='t=0.01')
    plt.plot(xaxis1, estimatedy1, label='t=0.1')
    plt.plot(xaxis2, estimatedy2, label='t=1.0')
    plt.plot(xaxisanalitic, analiticanalitic, label='analitic')
    plt.legend()
    plt.show()
    xaxisanal0 = np.linspace(0, 5, len(analiticvalue0))
    xaxisanal1 = np.linspace(0, 5, len(analiticvalue1))
    xaxisanal2 = np.linspace(0, 5, len(analiticvalue2))
    plt.plot(xaxisanal0, analiticvalue0, label='t=0.01')
    plt.plot(xaxisanal1, analiticvalue1, label='t=0.1')
    plt.plot(xaxisanal2, analiticvalue2, label='t=1.0')
    plt.legend()
    plt.show()


def RK4():
    y0 = 1
    lam = -1
    t = [0.01, 0.1, 1.0]
    estimatedy1 = []
    estimatedy2 = []
    estimatedy0 = []
    analiticvalue0 = []
    analiticvalue1 = []
    analiticvalue2 = []
    estimatedy1.append(y0)
    estimatedy2.append(y0)
    estimatedy0.append(y0)
    _i = 0.0

    while _i <= 5.0:
        k1 = lam*estimatedy0[-1]
        k2 = lam*(estimatedy0[-1] + k1 * t[0]/2.0)
        k3 = lam * (estimatedy0[-1] + k2 * t[0]/2.0)
        k4 = lam * (estimatedy0[-1] + k3 * t[0])
        analiticvalue0.append(estimatedy0[-1] - math.exp(lam * _i))
        estimatedy0.append(estimatedy0[-1] + t[0] / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4))
        _i = _i + t[0]
    _i = 0.0

    while _i <= 5.0:
        k1 = lam * estimatedy1[-1]
        k2 = lam * (estimatedy1[-1] + k1 * t[1]/2.0)
        k3 = lam * (estimatedy1[-1] + k2 * t[1] / 2.0)
        k4 = lam * (estimatedy1[-1] + k3 * t[1])
        analiticvalue1.append(estimatedy1[-1] - math.exp(lam * _i))
        estimatedy1.append(estimatedy1[-1] + t[1] / 6.0 * (k1 + 2*k2 + 2*k3 + k4))
        _i = _i + t[1]
    _i = 0.0

    while _i <= 5.0:
        k1 = lam * estimatedy2[-1]
        k2 = lam * (estimatedy2[-1] + k1 * t[2] / 2.0)
        k3 = lam * (estimatedy2[-1] + k2 * t[2] / 2.0)
        k4 = lam * (estimatedy2[-1] + k3 * t[2])
        analiticvalue2.append(estimatedy2[-1] - math.exp(lam * _i))
        estimatedy2.append(estimatedy2[-1] + t[2] / 6.0 * (k1 + 2*k2 + 2*k3 + k4))
        _i = _i + t[2]


    # wypisywanie
    xaxis0 = np.linspace(0, 5, len(estimatedy0))
    xaxis1 = np.linspace(0, 5, len(estimatedy1))
    xaxis2 = np.linspace(0, 5, len(estimatedy2))
    xaxisanal0 = np.linspace(0, 5, len(analiticvalue0))
    xaxisanal1 = np.linspace(0, 5, len(analiticvalue1))
    xaxisanal2 = np.linspace(0, 5, len(analiticvalue2))
    _i = 0.0
    analiticanalitic = []
    while _i <= 5.0:
        analiticanalitic.append(math.exp(lam * _i))
        _i += t[0]
    xaxisanalitic = np.linspace(0, 5, len(analiticanalitic))
    plt.plot(xaxis0, estimatedy0, label='t=0.01')
    plt.plot(xaxis1, estimatedy1, label='t=0.1')
    plt.plot(xaxis2, estimatedy2, label='t=1.0')
    plt.plot(xaxisanalitic, analiticanalitic, label='analitic')
    plt.legend()
    plt.show()

    plt.plot(xaxisanal0, analiticvalue0, label='t=0.01')
    plt.plot(xaxisanal1, analiticvalue1, label='t=0.1')
    plt.plot(xaxisanal2, analiticvalue2, label='t=1.0')
    plt.legend()
    plt.show()


def RK4v2(wv):
    deltat = 10 ** -4
    R = 100
    L = 0.1
    C = 0.001
    w0 = 1 / (math.sqrt(L * C))
    T0 = 2 * math.pi / w0
    tmax = 4 * T0
    Q0 = 0
    I0 = 0


    def Vfunc(time):
        return 10 * math.sin(wv*w0 * time)


    def gfunc(time, Q, I):
        return Vfunc(time)/L - (R/L)*I - (1/(L*C))*Q


    def ffunc(I):
        return I


    estimatedV = []
    estimatedQ = []
    estimatedI = []
    analiticvalueV = []
    estimatedI.append(I0)
    estimatedQ.append(Q0)
    t = 0.0

    while t <= tmax:
        k1Q = estimatedI[-1]
        k1I = gfunc(t, estimatedQ[-1], estimatedI[-1])
        k2Q = ffunc(estimatedI[-1] + (deltat/2)*k1I)
        k2I = gfunc(t+deltat/2, estimatedQ[-1] + (deltat/2)*k1Q, estimatedI[-1] + (deltat/2)*k1I)
        k3Q = ffunc(estimatedI[-1] + (deltat/2)*k2I)
        k3I = gfunc(t + deltat/2, estimatedQ[-1] + (deltat/2)*k2Q, estimatedI[-1]+ (deltat/2)*k1I)
        k4Q = ffunc(estimatedI[-1] + deltat*k3I)
        k4I = gfunc(t+deltat, estimatedQ[-1]+deltat*k3Q, estimatedI[-1] + deltat*k3I)
        estimatedI.append(estimatedI[-1] + (deltat/6)*(k1I + 2 * k2I + 2 * k3I + k4I))
        estimatedQ.append(estimatedQ[-1] + (deltat/6)*(k1Q + 2 * k2Q + 2 * k3Q + k4Q))
        t = t + deltat

    xaxis = np.linspace(0, tmax, len(estimatedQ))
    plt.plot(xaxis, estimatedQ, label='wv={}'.format(wv))


#zad1

metodaJawnaEulera()

#zad2

RK2()

#zad3

RK4()

#zad4

RK4v2(0.5)
RK4v2(0.8)
RK4v2(1.0)
RK4v2(1.2)
plt.legend()
plt.show()

import matplotlib.pyplot as plt
from math import sin
from math import pi
from math import fabs
from numpy import linspace

delta = 0.2
nx = 128
ny = 128
xmax = delta * nx
ymax = delta * ny
TOL = 10**-8


def stopintegral(v, nx, ny, delta, k):
    i = 0
    result = 0.0
    while i <= nx - k:
        j = 0
        while j <= ny - k:
            result += (k*delta)**2 / 2 * (((v[i+k][j] - v[i][j])/(2*k*delta) + ((v[i+k][j+k] - v[i][j+k])/(2*k*delta)))**2 + ((v[i][j+k] - v[i][j])/(2*k*delta) + (v[i+k][j+k] - v[i+k][j])/(2*k*delta))**2)
            j += k
        i += k
    return result


def preperemap(nx, ny, xmax, ymax, delta):
    v = [[0.0 for j in range(ny + 1)] for i in range(nx + 1)]
    for j in range(ny+1):
        v[0][j] = sin(pi*delta*j/ymax)
        v[nx][j] = sin(pi*delta*j/ymax)
    for i in range(nx+1):
        v[i][ny] = -1*sin(2*pi*delta*i/xmax)
        v[i][0] = sin(2*pi*delta*i/xmax)
    return v


def poissonsol(nx, ny, TOL, xmax, ymax, delta):
    v = preperemap(nx, ny, xmax, ymax, delta)
    k = 16
    s = 0.0
    sarr16 = []
    sarr8 = []
    sarr4 = []
    sarr2 = []
    sarr1 = []
    iter1 = 0
    iter2 = 0
    iter4 = 0
    iter8 = 0
    iter16 = 0
    iterall = 0
    

    varr16 = [[0.0 for j in range(int(ny/16)+1)] for i in range(int(nx/16)+1)]
    varr8 = [[0.0 for j in range(int(ny/8)+1)] for i in range(int(nx/8)+1)]
    varr4 = [[0.0 for j in range(int(ny/4)+1)] for i in range(int(nx/4)+1)]
    varr2 = [[0.0 for j in range(int(ny/2)+1)] for i in range(int(nx/2)+1)]
    varr1 = [[0.0 for j in range(int(ny/1)+1)] for i in range(int(nx/1)+1)]
    while k > 0:
        s1 = stopintegral(v, nx, ny, delta, k)
        while True:
            i = k
            while i <= nx-k:
                j = k
                while j <= ny-k:
                    v[i][j] = 0.25 * (v[i+k][j] + v[i-k][j] + v[i][j+k] + v[i][j-k])
                    j += k
                i += k
            s = s1
            s1 = stopintegral(v, nx, ny, delta, k)
            if k == 16:
                sarr16.append(s)
                iter16 = iterall
            if k == 8:
                sarr8.append(s)
                iter8 = iterall
            if k == 4:
                sarr4.append(s)
                iter4 = iterall 
            if k == 2:
                sarr2.append(s)
                iter2 = iterall 
            if k == 1:
                sarr1.append(s)
                iter1 = iterall 
            if fabs((s1-s)/s) <= TOL:
                break
            iterall += 1
        i = 0
        n = 0
        m = 0
        while i <= nx:
            j = 0
            while j <= ny:
                if v[i][j] != 0:
                    if k == 16:
                        n = int(i/16)
                        m = int(j/16)
                        varr16[n][m] = v[i][j]

                    if k == 8:
                        n = int(i / 8)
                        m = int(j / 8)
                        varr8[n][m] = v[i][j]

                    if k == 4:
                        n = int(i / 4)
                        m = int(j / 4)
                        varr4[n][m] = v[i][j]

                    if k == 2:
                        n = int(i / 2)
                        m = int(j / 2)
                        varr2[n][m] = v[i][j]

                    if k == 1:
                        n = int(i / 1)
                        m = int(j / 1)
                        varr1[n][m] = v[i][j]

                j += k
            i += k
        i = 0
        while i <= nx-k:
            j = 0
            while j <= ny-k:
                v[i + int(k/2)][j + int(k/2)] = 0.25*(v[i][j] + v[i+k][j] + v[i][j+k] + v[i+k][j+k])
                if i != nx-k:
                    v[i+k][j + int(k/2)] = 0.5*(v[i+k][j] + v[i+k][j+k])
                if j != ny-k:
                    v[i + int(k/2)][j+k] = 0.5 * (v[i][j+k] + v[i+k][j+k])
                if j != 0:
                    v[i + int(k/2)][j] = 0.5 * (v[i][j] + v[i+k][j])
                if i != 0:
                    v[i][j + int(k/2)] = 0.5 * (v[i][j] + v[i][j+k])
                j += k
            i += k
        k = int(k/2)
    return sarr1, sarr2, sarr4, sarr8, sarr16, varr1, varr2, varr4, varr8, varr16, iter1, iter2, iter4, iter8, iter16


def swapaxis(v):
    vnew = [[0.0 for j in range(len(v))] for i in range(len(v[0]))]
    for i in range(len(v)):
        for j in range(len(v[0])):
            vnew[j][i] = v[i][j]
    return vnew

 
sarr1, sarr2, sarr4, sarr8, sarr16, varr1, varr2, varr4, varr8, varr16, iter1, iter2, iter4, iter8, iter16 = poissonsol(nx, ny, TOL, xmax, ymax, delta)
varr1 = swapaxis(varr1) 
varr2 = swapaxis(varr2) 
varr4 = swapaxis(varr4) 
varr8 = swapaxis(varr8) 
varr16 = swapaxis(varr16) 
s16 = linspace(0, iter16 , iter16 +1 )
s8 = linspace(iter16 + 1, iter8  , iter8 + 1 - (iter16))
s4 = linspace(iter8 , iter4 , iter4 + 1 - (iter8))
s2 = linspace(iter4, iter2 , iter2 +1 -(iter4))
s1 = linspace(iter2, iter1 , iter1 + 1 - (iter2))


plt.pcolor(varr1)
plt.xlabel("x")
plt.ylabel("y")
plt.title("k = 1")
plt.colorbar()
plt.show()
plt.pcolor(varr2)
plt.xlabel("x")
plt.ylabel("y")
plt.title("k = 2")
plt.colorbar()
plt.show()
plt.pcolor(varr4)
plt.xlabel("x")
plt.ylabel("y")
plt.title("k = 4")
plt.colorbar()
plt.show()
plt.pcolor(varr8)
plt.xlabel("x")
plt.ylabel("y")
plt.title("k = 8")
plt.colorbar()
plt.show()
plt.pcolor(varr16)
plt.xlabel("x")
plt.ylabel("y")
plt.colorbar()
plt.title("k = 16")
plt.show()
plt.plot(s1+1, sarr1, label='k = 1')
plt.plot(s2+1 , sarr2, label='k = 2')
plt.plot(s4+1, sarr4, label='k = 4')
plt.plot(s8+1, sarr8, label='k = 8')
plt.plot(s16+1, sarr16, label='k = 16')
plt.xscale('log')
plt.legend()
plt.title("S(it)")
plt.show()





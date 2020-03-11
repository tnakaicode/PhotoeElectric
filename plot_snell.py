import numpy as np
import matplotlib.pyplot as plt
import numdifftools as nd
import sys
import os
import time
from scipy.integrate import odeint

sys.path.append(os.path.join('./'))
from base import plot2d

obj = plot2d()


def refract_1d(x=0, cff=2, rng=0.8):
    val = cff - (cff - 1) / (1 + np.exp(-(x - cff**2) / rng))
    return val


def n(r):
    # Refractive index function
    return refract_1d(r[0])


def diff_y(y, t):
    # Compute the differential
    grd = nd.Gradient(n)([y[0], y[1]])
    n_t = n([y[0], y[1]])
    return [y[2], y[3], grd[0] * n_t, grd[1] * n_t]


# initial position, angle
r_0 = [0, 1]
theta_0 = np.pi / 6

v_0 = [refract_1d(r_0[0]) * np.cos(theta_0),
       refract_1d(r_0[0]) * np.sin(theta_0)]

# Integration
t_range = np.arange(0, 20, 0.01)
sol = odeint(diff_y, r_0 + v_0, t_range)

X, Y = np.mgrid[0:20:1000j, 0:20:1000j]
obj.axs.plot(sol[:, 0], sol[:, 1], 'y', linewidth=2)
#pcm = obj.axs.pcolormesh(X, Y, n([X, Y]), cmap='jet', vmin=1)
pcm = obj.axs.contourf(X, Y, n([X, Y]), cmap='jet')
#bar = plt.colorbar(pcm)
#bar.ax.set_ylabel("n  Refractive index")

obj.axs.set_xlabel("cm")
obj.axs.set_ylabel("cm")
obj.SavePng("img/" + obj.rootname + ".png")

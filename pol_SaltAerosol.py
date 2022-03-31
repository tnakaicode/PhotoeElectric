import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps
import numdifftools as nd
import sys
import os
import time
# import a single function for integration using trapezoidal rule
from scipy.integrate import trapz
import argparse

sys.path.append(os.path.join("./"))
from src.base import plot2d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--pxyz", dest="point",
                      default=[0.0, 0.0, 0.0], type=float, nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    m = 1.536  # refractive index of NaCl
    wavelength = 405  # replace with the laser wavelength (nm)

    dp_g = 85  # geometric mean diameter - replace with your own (nm)
    # geometric standard deviation - replace with your own (unitless)
    sigma_g = 1.5
    N = 1e5  # total number of particles - replace with your own (cm^-3)

    # Calculate optical properties
    B = ps.Mie_Lognormal(m, wavelength, sigma_g, dp_g,
                         N, returnDistribution=True)
    S = ps.SF_SD(m, wavelength, B[7], B[8])

    # Make graphs - lots of this is really unnecessary decoration for a pretty graph.

    obj = plot2d(aspect="auto")
    obj.axs.plot(S[0], S[1], 'b', ls='dashdot',
                 lw=1, label="Parallel Polarization")
    obj.axs.plot(S[0], S[2], 'r', ls='dashed', lw=1,
                 label="Perpendicular Polarization")
    obj.axs.plot(S[0], S[3], 'k', lw=1, label="Unpolarized")

    x_label = ["0", r"$\mathregular{\frac{\pi}{4}}$", r"$\mathregular{\frac{\pi}{2}}$",
               r"$\mathregular{\frac{3\pi}{4}}$", r"$\mathregular{\pi}$"]
    x_tick = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
    obj.axs.set_xticks(x_tick)
    obj.axs.set_xticklabels(x_label)
    obj.axs.tick_params(which='both', direction='in')
    obj.axs.set_xlabel("Scattering Angle ϴ")
    obj.axs.set_ylabel(r"Intensity ($\mathregular{|S|^2}$)", labelpad=10)
    handles, labels = obj.axs.get_legend_handles_labels()
    obj.fig.legend(handles, labels, ncol=3, loc=8)
    obj.fig.suptitle("Scattering Intensity Functions")

    # Highlight certain angles and compute integral
    sixty = [0.96 < x < 1.13 for x in S[0]]
    ninety = [1.48 < x < 1.67 for x in S[0]]
    obj.axs.fill_between(S[0], 0, S[3], where=sixty, color='g', alpha=0.15)
    obj.axs.fill_between(S[0], 0, S[3], where=ninety, color='g', alpha=0.15)
    obj.axs.set_yscale('log')

    int_sixty = trapz(S[3][110:130], S[0][110:130])
    int_ninety = trapz(S[3][169:191], S[0][169:191])

    # Annotate plot with integral results
    obj.axs.annotate("Integrated value = {i:1.3f}".format(i=int_sixty),
                     xy=(np.pi / 3, S[3][120]), xycoords='data',
                     xytext=(np.pi / 6, 0.8), textcoords='data',
                     arrowprops=dict(arrowstyle="->",
                                     connectionstyle="arc3"),
                     )
    obj.axs.annotate("Integrated value = {i:1.3f}".format(i=int_ninety),
                     xy=(np.pi / 2, S[3][180]), xycoords='data',
                     xytext=(2 * np.pi / 5, 2), textcoords='data',
                     arrowprops=dict(arrowstyle="->",
                                     connectionstyle="arc3"),
                     )
    obj.SavePng()

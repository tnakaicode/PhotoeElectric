import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps
import sys
import os
import time
import glob
import imageio
import scipy.constants as cnt
import argparse
from scipy.integrate import trapz

"""
Recently, 
a colleague needed to know how much light a distribution of salt aerosol would scatter into two detectors, 
one at 60° and one at 90°. 
We modeled a lognormal distribution of NaCl particles based on laboratory measurements 
and then tried to figure out how much light we'd see at various angles.
"""

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default="./")
    parser.add_argument("--file", dest="file", default="plot_prop_mat.txt")
    parser.add_argument("--freq", dest="freq", default=100.0, type="float")
    parser.add_argument("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt = parser.parse_args()
    print(opt, argvs)

    # refractive index of NaCl
    m = 1.536

    # replace with the laser wavelength (nm)
    wavelength = 405

    # geometric mean diameter - replace with your own (nm)
    dp_g = 85

    # geometric standard deviation - replace with your own (unitless)
    sigma_g = 1.5

    # total number of particles - replace with your own (cm^-3)
    N = 1e5

    # Calculate optical properties
    B = ps.Mie_Lognormal(m, wavelength, sigma_g, dp_g, N,
                         returnDistribution=True)
    print(B)
    S = ps.SF_SD(m, wavelength, B[7], B[8])

    # Make graphs - lots of this is really unnecessary decoration for a pretty graph.
    plt.close('all')

    fig1 = plt.figure(figsize=(10.08, 6.08))
    ax1 = fig1.add_subplot(1, 1, 1)

    ax1.plot(S[0], S[1], 'b', ls='dashdot',
             lw=1, label="Parallel Polarization")
    ax1.plot(S[0], S[2], 'r', ls='dashed', lw=1,
             label="Perpendicular Polarization")
    ax1.plot(S[0], S[3], 'k', lw=1, label="Unpolarized")

    x_label = ["0", r"$\mathregular{\frac{\pi}{4}}$", r"$\mathregular{\frac{\pi}{2}}$",
               r"$\mathregular{\frac{3\pi}{4}}$", r"$\mathregular{\pi}$"]
    x_tick = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
    ax1.set_xticks(x_tick)
    ax1.set_xticklabels(x_label, fontsize=14)
    ax1.tick_params(which='both', direction='in')
    ax1.set_xlabel("Scattering Angle ϴ", fontsize=16)
    ax1.set_ylabel(
        r"Intensity ($\mathregular{|S|^2}$)", fontsize=16, labelpad=10)
    handles, labels = ax1.get_legend_handles_labels()
    fig1.legend(handles, labels, fontsize=14, ncol=3, loc=8)

    fig1.suptitle("Scattering Intensity Functions", fontsize=18)
    fig1.show()
    plt.tight_layout(rect=[0.01, 0.05, 0.915, 0.95])

    # Highlight certain angles and compute integral
    sixty = [0.96 < x < 1.13 for x in S[0]]
    ninety = [1.48 < x < 1.67 for x in S[0]]
    ax1.fill_between(S[0], 0, S[3], where=sixty, color='g', alpha=0.15)
    ax1.fill_between(S[0], 0, S[3], where=ninety, color='g', alpha=0.15)
    ax1.set_yscale('log')

    int_sixty = trapz(S[3][110:130], S[0][110:130])
    int_ninety = trapz(S[3][169:191], S[0][169:191])

    # Annotate plot with integral results
    ax1.annotate("Integrated value = {i:1.3f}".format(i=int_sixty),
                 xy=(np.pi / 3, S[3][120]), xycoords='data',
                 xytext=(np.pi / 6, 0.8), textcoords='data',
                 arrowprops=dict(arrowstyle="->",
                                 connectionstyle="arc3"),
                 )
    ax1.annotate("Integrated value = {i:1.3f}".format(i=int_ninety),
                 xy=(np.pi / 2, S[3][180]), xycoords='data',
                 xytext=(2 * np.pi / 5, 2), textcoords='data',
                 arrowprops=dict(arrowstyle="->",
                                 connectionstyle="arc3"),
                 )
    plt.savefig("./mie_aerosol.png")

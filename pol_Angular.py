import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps
import numdifftools as nd
import sys
import os
import time
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from src.base import plot2d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="point",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    m = 1.7 + 0.5j
    w = 532
    d = 5000

    theta, SL, SR, SU = ps.ScatteringFunction(m, w, d)
    qR, SLQ, SRQ, SUQ = ps.ScatteringFunction(m, w, d, space='qspace')

    obj = plot2d()

    ax1 = obj.add_axs(1, 2, 1)
    ax2 = obj.add_axs(1, 2, 2)

    ax1.semilogy(theta, SL, 'b', ls='dashdot', lw=1,
                 label="Parallel Polarization")
    ax1.semilogy(theta, SR, 'r', ls='dashed', lw=1,
                 label="Perpendicular Polarization")
    ax1.semilogy(theta, SU, 'k', lw=1, label="Unpolarized")

    x_label = ["0", r"$\mathregular{\frac{\pi}{4}}$", r"$\mathregular{\frac{\pi}{2}}$",
               r"$\mathregular{\frac{3\pi}{4}}$", r"$\mathregular{\pi}$"]
    x_tick = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
    ax1.set_xticks(x_tick)
    ax1.set_xticklabels(x_label, fontsize=14)
    ax1.tick_params(which='both', direction='in')
    ax1.set_xlabel("ϴ", fontsize=16)
    ax1.set_ylabel(
        r"Intensity ($\mathregular{|S|^2}$)", fontsize=16, labelpad=10)

    ax2.loglog(qR, SLQ, 'b', ls='dashdot', lw=1, label="Parallel Polarization")
    ax2.loglog(qR, SRQ, 'r', ls='dashed', lw=1,
               label="Perpendicular Polarization")
    ax2.loglog(qR, SUQ, 'k', lw=1, label="Unpolarized")

    ax2.tick_params(which='both', direction='in')
    ax2.set_xlabel("qR", fontsize=14)
    handles, labels = ax1.get_legend_handles_labels()

    obj.fig.legend(handles, labels, fontsize=14, ncol=3, loc=8)

    obj.fig.suptitle("Scattering Intensity Functions", fontsize=18)
    plt.tight_layout(rect=[0.01, 0.05, 0.915, 0.95])
    obj.SavePng()

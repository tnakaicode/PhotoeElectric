import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps
import numdifftools as nd
import sys
import os
import time
import glob
import imageio
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from base import plot2d

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

    obj = plot2d()
    obj.create_tempdir(flag=-1)

    wavelength = 450.0
    m = 1.536 + 0.0015j
    drange = np.logspace(1, np.log10(500 * 405 / np.pi), 250)
    for i, d in enumerate(drange):
        if 250 % (i + 1) == 0:
            print("Working on image " + str(i) + "...", flush=True)
        theta, SL, SR, SU = ps.ScatteringFunction(
            m, wavelength, d, space='theta')

        #ax2 = fig1.add_subplot(1,2,2)
        obj.new_2Dfig(aspect="auto")
        obj.axs.semilogy(theta, SL, 'b', ls='dashdot', lw=1,
                         label="Parallel Polarization")
        obj.axs.semilogy(theta, SR, 'r', ls='dashed', lw=1,
                         label="Perpendicular Polarization")
        obj.axs.semilogy(theta, SU, 'k', lw=1, label="Unpolarized")

        x_label = ["0", r"$\mathregular{\frac{\pi}{4}}$", r"$\mathregular{\frac{\pi}{2}}$",
                   r"$\mathregular{\frac{3\pi}{4}}$", r"$\mathregular{\pi}$"]
        x_tick = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
        obj.axs.set_xticks(x_tick)
        obj.axs.set_xticklabels(x_label, fontsize=14)
        obj.axs.tick_params(which='both', direction='in')
        obj.axs.set_xlabel("ϴ", fontsize=16)
        obj.axs.set_ylabel(
            r"Intensity ($\mathregular{|S|^2}$)", fontsize=16, labelpad=10)
        #obj.axs.set_ylim([1e-9, 10])
        obj.axs.set_xlim([1e-3, theta[-1]])
        obj.axs.annotate("x = πd/λ = {dd:1.2f}".format(dd=np.round(np.pi * d / 405, 2)), xy=(3, 1e-6), xycoords='data',
                         xytext=(0.05, 0.1), textcoords='axes fraction',
                         horizontalalignment='left', verticalalignment='top',
                         fontsize=18
                         )
        handles, labels = obj.axs.get_legend_handles_labels()
        obj.fig.legend(handles, labels, fontsize=14, ncol=3, loc=8)

        obj.fig.suptitle("Scattering Intensity Functions", fontsize=18)
        # plt.tight_layout(rect=[0.01,0.05,0.915,0.95])
        obj.SavePng_Serial()
        plt.close()
        #plt.savefig('output\\' + str(i).rjust(3,'0') + '.png')

    filenames = glob.glob(obj.tmpdir + '*.png')
    dur = [0.1 for x in range(250)]
    dur[249] = 10
    with imageio.get_writer(obj.tmpdir + 'mie_ripples.mp4', mode='I', fps=10) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

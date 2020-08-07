import PyMieScatt as ps
import numpy as np
import matplotlib.pyplot as plt
import imageio
import os

wavelength = 450.0
m = 1.536 + 0.0015j
drange = np.logspace(1, np.log10(500 * 405 / np.pi), 250)
for i, d in enumerate(drange):
    if 250 % (i + 1) == 0:
        print("Working on image " + str(i) + "...", flush=True)
    theta, SL, SR, SU = ps.ScatteringFunction(
        m, wavelength, d, space='theta')

    plt.close('all')

    fig1 = plt.figure(figsize=(10.08, 6.08))
    ax1 = fig1.add_subplot(1, 1, 1)
    #ax2 = fig1.add_subplot(1,2,2)

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
    ax1.set_ylim([1e-9, 1])
    ax1.set_xlim([1e-3, theta[-1]])
    ax1.annotate("x = πd/λ = {dd:1.2f}".format(dd=np.round(np.pi * d / 405, 2)), xy=(3, 1e-6), xycoords='data',
                 xytext=(0.05, 0.1), textcoords='axes fraction',
                 horizontalalignment='left', verticalalignment='top',
                 fontsize=18
                 )
    handles, labels = ax1.get_legend_handles_labels()
    fig1.legend(handles, labels, fontsize=14, ncol=3, loc=8)

    fig1.suptitle("Scattering Intensity Functions", fontsize=18)
    plt.tight_layout(rect=[0.01, 0.05, 0.915, 0.95])

    plt.savefig("./temp_output/" + str(i).rjust(3, '0') + '.png')

filenames = os.listdir('./temp_output/')
dur = [0.1 for x in range(250)]
dur[249] = 10
with imageio.get_writer('mie_ripples.mp4', mode='I', fps=10) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)

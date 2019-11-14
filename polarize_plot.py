import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import pi, sin, cos, tan, arcsin, linspace


def pol_plot(tdeg, Rp, Tp, Rs, Ts):
    plt.figure(figsize=(8, 6))
    plt.plot(tdeg, Rp, label=r"$R_{12}^{\rm{p}}$",
             linewidth=3.0, color='black', linestyle='dashed')
    plt.plot(tdeg, Tp, label=r"$T_{12}^{\rm{p}}$",
             linewidth=3.0, color='black')
    plt.plot(tdeg, Rs, label=r"$R_{12}^{\rm{s}}$",
             linewidth=3.0, color='gray', linestyle='dashed')
    plt.plot(tdeg, Ts, label=r"$T_{12}^{\rm{s}}$", linewidth=3.0, color='gray')

    plt.xlabel("Inject (deg.)", fontsize=20)
    plt.ylabel("Reflect", fontsize=20)
    plt.title("refelect and absorb", fontsize=18)

    plt.grid(True)
    plt.axis([0.0, 90, -1, 1])
    plt.legend(fontsize=20, loc='lower right')
    plt.tick_params(labelsize=20)
    plt.tight_layout()


def pol_plot_ref(tdeg, Rp, Rs):
    plt.figure(figsize=(8, 6))
    plt.plot(tdeg, Rp, label=r"$R_{12}^{\rm{p}}$",
             linewidth=3.0, color='black')
    plt.plot(tdeg, Rs, label=r"$R_{12}^{\rm{s}}$", linewidth=3.0, color='gray')

    plt.xlabel("Inject (deg.)", fontsize=20)
    plt.ylabel("Reflect", fontsize=20)
    plt.title("refelect rate", fontsize=18)

    plt.grid(True)
    #plt.axis([0.0, 90, -1, 1])
    plt.legend(fontsize=20, loc='lower right')
    plt.tick_params(labelsize=20)
    plt.tight_layout()

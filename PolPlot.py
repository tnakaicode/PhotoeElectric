from operator import imod
import numpy as np
import matplotlib.pyplot as plt

from src.base import plot2d
from src.RI import WLx, epAg, RIAg


class plot1d (plot2d):

    def __init__(self):
        plot2d.__init__(self)

    def plot_silver(self):
        self.new_2Dfig(aspect="auto")
        self.axs.axis([300, 700, -10, 5])
        self.axs.plot(WLx, np.real(epAg), label="real",
                      linewidth=3.0, color='black')
        self.axs.plot(WLx, np.imag(epAg), label="imag",
                      linewidth=3.0, color='gray')
        self.axs.set_xlabel("wave (nm)", fontsize=22)       # x-axis label
        self.axs.set_ylabel("dielective", fontsize=22)          # y-axis label
        self.axs.set_title("Dierective of SIlver", fontsize=22)
        self.axs.tick_params(labelsize=20)  # scale fontsize=18pt
        self.axs.legend(fontsize=20, loc='lower right')
        self.fig.tight_layout()
        self.SavePng()

    def plot_silver_scatter(self):
        n1 = 1
        n2 = RIAg
        r = 25
        k = 2 * np.pi / WLx
        alpha = 4 * np.pi * (r**3) * (n1**2) * \
            (n2**2 - n1**2) / (n2**2 + 2 * n1**2)

        Csca = k**4 / (6 * np.pi) * abs(alpha)**2
        Cabs = k * np.imag(alpha)
        Qsca = Csca / ((r**2) * np.pi)
        Qabs = Cabs / ((r**2) * np.pi)

        self.new_2Dfig(aspect="auto")
        self.axs.axis([300, 700, -10, 5])
        self.axs.plot(WLx, Csca / 10**4,
                      label=r"$C_{{\rm sca}}$",
                      linewidth=3.0, color='black')
        self.axs.plot(WLx, Cabs / 10**4,
                      label=r"$C_{{\rm abs}}$",
                      linewidth=3.0, color='gray')
        self.axs.set_xlabel("wave (nm)", fontsize=22)
        self.axs.set_ylabel("dielective", fontsize=22)
        self.axs.set_title("dierective of SIlver", fontsize=22)
        self.axs.tick_params(labelsize=20)
        self.axs.legend(fontsize=20, loc='lower right')
        self.fig.tight_layout()
        self.SavePng(self.tempname + "_scatter.png")


if __name__ == '__main__':
    obj = plot1d()
    obj.plot_silver()
    obj.plot_silver_scatter()

import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps
import numdifftools as nd
import sys
import os
import time
import glob
import imageio
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.ndimage import zoom
from optparse import OptionParser


sys.path.append(os.path.join("./"))
from base import plot2d

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(
        n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def ithPart(gammai, dp, dpgi, sigmagi):
    return (gammai / (np.sqrt(2 * np.pi) * np.log(sigmagi)
                      * dp)) * np.exp(-(np.log(dp) - np.log(dpgi))**2 / (2 * np.log(sigmagi)**2))


if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="point",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    N = 1e6
    w = 405
    maxDiameter = 3500
    numDiams = 1200

    dp = np.logspace(np.log10(1), np.log10(maxDiameter), numDiams)
    sigmaList = np.logspace(np.log10(1.005), np.log10(2), 49)
    mu = 300
    ndp = [N * ithPart(1, dp, mu, s) for s in sigmaList]

    deltaD = np.zeros(numDiams)
    deltaD[838] = N

    lognormalList = [deltaD] + ndp
    sigmaList = np.insert(sigmaList, 0, 1)

    # Test region - uncomment for a single graph
    #testCase = np.random.randint(1,49)
    #lognormalList = [lognormalList[testCase]]
    #sigmaList = [sigmaList[testCase]]

    BscaSolution = []
    BabsSolution = []

    for l in lognormalList:
        _, _s, _a, *rest = ps.Mie_SD(1.6 + 0.36j, w, dp, l)
        BscaSolution.append(_s)
        BabsSolution.append(_a)

    nMin = 1.3
    nMax = 3
    kMin = 0
    kMax = 2

    points = 40
    interpolationFactor = 2

    nRange = np.linspace(nMin, nMax, points)
    kRange = np.linspace(kMin, kMax, points)

    obj = plot2d()

    for i, (sigma, l, ssol, asol) in enumerate(zip(sigmaList, lognormalList, BscaSolution, BabsSolution)):
        start = time()
        BscaList = []
        BabsList = []
        nList = []
        kList = []
        for n in nRange:
            s = []
            a = []
            for k in kRange:
                m = n + k * 1.0j
                _, Bsca, Babs, *rest = ps.Mie_SD(m, w, dp, l)
                s.append(Bsca)
                a.append(Babs)
            BscaList.append(s)
            BabsList.append(a)
        n = zoom(nRange, interpolationFactor)
        k = zoom(kRange, interpolationFactor)
        BscaSurf = zoom(np.transpose(np.array(BscaList)), interpolationFactor)
        BabsSurf = zoom(np.transpose(np.array(BabsList)), interpolationFactor)
        nSurf, kSurf = np.meshgrid(n, k)

        c1 = truncate_colormap(cm.Reds, 0.2, 1, n=256)
        c2 = truncate_colormap(cm.Blues, 0.2, 1, n=256)

        xMin, xMax = nMin, nMax
        yMin, yMax = kMin, kMax

        plt.close('all')
        obj.new_2Dfig()
        plt.suptitle("σ={ww:1.3f}".format(ww=sigma), fontsize=24)

        ax1 = plt.subplot2grid(
            (3, 4), (0, 0), projection='3d', rowspan=2, colspan=2)
        ax2 = plt.subplot2grid(
            (3, 4), (0, 2), projection='3d', rowspan=2, colspan=2)
        ax3 = plt.subplot2grid((3, 4), (2, 0), colspan=3)
        ax4 = plt.subplot2grid((3, 4), (2, 3))

        ax1.plot_surface(nSurf, kSurf, BscaSurf, rstride=1,
                         cstride=1, cmap=c1, alpha=0.5)
        ax1.contour(nSurf, kSurf, BscaSurf, [
            ssol], lw=2, colors='r', linestyles='dashdot')
        ax1.contour(nSurf, kSurf, BscaSurf, [
            ssol], colors='r', linestyles='dashdot', offset=0)

        ax2.plot_surface(nSurf, kSurf, BabsSurf, rstride=1,
                         cstride=1, cmap=c2, alpha=0.5, zorder=-1)
        ax2.contour(nSurf, kSurf, BabsSurf, [
            asol], lw=2, colors='b', linestyles='solid', zorder=3)
        ax2.contour(nSurf, kSurf, BabsSurf, [
            asol], colors='b', linestyles='solid', offset=0)

        boxLabels = ["βsca", "βabs"]

        yticks = [2, 1.5, 1, 0.5, 0]
        xticks = [3, 2.5, 2, 1.5]

        for a, t in zip([ax1, ax2], boxLabels):
            lims = a.get_zlim3d()
            a.set_zlim3d(0, lims[1])
            a.text(1.5, 0, (a.get_zlim3d()[1]) * 1.15, t,
                   ha="center", va="center", size=18, zorder=5)
            a.set_ylim(2, 0)
            a.set_xlim(3, 1.3)
            a.set_xticks(xticks)
            a.set_xticklabels(xticks, rotation=-15, va='center', ha='left')
            a.set_yticks(yticks)
            a.set_yticklabels(yticks, rotation=-15, va='center', ha='left')
            a.set_zticklabels([])
            a.view_init(20, 120)
            a.tick_params(axis='both', which='major', labelsize=12, pad=0)
            a.tick_params(axis='y', pad=-2)
            a.set_xlabel("n", fontsize=18, labelpad=2)
            a.set_ylabel("k", fontsize=18, labelpad=3)

        ax3.semilogx(dp, l, c='g')
        ax3.set_xlabel('Diameter', fontsize=16)
        ax3.get_yaxis().set_ticks([])
        ax3.tick_params(which='both', direction='in')
        ax3.grid(color='#dddddd')

        giv = ps.ContourIntersection_SD(
            ssol, asol, w, dp, l, gridPoints=points * 1.5, kMin=0.001, kMax=2, axisOption=10, fig=fig1, ax=ax4)
        ax4.set_xlim(1.3, 3)
        ax4.yaxis.tick_right()
        ax4.yaxis.set_label_position("right")
        ax4.legend_.remove()
        ax4.set_title("")
        ax4.set_yscale('linear')

        plt.tight_layout()
        obj.SavePng_Serial()

        end = time()
        print(
            "Frame {n:1d}/30 done in {t:1.2f} seconds.".format(n=i + 1, t=end - start))

    #filenames = os.listdir('Distro\\')
    # with imageio.get_writer('SD.mp4', mode='I', fps=5) as writer:
    #    for filename in filenames:
    #        image = imageio.imread('Distro\\' + filename)
    #        writer.append_data(image)

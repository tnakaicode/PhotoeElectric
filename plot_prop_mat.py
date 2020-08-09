import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps
import numdifftools as nd
import sys
import os
import time
import glob
import imageio
import scipy.constants as cnt
from optparse import OptionParser

sys.path.append(os.path.join("./"))
from src.base import plot2d
from src.Unit import convert, convert_freq_to_wave, convert_wave_to_freq
from src.index import matFAI, mMATp, mMATs

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--freq", dest="freq", default=100.0, type="float")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    pt = np.linspace(-1, 1, 100) * 90
    data = np.loadtxt("plot_prop_mat.txt", comments="#")
    wave = 500.0  # [nm]
    freq = convert_wave_to_freq(wave, "nm", "THz")
    knum = 2 * np.pi / wave
    print(freq)
    print(data)
    print(data.shape)
    n0 = data[0, 0]
    matGs = np.empty((2, 2), dtype=complex)
    matGp = np.empty((2, 2), dtype=complex)
    # for t in pt:
    t = 45.0
    for idx, dat in enumerate(data[1:]):
        print(idx, *dat)
        i1, i2 = idx, idx + 1
        i, n1_real, n1_imag, d1 = data[i1]
        i, n2_real, n2_imag, d2 = data[i2]
        n1 = n1_real + 1j * n1_imag
        n2 = n2_real + 1j * n2_imag
        mPhi = matFAI(n2, d2, t, knum)
        matGs *= mPhi @ mMATs(n1, n2, t, knum)
        matGp *= mPhi @ mMATp(n1, n2, t, knum)
        rs = -1 * matGs[1, 0] / matGs[1, 1]
        rp = -1 * matGp[0, 1] / matGp[1, 1]
        ts = matGs[0, 0] - matGs[0, 1] * matGs[1, 0] / matGs[1, 1]
        tp = matGp[0, 0] - matGp[0, 1] * matGp[1, 0] / matGp[1, 1]
        print(i2, n2, d2)
        print(rs, ts)
        print(rp, tp)

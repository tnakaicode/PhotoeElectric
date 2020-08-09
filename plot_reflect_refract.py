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
from src.base import plot2d
from src.index import ref_index

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = OptionParser()
    parser.add_option("--dir", dest="dir", default="./")
    parser.add_option("--pxyz", dest="pxyz",
                      default=[0.0, 0.0, 0.0], type="float", nargs=3)
    opt, argc = parser.parse_args(argvs)
    print(opt, argc)

    obj = plot2d(aspect="auto")
    pt = np.linspace(0, 1, 100) * 90
    tp, rp, ts, rs = ref_index(n1=1.0, n2=1.0, t0=pt)
    obj.axs.plot(pt, tp, label="Transmit P-pol")
    obj.axs.plot(pt, rp, label="Reflect P-pol")
    obj.axs.plot(pt, ts, label="Transmit S-pol")
    obj.axs.plot(pt, rs, label="Reflect S-pol")
    obj.axs.legend()
    obj.SavePng()

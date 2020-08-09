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
from src.index import matFAI, mMATp, mMATs

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

    pt = np.linspace(-1, 1, 100) * 90
    data = np.loadtxt("plot_prop_mat.txt", comments="#")
    n0 = data[0, 0]
    for idx, dat in enumerate(data[1:]):
        i, j = idx, idx + 1
        print(i, *dat)

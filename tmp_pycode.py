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

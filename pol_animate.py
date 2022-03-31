import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt as ps
import numdifftools as nd
import sys
import os
import time
import glob
import imageio
import argparse

sys.path.append(os.path.join("./"))
from base import plot2d, create_tempdir

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('PIL').setLevel(logging.ERROR)

if __name__ == '__main__':
    argvs = sys.argv
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", dest="dir", default=None)
    opt = parser.parse_args()
    print(opt, argvs)

    if opt.dir == None:
        tmpdir = create_tempdir()
    else:
        tmpdir = opt.dir

    filenames = glob.glob(tmpdir + '*.png')
    dur = [0.1 for x in range(250)]
    dur[249] = 10
    with imageio.get_writer(tmpdir + 'mie_ripples.mp4', mode='I', fps=10) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

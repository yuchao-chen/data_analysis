# -*- coding: utf-8 -*-
import os
import shutil

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation as animation

import astropy.io.fits as fits
from astropy.visualization import ZScaleInterval

def Demod():
    src_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180705.nosync/SYNC_RECON_B1_BY_PHASE_AND_MODE_SEPARATELY'
    result_path = '/tmp/0/'
    if os.path.exists(result_path):
        shutil.rmtree(result_path)
    os.mkdir(result_path)

    files = [os.path.join(src_path, f) for f in os.listdir(src_path) if os.path.isfile(os.path.join(src_path, f)) and f.endswith('.fits')]
    files.sort()

    interval = ZScaleInterval()

    for j in range(0, len(files)//4):
        subfiles = files[j*4:j*4+4]
        data = []
        for f in subfiles:
            print(f)
            tmp = fits.getdata(f)
            data.append(tmp)
        print(len(data))
        i = 0.419 * data[0] + 0.581 * data[1] + 0.581 * data[2] + 0.418 * data[3]
        q = 0.889 * data[0] - 0.887 * data[1] - 0.888 * data[2] + 0.886 * data[3]
        u = 1.131 * data[0] - 0.547 * data[1] + 0.547 * data[2] - 1.133 * data[3]
        v = -0.687 * data[0] - 0.995 * data[1] + 0.994 * data[2] + 0.688 * data[3]

        plt.figure()
        plt.imshow(interval(i), cmap = 'gray', origin = 'lower')
        plt.colorbar()
        plt.figure()
        plt.imshow(interval(q/i), cmap = 'gray', origin = 'lower')
        plt.figure()
        plt.imshow(interval(u/i), cmap = 'gray', origin = 'lower')
        plt.figure()
        plt.imshow(interval(v/i), cmap = 'gray', origin = 'lower')
        plt.show()
        return

def main():
    Demod()

if __name__ == '__main__':
    main()

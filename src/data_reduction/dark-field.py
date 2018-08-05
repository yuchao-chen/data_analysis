# -*- coding: utf-8 -*-
import os

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import astropy.io.fits as fits
from astropy.visualization import ZScaleInterval

def difference_of_two_dark_fields():
    path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_PINHOLES&FOCUS/DARK.nosync/'

    d0 = fits.getdata(path + '200frames@201807018_10SECONDS_LATER.fits')
    d1 = fits.getdata(path + '200frames@201807018.fits')

    dif = d0 - d1
    std = np.std(dif)
    mean = np.mean(dif)

    fig, ax = plt.subplots()

    interval = ZScaleInterval()

    im = ax.imshow(interval(dif), cmap = 'gray')
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    ax.text(100, 1800, 'std: ' + "{:.5f}".format(std), fontsize = 10)
    ax.text(100, 1700, 'mean: ' + "{:.5f}".format(mean), fontsize = 10)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im, cax = cax)

    fig.savefig('/tmp/dark.eps')
    plt.show()

def main():
    difference_of_two_dark_fields()

if __name__ == '__main__':
    main()

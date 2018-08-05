# -*- coding: utf-8 -*-
import os

import matplotlib.pyplot as plt
import numpy as np

import astropy.io.fits as fits
from astropy.visualization import ZScaleInterval

def focus_detect(path):
    interval = ZScaleInterval()

    files = [os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.endswith('.fits')]
    files.sort()

    avg_energy_spectrum = 0
    energy_spectrums = []

    for f in files:
        tmp = fits.getdata(f)
        tmp = tmp[500:1500, 500:1500]
        tmp = tmp.astype(float)
        freq = np.fft.fftshift(np.fft.fftn(tmp))
        energy_spectrum = freq * freq.conj()
        energy_spectrums.append(energy_spectrum)
        avg_energy_spectrum += energy_spectrum / float(len(files))

    m = []
    for e in energy_spectrums:
        e /= avg_energy_spectrum
        m.append(np.mean(np.real(e)))
        #plt.imshow(interval(np.real(e)))
        #plt.show()

    #plt.imshow(interval(np.real(avg_energy_spectrum)))
    plt.figure()
    plt.plot(m)
    plt.show()

def main():
    path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_FLAT_PINHOLES&FOCUS/DEFOCUS.nosync/20180722/022624/CENT'
    focus_detect(path)

if __name__ == '__main__':
    main()

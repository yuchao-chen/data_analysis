from os import listdir
from os.path import isfile, join

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits

def FeIProfileAt5324A():
    # FTS profile
    fts_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/5TH_EXPERIMENT/20180118/solar_spectrum_at_5324A.txt'
    fts = np.genfromtxt(fts_path, delimiter=' ')
    center = 5324.184
    x0 = fts[:, 0:1]
    y0 = fts[:, 1:2]
    x0 = x0[850:1300]
    y0 = y0[850:1300]
    plt.plot(x0, y0)
    plt.axvline(x = center)

    # Measured profile by 5324A Lyot Filter @20180509
    measured_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180509/LYOT_FILTER_PROFILE.nosync'
    measured_files = [join(measured_path, f) for f in listdir(measured_path) if isfile(join(measured_path, f))]
    measured_files.sort()
    measured_profile_0 = []
    measured_profile_1 = []
    for f in measured_files:
        tmp = fits.getdata(f)
        tmp0 = tmp[1024:1500, 1024:1500]
        measured_profile_0.append(np.mean(tmp0))
        tmp1 = tmp[0:500, 1024:1500]
        measured_profile_1.append(np.mean(tmp1))

    x1 = center + 0.4 - np.arange(len(measured_profile_0)) * 0.01 - 0.07
    y1 = np.array(measured_profile_0)
    pos = np.where(x0 == x1[len(x1)-1])
    ratio = y0[pos]/y1[len(y1)-1]
    plt.plot(x1, y1*ratio)

    x1 = center + 0.4 - np.arange(len(measured_profile_1)) * 0.01 - 0.07
    y1 = np.array(measured_profile_1)
    pos = np.where(x0 == x1[len(x1)-1])
    ratio = y0[pos]/y1[len(y1)-1]
    plt.plot(x1, y1*ratio)

    # Measured profile by 5324A Lyot Filter @20180118
    data = np.genfromtxt('/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/5TH_EXPERIMENT/20180118/profile0.txt')
    x2 = data[:, 0:1]
    x2 = x2 / 1000.0 + center
    y2 = data[:, 1:2]
    ratio = y0[pos] / y2[7]
    y2 = y2 * ratio
    plt.plot(x2, y2, 'r.')
    plt.show()

if __name__ == '__main__':
    FeIProfileAt5324A()


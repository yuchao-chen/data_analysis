from os import listdir
from os.path import isfile, join

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits

def FeIProfileAt5324A():
    plt.figure(figsize = (8, 7))
    # FTS profile
    fts_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/5TH_EXPERIMENT/20180118/solar_spectrum_at_5324A.txt'
    fts = np.genfromtxt(fts_path, delimiter=' ')
    center = 5324.184
    x0 = fts[:, 0:1]
    y0 = fts[:, 1:2]
    x0 = x0[820:1330]
    y0 = y0[820:1330]
    p0 = plt.plot(x0, y0, label = 'FTS profile')
    plt.axvline(x = center)

    # measured profile by 5324A Lyot Filter @20180509
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
    p1 = plt.plot(x1, y1*ratio, '.', label = '@20180509')

    '''
    x1 = center + 0.4 - np.arange(len(measured_profile_1)) * 0.01 - 0.07
    y1 = np.array(measured_profile_1)
    pos = np.where(x0 == x1[len(x1)-1])
    ratio = y0[pos]/y1[len(y1)-1]
    plt.plot(x1, y1*ratio)
    '''

    # measured profile by 5324A Lyot Filter @20180118
    data = np.genfromtxt('/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/5TH_EXPERIMENT/20180118/profile0.txt')
    x2 = data[:, 0:1]
    x2 = x2 / 1000.0 + center
    y2 = data[:, 1:2]
    ratio = y0[pos] / y2[7]
    y2 = y2 * ratio
    p2 = plt.plot(x2, y2, '.', label = '@20180118')

    # scan from 0.4 to -0.4 @20180719
    measured_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_FLAT_PINHOLES&FOCUS/LYOT_FILTER_PROFILE.nosync/20180718/033202'

    measured_files = [join(measured_path, f) for f in listdir(measured_path) if isfile(join(measured_path, f)) and f.endswith('.fits')]
    measured_files.sort()
    x3 = []
    y3 = []
    for f in measured_files:
        header = fits.getheader(f)
        x3.append(float(header['OFFBAND']) + center - 0.09)

        tmp = fits.getdata(f)
        shape = tmp.shape
        tmp = tmp[shape[0]//4:shape[0]*3//4, shape[1]//4:shape[1]*3//4]
        y3.append(np.mean(tmp))

    pos = np.where(abs(x0 - x3[-1]) < 0.001)
    ratio = y0[pos]/y3[-1]
    p3 = plt.plot(np.array(x3), np.array(y3) * ratio, label = 'from 0.4Å to -0.4Å @20180719')

    # scan from -0.4 to 0.4 @20180719
    measured_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_FLAT_PINHOLES&FOCUS/LYOT_FILTER_PROFILE.nosync/20180718/032811'

    measured_files = [join(measured_path, f) for f in listdir(measured_path) if isfile(join(measured_path, f)) and f.endswith('.fits')]
    measured_files.sort()
    x3 = []
    y3 = []
    for f in measured_files:
        header = fits.getheader(f)
        x3.append(float(header['OFFBAND']) + center - 0.09)

        tmp = fits.getdata(f)
        shape = tmp.shape
        tmp = tmp[shape[0]//4:shape[0]*3//4, shape[1]//4:shape[1]*3//4]
        y3.append(np.mean(tmp))

    pos = np.where(abs(x0 - x3[4]) < 0.001)
    ratio = y0[pos]/y3[4]
    p4 = plt.plot(np.array(x3), np.array(y3) * ratio, label = 'from -0.4Å to 0.4Å @20180719')

    # scan from 0.4 to -0.4 @20180723
    measured_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_FLAT_PINHOLES&FOCUS/LYOT_FILTER_PROFILE.nosync/20180723/032304'

    measured_files = [join(measured_path, f) for f in listdir(measured_path) if isfile(join(measured_path, f)) and f.endswith('.fits')]
    measured_files.sort()
    x4 = []
    y4 = []
    for f in measured_files:
        header = fits.getheader(f)
        x4.append(float(header['OFFBAND']) + center - 0.07)

        tmp = fits.getdata(f)
        shape = tmp.shape
        tmp = tmp[shape[0]//4:shape[0]*3//4, shape[1]//4:shape[1]*3//4]
        y4.append(np.mean(tmp))

    pos = np.where(abs(x0 - x3[2]) < 0.001)
    ratio = y0[pos]/y3[2]
    p4 = plt.plot(np.array(x3), np.array(y3) * ratio, label = 'from 0.4Å to -0.4Å @20180723')

    # scan from -0.4 to 0.4 @20180723
    measured_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_FLAT_PINHOLES&FOCUS/LYOT_FILTER_PROFILE.nosync/20180723/032745'

    measured_files = [join(measured_path, f) for f in listdir(measured_path) if isfile(join(measured_path, f)) and f.endswith('.fits')]
    measured_files.sort()
    x5 = []
    y5 = []
    for f in measured_files:
        header = fits.getheader(f)
        x5.append(float(header['OFFBAND']) + center - 0.07)

        tmp = fits.getdata(f)
        shape = tmp.shape
        tmp = tmp[shape[0]//4:shape[0]*3//4, shape[1]//4:shape[1]*3//4]
        y5.append(np.mean(tmp))

    pos = np.where(abs(x0 - x3[0]) < 0.001)
    ratio = y0[pos]/y3[0]
    p5 = plt.plot(np.array(x3), np.array(y3) * ratio, label = 'from -0.4Å to 0.4Å @20180723')

    plt.legend()

    plt.savefig('/tmp/lyot_filter_profile.eps')

    plt.show()

def ProfileOfUTF(path):
    data = fits.getdata(path)
    size = data.shape
    profile = []
    for i in range(size[0]):
        tmp = data[i, :, :]
    plt.legend()

    plt.savefig('/tmp/lyot_filter_profile.eps')

    plt.show()

def ProfileOfUTF(path):
    data = fits.getdata(path)
    size = data.shape
    profile = []
    for i in range(size[0]):
        tmp = data[i, :, :]
        profile.append(np.mean(tmp))
        fig = plt.figure()
        plt.imshow(tmp, cmap='gray', origin='lower')
        fig.savefig('/tmp/2/' + f'{i:03}' + '.jpeg')

    #plt.imshow(data[0,:,:], origin='lower', cmap='gray')
    '''
    for i in range(size[0]):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(7, 3))
        axes[0].set_xticks([])
        axes[0].set_xticklabels([])
        axes[0].set_yticks([])
        axes[0].set_yticklabels([])
        axes[0].imshow(data[i, :, :], cmap='gray', origin='lower')
        axes[1].plot(profile)
        plt.tight_layout()
        plt.show()
        return
    '''
    return
def Temperature():
    path = '/Users/yuchao/Documents/memo/20180526_co_observation_with_hida_III/TEMPERATURE.NOSYNC/utfTlog_20180607.csv'
    data = np.genfromtxt(path, skip_header = 4800, delimiter=',')

    x = np.linspace(14, 17, 1689)
    print(x)
    plt.figure()
    plt.plot(x, data[:, 2])
    plt.title('1/8')
    plt.ylabel('Temperature')
    plt.xlabel('Time')
    plt.figure()
    plt.plot(x, data[:, 1])
    plt.title('UTF32')
    plt.xlabel('Time')
    plt.ylabel('Temperature')
    plt.show()
    return
if __name__ == '__main__':
    FeIProfileAt5324A()
    '''
    path = '/Users/yuchao/Documents/memo/20180526_HIDA/PROFILE.NOSYNC/WITH_SPIDER/'
    path = '/Users/yuchao/Documents/memo/20180526_co_observation_with_hida_III/PROFILE.NOSYNC/WITH_SPIDER/'
    path += 'Ca8542_20180607_141544.991.fits'
    path += 'Ha6563_20180607_140821.385.fits'
    #path += 'Mg5172_20180607_142838.729.fits'
    ProfileOfUTF(path)
    '''
    #Temperature()


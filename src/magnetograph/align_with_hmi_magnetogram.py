# -*- coding: utf-8 -*-

import astropy.io.fits as fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from astropy.visualization import ZScaleInterval
import scipy.ndimage.filters as scipyfilters
import scipy.ndimage as scipyimage
import scipy.misc as scipymisc
import drms
import sunpy.map as mp
import os
import math

def Align(im0, im1):
    return
# segment: continuum, magnetogram
def WriteFits(filename, data):
    if os.path.exists(filename):
        os.remove(filename)
    fits.writeto(filename, data)
    return
def RetrieveHMIFile(year, month, day, hour, minute, second, segment):
    query = 'hmi.'
    series = ''
    if segment == 'continuum':
        series = 'Ic_45s'
    else:
        series = 'M_45s'
    query += series + '[' + year + '.' + month + '.' + day + '_' + hour + ':' + minute + ':' + second + '_TAI]'
    #query = 'hmi.M_45s[2018.05.09_01:30:00_TAI]'

    client = drms.Client()
    hmi_keys, segments = client.query(query,
            key = drms.const.all, seg = segment)
    hmi_url = 'http://jsoc.stanford.edu'
    if segment == 'continuum':
        hmi_url += segments.continuum[0]
    else:
        hmi_url += segments.magnetogram[0]
    data = fits.getdata(hmi_url)
    header = dict(hmi_keys.iloc[0])
    path = '/tmp/'
    path += 'hmi.' + series + '.' + year + month + day + '_' + hour + minute + second + '.fits'
    fits.writeto(path, data)
    return data, header

def AlignWithHMI():
    #continuum, header = RetrieveHMIFile('2018', '05', '09', '01', '30', '00', 'continuum')
    continuum, header = RetrieveHMIFile('2018', '07', '22', '05', '57', '06', 'continuum')
    WriteFits('/tmp/continuum.fits', continuum)
    magnetogram, header = RetrieveHMIFile('2018', '07', '22', '05', '57', '06', 'magnetogram')
    WriteFits('/tmp/magnetogram.fits', magnetogram)
    '''
    magnetogram, header = RetrieveHMIFile('2018', '05', '09', '01', '30', '00', 'magnetogram')
    WriteFits('/tmp/magnetogram.fits', magnetogram)
    WriteFits('/tmp/continuum.fits', continuum)
    hmi_magneotgram = fits.getdata('/tmp/hmi.M_45s.20180509_013000.fits')
    hmi_continuum = fits.getdata('/tmp/hmi.Ic_45s.20180509_013000.fits')
    hmi_magneotgram = hmi_magneotgram[2020:2580, 1800:2300]
    hmi_continuum = hmi_continuum[2020:2580, 1800:2300]

    hmi_continuum = scipyimage.rotate(hmi_continuum, -22)
    hmi_magneotgram = scipyimage.rotate(hmi_magneotgram, -22)
    #fits.writeto('/tmp/tmp.fits', hmi_continuum)

    new_width = int(hmi_continuum.shape[0] * 6.4)
    new_height = int(hmi_continuum.shape[1] * 6.4)
    hmi_continuum = scipymisc.imresize(hmi_continuum, [new_width, new_height], interp = 'cubic')
    hmi_magneotgram = scipymisc.imresize(hmi_magneotgram, [new_width, new_height], interp = 'cubic')

    #WriteFits('/tmp/tmp.fits', hmi_continuum)

    path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180509/STOKES.nosync/STOKES_1/'
    continuum = fits.getdata(path + '00000_I.fits')
    q = fits.getdata(path + '00000_Q.fits')
    u = fits.getdata(path + '00000_U.fits')
    continuum = continuum[200:1600, 200:1600]
    q = q[200:1600, 200:1600]
    u = u[200:1600, 200:1600]
    u = scipyfilters.gaussian_filter(u, 0.8)
    u0 = fits.getdata('/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180509/STOKES.nosync/STOKES_BY_SPECKLE/00000_U.fits')
    u0 = u0[200:1600, 200:1600]
    i0 = fits.getdata('/tmp/F_00000_00.fits')
    i1 = fits.getdata('/tmp/F_00001_00.fits')
    i2 = fits.getdata('/tmp/F_00002_00.fits')
    i3 = fits.getdata('/tmp/F_00003_00.fits')
    i0 = i0[200:1600, 200:1600]
    i1 = i1[200:1600, 200:1600]
    i2 = i2[200:1600, 200:1600]
    i3 = i3[200:1600, 200:1600]
    #continuum = continuum[800:1200, 800:1200]

    #new_x0 = int(hmi_continuum.shape[0] / 2 - continuum.shape[0] / 2)
    #new_y0 = int(hmi_continuum.shape[1] / 2 - continuum.shape[1] / 2)
    #new_x1 = int(new_x0 + continuum.shape[0] / 2)
    #new_y1 = int(new_y0 + continuum.shape[1] / 2)

    new_x0 = int(hmi_continuum.shape[0] / 2 - continuum.shape[0] / 2)
    new_y0 = int(hmi_continuum.shape[1] / 2 - continuum.shape[1] / 2)
    new_x1 = int(new_x0 + continuum.shape[0])
    new_y1 = int(new_y0 + continuum.shape[1])

    hmi_continuum = hmi_continuum[new_x0:new_x1, new_y0:new_y1]
    hmi_magneotgram = hmi_magneotgram[new_x0:new_x1, new_y0:new_y1]

    continuum = continuum[500-340:, 500-440:]
    q = q[500-340:, 500-440:]
    u = u[500-340:, 500-440:]
    u0 = u0[500-340:, 500-440:]
    i0 = i0[500-340:, 500-440:]
    i1 = i1[500-340:, 500-440:]
    i2 = i2[500-340:, 500-440:]
    i3 = i3[500-340:, 500-440:]
    u = u * -1
    u0 = u0 * -1
    hmi_continuum = hmi_continuum[0:continuum.shape[0], 0:continuum.shape[1]]
    hmi_magneotgram = hmi_magneotgram[0:continuum.shape[0], 0:continuum.shape[1]]

    continuum = continuum / np.max(continuum)
    hmi_continuum = hmi_continuum / np.max(hmi_continuum)


    ims = []
    fig, axes = plt.subplots()
    ims.append([plt.imshow(continuum, cmap='gray', origin='lower', animated = True)])
    ims.append([plt.imshow(hmi_continuum, cmap='gray', origin='lower', animated = True)])
    ims.append([plt.imshow(hmi_magneotgram, cmap='gray', origin='lower', animated = True)])
    ims.append([plt.imshow(q, cmap='gray', origin='lower', animated = True)])
    ims.append([plt.imshow(u, cmap='gray', origin='lower', animated = True)])

    ani = animation.ArtistAnimation(fig, ims, interval=1000, repeat_delay=1000)
    #fits.writeto('/tmp/u0.fits', u0)
    #fits.writeto('/tmp/hmi.fits', hmi_magneotgram)
    #fits.writeto('/tmp/q.fits', q)
    #fits.writeto('/tmp/u.fits', u)
    fits.writeto('/tmp/i0.fits', i0)
    fits.writeto('/tmp/i1.fits', i1)
    fits.writeto('/tmp/i2.fits', i2)
    fits.writeto('/tmp/i3.fits', i3)
    #assert "timg" in result
    #ird.imshow(continuum, hmi_continuum, result['timg'], cmap='gray')
    #plt.figure()
    #plt.imshow(hmi_continuum, cmap = 'gray', origin = 'lower')
    #plt.figure()
    #plt.imshow(continuum, cmap = 'gray', origin = 'lower')
    #path = '/tmp/hmi.M_45s.20180509_013000_TAI.2.magnetogram.fits'
    #hmi_magnetograph = fits.getdata(path)
    #interval = ZScaleInterval()
    #hmi_magnetograph = interval(hmi_magnetograph)

    #plt.imshow(hmi_magnetograph, cmap = 'gray')

    #path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180509/STOKES.nosync/STOKES_1/'
    #filenames = [path + '00000_I.fits', path + '00000_Q.fits']
    #fig = plt.figure()
    #ims = []
    #for filename in filenames:
    #    data = fits.getdata(filename)
    #    interval = ZScaleInterval()
    #    data = interval(data)
    #    im = plt.imshow(data, cmap='gray', animated = True)
    #    ims.append([im])
    #ani = animation.ArtistAnimation(fig, ims, interval=300,
    #                            repeat_delay=1000)
    #hmimag = plt.get_cmap('hmimag')
    #plt.imshow(photosphere_data[1].data, cmap = hmimag, origin = 'lower')
    #interval = ZScaleInterval()
    #data = interval(data)

    #client = drms.Client()
    #hmi_keys, segments = client.query('hmi.M_45s[2018.05.09_01:30:00_TAI]',
    #        key = drms.const.all, seg = 'magnetogram')
    #hmi_url = 'http://jsoc.stanford.edu' + segments.magnetogram[0]
    #photosphere_full_image = fits.getdata(hmi_url)
    #header = dict(hmi_keys.iloc[0])
    #header['DATE-OBS'] = hmi_keys.DATE__OBS[0]
    #header['HGLN_OBS'] = 0.0
    #fits.writeto('/tmp/hmi.M_45s.20180509_013000_TAI.2.magnetogram.fits', photosphere_full_image)
    #mp.Map(photosphere_full_image, header).peek()
    plt.show()
    '''
    return

if __name__ == '__main__':
    #DownloadHMIFile('2018', '05', '09', '01', '30', '00', 'continuum')
    AlignWithHMI()

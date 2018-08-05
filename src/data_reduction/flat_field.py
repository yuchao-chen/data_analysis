# -*- coding: utf-8 -*-
import os
import shutil

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

import astropy.io.fits as fits
from astropy.visualization import ZScaleInterval

def histogram_qualize(data, nbin = 256):
    hist, bins = np.histogram(data.flatten(), nbin, normed = True)
    cdf = hist.cumsum()
    cdf = 255 * cdf / cdf[-1]
    #normalized_cdf = cdf * hist.max() / cdf.max()
    '''
    masked_cdf = np.ma.masked_equal(cdf, 0)
    masked_cdf = (masked_cdf - masked_cdf.min()) * 255 / (masked_cdf.max() - masked_cdf.min())
    cdf = np.ma.filled(masked_cdf, 0).astype('uint8')
    '''
    new_data = np.interp(data.flatten(), bins[:-1], cdf)
    return new_data.reshape(data.shape)

def difference_of_two_flat_fields0():
    path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_PINHOLES&FOCUS/20180718.nosync/B1/FLAT/'
    label = ['I0', 'I1', 'I2', 'I3']
    for l in label:
        '''
        dirs = [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and d.endswith(l)]
        '''
        dirs = ['/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_PINHOLES&FOCUS/FLATFIELD.nosync/20180718/B1/FLAT_' + label[0],
                '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_PINHOLES&FOCUS/FLATFIELD.nosync/20180718/B1/FLAT_' + label[3]]
        d0 = fits.getdata(dirs[0]+'/0015.fits')
        d1 = fits.getdata(dirs[1]+'/0015.fits')

        interval = ZScaleInterval()

        dif = d0 / d1
        std = np.std(dif[300:1500, 300:1600])

        fig, ax = plt.subplots(figsize = (5, 5))
        im = ax.imshow(histogram_qualize(dif) / 255.0, cmap = 'gray', origin = 'lower')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.text(100, 1800, 'std: ' + "{:.5f}".format(std), fontsize = 10)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax = cax)
        fig.savefig('/tmp/' + l + '.eps')
        plt.show()
        return

def difference_of_two_flat_fields1():
    path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_FLAT_PINHOLES&FOCUS/FLAT_FIELD.nosync/20180721'
    dirs = [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    dirs.sort(reverse = True)

    interval = ZScaleInterval()

    for root in dirs:
        subdirs = [os.path.join(root, d) for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))]
        subdirs.sort()
        files = []
        for subdir in subdirs:
            f = [os.path.join(subdir, f) for f in os.listdir(subdir) if os.path.isfile(os.path.join(subdir, f)) and f.endswith('.fits')]
            f.sort()
            files.append(f)
        files = list(map(list, zip(*files)))

        for index in range(1, len(files) - 1):
            ref = files[index - 1]
            cur = files[index]
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (7, 6))
            for r, c, ax in zip(ref, cur, axes.flat):
                d0 = fits.getdata(r)
                d1 = fits.getdata(c)
                dif = d1 / d0
                std = np.std(dif[300:1500, 300:1600])
                im = ax.imshow(interval(dif), cmap = 'gray', origin = 'lower')
                ax.text(100, 1800, 'std: ' + "{:.5f}".format(std), fontsize = 10)
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                if 'I00' in c:
                    ax.set_title('I0 (' + str(index * 50) + ' FRAMES)')
                elif 'I01' in c:
                    ax.set_title('I1 (' + str(index * 50) + ' FRAMES)')
                elif 'I02' in c:
                    ax.set_title('I2 (' + str(index * 50) + ' FRAMES)')
                elif 'I03' in c:
                    ax.set_title('I3 (' + str(index * 50) + ' FRAMES)')
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(im, cax = cax)
            fig.tight_layout()
            fig.savefig('/tmp/0/' + str(index).zfill(3) + '.eps')
        return
    return

def difference_of_two_flat_fields2():
    path = '/tmp/B1/'

    dirs = [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    dirs.sort(reverse = True)

    interval = ZScaleInterval()

    for root in dirs:
        subdirs = [os.path.join(root, d) for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))]
        subdirs.sort()
        files = []
        for subdir in subdirs:
            f = [os.path.join(subdir, f) for f in os.listdir(subdir) if os.path.isfile(os.path.join(subdir, f)) and f.endswith('.fits')]
            f.sort()
            files.append(f)
        files = list(map(list, zip(*files)))

        fig_path = '/tmp/' + os.path.basename(root) + '/'
        if os.path.exists(fig_path):
            shutil.rmtree(fig_path)
        os.mkdir(fig_path)

        titles = ['20180715', '20180718', '20180721', '20180722']
        for index in range(1, len(files)):
            ref = files[index - 1]
            cur = files[index]
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (7, 6))
            t = titles[index - 1] + '÷' + titles[index]
            for r, c, ax in zip(ref, cur, axes.flat):
                d0 = fits.getdata(r)
                d1 = fits.getdata(c)
                #dif = d1 * (index + 1) / (d0 * index)
                dif = d1 / d0
                std = np.std(dif[300:1500, 300:1600])
                im = ax.imshow(histogram_qualize(dif), cmap = 'gray', origin = 'lower')
                ax.text(100, 1800, 'std: ' + "{:.5f}".format(std), fontsize = 10)
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                title = 'I'
                if 'I00' in c:
                    title = 'I0'
                elif 'I01' in c:
                    title = 'I1'
                elif 'I02' in c:
                    title = 'I2'
                elif 'I03' in c:
                    title = 'I3'
                #title += ' (' +  str((index+1) * 1500) + ' FRAMES)'
                ax.set_title(title)
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(im, cax = cax)
            fig.suptitle(t)
            fig.tight_layout()
            fig.savefig(fig_path + str(index).zfill(3) + '.eps')
            plt.show()
    return

def main():
    difference_of_two_flat_fields2()

if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
import os
from functools import reduce

import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage

import astropy.io.fits as fits
from astropy.visualization import ZScaleInterval

from register_translation import register_translation

def write_fits(filename, data):
    if os.path.exists(filename):
        os.remove(filename)
    fits.writeto(filename, data)

def create_diectory(path):
    try:
        os.stat(path)
    except:
        os.makedirs(path)

def list_directories(path, suffix = None):
    return [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d)) and (suffix is None or d.endswith(suffix))]

def list_files_recursively(path, suffix = None):
    return [os.path.join(root, f) for root, dirs, files in os.walk(path) for f in files if suffix is None or f.endswith(suffix)]

def subpixel_shift(filename, shift, shifted_file_root):
    data = fits.getdata(filename)
    freq = np.fft.fftn(data)
    shifted_freq = ndimage.fourier_shift(freq, shift)
    shifted_data = np.real(np.fft.ifftn(shifted_freq))
    shifted_data = shifted_data.astype(np.uint16)
    # plt.imshow(shifted_data)
    # plt.show()
    subd0 = os.path.basename(os.path.dirname(filename))
    subd1 = os.path.basename(os.path.dirname(os.path.dirname(filename)))
    dir = os.path.join(shifted_file_root, subd1)
    dir = os.path.join(dir, subd0)
    create_diectory(dir)

    shifted_filename = os.path.join(dir, os.path.basename(filename))
    write_fits(shifted_filename, shifted_data)

def register(path0, path1, registered_path0, registered_path1, cropped_area = [1024-256, 1024-256, 512, 512]):
    # create_diectory(registered_path0)
    # create_diectory(registered_path1)

    files0 = list_files_recursively(path0, '.fits')
    files1 = list_files_recursively(path1, '.fits')
    files0.sort()
    files1.sort()
    # the first image is used as reference
    ref = fits.getdata(files0[0])
    ref = ref[cropped_area[0]: cropped_area[0]+cropped_area[2]-1,
          cropped_area[1]: cropped_area[1]+cropped_area[3]-1]
    ref = np.subtract(ref, np.mean(ref))
    ref_freq = np.fft.fftn(ref)

    for i in range(len(files0) - 1):
        target = fits.getdata(files0[i])
        target = target[cropped_area[0]: cropped_area[0] + cropped_area[2] - 1,
              cropped_area[1]: cropped_area[1] + cropped_area[3] - 1]
        target = np.subtract(target, np.mean(target))
        target_freq = np.fft.fftn(target)
        shift, error, diffphase = register_translation(ref_freq, target_freq, upsample_factor = 100, space = 'fourier')
        print(shift)
        subpixel_shift(files0[i], shift, registered_path0)
        subpixel_shift(files1[i], shift, registered_path1)

def sum(path, result_path):
    files0 = list_files_recursively(path, '.fits')
    sum0 = 0.0
    for f in files0:
        tmp = fits.getdata(f)
        tmp = tmp.astype(np.double)
        sum0 += tmp

    write_fits(result_path, sum0.astype(np.float32))

def deconvolution(b1_path, b2_path, b2_sum_filename, deconv_root, cropped_area = [1024-256, 1024-256, 512, 512]):
    b2_sum = fits.getdata(b2_sum_filename)
    b2_sum = b2_sum[cropped_area[0]: cropped_area[0]+cropped_area[2]-1,
          cropped_area[1]: cropped_area[1]+cropped_area[3]-1]
    # plt.imshow(b2_sum)
    # plt.show()
    b2_sum_freq = np.fft.fftn(b2_sum)
    b2_sum_phase = b2_sum_freq / np.abs(b2_sum_freq)
    files0 = list_files_recursively(b1_path, '.fits')
    files1 = list_files_recursively(b2_path, '.fits')
    files0.sort()
    files1.sort()
    tmp0 = []
    tmp1 = []
    for i in range(10):
        tmp0.append(files0[i::300])
        tmp1.append(files1[i::300])
    files0 = reduce(lambda x,y :x+y ,tmp0)
    files1 = reduce(lambda x,y :x+y ,tmp1)
    files0.sort()
    files1.sort()
    for f in files0:
        print(f)
    for i in range(len(files0) - 1):
        d0 = fits.getdata(files0[i])
        d1 = fits.getdata(files1[i])
        # interval = ZScaleInterval()
        # plt.imshow(interval(d0), cmap = 'gray', origin = 'lower')
        # plt.colorbar()
        # plt.show()
        d0 = d0[cropped_area[0]: cropped_area[0] + cropped_area[2] - 1,
                 cropped_area[1]: cropped_area[1] + cropped_area[3] - 1]
        d1 = d1[cropped_area[0]: cropped_area[0] + cropped_area[2] - 1,
                 cropped_area[1]: cropped_area[1] + cropped_area[3] - 1]

        d0_freq = np.fft.fftn(d0)
        d1_freq = np.fft.fftn(d1)
        d1_phase = d1_freq / np.abs(d1_freq)
        deconv_d0_freq = d0_freq * b2_sum_phase * d1_phase.conj() / np.abs(d1_phase)
        # deconv_d0_freq = d0_freq * b2_sum_freq * d1_freq.conj() / np.abs(d1_freq)
        deconv_d0 = np.real(np.fft.ifftn(deconv_d0_freq))
        smoothed_deconv_d0 = ndimage.filters.gaussian_filter(deconv_d0, 0.4)

        smoothed_deconv_d0 = smoothed_deconv_d0.astype(np.uint16)
        dir = os.path.join(deconv_root, os.path.basename(os.path.dirname(files0[i])))
        create_diectory(dir)

        deconv_filename = os.path.join(dir, os.path.basename(files0[i]))
        write_fits(deconv_filename, smoothed_deconv_d0)

def sum_all_files_in_current_directory(path, result_root):
    create_diectory(result_root)
    dirs = list_directories(path)
    for d in dirs:
        filename = os.path.join(result_root, os.path.basename(d)+'.fits')
        sum(d, filename)

def demod(path, result_dir):
    create_diectory(result_dir)
    files = list_files_recursively(path, '.fits')
    files.sort()
    for index in range(0, len(files)//4):
        d0 = fits.getdata(files[4 * index + 0])
        d1 = fits.getdata(files[4 * index + 1])
        d2 = fits.getdata(files[4 * index + 2])
        d3 = fits.getdata(files[4 * index + 3])

        i = 0.419 * d0 + 0.581 * d1 + 0.581 * d2 + 0.418 * d3
        q = 0.889 * d0 - 0.887 * d1 - 0.888 * d2 + 0.886 * d3
        u = 1.131 * d0 - 0.547 * d1 + 0.549 * d2 - 1.133 * d3
        v = -0.687 * d0 - 0.995 * d1 + 0.994 * d2 + 0.688 * d3
        qi = q/i
        ui = u/i
        vi = v/i
        write_fits(os.path.join(result_dir, str(index).zfill(4)+'_I.fits'), i)
        write_fits(os.path.join(result_dir, str(index).zfill(4) + '_Q.fits'), qi.astype(np.float32))
        write_fits(os.path.join(result_dir, str(index).zfill(4) + '_U.fits'), ui.astype(np.float32))
        write_fits(os.path.join(result_dir, str(index).zfill(4) + '_V.fits'), vi.astype(np.float32))
        plt.figure()
        plt.imshow(i, cmap='gray', origin='lower')
        plt.figure()
        plt.imshow(q/i, cmap='gray', origin='lower')
        plt.colorbar()

        plt.figure()
        plt.imshow(u/i, cmap='gray', origin='lower')
        plt.colorbar()

        plt.figure()
        plt.imshow(v/i, cmap='gray', origin='lower')
        plt.colorbar()
        plt.show()

def demod_with_deconv(b1_path, b2_path, result_path):
    registered_b1_path = 'G:\\20180807\\B1\\12717\\063007\\B0080_AFTER_ALIGNED_1'#os.path.join(result_path, 'REGISTERED\\B1')
    registered_b2_path = 'G:\\20180807\\B2\\12717\\063002\\CENT_AFTER_ALIGNED_1'#os.path.join(result_path, 'REGISTERED\\B2')

    register(b2_path, b1_path, registered_b2_path, registered_b1_path)
    # b2_sum_filename = os.path.join(result_path, 'SUM.fits')
    # # sum(registered_b2_path, b2_sum_filename)
    # b1_decon_root = os.path.join(result_path, 'DECONV')
    # # deconvolution(registered_b1_path, registered_b2_path, b2_sum_filename, b1_decon_root, cropped_area=[250, 200, 1600, 1600])
    # # deconvolution(registered_b1_path, registered_b2_path, b2_sum_filename, b1_decon_root, cropped_area=[0, 0, 0, 0])
    #
    # b1_sum_root = os.path.join(result_path, 'SUM1')
    # b2_sum_root = os.path.join(result_path, 'SUM2')
    # sum_all_files_in_current_directory(registered_b1_path, b1_sum_root)
    # sum_all_files_in_current_directory(registered_b2_path, b2_sum_root)
    #
    # deconvolution(b1_sum_root, b2_sum_root, b2_sum_filename, b1_decon_root,
    #               cropped_area=[250, 200, 1600, 1600])
    #
    # b1_stokes_root = os.path.join(result_path, 'STOKES')
    # demod(b1_sum_root, b1_stokes_root)

def main():
    b1_path = 'G:\\20180807\\B1\\12717\\063007\\B0080_AFTER_ALIGNED_0'
    b2_path = 'G:\\20180807\\B2\\12717\\063002\\CENT_AFTER_ALIGNED_0'
    result_path = 'G:\\20180722\\B1\\QS\\055706\\DECONVOLUTION'
    demod_with_deconv(b1_path, b2_path, result_path)

if __name__ == '__main__':
    main()

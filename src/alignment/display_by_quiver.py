# -*- coding: utf-8 -*-
# display pinholes by quiver
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
from pinholes_based_align import PinholesBasedAlign
import cv2
import math
import os
import pickle

def quiver_pinholes(ref_data, to_align_data, save_path = None, title = None, draw_circles = None):
    p = PinholesBasedAlign()
    ref_pinholes = p.count_pinholes(ref_data)
    to_align_pinholes = p.count_pinholes(to_align_data)

    ref_keys, to_align_keys = p._match(ref_pinholes, to_align_pinholes, 5)
    if ref_keys is None or to_align_keys is None:
        return 100

    x = []
    y = []
    u = []
    v = []

    distance = []
    for k0, k1 in zip(ref_keys, to_align_keys):
        x.append(k0[0])
        y.append(k0[1])
        u.append(k1[0] - k0[0])
        v.append(k1[1] - k0[1])

        distance.append(math.sqrt((k1[0] - k0[0])**2 + (k1[1] - k0[1])**2))

    print(x)
    print(y)
    print(u)

    m = np.hypot(u, v)

    std = np.sum(distance)

    fig, axe = plt.subplots(figsize = (5, 5))

    if draw_circles is not None:
        for cp in draw_circles:
            cir = Circle(cp, 50, color='lightskyblue')
            axe.add_patch(cir)

    axe.quiver(x, y, u, v, m, units = 'dots', pivot = 'tip', scale = 1)
    axe.text(100, 1900, '$\sum ||r||^2$: ' + str(std), fontsize = 9)
    axe.axis('equal')
    axe.set_xlim(0, 2048)
    axe.set_ylim(0, 2048)
    axe.set_xticks([0, 1024, 2048])
    axe.set_yticks([0, 1024, 2048])
    axe.set_xticklabels([])
    axe.set_yticklabels([])

    if title is not None:
        axe.set_title(title)

    if save_path is not None:
        fig.savefig(save_path)

    plt.show()
    return std

def display_two_channels():
    fig_path = '/tmp/0'
    try:
        os.stat(fig_path)
    except:
        os.mkdir(fig_path)

    ref_dark_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180704.nosync/B1/DARK.fits'
    ref_dark = fits.getdata(ref_dark_path)
    ref_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180704.nosync/B1/PINHOLES.fits'
    ref_image = fits.getdata(ref_path)[0,:,:]
    ref_image = ref_image.astype(float)
    ref_image -= ref_dark

    to_align_dark_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180704.nosync/B2/DARK.fits'
    to_align_dark = fits.getdata(to_align_dark_path)
    to_align_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180704.nosync/B2/PINHOLES.fits'
    to_align_image = fits.getdata(to_align_path)[0,:,:]
    to_align_image = to_align_image.astype(float)
    to_align_image -= to_align_dark

    shape = to_align_image.shape

    p = PinholesBasedAlign()
    '''
    # before align
    quiver_pinholes(ref_image, to_align_image, '/tmp/before_align.eps', title = 'BEFORE ALIGNMENT')

    h = p.find_homography(ref_image, to_align_image, estimator = cv2.RANSAC, rejection_threshold = 1)
    aligned_image = cv2.warpPerspective(to_align_image, h, (shape[0], shape[1]))
    quiver_pinholes(ref_image, aligned_image, save_path = '/tmp/homography_ransac.eps', title = 'PERSPECTIVE TRANSFORM WITH RANSAC')

    h = p.find_homography(ref_image, to_align_image, estimator = cv2.LMEDS)
    aligned_image = cv2.warpPerspective(to_align_image, h, (shape[0], shape[1]))
    quiver_pinholes(ref_image, aligned_image, save_path = '/tmp/homography_lmeds.eps', title = 'PERSPECTIVE TRANSFORM WITH LMEDS')

    # first rough register
    h = p.find_afffine_transform0(ref_image, to_align_image)
    aligned_image = cv2.warpAffine(to_align_image, h, (shape[0], shape[1]))
    quiver_pinholes(ref_image, aligned_image, save_path = '/tmp/shift_with_average_dx_dy.eps', title = 'AFFINE TRANSFORM')
    '''
    # second register
    selected_points = [0, 60, 11]
    h, selected_points = p.find_afffine_transform(ref_image, to_align_image, selected_points = selected_points)
    aligned_image = cv2.warpAffine(to_align_image, h, (shape[0], shape[1]))
    quiver_pinholes(ref_image, aligned_image, save_path = '/tmp/affine_transform.eps', title = 'AFFINE TRANSFORM')
    p._write_fits('/tmp/aligned.fits', aligned_image)
    p._write_fits('/tmp/ref.fits', ref_image)
    '''
    h, selected_points = p.find_afffine_transform(ref_image, to_align_image, selected_points = selected_points)
    aligned_image = cv2.warpAffine(ref_image, h, (shape[0], shape[1]))
    p._write_fits('/tmp/ref.fits', ref_image)
    quiver_pinholes(ref_image, aligned_image, save_path = '/tmp/homography_ransac.eps', title = 'PERSPECTIVE TRANSFORM WITH RANSAC')
    '''
    '''
    count = 94
    stds = np.zeros((count, count, count))
    #for i in range(0, count//2):
    #    for j in range(-1, i - count + 1, -1):
    #        for k in range(i + 1, count + j - 1):
    for i in [0]:#range(0,count-1):
        for j in range(1, count - 1):
            #for k in range(i + 1, count + j - 1):
            for k in [40]:
        #for j in range(-5, i - count + 1, -1):
        #    for k in range(26, count + j - 1):
                print(i, j, k)
                selected_points = [i, j, k]
                h, selected_points = p.find_afffine_transform(ref_image, to_align_image, selected_points = selected_points)
                if h is not None:
                    aligned_image = cv2.warpAffine(to_align_image, h, (shape[0], shape[1]))
                    std = quiver_pinholes(ref_image, aligned_image, save_path = fig_path + '/affine_transform_' + str(i).zfill(3) + str(j).zfill(3) + str(k).zfill(3) + '.eps', title = 'AFFINE TRANSFORM', draw_circles = selected_points)
                    stds[i, j, k] = std
    with open('/tmp/stds.dat', 'wb') as f:
        pickle.dump(stds, f, protocol=pickle.HIGHEST_PROTOCOL)
    '''
    '''
    store = []
    count = 94
    p = PinholesBasedAlign()
    for i in range(0, count - 1):
        selected_points = [i, 60, 11]
        h, selected_points = p.find_afffine_transform(ref_image, to_align_image, selected_points = selected_points)
        if h is not None:
            aligned_image = cv2.warpAffine(to_align_image, h, (shape[0], shape[1]))
            r = quiver_pinholes(ref_image, aligned_image, save_path = fig_path + '/affine_transform_' + str(i).zfill(3) + '.eps', title = 'AFFINE TRANSFORM', draw_circles = selected_points)
            store.append(r)
    for index, s in enumerate(store):
        print(index, s)
    '''
def display1():
    fig_path = '/tmp/1'
    try:
        os.stat(fig_path)
    except:
        os.mkdir(fig_path)

    dark_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180704.nosync/B2/DARK.fits'
    dark0 = fits.getdata(dark_path)

    dark_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180704.nosync/B1/DARK.fits'
    dark1 = fits.getdata(dark_path)

    path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_PINHOLES&FOCUS/PINHOLES.nosync/20180714/B2'
    dirs0 = [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

    path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_PINHOLES&FOCUS/PINHOLES.nosync/20180714/B1'
    dirs1 = [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

    dirs0.sort()
    dirs1.sort()

    p = PinholesBasedAlign()
    selected_points = [0, 59, 17]

    for d0, d1 in zip(dirs0, dirs1):
        t = []

        files0 = []
        for root, subdirs, fs in os.walk(d0):
            for f in fs:
                if os.path.isfile(os.path.join(root, f)) and f.endswith('.fits'):
                    files0.append(os.path.join(root, f))
        files1 = []
        for root, subdirs, fs in os.walk(d1):
            for f in fs:
                if os.path.isfile(os.path.join(root, f)) and f.endswith('.fits'):
                    files1.append(os.path.join(root, f))

        files0.sort()
        files1.sort()

        for f0, f1 in zip(files0, files1):
            data0 = fits.getdata(f0)
            data0 = data0.astype(float)
            data0 -= dark0

            data1 = fits.getdata(f1)
            data1 = data1.astype(float)
            data1 -= dark1
            h, dump = p.find_afffine_transform(data1, data0, selected_points = selected_points)
            #quiver_pinholes(ref, data, save_path = new_fig_path + '/' + str(index).zfill(3) + '.eps', title = 'Derotation OFF + Polarimeter ON')
            t.append(h)

        with open(fig_path + '/' + os.path.basename(d0) + 'transform.dat', 'wb') as f:
            pickle.dump(t, f, protocol=pickle.HIGHEST_PROTOCOL)

def display0():
    fig_path = '/tmp/0'
    try:
        os.stat(fig_path)
    except:
        os.mkdir(fig_path)

    dark_path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/6TH_EXPERIMENT/20180704.nosync/B2/DARK.fits'
    dark = fits.getdata(dark_path)
    path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_PINHOLES&FOCUS/PINHOLES.nosync/20180714/B2'
    dirs = [os.path.join(path, d) for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

    p = PinholesBasedAlign()
    selected_points = [0, 59, 17]

    ref = None
    for d in dirs:
        t = []
        files = []
        new_fig_path = fig_path + '/' + os.path.basename(d)
        try:
            os.stat(new_fig_path)
        except:
            os.mkdir(new_fig_path)

        for root, subdirs, fs in os.walk(d):
            for f in fs:
                if os.path.isfile(os.path.join(root, f)) and f.endswith('.fits'):
                    files.append(os.path.join(root, f))
        for index, f in enumerate(files):
            data = fits.getdata(f)
            data = data.astype(float)
            data -= dark
            if ref is None:
                ref = data
            h, dump = p.find_afffine_transform(ref, data, selected_points = selected_points)
            #quiver_pinholes(ref, data, save_path = new_fig_path + '/' + str(index).zfill(3) + '.eps', title = 'Derotation OFF + Polarimeter ON')
            t.append(h)

        with open(new_fig_path + 'transform.dat', 'wb') as f:
            pickle.dump(t, f, protocol=pickle.HIGHEST_PROTOCOL)


def main():
    display_two_channels()

if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
# alignment of two pinholearrays
#
import os
from collections import defaultdict
import math

import pickle

import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt

import astropy.io.fits as fits

import cv2

from thresholding import otsu_threshold


class PinholesBasedAlign():
    def __init__(self):
        pass

    @classmethod
    def _center_of_mass(self, data, found_regions):
        no = np.max(found_regions)
        center_of_found_regions = ndimage.measurements.center_of_mass(data, found_regions, range(1, no + 1))
#        return center_of_found_regions
        # dump code from here
        center_of_mass = []
        r = 20
        for p in center_of_found_regions:
            start_x = int(p[0] - r)
            end_x = int(p[0] + r)
            start_y = int(p[1] - r)
            end_y = int(p[1] + r)

            if start_x > 0 and start_y > 0:
                roi = data[start_x:end_x, start_y:end_y]
                #mc = ndimage.measurements.center_of_mass(roi)
                mc = ndimage.measurements.maximum_position(roi)
                center_of_mass.append((mc[0] + start_x, mc[1] + start_y))
                #xx = [1032.2545195328107, 710.9350263774473]
                #if abs(mc[0] + start_x - xx[0]) < 1 and abs(mc[1] + start_y - xx[1]) < 1:
                #    self._write_fits('/tmp/roi.fits', roi)
        return center_of_mass

    @classmethod
    def count_pinholes(self, data):
        # remove noise
        threshold = otsu_threshold(data)
        binary = data > threshold
        open_binary = ndimage.binary_opening(binary)
        eroded_binary = ndimage.binary_erosion(open_binary)
        recon_binary = ndimage.binary_propagation(eroded_binary, mask = binary)

        found_pinholes, found_pinholes_no = ndimage.label(recon_binary)
        #found_pinholes, found_pinholes_no = ndimage.label(open_binary)
        print("FOUND PINHOLES: " + str(found_pinholes_no))

    #    # remove the pinholes with small size
    #    sizes = ndimage.sum(recon_binary, found_pinholes, range(found_pinholes_no + 1))
    #    mean_size = ndimage.sum(data, found_pinholes, range(1, found_pinholes_no + 1))
    #    mask_size = sizes < 15#np.mean(mean_size)
    #    remove_pinholes = mask_size[found_pinholes]
    #    found_pinholes[remove_pinholes] = 0
    #    pinholes = np.unique(found_pinholes)
    #    found_pinholes = np.searchsorted(pinholes, found_pinholes)

        center_of_mass = self._center_of_mass(data, found_pinholes)

        return center_of_mass

    @classmethod
    def find_afffine_transform0(self, ref_image, to_align_image, selected_points = None, threshold = 5):
        ref_pinholes = self.count_pinholes(ref_image)
        to_align_pinholes = self.count_pinholes(to_align_image)

        ref_keys, to_align_keys = self._match(ref_pinholes, to_align_pinholes, threshold)

        if ref_keys is None or to_align_keys is None:
            return None, None

        tx = 0.0
        ty = 0.0
        for r, t in zip(ref_keys, to_align_keys):
            tx += t[0] - r[0]
            ty += t[1] - r[1]
        tx /= len(ref_keys)
        ty /= len(ref_keys)
        h = np.array([[1, 0, tx], [0, 1, ty]])
        print('============= Affine Transform ===============')
        print(h)
        print('=======================================')
        return h

    # available estimator: cv2.LMEDS, cv2.RANSAC
    @classmethod
    def find_afffine_transform(self, ref_image, to_align_image, selected_points = None, threshold = 5):
        ref_pinholes = self.count_pinholes(ref_image)
        to_align_pinholes = self.count_pinholes(to_align_image)

        ref_keys, to_align_keys = self._match(ref_pinholes, to_align_pinholes, threshold)

        if ref_keys is None or to_align_keys is None:
            return None, None

        if selected_points is None:
            selected_points = []
            selected_points.append(0)
            mid = len(ref_pinholes) // 2
            selected_points.append(mid)
            selected_points.append(-1)

        selected_ref_keys = []
        selected_to_align_keys = []
        for p in selected_points:
            selected_ref_keys.append(ref_keys[p])
            selected_to_align_keys.append(to_align_keys[p])

        h = cv2.getAffineTransform(np.float32(selected_ref_keys), np.float32(selected_to_align_keys))
        print('============= Affine Transform ===============')
        print(h)
        print('=======================================')
        return h, selected_ref_keys

    # available estimator: cv2.LMEDS, cv2.RANSAC
    @classmethod
    def find_homography(self, ref_image, to_align_image, threshold = 5, estimator = cv2.LMEDS, rejection_threshold = 0):
        ref_image_copy = ref_image
        to_align_image_copy = to_align_image

        ref_pinholes = self.count_pinholes(ref_image_copy)
        to_align_pinholes = self.count_pinholes(to_align_image_copy)

        ref_keys, to_align_keys = self._match(ref_pinholes, to_align_pinholes, threshold)

        h, mask = cv2.findHomography(np.array(ref_keys), np.array(to_align_keys), estimator, rejection_threshold)
        print('============= Homography ===============')
        print(h)
        print('=======================================')
        return h

    @classmethod
    def _distance(self, p0, p1):
        d = math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[0])**2)
        return d

    @classmethod
    def _cluster_pinholes(self, points, threshold):
        #print('======================= CLUTSERING ===========================')
        cluster = []
        row = []
        for p in points:
            if (len(row) > 0):
                previous = row[-1]
                if (abs(previous[0] - p[0]) > threshold):
                    cluster.append(row)
                    row = []
            row.append(p)
        return cluster

    @classmethod
    def _find_unique_point(self, points, threshold):
        key = []
        cluster = self._cluster_pinholes(points, threshold)
        pos = 0
        for index, c in enumerate(cluster):
            if len(c) == 1:
                key.append([c[0][0], c[0][1]])
                pos = index
        return key, cluster, pos

    @classmethod
    def _map_to_grid(self, points, threshold = 0.2):
        grid = []
        selected_angles = []

        points_copy = points
        # the pinholes are distributed in grid,
        # there are two main gradients, m and -1/m
        for index, p0 in enumerate(points_copy):

            angles = defaultdict(list)

            if index + 1 > len(points) - 1:
                break
            left_points = points[index+1:]

            # perpendicular point pairs
            pairs = []

            for p1 in left_points:
                angle = self._gradient(p0, p1)
                not_found = True
                for k in angles:
                    # check the parallel points
                    if abs(k - angle) > threshold:
                        if abs(abs(k - angle) - 90.0) < threshold:
                            pairs.append(k)
                            pairs.append(angle)
                        continue
                    angles[k].append(p1)
                    not_found = False
                    break
                if not_found:
                    angles[angle].append(p1)

            max_leng_angle = None
            for k in angles:
                if max_leng_angle is None:
                    max_leng_angle = k
                else:
                    if len(angles[k]) > len(angles[max_leng_angle]):
                        max_leng_angle = k
            if max_leng_angle is not None:
                selected_angles.append(max_leng_angle)
                selected_points = angles[max_leng_angle]
                selected_points.append(p0)
                selected_points.sort(key = lambda p: math.sqrt((p0[0] - p[0])**2 + (p0[1] - p[1])**2))

                grid.append(selected_points)
                if max_leng_angle in pairs:
                    index = pairs.index(max_leng_angle)
                    if index % 2 == 1:
                        index -= 1
                    else:
                        index += 1
                    pp = angles[pairs[index]]
                    selected_angles.append(pairs[index])
                    pp.append(p0)
                    pp.sort(key = lambda p: math.sqrt((p0[0] - p[0])**2 + (p0[1] - p[1])**2), reverse = True)
                    for index, p in enumerate(pp):
                        if p == p0:
                            continue
                        grid.insert(index, [p])

            if len(grid) > 1:
                break
        # remove the duplicate points
        for row in grid:
            for e in row:
                points_copy.remove(e)
        not_found_points = []
        for p0 in points_copy:
            not_found = True
            for index, s in enumerate(grid):
                angle = self._gradient(s[0], p0)
                if abs(angle - selected_angles[0]) < threshold:
                    not_found = False
                    grid[index].append(p0)
                    break
            if not_found:
                not_found_points.append(p0)

        longest_row = []
        for g in grid:
            if len(g) > len(longest_row):
                longest_row = g
        print(longest_row)

    @classmethod
    def _match(self, points0, points1, threshold = 5):

        key0, cluster0, index0 = self._find_unique_point(points0, threshold)
        key1, cluster1, index1 = self._find_unique_point(points1, threshold)

        if len(key0) == 0 or len(key1) == 0:
            return None, None
        print('THE RECOGNIZED UNIQUE POINT:')
        print(key0, index0)
        print(key1, index1)
        distance = key1[0][1] - key0[0][1]

        if (index0 != index1):
            print('Error: the position of the unique point is not same')
            return None, None

        sorted_key0 = []
        sorted_key1 = []

        for c0, c1 in zip(cluster0, cluster1):
            c1_duplicate = c1
            for c00 in c0:
                selected_p = None
                for c11 in c1_duplicate:
                    if abs(c11[1] - c00[1] - distance) < threshold:
                        selected_p = c11
                        c1_duplicate.remove(selected_p)
                        break
                if selected_p != None:
                    #if self._distance(c00, key0[0]) < 400:
                    sorted_key0.append([c00[0], c00[1]])
                    sorted_key1.append([selected_p[0], selected_p[1]])
        return sorted_key0, sorted_key1

    @classmethod
    def _write_fits(self, filename, data):
        if os.path.exists(filename):
            os.remove(filename)
        fits.writeto(filename, data)

    @classmethod
    def _gradient(self, p0, p1):
        a = p0[0] - p1[0]
        b = p0[1] - p1[1]
        if a == 0:
            return 90.0
        return math.atan(b/a) * 180.0 / 3.14

def main():
    path = '/Users/yuchao/Documents/memo/20161130_magnetogram/MEMO/7TH_EXPERIMENT_FLAT_PINHOLES&FOCUS/PINHOLES.nosync/20180714/B1/031915/CENT/F_031915558.fits'
    ref_image = fits.getdata(path)

    p = PinholesBasedAlign()
    pinholes = p.count_pinholes(ref_image)
    grid = p._map_to_grid(pinholes, 0.5)

if __name__ == '__main__':
    main()

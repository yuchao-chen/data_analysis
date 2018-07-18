# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
"""
# Return threshold value based on Otsu's method
# references:
# .. [1] Wikipedia, http://en.wikipedia.org/wiki/Otsu's_Method
"""
def otsu_threshold(data, nbins = 256, debug = False):
    if np.min(data) == np.max(data):
        return -1
    hist, bin_centers = np.histogram(data.ravel(), nbins)
    bin_centers = bin_centers[:-1]

    hist = hist.astype(float)
    # probabilities of all possible thresholds
    weight0 = np.cumsum(hist)
    weight1 = np.cumsum(hist[::-1])[::-1]
    # means for all possible thresholds
    mean0 = np.cumsum(hist * bin_centers) / weight0
    mean1 = (np.cumsum((hist * bin_centers)[::-1]) / weight1[::-1])[::-1]

    variance12 = weight0[:-1] * weight1[1:] * (mean0[:-1] - mean1[1:]) ** 2
    idx = np.argmax(variance12)
    threshold = bin_centers[:-1][idx]

    if debug:
        plt.figure()
        plt.plot(np.log(hist))
        plt.axvline(x = idx)
        plt.title('Histogram in \'Ostu Threshold\'')
    return threshold

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Auxiliary functions for processing and analyzing multiplexed images from project:
# https://biifsweden.github.io/projects/2023/03/15/CeciliaLindskog2023-1/

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import sklearn

# pre-processing filters
from skimage.restoration import rolling_ball
from skimage.filters import median
from skimage.morphology import disk

from skimage import io, filters, measure, segmentation, color, util, exposure, morphology
from skimage.filters import threshold_otsu, threshold_multiotsu

def preprocess(img):
    background = rolling_ball(img, radius=50)
    filtered_backgd = img - background
    filtered_median = median(filtered_backgd, disk(3))
    
    return filtered_backgd, filtered_median

def nonzero_intensity_mean(mask: np.ndarray, img: np.ndarray) -> float:
    data = img[mask]
    data = data[data != 0]
    
    if data.size != 0:
        return np.mean(data)
    else:
        return 0
    
def get_avg_intensity(ref_img, labels, cols, properties):
    images = [ref_img[(x),:,:] for x in range(ref_img.shape[0])]

    tables = [measure.regionprops_table(labels, image, properties=properties) for image in images]
    mean_intens = export_table_to_dataframe(tables, cols)
    
    return mean_intens
    
def export_table_to_dataframe(tables, cols):
    tables = [pd.DataFrame(table) for table in tables] # create dataframe for each table
    tables = [table.set_index('label') for table in tables] # reset segmentation label as table index
    
    mean_intens = pd.concat(tables, axis=1)
    mean_intens.columns = cols
    
    return mean_intens

def opal_quantification(ref_img, labels, bin_mask, cols, filter_area, filter_size, preproc):

    images = [ref_img[(x),:,:] for x in range(ref_img.shape[0])] # for each channel
    thresh_images = []
    intens_masks = []
    tables = []

    excl_cols = ['DAPI', 'Autofluorescence']

    i = 0
    for img in images:
        print('channel: ' + cols[i])

        if cols[i] in excl_cols:
            intensity_means = measure.regionprops_table(labels, img, properties=['label', 'intensity_mean'])
            tables.append(intensity_means)
        else:
            if preproc:
                print('processing...')
                # pre-processing
                background, filtered = preprocess(img)
            else:
                filtered = img

            filteredByCellMask = bin_mask * filtered # for each channel, exclude regions that are outside of the cellular region defined by the segmentation segmentation
            
            thresholds = threshold_multiotsu(filtered)
            thr = thresholds[1]
            thresholded = (img >= thr)

            if filter_area:
                thresholded = morphology.remove_small_objects(thresholded, min_size=filter_size, connectivity=2)
            thresholded = thresholded * 1

            final_mask = filteredByCellMask * thresholded

            thresh_images.append(thresholded)
            intens_masks.append(final_mask)

            nonzero_intensity_means = measure.regionprops_table(labels, final_mask, 
                properties=['label'], extra_properties=[nonzero_intensity_mean])
            tables.append(nonzero_intensity_means)
        i = i + 1

    mean_intens = export_table_to_dataframe(tables, cols)

    return mean_intens, thresh_images, intens_masks

def save_results_opal_quantification(cols, outpath, thresh_images, intens_masks):
    thresh = [io.imsave(outpath + '_thresh_mask_' + cols[x] + '.tif',util.img_as_ubyte(thresh_images[x]*255))
        for x in range(len(thresh_images))] # for each channel

    intens = [io.imsave(outpath + '_combined_mask_' + cols[x] + '.tif',intens_masks[x])
        for x in range(len(intens_masks))] # for each channel

def filter_columns(cols, df):
    newdf = df.copy()
    newdf = newdf.drop(columns=cols, axis=1)
    return newdf

def get_hist_pos_signal(mean_intens_thres_OPAL, mean_intens_thres_OPAL_nonzero):

    # counting of cols with values > 0 per row
    OPAL_counts = mean_intens_thres_OPAL_nonzero.copy()
    OPAL_counts['total_positive'] = mean_intens_thres_OPAL_nonzero[mean_intens_thres_OPAL_nonzero > 0].count(axis='columns')

    # histogram of the distribution of the postive signals
    #n, bins, patches = np.histogram(OPAL_counts['total_positive'], bins=[1,2,3,4,5,6])
    n, bins = np.histogram(OPAL_counts['total_positive'], bins=[1,2,3,4,5,6])

    # transform the data to count unique patterns
    mean_intens_unit = mean_intens_thres_OPAL.copy()
    mean_intens_unit = mean_intens_unit.loc[~(mean_intens_unit==0).all(axis=1)]
    mean_intens_unit[mean_intens_unit > 0] = 1

    signal_stats = mean_intens_unit.groupby(mean_intens_unit.columns.tolist(),as_index=False).size()

    return n, bins, signal_stats

def plot_clusters(fit_pca, labels, centers, plotcenters):
    # Show the segmentations.
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(12, 5),
        sharex=True,
        sharey=True,
    )

    axes[0].scatter(fit_pca[:, 0], fit_pca[:, 1], c=labels, s=50, cmap='Accent')
    axes[0].set_title("1st and 2nd")
    
    axes[1].scatter(fit_pca[:, 0], fit_pca[:, 2], c=labels, s=50, cmap='Accent')
    axes[1].set_title("1st and 3rd")

    axes[2].scatter(fit_pca[:, 1], fit_pca[:, 2], c=labels, s=50, cmap='Accent')
    axes[2].set_title("2nd and 3rd")
    
    if plotcenters:
        axes[0].scatter(centers[:, 0], centers[:, 1], c='black', s=200, alpha=0.5);
        axes[1].scatter(centers[:, 0], centers[:, 2], c='black', s=200, alpha=0.5);
        axes[2].scatter(centers[:, 1], centers[:, 2], c='black', s=200, alpha=0.5);

    fig.tight_layout()
    plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Auxiliary functions for processing and analyzing multiplexed images from project:
# https://biifsweden.github.io/projects/2023/03/15/CeciliaLindskog2023-1/

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sklearn

# pre-processing filters
from skimage.restoration import rolling_ball
from skimage.filters import median, threshold_otsu, threshold_multiotsu
from skimage.morphology import remove_small_objects, disk
from skimage.measure import regionprops_table
from skimage.util import img_as_ubyte
from skimage.io import imread, imsave

from stardist.models import StarDist2D, Config2D
from csbdeep.utils import normalize

from tiler import Tiler, Merger 

def segment_nuclei(img, split, method, param_stardist, param_cellpose):
    if split:
        image16 = img.astype('uint16')
        full_image = np.reshape(image16, [1, img.shape[0], img.shape[1]]) # reshape to [#channels, width, height]

        blockSize = 1000
        # Setup tiling parameters
        tiler = Tiler(data_shape=full_image.shape,tile_shape=(1, blockSize, blockSize),overlap=0,channel_dimension=0)
        merger = Merger(tiler)

        countLabels = 0

        # Process the image tile-by-tile
        for tile_id, tile in tiler(full_image):
            print(f"segmenting tile {tile_id}...")
            # run inference

            tile_new = tile[0,:,:]

            labels = segment_with_stardist(normalize(tile_new),param_stardist)

            maxLabelTile = np.amax(labels) # get max label of the current mask
            labels = np.where(labels > 0, labels+countLabels, labels) # sum countLabels to the current mask
            countLabels = countLabels + maxLabelTile # sum current max label to countLabels

            tile_res = np.reshape(labels, [1, labels.shape[0], labels.shape[1]])
            merger.add(tile_id, tile_res)

        # Merge the image
        merged_image = merger.merge(unpad=True)
        
        labels = stitchSegmentedTiles(merged_image, blockSize)
        labels = labels[:,:,0]
    else:
        print('segmenting...')
        labels = segment_with_stardist(normalize(img),param_stardist)

    return labels

def segment_with_stardist(image, param):
    model = StarDist2D.from_pretrained('2D_versatile_fluo') # load pretrained model
    labels, _ = model.predict_instances(normalize(image),nms_thresh=0.8, prob_thresh=0.7) # get predictions for nuclei
    return labels

def stitchSegmentedTiles(merged_image, blockSize):

    labeled_image = merged_image.transpose(1, 2, 0).astype('uint16')

    delta = 3

    h = labeled_image.shape[0]
    w = labeled_image.shape[1]

    i = blockSize - delta # initialize counter

    labels_to_replace = []

    # check labels along rows
    while i < w:
        rows = labeled_image[i:(i+delta*2),]
        
        for j in range(rows.shape[1]):
            labels = rows[:,j,0]
            labels = labels[labels > 0] # exclude background
            
            if(labels.size > 0):
                if max(labels) != min(labels):
                    labels_to_replace.append((min(labels), max(labels)))
            
        i = i + blockSize

        
    # check labels along columns
    i = blockSize - delta # initialize counter
    while i < h:
        cols = labeled_image[:,i:(i+delta*2),:]
        
        for j in range(cols.shape[0]):
            labels = cols[j,:,0]
            labels = labels[labels > 0] # exclude background
            
            if(labels.size > 0):
                if max(labels) != min(labels):
                    labels_to_replace.append((min(labels), max(labels)))
            
        i = i + blockSize

    labels_to_replace = list(set(labels_to_replace))

    for k in range(len(labels_to_replace)):
        labeled_image[labeled_image == labels_to_replace[k][0]] = labels_to_replace[k][1]
        
    return labeled_image

def preprocess(img):
    filtered_median = median(img, disk(3))
    
    return filtered_median

def nonzero_intensity_mean(mask: np.ndarray, img: np.ndarray) -> float:
    data = img[mask]
    data = data[data != 0]
    
    if data.size != 0:
        return np.mean(data)
    else:
        return 0
    
def get_avg_intensity(ref_img, labels, cols, properties):
    images = [ref_img[(x),:,:] for x in range(ref_img.shape[0])]

    tables = [regionprops_table(labels, image, properties=properties) for image in images]
    mean_intens = export_table_to_dataframe(tables, cols)
    
    return mean_intens
    
def export_table_to_dataframe(tables, cols):
    tables = [pd.DataFrame(table) for table in tables] # create dataframe for each table
    tables = [table.set_index('label') for table in tables] # reset segmentation label as table index
    
    mean_intens = pd.concat(tables, axis=1)
    mean_intens.columns = cols
    
    return mean_intens

def opal_quantification(ref_img, labels, bin_mask, ilastik_mask, cols, filter_area, minSize, maxSize, preproc, multi_otsu_levels):

    images = [ref_img[(x),:,:] for x in range(ref_img.shape[0])] # for each channel
    thresh_images = []
    intens_masks = []
    tables = []

    excl_cols = ['DAPI', 'Autofluorescence']

    i = 0
    for img in images:
        print('channel: ' + cols[i])

        if cols[i] in excl_cols:
            intensity_means = regionprops_table(labels, img, properties=['label', 'intensity_mean'])
            tables.append(intensity_means)

        else:
            if cols[i] == 'OPAL520':
                # load Ilastik mask
                label520 = imread(ilastik_mask)
                thresholded = (label520 == 1) # * img 
                filtered = img
            else:
                if preproc:
                    print('preprocessing...')
                    # pre-processing
                    background, filtered = preprocess(img)
                else:
                    filtered = img

                    level = multi_otsu_levels[i]
                    thresholds = threshold_multiotsu(filtered,classes=level)
                    
                    j = level - 2
                    thr = thresholds[j]
                    thresholded = (img >= thr)

            filteredByCellMask = bin_mask * filtered # for each channel, exclude regions that are outside of the cellular region defined by the segmentation segmentation
            
            if filter_area:
                #thresholded = morphology.remove_small_objects(thresholded, min_size=filter_size, connectivity=2)
                # filter by size - remove small objects
                thresholded_small = remove_small_objects(thresholded, min_size=minSize[i], connectivity=2)
                
                # filter by size - remove big objects
                thresholded_mid = remove_small_objects(thresholded_small, maxSize[i])
                thresholded = thresholded_small ^ thresholded_mid

            thresholded = thresholded * 1

            final_mask = filteredByCellMask * thresholded

            thresh_images.append(thresholded)
            intens_masks.append(final_mask)

            nonzero_intensity_means = regionprops_table(labels, final_mask, 
                properties=['label'], extra_properties=[nonzero_intensity_mean])
            tables.append(nonzero_intensity_means)
        i = i + 1

    mean_intens = export_table_to_dataframe(tables, cols)

    return mean_intens, thresh_images, intens_masks

def save_results_opal_quantification(cols, outpath, thresh_images, intens_masks):
    thresh = [imsave(outpath + '_thresh_mask_' + cols[x] + '.tif',img_as_ubyte(thresh_images[x]*255))
        for x in range(len(thresh_images))] # for each channel

    intens = [imsave(outpath + '_combined_mask_' + cols[x] + '.tif',intens_masks[x])
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
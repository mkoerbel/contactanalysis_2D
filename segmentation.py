''' 
Part of contactanalysis scripts: segmentation and tiff output

@author: Markus Koerbel
@email: mk988@cam.ac.uk

Input:  -folder: list with root directory as first element, and pandas dataframe with all files to be analysed
        -parameter: user input parameters for segmentation


'''

import numpy as np
import pandas as pd
import os.path
from scipy import ndimage, optimize
from skimage import filters, morphology, feature, io, measure, segmentation
from skimage.io._plugins import tifffile_plugin as tf


def remove_lin_background(image):
    """
    Removes a linearly increasing background from an image stack. Estimtates a correction factor c by calculating C(t,x,y) = I(t,x,y)/t and taking the median of all C as c.
    The output image is then Io = Ii - c*t
    """
    c = np.zeros(image.shape)
    for t in range(image.shape[0]):
        c[t,:,:] = image[t,:,:]/(t+1) 
    c_final = np.median(c)
    output = np.zeros(image.shape)
    for t in range(image.shape[0]):
        output[t,:,:] = image[t,:,:] - t*c_final
    return output

def remove_1storder_background(image):
    """
    Removes a background increasing acording to pseudo first order binding kinetics of the form F + S <=> Fb. S ist the free binding sites on the SLB, F the fluorophore at constant concentration. 
    This gives the integrated rate equation: Fb = S0*(1-exp(-kt))
    """

    def kinetics(t, S0, k):
        return S0*(1-np.exp(-k*t))

    I_bg = np.median(image, axis = (1,2))
    popt, _ = optimize.curve_fit(kinetics, np.linspace(0,len(I_bg)-1, len(I_bg)), I_bg, p0 = [10, 0.01], bounds = (0, [2**16, 1]))
    S0, k = popt
    print(S0, k)
    output = np.zeros(image.shape)
    for t in range(image.shape[0]):
        output[t,:,:] = image[t,:,:] - kinetics(t, S0, k)
    return output

def DoG_segmentation(image, s_low, s_high, thresh, min_area, remove_lin_bg):
    '''
    input image is 3 dimensional [t, x, y]. Bandpass filter applied in spatial (x,y) domain.
    s_low is sigma for Guassian blur of lower spatial frequency cutoff, s_high for higher (s_low > s_high).
    Use two thresholds, one for small contacts, and Otsu threshold for big, bright contacts
    '''

    if remove_lin_bg == "linear":
        image = remove_lin_background(image)
    elif remove_lin_bg == '1st_order':
        image = remove_1storder_background(image)
    
    # Detect presence of small features with DoG combined with Otsu threshold for large contacts
    image_low = filters.gaussian(image, sigma = (0,s_low,s_low), preserve_range=True)
    image_high = filters.gaussian(image, sigma = (0,s_high,s_high), preserve_range=True)
    image_bin = np.zeros(image.shape, dtype = np.bool)
    otsu_thresh = filters.threshold_otsu(image_high)
    image_bin = ((image_high - image_low) > thresh) | (image_high > otsu_thresh)
    for i_frame in range(image.shape[0]):
        image_bin[i_frame,:,:] = morphology.remove_small_objects(image_bin[i_frame,:,:], min_size = min_area)
    return image_bin

def get_cell_markers(image, frame, cell_radius, redo_seeds, i_file):
    '''
    Established the marker array to assign segmented regions to cells via Watershed. If redo_seeds = False, it will look for a previously generated markers file as input.
    If none is found, a message is displayed and a new one generated. 
    image: input image array, 3D
    frame: frame of image to be used to detect cells
    '''
    if (redo_seeds == False) & (os.path.isfile(i_file.timelapse[:-4] + '_cell_markers.tif')):
        markers = np.array(io.imread(i_file.timelapse[:-4] + '_cell_markers.tif'), dtype = np.uint8)
    elif redo_seeds == False:
        print('       - !!! No cell marker file found to load. Start cell detection ...')
        redo_seeds = True
    
    if redo_seeds == True:
        image_low = filters.gaussian(image[frame, :,:], sigma=4*cell_radius, preserve_range=True)
        image_high = filters.gaussian(image[frame, :,:], sigma=cell_radius, preserve_range=True)
        cell_pos = feature.peak_local_max(image_high - image_low, min_distance=cell_radius, indices = True, exclude_border=False, threshold_abs=0)
        markers = np.zeros(image.shape, dtype=np.int8)
        markers[frame, cell_pos[:,0], cell_pos[:,1]] = 1
        markers[frame, :,:] = morphology.dilation(markers[frame,:,:], morphology.disk(cell_radius))
        markers = ndimage.label(markers)[0]
        tf.imsave(i_file.timelapse[:-4] + '_cell_markers.tif', markers.astype(np.uint8))

    return markers

def corr_flatfield(flatfield, bitdepth, bias):
    '''
    For ratiometric flat field correction
    Input flat field image stack. Averages and normalizes flat field image stack and divides input stack. Removes overexposed images from flat field stack.
    '''
    flat_max = np.amax(flatfield, axis = (1,2))
    flat_index_del = np.array(np.where(flat_max == 2**bitdepth-1))
    flat_index = np.array([i for i in range(len(flat_max)) if i not in flat_index_del])
    corr_image = flatfield[flat_index,:,:]
    corr_image = flatfield - bias
    corr_image = np.mean(corr_image, axis = 0)
    corr_image = corr_image / np.amax(corr_image)
    return corr_image

def regions2cell_watershed(markers, regions, cell_radius):
    '''
    Takes a 3D image stack and positions for the second to last frame, then performs Watershed segmentation and returns labelles image stack of the same size.
    '''
    # use 3D watershed to segment regions into cells
    n_frames = regions.shape[0]
    distance = -ndimage.distance_transform_edt(regions)
    labels = segmentation.watershed(distance, markers, mask = regions)
    return labels

def small_regions2cell_watershed(regions_bin, labels, radius):
    '''
    Takes two input images: regions_bin a binary 3D array with segmented areas; labels a 3D array with regions labelled via a previous algorithm.
    Function iterates backwards over time (fist array coordinate) and assigns regions in 'regions_bin' that are not present in 'labels' to the closest regions labelled in 'labels' within distance 'radius'.
    Returns an array of the same size as 'labels' with the added regions labelled.
    '''
    #assign unconnected regions
    radius = np.ceil(radius).astype(np.int)
    n_frames = labels.shape[0]
    assigned = labels.astype(np.int8)
    assigned[assigned > 0] = 1
    to_assign = regions_bin - assigned
    for i_frame in range(n_frames-1, -1, -1):
        if np.sum(to_assign[i_frame,:,:]):
            regions_label = measure.label(to_assign[i_frame,:,:])
            regions_props = measure.regionprops_table(regions_label, properties = ('label', 'area', 'coords', 'centroid'))
            regions_props = pd.DataFrame(regions_props)
            for _, i_region in regions_props.iterrows():
                #define subarea to look at
                y_min = max(0, int(i_region['centroid-0']) - radius)
                y_max = min(labels.shape[0], int(i_region['centroid-0']) + radius + 1)
                x_min = max(0, int(i_region['centroid-1']) - radius)
                x_max = min(labels.shape[1], int(i_region['centroid-1']) + radius + 1) 
                t_min = max(0,i_frame-5)
                t_max = min(n_frames-1, i_frame + 5 + 1)
                i_subregion = labels[t_min:t_max, y_min:y_max, x_min:x_max]
                sub_centroid = [min(i_frame, 5), min(i_region['centroid-0'], radius), min(i_region['centroid-1'], radius)]
                #calculate distances in subregion
                calc_dist = np.array(np.where(i_subregion > 0))
                # 2 fold penalty on connections in t, i.e. t distance is counted twice #removed
                dist = np.sqrt(((calc_dist[0,:]-sub_centroid[0])) ** 2 + (calc_dist[1,:]-sub_centroid[1]) ** 2 + (calc_dist[2,:]-sub_centroid[2]) ** 2)
                # if no contacts were found, check in frame+1
                '''
                if len(dist) == 0 and i_frame < n_frames-1:
                    i_subregion = labels[i_frame+1, y_min:y_max, x_min:x_max]
                    sub_centroid = [min(i_region['centroid-0'], radius), min(i_region['centroid-1'], radius)]
                    calc_dist = np.array(np.where(i_subregion > 0))
                    dist = np.sqrt((calc_dist[0,:]-sub_centroid[0]) ** 2 + (calc_dist[1,:]-sub_centroid[1]) ** 2)
                '''
                if len(dist) > 0:
                    # if one contact at the minimal distance within radius was found
                    dist_min = np.min(dist)
                    if np.sum(dist == dist_min) == 1 and dist_min < radius:
                        region_n = i_subregion[calc_dist[0,dist == dist_min], calc_dist[1,dist == dist_min], calc_dist[2,dist == dist_min]]
                        labels[i_frame, i_region.coords[:,0], i_region.coords[:,1]] = region_n
                        #print(i_frame, i_region['centroid-0'], i_region['centroid-1'], region_n)
    return labels

def rollingball(image, size):
    '''
    Rolling ball average with a window of size frames. size has to be uneven; if even, size is set to size+1. 
    '''
    image_rb = np.zeros(image.shape)
    r = int(np.ceil((size-1)/2))
    n_frames = image.shape[0]
    for i_frame in range(n_frames):
        t_min = max(0, i_frame-r)
        t_max = min(n_frames-1, i_frame+r)
        image_rb[i_frame,:,:] = np.mean(image[t_min:t_max,:,:], axis = 0)
    return image_rb

def remove_small_zones(image, size):
    '''
    Remove objects that have a smaller size (volume) than size, including diagonal connectivity. Input image is a 3D array with segmented regions.
    '''
    image_labels = measure.label(image, connectivity=3)
    image_props = measure.regionprops(image_labels)
    for i_prop in image_props:
        if i_prop.area < size:
            i_coords = np.array(i_prop.coords)
            image[i_coords[:,0], i_coords[:,1], i_coords[:,2]] = 0
    return image

def apply_hysteresis_threshold_size(image, low, high, size):
    """Apply hysteresis thresholding to ``image``.
    Adapted from skimage package

    This algorithm finds regions where ``image`` is greater than ``high``
    OR ``image`` is greater than ``low`` *and* that region is connected to
    a region greater than ``high``.
    Parameters
    ----------
    image : array, shape (M,[ N, ..., P])
        Grayscale input image.
    low : float, or array of same shape as ``image``
        Lower threshold.
    high : float, or array of same shape as ``image``
        Higher threshold.
    Returns
    -------
    thresholded : array of bool, same shape as ``image``
        Array in which ``True`` indicates the locations where ``image``
        was above the hysteresis threshold.
    Examples
    --------
    >>> image = np.array([1, 2, 3, 2, 1, 2, 1, 3, 2])
    >>> apply_hysteresis_threshold(image, 1.5, 2.5).astype(int)
    array([0, 1, 1, 1, 0, 0, 0, 1, 1])
    References
    ----------
    .. [1] J. Canny. A computational approach to edge detection.
           IEEE Transactions on Pattern Analysis and Machine Intelligence.
           1986; vol. 8, pp.679-698.
           :DOI:`10.1109/TPAMI.1986.4767851`
    """
    low = np.clip(low, a_min=None, a_max=high)  # ensure low always below high
    mask_low = image > low
    mask_high = image > high
    # Connected components of mask_low
    labels_low, num_labels = ndimage.label(mask_low)
    mask_high = remove_small_zones(mask_high, size)
    # Check which connected components contain pixels from mask_high
    sums = ndimage.sum(mask_high, labels_low, np.arange(num_labels + 1))
    connected_to_high = sums > 0
    thresholded = connected_to_high[labels_low]
    return thresholded

#close contact analsysi function to include big contacts
def contactanalysis_big(image, labels, sigma, n_rb, thresh_std, gauss_thresh, confine, contact_min_volume):
    '''
    Uses LoG to detect small close contact zones in bilayer, combined with simple thresholding for bigger contacts. 
    Input image stack to analyse, image with labelled cell contacts, sigma for Gaussian filtering, thresh_std for hysteresis filtering of LoG image with [s_high, s_low], n_rb number of frames to average for rolling ball.
    The image Laplacian is calculated with Sobel operators. 
    '''
    n_frames = image.shape[0]
    exclusion_labels = np.zeros(image.shape, dtype=np.int16)
    accummulation_labels = np.zeros(image.shape, dtype=np.int16)
    image_div = np.zeros(image.shape)                 
    image_gauss = np.zeros(image.shape)
    
    image_rb = rollingball(image, n_rb)
    for i_frame in range(n_frames):
        image_gauss[i_frame,:,:] = filters.gaussian(image_rb[i_frame,:,:,], sigma, preserve_range=True)
        image_div[i_frame,:,:] = filters.laplace(image_gauss[i_frame,:,:])
    div_outside_mean = 0 #np.mean(image_div[labels == 0])
    div_outside_std = np.std(image_div[labels == 0])
    thresh_high = div_outside_mean + thresh_std[0]*div_outside_std
    thresh_low = div_outside_mean + thresh_std[1]*div_outside_std
    #print(thresh_high, thresh_low, np.mean(image_div[labels == 0]), outside_std)
    gauss_outside_mean = np.mean(image_gauss[labels == 0])
    gauss_outside_std = np.std(image_gauss[labels == 0])
    snr = gauss_outside_mean/gauss_outside_std
    print('         SNR of Gaussian filtered bilayer: {}'.format(snr))

    def thresh_label(thresh_sign):
        image_thresh = filters.apply_hysteresis_threshold(thresh_sign*image_div, thresh_low, thresh_high)
        image_thresh = image_thresh | exclusion_big
        image_distance = -ndimage.distance_transform_edt(image_thresh)
        image_markers = np.multiply(labels, image_thresh)
        image_markers[0,:,:] = np.amax(labels) + 10 #to remove anything static present in frame 1
        image_labels = segmentation.watershed(image_distance, markers = image_markers, mask = image_thresh)
        image_labels[image_labels == np.amax(labels)+10] = 0
        return image_labels
    
    # calculate simple threshold
    image_enh = image_gauss + image_div
    exclusion_big = image_enh < gauss_outside_mean - gauss_thresh * gauss_outside_std
    exclusion_big = morphology.erosion(exclusion_big, morphology.ball(np.ceil(sigma)))
    accummulation_big = image_enh > gauss_outside_mean + gauss_thresh * gauss_outside_std
    accummulation_big = morphology.erosion(accummulation_big, morphology.ball(np.ceil(sigma)))

    if confine:
        exclusion_labels = filters.apply_hysteresis_threshold(-image_div, thresh_low, thresh_high)
        exclusion_labels = np.maximum(exclusion_labels, exclusion_big)
        exclusion_labels = np.multiply(exclusion_labels, labels)
        accummulation_labels = filters.apply_hysteresis_threshold(image_div, thresh_low, thresh_high)
        accummulation_labels = np.maximum(accummulation_labels, accummulation_big)
        accummulation_labels = np.multiply(accummulation_labels, labels)
    else:
        exclusion_labels = thresh_label(-1)
        #exclusion_labels = np.maximum(exclusion_labels, exclusion_big) ### I think here is mistake
        accummulation_labels = thresh_label(1)
        #accummulation_labels = np.maximum(accummulation_labels, accummulation_big)

    exclusion_labels = remove_small_zones(exclusion_labels, contact_min_volume)
    accummulation_labels = remove_small_zones(accummulation_labels, contact_min_volume)

    return exclusion_labels, accummulation_labels, image_div

def CA2D_segmentation(i_file, channel, bias, bitdepth, cell_radius, psf_radius, cell_thresh, cell_min_area, LoG_sigma, n_rb, contact_thresh_high, contact_thresh_low, gauss_thresh, confine_CCZ, contact_min_volume, redo_seeds, remove_lin_bg):
    # load tif file
    print('    - Starting image segmentation')
    image_raw = np.array(io.imread(i_file.timelapse), dtype = np.int32)
    [n_frames, n_X, n_Y, n_channels] = np.array(image_raw.shape)
    image_raw = image_raw - bias
    image_raw[image_raw < 0] = 1E-6 #set values below 0 to 1E-6, as they results from read noise and do not contain signal
    # Segment Contact Zones (CZ)
    print('       - Detecting Contact Zones to define cells') 
    image_CZ_bin = DoG_segmentation(image_raw[:,:,:,channel['cellmask']], 0.5*cell_radius, psf_radius, cell_thresh, cell_min_area, remove_lin_bg)    
    # Find cell positions for Watershed seeds
    image_CZ_markers = get_cell_markers(image_raw[:,:,:,channel['cellmask']], n_frames-1, cell_radius, redo_seeds, i_file)
    # Cell segmentation using Watershed 
    image_CZ_labels = regions2cell_watershed(image_CZ_markers, image_CZ_bin, cell_radius)
    image_CZ_labels = small_regions2cell_watershed(image_CZ_bin, image_CZ_labels, 2*cell_radius)
    # save CZ image stacks
    print('       - Saving Contact Zone images') 
    #tf.imsave(i_file.timelapse[:-4] + '_CZ_bin.tif', image_CZ_bin.astype(np.uint8))
    tf.imsave(i_file.timelapse[:-4] + '_CZ.tif', image_CZ_labels.astype(np.uint8))
    
    # Analyse Close Contact Zones (CCZ)
    print('       - Detecting Close Contacts Zones in bilayer')
    # Flat field correction if ilumination file provided
    if i_file.illumination == '':
        print('         NO ILLUMINATION FILE FOR FLAT-FIELD CORRECTION PROVIDED!')
        flatfield_CZ = np.ones((n_Y,n_X))
    else:
        image_flatfield = io.imread(i_file.illumination)
        flatfield_CZ = corr_flatfield(image_flatfield[:,:,:,channel['bilayer']], bitdepth, bias)
    # Detect Close Contact Zones
    image_CCZ_corr = image_raw[:,:,:,channel['bilayer']]/flatfield_CZ
    image_CCZ_exclusion, image_CCZ_accummulation, image_CCZ_filt = contactanalysis_big(image_CCZ_corr, image_CZ_labels, LoG_sigma, n_rb, [contact_thresh_high, contact_thresh_low], gauss_thresh, confine_CCZ, contact_min_volume)
    # save CCZ images
    print('       - Saving Close Contacts Zones detected')
    tf.imsave(i_file.timelapse[:-4] + '_SLB_flatfield-corrected.tif', image_CCZ_corr.astype(np.uint32))
    tf.imsave(i_file.timelapse[:-4] + '_CCZ_exc.tif', image_CCZ_exclusion.astype(np.uint8))
    #tf.imsave(i_file.timelapse[:-4] + '_CCZ_acc.tif', image_CCZ_accummulation.astype(np.uint8))
    tf.imsave(i_file.timelapse[:-4] + '_CCZ_filt.tif', image_CCZ_filt.astype(np.uint32))

    return image_raw, image_CZ_labels, image_CCZ_exclusion, image_CCZ_accummulation, image_CCZ_corr
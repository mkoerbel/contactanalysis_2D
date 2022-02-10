'''
Part of contactanalysis scripts: analysis and feature detection

@author: Markus Koerbel
@email: mk988@cam.ac.uk

'''

import numpy as np 
import pandas as pd
from scipy import ndimage
from scipy import signal
from skimage import morphology
from skimage import measure
from skimage.io._plugins import tifffile_plugin as tf


def regions2tracks(regions, smoothing = True):
    '''
    Get cell position track from labelled areas.
    Input image stack with labelled regions. Returns a dataframe with tracks of the centroid for each region. Missing positions are interpolated. 
    If smoothing = True, a median filter is applied to the coordinates to avoid jumps. 
    '''
    region_labels = np.unique(regions)
    tracks_df = pd.DataFrame() #('frame', 'particle', 'centroid-x', 'centroid-y')
    for i_label in region_labels:
        if i_label == 0:
                continue
        i_coords = np.array(np.where(regions == i_label)) # [[frame, ...],[y,...],[x,...]]
        frames = np.unique(i_coords[0,:])
        tracks = np.zeros((4,frames[-1]-frames[0]+1))
        tracks[0,:] = frames[0] + np.array(range(frames[-1]-frames[0]+1))
        tracks[3,:] = i_label
        tracks[1:3,0] = np.sum(i_coords[1:3,i_coords[0,:] == frames[0]], axis = 1)/np.sum(i_coords[0,:] == frames[0])
        next_track = 1 # counter to go through tracks (counts from frame[0] to frame[-1])
        for i_ctr in range(len(frames)-1): # counter to go through frames (if timepoints are missing these are different)
            frame_diff = frames[i_ctr+1]-frames[i_ctr]
            # calculate difference vector between centroids at frame[i_ctr] and frame[i_ctr+1]
            diff = np.sum(i_coords[1:3, i_coords[0,:] == frames[i_ctr+1]], axis = 1)/np.sum(i_coords[0,:] == frames[i_ctr+1]) - tracks[1:3,i_ctr]
            # interpolate missing timepoints and arrange into array
            tracks[1:3,next_track:next_track+frame_diff] = np.transpose((diff/frame_diff) * np.transpose(np.array([range(frame_diff), range(frame_diff)])+1) + tracks[1:3,i_ctr])
            next_track += frame_diff
        tracks_df = tracks_df.append(pd.DataFrame(data = tracks.T, columns = ['frame','y','x','particle']))

    return tracks_df

def signal_around_point2D(image, point, radius):
    # calculate the mean intensity of the pixel values around point within radius 
    # image is 2D 

    if len(image.shape) > 2:
        raise Exception('Image handed to signal_around_point has too many dimensions')

    point = point.astype(np.int)
    y_min = max(0, point[0] - radius)
    y_max = min(image.shape[0], point[0] + radius + 1)
    x_min = max(0, point[1] - radius)
    x_max = min(image.shape[1], point[1] + radius + 1) 
    image_masked = image[y_min:y_max, x_min:x_max]

    mask = morphology.disk(radius)
    #m0, m1 = np.array(mask.shape)
    ym_min = max(0, radius - point[0])
    ym_max = min(mask.shape[0], mask.shape[0] + image.shape[0] - point[0] - radius - 1)
    xm_min = max(0, radius - point[1])
    xm_max = min(mask.shape[1], mask.shape[1] + image.shape[1] - point[1] - radius - 1)
    mask = mask[ym_min:ym_max, xm_min:xm_max]

    return np.mean(image_masked[mask==1])

def analyse_calcium(trace, min_gradient, smoothing, time_interval, min_height, min_width):
    '''
    Input: Dataframe with 'time' column and 'signal' column to be analysed. The same approach as in Jane's calcium code is taken.
    output: time of first Ca spike. 
    '''

    #18/6 smooth_trace = ndimage.gaussian_filter1d(trace['signal'], smoothing/time_interval)
    smooth_trace = signal.wiener(trace['signal'], smoothing//time_interval)
    sorted_trace = np.sort(smooth_trace)
    baseline_est = np.mean(sorted_trace[1:np.int(len(sorted_trace)/4)])
    normTrace = smooth_trace/baseline_est
    gradient = np.gradient(normTrace)
    peaks, properties = signal.find_peaks(gradient, height = min_gradient, prominence = 0)
    spikeStart = peaks #18/6 peaks[properties['prominences'] > 0.75*properties['peak_heights']]
    if spikeStart.size == 0:
        return -1, 0, 0, 0 
    else:
        # simple version without peak control
        #peak_pos = np.amin(spikeStart)
        #return trace['time'].iloc[peak_pos]

        nSpikes = len(spikeStart)
        overlapping = np.zeros(nSpikes, dtype = np.bool)
        spikeEnd = np.zeros(nSpikes)
        spikeHeight = np.zeros(nSpikes)
        spikeWidth = np.zeros(nSpikes)
        for i_spike in range(len(spikeStart)):
            heightAtStart = normTrace[spikeStart[i_spike]]
            if np.any(spikeStart[i_spike] < spikeEnd):
                overlapping[i_spike] = True
            afterSpikeStart = trace['time'] > trace['time'].iloc[spikeStart[i_spike]]
            belowStartHeight = normTrace < heightAtStart
            findEnd = np.where(afterSpikeStart & belowStartHeight)[0]
            if np.any(findEnd):
                spikeEnd[i_spike] = findEnd[0]
            else:
                spikeEnd[i_spike] = len(trace['time']) - 1
            spikeHeight[i_spike] = np.amax(normTrace[int(spikeStart[i_spike]):int(spikeEnd[i_spike])])
            spikeWidth[i_spike] = (spikeEnd[i_spike] - spikeStart[i_spike]) * time_interval
        spikeRawHeight = smooth_trace[spikeStart] > baseline_est #18/6 np.amax(smooth_trace)*0.25
        selSpikes = (spikeHeight > min_height) & (spikeWidth > min_width) & (spikeRawHeight) #18/6 & (~overlapping)
        #print(trace['time'].iloc[spikeStart].values, spikeHeight, spikeWidth, overlapping, trace['time'].iloc[peaks].values, properties['peak_heights'])
        if np.any(selSpikes):
            retSpike = np.where((spikeHeight == np.max(spikeHeight[selSpikes])) & selSpikes)[0]
            return trace['time'].iloc[np.amin(spikeStart[retSpike])], spikeWidth[np.min(np.where(selSpikes))], spikeHeight[np.min(np.where(selSpikes))], nSpikes
            
            #18/6 return trace['time'].iloc[np.amin(spikeStart[selSpikes])], spikeWidth[np.min(np.where(selSpikes))], spikeHeight[np.min(np.where(selSpikes))], nSpikes
        else:
            return -1, 0, 0, nSpikes

def analyse_speed(x_coords, y_coords, smoothing, time_interval):
    '''
    Input: coordinates of cell trace in pixel space, smoothing parameter as sigma for Gaussian smoothing of coordinates in s
        time_interval the ttime bewtween frames in s
    Output: The cellular speed per frame
    '''
    x_smooth = ndimage.gaussian_filter1d(x_coords, smoothing/time_interval)
    y_smooth = ndimage.gaussian_filter1d(y_coords, smoothing/time_interval)
    x_diff = np.diff(x_smooth)
    y_diff = np.diff(y_smooth)
    displacement = np.sqrt(x_diff**2 + y_diff**2)
    return displacement

def analyse_contacts(contact_results, i_file, i_cell, i_frame, i_CCZ_start, i_CCZ_end, i_image_CZ_labels, i_image_CCZ_labels, image_CCZ_corr, image_CCZ_times, image_CCZ_center, time_interval, time_Ca, frame_adhesion, pixel_size, QC):
    '''
    Input binary images, outputs contacts zone features. Input is just one frame
    '''
    [n_frames, n_X, n_Y] = np.shape(i_image_CZ_labels)
    i_CZ_coords = np.array(np.where(i_image_CZ_labels[i_frame,:,:] == 1))
    y_min = max(np.min(i_CZ_coords[0]) - 1, 0)
    y_max = min(np.max(i_CZ_coords[0]) + 1, n_X-1)
    x_min = max(np.min(i_CZ_coords[1]) - 1, 0)
    x_max = min(np.max(i_CZ_coords[1]) + 1, n_Y-1)
    #i_CZ_image = np.zeros((np.max(i_CZ_coords[0]) - np.min(i_CZ_coords[0]) + 1, np.max(i_CZ_coords[1]) - np.min(i_CZ_coords[1]) + 1))
    #i_CZ_image[i_CZ_coords[0] - np.min(i_CZ_coords[0]), i_CZ_coords[1] - np.min(i_CZ_coords[1])] = 1
    i_contact_CZ = measure.label(i_image_CZ_labels[i_frame,y_min:y_max, x_min:x_max])
    i_contact_CZ_props = measure.regionprops_table(i_contact_CZ, properties = ('area', 'major_axis_length', 'minor_axis_length', 'coords'))
    i_contact_CZ_props = pd.DataFrame(i_contact_CZ_props)
    # add contact results
    for _, i_CZ in i_contact_CZ_props.iterrows():
        i_CZ_center = np.mean(i_CZ.coords, axis = 0)
        contact_results = contact_results.append({
            'file': i_file.folder, 
            'cell': i_cell,
            'contact': 'CZ',
            'frame': i_frame,
            'time [s]': i_frame*time_interval,
            'time_to_Ca [s]': i_frame*time_interval - time_Ca,
            'time_to_adhesion [s]': (i_frame - frame_adhesion)*time_interval,
            'contact_time [s]': -1,
            'area [um2]': i_CZ.area * pixel_size**2,
            'exclusion': 0,
            'time_evol': 'NA',
            'QC': QC,
            'major_length [um]': i_CZ.major_axis_length,
            'minor_length [um]': i_CZ.minor_axis_length, 
            'x-coord [um]': (i_CZ_center[1] + x_min) * pixel_size,  # absolute position in whole image frame. 
            'y-coord [um]': (i_CZ_center[0] + y_min) * pixel_size, 
            'x-displ [um]': 0,
            'y-displ [um]': 0}, ignore_index=True)
    i_CZ_outline = morphology.dilation(i_image_CZ_labels[i_frame,y_min:y_max, x_min:x_max], morphology.disk(1)) - i_image_CZ_labels[i_frame,y_min:y_max, x_min:x_max]
    i_image_CCZ_corr = image_CCZ_corr[i_frame, y_min:y_max, x_min:x_max]
    i_SLB_baseline = np.mean(i_image_CCZ_corr[i_CZ_outline > 0])
    i_contact_CCZ = measure.label(i_image_CCZ_labels[i_frame,y_min:y_max, x_min:x_max])
    i_contact_CCZ_props = measure.regionprops_table(i_contact_CCZ, intensity_image=i_image_CCZ_corr, properties = ('area', 'major_axis_length', 'minor_axis_length', 'coords', 'min_intensity'))
    i_contact_CCZ_props = pd.DataFrame(i_contact_CCZ_props)
    for _,i_CCZ in i_contact_CCZ_props.iterrows():
        exclusion_min = 1 - i_CCZ.min_intensity/i_SLB_baseline
        int_values = i_image_CCZ_corr[i_CCZ.coords[:,0], i_CCZ.coords[:,1]]
        exclusion_10 = 1-np.percentile(int_values, 10)/i_SLB_baseline
        exclusion_mean = 1-np.mean(int_values)/i_SLB_baseline
        # get contact times and fill image_CCZ_times
        if i_frame > 0:
            temp_contact_times = np.unique(image_CCZ_times[i_frame-1, i_CCZ.coords[:,0]+y_min, i_CCZ.coords[:,1]+x_min])
            temp_contact_times = np.max(temp_contact_times)
            x_displ = 0
            y_displ = 0
            i_CCZ_center = np.mean(i_CCZ.coords, axis = 0)
            image_CCZ_center[0,i_frame, i_CCZ.coords[:,0]+y_min, i_CCZ.coords[:,1]+x_min] = i_CCZ_center[0] + y_min # y
            image_CCZ_center[1,i_frame, i_CCZ.coords[:,0]+y_min, i_CCZ.coords[:,1]+x_min] = i_CCZ_center[1] + x_min # x
            if temp_contact_times == 0:
                i_contact_time = 1
                i_CCZ_start += 1
                CCZ_type = 'start'
            else:
                i_contact_time = temp_contact_times + 1
                CCZ_type = 'continuing'
                prev_coords = np.unique(image_CCZ_center[:, i_frame-1, i_CCZ.coords[:,0]+y_min, i_CCZ.coords[:,1]+x_min], axis = 1)
                prev_x = prev_coords[1]
                prev_y = prev_coords[0]
                if (len(prev_x) > 1) & (len(prev_y) > 1) & (prev_x[0] == 0) & (prev_y[0] == 0):
                    prev_x = prev_x[1:]
                    prev_y = prev_y[1:]
                x_displ = i_CCZ_center[1] + x_min - prev_x
                y_displ = i_CCZ_center[0] + y_min - prev_y
                #x_displ = prev_coords[1,:] - i_CCZ_center[1] - x_min
                #y_displ = prev_coords[0,:] - i_CCZ_center[0] - y_min
                if len(x_displ) > 1:
                    displ2 = x_displ*x_displ + y_displ*y_displ
                    displ_pos = np.where(displ2 == np.min(displ2))[0]
                    if len(displ_pos) > 1:
                        displ_pos = displ_pos[0]
                    x_displ = x_displ[displ_pos]
                    y_displ = y_displ[displ_pos]
            if (i_frame < n_frames-1) and not np.any(i_image_CCZ_labels[i_frame + 1, i_CCZ.coords[:,0]+y_min, i_CCZ.coords[:,1]+x_min]):
                i_CCZ_end += 1
                if CCZ_type == 'start':
                    CCZ_type = 'startend'
                else:
                    CCZ_type = 'end'                
            image_CCZ_times[i_frame, i_CCZ.coords[:,0]+y_min, i_CCZ.coords[:,1]+x_min] = i_contact_time
            contact_results = contact_results.append({
                'file': i_file.folder, 
                'cell': i_cell,
                'contact': 'CCZ',
                'frame': i_frame,
                'time [s]': i_frame*time_interval,
                'time_to_Ca [s]': i_frame*time_interval - time_Ca,
                'time_to_adhesion [s]': (i_frame - frame_adhesion)*time_interval,
                'contact_time [s]': i_contact_time*time_interval,
                'area [um2]': i_CCZ.area * pixel_size**2,
                'exclusion': exclusion_10,
                'exclusion_min': exclusion_min,
                'exclusion_10': exclusion_10,
                'exclusion_mean': exclusion_mean,
                'time_evol': CCZ_type,
                'QC': QC,
                'major_length [um]': i_CCZ.major_axis_length * pixel_size,
                'minor_length [um]': i_CCZ.minor_axis_length * pixel_size,
                'x-coord [um]': (i_CCZ_center[1] + x_min) * pixel_size, 
                'y-coord [um]': (i_CCZ_center[0] + y_min) * pixel_size,
                'x-displ [um]': np.float(x_displ) * pixel_size,
                'y-displ [um]': np.float(y_displ) * pixel_size
                }, ignore_index=True)

    return contact_results, i_CCZ_start, i_CCZ_end, image_CCZ_times

def analysis(i_file, cell_summary, cell_results, contact_results, image_CZ_labels, image_CCZ_exclusion, image_CCZ_accummulation, image_raw_signal, image_CCZ_corr, parameters):
    '''
    1 loop over cell
    2 analyse Ca and adhesion 
    3 do quality control
    4 loop over frames and fill cell_results and contact_results

    '''
    print('    - Starting Analysis')
    [cell_radius, min_gradient, trace_smoothing, time_interval, min_height, min_width, pixel_size, max_speed, min_track_length, tmin_noCa] = parameters
    # track CZ regions
    cell_tracks = regions2tracks(image_CZ_labels)

    [n_frames, n_X, n_Y] = np.shape(image_CZ_labels)
    cell_labels = np.unique(image_CZ_labels)[1:]
    image_CCZ_times = np.zeros(image_CZ_labels.shape, dtype = np.int)
    image_CCZ_center = np.zeros((2,n_frames,n_X,n_Y))

    for i_cell in cell_labels:
        # get Ca
        # get adhesion
        # check quality
        # loop over frames and 
            # calculate contact results
            # infer cell results
        print('      - Analysing cell number ' + str(i_cell) + ' (' + str(np.where(cell_labels == i_cell)[0][0] + 1) + '/' + str(len(cell_labels)) + ')')
        i_cell_frames = np.unique(np.where(image_CZ_labels == i_cell)[0])
        track_length = i_cell_frames[-1] - i_cell_frames[0] + 1
        if track_length < min_track_length:
            continue
        # get signalling trace
        signal = np.zeros(track_length)
        for i_frame in range(i_cell_frames[0], i_cell_frames[-1] + 1):
            cell_center = np.array([cell_tracks[(cell_tracks['frame']==i_frame) & (cell_tracks['particle']==i_cell)]['y'].iloc[0], cell_tracks[(cell_tracks['frame']==i_frame) & (cell_tracks['particle']==i_cell)]['x'].iloc[0]])
            signal[i_frame-i_cell_frames[0]] = signal_around_point2D(image_raw_signal[i_frame,:,:], cell_center, cell_radius)
        # get time of signalling (Ca)
        trace = pd.DataFrame({'time': (np.array(range(i_cell_frames[-1]-i_cell_frames[0]+1))+i_cell_frames[0])*time_interval, 'signal': signal}) # trace starts with time 0
        # trace.to_csv(i_file.timelapse[:-4] + '_' + str(i_cell) + '_Ca-trace.csv')
        time_Ca, signal_width, signal_height, signal_cnt = analyse_calcium(trace, min_gradient, trace_smoothing, time_interval, min_height, min_width)
        # get time of adhesion
        cell_speed = analyse_speed(np.array(cell_tracks[cell_tracks['particle'] == i_cell]['x']), np.array(cell_tracks[cell_tracks['particle'] == i_cell]['y']), trace_smoothing, time_interval)*pixel_size/time_interval #um/s
        change_speed = np.diff(np.array(cell_speed < max_speed).astype(np.int)) # returns 1 if the cell speed changes below max_speed (negative slope), 0 if it stays above/below max_speed, and -1 if the cell starts to move faster than max_speed (positive slope)
        change_pos = np.where(change_speed == 1)[0]
        if len(change_pos) > 0:
            frame_adhesion = np.amax(change_pos) + i_cell_frames[0]
        elif cell_speed[0] < max_speed/pixel_size:  
            frame_adhesion = i_cell_frames[0]
        else:
            frame_adhesion = -1
        # get time of first CCZ
        CCZ_frames = np.where(image_CCZ_exclusion == i_cell)[0]
        if len(CCZ_frames) > 0:
            frame_CCZ_first = CCZ_frames[0]
        else:
            frame_CCZ_first = -1
        # check quality
        # check if cell touches image border
        i_cell_CZ_border = image_CZ_labels == i_cell
        i_cell_CZ_border[:,1:-1,1:-1] = 0
        i_cell_CZ_border_frames = np.array(np.where(i_cell_CZ_border == 1)[0])
        if i_cell_CZ_border_frames.size > 0:
            i_frame_border = i_cell_CZ_border_frames[0]
        else:
            i_frame_border = n_frames + 1
        if time_Ca > 0:
            # cell is triggering
            if (i_frame_border > time_Ca/time_interval) or ((i_frame_border < frame_CCZ_first) and (int(time_Ca/time_interval) not in i_cell_CZ_border_frames)) :
                QC = 'good'
                print('         Quality Control PASSED')
            else:  
                QC = 'border_before_Ca'
                print('         Quality Control FAILED')
        else:
            # cell is not triggering
            # time_Ca = track_length*time_interval #is in s
            if (i_frame_border > frame_CCZ_first + tmin_noCa/time_interval) and (frame_CCZ_first > 0):
                QC = 'good_noCa'
                print('         Quality Control PASSED')
            else:
                QC = 'border_before_CCZ'
                print('         Quality Control FAILED')
        cell_summary = cell_summary.append({
            'file': i_file.folder,
            'cell': i_cell,
            'QC': QC,
            'time_Ca [s]': time_Ca,
            'time_adhesion [s]': frame_adhesion*time_interval,
            'time_CCZ_first [s]': frame_CCZ_first*time_interval, 
            'time_CZ_first [s]': i_cell_frames[0]*time_interval,
            'Ca_signal_width': signal_width, 
            'Ca_signal_height': signal_height, 
            'Ca_signal_cnt': signal_cnt}, ignore_index=True)
        if (QC == 'border_before_Ca') or (QC == 'border_before_CCZ'):           
            continue

        # prepare arrays and non-framewise properties
        print('         Calculating contact features:')
        i_image_CZ_labels = (image_CZ_labels == i_cell).astype(np.uint8)
        i_image_CCZ_labels = (image_CCZ_exclusion == i_cell).astype(np.uint8)
        image_CZ_sum = np.zeros((n_X, n_Y), dtype = np.uint8)
        image_CCZ_sum = np.zeros((n_X, n_Y), dtype = np.uint8)
        i_CCZ_start = 0 #unused?
        i_CCZ_end = 0 #unused?
        for i_frame in i_cell_frames:
            print('          Frame ' + str(i_frame - i_cell_frames[0] + 1) + '/' + str(i_cell_frames[-1] - i_cell_frames[0] + 1), end = '\r')
            # get cell stats
            # get contact stats
            contact_results, i_CCZ_start, i_CCZ_end, image_CCZ_times = analyse_contacts(contact_results, i_file, i_cell, i_frame, i_CCZ_start, i_CCZ_end, i_image_CZ_labels, i_image_CCZ_labels, image_CCZ_corr, image_CCZ_times, image_CCZ_center, time_interval, time_Ca, frame_adhesion, pixel_size, QC)
            
            if i_file.illumination == '':
                SLB_mean_signal = 0
            else:
                SLB_mean_signal = np.mean(image_CCZ_corr[0,:,:])
            image_CZ_sum += i_image_CZ_labels[i_frame,:,:]
            image_CCZ_sum += i_image_CCZ_labels[i_frame,:,:]
            cell_results = cell_results.append({
            'cell':                 i_cell,
            'file':                 i_file.folder,
            'frame':                i_frame,
            'QC':                   QC,
            'time [s]':             i_frame * time_interval,
            'time_to_CZ [s]':          (i_frame - i_cell_frames[0]) * time_interval,
            'time_to_Ca [s]':       i_frame*time_interval - time_Ca,
            'time_to_adhesion [s]': (i_frame - frame_adhesion)*time_interval,
            'time_to_CCZ_first [s]': (i_frame - frame_CCZ_first) * time_interval,
            'Ca_signal [AU]':       trace[trace['time'] == i_frame * time_interval]['signal'].iloc[0],
            'cell_center_x [px]':   cell_tracks[(cell_tracks['frame']==i_frame) & (cell_tracks['particle']==i_cell)]['x'].iloc[0],
            'cell_center_y [px]':   cell_tracks[(cell_tracks['frame']==i_frame) & (cell_tracks['particle']==i_cell)]['y'].iloc[0],
            'CZ_area_total [um2]':  np.sum(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CZ')]['area [um2]']),
            'CZ_area_cnt':          len(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CZ')]['area [um2]']),
            'CZ_area_max [um2]':    np.max(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CZ')]['area [um2]']),
            'CCZ_area_total [um2]': np.sum(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['area [um2]']),
            'CCZ_area_cnt':         len(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['area [um2]']),
            'CCZ_area_mean [um2]': np.mean(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['area [um2]']),
            'CCZ_area_std [um2]':   np.std(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['area [um2]']),
            'CCZ_area_min [um2]':   np.min(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['area [um2]']),
            'CCZ_area_max [um2]':   np.max(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['area [um2]']),
            'CCZ_excl_mean': np.mean(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['exclusion']),
            'CCZ_excl_std':   np.std(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['exclusion']),
            'CCZ_excl_min':   np.min(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['exclusion']),
            'CCZ_excl_max':   np.max(contact_results[(contact_results['frame'] == i_frame) & (contact_results['cell'] == i_cell) & (contact_results['file'] == i_file.folder) & (contact_results['contact'] == 'CCZ')]['exclusion']),
            'SLB_mean_signal':      SLB_mean_signal,
            'CZ_area_sampled [um2]': np.sum(image_CZ_sum > 0) * pixel_size ** 2,
            'CCZ_area_sampled [um2]': np.sum(image_CCZ_sum > 0) * pixel_size ** 2
            }, ignore_index=True)

    tf.imsave(i_file.timelapse[:-4] + '_CCZ_contact times.tif', image_CCZ_times.astype(np.uint16))
    return cell_summary, cell_results, contact_results
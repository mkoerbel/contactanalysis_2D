"""
Contactanalysis

@author: Markus Koerbel
@email: mk988@cam.ac.uk

Temporal contact zone analysis of the interaction of cells with a supported lipid bilayer (SLB). Base on three-colour TIRF microscopy
imaging a signaling readout, contact zones and a cell membrane marker. 
Consult README for more information.


INPUT 
    - Select a folder with follwing folder structure
        . Each tif stack to be analysed needs to be in a subfolder including id_stack in its name. 
        . Corresponding to each analysis subfolder there may be a folder with a tif stack for flat-field correction of the 'Bilayer signal'-channel. This needs to have the same number of channels as the analysis tif stack, and be in a subfolder having the same name but id_stack replaced by id_flatfield. If there is no such folder, no flat-field correction will be performed.
        . In each subfolder the .tif file with the shortest name will be taken as the raw input. Generated ouytput tif stacks always have a longer file name. The code can therefore be run on the same folder-tree and the oputput tif files will be replaced. 
        . The code will loop over all subfolders containing id_stack in the name and save the results in the selected folder. 
        . Exmaple structure (id_stack = 'timelapse, id_flatfield = 'illumination):
            Condition_1 [folder]
                |-- SLB1_timelapse [folder]
                    |-- Celltype_SLB1.tif [tif stack to analyse]
                |-- SLB1_illumination [folder]
                    |-- SLB1_illumination_stack.tif [tif stack for flatfield correction]
                |-- SLB2_timelapse [folder]
                    |-- Celltype_2_SLB2.tif [tif stack to analyse - since there is no SLB2_illumination folder, no flatfield correction will be performed on this stack] 

    - Input tif files need to be multi channel images and the order of the channels defined in "channel". 
        There needs to be a channel assigned 'cellmask', 'bilayer', and 'calcium' each. E.g.
        0: CellMask (in TIRF)
        1: Bilayer signal (ICAM or CD45/CD43, in TIRF)
        2: Calcium signal (in Epifluorescence)
        If channel order is different, change manually in USER settings

    - The analysis is performed in 3 steps:
        1) segmentation 
            Run analsyis only with this option selected first and make sure the segmentation parameters are set appropriately to give good labels. If the cell assignment contains errors, the file "*_cell_markers.tif" may be edited manually and then the code rerun. 
        2) analysis
            Analyse the segmented regions, detect signaling event and output information for each cell. 
        3) plotting
            Takes the "analysis" output, restructures it for easier visualisation and adds some additional interaction features. 
      Always the latest saved output will be used in subsequent runs.
      
OUTPUT
    - various tif stacks of the filtered and segmented channels. 
    - .csv files for a time-independent cell cummary ("*_cell_summary_plot.csv"), a time-dependent cell summary ("*_cell_results_plot.csv"), and a time-dependent contact summary ("*_contact_results_plot.csv")

Last edited 19/05/2022
"""

### USER SETTINGS

### Input Settings
channel = {'cellmask': 2, 'bilayer':1, 'calcium':0}  # Dictionary. Change numbers to match channel order in tif stack. Start counting at 0
id_stack = 'timelapse' 
id_flatfield = 'illumination'
confine_CCZ = True   # set True if only contacts within segmented cell membrane areas should be considered. Otherwise in time connected contacts will be added.
redo_seeds = False   # set to True if the cell seeds should be calculated again; otherwise an already present file will be used (same filename as input tif stack but ending "_cell_markers.tif"). 

### Acquisition Parameters
bitdepth = 16                   # camera bit depth. images in flat-field stack will be excluded if they are overexposed.
bias = 400                      # camera bias. Substraced from raw data.
pixel_size = 16/150             # um ; camera pixel size
time_interval = 2               # s ; time interval between frames in tif stack

### Segmentation Parameters
cell_radius = 25
psf_radius = 1.5                # px ; used as lowpass filter cutoff
LoG_sigma = 1.5                 # px ; sigma for Gaussian blur before Lalpace filter.
n_rb = 3                        # frames ; number of frames to average with rolling ball average in time
cell_thresh = 2000              # AU ; threshold of DoG filter for CellMask segmentation
contact_thresh_low = 2.5        # 2 # x-fold of signal standard deviation outside cellmask ; Will be used as 0 + contact_thresh_low*std for the lower threshold of hysteresis thresholding of close contact zones. 
contact_thresh_high = 4         # 4 # x-fold of signal standard deviation outside cellmask ; Will be used as 0 + contact_thresh_low*std for the upper threshold of hysteresis thresholding of close contact zones.
gauss_thresh = 2.5              # x-fold of signal standard deviation outside cellmask ; Will be used as mean + gauss_thresh*std for thresholding of big close contact zones.
cell_min_area = 6               # px ; minimum size of CellMask signal per frame. 6 px for FWHM of diffraction limited microvillus (70 - 150 nm).
contact_min_volume = 12         # px ; minimum (x,y,t) volume size of detected close contact zones. E.g. for contacts with a min area of 6 px to be present in two consecutive frames set to 12
remove_lin_bg = "linear"        # CellMask can cause the background to increase over time. Choose from "None", "linear", and "1st_order" for background correction.

### Analysis Parameters
tmin_noCa = 30                  # s ; time for a cell that has not Ca triggerred to have CCZ formed to be "good" and valid for further analysis. 
max_speed = 0.15                # um/s ; threshold of cell movement to define adhered cell
min_gradient = 0.3              # /s ; threshold for intensity gradient to define Ca spike
min_height = 2.5                # min height of calcium spike relative to baseline estimate
min_width =  5                  # s ; min width of calcium spike in s
trace_smoothing = 20            # s ; gaussian filter for smoothing Ca and displacement traces
min_track_length = 10           # s ; minimum number of time a cell has to be detected (=min_track_length/time_interval frames)
QC = ['good', 'good_noCa']      # list of QC conditions that will be further investigated when running "plotting". Possible options: 'good', 'good_noCa' 
plot_smoothing = 20             # s ; smoothing windows for lineplots of timetraces






###########################################
# IMPORT MODULES
# general modules
import numpy as np 
import pandas as pd
from skimage import io
import datetime
import os

# Custom functions
from gui import CA2D_gui
from segmentation import CA2D_segmentation
from analysis import CA2D_analysis
from plotting import CA2D_plotting, CA2D_plot_Ca_traces

# MAIN

# have user input list of folder and tick which part of the program to run
[folders, programs] = CA2D_gui([id_stack, id_flatfield])
print(str(len(folders)) + ' condition(s) selected. Lets get to work!\n')

for i_folder in folders:        # one folder is one condition
    dirname = i_folder[0]       # root directory name for given condition(i_folder)
    sel_files = i_folder[1]     # pandas dataframe with paths to files

    cell_summary = pd.DataFrame()
    cell_results = pd.DataFrame()
    contact_results = pd.DataFrame()

    if programs['segmentation'] and programs['plotting'] and (not programs['analysis']):
        print('!! CANNOT DO PLOTTING WITHOUT ANALYSIS! Analysis added to selection. !!')
        programs['analysis'] = True

    for index, i_file in sel_files.iterrows():
        if programs['segmentation']:
            print(' Working on file ' + str(i_file.folder))
            image_raw, image_CZ_labels, image_CCZ_exclusion, image_CCZ_accummulation, image_CCZ_corr = CA2D_segmentation(i_file, channel, bias, bitdepth, cell_radius, psf_radius, cell_thresh, cell_min_area, LoG_sigma, n_rb, contact_thresh_high, contact_thresh_low, gauss_thresh, confine_CCZ, contact_min_volume, redo_seeds, remove_lin_bg)
            
        if programs['analysis']:
            if not programs['segmentation']:
                print(' Working on file ' + str(i_file.folder))
                #load segmentation
                print('    - Loading segmented images')
                image_raw = np.array(io.imread(i_file[id_stack]), dtype = np.int32)
                image_raw = image_raw - bias
                image_raw[image_raw < 0] = 1E-6
                image_CZ_labels = np.array(io.imread(i_file[id_stack][:-4] + '_CZ.tif'), dtype = np.uint8)
                image_CCZ_exclusion = np.array(io.imread(i_file[id_stack][:-4] + '_CCZ_exc.tif'), dtype = np.uint8)
                image_CCZ_accummulation = np.array(io.imread(i_file[id_stack][:-4] + '_CCZ_acc.tif'), dtype = np.uint8)
                image_CCZ_corr = np.array(io.imread(i_file[id_stack][:-4] + '_SLB_flatfield-corrected.tif'), dtype = np.int32)

            analysis_parameters = [cell_radius, min_gradient, trace_smoothing, time_interval, min_height, min_width, pixel_size, max_speed, min_track_length, tmin_noCa]
            cell_summary, cell_results, contact_results = CA2D_analysis(i_file, cell_summary, cell_results, contact_results, image_CZ_labels, image_CCZ_exclusion, image_CCZ_accummulation, image_raw[:,:,:,channel['calcium']], image_CCZ_corr, analysis_parameters)
            contact_results.to_csv(dirname + '/contact_results.csv')
            cell_results.to_csv(dirname + '/cell_results.csv')
            cell_summary.to_csv(dirname + '/cell_summary.csv')
            contact_results.to_csv(dirname +'/' + datetime.datetime.now().strftime("%Y%m%d") + '_contact_results.csv')
            cell_results.to_csv(dirname +'/' + datetime.datetime.now().strftime("%Y%m%d") + '_cell_results.csv')
            cell_summary.to_csv(dirname +'/' + datetime.datetime.now().strftime("%Y%m%d") + '_cell_summary.csv')
            
    if programs['plotting']: 
        condition = os.path.basename(dirname)
        print(' Plotting condition ' + str(condition))
        if not programs['analysis']:
            #load analysis results
            contact_results = pd.read_csv(dirname + '/contact_results.csv')
            cell_results = pd.read_csv(dirname + '/cell_results.csv')
            cell_summary = pd.read_csv(dirname + '/cell_summary.csv')
        cell_summary_plot = pd.DataFrame()
        cell_results_plot = pd.DataFrame()
        contact_results_plot = pd.DataFrame()
        for iq in QC:
            cell_summary_plot_t, cell_results_plot_t, contact_results_plot_t = CA2D_plotting(cell_results, cell_summary, contact_results, iq, plot_smoothing/time_interval, pixel_size)
            cell_summary_plot_t['condition'] = condition
            cell_results_plot_t['condition'] = condition
            contact_results_plot_t['condition'] = condition
            cell_summary_plot = cell_summary_plot.append(cell_summary_plot_t)
            cell_results_plot = cell_results_plot.append(cell_results_plot_t)
            contact_results_plot = contact_results_plot.append(contact_results_plot_t)
        cell_summary_plot.to_csv(dirname + '/cell_summary_plot.csv')
        cell_results_plot.to_csv(dirname + '/cell_results_plot.csv')
        contact_results_plot.to_csv(dirname + '/contact_results_plot.csv')
        cell_summary_plot.to_csv(dirname +'/' + datetime.datetime.now().strftime("%Y%m%d") + '_cell_summary_plot.csv')
        cell_results_plot.to_csv(dirname +'/' + datetime.datetime.now().strftime("%Y%m%d") + '_cell_results_plot.csv')
        contact_results_plot.to_csv(dirname +'/' + datetime.datetime.now().strftime("%Y%m%d") + '_contact_results_plot.csv')
        print('    - Generating signaling plots')
        Ca_fig = CA2D_plot_Ca_traces(cell_summary_plot, cell_results_plot, trace_smoothing//time_interval)
        Ca_fig.savefig(dirname + '/Signaling_traces.pdf', dpi = 300)

    # save parameters and date
    save_param = {
        'channel':[channel],
        'programs (segmentation - analysis - plotting': ['{} - {} - {}'.format(programs['segmentation'], programs['analysis'], programs['plotting'])],
        'id_stack': [id_stack],
        'id_flatfield': [id_flatfield],
        'confine_CCZ': confine_CCZ,
        'bitdepth':[bitdepth],
        'bias [AU]':[bias],
        'pixel_size [um]':[pixel_size],
        'time_interval [s]':[time_interval],
        'cell_radius [px]':[cell_radius],
        'thresh [AU]':[cell_thresh],
        'contact_thresh_low': [contact_thresh_low],
        'contact_thresh_high': [contact_thresh_high],
        'psf_radius [px]':[psf_radius],
        'LoG_sigma': [LoG_sigma], 
        'n_rb': [n_rb],
        'cell_min_area [px]':[cell_min_area],
        'contact_min_volume [px]':[contact_min_volume],
        'tmin_noCa [s]': [tmin_noCa],
        'max_speed [um/s]':[max_speed],
        'min_gradient [s-1]':[min_gradient],
        'min_height':[min_height],
        'min_width [s]':[min_width],
        'trace_smoothing':[trace_smoothing],
        'min_track_length':[min_track_length],
        'plot_smoothing': [plot_smoothing],
        'date_analysed':[datetime.datetime.now().strftime("%Y-%B-%d %H:%M:%S")],
        'contactgrowth_version': ["1.4"]}
    save_param = pd.DataFrame(data=save_param)
    save_param.to_csv(dirname+'/' + datetime.datetime.now().strftime("%Y%m%d") + '_contactgrowth_parameters.csv', index = False)

print('\n Contactgrowth Analysis complete. ENJOY YOUR DATA!')
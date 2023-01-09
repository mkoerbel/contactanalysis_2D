# contactanalysis_2D

Analysis of three-color total internal reflection fluorescence (TIRF) microscopy images of T cells interacting with supported lipid bilayers (SLBs). 

Written by Markus Koerbel. Email mk988@cam.ac.uk

The input TIRF images (.tif) are timelapse data with one channel capturing the T-cell membrane, one the exclusion of glycocalyx proteins on the SLB, and the third intracellular calcium of the T cell (acquired in pseudo-TIRF).

![](https://github.com/mkoerbel/contactanalysis_2D/raw/main/docs/contactanalysis_2D_overview.png)
**Overview of the image analysis** 

---

## 1. Setup

System Requirements:

The code has been tested on Windows 10 with an Intel i7-7700K CPU and 32 GB RAM. Depending on the size of the input files at least 16 GB RAM are recommended. 

Installation (on Windwos, about 10 min):

If you are familiar with Python and have a working installation, jump to step 5 or make your own virtual environment and install the dependencies listed below. 
1. Download the Python 3 installer (.exe, www.python.org, tested with Python 3.8 and 3.9. 3.11 is not recommended, installation of scikit-image fails). 
2. Install Python. Check "Add to PATH". 
3. Open a command prompt window. (To test the Python installation, type `python`. You should see the command line starting with `>>>` and you could start typing Python code. Leave this mode with `exit()`.)
4. Navigate to the folder you want the new virtual environment to be built in, e.g. `cd C:\Users\Me\Documents\Contactanalysis`. Change `Me` to your user name.
5. Create a virtual environment called `env` in the folder by running `python -m venv \env`. The `python` commmad might have to be adjusted based on your Python installation, often to `python3`. You should now have a folder `C:\Users\Me\Documents\Contactanalysis\env`. 
6. Activate your new virtual environment by typing `\env\Scripts\activate.bat` in the terminal. 
7. Once activated (If you have done steps 1-4 the command line will start with `(env) C:\Users\Me\Documents\Contactanalysis>`), install the required packages via pip (21.3.1) to your virtual environment via `pip install numpy pandas matplotlib scipy scikit-image seaborn`:

Package dependencies (tested with the versions in brackets):

- numpy (1.21.4) 
- pandas (1.3.4)
- matplotlib (3.5.0)
- scipy (1.7.3)
- skimage (scikit-image 0.18.3)
- seaborn (0.11.2)

The following files are needed to run the analysis:

- main.py
- gui.py
- segmentation.py
- analysis.py 
- plotting.py

Runing the code:

1. Download and extract all the files in contactanalysis_2D. 
2. In the command prompt, with the virtual environment active, navigate to the folder where you extracted the files and run `python main.py`. The GUI will appear to select data input folders (see below). 

## 2. Using the GUI and output

### User parameter settings and input structure

- User settings can be adjusted in the first part of main.py using any editor. 

- In the GUI you can add several input folders, which should contain all the files for one experimental condition. Every input folder should have the tif stacks to analyse in sub-folders whose name includes the 'id_stack' keyword. You can additionally provide a file for flatfield correction of the SLB channel (needs to have the same channels as the main tif stack). This should be present in a separate folder within the input folder and named the same as the folder for the tif stack, but in which 'id-stack' is exactly replaced by 'id_illumination'. Only the tif file with the shortest name in each subfolder will be used.
	e.g. select "SLB_Condition_1" in the GUI when running the program. All folders except marked will be included and their respective tif files analysed. (`id_stack = 'timelapse'`, `id_illumination = 'illumination'`)
	```
		-- SLB_Condition_1		               	# select this folder
		   |-- 2050_01_01_SLB_timelapse_1		 # This folder will be included 
		   |   |-- stack.tif
		   |-- 2050_01_01_SLB_timelapse_2		 # this folder will be included
		   |   |-- stack.tif
		   |   |-- stack_2.tif  		       	# not included
		   |-- 2050_01_01_SLB_illumination_1
		   |   |-- file.tif			           	# stack used to flatfield correct stack.tif in 2020_00_00_SLB_timelapse_1
		   `-- 2050_01_01_SLB_snapshot_1       	# not included, wrong 'id_stack'
			   |-- image_time.tif   
	```
- The code was designed to work with three-color TIRF imaging data. Thus, the tif stack has to contain three channels: One of the SLB in which reduced signal reports close contact zones, one of the cell membrane in TIRF, one of the intracellular Ca signal of the T-cell. The order of these has to be set in main.py.
- In the GUI select which part(s) of the code to run. 'Contact Segmentaton' processes the image stacks to generate segmented image stacks with cell numbers assigned. 'Contact Analysis' will then extract quantitative data from the segmented images. 'Contact Plotting' will adjust the analysis output slighty, calculate a few more metrics and prepare Pandas DataFrames for easier data visualisation.  
Quantitative features are reported on a cell- and contact-based level. Note that only cells that pass the internal quality control ("QC", excluding cells that cut off at the edge) are analysed. 'Contact Analysis' cannot be run unless 'Contact Segmentation' has been selected as well or has been performed beforehand. 

### Output

The code will generate a series of tif stack outputs that allow inspection of the segmentation and filtering:
  - *_SLB_flatfield-corrected.tif: Flatfield corrected SLB channel
  - *_CCZ_filt.tif: Filtered SLB channel 
  - *_CCZ_exc.tif: Segmented SLB channel/detected close contacts. Labelled by cell number.
  - *_CZ.tif: Segmented cell membrane channel/detected cell footprint. Labelled by cell number.
  - *_CCZ_contact_times.tif: Image stack where the pixel values correspond to the contact time of the respective close contact. 
  - *_cell_markers.tif: If not supplied and `redo_seeds = True` in main.py this file will be generated. Once generated, it can be edited to adjust the seed points for assigning cell labels. It is recommended to first run the analysis with only 'Contact Segmentaton' active and then verify the segmentation. If necessary, manually adjust *_cell_markers.tif in an image processing software (Fiji, Napari, ...).

Main outputs in each selected folder:
  - contactgrowth_parameters: saves the input parameters when running the code. 
  - cell_sumary: Gives an overview of the cells dected. 
  - cell_results: Reports for each frame some key characteristics of the contacts formed. Statistical measures are calculated from contact_results. 
  - contact_results: Reports characteristics for each contact at each timepoint. Most detailed output. 
  - Signaling_traces.pdf: an overview of the calcium and cell footprint area for each cell over time

## 3. Accesory scripts

### Triggering probability

triggering_probability.py uses the output of contactanalysis_2D to calculate the triggering probability based on close contact size or duration. 

### Single frame enrichment in close contacts

enrichment.ipynb is a stand-alone Jupyter notebook to calculate local protein enrichment at close contacts. 

Needed is an input folder with multi-color tif stacks where one channel imaged the glycocalyx on the SLB and one channel the protein of interest for which the enrichment should be calculated. 

## 4. Citation

Please consider citing the corresponding paper when using contactanalysis_2D in your work. 

## 5. Troubleshooting

Example TIRF data have been deposited (DOI: 10.5281/zenodo.7509910) alongside the publication which include the expected output and can be used to test the analysis. 


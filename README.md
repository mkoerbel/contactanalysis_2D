# contactanalysis_2D

Analyse 3-colour TIRF image stacks of T cells interacting with supperted lipid bilayers. 



Written by Markus Koerbel. Email mk988@cam.ac.uk

## 1. Setup

System Requirements:

The code has been tested on Windows 10 with an Intel i7-7700K CPU and 32 GB RAM. Depending on the size of the input files at least 16 GB RAM are recommended. 

Installation (on Windwos, about 10 min):

Install Python 3 (www.python.org, tested with Python 3.8). Create a virtual envirnment by running `python -m venv /path/to/new/virtual/environment` (the `python` commnad might have to be adjusted based on your Python installation, often to `python3`) and activate it by `/path/to/new/virtual/environment/Scripts/activate.bat` in the terminal. Once activated, install the required packages via pip (21.3.1) to you virtual environment via `pip install numpy pandas matplotlib scipy scikit-image seaborn`:

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

With the virtual environment active, `python main.py` to run. 

## 2. Using the GUI and output

### User parameter settings and input structure

User Setting can be made in the first part of main.py using any editor. 

- In the GUI you can add several input folders, which normally contain all the files for one experimental condition. Every input folder should have the tif stacks to analyse in sub-folders whose name includes the 'id_stack' word. You can additionally provide a file for flatfield correction of the SLB channel (needs to have the same channels as the main tif stack) which should be present in a folder (for each condition-folder) containing 'id_stack', in which 'id-stack' is exactly replaced by 'id_illumination'. Only the tif file with the shortest name in each subfolder will be used.
	E.g. Select "SLB_Condition_1" in the GUI when running the program. All folders except marked will be included and their respective tif files analysed. (`id_stack = 'timelapse'`, `id_illumination = 'illumination'`)
	```
		-- SLB_Condition_1		               # select this folder
		   |-- 2050_01_01_SLB_timelapse_1
		   |   |-- stack.tif
		   |-- 2050_01_01_SLB_timelapse_2
		   |   |-- stack.tif
		   |   |-- stack_2.tif  		       # not included
		   |-- 2050_01_01_SLB_illumination_1
		   |   |-- file.tif			           # stack used to flatfield correct stack.tif in 2020_00_00_SLB_timelapse_1
		   `-- 2050_01_01_SLB_snapshot_1       # not included
			   |-- image_time.tif   
	```
- The code was designed to work with three color TIRF imaging data (refer to the reference above). Thus, the tiff stack has to contain three channels: One of the SLB in which reduced signal reports close contact zones, one of the cell membrane in TIRF, one of the intracellular Ca signal of the T-cell. The orrder of these has to be set in main.py.
- In the GUI select which part(s) of the code to run. 'Contact Segmentaton' processes the image stacks to generate segmented image stacks with cell numbers assigned. 'Contact Analysis' will then extract quantitative data from the segmented images. 'Contact Plotting' will adjust the analysis output slighty, calculate a few more metrics and prepare Pandas DataFrames for easier data visualisation. 
Quantitative features are reported on a cell and contact based level. Note that only cells that pass the internal quality control ("QC", excluding cells that cut off at the edge) are analysed. 'Contact Analysis' cannot be run unless 'Contact Segmentation' has been selected as well or has been performed beforehand. 

### Output

The code will generate a series of tif stack outputs that allow inspection of the segmentation and filtering:
  - *_SLB_flatfield-corrected.tif: Flatfield corrected SLB channel
  - *_CCZ_filt.tif: Filtered SLB channel 
  - *_CCZ_exc.tif: Segmented SLB channel/detected close contacts. Labelled by cell number.
  - *_CZ.tif: Segmented cell membrane channel/detected cell footprint. Labelled by cell number.
  - *_CCZ_contact_times.tif: Image stack where the pixel values correspond to the contact time of the respective close contact. 
  - *_cell_markers.tif: If not supplied and `redo_seeds = True` in main.py this file will be generated. Once generated, it can be edited to adjust the seed points for assigning cell labels.

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

## 4. Citation

Contactanalysis_2D has been used and is also described in the following publication:

//CITE

Please cite this reference when using contactanalysis_2D in your work. 


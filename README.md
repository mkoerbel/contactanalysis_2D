# contactanalysis_2D

Analyse 3-colour TIRF image stacks of T cells interacting with supperted lipid bilayers. Written by Markus Koerbel. Email mk988@cam.ac.uk

The presented analysis scripts have been used and are also described in the following publication:

//CITE

Please cite this reference when using all or part of this analysis. 

## 1. Installation

Python 3 installation needed. In addition the packages (tested with the versions in brackets)

- numpy (1.21.4) 
- pandas (1.3.4)
- matplotlib (3.5.0)
- scipy (1.7.3)
- skimage (scikit-image 0.18.3)
- seaborn (0.11.2)

Included files:

- main.py
- gui.py
- segmentation.py
- analysis.py 
- plotting.py

## 2. Using the GUI and output

Run main.py. User Setting can be made in the first part of main.py using any editor. 

- In the GUI you can add several input folders, which normally contain all the files for one experimental condition. Every input folder should have the tif stacks to analyse in sub-folders whose name includes the 'id_stack' word. For the SLB channel you can additionally provide a file for flatfield correction (needs to have the same channels as the main tif stack). The code will look for a folder (for each condition-folder) containing 'id_stack', in which 'id-stack' is exactly replaced by 'id_illumination'.
	E.g. All folder except marked will be included. Only the tif file with the shortest name in each subfolder will be used. 'id_stack' = 'timelapse'. 
		- SLB_Condition_1		# select this folder
			|-- 2020_00_00_SLB_timelapse_1
				|-- stack.tif
			|-- 2020_00_00_SLB_timelapse_2
				|-- stack.tif
				|-- stack_2.tif  		# not included
			|-- 2020_00_000_SLB_illumination_1
				|-- file.tif			# stack used to flatfield correct stack.tif in 2020_00_00_SLB_timelapse_1
			|-- 2020_00_00_SLB_snapshot_1  # not included
				|-- image_time.tif    
- Select which part(s) of the code to run. 'Contact Segmentaton' processes the iamge stacks to generate segmented image stacks with cell numbers assigned. 'Contact Analysis' will then extract quantitative data from the segmented images. 'Contact Plotting' will adjust the analysis output slighty, calculate a few more metric and prepare DataFrames for data visualisation. Quantitative features are reported on a cell and contact based level. Note that only cells that pass the internal quality control (excluding cells that cut off at the edge) are analysed. 'Contact Analysis' cannot be run unless 'Contact Segmentation' has been selected as well or has been performed beforehand. 
- Main outputs in each selected folder:
  - contactgrowth_parameters: saves the input parameters when running the code. 
  - cell_sumary: Gives an overview of the cells dected. 
  - cell_results: Reports for each frame some key characteristics of the contacts formed. Statistical measures are calculated from contact_results. 
  - contact_results: Reports characteristics for each contact at each timepoint. Most detailed output. 


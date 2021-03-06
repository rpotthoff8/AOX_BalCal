******* AOX EXAMPLE FILES README *******
******* 29 April 2020, JRP

This folder contains example data files to begin using AOX_BalCal and AOX_Approx. 

Two versions of artificial data are provided in individual folders.
The 'LargeArtificial' folder contains artificial data that was constructed mapping loads to voltages. 
	The full AIAA recommended 6x96 iterative approach polynomial term set will perfectly map the relationship with the exception of noise. 
	As AOX_BalCal uses a non-iterative approach, its polynomial terms will not exactly match this data.
	Normally-distributed artificial noise with a mean of 0 and a standard deviation of 0.7 microV/V was added to the voltages (gage outputs).
The 'LargeArtificial-NoisyPerfect' folder contains artificial data that was constructed mapping voltages to loads.
	AOX_BalCal's full polynomial model will perfectly map the relationship with the exception of noise.
	Normally-distributed artificial noise with a mean of 0 and a standard deviation of 0.01% load capacity was added to the loads.

Each folder contains .csv files for the full dataset and the dataset split into calibration and validation subsets. 
	Calibration and validation are labeled '_cal_rand' and '_val_rand', respectively.
	The subsets were generated by randomly splitting the full dataset for 70% calibration, 30% validation.

Each folder also contains .cal and .val files for each of the .csv files. These may be loaded in AOX_BalCal to quickly load the data without providing the spreadsheet ranges.

Each folder also contains a .ini file. Use the 'Load Settings' button in AOX_BalCal to load this file and the settings for the performing calibration and validation on the example dataset.

Both datasets include tare loads. The true tares are provided in this folder in the file 'Artificial Data Tares (Cheat Sheet).xlsx'.


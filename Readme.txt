Icy automated vector value extraction for foci and nuclei
by Anthony Lagnado and Abbas Ishaq

General automated vector-based foci colocalisation script (Python v3.6.6, Pandas v0.23.3), v1.0, 10/10/2018 - initially written for:
Telomere-associated gamma H2A.X foci (TAF) automated vector colocalisation
by Abbas Ishaq

--------

Table of contents:
1 - Background
2 - Icy image analysis method and parameters
3 - Spreadsheet format:
4 - Python method and parameters:
5 - Citation / Methods:
6 - Other notes:

--------

1 - Background

Colocalisation of telomere foci with H2A.X DNA damage foci is the basis for the TAF senescence biomarker. TAF is used as the primary example in this document.

This script colocalises nuclear telomere and H2A.X foci based on vectors determined by ICY. For TAF colocalisation v1.0, quantification by image has not been integrated.
This should not be a problem for sample-based quantitative statistics, since SD and SEM would be between average for samples. As such, it is a requirement that each sample
to be analysed be placed on separate spreadsheets. 

v1.0 is a framework for quantification by image which is in progress; v1.0 only performs quantification by sample. Completion of v2.0 will promise iteration through 
multiple spreadsheets to automatically analyse multiple samples and images, and collate them into a single spreadsheet.

--------

2 - Icy method and parameters:

Download and install Icy (www.icy.bioimageanalysis.org). Run Icy. If you are running Icy on a more powerful system, tune up the RAM settings under preferences. 
	a - Load your image using Icy. 
	b - Under "Tools", select "Protocols". Load the protocols template file from the "TAF colocalisation" folder. The individual processes are called blocks,
	    and the function of each block is in the respective header.
	c - Ensure that each "Extract channel" block extracts the correct channel. "HK-Means" and "Active Contours" are used to determine the outline
	    of your nuclei. There are two spot detector channels which are interchangeable in the template. Specify one for your telomeres, and another for H2A.X.
	d - Ensure that the correct input files and folder paths are selected.
	e - Modify the save files for each "Workbook to file" block. Ensure that each stream of method leads to the appropriately named save-path block.
	f - For "if file exists" setting, select "Merge sheets". Ensure that first row is not excluded. This is critical for accurate image analysis.

Each method of image acquisition will result in slight variations to the settings in blocks "HK-Means", "Active Contours", and "Spot detector". Thus they may require manual
tuning if your samples have changed, or if your image acquisition methods and exposures are different. Manual activation for "Extract channel" is under "Image/Sequence" tab.
Manual activation for "HK-Means", "Active Contours", and "Spot detector" are under "Detection and Tracking". The default "HK-Means" settings in the protocol is only sensitive
enough for maximum projection of nuclei. Non-maximum projection of nuclei may result in errors with the Python script.

	f - Run the analysis protocol. Go for a coffee.

--------

3 - Spreadsheet format:

- The generated spreadsheet will contain N/A in certain columns. Replace all with "1" (select the column, ctrl+F, replace, "N/A" -> "1").

- Icy can be foregone for a different image analysis software, however ensure that the correct spreadsheet formats are observed (see sample data).

--------

4 - Python method and parameters:

Download and install Anaconda. Launch Spyder through the Anaconda GUI or command line. Ensure the appropriate libraries are installed and up-to-date ("import etc").
The numba library requires an NVIDIA graphics card and cudatoolkit. If cudatoolkit has not been installed, the quickest way would be to launch the Anaconda Prompt and enter
"pip install cudatoolkit". For v1.0, the numba library only improves analysis time by 1%, so numba is not currently required. As datasets increase, it is expected that
efficiency provided by numba will increase. Numba and the parts using it can be commented out (lines 4, 262, 266, 270, 274, 278).

	a - Ensure that the correct paths are selected for input and output files (lines 12-16).
	b - Ensure that the correct values for pixel size etc are entered, removing the constant need to input these values(lines 19-25). 
	    Alternatively, comment out the flat values using # and allow input through the console. Ensure that either method of input is prefaced by the float() function.
	c - Ensure that the correct columns in each file is selected for x, y, and z vector values (lines 32-49).
	d - Run the script. Don't go for a coffee, there is no time.

--------

5 - Citation / Methods:
	a - Cite the following paper - 
	b - For methods, refer to my github repository for this script package . 

--------

6 - Other notes: 

List composition for editing:

all_H2AX = list(zip(x_H2AX, y_H2AX, z_H2AX, width_H2AX, height_H2AX, 
                    depth_H2AX, xmicron_H2AX, ymicron_H2AX, zmicron_H2AX, 
                    sxmicron_H2AX, symicron_H2AX, szmicron_H2AX, 
                    xmicron_H2AX_end, ymicron_H2AX_end, zmicron_H2AX_end))

0  =  x_H2AX
1 =  y_H2AX
2 = z_H2AX
3 = width_H2AX
4 = height_H2AX
5 = depth_H2AX
6 = xmicron_H2AX
7 = ymicron_H2AX
8 = zmicron_H2AX
9 = sxmicron_H2AX
10 = symicron_H2AX
11 = szmicron_H2AX
12 = xmicron_H2AX_end
13 = ymicron_H2AX_end
14 = zmicron_H2AX_end

--------
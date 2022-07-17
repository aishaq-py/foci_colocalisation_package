# Foci_colocalisation_package
Introduction

This script was designed to colocalise 3D dots on multi-channel Z-stack microscopy images, particularly γH2A.X foci and telomere foci. The Python script has been coded to filter all γH2A.X foci and telomere foci in their respective spreadsheets by the nuclear regions in the DAPI spreadsheet, then determine if they are overlapping in the x,y,z axes. 

I - Data generation:

This workflow begins with and requires capture of z-stack images with high resolution at approximately 0.1 µm for pixel length, width, and z-step. Using the 40x oil or 63x oil objectives on our Leica DMi8 microscope are both acceptable. The Leica DM5500 can also be used for these image captures. Considering that the Leica DMi8 costs twice as much to use but captures z-stacks in far less time than the DM5500, save some time and just use the DMi8. At 5 z-stacks for 8 samples per hour for the DMi8, one would be able to work through 64 samples in an 8-hour day. 
The main drawback for generating this amount of data is the sheer requirement of data storage. Each file (sample) would save at between 4 to 8 GB. At this point in time (01/2020), we do not have a sustainable data storage solution. Some suggestions are to acquire 4 TB refurbished hard drives from Western Digital for £50 each and store data sets sorted by treatment or cohort or study, then file them away. Alternatively, we could be allocated more server space by the university, but they seem to be rather hesitant to do this. Alternatively, we could rent a cloud service and code up an SQL database to sort and store the data. Rather than downloading files constantly, we could just pull the data by calling them, and run meta-analyses easily.
The image files will be stored as LIF. This is a series of TIFFs that have been concatenated together into a Leica format container. To access these images, you require a program that uses the Bioformats Java plugin, such as ImageJ or Icy.

-------------------------

II – Icy Image analysis
http://icy.bioimageanalysis.org/
Icy is an open source image analysis software designed and written by a French group. Icy is written in Java and contains ImageJ within it. Thus, while it may look like an intimidating software, Icy is actually quite relatable for most analysts. Icy contains many powerful functionalities far beyond what we will use in the TAF analysis workflow. As someone had previously described though, documentation for Icy functions is absolutely abysmal. 
To begin with, try to acquire Icy version 1.9.7. It seems that the later versions have Bioformats modules that will not open LIF files. There are copies of this version around. Ed, Tengfei, and Abbas have one each.
There are two facets to Icy. The first is that every module, plugin, and function in Icy can be used by themselves (Fig 2A). The second is that all of these modules can be chained into a workflow with settings that are saved, then run over and over for different images. This second function is called a “Protocol” in Icy (Fig 2B). A protocol is composed of various blocks that perform these functionalities. Henceforth, when blocks are referred to, the modules referred to are in the protocol window. If block is not specified, then the document refers to non-protocol variants of the modules. Folders can also be looped through a protocol so that an entire folder of images can be analysed after initialising the protocol without further input from the user. Analysis time for 5 z-stacks (1 sample) takes approximately 5 minutes.
Fig 2. A) Location of all the plugins and functions. B) The protocols button, under tools tab.
To begin the analysis process, load up Icy. Icy can be loaded in a whole window format or separate windows. Toggle this under the ImageJ tab. Next, the “Inspector” window. This window informs you of everything that is being manually done in the workspace, as well as what physical resources are currently being used by Icy (CPU or RAM consumption) shown at the bottom as a running graph for memory and CPU. Image analysis is most reliant on RAM usage. Depending on the size of the image opened, up to 10 GB of RAM can be tied up by Icy or ImageJ at once. Including background programmes and other functions, this means that a computer with 16 GB of RAM will be bottlenecked for resources and will not be usable for anything else during the analysis. If you think that there is a higher than expected RAM consumption by Icy, click on the black box with “Memory” and “CPU” in it. This will attempt to clear up any resources that have been tied up. To prevent any issues with bottlenecking of RAM, use Image analysis computers 1 or 2 in the Image analysis room in CAV. If the maximum RAM value appears to be too low (around 10 to 20 GB), increase this value in the “Preferences”.
Next, open a protocols window. From this window, open one of the TAF analysis protocols. This should serve as a template. The protocol can be thought of to do three separate jobs for each image file. It extracts the blue (Fig 3 block 5), green (Fig 3 block 11) and red (Fig 3 block 2) channels, then performs various tasks with them, and output your data in an Excel spreadsheet. The values for each individual block within the protocol should be manually measured and tailored to your images.
To begin with, determine the orders of your blue, green and red channels, starting with 0 and ending with 2. Put the correct numbers into the extract channel blocks. 

DAPI channel parameters
Then, begin by tuning the values for DAPI. In figure 3, 2 lines come off block 5, joining to intensity projection (block 6) and a later block called ROI statistics. These blocks do not need modification. Intensity projection then leads into the median filter (block 7) which also does not require modification. Median filter then leads into the HK-means for nucleus (block 8) block. This block is for detection of nuclei. To tune this block, open one z-stack image in Icy by dragging and dropping a file onto the main Icy toolbar. Navigate to “Detection and Tracking”, then locate the HK-means button. Navigate to “Image/Sequence” and select “Extract”. Extract your blue channel. Then navigate to “Processing” and select “Projection”. Ensure that your extracted blue channel is selected. Press start on the projection window. You now have a maximum z-projection for blue. Navigate to “Region of Interest” and select “Polygon”. Draw a ROI around a large nucleus to get an idea of what nuclear size range you are expecting on your image. Put this value + about 20% into the HK-means window that was opened earlier. Press start. If done correctly, all nuclei within the size specified will be detected. Put this value into the HK-means block on the protocols window. Then navigate to “Detection and Tracking” and select “Active Contours”. Apply the parameters in the “Active Contours” block to the “Active Contours” window. Press start and determine if these parameters are correctly defining the nuclei. Tweak the parameters as necessary and update the “Active Contours” block. Close the HK-means window and extracted blue channel.

γH2A.X and telomere channel parameters
The steps for γH2A.X and telomere parameters optimisation are the same, just performed on different extracted channels. Navigate to “Image/Sequence” and select “Extract” to extract the green channel. Navigate to “Detection and Tracking” and select “Spot Detector”. Select “Detection” from the left menu. Input the detection parameters from the “Spot Detector” block into the “Spot Detector” window and run. Tweak the spot detector settings as necessary. Apply the settings to the “Spot Detector” block.

Saving and loading
There are three blocks for saving, labelled accordingly for DAPI, H2AX and telomeres. Set your save paths on the blocks, then ensure that “merge sheets” is selected, not one of the other selections. Navigate to the top left of the protocols page. Select the folder in which you have your LIF files. Once these parameters have been set, save your protocol and start run. 

-------------------------

III – Data Clean-up
Data clean-up will be the same for all three Excel spreadsheets. Open one spreadsheet. Copy column A and paste into column B. Highlight column B, press CTRL + F, then select replace. Replace all “Dataset” in column B with “Image”. Highlight column A, press CTRL + F, then select replace. Replace all “ - Mark_and_Find 00*/Position00* - channel: *” with “ “. Highlight column F, press CTRL+F, then select replace. Replace all “ALL” with “1”. Highlight column I, press CTRL+F, then select replace. Replace all “N/A” with “1”. Save and close. Repeat for all three spreadsheets. There is also a VBA macro that does all this for each spreadsheet.

-------------------------

IV - Script execution
*At the time of initiating this learning project, I had not learned to execute python scripts from the command line. Thus parameters had to be changed from within the script. This was normally okay because the parameters of the regions of interest did not change very much between experiments.

To prepare the interface for working with the script, download and install Anaconda Navigator. Then, from either the Anaconda Navigator UI or the Anaconda Prompt command line, launch the Spyder IDE. Load in the Python script. Modify lines 15 to 18 for your input directory. Then, modify lines 25, 27, and 31 for the values you have. 25 for pixel length or width, 27 for z step size, and 31 for the TAF count above which you wish to consider your nuclei senescent. Normally, 1 is a good value. 

-------------------------

V – Checking your data
Two files will be generated for each analysis. The diagnostics file contains all raw data for the summary analysis in the main output file. To determine if a TAF is actually colocalising, open the diagnostics file. Find the sample, then the nucleus you are looking for. Find the coordinates of the nucleus. Open the corresponding z-stack image, then locate the nucleus of interest on the image.
Another way to confirm that the data output is accurate is to analyse the distribution of TAF. There should be a decreasing distribution from 1 to 10+, with the highest category being the no-TAF category. This function is not yet coded in, but there is an Excel template for doing this manually. 

-------------------------

Publication
This script was initially targetted to be submitted to Scientific Reports. However, no funding bodies are tied to this work, and the idea seems moot at this point considering the progress that has been made in image analysis research since then. Also there are some bugs in the code that have not been worked out.

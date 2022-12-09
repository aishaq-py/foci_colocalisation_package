from skimage import io
import numpy as np

# Read the immunofluorescence microscopy images
red_channel = io.imread('red_channel.tif')
green_channel = io.imread('green_channel.tif')
blue_channel = io.imread('blue_channel.tif')

# Compute the colocalization of red and green foci in 3D space
colocalization = np.logical_and(red_channel, green_channel)

# Compute the number of colocalizing foci per nuclei
nuclei_labels, num_nuclei = ndimage.label(blue_channel)
coloc_per_nuclei = ndimage.sum(colocalization, nuclei_labels, range(1, num_nuclei + 1))

# Get the coordinates of the colocalizing foci
coordinates = np.where(colocalization)

# Print the results
print("Number of colocalizing foci per nuclei:", coloc_per_nuclei)
print("Coordinates of colocalizing foci:", coordinates)

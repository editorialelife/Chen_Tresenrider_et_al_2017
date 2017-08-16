FISH-FINDER: Version 4
David McSwiggen, Dec 2016

This script uses gaussian detection/localization as applied by
the SLIMfast implementation of MTT in order to detect spots coming
from a maximum projection image of FISH data in multiple colors.
Then it will attempt to pair up spots from the different images
in order to determine colocalization.
 
To effectively use this code, one should have images in a .TIFF
format, and have at least two separate FISH color channels as well
as DAPI (or similar) to allow for the segmentation of cells.
Filepaths to each image should be provided, and the code can take
into account different replicates uding the "Position_number" and
"Replicate_number" variables.

After detection of all spots, a threshold must be used to determine
which are truly mRNA and which are background noise. This threshold
will depend on the parameters used to aquire the images, and may be
very different from those used in the Chen, Tresenrider, et. al.
manuscript.

For more information, or to address specific questions, please
email dmcswiggen@berkeley.edu or direct them to the corresponding
authors of the Chen, Tresenrider, et. al. manuscript.

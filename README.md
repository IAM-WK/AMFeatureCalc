# AMFeatureCalc

This is a tool to determine features quasi-continuously over a large microscope image of an Additively Manufactured Foam Structure (AMF). Quasi-continuous refers to the fact that not only one value for each feature is determined on the entire microscope image, instead many sub-images are extracted from the microscope image, on which all features are determined in each case. This results in a quasi-continuous curve for each feature. The distance between the sub-images can be changed by the user. If required, all features can also be calculated on a single image. The features are the following:
- Total porosity
- Texture features from the publication by Haralick [1]
  - ASM (Output not available in this version - must be added on your own)
  - DifferenceEntropy (Output not available in this version - must be added on your own)
  - TextureHomogeneity (Output not available in this version - must be added on your own)
  - Entropy
- Features extracted by MorphoLibJ [2]
  - Pore size in mm^2
  - Perimeter
  - Circularity
  - Ellipse radius 1
  - Ellipse radius 2
  - Ellipse orientation
  - Ellipse elongation
- Pore amount in sub-image (1/mm^2)

All features from MorphoLibJ are calculated by the list of pores within the sub-image. You can choose between median and average. The process of the algorithm is as follows:
1) Binaryise image with local threshold.
2) Label pores with MorphoLibJ. Then filter the list (by pore size).
3) Filter between different shaped pores as desired.
4) Start with first sub-image and calculate features.
5) Rasterise the sub-image over the entire image until the limit is reached.
6) Save the data in a CSV file.

A more detailed description of the algorithm can be found in the publication (see Cite section).

[1] R. M. Haralick, K. Shanmugam, and I. Dinstein. “Textural features for image classification”. IEEE Transactions on systems, man, and cybernetics SMC-3.6 (1973), pp. 610–621.\
[2] D. Legland, I. Arganda-Carreras, and P. Andrey. “MorphoLibJ: integrated library and plugins for mathematical morphology with ImageJ”. Bioinformatics 32.22 (July 2016), pp. 3532–3534.

## Citation

... To be done

## Prerequisites

- Fiji software (https://imagej.net/software/fiji/downloads) is required.
- Within Fiji the plugin MorphoLibJ must be active.
- Create microscope image(s) and save in usual format as greyscale image (.png, .tif, ...).
- Create a folder with all microscope images if batch processing is desired.


## Usage

The script can be executed in Fiji. Beforehand, the parameters in the header of the script must be adjusted. Short description of use:
- Set all parameters for calculating the features according to the preferences. Suggested values for parameters can be found in the publication (see section "Cite").
- Specify image resolution of the microscope image.
- Specify layer thickness of the AMF.
- Select if a parameter and which parameter should be varied (a separate CSV file is created for each parameter value).
- Select calculation mode (quasi-continuous or on the entire image).
- Specify file path to folder with images, file names of images, storage path and name of output CSV file.

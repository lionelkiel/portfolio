# Extracts from [Report](Report.pdf)
## Introduction and Method
To find out about the structure of some sample of interest one can use electron diffraction to image it. Electron diffraction is done in a Low-Energy Electron Microscope (LEEM for short). It is done by accelerating electrons towards the sample where they will scatter and create a diffraction image. The LEEM scheme used is shown here: ![Scheme](In-place correction of diffraction distortion in Low-Energy Electron Microscopy/Figures/Scheme.png)

Due to the complexity of this LEEM there is often a strong distortion in the diffraction images. This is problematic as it prevents the measuring of angles and distances. A method that corrects this distortion was successfully implemented in similar machines, however they require replacing of the sample. Replacing a sample causes the LEEM setup to be decalibrated. Therefore, an in-place correction of the diffraction distortion needs to be done.

To do the in-place correction of the diffraction distortion the following will be done. Firstly, the transformation function for the diffraction distortion will be closely approximated. Then, the inverse transformation will be applied to the diffraction images. Finding the transformation function consists of a couple of steps explained in more detail in the report. The main idea is to generate a calibration image of which its corrected and diffraction distorted image is known. Using this the full transformation function can then be derived using interpolation.

## Conclusion
The method presented in this thesis allows to accurately measure angles and distances, even when the distortion in the original data is so bad it is not possible at all. This is very clearly seen by the fact that measuring an angle in the uncorrected data has a variance of 64.1 degrees (squared) while it only has a variance of 5.5 degrees (squared) when corrected. The results suggest that if less of the data points are used, with a wide spread, the variance in both the distance and angles only increase slightly. With 25% of the data the variance in both the distance and angles doubles. However, with an even wider spread of data this could be achieved with even less data points and less of an increase in the variance. More research into this can thus be done. Other suggestions for future research include a better implementation of the nearest peak detection and cluster detection algorithms.

# Run notebook ([BRP](BRP.ipynb))
- The jupyter notebook should be executed from top to bottom.
- The first cell shows all the dependencies for the notebook
- The second cell contains the functions used for the distortion correction
- The third cell contains the functions used for the various metrics used to check the quality of the correction
- Cells 4-6 applies the correction on some data (this data is not given in the repository)
- The rest of the cells analyze how well the correction does when only using certain subsets of the given data

[![Build Status](https://travis-ci.org/CMET-UGent/MicroRaman.svg?branch=master)](https://travis-ci.org/CMET-UGent/MicroRaman)

# MicroRaman
*******************
- **Authors**: [Frederiek-Maarten Kerchkof](mailto:FrederiekMaarten.Kerckhof@UGent.be), [Benjamin Buysschaert](mailto:Benjamin.Buysschaert@Ugent.be), [Dmitry Khalenkow](mailto:Dmitry.Khalenkow@Ugent.be), [Jasmine Heyse](mailto:Jasmine.Heyse@ugent.be), [Ruben Props](mailto:Ruben.Props@ugent.be) and [Cristina Garcia Timermans](mailto:Cristina.GarciaTimermans@ugent.be).

The goal of this package is to provide a standardized and automated workflow for spectral analysis in a fast and objective way, removing the post-processing bias of different operators. 

If you use this package, please consider citing the original publication:  

García-Timermans, C., Rubbens, P., Kerckhof, F. M., Buysschaert, B., Khalenkow, D., Waegeman, W., Skirtach, A. G. & Boon, N. (2018). Label-free Raman characterization of bacteria calls for standardized procedures. Journal of microbiological methods, 151, 69-75.

Install the package:
```R
library("devtools")
install_github("CMET-UGent/MicroRaman")
```

## Available functions

Functions  | Actions
------------| -----------
baselinecorrplot | plot the result of optimal baseline correction
hs2mq | converts a `hyperSpec::hyperSpec` object directly to a `MALDIquant::MassSpectrum` object
intervalplot | plotting function that shows if (trimmed) wavelengths are evenly spaced
iterationsplot | plot the effect of iterations on baseline correction
grid_arrange_shared_legend | Function for creating a common legend for 2 ggplot2 figures.
iterationsplot | plots baseline correction for several iterations on the spectral data to determine the optimal number of iterations
mq2hs | Convert a MassSpectrum object to a hyperSpec object
ramplot | Convenience function to plot Raman spectra (from `MALDIquant::MassSpectrum`)
SCA | calculates the spectral contracst angle between two vectors
SCA.diss | calculates the spectral contrast angle (SCA) between all cells in a hyperSpec object
wlcutter | Set intensities of a given bandpass to a defined value (anomalous peak removal)


## Avaiable datasets

Some datasets are included in the package. They allow the examples and vigenttes
to be run. They can be loaded using:
```R
library("MicroRaman")
data("<name of dataset>")
```

Dataset name | Data contents
-------------| ----------------
mass.spectra.baseline.corr |
mdqs | Raman spectrometry test single-cell data of Escherichia coli LMG2092 biological replicates (list of MALDIquant objects)
spx.all | Hyperspec dataset (list of hyperSpec class objects) from single-cell data of *E. coli* LMG 2092 in triplicate in NB and LB test.


[![Build Status](https://travis-ci.org/CMET-UGent/MicroRaman.svg?branch=master)](https://travis-ci.org/CMET-UGent/MicroRaman)

# MicroRaman
*******************
- **Authors**: [Frederiek-Maarten Kerchkof](mailto:FrederiekMaarten.Kerckhof@UGent.be), [Cristina Garcia Timermans](mailto:Cristina.GarciaTimermans@ugent.be), and [Ruben Props](mailto:Ruben.Props@ugent.be).

- **Contributors**:
[Benjamin Buysschaert](mailto:Benjamin.Buysschaert@Ugent.be) and [Jasmine Heyse](mailto:Jasmine.Heyse@ugent.be) 

The goal of this package is to provide a standardized and automated workflow for Raman spectra analysis. 

If you use this package, please consider citing the original publication in which is was first used:  

García-Timermans, C., Rubbens, P., Kerckhof, F. M., Buysschaert, B., Khalenkow, D., Waegeman, W., Skirtach, A. G. & Boon, N. (2018). Label-free Raman characterization of bacteria calls for standardized procedures. Journal of microbiological methods, 151, 69-75.

Install the package:
```R
library("devtools")
install_github("CMET-UGent/MicroRaman")
```

## Core functions

Functions  | Actions
------------| -----------
hs_import | Import Thermo Galactic's spc file format data into the R environment
hs_preprocess | Preprocesses the data using the Garcia-Timermans et al. (2020) workflow
hs_contrast | Calculate contrasts between spectra of specified groups of cells
hs_hclust | Hierarchical clustering of Raman spectra (with or without bootstrap support)
hs_hclust_cutoff | Visualization of distance cut-off in hclust plots
hs_PCA | Principal Component Analysis of Raman spectra
hs_tsne | t-distributed stochastic neighbor embedding of Raman spectra 
hs_phenoRam | Calculation of Hill diversity numbers for each individual Raman spectrum
hs_SCAdiss | Calculates the spectral contrast angle (SCA) between all cells in a hyperSpec object
hs_SCA_conv_itol | Convert SCA dissimilarity matrix to itol-compatible object

## Convenience functions
Functions  | Actions
------------| -----------
hs_conv_mq | converts a `hyperSpec::hyperSpec` object directly to a `MALDIquant::MassSpectrum` object
mq_conv_hs | converts a `MALDIquant::MassSpectrum` object directly to a `hyperSpec::hyperSpec` object
mq_plot |
mq_baseline_plot |
mq_iter_plot |
intervalplot |
model_fit_stats |
pred_r_squared |
PRESS |
SCA | Calculates the spectral contrast angle between two vectors
wlcutter |

## Avaiable datasets

Some datasets are included in the package. They allow the examples and vigenttes
to be run. They can be loaded using:
```R
library("MicroRaman")
data("<name of dataset>")
```

Dataset name | Data contents
-------------| ----------------
mass.spectra.baseline.corr |baseline-corrected raman spectra
mdqs | Raman spectrometry test single-cell data of Escherichia coli LMG2092 biological replicates (list of MALDIquant objects)
spx.all | Hyperspec dataset (list of hyperSpec class objects) from single-cell data of *E. coli* LMG 2092 in triplicate in NB and LB test.
SCA.similarity | Spectral contrast angle similarities calculated on the spx.all/mdqs data

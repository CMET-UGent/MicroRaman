[![Build Status](https://travis-ci.org/CMET-UGent/MicroRaman.svg?branch=master)](https://travis-ci.org/CMET-UGent/MicroRaman)

# MicroRaman
*******************
- **Authors**: [Frederiek-Maarten Kerchkof](mailto:FrederiekMaarten.Kerckhof@UGent.be), [Cristina Garcia Timermans](mailto:Cristina.GarciaTimermans@ugent.be), and [Ruben Props](mailto:Ruben.Props@ugent.be).

- **Contributors**:
[Benjamin Buysschaert](mailto:Benjamin.Buysschaert@Ugent.be) and [Jasmine Heyse](mailto:Jasmine.Heyse@ugent.be) 

The goal of this package is to provide a standardized and automated workflow for Raman spectra analysis. 

If you use this package, please consider citing the original publication in which is was first used:  

García-Timermans, C., Rubbens, P., Kerckhof, F. M., Buysschaert, B., Khalenkow, D., Waegeman, W., Skirtach, A. G. & Boon, N. (2018). Label-free Raman characterization of bacteria calls for standardized procedures. Journal of microbiological methods, 151, 69-75.

García‐Timermans, C., Rubbens, P., Heyse, J., Kerckhof, F.‐M., Props, R., Skirtach, A.G., Waegeman, W. and Boon, N. (2020), Discriminating Bacterial Phenotypes at the Population and Single‐Cell Level: A Comparison of Flow Cytometry and Raman Spectroscopy Fingerprinting. Cytometry. doi:10.1002/cyto.a.23952

Install the package:
```R
library("devtools")
install_github("CMET-UGent/MicroRaman")
```

## Core functions

Functions  | Actions | Functional
------------ | ----------- | -----------
hs_import | Import Thermo Galactic's spc file format data into the R environment | YES
hs_preprocess | Preprocesses the data using the Garcia-Timermans et al. (2020) workflow | YES
hs_contrast | Calculate contrasts between spectra of specified groups of cells | NO
hs_hclust | Hierarchical clustering of Raman spectra (with or without bootstrap support) | YES
hs_hclust_cutoff | Visualization of distance cut-off in hclust plots | YES
hs_type | Clusters spectra using partitioning around medoids  | NO
hs_PCA | Principal Component Analysis of Raman spectra | YES
hs_tsne | t-distributed stochastic neighbor embedding of Raman spectra  | YES
hs_phenoRam | Calculation of Hill diversity numbers for each individual Raman spectrum | YES
hs_RF | Train Random Forest classifier to distinguish between groups of cells | NO
hs_RF_pred | Predict using Random Forest classifier on new data | NO
hs_SCAdiss | Calculates the spectral contrast angle (SCA) between all cells in a hyperSpec object | YES
hs_SCA_conv_itol | Convert SCA dissimilarity matrix to itol-compatible object | NO

## Convenience functions
Functions  | Actions | Functional
------------| ----------- | -----------
hs_conv_mq | Converts a `hyperSpec::hyperSpec` object directly to a `MALDIquant::MassSpectrum` object | YES
mq_conv_hs | Converts a `MALDIquant::MassSpectrum` object directly to a `hyperSpec::hyperSpec` object | YES
hs_tidy_filenames | Tidies up hyperspec spectral IDs | YES
mq_plot | | NO
mq_baseline_plot | | NO
mq_iter_plot | | NO
intervalplot | | NO
model_fit_stats | | NO
pred_r_squared | | NO
PRESS | | NO
SCA | Calculates the spectral contrast angle between two vectors | YES
wlcutter | | NO

## Available datasets

Some datasets are included in the package. They allow the examples and vigenttes
to be run. They can be loaded using:
```R
library("MicroRaman")
data("<name of dataset>")
```

Dataset name | Data contents
-------------| ----------------
mdqs | Raman spectrometry test single-cell data of Escherichia coli LMG2092 biological replicates (list of MALDIquant objects)
spx.all | Hyperspec dataset (list of hyperSpec class objects) from single-cell data of *E. coli* LMG 2092 in triplicate in NB and LB test.

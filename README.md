# MicroRaman
*******************
- **Authors**: Frederiek-Maarten Kerchkof [FrederiekMaarten.Kerckhof@UGent.be], Benjamin Buysschaert [Benjamin.Buysschaert@Ugent.be], Dmitry Khalenkow [Dmitry.Khalenkow@Ugent.be] and Cristina Garcia Timermans [Cristina.GarciaTimermans@ugent.be]

The goal of this package is to provide a standardized automated workflow for spectral analysis in a fast and objective way, removing the post-processing bias of different operators. 



Install the package:
```R
library("devtools")
install_github("CMET-UGent/MicroRaman")
```

## Available functions

Functions  | Actions
------------| -----------
hs2mq | converts a `hyperSpec::hyperSpec` object directly to a `MALDIquant::MassSpectrum` object
hyperSpec2MALDI | converts the extracted information from  `hyperSpec::hyperSpec` object  to a `MALDIquant::MassSpectrum` object
intervalplot | plotting function that shows if (trimmed) wavelengths are evenly spaced


## Avaiable datasets

Some datasets are included in the package. They allow the examples and vigenttes
to be run. They can be loaded using:
```R
library("MicroRaman")
data("<name of dataset>")
```

Dataset name | Data contents
-------------| ----------------
cell.names   | simple list of filenames of spectra of E. coli LMG 2092 in triplicate in NB and LB
hs           | a hyperSpec object (empty) with our Raman wavelengths
hyperspecobj | a hyperSpec object (empty) with our Raman wavelengths
LMG2092_Reproducibility | large matrix with raw data extracted from the . coli LMG 2092 in triplicate in NB and LB test.
mdqs | a list of `MALDIquant::MassSpectrum` objects of single-cell data of Escherichia coli LMG2092 biological replicates (see Vignettes)
raw.data | large matrix with raw data extracted from the . coli LMG 2092 in triplicate in NB and LB test.

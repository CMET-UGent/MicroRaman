---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
---

# Load library & data


```r
library("MicroRaman")
# hs.x <- hs_import(path = "./test_data/")
# summary(hs.x)

# Load data from package
data("hs_example")
hs.x <- hs_example
summary(hs.x)
```

```
## hyperSpec object
##    64 spectra
##    4 data columns
##    1024 data points / spectrum
## wavelength: Delta * tilde(nu)/cm^-1 [numeric] -410.3763 -405.6340 ... 3096.027
## data:  (64 rows x 4 columns)
##    1. z: T/(degree * C) [numeric] 0 0 ... 0
##    2. z.end: T/(degree * C) [numeric] 8.437222e-19 8.437222e-19 ... 8.437222e-19
##    3. spc: I/"a. u." [matrix, array1024] 18326 18140 ... 16376
##    4. filename: filename [character] ./test_data/gfp+_001_spec.data 1.spc ./test_data/gfp+_002_spec.data 1.spc ... ./test_data/gfp+_065_spec.data 1.spc
```

# Tidy data


```r
# Spectra IDs contain full/absolute path
# This is annoying for downstream analysis
# These can be shortened as such:
head(hs.x@data$filename, 5)
```

```
## [1] "./test_data/gfp+_001_spec.data 1.spc" "./test_data/gfp+_002_spec.data 1.spc"
## [3] "./test_data/gfp+_004_spec.data 1.spc" "./test_data/gfp+_005_spec.data 1.spc"
## [5] "./test_data/gfp+_006_spec.data 1.spc"
```

```r
hs.x <- hs_tidy_filenames(hs.x, remove_pattern = "path")
head(hs.x@data$filename, 5)
```

```
## [1] "gfp+_001_spec.data 1.spc" "gfp+_002_spec.data 1.spc" "gfp+_004_spec.data 1.spc"
## [4] "gfp+_005_spec.data 1.spc" "gfp+_006_spec.data 1.spc"
```

# Preprocess data


```r
# Denoise/clean
hs.x.proc <- hs_preprocess(hs.x, smooth = FALSE)
```

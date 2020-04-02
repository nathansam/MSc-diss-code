# MSC dissertation and supporting material

The dissertation document can be found in a pdf format and is called
"nsc_msc_diss"

Included in this repo is the version of `CircadianTools`, an R package developed
for the dissertation, used for analysis in the dissertation. This package can be
installed using the instructions included in the CircadianTools appendix found
in the main document or by following the instructions in the below section.

For ease of use, an uncompressed folder containing all of the source code for
CircadianTools can be found in the R directory.

The R markdown directory contains R markdown files which can be run to reproduce
the results presented in the dissertation. Also included in this directory is
any saved file which is produced by the R markdown files.

The Raw_Data folder contains the raw *T. saltator* transcriptomics data which
was used for this project which has been acquired from
[*Identification and temporal expression of putative circadian clock transcripts in the amphipod crustacean Talitrus saltator*](https://doi.org/10.7717/peerj.2555)
by O'Grady *et al.*

## CircadianTools Installation

Please note either a MacOS or Linux/GNU system is recommended for CircadianTools
due to the use of Forked clusters for parallel calculations.

Installing the CircadianTools R package is best done through the use of the 
```devtools``` package as this will ensure  all of the packages which
CircadianTools depends on will also be installed.
Running the below code in R will install the devtools package. 

```{R}
install.packages("devtools")
```

Before CircadianTools can be installed, dependencies from Bioconductor must also
be installed. This can be done by using the following R command:

```{R}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rain", "seqinr")
```

The most up-to-date version of CircadianTools can be installed from Github by
entering the below line in to the R console:

```{R}
devtools::install_github("nathansam/CircadianTools")
```

Included in this repo is the most recent version of the CircadianTools at the
time of submission. This version of the package can be installed by the
`devtools::install_local("path")` function where "path" is the path to the
included "CircadianTools_1.0.0.tar.gz" file  For example, if the
"CircadianTools_1.0.0.tar.gz" file is moved to my linux desktop then
CircadianTools can be installed by the following command:

```{R}
devtools::install_local("/home/nathan/desktop/CircadianTools_1.0.0.tar.gz")
```
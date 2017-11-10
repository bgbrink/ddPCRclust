# ddPCRclust
## A package for automated quantification of multiplexed ddPCR data.

This R package was designed to automatically quantify the events of a multiplexed ddPCR reaction.
During a ddPCR run, each marker gene is fluorescently labeled with a combination of FAM and/or HEX fluorophore, giving it a unique footprint in the two-dimensional space represented by the intensities for each color channel. The position of each droplet within this space reveals, how many and, more importantly, which marker genes it contains. Thus, droplets that belong to the same marker cluster together. However, correctly identifying and labelling these clusters is not trivial, since one droplet can contain more than one marker. To accurately quantify the droplet count for each marker, it is crucial to both identify the clusters and label them correctly, based on their position. 
<p align="center">
<img 
src="https://cloud.githubusercontent.com/assets/11661112/25387153/1b9928d2-29ca-11e7-97b3-ce67694eed0f.png"  
alt="Example B1"
width="400">
<img 
src="https://cloud.githubusercontent.com/assets/11661112/25387163/20f98c90-29ca-11e7-8a9e-a7efaaf82fd9.png" 
alt="Example G1"
width="400">
</p>

For robustness, ddPCRclust incorporates adapted versions of three established, independent clustering algorithms: the [flowDensity](https://bioconductor.org/packages/release/bioc/html/flowDensity.html) algorithm, published in 2012 by M. Jafar Taghiyar and Mehrnoush Malek, the clustering algorithm [SamSPECTRAL](https://bioconductor.org/packages/release/bioc/html/SamSPECTRAL.html), a version of spectral clustering adapted to flow cytometry data and the [flowPeaks](https://bioconductor.org/packages/release/bioc/html/flowPeaks.html) package, developed by Yongchao Ge and Stuart C. Sealfon. The results are combined into a cluster ensemble. This enhances the precision and the agreement between the three approaches provides a measure of confidence for clustering results.
<p align="center">
<img 
src="https://cloud.githubusercontent.com/assets/11661112/25387160/1e5eea02-29ca-11e7-871b-e2e3cd2639ec.png"  
alt="Result B1"
width="400">
<img 
src="https://cloud.githubusercontent.com/assets/11661112/25387164/224c4876-29ca-11e7-9f8f-557b7e515f0a.png" 
alt="Example G1"
width="400">
</p>

## Installation
You can install this package like any other package from GitHub using devtools
```R
library(devtools)
install_github("bgbrink/dropclust")
```
**Disclaimer:** This method currently only works, when the GNU Scientific Library (GSL) is installed on your machine, because one of the dependencies (flowPeaks) needs it in order to compile. Otherwise you will have to install the packages from bioconductor manually before installing dropclust:
```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite(c('flowDensity', 'SamSPECTRAL', 'flowPeaks'))
```

## Usage
This package was written in close cooperation with the BC Cancer Agency in Vancouver, Canada. Please read their recently published [manuscript](https://doi.org/10.1371/journal.pone.0161274) for details on the background and how to produce the necessary data.

The raw data should be *csv* files. Each file represents a two-dimensional data frame. Each row within the data frame represents a single droplet, each column the respective intensities per colour channel:

Ch1 Amplitude | Ch2 Amplitude 
--- | --- 
2360.098 |	6119.26953
2396.3916 |	1415.31665
2445.838 |	6740.79639
2451.63867 |	1381.74683
2492.55884 |	1478.19617
2519.6355 |	7082.25049
&#8942; | &#8942;

Since one experiment most likely consists of many different files, naming them apropriately is important in order to keep things organized. We chose to use a unique identifier in each filename of the form `"^[[:upper:]][[:digit:]][[:digit:]]$"` (A01, A02, A03, B01, B02, ...), which is usually included automatically by the ddPCR machine. A set of eight example files is included in this package. 

In order to use all functions of this package, it is also necessary to create a template with more information about this experiment. The template has to be a *csv* file with a header, which contains information about each of the raw data files according to their unique identifier, as explained above. A template for the eight example files is also included in this package.

*> Name of your experiment, channel1=HEX, channel2=FAM, annotations(date, experimentor, etc)*

Well|Sample type|No of markers|Marker 1|Marker 2|Marker 3|Marker 4
---|---|---|---|---|---|---
B01|Blood|4|a|b|c|d
G01|FFPE|4|a|b|c|d
F02|Blood|3|a||c|d
D03|FFPE|3|a||c|d
A04|FFPE|4|a|b|c|d
G07|Cell line|3|a||c|d
G08|Cell line|3|a||c|d
E09|FFPE|2|||c|d

### Example Files
We provide eight examplary ddPCR files under ```ddPCRclust/inst/extdata/```

Run the algorithm using the provided examples with the following command:
```R
# Run ddPCRclust
exampleFiles <- list.files(paste0(find.package("ddPCRclust"), "/extdata"), full.names = TRUE)
result <- ddPCRclust(files = exampleFiles[1:8], template = exampleFiles[9])
```
All functions are documented. You can find additional information using the help function of R: `?ddPCRclust`

## License
    ddPCRclust
    Copyright (C) 2017  Benedikt G. Brink, Bielefeld University
                        BC Cancer Agency Branch (“BCCA”)
                        
    CAREFULLY READ THE FULL TERMS AND CONDITIONS INCLUDED IN THE FILE 'LICENSE'.  
    
    LICENSE TO USE: BCCA hereby grants to You a personal, non-exclusive, non-transferable, limited license 
    to use the Product solely for internal, non-commercial use for non-profit research or educational purposes
    only on the terms and conditions contained in this Agreement. The Product may be installed on a maximum of 
    one machine per Qualified User at Your premises only.  A copy of the Product installed on a single common 
    machine may be shared for internal use by Qualified Users only. In order to be a “Qualified User”, an 
    individual must be a student, researcher, professor, instructor or staff member of a non-profit educational
    institution or organization who uses the Product solely for non-profit research or educational purposes.


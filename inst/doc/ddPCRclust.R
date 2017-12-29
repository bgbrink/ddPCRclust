### R code from vignette source 'ddPCRclust.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: ddPCRclust.Rnw:36-38 (eval = FALSE)
###################################################
## library(devtools)
## install_github("bgbrink/ddPCRclust")


###################################################
### code chunk number 3: ddPCRclust.Rnw:43-46 (eval = FALSE)
###################################################
## ## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R")
## biocLite("ddPCRclust")


###################################################
### code chunk number 4: ddPCRclust.Rnw:62-73
###################################################
# Run ddPCRclust
library(ddPCRclust)
exampleFiles <- list.files(paste0(find.package("ddPCRclust"), "/extdata"), full.names = TRUE)
result <- ddPCRclust(files = exampleFiles[3], template = exampleFiles[9])

# Plot the results
library(ggplot2)
p <- ggplot(data = result$results$B01$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = TRUE) +
    ggtitle("ddPCRclust example")+theme_bw() + theme(legend.position="none")
p


###################################################
### code chunk number 5: ddPCRclust.Rnw:94-106
###################################################
# Run the flowDensity based approach
library(ddPCRclust)
exampleFiles <- list.files(paste0(find.package("ddPCRclust"), "/extdata"), full.names = TRUE)
file <- read.csv(exampleFiles[3])
densResult <- runDensity(file = file, numOfMarkers = 4)

# Plot the results
library(ggplot2)
p <- ggplot(data = densResult$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = TRUE) +
    ggtitle("flowDensity example")+theme_bw() + theme(legend.position="none")
p


###################################################
### code chunk number 6: ddPCRclust.Rnw:114-126
###################################################
# Run the SamSPECTRAL based approach
library(ddPCRclust)
exampleFiles <- list.files(paste0(find.package("ddPCRclust"), "/extdata"), full.names = TRUE)
file <- read.csv(exampleFiles[3])
samResult <- runSam(file = file, numOfMarkers = 4)

# Plot the results
library(ggplot2)
p <- ggplot(data = samResult$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = TRUE) +
     ggtitle("SamSPECTRAL example")+theme_bw() + theme(legend.position="none")
p


###################################################
### code chunk number 7: ddPCRclust.Rnw:136-148
###################################################
# Run the flowPeaks based approach
library(ddPCRclust)
exampleFiles <- list.files(paste0(find.package("ddPCRclust"), "/extdata"), full.names = TRUE)
file <- read.csv(exampleFiles[3])
peaksResult <- runPeaks(file = file, numOfMarkers = 4)

# Plot the results
library(ggplot2)
p <- ggplot(data = peaksResult$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = TRUE) +
     ggtitle("flowPeaks example")+theme_bw() + theme(legend.position="none")
p



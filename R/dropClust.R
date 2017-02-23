#' dropClust: A package for automated quantification of multiplexed ddPCR data.
#'
#' The dropClust algorithm is used to automatically quantify the events of a multiplexed ddPCR reaction.
#' In order to determine the correct droplet count for each marker, it is crucial to both identify all clusters and label them correctly based on their position.
#' To be more robust, three different algorithms for flow cytometry data were adapted and extended in order to analyze ddPCR data.
#'
#' @section Foo functions:
#' The foo functions ...
#'
#' @docType package
#' @name dropClust
NULL
#> NULL

#' Find the clusters 
#'
#' Description.
#'
#' @param file The input files. More specifically, csv files with two dimensions, each dimension representing the intensity for one color.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up.
#' @param sensitivity An integer between 0.1 and 2 determining sensitivity of the initial clustering, e.g. the number of clusters. A higher value means more clusters are being found. Standard is 1.
#' @return
#' \item{data}{The original input data minus the removed events (for plotting)}
#' \item{confidence} The agreement between the different clustering results in percent. If all algorithms calculated the same result, the clustering is likely to be correct, thus the confidence is high.
#' \item{counts}{The droplet count for each cluster.}
#' @export
#' @import parallel
#' @examples
#' # Run dropClust:
#' file <- "example.csv"
#' result <- runDropClust(file)
#' 
#' # Plot the result:
#' p <- ggplot(data = result[[1]]$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
#' p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = T) + ggtitle("example") + theme_bw()+ theme(legend.position="none")
#' p
#' 
#' # Run dropClust with multiple files:
#' files <- c("example.csv", "example2.csv", "example3.csv", "example4.csv")
#' result <- runDropClust(files)
runDropClust <- function(files, numOfMarkers = 4, sensitivity = 1) {
  library(parallel)
  
  if(length(files>1)) {
    csvFiles <- mclapply(files, read.csv)
    dens_result <- mcmapply(runDensity, file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity, SIMPLIFY = F)
    sam_result <- mcmapply(runSam, file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity, SIMPLIFY = F)
    peaks_result <- mcmapply(runPeaks, file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity, SIMPLIFY = F)
    superResults <- mcmapply(createEnsemble, dens_result, sam_result, peaks_result, plex, csvFiles, SIMPLIFY = F)
  } else {
    csvFiles <- read.csv(files)
    dens_result <- runDensity(file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity)
    sam_result <- runSam(file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity)
    peaks_result <- runPeaks(file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity)
    superResults <- createEnsemble(dens_result, sam_result, peaks_result, plex, csvFiles)
  }
}




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
#' \item{confidence}{The agreement between the different clustering results in percent. If all algorithms calculated the same result, the clustering is likely to be correct, thus the confidence is high.}
#' \item{counts}{The droplet count for each cluster.}
#' @export
#' @import parallel
#' @examples
#' # Run dropClust
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = T)
#' result <- runDropClust(exampleFiles)
#' # Plot the results
#' exportPlots(result, "./results/", "examples")
runDropClust <- function(files, numOfMarkers = 4, sensitivity = 1) {
  library(parallel)
  
  if(length(files>1)) {
    csvFiles <- mclapply(files, read.csv)
    dens_result <- mcmapply(dens_wrapper, file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity, SIMPLIFY = F)
    sam_result <- mcmapply(sam_wrapper, file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity, SIMPLIFY = F)
    peaks_result <- mcmapply(peaks_wrapper, file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity, SIMPLIFY = F)
    superResults <- mcmapply(ensemble_wrapper, dens_result, sam_result, peaks_result, numOfMarkers, csvFiles, SIMPLIFY = F)
  } else {
    csvFiles <- read.csv(files)
    dens_result <- dens_wrapper(file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity)
    sam_result <- sam_wrapper(file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity)
    peaks_result <- peaks_wrapper(file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity)
    superResults <- ensemble_wrapper(dens_result, sam_result, peaks_result, numOfMarkers, csvFiles)
  }
}


#' Plot the algorithm's results with ggplot2
#'
#' A convinience function that takes the results of the droplClust algorithm and plots them using the ggplot2 library and a custom colour palette.
#'
#' @param data The result of the dropClust algorithm
#' @param directory The parent directory where the plots should saved
#' @param sensitivity The name of this experiment. A subfolder in the parent directory will be created with its name to keep things organized.
#' @export
#' @import ggplot2
#' @examples
#' # Run dropClust
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = T)
#' result <- runDropClust(exampleFiles)
#' # Plot the results
#' exportPlots(result, "./results/", "examples")
exportPlots <- function(data, directory, experiment) {
  library(ggplot2)
  ifelse(!dir.exists(file.path(directory, experiment)), dir.create(file.path(directory, experiment)), FALSE)
  
  for (i in 1:length(data)) {
    result <- data[[i]]
    if (is.null(result$data)) {
      next
    }
    p <- ggplot(data = result$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
    if (length(unique(result$data$Cluster)) > 8) {
      cbPalette <- c("#999999", "#f272e6","#e5bdbe","#bf0072","#cd93c5", "#1fba00","#5e7f65","#bdef00","#2c5d26","#ffe789","#4a8c00", "#575aef","#a3b0fa","#005caa","#019df8", "#bc8775")
    } else if (length(unique(result$data$Cluster)) > 4) {
      cbPalette <- c("#999999", "#d800c4","#fca3a7","#bb004e", "#70cf56","#8b9d61","#ccd451", "#bc8775")
    } else if (length(unique(result$data$Cluster)) > 2) {
      cbPalette <- c("#999999", "#8d5286","#b42842", "#c2ff79","#076633", "#bc8775")
    } else {
      cbPalette <- c("#999999", "#bc8775")
    }
    p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = T) + ggtitle(paste(experiment, i)) + theme_bw()+ theme(legend.position="none") + scale_colour_manual(values=cbPalette)
    ggsave(paste0(directory,"/",experiment,"/",experiment,"_", i, ".png"), p)
  }
}


dens_wrapper <- function(file, sensitivity=1, numOfMarkers) {
  result <- tryCatch(runDensity(file, sensitivity, numOfMarkers), error = function(e) {
    print(e)
  })
}

sam_wrapper <- function(file, sensitivity=1, numOfMarkers) {
  result <- tryCatch(runSam(file, sensitivity, numOfMarkers), error = function(e) {
    print(e)
  })
}

peaks_wrapper <- function(file, sensitivity=1, numOfMarkers) {
  result <- tryCatch(runPeaks(file, sensitivity, numOfMarkers), error = function(e) {
    print(e)
  })
}

ensemble_wrapper <- function(dens_result, sam_result, peaks_result, plex, csvFiles) {
  result <- tryCatch(createEnsemble(dens_result, sam_result, peaks_result, plex, csvFiles), error = function(e) {
    print(e)
  })
}


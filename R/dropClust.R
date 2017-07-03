# ##############################   dropClust   ################################# #  
# Copyright (C) 2017  Benedikt G. Brink, Bielefeld University                    #
#                                                                                #
# This program is free software; you can redistribute it and/or                  #
# modify it under the terms of the GNU General Public License                    #
# as published by the Free Software Foundation; either version 2                 #
# of the License, or (at your option) any later version.                         #
#                                                                                #
# This program is distributed in the hope that it will be useful,                #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                 #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                   #
# GNU General Public License for more details.                                   #
#                                                                                #
# You should have received a copy of the GNU General Public License              #
# along with this program; if not, write to the Free Software                    #
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. #
# ############################################################################## #

#' dropClust 
#' A package for automated quantification of multiplexed ddPCR data
#'
#' The dropClust algorithm can automatically quantify the events of ddPCR reaction with up to four markers.
#' In order to determine the correct droplet count for each marker, it is crucial to both identify all clusters and label them correctly based on their position.
#' For more information on what data can be analyzed and how a template needs to be formatted, please check the project repository on github.
#'
#' @section Usage:
#' The main function of the package is \code{\link{runDropClust}}. This function runs the algorithm with one or multiple files, automatically distributing them amongst all cpu cores using the \link[parallel]{parallel} package (parallelization does not work on windows).
#' Afterwards, the results can be exported in different ways, using \code{\link{exportPlots}}, \code{\link{exportToExcel}} and \code{\link{exportToCSV}}.
#' Once the clustering is finished, copies per droplet (CPD) for each marker can be calculated using \code{\link{calculateCPDs}}.
#' 
#' These functions provide access to all functionalities of the dropClust package. However, expert users can directly call some internal functions of the algorithm, if they find it necessary. 
#' Here is a list of all available supplemental functions: \cr
#' \code{\link{runDensity}} \cr
#' \code{\link{runSam}} \cr
#' \code{\link{runPeaks}} \cr
#' \code{\link{createEnsemble}} 
#'
#' @docType package
#' @name dropClust
"_PACKAGE"
#> [1] "_PACKAGE"
library(parallel)
library(ggplot2)
library(openxlsx)

#' Run the dropClust algorithm
#'
#' This is the main function of this package. It automatically runs the dropClust algorithm on one or multiple csv files containing the raw data from a ddPCR run with up to 4 markers. 
#'
#' @param files The input file(s). More specifically, csv files with two dimensions, each dimension representing one color choannel of the ddPCR reaction.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up. Can be ignored if a template is provided. 
#' Else, a vector with length equal to \code{length(files)} should be provided, containing the number of markers used for the respective reaction.
#' @param sensitivity An integer between 0.1 and 2 determining sensitivity of the initial clustering, e.g. the number of clusters. A higher value means more clusters are being found. Standard is 1.
#' @param template A csv file containing information about the individual ddPCR runs. An example template is provided with this package. For more information, please check the repository on github. 
#' @param fast Run a simple version of the algorithm that is about 10x faster. For clean data, this can already deliver very good results. In any case useful to get a quick overview over the data.
#' @return
#' \item{results}{The results of the dropClust algorithm. It contains three fields: \cr
#' \code{data} The original input data minus the removed events (for plotting) 
#' \code{confidence} The agreement between the different clustering results in percent
#' If all parts of the algorithm calculated the same result, the clustering is likely to be correct, thus the confidence is high\cr
#' \code{counts} The droplet count for each cluster
#' }
#' \item{annotations}{The metatdata provided in the header of the template. It contains four fields: \cr
#' \code{Name} The name given to this ddPCR experiment \cr
#' \code{Ch1} Color channel 1 (usually HEX) \cr
#' \code{Ch2} Color channel 2 (usually FAM) \cr
#' \code{descriptions} Additional descriptions about this ddPCR experiment (e.g. date, exprimentor, etc.)
#' }
#' \item{template}{A parsed dataframe containing the template, if one was provided.}
#' @export
#' @import parallel
#' @examples
#' # Run dropClust
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = TRUE)
#' result <- runDropClust(files = exampleFiles[1:8], template = exampleFiles[9])
#'
runDropClust <- function(files, numOfMarkers = 4, sensitivity = 1, template = NULL, fast = FALSE) {
  time <- proc.time()
  ids <- annotations <- vector()
  markerNames <- list(c("M1", "M2", "M3", "M4"))
  ids <- unlist(lapply(files, function(x) {grep("^[[:upper:]][[:digit:]][[:digit:]]$", unlist(strsplit(x, "_")), value = T)}))
  if (!is.null(template)) {
    header <- readLines(template, n = 1)
    header <- gsub("[\\\"]", "", header)
    if (substr(header, start = 1, stop = 1) == ">") {
      annotations <- unlist(strsplit(substring(header, 2), ","))
      annotations <- trimws(annotations)
      if (length(annotations) < 4) {
        warning("Missing fields in header. Please use the following layout: experiment_name, ch1, ch2, descriptions")
        annotations <- NULL
      } else {
        ch1 <- strsplit(annotations[2], "1=")[[1]][2]
        ch2 <- strsplit(annotations[3], "2=")[[1]][2]
        if (any(is.na(c(ch1, ch2)))) {
          warning("Unknown fields in header. Please use the following layout: experiment_name, ch1, ch2, descriptions")
          annotations <- NULL
        } else {
          annotations <- c(annotations[1], ch1, ch2, paste(annotations[4:length(annotations)], collapse = ","))
          names(annotations) <- c("Name", "Ch1", "Ch2", "Descriptions")
        }
      }
      template <- read.csv(template, skip = 1)
      numOfMarkers <- lapply(ids, function(x) {unlist(template[which(template[,1] == x), 3])})
      markerNames <- lapply(ids, function(x) {unlist(template[which(template[,1] == x), 4:7])})
    } else {
      stop(paste("Invalid Template file! This file starts with:\n", substr(header, start = 1, stop = 10)))
    }
  }
  
  if (Sys.info()['sysname'] == "Windows") {
    nrOfCores <- 1
  } else {
    nrOfCores <- detectCores()
  }
  dens_result <- sam_result <- peaks_result <- rep(0, length(files))
  names(dens_result) <- names(sam_result) <- names(peaks_result) <- ids
  if(length(files)>1) {
    csvFiles <- mclapply(files, read.csv, mc.cores = nrOfCores)
    names(csvFiles) <- ids
    dens_result <- mcmapply(dens_wrapper, file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity, markerNames=markerNames, SIMPLIFY = F, mc.cores = nrOfCores)
    if (!fast) {
      sam_result <- mcmapply(sam_wrapper, file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity, markerNames=markerNames, SIMPLIFY = F, mc.cores = nrOfCores)
      peaks_result <- mcmapply(peaks_wrapper, file=csvFiles, numOfMarkers=numOfMarkers, sensitivity=sensitivity, markerNames=markerNames, SIMPLIFY = F, mc.cores = nrOfCores)
    }
    superResults <- mcmapply(ensemble_wrapper, dens_result, sam_result, peaks_result, csvFiles, SIMPLIFY = F, mc.cores = nrOfCores)
  } else {
    csvFiles <- read.csv(files)
    dens_result <- dens_wrapper(file=csvFiles, numOfMarkers=numOfMarkers[[1]], sensitivity=sensitivity, markerNames = markerNames[[1]])
    if (!fast) {
      sam_result <- sam_wrapper(file=csvFiles, numOfMarkers=numOfMarkers[[1]], sensitivity=sensitivity, markerNames = markerNames[[1]])
      peaks_result <- peaks_wrapper(file=csvFiles, numOfMarkers=numOfMarkers[[1]], sensitivity=sensitivity, markerNames = markerNames[[1]])
    }
    superResults <- list()
    superResults[[ids]] <- ensemble_wrapper(dens_result, sam_result, peaks_result, csvFiles)
  }    
  time <- (proc.time()-time)[3]
  return(list(results=superResults, annotations=annotations, template=template, runtime=time))
}


#' Plot the algorithms results with ggplot2
#'
#' A convinience function that takes the results of the dropClust algorithm and plots them using the ggplot2 library and a custom colour palette.
#'
#' @param data The result of the dropClust algorithm
#' @param directory The parent directory where the files should saved. A new folder with the experiment name will be created (see below).
#' @param annotations Some basic metadata about the ddPCR reaction. If you provided \code{\link{runDropClust}} a template, this paramater can be filled with the corresponding field in the result. 
#' Otherwise, you have to provide a character vector containing a name and the the color channels, e.g. \code{c(Name="ddPCR_01-04-2017", Ch1="HEX", Ch2="FAM")}
#' @param invert Invert the axis, e.g. x = Ch2.Amplitude, y = Ch1.Amplitude
#' @export
#' @import ggplot2
#' @examples
#' # Run dropClust
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = TRUE)
#' result <- runDropClust(files = exampleFiles[1:8], template = exampleFiles[9])
#' 
#' # Plot the results
#' dir.create("./Results")
#' exportPlots(data = result$results, directory = "./Results/", annotations = result$annotations)
#'
exportPlots <- function(data, directory, annotations, format = ".png", invert = FALSE) {

  directory <- normalizePath(directory, mustWork = T)
  ifelse(!dir.exists(paste0(directory,"/",annotations[1])), dir.create(paste0(directory,"/",annotations[1])), FALSE)
  
  for (i in 1:length(data)) {
    id <- names(data[i])
    result <- data[[i]]
    if (is.null(result$data)) {
      next
    }
    if (invert) {
      p <- ggplot(data = result$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
    } else {
      p <- ggplot(data = result$data, mapping = aes(x = Ch1.Amplitude, y = Ch2.Amplitude))
    }
    if (length(unique(result$data$Cluster)) > 8) {
      cbPalette <- c("#999999", "#f272e6","#e5bdbe","#bf0072","#cd93c5", "#1fba00","#5e7f65","#bdef00","#2c5d26","#ffe789","#4a8c00", "#575aef","#a3b0fa","#005caa","#01c8fe", "#bc8775")
    } else if (length(unique(result$data$Cluster)) > 4) {
      cbPalette <- c("#999999", "#d800c4","#fca3a7","#bb004e", "#70cf56","#8b9d61","#ccd451", "#bc8775")
    } else if (length(unique(result$data$Cluster)) > 2) {
      cbPalette <- c("#999999", "#8d5286","#b42842", "#c2ff79","#076633", "#bc8775")
    } else {
      cbPalette <- c("#999999", "#bc8775")
    }
    if (invert) {
      p <- p + geom_point(aes(color = factor(Cluster)), size = .4, na.rm = T) + ggtitle(id) + theme_bw()+ theme(legend.position="none") + 
        scale_colour_manual(values=cbPalette) + labs(x = paste(annotations[3], "Amplitude"), y = paste(annotations[2], "Amplitude")) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
    } else {
      p <- p + geom_point(aes(color = factor(Cluster)), size = .4, na.rm = T) + ggtitle(id) + theme_bw()+ theme(legend.position="none") + 
        scale_colour_manual(values=cbPalette) + labs(x = paste(annotations[2], "Amplitude"), y = paste(annotations[3], "Amplitude")) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
    }
    ggsave(paste0(directory,"/",annotations[1],"/", id, format), p, dpi = 1200)
  }
}

#' Export the algorithms results to an Excel file
#'
#' A convinience function that takes the results of the droplClust algorithm and exports them to an Excel file.
#'
#' @param data The result of the dropClust algorithm
#' @param directory The parent directory where the files should saved. A new folder with the experiment name will be created (see below).
#' @param annotations Some basic metadata about the ddPCR reaction. If you provided \code{\link{runDropClust}} a template, this paramater can be filled with the corresponding field in the result. 
#' Otherwise, you have to provide a character vector containing a name and the the color channels, e.g. \code{c(Name="ddPCR_01-04-2017", Ch1="HEX", Ch2="FAM")}
#' @param raw Boolean which determines if the annotated raw data should be exported along with the final counts. Basically, a third column will be added to the original data, which contains the cluster number to which this point was assigned to.
#' Useful for example to visualize the clustering later on. (Warning: this can take a while!)
#' @export
#' @import openxlsx
#' @examples
#' # Run dropClust
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = TRUE)
#' result <- runDropClust(files = exampleFiles[1:8], template = exampleFiles[9])
#' 
#' # Export the results
#' dir.create("./Results")
#' exportToExcel(data = result$results, directory = "./Results/", annotations = result$annotations)
#'
exportToExcel <- function(data, directory, annotations, raw = FALSE) {

  directory <- normalizePath(directory, mustWork = T)
  ifelse(!dir.exists(paste0(directory,"/", annotations[1])), dir.create(paste0(directory,"/",annotations[1])), FALSE)
  dataToWrite <- NULL
  for (i in 1:length(data)) {
    id <- names(data[i])
    result <- data[[i]]
    if (is.null(result$data)) {
      next
    }
    conf <- result$confidence
    names(conf) <- "Confidence"
    tmp <- as.data.frame(t(c(result$counts, conf)))
    if (is.null(dataToWrite)) {
      dataToWrite <- rbind(dataToWrite, tmp)
    } else if (ncol(dataToWrite) < ncol(tmp)) {
      dataToWrite <- merge(dataToWrite, tmp, all.x=T)
      dataToWrite <- rbind(dataToWrite, tmp[match(colnames(dataToWrite), colnames(tmp))])
    } else {
      tmp <- merge(dataToWrite, tmp, all.y=T)
      dataToWrite <- rbind(dataToWrite, tmp[match(colnames(dataToWrite), colnames(tmp))])
    }
    if (raw) {
      file <- paste0(directory,"/",annotations[1],"/", id, "_annotated_raw.xlsx")
      write.xlsx(result$data, file = file)
    }
  }
  mynames <- c("Empties","1","2","3","4","1+2","1+3","1+4","2+3","2+4","3+4","1+2+3","1+2+4","1+3+4","2+3+4","1+2+3+4","Removed","Total","Confidence")
  dataToWrite <- t(dataToWrite[,match(mynames, colnames(dataToWrite))])
  colnames(dataToWrite) <- names(data)
  file <- paste0(directory,"/",annotations[1],"/", annotations[1], "_results.xlsx")
  write.xlsx(dataToWrite, file = file, colNames = T, rowNames = T)
}

#' Export the algorithms results to a csv file
#'
#' A convinience function that takes the results of the droplClust algorithm and exports them to a csv file.
#'
#' @param data The result of the dropClust algorithm
#' @param directory The parent directory where the files should saved. A new folder with the experiment name will be created (see below).
#' @param annotations Some basic metadata about the ddPCR reaction. If you provided \code{\link{runDropClust}} a template, this paramater can be filled with the corresponding field in the result. 
#' Otherwise, you have to provide a character vector containing a name and the the color channels, e.g. \code{c(Name="ddPCR_01-04-2017", Ch1="HEX", Ch2="FAM")}
#' @param raw Boolean which determines if the annotated raw data should be exported along with the final counts. Basically, a third column will be added to the original data, which contains the cluster number to which this point was assigned to.
#' Useful for example to visualize the clustering later on. (Warning: this can take a while!)
#' @export
#' @examples
#' # Run dropClust
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = T)
#' result <- runDropClust(exampleFiles)
#' 
#' # Export the results
#' dir.create("./Results")
#' exportToCSV(data = result$results, directory = "./Results/", annotations = result$annotations)
#'
exportToCSV <- function(data, directory, annotations, raw = FALSE) {
  directory <- normalizePath(directory, mustWork = T)
  ifelse(!dir.exists(paste0(directory,"/", annotations[1])), dir.create(paste0(directory,"/",annotations[1])), FALSE)
  dataToWrite <- NULL
  for (i in 1:length(data)) {
    id <- names(data[i])
    result <- data[[i]]
    if (is.null(result$data)) {
      next
    }
    conf <- result$confidence
    names(conf) <- "Confidence"
    tmp <- as.data.frame(t(c(result$counts, conf)))
    if (is.null(dataToWrite)) {
      dataToWrite <- rbind(dataToWrite, tmp)
    } else if (ncol(dataToWrite) < ncol(tmp)) {
      dataToWrite <- merge(dataToWrite, tmp, all.x=T)
      dataToWrite <- rbind(dataToWrite, tmp[match(colnames(dataToWrite), colnames(tmp))])
    } else {
      tmp <- merge(dataToWrite, tmp, all.y=T)
      dataToWrite <- rbind(dataToWrite, tmp[match(colnames(dataToWrite), colnames(tmp))])
    }
    if (raw) {
      file <- paste0(directory,"/",annotations[1],"/", id, "_annotated_raw.csv")
      write.csv(result$data, file = file)
    }
  }
  mynames <- c("Empties","1","2","3","4","1+2","1+3","1+4","2+3","2+4","3+4","1+2+3","1+2+4","1+3+4","2+3+4","1+2+3+4","Removed","Total","Confidence")
  dataToWrite <- t(dataToWrite[,match(mynames, colnames(dataToWrite))])
  colnames(dataToWrite) <- names(data)
  file <- paste0(directory,"/",annotations[1],"/", annotations[1], "_results.csv")
  write.csv(dataToWrite, file = file)
}

# wrapper function for exception handling in mcmapply
dens_wrapper <- function(file, sensitivity=1, numOfMarkers, markerNames) {
  missingClusters <- which(markerNames == "")
  result <- tryCatch(runDensity(file[,c(2,1)], sensitivity, numOfMarkers, missingClusters), error = function(e) {
    print(e)
  })
}

# wrapper function for exception handling in mcmapply
sam_wrapper <- function(file, sensitivity=1, numOfMarkers, markerNames) {
  missingClusters <- which(markerNames == "")
  result <- tryCatch(runSam(file[,c(2,1)], sensitivity, numOfMarkers, missingClusters), error = function(e) {
    print(e)
  })
}

# wrapper function for exception handling in mcmapply
peaks_wrapper <- function(file, sensitivity=1, numOfMarkers, markerNames) {
  missingClusters <- which(markerNames == "")
  result <- tryCatch(runPeaks(file[,c(2,1)], sensitivity, numOfMarkers, missingClusters), error = function(e) {
    print(e)
  })
}

# wrapper function for exception handling in mcmapply
ensemble_wrapper <- function(dens_result, sam_result, peaks_result, file) {
  if (is.numeric(dens_result)) {
    dens_result <- NULL
  }
  if (is.numeric(sam_result)) {
    sam_result <- NULL
  }
  if (is.numeric(peaks_result)) {
    peaks_result <- NULL
  }
  result <- tryCatch(createEnsemble(dens_result, sam_result, peaks_result, file[,c(2,1)]), error = function(e) {
    print(e)
  })
}

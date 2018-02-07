# ###################### ddPCRclust ######################## #
# Copyright (C) 2017 Benedikt G. Brink, Bielefeld University # 
# ########################################################## #

#' ddPCRclust
#' A package for automated quantification of multiplexed ddPCR data
#'
#' The ddPCRclust algorithm can automatically quantify the events of ddPCR reaction with up to four markers.
#' In order to determine the correct droplet count for each marker, 
#' it is crucial to both identify all clusters and label them correctly based on their position.
#' For more information on what data can be analyzed and how a template needs to be formatted, 
#' please check the project repository on github.
#'
#' @section Usage:
#' The main function of the package is \code{\link{ddPCRclust}}. This function runs the algorithm with one or multiple files, 
#' automatically distributing them amongst all cpu cores using the \link[parallel]{parallel} package 
#' (parallelization does not work on windows). Afterwards, the results can be exported in different ways, 
#' using \code{\link{exportPlots}}, \code{\link{exportToExcel}} and \code{\link{exportToCSV}}.
#' Once the clustering is finished, copies per droplet (CPD) for each marker can be calculated using \code{\link{calculateCPDs}}.
#'
#' These functions provide access to all functionalities of the ddPCRclust package. 
#' However, expert users can directly call some internal functions of the algorithm, if they find it necessary.
#' Here is a list of all available supplemental functions: \cr
#' \code{\link{runDensity}} \cr
#' \code{\link{runSam}} \cr
#' \code{\link{runPeaks}} \cr
#' \code{\link{createEnsemble}}
#'
#' @docType package
#' @name ddPCRclust
#' @import clue ggplot2 plotrix
#' @importFrom parallel mclapply mcmapply detectCores
#' @importFrom flowDensity deGate
#' @importFrom SamSPECTRAL SamSPECTRAL
#' @importFrom flowPeaks flowPeaks
#' @include cluster_functions.R
#' @include functions.R
"_PACKAGE"
# > [1] '_PACKAGE'


#' Read the csv files from your disk
#'
#' This function reads the raw csv files for ddPCRclust from disk and returns the experiment data.
#' Please refer to the vignette for more information on how these files need to be formatted.
#'
#' @param files The input file(s), specifically csv files. Each file represents a two-dimensional data frame.
#' Each row within the data frame represents a single droplet, each column the respective intensities per colour channel.
#' @return
#' \item{files}{A data frame composed of the experiment data}
#' \item{ids}{The file ids, e.g. A01, A02, etc.}
#' @export
#' @examples
#' # Read files
#' files <- readFiles(exampleFiles[1:8])
#'
readFiles <- function(files) {
  ids <- unlist(lapply(files, function(x) {
    grep("^[[:upper:]][[:digit:]][[:digit:]]$", unlist(strsplit(x, "_")), value = TRUE)
  }))
  csvFiles <- lapply(files, utils::read.csv)
  names(csvFiles) <- ids
  return(list(ids = ids, data = csvFiles))
}

#' Read a template file from disk
#'
#' This function reads a template file for ddPCRclust from disk and returns a run template and annotations.
#' Please refer to the vignette for information on how this file need to be formatted.
#'
#' @param template A csv file containing information about the individual ddPCR runs. 
#' An example template is provided with this package. For more information, please check the vignette or the repository on github.
#' @return 
#' \item{annotations}{The metatdata provided in the header of the template. It contains four fields: \cr
#' \code{Name} The name given to this ddPCR experiment \cr
#' \code{Ch1} Color channel 1 (usually HEX) \cr
#' \code{Ch2} Color channel 2 (usually FAM) \cr
#' \code{descriptions} Additional descriptions about this ddPCR experiment (e.g. date, exprimentor, etc.)}
#' \item{template}{A parsed dataframe containing the template.}
#' @export
#' @examples
#' # Read template
#' template <- readTemplate(exampleFiles[9])
#'
readTemplate <- function(template) {
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
        annotations <- c(annotations[1], ch1, ch2, paste(annotations[4:length(annotations)], 
                                                         collapse = ","))
        names(annotations) <- c("Name", "Ch1", "Ch2", "Descriptions")
      }
    }
    template <- utils::read.csv(template, skip = 1)
  } else {
    stop(paste("Invalid Template file! This file starts with:\n", substr(header, start = 1, stop = 10)))
  }
  return(list(annotations = annotations, template = template))
}


#' Run the ddPCRclust algorithm
#'
#' This is the main function of this package. 
#' It automatically runs the ddPCRclust algorithm on one or multiple csv files 
#' containing the raw data from a ddPCR run with up to 4 markers.
#'
#' @param files The input data obtained from the csv files. For more information, please see \code{\link{readFiles}}.
#' @param template A data frame containing information about the individual ddPCR runs. 
#' An example template is provided with this package. For more information, please see \code{\link{readTemplate}}.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up. 
#' Can be ignored if a template is provided. Else, a vector with length equal to \code{length(files)} should be provided, 
#' containing the number of markers used for the respective reaction.
#' @param sensitivity A number between 0.1 and 2 determining sensitivity of the initial clustering, 
#' e.g. the number of clusters. A higher value means more clusters are being found. Standard is 1.
#' @param similarityParam If the distance of a droplet between two or more clusters is very similar, it will not be counted for either. 
#' The standard it 0.95, i.e. at least 95\% similarity. A sensible value lies between 0 and 1, 
#' where 0 means none of the 'rain' droplets will be counted and 1 means all droplets will be counted.
#' @param distanceParam When assigning rain between to clusters, typically the bottom 20\% are assigned to the lower cluster 
#' and remaining 80\% to the higher cluster. This parameter changes the ratio.
#' @param fast Run a simpler version of the algorithm that is about 10x faster. For clean data, 
#' this can already deliver very good results. In any case useful to get a quick overview over the data.
#' @param multithread Distribute the algorithm amongst all CPU cores to speed up the computation.
#' @return
#' \item{results}{The results of the ddPCRclust algorithm. It contains three fields: \cr
#' \code{data} The original input data minus the removed events (for plotting) \cr
#' \code{confidence} The agreement between the different clustering results in percent
#' If all parts of the algorithm calculated the same result, the clustering is likely to be correct, thus the confidence is high\cr
#' \code{counts} The droplet count for each cluster
#' }
#' @export
#' @examples
#' # Read files
#' files <- readFiles(exampleFiles[1:8])
#' 
#' # Read template
#' template <- readTemplate(exampleFiles[9])
#' 
#' # Run ddPCRclust
#' result <- ddPCRclust(files, template)
#'
#' # Plot the results
#' library(ggplot2)
#' p <- ggplot(data = result$B01$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
#' p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = TRUE) +
#'   ggtitle('B01 example')+theme_bw() + theme(legend.position='none')
#' p
#'
ddPCRclust <- function(files, template, numOfMarkers = 4, sensitivity = 1, similarityParam = 0.95, 
                       distanceParam = 0.2, fast = FALSE, multithread = FALSE) {
    if (!is.numeric(numOfMarkers) || numOfMarkers > 4 || numOfMarkers < 1) {
        stop("Invalid argument for numOfMarkers. Currently only the detection of 1-4 markers is supported.")
    } else if (!is.numeric(similarityParam) || similarityParam > 1 || similarityParam < 0) {
        stop("Invalid argument for similarityParam. Only values between 0 and 1 are supported.")
    } else if (!is.numeric(distanceParam) || distanceParam > 1 || distanceParam < 0) {
        stop("Invalid argument for distanceParam. Only values between 0 and 1 are supported.")
    } else if (!is.numeric(sensitivity) || sensitivity > 2 || sensitivity < 0.1) {
        stop("Invalid argument for sensitivity. Only values between 0.1 and 2 are supported.")
    }
    
    time <- proc.time()
    if (!is.null(template)) {
      numOfMarkers <- lapply(files$ids, function(x) {
        unlist(template$template[which(template$template[, 1] == x), 3])
      })
      markerNames <- lapply(files$ids, function(x) {
        unlist(template$template[which(template$template[, 1] == x), 4:7])
      })
      numOfMarkers <- as.numeric(numOfMarkers)
    } else {
      numOfMarkers <- 4
      markerNames <- list(c("M1", "M2", "M3", "M4"))
    }
    if (Sys.info()["sysname"] == "Windows" || !multithread) {
        nrOfCores <- 1
    } else {
        nrOfCores <- detectCores()
    }
    dens_result <- sam_result <- peaks_result <- rep(0, length(files$data))
    names(dens_result) <- names(sam_result) <- names(peaks_result) <- files$ids
    if (length(files$data) > 1) {
        dens_result <- mcmapply(dens_wrapper, file = files$data, numOfMarkers = numOfMarkers, 
            sensitivity = sensitivity, markerNames = markerNames, similarityParam = similarityParam, 
            distanceParam = distanceParam, SIMPLIFY = FALSE, mc.cores = nrOfCores)
        if (!fast) {
            sam_result <- mcmapply(sam_wrapper, file = files$data, numOfMarkers = numOfMarkers, 
                sensitivity = sensitivity, markerNames = markerNames, similarityParam = similarityParam, 
                distanceParam = distanceParam, SIMPLIFY = FALSE, mc.cores = nrOfCores)
            peaks_result <- mcmapply(peaks_wrapper, file = files$data, numOfMarkers = numOfMarkers, 
                sensitivity = sensitivity, markerNames = markerNames, similarityParam = similarityParam, 
                distanceParam = distanceParam, SIMPLIFY = FALSE, mc.cores = nrOfCores)
        }
        superResults <- mcmapply(ensemble_wrapper, dens_result, sam_result, peaks_result, 
            files$data, SIMPLIFY = FALSE, mc.cores = nrOfCores)
    } else {
        dens_result <- dens_wrapper(file = files$data, numOfMarkers = numOfMarkers[[1]], 
            sensitivity = sensitivity, markerNames = markerNames[[1]], similarityParam = similarityParam, 
            distanceParam = distanceParam)
        if (!fast) {
            sam_result <- sam_wrapper(file = files$data, numOfMarkers = numOfMarkers[[1]], 
                sensitivity = sensitivity, markerNames = markerNames[[1]], similarityParam = similarityParam, 
                distanceParam = distanceParam)
            peaks_result <- peaks_wrapper(file = files$data, numOfMarkers = numOfMarkers[[1]], 
                sensitivity = sensitivity, markerNames = markerNames[[1]], similarityParam = similarityParam, 
                distanceParam = distanceParam)
        }
        superResults <- list()
        superResults[[files$ids]] <- ensemble_wrapper(dens_result, sam_result, peaks_result, 
            files$data)
    }
    time <- (proc.time() - time)[3]
    return(superResults)
}


#' Plot the algorithms results with ggplot2
#'
#' A convinience function that takes the results of the ddPCRclust algorithm,
#' plots them using the ggplot2 library and a custom colour palette and saves the plots to a folder.
#'
#' @param data The result of the ddPCRclust algorithm
#' @param directory The parent directory where the files should saved.
#'  A new folder with the experiment name will be created (see below).
#' @param annotations Some basic metadata about the ddPCR reaction. 
#' If you provided \code{\link{ddPCRclust}} a template, 
#' this paramater can be filled with the corresponding field in the result.
#' Otherwise, you have to provide a character vector containing a name and the the color channels, 
#' e.g. \code{c(Name='ddPCR_01-04-2017', Ch1='HEX', Ch2='FAM')}
#' @param format Which file format to use. Can be either be a device function (e.g. png), 
#' or one of 'eps', 'ps', 'tex' (pictex), 'pdf', 'jpeg', 'tiff', 'png', 'bmp', 'svg' or 'wmf' (windows only). 
#' See also \code{\link{ggsave}}
#' @param invert Invert the axis, e.g. x = Ch2.Amplitude, y = Ch1.Amplitude
#' @return None
#' @examples
#' # Read files
#' files <- readFiles(exampleFiles[1:8])
#' 
#' # Read template
#' template <- readTemplate(exampleFiles[9])
#' 
#' # Run ddPCRclust
#' result <- ddPCRclust(files, template)
#'
#' # Export the plots
#' dir.create('./Results')
#' exportPlots(data = result, directory = './Results/', annotations = result$annotations)
#'
exportPlots <- function(data, directory, annotations, format = "png", invert = FALSE) {
    directory <- normalizePath(directory, mustWork = TRUE)
    ifelse(!dir.exists(paste0(directory, "/", annotations[1])), dir.create(paste0(directory, 
        "/", annotations[1])), FALSE)
    
    for (i in seq_along(data)) {
        id <- names(data[i])
        result <- data[[i]]
        if (is.null(result$data)) {
            next
        }
        if (invert) {
            p <- ggplot(data = result$data, mapping = aes_(x = ~Ch2.Amplitude, y = ~Ch1.Amplitude))
        } else {
            p <- ggplot(data = result$data, mapping = aes_(x = ~Ch1.Amplitude, y = ~Ch2.Amplitude))
        }
        if (length(unique(result$data$Cluster)) > 8) {
            cbPalette <- c("#999999", "#f272e6", "#e5bdbe", "#bf0072", "#cd93c5", 
                "#1fba00", "#5e7f65", "#bdef00", "#2c5d26", "#ffe789", "#4a8c00", 
                "#575aef", "#a3b0fa", "#005caa", "#01c8fe", "#bc8775")
        } else if (length(unique(result$data$Cluster)) > 4) {
            cbPalette <- c("#999999", "#d800c4", "#fca3a7", "#bb004e", "#70cf56", 
                "#8b9d61", "#ccd451", "#bc8775")
        } else if (length(unique(result$data$Cluster)) > 2) {
            cbPalette <- c("#999999", "#8d5286", "#b42842", "#c2ff79", "#076633", 
                "#bc8775")
        } else {
            cbPalette <- c("#999999", "#bc8775")
        }
        if (invert) {
            p <- p + geom_point(aes(color = factor(Cluster)), size = 0.4, na.rm = TRUE) + 
                ggtitle(id) + theme_bw() + theme(legend.position = "none") + scale_colour_manual(values = cbPalette) + 
                labs(x = paste(annotations[3], "Amplitude"), y = paste(annotations[2], 
                  "Amplitude")) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
        } else {
            p <- p + geom_point(aes(color = factor(Cluster)), size = 0.4, na.rm = TRUE) + 
                ggtitle(id) + theme_bw() + theme(legend.position = "none") + scale_colour_manual(values = cbPalette) + 
                labs(x = paste(annotations[2], "Amplitude"), y = paste(annotations[3], 
                  "Amplitude")) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
        }
        ggsave(paste0(directory, "/", annotations[1], "/", id, ".", format), p, dpi = 1200)
    }
}

#' Export the algorithms results to an Excel file
#'
#' A convinience function that takes the results of the droplClust algorithm and exports them to an Excel file.
#'
#' @param data The result of the ddPCRclust algorithm
#' @param directory The parent directory where the files should saved. 
#' A new folder with the experiment name will be created (see below).
#' @param annotations Some basic metadata about the ddPCR reaction. 
#' If you provided \code{\link{ddPCRclust}} a template, 
#' this paramater can be filled with the corresponding field in the result.
#' Otherwise, you have to provide a character vector containing a name and the the color channels, 
#' e.g. \code{c(Name='ddPCR_01-04-2017', Ch1='HEX', Ch2='FAM')}
#' @param raw Boolean which determines if the annotated raw data should be exported along with the final counts. 
#' Basically, a third column will be added to the original data, 
#' which contains the cluster number to which this point was assigned to.
#' Useful for example to visualize the clustering later on. (Warning: this can take a while!)
#' @return None
#' @examples
#' # Read files
#' files <- readFiles(exampleFiles[1:8])
#' 
#' # Read template
#' template <- readTemplate(exampleFiles[9])
#' 
#' # Run ddPCRclust
#' result <- ddPCRclust(files, template)
#'
#' # Export the results
#' dir.create('./Results')
#' exportToExcel(data = result, directory = './Results/', annotations = result$annotations)
#'
exportToExcel <- function(data, directory, annotations, raw = FALSE) {
    directory <- normalizePath(directory, mustWork = TRUE)
    ifelse(!dir.exists(paste0(directory, "/", annotations[1])), dir.create(paste0(directory, 
        "/", annotations[1])), FALSE)
    dataToWrite <- NULL
    for (i in seq_along(data)) {
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
            dataToWrite <- merge(dataToWrite, tmp, all.x = TRUE)
            dataToWrite <- rbind(dataToWrite, tmp[match(colnames(dataToWrite), colnames(tmp))])
        } else {
            tmp <- merge(dataToWrite, tmp, all.y = TRUE)
            dataToWrite <- rbind(dataToWrite, tmp[match(colnames(dataToWrite), colnames(tmp))])
        }
        if (raw) {
            file <- paste0(directory, "/", annotations[1], "/", id, "_annotated_raw.xlsx")
            openxlsx::write.xlsx(result$data, file = file)
        }
    }
    mynames <- c("Empties", "1", "2", "3", "4", "1+2", "1+3", "1+4", "2+3", "2+4", 
        "3+4", "1+2+3", "1+2+4", "1+3+4", "2+3+4", "1+2+3+4", "Removed", "Total", 
        "Confidence")
    dataToWrite <- t(dataToWrite[, match(mynames, colnames(dataToWrite))])
    colnames(dataToWrite) <- names(data)
    file <- paste0(directory, "/", annotations[1], "/", annotations[1], "_results.xlsx")
    openxlsx::write.xlsx(dataToWrite, file = file, colNames = TRUE, rowNames = TRUE)
}

#' Export the algorithms results to a csv file
#'
#' A convinience function that takes the results of the droplClust algorithm and exports them to a csv file.
#'
#' @param data The result of the ddPCRclust algorithm
#' @param directory The parent directory where the files should saved. 
#' A new folder with the experiment name will be created (see below).
#' @param annotations Some basic metadata about the ddPCR reaction. 
#' If you provided \code{\link{ddPCRclust}} a template, 
#' this paramater can be filled with the corresponding field in the result.
#' Otherwise, you have to provide a character vector containing a name and the the color channels, 
#' e.g. \code{c(Name='ddPCR_01-04-2017', Ch1='HEX', Ch2='FAM')}
#' @param raw Boolean which determines if the annotated raw data should be exported along with the final counts. 
#' Basically, a third column will be added to the original data, 
#' which contains the cluster number to which this point was assigned to.
#' Useful for example to visualize the clustering later on. (Warning: this can take a while!)
#' @return None
#' @examples
#' # Read files
#' files <- readFiles(exampleFiles[1:8])
#' 
#' # Read template
#' template <- readTemplate(exampleFiles[9])
#' 
#' # Run ddPCRclust
#' result <- ddPCRclust(files, template)
#'
#' # Export the results
#' dir.create('./Results')
#' exportToCSV(data = result, directory = './Results/', annotations = result$annotations)
#'
exportToCSV <- function(data, directory, annotations, raw = FALSE) {
    directory <- normalizePath(directory, mustWork = TRUE)
    ifelse(!dir.exists(paste0(directory, "/", annotations[1])), dir.create(paste0(directory, 
        "/", annotations[1])), FALSE)
    dataToWrite <- NULL
    for (i in seq_along(data)) {
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
            dataToWrite <- merge(dataToWrite, tmp, all.x = TRUE)
            dataToWrite <- rbind(dataToWrite, tmp[match(colnames(dataToWrite), colnames(tmp))])
        } else {
            tmp <- merge(dataToWrite, tmp, all.y = TRUE)
            dataToWrite <- rbind(dataToWrite, tmp[match(colnames(dataToWrite), colnames(tmp))])
        }
        if (raw) {
            file <- paste0(directory, "/", annotations[1], "/", id, "_annotated_raw.csv")
            utils::write.csv(result$data, file = file)
        }
    }
    mynames <- c("Empties", "1", "2", "3", "4", "1+2", "1+3", "1+4", "2+3", "2+4", 
        "3+4", "1+2+3", "1+2+4", "1+3+4", "2+3+4", "1+2+3+4", "Removed", "Total", 
        "Confidence")
    dataToWrite <- t(dataToWrite[, match(mynames, colnames(dataToWrite))])
    colnames(dataToWrite) <- names(data)
    file <- paste0(directory, "/", annotations[1], "/", annotations[1], "_results.csv")
    utils::write.csv(dataToWrite, file = file)
}

# wrapper function for exception handling in mcmapply
dens_wrapper <- function(file, sensitivity = 1, numOfMarkers, markerNames, similarityParam, 
    distanceParam) {
    missingClusters <- which(markerNames == "")
    result <- tryCatch(expr = R.utils::withTimeout(runDensity(file[, c(2, 1)], sensitivity, 
        numOfMarkers, missingClusters, similarityParam, distanceParam), timeout = 60), 
        TimeoutException = function(ex) "TimedOut", error = function(e) print(e))
}

# wrapper function for exception handling in mcmapply
sam_wrapper <- function(file, sensitivity = 1, numOfMarkers, markerNames, similarityParam, 
    distanceParam) {
    missingClusters <- which(markerNames == "")
    result <- tryCatch(expr = R.utils::withTimeout(runSam(file[, c(2, 1)], sensitivity, 
        numOfMarkers, missingClusters, similarityParam, distanceParam), timeout = 60), 
        TimeoutException = function(ex) "TimedOut", error = function(e) print(e))
}

# wrapper function for exception handling in mcmapply
peaks_wrapper <- function(file, sensitivity = 1, numOfMarkers, markerNames, similarityParam, 
    distanceParam) {
    missingClusters <- which(markerNames == "")
    result <- tryCatch(expr = R.utils::withTimeout(runPeaks(file[, c(2, 1)], sensitivity, 
        numOfMarkers, missingClusters, similarityParam, distanceParam), timeout = 60), 
        TimeoutException = function(ex) "TimedOut", error = function(e) print(e))
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
    result <- tryCatch(createEnsemble(dens_result, sam_result, peaks_result, file[, 
        c(2, 1)]), error = function(e) {
        print(e)
    })
}


#' Find the clusters using flowDensity
#'
#' Use the local density function of the flowDensity package to find the cluster centres of the ddPCR reaction. 
#' Clusters are then labelled based on their rotated position and lastly the rain is assigned.
#'
#' @param file The input data. More specifically, a data frame with two dimensions, 
#' each dimension representing the intensity for one color channel.
#' @param sensitivity A number between 0.1 and 2 determining sensitivity of the initial clustering, 
#' e.g. the number of clusters. A higher value means more clusters are being found. Standard is 1.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up.
#' @param missingClusters A vector containing the number of primary clusters, 
#' which are missing in this dataset according to the template.
#' @param similarityParam If the distance of a droplet between two or more clusters is very similar, 
#' it will not be counted for either. The standard it 0.95, i.e. at least 95\% similarity. 
#' A sensible value lies between 0 and 1, where 0 means none of the 'rain' droplets will be counted and 
#' 1 means all droplets will be counted.
#' @param distanceParam When assigning rain between to clusters, 
#' typically the bottom 20\% are assigned to the lower cluster and remaining 80\% to the higher cluster. 
#' This parameter changes the ratio.
#' @return
#' \item{data}{The original input data minus the removed events (for plotting)}
#' \item{counts}{The droplet count for each cluster.}
#' \item{firstClusters}{The position of the primary clusters.}
#' \item{partition}{The cluster numbers as a CLUE partition (see clue package for more information).}
#' @export
#' @examples
#' # Run the flowDensity based approach
#' exampleFiles <- list.files(paste0(find.package('ddPCRclust'), '/extdata'), full.names = TRUE)
#' file <- read.csv(exampleFiles[3])
#' densResult <- runDensity(file = file, numOfMarkers = 4)
#'
#' # Plot the results
#' library(ggplot2)
#' p <- ggplot(data = densResult$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
#' p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = TRUE) +
#'      ggtitle('flowDensity example')+theme_bw() + theme(legend.position='none')
#' p
#'
runDensity <- function(file, sensitivity = 1, numOfMarkers, missingClusters = NULL, 
    similarityParam = 0.95, distanceParam = 0.2) {
    # ****** Parameters *******
    scalingParam <- c(max(file[, 1])/25, max(file[, 2])/25)
    epsilon <- 0.02/sensitivity^3
    # *************************
    
    if (!is.numeric(numOfMarkers) || numOfMarkers > 4 || numOfMarkers < 1) {
        stop("Invalid argument for numOfMarkers. Currently only the detection of 1-4 markers is supported.")
    } else if (!is.numeric(similarityParam) || similarityParam > 1 || similarityParam < 
        0) {
        stop("Invalid argument for similarityParam. Only values between 0 and 1 are supported.")
    } else if (!is.numeric(distanceParam) || distanceParam > 1 || distanceParam < 0) {
        stop("Invalid argument for distanceParam. Only values between 0 and 1 are supported.")
    } else if (!is.numeric(sensitivity) || sensitivity > 2 || sensitivity < 0.1) {
        stop("Invalid argument for sensitivity. Only values between 0.1 and 2 are supported.")
    }
    
    f <- flowCore::flowFrame(exprs = as.matrix(file))
    
    DataRemoved <- FinalResults <- NULL
    
    up1max <- deGate(f, c(1), percentile = 0.999, use.percentile = TRUE)
    up1min <- deGate(f, c(1), percentile = 0.001, use.percentile = TRUE)
    up2max <- deGate(f, c(2), percentile = 0.999, use.percentile = TRUE)
    up2min <- deGate(f, c(2), percentile = 0.001, use.percentile = TRUE)
    
    indices <- unique(c(which(flowCore::exprs(f)[, 1] >= 0.15 * (up1max - up1min) + 
        up1min), which(flowCore::exprs(f)[, 2] >= 0.15 * (up2max - up2min) + up2min)))
    
    f_remNeg <- f
    flowCore::exprs(f_remNeg) <- flowCore::exprs(f_remNeg)[indices, ]  
    # keep the non 15% bottom left corner (rn for removed neg)
    
    f_onlyNeg <- f
    flowCore::exprs(f_onlyNeg) <- flowCore::exprs(f_onlyNeg)[-indices, ]  
    # keep the     15% bottom left corner
    
    # coordinates of negative populations
    XcN <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_onlyNeg)[, 1], 
        width = 1000), tinypeak.removal = 0.2)$Peaks[1]
    YcN <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_onlyNeg)[, 2], 
        width = 1000), tinypeak.removal = 0.2)$Peaks[1]
    
    emptyDroplets <- c(XcN, YcN)
    firstClusters <- secondClusters <- tertClusters <- quadCluster <- NULL
    
    #---- find 1st gen clusters---------------------------------------------#
    
    firstClusters <- findPrimaryClustersDensity(f, file, f_remNeg, numOfMarkers, 
        scalingParam, epsilon)
    
    NumberOfSinglePos <- nrow(firstClusters$clusters)
    NumOfClusters <- 2^NumberOfSinglePos
    
    f_findExtremes_temp <- f
    
    x_leftPrim <- firstClusters$clusters[1, 1]
    y_leftPrim <- firstClusters$clusters[1, 2]
    x_rightPrim <- firstClusters$clusters[NumberOfSinglePos, 1]
    y_rightPrim <- firstClusters$clusters[NumberOfSinglePos, 2]
    
    mSlope <- (y_leftPrim - y_rightPrim)/(x_leftPrim - x_rightPrim)
    theta <- abs(atan(mSlope))
    rotate <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
    
    Rot_xy_leftPrim <- rotate %*% c(x_leftPrim, y_leftPrim)  # coordinates of rotated left  primary cluster
    Rot_xy_rightPrim <- rotate %*% c(x_rightPrim, y_rightPrim)  # coordinates of rotated right primary cluster
    
    flowCore::exprs(f_findExtremes_temp)[, c(1, 2)] <- t(rotate %*% t(flowCore::exprs(f_findExtremes_temp)[, 
        c(1, 2)]))
    upSlantmax <- deGate(f_findExtremes_temp, c(2), percentile = 0.999, use.percentile = TRUE)
    upSlantmin <- deGate(f_findExtremes_temp, c(2), percentile = 0.001, use.percentile = TRUE)
    ScaleChop <- (upSlantmax - upSlantmin)/max(file)
    
    #---- remove 1st gen clusters------------------------------------------------------------#
    
    f_temp <- f_remNeg
    
    
    for (o1 in seq_len(NumberOfSinglePos)) {
        indices <- union(which(flowCore::exprs(f_temp)[, 1] >= firstClusters$clusters[o1, 
            1] + ScaleChop * (scalingParam[1] + firstClusters$deviation[o1, 1])), 
            which(flowCore::exprs(f_temp)[, 2] >= firstClusters$clusters[o1, 2] + 
                ScaleChop * (scalingParam[2] + firstClusters$deviation[o1, 2])))
        flowCore::exprs(f_temp) <- flowCore::exprs(f_temp)[indices, ]
    }
    
    if (NumberOfSinglePos > 2) {
        #---- find 2nd gen clusters-----------------------------------------------------------#
        
        secondClusters <- findSecondaryClustersDensity(f, file, f_temp, emptyDroplets, 
            firstClusters, scalingParam, epsilon)
        
        #---- remove 2nd gen clusters---------------------------------------------------------#
        
        for (o1 in seq_len(nrow(secondClusters$clusters))) {
            indices <- union(which(flowCore::exprs(f_temp)[, 1] >= secondClusters$clusters[o1, 
                1] + ScaleChop * (scalingParam[1] + secondClusters$deviation[o1, 
                1])), which(flowCore::exprs(f_temp)[, 2] >= secondClusters$clusters[o1, 
                2] + ScaleChop * (scalingParam[2] + secondClusters$deviation[o1, 
                2])))
            if (length(indices) <= 1) {
                flowCore::exprs(f_temp) <- flowCore::exprs(f_temp)[c(indices, indices), 
                  ]
                next
            }
            flowCore::exprs(f_temp) <- flowCore::exprs(f_temp)[indices, ]
        }
        flowCore::exprs(f_temp)[, c(1, 2)] <- t(rotate %*% t(flowCore::exprs(f_temp)[, 
            c(1, 2)]))
    }
    
    if (NumberOfSinglePos > 3) {
        #---- find 3rd gen clusters-------------------------------------------------------#
        
        tertClusters <- findTertiaryClustersDensity(f, f_temp, emptyDroplets, firstClusters, 
            secondClusters, scalingParam, epsilon)
        
        #---- remove 3rd gen clusters-----------------------------------------------------#
        
        flowCore::exprs(f_temp)[, c(1, 2)] <- t(t(rotate) %*% t(flowCore::exprs(f_temp)[, 
            c(1, 2)]))
        for (o1 in seq_len(4)) {
            indices <- union(which(flowCore::exprs(f_temp)[, 1] >= tertClusters$clusters[o1, 
                1] + ScaleChop * (scalingParam[1] + tertClusters$deviation[o1, 1])), 
                which(flowCore::exprs(f_temp)[, 2] >= tertClusters$clusters[o1, 2] + 
                  ScaleChop * (scalingParam[2] + tertClusters$deviation[o1, 2])))
            if (length(indices) <= 1) {
                flowCore::exprs(f_temp) <- flowCore::exprs(f_temp)[c(indices, indices), 
                  ]
                next
            }
            flowCore::exprs(f_temp) <- flowCore::exprs(f_temp)[indices, ]
        }
    }
    
    if (NumberOfSinglePos > 1) {
        #---- find the 4th gen cluster----------------------------------------------------------#
        
        quadCluster <- findQuaternaryClusterDensity(f, f_temp, emptyDroplets, firstClusters, 
            secondClusters, tertClusters)
    }
    
    #-------------------------------------------------------------------------------------------#
    
    # 
    posOfFirsts <- 2:(1 + nrow(firstClusters$clusters))
    ClusterCentres <- rbind(emptyDroplets, firstClusters$clusters)
    posOfSeconds <- posOfThirds <- posOfFourth <- NULL
    if (length(secondClusters$clusters)) {
        posOfSeconds <- (utils::tail(posOfFirsts, 1) + 1):(utils::tail(posOfFirsts, 
            1) + nrow(secondClusters$clusters))
        ClusterCentres <- rbind(ClusterCentres, secondClusters$clusters)
    }
    if (length(tertClusters$clusters)) {
        posOfThirds <- (utils::tail(posOfSeconds, 1) + 1):(utils::tail(posOfSeconds, 
            1) + nrow(tertClusters$clusters))
        ClusterCentres <- rbind(ClusterCentres, tertClusters$clusters)
    }
    if (length(quadCluster$clusters)) {
        posOfFourth <- max((utils::tail(posOfFirsts, 1) + 1), (utils::tail(posOfSeconds, 
            1) + 1), (utils::tail(posOfThirds, 1) + 1), na.rm = TRUE)
        ClusterCentres <- rbind(ClusterCentres, abs(quadCluster$clusters))
    }
    
    angles <- vapply(seq_len(nrow(firstClusters$clusters)), function(x) return(atan2(firstClusters$clusters[x, 
        1] - emptyDroplets[1], firstClusters$clusters[x, 2] - emptyDroplets[2])), 
        double(length = 1))
    cuts <- c(0, 0.5 * pi/4, 1.5 * pi/4, pi/2)
    for (i in missingClusters) {
        angles <- c(angles, cuts[i])
    }
    distMatrix <- t(vapply(seq_along(angles), function(x) return(abs(angles[x] - 
        cuts)), double(length = length(cuts))))
    order <- solve_LSAP(distMatrix)
    deletions <- which(!seq_len(4) %in% order)
    deletions <- c(deletions, missingClusters)
    # deletions <- deletions[!deletions %in% missingClusters] order <-
    # order[1:nrow(firstClusters$clusters)]
    names <- c("Empties", "1", "2", "3", "4", "1+2", "1+3", "1+4", "2+3", "2+4", 
        "3+4", "1+2+3", "1+2+4", "1+3+4", "2+3+4", "1+2+3+4", "Removed", "Total")
    names_indices <- seq_along(names)
    indices <- vector()
    for (cl in deletions) {
        indices <- c(indices, grep(cl, names))
    }
    if (length(indices)) 
        names_indices <- names_indices[-indices]
    
    result <- rep(0, nrow(f))
    newData <- flowCore::exprs(f)
    rownames(newData) <- seq_len(nrow(f))
    for (cl in seq_len(nrow(ClusterCentres))) {
        result[as.numeric(rownames(newData[(newData[, 1] < ClusterCentres[cl, 1] + 
            scalingParam[1] & newData[, 1] > ClusterCentres[cl, 1] - scalingParam[1] & 
            newData[, 2] < ClusterCentres[cl, 2] + scalingParam[2] & newData[, 2] > 
            ClusterCentres[cl, 2] - scalingParam[2]), ]))] <- cl
        newData <- subset(newData, !(newData[, 1] < ClusterCentres[cl, 1] + scalingParam[1] & 
            newData[, 1] > ClusterCentres[cl, 1] - scalingParam[1] & newData[, 2] < 
            ClusterCentres[cl, 2] + scalingParam[2] & newData[, 2] > ClusterCentres[cl, 
            2] - scalingParam[2]))
    }
    
    for (i in seq_len(nrow(newData))) {
        temp <- apply(ClusterCentres, 1, function(x) {
            euc.dist(x, newData[i, ])
        })
        result[as.numeric(rownames(newData)[i])] <- which.min(temp)
    }
    
    ClusterCentresNew <- t(vapply(seq_len(NumOfClusters), function(x) return(colMeans(flowCore::exprs(f)[result == 
        x, , drop = FALSE])), double(length = 2)))
    ClusterCentresNew[which(is.nan(ClusterCentresNew))] <- ClusterCentres[which(is.nan(ClusterCentresNew))]
    
    fDensResult <- assignRain(clusterMeans = ClusterCentresNew, data = flowCore::exprs(f), 
        result = result, emptyDroplets = 1, firstClusters = posOfFirsts, secondClusters = posOfSeconds, 
        thirdClusters = posOfThirds, fourthCluster = posOfFourth, scalingParam = scalingParam, 
        similarityParam = similarityParam, distanceParam = distanceParam)
    
    fDensResult$result[fDensResult$result == 0] <- 0/0
    if (NumberOfSinglePos < 4) {
        tempResult <- fDensResult$result
        for (i in seq_len(nrow(ClusterCentres))) {
            fDensResult$result[which(tempResult == i)] <- names_indices[i]
        }
    } else {
        tempResult <- fDensResult$result
        fDensResult$result[which(tempResult == 8)] <- 9
        fDensResult$result[which(tempResult == 9)] <- 8
    }
    
    removed <- c(fDensResult$removed, which(is.nan(fDensResult$result)))
    
    NumOfEventsClust <- table(c(fDensResult$result, seq_len(length(names) - 2))) - 
        1
    NumOfEventsClust <- c(NumOfEventsClust, length(removed))  # add on removed
    NumOfEventsClust <- c(NumOfEventsClust, sum(NumOfEventsClust))  # add on total
    names(NumOfEventsClust) = names
    
    if (length(removed)) {
        fDensResult$result[removed] <- length(names) - 1  # remove the removed ones
    }
    result <- cbind(file[, c(2, 1)], Cluster = fDensResult$result)
    partition <- as.cl_partition(c(fDensResult$result, seq_len(length(names) - 1)))
    return(list(data = result, counts = NumOfEventsClust, firstClusters = firstClusters$clusters, 
        partition = partition))
}


#' Find the clusters using SamSPECTRAL
#'
#' Find the rain and assign it based on the distance to vector lines connecting the cluster centres.
#'
#' @param file The input data. More specifically, a data frame with two dimensions, 
#' each dimension representing the intensity for one color channel.
#' @param sensitivity A number between 0.1 and 2 determining sensitivity of the initial clustering, 
#' e.g. the number of clusters. A higher value means more clusters are being found. Standard is 1.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up.
#' @param missingClusters A vector containing the number of primary clusters, 
#' which are missing in this dataset according to the template.
#' @param similarityParam If the distance of a droplet between two or more clusters is very similar, 
#' it will not be counted for either. The standard it 0.95, i.e. at least 95\% similarity. 
#' A sensible value lies between 0 and 1, where 0 means none of the 'rain' droplets will be counted and 
#' 1 means all droplets will be counted.
#' @param distanceParam When assigning rain between to clusters, 
#' typically the bottom 20\% are assigned to the lower cluster and remaining 80\% to the higher cluster. 
#' This parameter changes the ratio.
#' @return
#' \item{data}{The original input data minus the removed events (for plotting)}
#' \item{counts}{The droplet count for each cluster.}
#' \item{firstClusters}{The position of the primary clusters.}
#' \item{partition}{The cluster numbers as a CLUE partition (see clue package for more information).}
#' @export
#' @examples
#' # Run the SamSPECTRAL based approach
#' exampleFiles <- list.files(paste0(find.package('ddPCRclust'), '/extdata'), full.names = TRUE)
#' file <- read.csv(exampleFiles[3])
#' samResult <- runSam(file = file, numOfMarkers = 4)
#'
#' # Plot the results
#' library(ggplot2)
#' p <- ggplot(data = samResult$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
#' p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = TRUE) +
#'      ggtitle('SamSPECTRAL example')+theme_bw() + theme(legend.position='none')
#' p
#'
runSam <- function(file, sensitivity = 1, numOfMarkers, missingClusters = NULL, similarityParam = 0.95, 
    distanceParam = 0.2) {
    # ****** Parameters *******
    scalingParam <- c(max(file[, 1])/25, max(file[, 2])/25)
    epsilon <- 0.02/sensitivity^3
    m <- trunc(nrow(file)/20)
    # *************************
    
    if (!is.numeric(numOfMarkers) || numOfMarkers > 4 || numOfMarkers < 1) {
        stop("Invalid argument for numOfMarkers. Currently only the detection of 1-4 markers is supported.")
    } else if (!is.numeric(similarityParam) || similarityParam > 1 || similarityParam < 
        0) {
        stop("Invalid argument for similarityParam. Only values between 0 and 1 are supported.")
    } else if (!is.numeric(distanceParam) || distanceParam > 1 || distanceParam < 0) {
        stop("Invalid argument for distanceParam. Only values between 0 and 1 are supported.")
    } else if (!is.numeric(sensitivity) || sensitivity > 2 || sensitivity < 0.1) {
        stop("Invalid argument for sensitivity. Only values between 0.1 and 2 are supported.")
    }
    
    f <- flowCore::flowFrame(exprs = as.matrix(file))
    
    #### Start of algorithm ####
    samRes <- SamSPECTRAL(data.points = as.matrix(file), dimensions = c(1, 2), normal.sigma = (200 * 
        sensitivity^2), separation.factor = (0.88 * sensitivity), m = m, talk = FALSE)
    data <- file
    temp <- lapply(min(samRes, na.rm = TRUE):max(samRes, na.rm = TRUE), function(x) return(apply(data[samRes == 
        x, ], 2, stats::median, na.rm = TRUE)))
    clusterMeans <- t(do.call(cbind, temp))
    rowSums <- vapply(seq_len(nrow(clusterMeans)), function(x) return(sum(clusterMeans[x, 
        ])), double(length = 1))
    emptyDroplets <- match(min(rowSums), rowSums)
    dimensions <- c(max(data[1]), max(data[2]))
    temp <- vapply(min(samRes, na.rm = TRUE):max(samRes, na.rm = TRUE), function(x) return(abs(stats::var(data[samRes == 
        x, 1], data[samRes == x, 2], na.rm = TRUE))), double(length = 1))
    badClusters <- match(temp[temp > sum(dimensions) * 25], temp)
    samTable <- table(samRes)
    secondaryClusters <- tertiaryClusters <- quaternaryCluster <- NULL
    firstClusters <- findPrimaryClusters(samRes, clusterMeans, emptyDroplets, badClusters, 
        dimensions, file, f, numOfMarkers, scalingParam, epsilon)
    ## estimate missing clusters based on angle:
    angles <- vapply(seq_along(firstClusters), function(x) return(atan2(clusterMeans[firstClusters[x], 
        1] - clusterMeans[emptyDroplets, 1], clusterMeans[firstClusters[x], 2] - 
        clusterMeans[emptyDroplets, 2])), double(length = 1))
    cuts <- c(0, 0.5 * pi/4, 1.5 * pi/4, pi/2)
    for (i in missingClusters) {
        angles <- c(angles, cuts[i])
    }
    distMatrix <- t(vapply(seq_along(angles), function(x) return(abs(angles[x] - 
        cuts)), double(length = length(cuts))))
    order <- solve_LSAP(distMatrix)
    deletions <- which(!seq_len(4) %in% order)
    deletions <- c(deletions, missingClusters)
    # deletions <- deletions[!deletions %in% missingClusters] order <-
    # order[1:nrow(firstClusters$clusters)]
    names <- c("Empties", "1", "2", "3", "4", "1+2", "1+3", "1+4", "2+3", "2+4", 
        "3+4", "1+2+3", "1+2+4", "1+3+4", "2+3+4", "1+2+3+4", "Removed", "Total")
    names_indices <- seq_along(names)
    indices <- vector()
    for (cl in deletions) {
        indices <- c(indices, grep(cl, names))
    }
    if (length(indices)) 
        names_indices <- names_indices[-indices]
    
    if (length(firstClusters) == 1) {
        samResult <- c(emptyDroplets, firstClusters)
    }
    if (length(firstClusters) == 2) {
        quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, 
            firstClusters)
        samResult <- c(emptyDroplets, firstClusters, quaternaryCluster)
    }
    if (length(firstClusters) == 3) {
        secondaryClusters <- findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, 
            badClusters, sum(dimensions), samTable)
        quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, 
            firstClusters, secondaryClusters$clusters)
        samResult <- c(emptyDroplets, firstClusters, secondaryClusters$clusters, 
            quaternaryCluster)
        # if (length(firstClusters) == numOfMarkers) { names <-
        # c('1','2','3','1+2','1+3','2+3','1+2+3','Empties','Removed','Total') } else {
        # missingCluster <- findDeletion(clusterMeans, firstClusters, emptyDroplets,
        # numOfMarkers-length(firstClusters), dimensions) names <-
        # c('1','2','3','4','1+2','1+3','1+4','2+3','2+4','3+4','1+2+3','1+2+4','1+3+4','2+3+4','1+2+3+4','Empties','Removed','Total')
        # pos <- grep(missingCluster, names) for (foo in pos) { newsamResult <-
        # c(samResult, 0) id <- c( seq_along(samResult), foo-0.5 ) samResult <-
        # newsamResult[order(id)] } }
    }
    if (length(firstClusters) == 4) {
        secondaryClusters <- findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, 
            badClusters, sum(dimensions), samTable)
        tertiaryClusters <- findTertiaryClusters(emptyDroplets, firstClusters, secondaryClusters$clusters, 
            badClusters, clusterMeans, secondaryClusters$correctionFactor, sum(dimensions), 
            samTable)
        quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, 
            firstClusters, secondaryClusters$clusters, tertiaryClusters$clusters)
        samResult <- c(emptyDroplets, firstClusters, secondaryClusters$clusters, 
            tertiaryClusters$clusters, quaternaryCluster)
    }
    samRes <- mergeClusters(samRes, clusterMeans, samResult, badClusters)
    rain <- assignRain(clusterMeans, data, samRes, emptyDroplets, firstClusters, 
        secondaryClusters$clusters, tertiaryClusters$clusters, quaternaryCluster, 
        scalingParam, similarityParam, distanceParam)
    samRes <- rain$result
    firstClusters <- clusterMeans[firstClusters, ]
    
    for (missingCluster in deletions) {
        firstClusters <- insertRow(firstClusters, cbind(0, 0), missingCluster)
    }
    finalSamRes <- rep(0/0, length(samRes))
    for (i in seq_along(samResult)) {
        finalSamRes[which(samRes == samResult[i])] <- names_indices[i]
    }
    clusterCount <- table(c(finalSamRes, seq_len(length(names) - 2))) - 1
    removed <- c(rain$removed, which(is.nan(finalSamRes)))
    if (length(removed)) {
        finalSamRes[removed] <- length(names) - 1
    }
    
    samCount <- c(clusterCount, length(removed))
    samCount <- c(samCount, sum(samCount))
    names(samCount) = names
    result <- cbind(data[, c(2, 1)], Cluster = finalSamRes)
    partition <- as.cl_partition(c(finalSamRes, seq_len(length(names) - 1)))
    return(list(data = result, counts = samCount, firstClusters = firstClusters, 
        partition = partition))
}


#' Find the clusters using flowPeaks
#'
#' Find the rain and assign it based on the distance to vector lines connecting the cluster centres.
#'
#' @param file The input data. More specifically, a data frame with two dimensions, 
#' each dimension representing the intensity for one color channel.
#' @param sensitivity A number between 0.1 and 2 determining sensitivity of the initial clustering, 
#' e.g. the number of clusters. A higher value means more clusters are being found. Standard is 1.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up.
#' @param missingClusters A vector containing the number of primary clusters, 
#' which are missing in this dataset according to the template.
#' @param similarityParam If the distance of a droplet between two or more clusters is very similar, 
#' it will not be counted for either. The standard it 0.95, i.e. at least 95\% similarity. 
#' A sensible value lies between 0 and 1, where 0 means none of the 'rain' droplets will be counted 
#' and 1 means all droplets will be counted.
#' @param distanceParam When assigning rain between to clusters, 
#' typically the bottom 20\% are assigned to the lower cluster and remaining 80\% to the higher cluster. 
#' This parameter changes the ratio.
#' @return
#' \item{data}{The original input data minus the removed events (for plotting)}
#' \item{counts}{The droplet count for each cluster.}
#' \item{firstClusters}{The position of the primary clusters.}
#' \item{partition}{The cluster numbers as a CLUE partition (see clue package for more information).}
#' @export
#' @examples
#' # Run the flowPeaks based approach
#' exampleFiles <- list.files(paste0(find.package('ddPCRclust'), '/extdata'), full.names = TRUE)
#' file <- read.csv(exampleFiles[3])
#' peaksResult <- runPeaks(file = file, numOfMarkers = 4)
#'
#' # Plot the results
#' library(ggplot2)
#' p <- ggplot(data = peaksResult$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
#' p <- p + geom_point(aes(color = factor(Cluster)), size = .5, na.rm = TRUE) +
#'      ggtitle('flowPeaks example')+theme_bw() + theme(legend.position='none')
#' p
#'
runPeaks <- function(file, sensitivity = 1, numOfMarkers, missingClusters = NULL, 
    similarityParam = 0.95, distanceParam = 0.2) {
    # ****** Parameters *******
    scalingParam <- c(max(file[, 1])/25, max(file[, 2])/25)
    epsilon <- 0.02/sensitivity^3
    # *************************
    
    if (!is.numeric(numOfMarkers) || numOfMarkers > 4 || numOfMarkers < 1) {
        stop("Invalid argument for numOfMarkers. Currently only the detection of 1-4 markers is supported.")
    } else if (!is.numeric(similarityParam) || similarityParam > 1 || similarityParam < 
        0) {
        stop("Invalid argument for similarityParam. Only values between 0 and 1 are supported.")
    } else if (!is.numeric(distanceParam) || distanceParam > 1 || distanceParam < 0) {
        stop("Invalid argument for distanceParam. Only values between 0 and 1 are supported.")
    } else if (!is.numeric(sensitivity) || sensitivity > 2 || sensitivity < 0.1) {
        stop("Invalid argument for sensitivity. Only values between 0.1 and 2 are supported.")
    }
    
    f <- flowCore::flowFrame(exprs = as.matrix(file))
    
    #### Start of algorithm ####
    fPeaksRes <- flowPeaks(file, tol = 0, h0 = (0.3/sensitivity^1.5), h = (0.5/sensitivity^1.5))
    data <- file
    clusterMeans <- fPeaksRes$peaks$mu
    rowSums <- vapply(seq_len(nrow(clusterMeans)), function(x) return(sum(clusterMeans[x, 
        ])), double(length = 1))
    emptyDroplets <- match(min(rowSums), rowSums)
    dimensions <- c(max(data[1]), max(data[2]))
    temp <- vapply(min(fPeaksRes$peaks.cluster, na.rm = TRUE):max(fPeaksRes$peaks.cluster, 
        na.rm = TRUE), function(x) return(abs(stats::var(data[fPeaksRes$peaks.cluster == 
        x, 1], data[fPeaksRes$peaks.cluster == x, 2], na.rm = TRUE))), double(length = 1))
    badClusters <- match(temp[temp > sum(dimensions) * 25], temp)
    fPeaksTable <- table(fPeaksRes$peaks.cluster)
    secondaryClusters <- tertiaryClusters <- quaternaryCluster <- NULL
    firstClusters <- findPrimaryClusters(fPeaksRes$peaks.cluster, clusterMeans, emptyDroplets, 
        badClusters, dimensions, file, f, numOfMarkers, scalingParam, epsilon)
    ## estimate missing clusters based on angle:
    angles <- vapply(seq_along(firstClusters), function(x) return(atan2(clusterMeans[firstClusters[x], 
        1] - clusterMeans[emptyDroplets, 1], clusterMeans[firstClusters[x], 2] - 
        clusterMeans[emptyDroplets, 2])), double(length = 1))
    cuts <- c(0, 0.5 * pi/4, 1.5 * pi/4, pi/2)
    for (i in missingClusters) {
        angles <- c(angles, cuts[i])
    }
    distMatrix <- t(vapply(seq_along(angles), function(x) return(abs(angles[x] - 
        cuts)), double(length = length(cuts))))
    order <- solve_LSAP(distMatrix)
    deletions <- which(!seq_len(4) %in% order)
    deletions <- c(deletions, missingClusters)
    # deletions <- deletions[!deletions %in% missingClusters] order <-
    # order[1:nrow(firstClusters$clusters)]
    names <- c("Empties", "1", "2", "3", "4", "1+2", "1+3", "1+4", "2+3", "2+4", 
        "3+4", "1+2+3", "1+2+4", "1+3+4", "2+3+4", "1+2+3+4", "Removed", "Total")
    names_indices <- seq_along(names)
    indices <- vector()
    for (cl in deletions) {
        indices <- c(indices, grep(cl, names))
    }
    if (length(indices)) 
        names_indices <- names_indices[-indices]
    
    if (length(firstClusters) == 1) {
        fPeaksResult <- c(emptyDroplets, firstClusters)
    } else if (length(firstClusters) == 2) {
        quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, 
            firstClusters)
        fPeaksResult <- c(emptyDroplets, firstClusters, quaternaryCluster)
    } else if (length(firstClusters) == 3) {
        secondaryClusters <- findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, 
            badClusters, sum(dimensions), fPeaksTable)
        quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, 
            firstClusters, secondaryClusters$clusters)
        fPeaksResult <- c(emptyDroplets, firstClusters, secondaryClusters$clusters, 
            quaternaryCluster)
    } else if (length(firstClusters) == 4) {
        secondaryClusters <- findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, 
            badClusters, sum(dimensions), fPeaksTable)
        tertiaryClusters <- findTertiaryClusters(emptyDroplets, firstClusters, secondaryClusters$clusters, 
            badClusters, clusterMeans, secondaryClusters$correctionFactor, sum(dimensions), 
            fPeaksTable)
        quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, 
            firstClusters, secondaryClusters$clusters, tertiaryClusters$clusters)
        fPeaksResult <- c(emptyDroplets, firstClusters, secondaryClusters$clusters, 
            tertiaryClusters$clusters, quaternaryCluster)
    }
    
    fPeaksRes$peaks.cluster <- mergeClusters(fPeaksRes$peaks.cluster, clusterMeans, 
        fPeaksResult, badClusters)
    rain <- assignRain(clusterMeans, data, fPeaksRes$peaks.cluster, emptyDroplets, 
        firstClusters, secondaryClusters$clusters, tertiaryClusters$clusters, quaternaryCluster, 
        scalingParam, similarityParam, distanceParam)
    fPeaksRes$peaks.cluster <- rain$result
    firstClusters <- clusterMeans[firstClusters, ]
    
    for (missingCluster in deletions) {
        firstClusters <- insertRow(firstClusters, cbind(0, 0), missingCluster)
    }
    finalPeaksRes <- rep(0/0, length(fPeaksRes$peaks.cluster))
    for (i in seq_along(fPeaksResult)) {
        finalPeaksRes[which(fPeaksRes$peaks.cluster == fPeaksResult[i])] <- names_indices[i]
    }
    clusterCount <- table(c(finalPeaksRes, seq_len(length(names) - 2))) - 1
    removed <- c(rain$removed, which(is.nan(finalPeaksRes)))
    if (length(removed)) {
        finalPeaksRes[removed] <- length(names) - 1
    }
    
    fPeaksCount <- c(clusterCount, length(removed))
    fPeaksCount <- c(fPeaksCount, sum(fPeaksCount))
    names(fPeaksCount) = names
    result <- cbind(data[, c(2, 1)], Cluster = finalPeaksRes)
    partition <- as.cl_partition(c(finalPeaksRes, seq_len(length(names) - 1)))
    return(list(data = result, counts = fPeaksCount, firstClusters = firstClusters, 
        partition = partition))
}

#' Calculates the copies per droplet
#'
#' This function takes the results of the clustering and calculates the actual counts per target, 
#' as well as the counts per droplet (CPD) for each marker.
#'
#' @param results The result of the ddPCRclust algorithm.
#' @param template The parsed dataframe containing the template.
#' @param constantControl The constant refrence control, which should be present in each reaction. 
#' It is used to normalize the data.
#' @return
#' A list of lists, containing the counts for empty droplets, each marker with both total droplet count and CPD, 
#' and total number of droplets, for each element of the input list respectively.
#' @export
#' @examples
#' # Read files
#' files <- readFiles(exampleFiles[1:8])
#' 
#' # Read template
#' template <- readTemplate(exampleFiles[9])
#' 
#' # Run ddPCRclust
#' result <- ddPCRclust(files, template)
#'
#' # Calculate the CPDs
#' markerCPDs <- calculateCPDs(result, result$template)
#'
calculateCPDs <- function(results, template = NULL, constantControl = NULL) {
    countedResult <- list()
    maxctrl <- 0
    for (i in seq_along(results)) {
        id <- names(results[i])
        result <- results[[i]]$counts
        if (is.null(result)) 
            next
        markers <- as.character(unlist(template[which(template[, 1] == id), 4:7]))
        if (!is.null(constantControl)) {
            sctrl <- which(markers == constantControl)
            if (length(sctrl) != 1) {
                stop("The constant control does not seem to be valid!")
            }
        }
        total <- as.integer(result[grep("Total", names(result))])
        empties <- as.integer(result[grep("Empties", names(result))])
        for (j in seq_len(4)) {
            if (is.null(template)) {
                marker <- paste0("M", j)
            } else {
                marker <- markers[j]
            }
            if (marker == "") 
                next
            counts <- as.integer(result[grep(j, names(result))])
            counts <- sum(counts, na.rm = TRUE)
            if (total == 0) {
                cpd <- 0
            } else {
                cpd <- -log(1 - (counts/total))
            }
            if (!is.null(constantControl)) {
                if (j == sctrl && cpd > maxctrl) {
                  maxctrl <- cpd
                }
            }
            countedResult[[id]][[marker]] <- list(counts = counts, cpd = cpd)
        }
        countedResult[[id]][["Empties"]] <- empties
        countedResult[[id]][["Total"]] <- total
    }
    # Normalize to constant control
    if (!is.null(constantControl)) {
        for (i in seq_along(countedResult)) {
            id <- names(countedResult[i])
            modFactor <- maxctrl/countedResult[[id]][[constantControl]]$cpd
            markers <- as.character(unlist(template[which(template[, 1] == id), 4:7]))
            for (j in markers) {
                if (j == "") 
                  next
                countedResult[[id]][[j]]$cpd <- countedResult[[id]][[j]]$cpd * modFactor
            }
        }
    }
    return(countedResult)
}

#' Create a cluster ensemble
#'
#' This function takes the three (or less) clustering approaches of the ddPCRclust package 
#' and combines them to one cluster ensemble. See \link{cl_medoid} for more information.
#'
#' @param dens The result of the flowDensity algorithm as a CLUE partition.
#' @param sam The result of the samSPECTRAL algorithm as a CLUE partition.
#' @param peaks The result of the flowPeaks algorithm as a CLUE partition.
#' @param file The input data. More specifically, a data frame with two dimensions, 
#' each dimension representing the intensity for one color.
#' @return
#' \item{data}{The original input data minus the removed events (for plotting)}
#' \item{confidence}{The agreement between the different clustering results in percent. 
#' If all algorithms calculated the same result, the clustering is likely to be correct, thus the confidence is high.}
#' \item{counts}{The droplet count for each cluster.}
#' @export
#' @examples
#' exampleFiles <- list.files(paste0(find.package('ddPCRclust'), '/extdata'), full.names = TRUE)
#' file <- read.csv(exampleFiles[3])
#' densResult <- runDensity(file = file, numOfMarkers = 4)
#' samResult <- runSam(file = file, numOfMarkers = 4)
#' peaksResult <- runPeaks(file = file, numOfMarkers = 4)
#'
#' superResult <- createEnsemble(densResult, samResult, peaksResult, file)
#'
createEnsemble <- function(dens = NULL, sam = NULL, peaks = NULL, file) {
    listResults <- list()
    if (!is.null(dens$partition)) {
        names <- names(dens$counts)
        listResults <- c(listResults, list(dens$partition))
    }
    if (!is.null(sam$partition)) {
        names <- names(sam$counts)
        listResults <- c(listResults, list(sam$partition))
    }
    if (!is.null(peaks$partition)) {
        names <- names(peaks$counts)
        listResults <- c(listResults, list(peaks$partition))
    }
    
    if (!length(listResults)) {
        next
    } else if (length(listResults) == 1) {
        comb_ids <- cl_class_ids(listResults[[1]])
        conf <- 1
    } else {
        cens <- cl_ensemble(list = listResults)
        conf <- mean(cl_agreement(cens, method = "cRand"))
        comb <- cl_medoid(cens)
        comb_ids <- cl_class_ids(comb)
    }
    
    superCounts <- table(comb_ids) - 1
    superCounts <- c(superCounts, sum(superCounts))
    names(superCounts) <- names
    comb_ids[comb_ids == length(names) - 1] <- 0/0  ## Remove the removed ones
    superResult <- cbind(file, Cluster = as.integer(comb_ids[seq_len(nrow(file))]))
    return(list(data = superResult, confidence = conf, counts = superCounts, nrOfAlgorithms = length(listResults)))
}

#' Correct for DNA shearing
#'
#' Longer DNA templates produce a lower droplet count due to DNA shearing.
#' This function normalizes the ddPCRclust result based on a stable marker of different lengths 
#' to negate the effect of differences in the lengths of the actual markers of interest.
#' (Work in progress)
#'
#' @param counts The counts per marker as provided by \link{calculateCPDs}.
#' @param lengthControl The name of the length Control. If the template name is for example CPT2, 
#' the name in the template should be CPT2-125, where 125 represents the number of basepairs.
#' @param stableControl The name of the stable Control used as a reference for this experiment.
#' @return
#' A linear regression model fitting the length vs ln(ratio) (see \link{lm} for details on linear regression).
#'
shearCorrection <- function(counts, lengthControl, stableControl) {
    controlRatios <- data.frame()
    for (well in counts) {
        markerPos = grep(paste0("^", lengthControl), names(well))
        controlPos = grep(paste0("^", stableControl), names(well))
        if (!length(markerPos) || !length(controlPos) || 
            !length(well[[controlPos]]$cpd) || well[[controlPos]]$cpd == 0) 
            next
        ratio <- log(well[[markerPos]]$cpd/well[[controlPos]]$cpd)
        len <- unlist(strsplit(names(well[markerPos]), "-"))[2]
        if (is.na(len)) {
            warning(paste("Bad marker name", names(well[markerPos])))
            next
        }
        controlRatios <- rbind(controlRatios, as.numeric(c(len, ratio)))
    }
    colnames(controlRatios) <- c("Length", "Ratio")
    stats::lm(Ratio ~ Length, controlRatios)
}

## Part of the dropClust algorithm
## Author: Benedikt G Brink, Bielefeld University
## April 2017

#' @include cluster_functions.R
#' @include functions.R

source("R/cluster_functions.R")
source("R/functions.R")

library(flowDensity)
library(SamSPECTRAL)
library(flowPeaks)
library(plotrix)
library(clue)


#' Find the clusters using flowDensity
#'
#' Use the local density function of the flowDensity package to find the cluster centres of the ddPCR reaction. Clusters are then labelled based on their rotated position and lastly the rain is assigned.
#'
#' @param file The input data. More specifically, a data frame with two dimensions, each dimension representing the intensity for one color channel.
#' @param sensitivity An integer between 0.1 and 2 determining sensitivity of the initial clustering, e.g. the number of clusters. A higher value means more clusters are being found. Standard is 1.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up.
#' @return
#' \item{data}{The original input data minus the removed events (for plotting)}
#' \item{counts}{The droplet count for each cluster.}
#' \item{firstClusters}{The position of the primary clusters.}
#' \item{partition}{The cluster numbers as a CLUE partition (see clue package for more information).}
#' @export
#' @examples
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = TRUE)
#' file <- read.csv(exampleFiles[3])
#' densResult <- runDensity(file = file, numOfMarkers = 4)
#' 
#' library(ggplot2)
#' p <- ggplot(data = densResult$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
#' p <- p+geom_point(aes(color = factor(Cluster)), size = .5, na.rm = T)+ggtitle("flowDensity example")+theme_bw()+theme(legend.position="none")
#' p
#' 
runDensity <- function(file, sensitivity=1, numOfMarkers) {
  
  # ****** Parameters *******
  scalingParam <<- c(max(file[,1])/25, max(file[,2])/25)
  CutAbovePrimary <<- mean(scalingParam)*2
  epsilon <<- 0.02/sensitivity^3
  threshold <<- 0.1/sensitivity^2
  # *************************
  
  data_dir <- system.file("extdata", package = "flowDensity")
  load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE)) # load f to copy over later so we have an FCS file to use flowDensity
  
  f@exprs <- as.matrix(file) # overide the FCS file. This allows us to use flowDensity on data that is not truely a FCS data.
  #file <- file[,c(2,1)] # switch axis to be consistent
  
  DataRemoved <- FinalResults <- NULL
  
  up1max <- deGate(f, c(1), percentile=0.999, use.percentile=T)
  up1min <- deGate(f, c(1), percentile=0.001, use.percentile=T)
  up2max <- deGate(f, c(2), percentile=0.999, use.percentile=T)
  up2min <- deGate(f, c(2), percentile=0.001, use.percentile=T)
  
  indices <- unique( c( which(f@exprs[,1] >= 0.15 * (up1max - up1min) + up1min),
                        which(f@exprs[,2] >= 0.15 * (up2max - up2min) + up2min) ) )
  
  f_remNeg  <- f; f_remNeg @exprs <- f_remNeg @exprs[ indices,] # keep the non 15% bottom left corner (rn for removed neg)
  f_onlyNeg <- f; f_onlyNeg@exprs <- f_onlyNeg@exprs[-indices,] # keep the     15% bottom left corner
  
  # coordinates of negative populations
  XcN <- flowDensity:::.getPeaks(density(f_onlyNeg@exprs[,1], width=1000), tinypeak.removal=0.2)$Peaks[1]
  YcN <- flowDensity:::.getPeaks(density(f_onlyNeg@exprs[,2], width=1000), tinypeak.removal=0.2)$Peaks[1]
  
  emptyDroplets <- c(XcN, YcN)
  firstClusters <- secondClusters <- tertClusters <- quadCluster <- NULL
  
  #---- find 1st gen clusters------------------------------------------------------------------------------------------------------------------------#
  
  firstClusters <- findPrimaryClustersDensity(f, file, f_remNeg, numOfMarkers)
  
  NumberOfSinglePos <- nrow(firstClusters$clusters)
  NumOfClusters <- 2^NumberOfSinglePos
  
  #---- remove 1st gen clusters------------------------------------------------------------------------------------------------------------------------#
  
  f_temp <- f_remNeg;
  
  for ( o1 in 1:NumberOfSinglePos ) {
    indices <- union(which(f_temp@exprs[,1] >= firstClusters$clusters[o1,1]+ScaleChop*(scalingParam[1]+firstClusters$deviation[o1,1])), 
                     which(f_temp@exprs[,2] >= firstClusters$clusters[o1,2]+ScaleChop*(scalingParam[2]+firstClusters$deviation[o1,2])))
    f_temp@exprs <- f_temp@exprs[indices,]
  }
  
  if (NumberOfSinglePos > 2) {
    #---- find 2nd gen clusters------------------------------------------------------------------------------------------------------------------------#
    
    secondClusters <- findSecondaryClustersDensity(f, f_temp, emptyDroplets, firstClusters)
    
    #---- remove 2nd gen clusters----------------------------------------------------------------------------------------------------------------------#
    
    for ( o1 in 1:nrow(secondClusters$clusters) ) {
      indices <- union(which(f_temp@exprs[,1] >= secondClusters$clusters[o1, 1]+ScaleChop*(scalingParam[1]+secondClusters$deviation[o1, 1])), 
                       which(f_temp@exprs[,2] >= secondClusters$clusters[o1, 2]+ScaleChop*(scalingParam[2]+secondClusters$deviation[o1, 2])))
      if (length(indices) <= 1) {
        f_temp@exprs <- f_temp@exprs[c(indices, indices),]
        next
      }
      f_temp@exprs <- f_temp@exprs[indices,]
    }
    f_temp@exprs[,c(1,2)] <- t(R %*% t(f_temp@exprs[,c(1,2)]))
  }
  
  if (NumberOfSinglePos > 3) {
    #---- find 3rd gen clusters------------------------------------------------------------------------------------------------------------------------#
    
    tertClusters <- findTertiaryClustersDensity(f, f_temp, emptyDroplets, firstClusters, secondClusters)
    
    #---- remove 3rd gen clusters----------------------------------------------------------------------------------------------------------------------#
    
    f_temp@exprs[,c(1,2)] <- t(t(R) %*% t(f_temp@exprs[,c(1,2)]))
    for ( o1 in 1:4 ) {
      indices <- union(which(f_temp@exprs[,1] >= tertClusters$clusters[o1,1]+ScaleChop*(scalingParam[1]+tertClusters$deviation[o1, 1])), 
                       which(f_temp@exprs[,2] >= tertClusters$clusters[o1,2]+ScaleChop*(scalingParam[2]+tertClusters$deviation[o1, 2])))
      if (length(indices) <= 1) {
        f_temp@exprs <- f_temp@exprs[c(indices, indices),]
        next
      }
      f_temp@exprs <- f_temp@exprs[indices,]
    }
  }
  
  if (NumberOfSinglePos > 1) {
    #---- find the 4th gen cluster------------------------------------------------------------------------------------------------------------------------#
    
    quadCluster <- findQuaternaryClusterDensity(f, f_temp, emptyDroplets, firstClusters, secondClusters, tertClusters)
  }
  
  #----------------------------------------------------------------------------------------------------------------------------------------------------------#
  
  #
  posOfFirsts <- 2:(1+nrow(firstClusters$clusters))
  ClusterCentres <- rbind(emptyDroplets, firstClusters$clusters)
  posOfSeconds <- posOfThirds <- posOfFourth <- NULL
  if (length(secondClusters$clusters) > 0) {
    posOfSeconds <- (tail(posOfFirsts, 1)+1):(tail(posOfFirsts, 1)+nrow(secondClusters$clusters))
    ClusterCentres <- rbind(ClusterCentres, secondClusters$clusters)
  }
  if (length(tertClusters$clusters) > 0) {
    posOfThirds <- (tail(posOfSeconds, 1)+1):(tail(posOfSeconds, 1)+nrow(tertClusters$clusters))
    ClusterCentres <- rbind(ClusterCentres, tertClusters$clusters)
  }
  if (length(quadCluster$clusters) > 0) {
    posOfFourth <- max((tail(posOfFirsts, 1)+1), (tail(posOfSeconds, 1)+1), (tail(posOfThirds, 1)+1), na.rm = T)
    ClusterCentres <- rbind(ClusterCentres, abs(quadCluster$clusters))
  }
  
  result <- rep(0, nrow(f))
  newData <- f@exprs
  rownames(newData) <- 1:nrow(f)
  for (c in 1:nrow(ClusterCentres)) {
    result[as.numeric(rownames(newData[(newData[,1] < ClusterCentres[c,1]+scalingParam[1] & newData[,1] > ClusterCentres[c,1]-scalingParam[1]
                                        & newData[,2] < ClusterCentres[c,2]+scalingParam[2] & newData[,2] > ClusterCentres[c,2]-scalingParam[2]),]))] <- c
    newData <- subset(newData, !(newData[,1] < ClusterCentres[c,1]+scalingParam[1] & newData[,1] > ClusterCentres[c,1]-scalingParam[1]
                                 & newData[,2] < ClusterCentres[c,2]+scalingParam[2] & newData[,2] > ClusterCentres[c,2]-scalingParam[2]))
  }
  
  for (i in 1:nrow(newData)) {
    temp <- apply(ClusterCentres,1, function(x){euc.dist(x, newData[i,])})
    result[as.numeric(rownames(newData)[i])] <- which.min(temp)
  }
  
  ClusterCentresNew <- t(sapply(1:NumOfClusters, function(x) return(colMeans(f@exprs[result == x, , drop=F]))))
  ClusterCentresNew[which(is.nan(ClusterCentresNew))] <- ClusterCentres[which(is.nan(ClusterCentresNew))]
  
  
  fDensResult <- assignRain(clusterMeans = ClusterCentresNew, data = f@exprs, result = result, emptyDroplets = 1, firstClusters = posOfFirsts, secondClusters = posOfSeconds, thirdClusters = posOfThirds, fourthCluster = posOfFourth, flowDensity = T)
  
  if (NumberOfSinglePos == 1) {
    names <- c("1","Empties","Removed","Total")
  } else if (NumberOfSinglePos == 2) {
    names <- c("1","2","1+2","Empties","Removed","Total")
  } else if (NumberOfSinglePos == 3) {
    names <- c("1","2","3","1+2","1+3","2+3","1+2+3","Empties","Removed","Total")
  } else if (NumberOfSinglePos == 4) {
    tempResult <- fDensResult$result
    fDensResult$result[which(tempResult == 8)] <- 9
    fDensResult$result[which(tempResult == 9)] <- 8
    names <- c("1","2","3","4","1+2","1+3","1+4","2+3","2+4","3+4","1+2+3","1+2+4","1+3+4","2+3+4","1+2+3+4","Empties","Removed","Total")
  }
  
  NumOfEventsClust <- matrix(0,NumOfClusters,1)
  for ( t1 in 1:NumOfClusters) {
    NumOfEventsClust[t1] <- length(which(fDensResult$result == t1 ))
  }
  
  NumOfEventsClust <- c(NumOfEventsClust[2:length(NumOfEventsClust)],NumOfEventsClust[1],length(fDensResult$remove))
  
  NumOfEventsClust <- c(NumOfEventsClust, sum(NumOfEventsClust)) # add on total
  #   NumOfEventsClust[c(7,8)] <- NumOfEventsClust[c(8,7)]
  
  names(NumOfEventsClust) = names
  
  #     newResult <- rep(0/0, length(result))
  # 
  #     for (i in 1:16) {
  #       newResult[which(result == clusterOrder[i])] <- i
  #     }
  
  if ( length(fDensResult$removed) > 0 ) {
    fDensResult$result[fDensResult$remove] <- NumOfClusters+1# remove the removed ones
  }
  result <- cbind(file, "Cluster" = fDensResult$result)
  return(list(data=result, counts=NumOfEventsClust, firstClusters=firstClusters$clusters, partition=as.cl_partition(fDensResult$result)))
}


#' Find the clusters using SamSPECTRAL
#'
#' Find the rain and assign it based on the distance to vector lines connecting the cluster centres.
#'
#' @param file The input data. More specifically, a data frame with two dimensions, each dimension representing the intensity for one color channel.
#' @param sensitivity An integer between 0.1 and 2 determining sensitivity of the initial clustering, e.g. the number of clusters. A higher value means more clusters are being found. Standard is 1.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up.
#' @return
#' \item{data}{The original input data minus the removed events (for plotting)}
#' \item{counts}{The droplet count for each cluster.}
#' \item{firstClusters}{The position of the primary clusters.}
#' \item{partition}{The cluster numbers as a CLUE partition (see clue package for more information).}
#' @export
#' @import SamSPECTRAL flowDensity clue
#' @examples
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = TRUE)
#' file <- read.csv(exampleFiles[3])
#' samResult <- runSam(file = file, numOfMarkers = 4)
#' 
#' library(ggplot2)
#' p <- ggplot(data = samResult$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
#' p <- p+geom_point(aes(color = factor(Cluster)), size = .5, na.rm = T)+ggtitle("SamSPECTRAL example")+theme_bw()+theme(legend.position="none")
#' p
#' 
runSam <- function(file, sensitivity = 1, numOfMarkers) {
  
  # ****** Parameters *******
  scalingParam <<- c(max(file[,1])/25, max(file[,2])/25)
  CutAbovePrimary <<- mean(scalingParam)*2
  epsilon <<- 0.02/sensitivity^3
  threshold <<- 0.1/sensitivity^2
  m <- trunc(nrow(file)/20)
  # *************************
  
  data_dir <- system.file("extdata", package = "flowDensity")
  load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE)) # load f to copy over later so we have an FCS file to use flowDensity
  
  f@exprs <- as.matrix(file) # overide the FCS file. This allows us to use flowDensity on data that is not truely a FCS data.

  #### Start of algorithm ####
  samRes <- SamSPECTRAL(data.points=as.matrix(file),dimensions=c(1,2), normal.sigma = (400*sensitivity^2), separation.factor = (0.88*sensitivity), m = m, talk=F)
  data <- file
  temp <- lapply(min(samRes, na.rm = T):max(samRes, na.rm = T), function(x) return(apply(data[samRes==x,], 2, median, na.rm = T)))
  clusterMeans <- t(do.call(cbind, temp))
  rowSums <- sapply(1:nrow(clusterMeans), function(x) return(sum(clusterMeans[x,])))
  emptyDroplets <- match(min(rowSums), rowSums)
  dimensions <- c(max(data[1]), max(data[2]))
  temp <- sapply(min(samRes, na.rm = T):max(samRes, na.rm = T), function(x) return(abs(var(data[samRes==x,1], data[samRes==x,2], na.rm = T))))
  badClusters <- match(temp[temp>sum(dimensions)*25], temp)
  samTable <- table(samRes)
  #firstClusters <- findPrimaryClusters(samRes, clusterMeans, rowSums, emptyDroplets, badClusters, dimensions)
  missingCluster <- 0
  secondaryClusters <- tertiaryClusters <- quaternaryCluster <- NULL
  firstClusters <- findPrimaryClusters(samRes, clusterMeans, emptyDroplets, badClusters, dimensions, file, f, numOfMarkers)
  if(length(firstClusters) == 1) {
    samResult <- c(emptyDroplets, firstClusters)
    names <- c("1","Empties","Removed","Total")
  }
  if(length(firstClusters) == 2) {
    quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, firstClusters)
    samResult <- c(emptyDroplets, firstClusters, quaternaryCluster)
    names <- c("1","2","1+2","Empties","Removed","Total")
  }
  if(length(firstClusters) == 3) {
    secondaryClusters <- findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, badClusters, sum(dimensions), samTable)
    quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, firstClusters, secondaryClusters$clusters)
    samResult <- c(emptyDroplets, firstClusters, secondaryClusters$clusters, quaternaryCluster)
    if (length(firstClusters) == numOfMarkers) {
      names <- c("1","2","3","1+2","1+3","2+3","1+2+3","Empties","Removed","Total")
    } else {
      missingCluster <- findDeletion(clusterMeans, firstClusters, emptyDroplets, numOfMarkers-length(firstClusters), dimensions)
      names <- c("1","2","3","4","1+2","1+3","1+4","2+3","2+4","3+4","1+2+3","1+2+4","1+3+4","2+3+4","1+2+3+4","Empties","Removed","Total")
      pos <- grep(missingCluster, names)
      for (foo in pos) {
        newsamResult <- c(samResult, 0)
        id  <- c( seq_along(samResult), foo-0.5 )
        samResult <- newsamResult[order(id)]
      }
    }
  }
  if(length(firstClusters) == 4) {
    secondaryClusters <- findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, badClusters, sum(dimensions), samTable)
    tertiaryClusters <- findTertiaryClusters(emptyDroplets, firstClusters, secondaryClusters$clusters, badClusters, clusterMeans, secondaryClusters$correctionFactor, sum(dimensions), samTable)
    quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, firstClusters, secondaryClusters$clusters, tertiaryClusters$clusters)
    samResult <- c(emptyDroplets, firstClusters, secondaryClusters$clusters, tertiaryClusters$clusters, quaternaryCluster)
    names <- c("1","2","3","4","1+2","1+3","1+4","2+3","2+4","3+4","1+2+3","1+2+4","1+3+4","2+3+4","1+2+3+4","Empties","Removed","Total")
    
  }
  samRes <- mergeClusters(samRes, clusterMeans, samResult, badClusters)
  rain <- assignRain(clusterMeans, data, samRes, emptyDroplets, firstClusters, secondaryClusters$clusters, tertiaryClusters$clusters, quaternaryCluster, F)
  samRes <- rain$result
  samTable <- table(samRes)
  clusterCount <- samTable[as.character(samResult)]
  clusterCount[which(is.na(clusterCount))] = 0
  clusterCount <- c(clusterCount[2:length(clusterCount)], clusterCount[1])
  
  firstClusters <- clusterMeans[firstClusters,]
  if (missingCluster > 0) {
    firstClusters <- insertRow(firstClusters, cbind(0,0), missingCluster)
  }
  
  finalSamRes <- rep(0/0, length(samRes))
  for (i in 1:length(samResult)) {
    finalSamRes[which(samRes == samResult[i])] <- i
  }
  
  removed <- which(is.nan(finalSamRes))
  if (length(removed) > 0) {
    finalSamRes[removed] <- 2^numOfMarkers+1
  }
  

  samCount <- c(clusterCount, length(removed), length(file[,1]))
  names(samCount) = names
  result <- cbind(data, "Cluster" = finalSamRes)
  return(list(data=result, counts=samCount, firstClusters=firstClusters, partition=as.cl_partition(finalSamRes)))
}


#' Find the clusters using flowPeaks
#'
#' Find the rain and assign it based on the distance to vector lines connecting the cluster centres.
#'
#' @param file The input data. More specifically, a data frame with two dimensions, each dimension representing the intensity for one color channel.
#' @param sensitivity An integer between 0.1 and 2 determining sensitivity of the initial clustering, e.g. the number of clusters. A higher value means more clusters are being found. Standard is 1.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up.
#' @return
#' \item{data}{The original input data minus the removed events (for plotting)}
#' \item{counts}{The droplet count for each cluster.}
#' \item{firstClusters}{The position of the primary clusters.}
#' \item{partition}{The cluster numbers as a CLUE partition (see clue package for more information).}
#' @export
#' @examples
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = TRUE)
#' file <- read.csv(exampleFiles[3])
#' peaksResult <- runPeaks(file = file, numOfMarkers = 4)
#' 
#' library(ggplot2)
#' p <- ggplot(data = peaksResult$data, mapping = aes(x = Ch2.Amplitude, y = Ch1.Amplitude))
#' p <- p+geom_point(aes(color = factor(Cluster)), size = .5, na.rm = T)+ggtitle("flowPeaks example")+theme_bw()+theme(legend.position="none")
#' p
#' 
runPeaks <- function(file, sensitivity = 1, numOfMarkers) {
  
  # ****** Parameters *******
  scalingParam <<- c(max(file[,1])/25, max(file[,2])/25)
  CutAbovePrimary <<- mean(scalingParam)*2
  epsilon <<- 0.02/sensitivity^3
  threshold <<- 0.1/sensitivity^2
  # *************************
  
  data_dir <- system.file("extdata", package = "flowDensity")
  load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE)) # load f to copy over later so we have an FCS file to use flowDensity
  
  f@exprs <- as.matrix(file) # overide the FCS file. This allows us to use flowDensity on data that is not truely a FCS data.
  
  #### Start of algorithm ####
  fPeaksRes <- flowPeaks(file, tol=0, h0=(0.3/sensitivity^1.5), h=(0.4/sensitivity^1.5))
  data <- file
  clusterMeans <- fPeaksRes$peaks$mu
  rowSums <- sapply(1:nrow(clusterMeans), function(x) return(sum(clusterMeans[x,])))
  emptyDroplets <- match(min(rowSums), rowSums)
  dimensions <- c(max(data[1]), max(data[2]))
  temp <- sapply(min(fPeaksRes$peaks.cluster, na.rm = T):max(fPeaksRes$peaks.cluster, na.rm = T), function(x) return(abs(var(data[fPeaksRes$peaks.cluster==x,1], data[fPeaksRes$peaks.cluster==x,2], na.rm = T))))
  badClusters <- match(temp[temp>sum(dimensions)*25], temp)
  fPeaksTable <- table(fPeaksRes$peaks.cluster)
  #firstClusters <- findPrimaryClusters(fPeaksRes$peaks.cluster, clusterMeans, rowSums, emptyDroplets, badClusters, dimensions)
  missingCluster <- 0
  secondaryClusters <- tertiaryClusters <- quaternaryCluster <- NULL
  firstClusters <- findPrimaryClusters(fPeaksRes$peaks.cluster, clusterMeans, emptyDroplets, badClusters, dimensions, file, f, numOfMarkers)
  if(length(firstClusters) == 1) {
    fPeaksResult <- c(emptyDroplets, firstClusters)
    names <- c("1","Empties","Removed","Total")
  } else
    if(length(firstClusters) == 2) {
      quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, firstClusters)
      fPeaksResult <- c(emptyDroplets, firstClusters, quaternaryCluster)
      names <- c("1","2","1+2","Empties","Removed","Total")
    } else
      if(length(firstClusters) == 3) {
        secondaryClusters <- findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, badClusters, sum(dimensions), fPeaksTable)
        quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, firstClusters, secondaryClusters$clusters)
        fPeaksResult <- c(emptyDroplets, firstClusters, secondaryClusters$clusters, quaternaryCluster)
        if (length(firstClusters) == numOfMarkers) {
          names <- c("1","2","3","1+2","1+3","2+3","1+2+3","Empties","Removed","Total")
        } else {
          missingCluster <- findDeletion(clusterMeans, firstClusters, emptyDroplets, numOfMarkers-length(firstClusters), dimensions)
          # firstClusters <- append(firstClusters, 0, after=missingCluster-1)
          names <- c("1","2","3","4","1+2","1+3","1+4","2+3","2+4","3+4","1+2+3","1+2+4","1+3+4","2+3+4","1+2+3+4","Empties","Removed","Total")
          pos <- grep(missingCluster, names)
          for (foo in pos) {
            newfPeaksResult <- c(fPeaksResult, 0)
            id  <- c( seq_along(fPeaksResult), foo-0.5 )
            fPeaksResult <- newfPeaksResult[order(id)]
          }
        }
      } else
        if(length(firstClusters) == 4) {
          secondaryClusters <- findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, badClusters, sum(dimensions), fPeaksTable)
          tertiaryClusters <- findTertiaryClusters(emptyDroplets, firstClusters, secondaryClusters$clusters, badClusters, clusterMeans, secondaryClusters$correctionFactor, sum(dimensions), fPeaksTable)
          quaternaryCluster <- findQuaternaryCluster(clusterMeans, emptyDroplets, badClusters, firstClusters, secondaryClusters$clusters, tertiaryClusters$clusters)
          fPeaksResult <- c(emptyDroplets, firstClusters, secondaryClusters$clusters, tertiaryClusters$clusters, quaternaryCluster)
          names <- c("1","2","3","4","1+2","1+3","1+4","2+3","2+4","3+4","1+2+3","1+2+4","1+3+4","2+3+4","1+2+3+4","Empties","Removed","Total")
        }
  
  fPeaksRes$peaks.cluster <- mergeClusters(fPeaksRes$peaks.cluster, clusterMeans, fPeaksResult, badClusters)
  rain <- assignRain(clusterMeans, data, fPeaksRes$peaks.cluster, emptyDroplets, firstClusters, secondaryClusters$clusters, tertiaryClusters$clusters, quaternaryCluster, F)
  fPeaksRes$peaks.cluster <- rain$result
  fPeaksTable <- table(fPeaksRes$peaks.cluster)
  clusterCount <- fPeaksTable[as.character(fPeaksResult)]
  clusterCount[which(is.na(clusterCount))] = 0
  clusterCount <- c(clusterCount[2:length(clusterCount)], clusterCount[1])
  firstClusters <- clusterMeans[firstClusters,]
  
  if (missingCluster > 0) {
    firstClusters <- insertRow(firstClusters, cbind(0,0), missingCluster)
  }
  finalPeaksRes <- rep(0/0, length(fPeaksRes$peaks.cluster))
  for (i in 1:length(fPeaksResult)) {
    finalPeaksRes[which(fPeaksRes$peaks.cluster == fPeaksResult[i])] <- i
  }
  removed <- which(is.nan(finalPeaksRes))
  if (length(removed) > 0) {
    finalPeaksRes[removed] <- 2^numOfMarkers+1
  }
  
  fPeaksCount <- c(clusterCount, length(removed), length(file[,1]))
  names(fPeaksCount) = names
  result <- cbind(data, "Cluster" = finalPeaksRes)
  return(list(data=result, counts=fPeaksCount, firstClusters=firstClusters, partition=as.cl_partition(finalPeaksRes)))
}

#' Calculates the copies per droplet
#'
#' This function takes the results of the clustering and calculates the actual counts per marker, as well as the counts per droplet (CPD) for each marker. 
#'
#' @param results The result of the dropClust algorithm.
#' @param template A parsed dataframe containing the template.
#' @return
#' A list of lists, containing the counts for empty droplets, each marker with both total droplet count and CPD, and total number of droplets, for each element of the input list respectively. 
#' @export
#' @examples
#' # Run dropClust
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = TRUE)
#' result <- runDropClust(files = exampleFiles[1:8], template = exampleFiles[9])
#' 
#' # Calculate the CPDs
#' markerCPDs <- calculateCPDs(result$results)
#'
calculateCPDs <- function(results, template = NULL) {
  countedResult <- list()

  for (i in 1:length(results)) {
    id <- names(results[i])
    result <- results[[i]]$counts
    total <- sum(as.integer(result[grep('Total', names(result))]), na.rm = T)
    empties <- sum(as.integer(result[grep('Empties', names(result))]), na.rm = T)
    for (j in 1:4) {
      counts <- sum(as.integer(result[grep(j, names(result))]), na.rm = T)
      cpd <- -log(1-counts/total)
      if (is.null(template)) {
        marker <- paste0("M", j)
      } else {
        marker <- as.character(template[which(template[,1] == id), j+3])
      }
      countedResult[[id]][[marker]] <- list(counts=counts, cpd=cpd)
    }
    countedResult[[id]][["Empties"]] <- empties
    countedResult[[id]][["Total"]] <- total
  }
  return(countedResult)
}

#' Create a cluster ensemble.
#'
#' This function takes the three (or less) clustering approaches of the dropClust package and combines them to one cluster ensemble. See \link{cl_medoid} for more information.
#'
#' @param dens The result of the flowDensity algorithm as a CLUE partition.
#' @param sam The result of the samSPECTRAL algorithm as a CLUE partition.
#' @param peaks The result of the flowPeaks algorithm as a CLUE partition.
#' @param numOfMarkers The number of primary clusters that are expected according the experiment set up.
#' @param file The input data. More specifically, a data frame with two dimensions, each dimension representing the intensity for one color.
#' @return
#' \item{data}{The original input data minus the removed events (for plotting)}
#' \item{confidence}{The agreement between the different clustering results in percent. If all algorithms calculated the same result, the clustering is likely to be correct, thus the confidence is high.}
#' \item{counts}{The droplet count for each cluster.}
#' @export
#' @import clue
#' @examples
#' exampleFiles <- list.files(paste0(find.package("dropClust"), "/extdata"), full.names = TRUE)
#' file <- read.csv(exampleFiles[3])
#' densResult <- runDensity(file = file, numOfMarkers = 4)
#' samResult <- runSam(file = file, numOfMarkers = 4)
#' peaksResult <- runPeaks(file = file, numOfMarkers = 4)
#' 
#' superResult <- createEnsemble(densResult, samResult, peaksResult, 4, file)
#' 
createEnsemble <- function(dens = NULL, sam = NULL, peaks = NULL, numOfMarkers, file) {

  listResults <- list()
  if (!is.null(dens$partition)) {
    listResults <- c(listResults, list(dens$partition))
  }
  if (!is.null(sam$partition)) {
    listResults <- c(listResults, list(sam$partition))
  }  
  if (!is.null(peaks$partition)) {
    listResults <- c(listResults, list(peaks$partition))
  }
  
  if (length(listResults) == 0) {
    next
  } else if (length(listResults) == 1) {
    comb_ids <- cl_class_ids(listResults[[1]])
    conf <- 1
  } else {
    cens <- cl_ensemble(list = listResults)
    conf <- mean(cl_agreement(cens))
    comb <- cl_medoid(cens)
    comb_ids <- cl_class_ids(comb)
  }
  temp_comb_ids <- comb_ids
  for (i in 1:2^numOfMarkers+1) {
    temp_comb_ids <- c(temp_comb_ids, i)
  }
  superCounts <- table(temp_comb_ids)-1
  superCounts <- c(superCounts, sum(superCounts))
  switch(as.character(numOfMarkers),
         '1' = names(superCounts) <- c("Empties","1","Removed","Total"),
         '2' = names(superCounts) <- c("Empties","1","2","1+2","Removed","Total"),
         '3' = names(superCounts) <- c("Empties","1","2","3","1+2","1+3","2+3","1+2+3","Removed","Total"),
         '4' = names(superCounts) <- c("Empties","1","2","3","4","1+2","1+3","1+4","2+3","2+4","3+4","1+2+3","1+2+4","1+3+4","2+3+4","1+2+3+4","Removed","Total"))
  
  comb_ids[comb_ids==2^numOfMarkers+1] <- 0/0
  superResult <- cbind(file, "Cluster" = as.integer(comb_ids))
  return(list(data=superResult, confidence=conf, counts=superCounts))
}


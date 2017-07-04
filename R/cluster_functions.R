## Part of the dropClust algorithm
## Author: Benedikt G Brink, Bielefeld University
## July 2017

# Find the primary clusters based on their density peaks found by flowDensity.
findPrimaryClustersDensity <- function(f, File, f_remNeg, NumOfMarkers) {
  
  # Calculate the threshold to remove events that are too positive. This makes it easier to find single positives.
  f_findExtremes_temp <- f
  f_remNeg_temp <- f_remNeg
  
  
  # find left and right primary
  tinyP <- 0.6
  repeat {
    tinyP <- tinyP/2
    if (tinyP < epsilon) break
    leftPrim <-  deGate(f_remNeg_temp, 1, all.cut=T, tinypeak.removal=tinyP, adjust=0.1)[1]
    indices <- which(f_remNeg_temp@exprs[,1] <= leftPrim)
    f_leftPrim <- f_remNeg_temp; f_leftPrim@exprs <- f_leftPrim@exprs[indices,]
    x_leftPrim <- max(flowDensity:::.getPeaks(density(f_leftPrim@exprs[,1], width=1000), tinypeak.removal=tinyP)$Peaks)
    y_leftPrim <- max(flowDensity:::.getPeaks(density(f_leftPrim@exprs[,2], width=1000), tinypeak.removal=tinyP)$Peaks)
    if (x_leftPrim < (min(File[,1]) + (max(File[,1])-min(File[,1]))/6)) {
      threshold <- tinyP
      break
    }
  }
  
  tinyP <- 0.6
  repeat{
    tinyP <- tinyP/2
    if (tinyP < epsilon) break
    rightPrim <- deGate(f_remNeg_temp, 2, all.cut=T, tinypeak.removal=tinyP, adjust=0.1)[1]
    indices <- which(f_remNeg_temp@exprs[,2] <= rightPrim)
    f_rightPrim <- f_remNeg_temp; f_rightPrim@exprs <- f_rightPrim@exprs[indices,]
    x_rightPrim <- max(flowDensity:::.getPeaks(density(f_rightPrim@exprs[,1], width=1000), tinypeak.removal=tinyP)$Peaks)
    y_rightPrim <- max(flowDensity:::.getPeaks(density(f_rightPrim@exprs[,2], width=1000), tinypeak.removal=tinyP)$Peaks)
    if (y_rightPrim < (min(File[,2]) + (max(File[,2])-min(File[,2]))/6)) {
      threshold <- min(threshold, tinyP)
      break
    }
  }
  
  mSlope <<- (y_leftPrim - y_rightPrim)/(x_leftPrim - x_rightPrim)
  theta <<- abs(atan(mSlope))
  R <<- matrix( c(cos(theta), sin(theta), -sin(theta), cos(theta)) ,2 ,2)
  
  Rot_xy_leftPrim <<- R %*% c(x_leftPrim, y_leftPrim) # coordinates of rotated left  primary cluster
  Rot_xy_rightPrim <<- R %*% c(x_rightPrim, y_rightPrim) # coordinates of rotated right primary cluster
  
  f_findExtremes_temp@exprs[,c(1,2)] <- t(R %*% t(f_findExtremes_temp@exprs[,c(1,2)]))
  f_remNeg_temp@exprs    [,c(1,2)] <- t(R %*% t(f_remNeg_temp@exprs    [,c(1,2)]))
  upSlantmax <<- deGate(f_findExtremes_temp, c(2), percentile=0.999, use.percentile=T)
  upSlantmin <<- deGate(f_findExtremes_temp, c(2), percentile=0.001, use.percentile=T)
  ScaleChop <<-  (upSlantmax - upSlantmin) / max(File)
  # Chop off the very positive events
  indices <- which(f_remNeg_temp@exprs[,2] <= (max(Rot_xy_leftPrim[2], Rot_xy_rightPrim[2]) + CutAbovePrimary*ScaleChop))
  f_remNeg_temp@exprs <- f_remNeg_temp@exprs[indices,]
  
  # find thresholds to divide the primary clusters
  highP <- 1
  lowP <- 0
  repeat { # newton iteration
    newP <- (highP+lowP)/2
    if (newP < 2*epsilon) break
    Cuts <- deGate(f_remNeg_temp, 1, all.cut=T, tinypeak.removal=newP, adjust=0.1)
    if (length(Cuts) < (NumOfMarkers - 1)){
      highP <- newP
    } else if (length(Cuts) > (NumOfMarkers - 1)) {
      lowP <- newP
    } else {
      break
    }
  }
  Xc <- Yc <- NULL
  Cuts <- c(min(f_remNeg_temp@exprs[,1]), Cuts, max(f_remNeg_temp@exprs[,1]))
  deviationX <- vector()
  deviationY <- vector()
  
  # find the coordinates of the primary clusters
  for ( r1 in 1:(length(Cuts)-1) ) {
    indices <- intersect(which(f_remNeg_temp@exprs[,1] >= Cuts[r1]), which(f_remNeg_temp@exprs[,1] < Cuts[r1+1]))
    f_Cuts_temp <- f_remNeg_temp; f_Cuts_temp@exprs <- f_Cuts_temp@exprs[indices,]
    f_Cuts_temp@exprs[,c(1,2)] <- t(t(R) %*% t(f_Cuts_temp@exprs[,c(1,2)]))
    Xx <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,1], width=1000), tinypeak.removal=newP)$Peaks
    Yy <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,2], width=1000), tinypeak.removal=newP)$Peaks
    if (length(Xx) > 1 || length(Yy) > 1) {
      highP <- 1
      lowP <- 0
      repeat { # newton iteration
        newP <- (highP+lowP)/2
        if (newP < epsilon) break
        Xx <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,1], width=1000), tinypeak.removal=newP)$Peaks
        Yy <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,2], width=1000), tinypeak.removal=newP)$Peaks
        if (length(Xx) > 1 || length(Yy) > 1){
          highP <- newP
        } else if (length(Xx) < 1 || length(Yy) < 1) {
          lowP <- newP
        } else {
          break
        }
      }
    }
    Xx <- Xx[1]
    Yy <- Yy[1]
    deviationX <- c(deviationX, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),1]))
    deviationY <- c(deviationY, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),2]))
    Xc <- c(Xc, Xx)
    Yc <- c(Yc, Yy)
  }
  if (length(Xc) == 4 && length(Yc) == 4) {
    Rot_xy_midLeftPrim <- R %*% c(Xc[2], Yc[2]) # coordinates of rotated middle left  primary cluster
    Rot_xy_midRightPrim <- R %*% c(Xc[3], Yc[3]) # coordinates of rotated middle right primary cluster
    
    test_left <- abs(det(rbind(cbind(1,1,1), cbind(Rot_xy_leftPrim, Rot_xy_midLeftPrim, rbind(0,0)))))/100
    test_right <- abs(det(rbind(cbind(1,1,1), cbind(Rot_xy_rightPrim, Rot_xy_midRightPrim, rbind(0,0)))))/100
    
    if (!is.na(test_left) && test_left < max(File)/2 && test_left > max(File)/20) {
      Xc[1] <- mean(Xc[1:2])
      Yc[1] <- mean(Yc[1:2])
      Xc <- Xc[-2]
      Yc <- Yc[-2]
    }
    if (!is.na(test_right) && test_right < max(File)/2 && test_right > max(File)/20) {
      Xc[3] <- mean(Xc[3:4])
      Yc[3] <- mean(Yc[3:4])
      Xc <- Xc[-4]
      Yc <- Yc[-4]
    }
  }
  return(list(clusters=cbind(Xc, Yc), deviation = cbind(deviationX, deviationY)))
}

# Find the secondary clusters based on their density peaks found by flowDensity.
findSecondaryClustersDensity <- function(f, f_remNegPrim, emptyDroplets, firstClusters) {
  f_remNegPrim_chopTertQuat <- f_remNegPrim
  NumberOfSinglePos <- nrow(firstClusters$clusters)
  f_remNegPrim_chopTertQuat@exprs    [,c(1,2)] <- t(R %*% t(f_remNegPrim_chopTertQuat@exprs    [,c(1,2)]))
  indices <- which(f_remNegPrim_chopTertQuat@exprs[,2] <= (upSlantmin + 2*((max(Rot_xy_leftPrim[2], Rot_xy_rightPrim[2]) + CutAbovePrimary*ScaleChop)-upSlantmin)) )
  f_remNegPrim_chopTertQuat@exprs <- f_remNegPrim_chopTertQuat@exprs[indices,]
  f_remNegPrim_chopTertQuat@exprs[,c(1,2)] <- t(t(R) %*% t(f_remNegPrim_chopTertQuat@exprs[,c(1,2)]))
  
  deviation2X <- deviation2Y <- vector()
  Xc2g <- Yc2g <- NULL
  
  if ( nrow(f_remNegPrim) >= NumberOfSinglePos^2*6 ) { # If less than 96 (4-plex) points for doub, trip and quad pops together, then use vector addition
    
    adjustDens <- 0.3
    repeat {
      
      # find thresholds to divide the secondary clusters
      highP <- 0.4
      lowP <- 0
      repeat{ # tinyP from 0.4 down to 0.04
        tinyP <- (highP+lowP)/2
        if ((tinyP + epsilon) > highP) break
        theta <- 0
        repeat{ # theta from 0 to pi/2
          f_rotate_tinyP_temp <- f_remNegPrim_chopTertQuat
          R <- matrix( c(cos(theta), sin(theta), -sin(theta), cos(theta)) ,2 ,2)
          f_rotate_tinyP_temp@exprs[,c(1,2)] <- t(R %*% t(f_rotate_tinyP_temp@exprs[,c(1,2)]))
          Cuts <- deGate(f_rotate_tinyP_temp, 1, all.cut=T, tinypeak.removal=tinyP, adjust.dens=adjustDens)
          if (length(Cuts) == (choose(NumberOfSinglePos, 2) - 1) ) break;
          if (theta >= pi/2) break;
          theta <- theta + pi/16
        }
        if (length(Cuts) < (choose(NumberOfSinglePos, 2) - 1)){
          highP <- tinyP
        } else if (length(Cuts) > (choose(NumberOfSinglePos, 2) - 1)) {
          lowP <- tinyP
        } else {
          break
        }
      }
      if ( length(Cuts) == (choose(NumberOfSinglePos, 2) - 1) ) break;
      if ( adjustDens > 0.5) {
        message("double positive detection failed.")
        break;
      }
      adjustDens <- adjustDens + 0.05
    }
    
    # find coordiates of secondary populations
    Cuts <- c(min(f_rotate_tinyP_temp@exprs[,1]), Cuts, max(f_rotate_tinyP_temp@exprs[,1]))
    
    for ( r1 in 1:(length(Cuts)-1) ) {
      indices <- intersect(which(f_rotate_tinyP_temp@exprs[,1] >= Cuts[r1]), which(f_rotate_tinyP_temp@exprs[,1] < Cuts[r1+1]))
      switch(r1,
             {
               XxA <- (firstClusters$clusters[1,1]+firstClusters$clusters[2,1] - emptyDroplets[1])
               YyA <- (firstClusters$clusters[1,2]+firstClusters$clusters[2,2] - emptyDroplets[2])
             },
             {
               XxA <- (firstClusters$clusters[1,1]+firstClusters$clusters[3,1] - emptyDroplets[1])
               YyA <- (firstClusters$clusters[1,2]+firstClusters$clusters[3,2] - emptyDroplets[2])
             },
             {
               XxA <- (firstClusters$clusters[2,1]+firstClusters$clusters[3,1] - emptyDroplets[1])
               YyA <- (firstClusters$clusters[2,2]+firstClusters$clusters[3,2] - emptyDroplets[2])
             },
             {
               XxA <- (firstClusters$clusters[1,1]+firstClusters$clusters[4,1] - emptyDroplets[1])
               YyA <- (firstClusters$clusters[1,2]+firstClusters$clusters[4,2] - emptyDroplets[2])
             },
             {
               XxA <- (firstClusters$clusters[2,1]+firstClusters$clusters[4,1] - emptyDroplets[1])
               YyA <- (firstClusters$clusters[2,2]+firstClusters$clusters[4,2] - emptyDroplets[2])
             },
             {
               XxA <- (firstClusters$clusters[3,1]+firstClusters$clusters[4,1] - emptyDroplets[1])
               YyA <- (firstClusters$clusters[3,2]+firstClusters$clusters[4,2] - emptyDroplets[2])
             }
      )
      if (length(indices) < 2) {
        deviation2X <- c(deviation2X, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= XxA+1000), which(f@exprs[,1] >= XxA-1000)), intersect(which(f@exprs[,2] <= YyA+1000), which(f@exprs[,2] >= YyA-1000))),1]))
        deviation2Y <- c(deviation2Y, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= XxA+1000), which(f@exprs[,1] >= XxA-1000)), intersect(which(f@exprs[,2] <= YyA+1000), which(f@exprs[,2] >= YyA-1000))),2]))
        Xc2g <- c(Xc2g, XxA)
        Yc2g <- c(Yc2g, YyA)
        next
      }
      f_Cuts_temp <- f_rotate_tinyP_temp; f_Cuts_temp@exprs <- f_Cuts_temp@exprs[indices,]
      f_Cuts_temp@exprs[,c(1,2)] <- t(t(R) %*% t(f_Cuts_temp@exprs[,c(1,2)]))
      
      highPXx <- 1
      lowPXx <- 0
      repeat{
        tinyPXx <- (highPXx + lowPXx)/2
        if ((tinyPXx + epsilon) > highPXx) {
          Xx <- Xx[1]
          message(paste0("Error finding X value for peak ", r1, " of ", choose(NumberOfSinglePos, 2)))
          break
        }
        Xx <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,1], width=1000), tinypeak.removal=tinyPXx)$Peaks
        if ( length (Xx) > 1 ) {
          lowPXx <- tinyPXx
        } else if (length (Xx) < 1) {
          highPXx <- tinyPXx
        } else {
          break
        }
      }
      highPYy <- 1
      lowPYy <- 0
      repeat{
        tinyPYy <- (highPYy + lowPYy)/2
        if ((tinyPYy + epsilon) > highPXx) {
          Yy <- Yy[1]
          message(paste0("Error finding Y value for peak ", r1, " of ", choose(NumberOfSinglePos, 2)))
          break
        }
        Yy <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,2], width=1000), tinypeak.removal=tinyPYy)$Peaks
        if ( length (Yy) > 1 ) {
          lowPYy <- tinyPYy
        } else if (length (Yy) < 1) {
          highPYy <- tinyPYy
        } else {
          break
        }
      }
      if (abs(Xx - XxA) > 3*scalingParam[1] || abs(Yy - YyA) > 3*scalingParam[2]) {
        deviation2X <- c(deviation2X, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= XxA+1000), which(f@exprs[,1] >= XxA-1000)), intersect(which(f@exprs[,2] <= YyA+1000), which(f@exprs[,2] >= YyA-1000))),1]))
        deviation2Y <- c(deviation2Y, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= XxA+1000), which(f@exprs[,1] >= XxA-1000)), intersect(which(f@exprs[,2] <= YyA+1000), which(f@exprs[,2] >= YyA-1000))),2]))
        Xc2g <- c(Xc2g, XxA)
        Yc2g <- c(Yc2g, YyA)
      } else {
        deviation2X <- c(deviation2X, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),1]))
        deviation2Y <- c(deviation2Y, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),2]))
        Xc2g <- c(Xc2g, Xx)
        Yc2g <- c(Yc2g, Yy)
      }
    }
  } else {
    # vector additon
    for (r1 in 1:NumberOfSinglePos) {
      for (r2 in 1:NumberOfSinglePos) {
        if (r1 >= r2) next
        Xx <- firstClusters$clusters[r1,1] + firstClusters$clusters[r2,1] - emptyDroplets[1]
        Yy <- firstClusters$clusters[r1,2] + firstClusters$clusters[r2,2] - emptyDroplets[2]
        deviation2X <- c(deviation2X, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),1]))
        deviation2Y <- c(deviation2Y, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),2]))
        Xc2g <- c(Xc2g, Xx)
        Yc2g <- c(Yc2g, Yy)
      }
    }
  }
  # switch ordering of secondary populations if if staement is true
  if( length(Xc2g) > 3 && (Xc2g[3]^2 + Yc2g[3]^2) >= 1.2*(Xc2g[4]^2 + Yc2g[4]^2) ) {
    message("Switching the middle two double populations.")
    Xc2g[c(3,4)] <- Xc2g[c(4,3)]
    Yc2g[c(3,4)] <- Yc2g[c(4,3)]
  }
  if( length(Xc2g) > 1 && Xc2g[1] > Xc2g[2] ) {
    message("Switching the top left two double populations.")
    Xc2g[c(1,2)] <- Xc2g[c(2,1)]
    Yc2g[c(1,2)] <- Yc2g[c(2,1)]
  }
  if( length(Yc2g) > 5 && Yc2g[5] < Yc2g[6] ) {
    message("Switching the bottom right two double populations.")
    Xc2g[c(5,6)] <- Xc2g[c(6,5)]
    Yc2g[c(5,6)] <- Yc2g[c(6,5)]
  }
  
  return(list(clusters=cbind(Xc2g, Yc2g), deviation = cbind(deviation2X, deviation2Y)))
}

# Find the tertiary clusters based on their density peaks found by flowDensity.
findTertiaryClustersDensity <- function(f, f_remNegPrimSec, emptyDroplets, firstClusters, secondaryClusters) {
  
  NumberOfSinglePos <- nrow(firstClusters$clusters)
  XcTert <- YcTert <- NULL
  deviation2X <- deviation2Y <- vector()
  
  if ( nrow(f_remNegPrimSec) >= NumberOfSinglePos^2*3 ) { # If less than 50 points for trip and quad pops together, then use vector addition
    
    # find thresholds to divide the tertiary clusters
    highP <- 1
    lowP <- 0
    repeat { # newton iteration
      newP <- (highP+lowP)/2
      Cuts <- deGate(f_remNegPrimSec, 1, all.cut=T, tinypeak.removal=newP, adjust=0.1)
      if (length(Cuts) < (NumberOfSinglePos - 1)){
        highP <- newP
      } else if (length(Cuts) > (NumberOfSinglePos - 1)) {
        lowP <- newP
      } else {
        break
      }
    }
    
    Cuts <- c(min(f_remNegPrimSec@exprs[,1]), Cuts, max(f_remNegPrimSec@exprs[,1]))
    
    # find coordiates of tertiary populations
    for ( r1 in 1:(length(Cuts)-1) ) {
      
      indices <- intersect(which(f_remNegPrimSec@exprs[,1] >= Cuts[r1]), which(f_remNegPrimSec@exprs[,1] < Cuts[r1+1]))
      
      switch(r1,
             {
               XxA <- mean(secondaryClusters$clusters[1,1]+firstClusters$clusters[3,1] - emptyDroplets[1],
                           secondaryClusters$clusters[2,1]+firstClusters$clusters[2,1] - emptyDroplets[1],
                           secondaryClusters$clusters[3,1]+firstClusters$clusters[1,1] - emptyDroplets[1])
               YyA <- mean(secondaryClusters$clusters[1,2]+firstClusters$clusters[3,2] - emptyDroplets[2],
                           secondaryClusters$clusters[2,2]+firstClusters$clusters[2,2] - emptyDroplets[2],
                           secondaryClusters$clusters[3,2]+firstClusters$clusters[1,2] - emptyDroplets[2])
             },
             {
               XxA <- mean(secondaryClusters$clusters[1,1]+firstClusters$clusters[4,1] - emptyDroplets[1],
                           secondaryClusters$clusters[4,1]+firstClusters$clusters[2,1] - emptyDroplets[1],
                           secondaryClusters$clusters[5,1]+firstClusters$clusters[1,1] - emptyDroplets[1])
               YyA <- mean(secondaryClusters$clusters[1,2]+firstClusters$clusters[4,2] - emptyDroplets[2],
                           secondaryClusters$clusters[4,2]+firstClusters$clusters[2,2] - emptyDroplets[2],
                           secondaryClusters$clusters[5,2]+firstClusters$clusters[1,2] - emptyDroplets[2])
             },
             {
               XxA <- mean(secondaryClusters$clusters[2,1]+firstClusters$clusters[4,1] - emptyDroplets[1],
                           secondaryClusters$clusters[4,1]+firstClusters$clusters[3,1] - emptyDroplets[1],
                           secondaryClusters$clusters[6,1]+firstClusters$clusters[1,1] - emptyDroplets[1])
               YyA <- mean(secondaryClusters$clusters[2,2]+firstClusters$clusters[4,2] - emptyDroplets[2],
                           secondaryClusters$clusters[4,2]+firstClusters$clusters[3,2] - emptyDroplets[2],
                           secondaryClusters$clusters[6,2]+firstClusters$clusters[1,2] - emptyDroplets[2])
             },
             {
               XxA <- mean(secondaryClusters$clusters[3,1]+firstClusters$clusters[4,1] - emptyDroplets[1],
                           secondaryClusters$clusters[5,1]+firstClusters$clusters[3,1] - emptyDroplets[1],
                           secondaryClusters$clusters[6,1]+firstClusters$clusters[2,1] - emptyDroplets[1])
               YyA <- mean(secondaryClusters$clusters[3,2]+firstClusters$clusters[4,2] - emptyDroplets[2],
                           secondaryClusters$clusters[5,2]+firstClusters$clusters[3,2] - emptyDroplets[2],
                           secondaryClusters$clusters[6,2]+firstClusters$clusters[2,2] - emptyDroplets[2])
             }
      )
      if (length(indices) < 2) {
        deviation2X <- c(deviation2X, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= XxA+1000), which(f@exprs[,1] >= XxA-1000)), intersect(which(f@exprs[,2] <= YyA+1000), which(f@exprs[,2] >= YyA-1000))),1]))
        deviation2Y <- c(deviation2Y, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= XxA+1000), which(f@exprs[,1] >= XxA-1000)), intersect(which(f@exprs[,2] <= YyA+1000), which(f@exprs[,2] >= YyA-1000))),2]))
        XcTert <- c(XcTert, XxA)
        YcTert <- c(YcTert, YyA)
        next
      }
      f_remNegPrimSec_temp <- f_remNegPrimSec; f_remNegPrimSec_temp@exprs <- f_remNegPrimSec_temp@exprs[indices,]
      f_remNegPrimSec_temp@exprs[,c(1,2)] <- t(t(R) %*% t(f_remNegPrimSec_temp@exprs[,c(1,2)]))
      
      highPXx <- 1
      lowPXx <- 0
      repeat{
        tinyPXx <- (highPXx + lowPXx)/2
        if ((tinyPXx + epsilon) > highPXx) {
          Xx <- Xx[1]
          message(paste0("Default X value for tertiary peak ", r1, " of 4"))
          break
        }
        Xx <- flowDensity:::.getPeaks(density(f_remNegPrimSec_temp@exprs[,1], width=1000), tinypeak.removal=tinyPXx)$Peaks
        if ( length (Xx) > 1 ) {
          lowPXx <- tinyPXx
        } else if (length (Xx) < 1) {
          highPXx <- tinyPXx
        } else {
          break
        }
      }
      highPYy <- 1
      lowPYy <- 0
      repeat{
        tinyPYy <- (highPYy + lowPYy)/2
        if ((tinyPYy + epsilon) > highPYy) {
          Yy <- Yy[1]
          message(paste0("Default Y value for tertiary peak ", r1, " of 4"))
          break
        }
        Yy <- flowDensity:::.getPeaks(density(f_remNegPrimSec_temp@exprs[,2], width=1000), tinypeak.removal=tinyPYy)$Peaks
        if ( length (Yy) > 1 ) {
          lowPYy <- tinyPYy
        } else if (length (Yy) < 1) {
          highPYy <- tinyPYy
        } else {
          break
        }
      }
      if (abs(Xx - XxA) > 4*scalingParam[1] || abs(Yy - YyA) > 4*scalingParam[2]) {
        deviation2X <- c(deviation2X, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= XxA+1000), which(f@exprs[,1] >= XxA-1000)), intersect(which(f@exprs[,2] <= YyA+1000), which(f@exprs[,2] >= YyA-1000))),1]))
        deviation2Y <- c(deviation2Y, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= XxA+1000), which(f@exprs[,1] >= XxA-1000)), intersect(which(f@exprs[,2] <= YyA+1000), which(f@exprs[,2] >= YyA-1000))),2]))
        XcTert <- c(XcTert, XxA)
        YcTert <- c(YcTert, YyA)
      } else {
        deviation2X <- c(deviation2X, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),1]))
        deviation2Y <- c(deviation2Y, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),2]))
        XcTert <- c(XcTert, Xx)
        YcTert <- c(YcTert, Yy)
      }
    }
    
  } else {
    for (r1 in 1:4) {
      switch(r1,
             {
               Xx <- mean(secondaryClusters$clusters[1,1]+firstClusters$clusters[3,1] - emptyDroplets[1],
                          secondaryClusters$clusters[2,1]+firstClusters$clusters[2,1] - emptyDroplets[1],
                          secondaryClusters$clusters[3,1]+firstClusters$clusters[1,1] - emptyDroplets[1])
               Yy <- mean(secondaryClusters$clusters[1,2]+firstClusters$clusters[3,2] - emptyDroplets[2],
                          secondaryClusters$clusters[2,2]+firstClusters$clusters[2,2] - emptyDroplets[2],
                          secondaryClusters$clusters[3,2]+firstClusters$clusters[1,2] - emptyDroplets[2])
             },
             {
               Xx <- mean(secondaryClusters$clusters[1,1]+firstClusters$clusters[4,1] - emptyDroplets[1],
                          secondaryClusters$clusters[4,1]+firstClusters$clusters[2,1] - emptyDroplets[1],
                          secondaryClusters$clusters[5,1]+firstClusters$clusters[1,1] - emptyDroplets[1])
               Yy <- mean(secondaryClusters$clusters[1,2]+firstClusters$clusters[4,2] - emptyDroplets[2],
                          secondaryClusters$clusters[4,2]+firstClusters$clusters[2,2] - emptyDroplets[2],
                          secondaryClusters$clusters[5,2]+firstClusters$clusters[1,2] - emptyDroplets[2])
             },
             {
               Xx <- mean(secondaryClusters$clusters[2,1]+firstClusters$clusters[4,1] - emptyDroplets[1],
                          secondaryClusters$clusters[4,1]+firstClusters$clusters[3,1] - emptyDroplets[1],
                          secondaryClusters$clusters[6,1]+firstClusters$clusters[1,1] - emptyDroplets[1])
               Yy <- mean(secondaryClusters$clusters[2,2]+firstClusters$clusters[4,2] - emptyDroplets[2],
                          secondaryClusters$clusters[4,2]+firstClusters$clusters[3,2] - emptyDroplets[2],
                          secondaryClusters$clusters[6,2]+firstClusters$clusters[1,2] - emptyDroplets[2])
             },
             {
               Xx <- mean(secondaryClusters$clusters[3,1]+firstClusters$clusters[4,1] - emptyDroplets[1],
                          secondaryClusters$clusters[5,1]+firstClusters$clusters[3,1] - emptyDroplets[1],
                          secondaryClusters$clusters[6,1]+firstClusters$clusters[2,1] - emptyDroplets[1])
               Yy <- mean(secondaryClusters$clusters[3,2]+firstClusters$clusters[4,2] - emptyDroplets[2],
                          secondaryClusters$clusters[5,2]+firstClusters$clusters[3,2] - emptyDroplets[2],
                          secondaryClusters$clusters[6,2]+firstClusters$clusters[2,2] - emptyDroplets[2])
             }
      )
      deviation2X <- c(deviation2X, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),1]))
      deviation2Y <- c(deviation2Y, sd(f@exprs[intersect(intersect(which(f@exprs[,1] <= Xx+1000), which(f@exprs[,1] >= Xx-1000)), intersect(which(f@exprs[,2] <= Yy+1000), which(f@exprs[,2] >= Yy-1000))),2]))
      XcTert <- c(XcTert, Xx)
      YcTert <- c(YcTert, Yy)
    }
  }
  return(list(clusters=cbind(XcTert, YcTert), deviation = cbind(deviation2X, deviation2Y)))
}

# Find the quaternary clusters based on their density peaks found by flowDensity.
findQuaternaryClusterDensity <- function(f, f_onlyQuad, emptyDroplets, firstClusters, secondaryClusters, tertiaryClusters) {
  
  NumberOfSinglePos <- nrow(firstClusters$clusters)
  if ( nrow(f_onlyQuad) >= 4 ) {
    
    highP <- 1
    lowP <- 0
    repeat{
      tinyP <- (highP+lowP)/2
      XcQuad <- flowDensity:::.getPeaks(density(f_onlyQuad@exprs[,1], width=1000), tinypeak.removal=tinyP)$Peaks
      if (length(XcQuad) < 1 ) {
        highP <- tinyP
      } else if (length(XcQuad) > 1){
        lowP <- tinyP
      } else {
        break
      }
    }
    highP <- 1
    lowP <- 0
    repeat{
      tinyP <- (highP+lowP)/2
      YcQuad <- flowDensity:::.getPeaks(density(f_onlyQuad@exprs[,2], width=1000), tinypeak.removal=tinyP)$Peaks
      if (length(YcQuad) < 1 ) {
        highP <- tinyP
      } else if (length(YcQuad) > 1){
        lowP <- tinyP
      } else {
        break
      }
    }
  } else {
    if (is.null(tertiaryClusters)) {
      XcQuad <- mean(secondaryClusters$clusters[1,1]+firstClusters$clusters[3,1] - emptyDroplets[1],
                     secondaryClusters$clusters[2,1]+firstClusters$clusters[2,1] - emptyDroplets[1],
                     secondaryClusters$clusters[3,1]+firstClusters$clusters[1,1] - emptyDroplets[1])
      YcQuad <- mean(secondaryClusters$clusters[1,2]+firstClusters$clusters[3,2] - emptyDroplets[2],
                     secondaryClusters$clusters[2,2]+firstClusters$clusters[2,2] - emptyDroplets[2],
                     secondaryClusters$clusters[3,2]+firstClusters$clusters[1,2] - emptyDroplets[2])
    } else if (is.null(secondaryClusters)) {
      XcQuad <- firstClusters$clusters[1,1] + firstClusters[2,1] - emptyDroplets[1]
      YcQuad <- firstClusters$clusters[1,2] + firstClusters[2,2] - emptyDroplets[2] 
    } else {
      XcQuad <- mean(tertiaryClusters$clusters[1,1]+firstClusters$clusters[4,1] - emptyDroplets[1],
                     tertiaryClusters$clusters[2,1]+firstClusters$clusters[3,1] - emptyDroplets[1],
                     tertiaryClusters$clusters[3,1]+firstClusters$clusters[2,1] - emptyDroplets[1],
                     tertiaryClusters$clusters[4,1]+firstClusters$clusters[1,1] - emptyDroplets[1])
      YcQuad <- mean(tertiaryClusters$clusters[1,2]+firstClusters$clusters[4,2] - emptyDroplets[2],
                     tertiaryClusters$clusters[2,2]+firstClusters$clusters[3,2] - emptyDroplets[2],
                     tertiaryClusters$clusters[3,2]+firstClusters$clusters[2,2] - emptyDroplets[2],
                     tertiaryClusters$clusters[4,2]+firstClusters$clusters[1,2] - emptyDroplets[2])
    }
  }
  return(list(clusters=cbind(XcQuad, YcQuad)))
} 

# Find the primary clusters based on their density peaks found by flowDensity.
findPrimaryClusters <- function(data, clusterMeans, emptyDroplets, remove=0, dimensions, File, f, NumberOfSinglePos=4) {
  
  NumOfClusters <- NumberOfSinglePos^2
  DataRemoved <- FinalResults <- NULL
  
  up1max <- deGate(f, c(1), percentile=0.999, use.percentile=T)
  up1min <- deGate(f, c(1), percentile=0.001, use.percentile=T)
  up2max <- deGate(f, c(2), percentile=0.999, use.percentile=T)
  up2min <- deGate(f, c(2), percentile=0.001, use.percentile=T)
  
  indices <- unique( c( which(f@exprs[,1] >= 0.15 * (up1max - up1min) + up1min),
                        which(f@exprs[,2] >= 0.15 * (up2max - up2min) + up2min) ) )
  
  f_remNeg  <- f; f_remNeg @exprs <- f_remNeg @exprs[ indices,] # keep the non 15% bottom left corner (rn for removed neg)
  f_onlyNeg <- f; f_onlyNeg@exprs <- f_onlyNeg@exprs[-indices,] # keep the     15% bottom left corner
  
  ClusterCentres <- vector()
  
  #---- find 1st gen clusters------------------------------------------------------------------------------------------------------------------------#
  
  # Calculate the threshold to remove events that are too positive. This makes it easier to find single positives.
  f_findExtremes_temp <- f
  f_remNeg_temp <- f_remNeg
  threshold <- 0.1
  
  # find left and right primary
  dataTable <- table(data)
  clusterMeans2 <- clusterMeans[-c(emptyDroplets, remove), , drop=F]
  if (length(clusterMeans2) == 0) return(NULL)
  minimum <- match(min(clusterMeans2[,1]), clusterMeans[,1])
  a <- matrix(1, nrow=3)
  realFirstCluster1 <- vector()
  collinear1 <- minimum
  for (i in 1:nrow(clusterMeans)) {
    if (i == minimum || i == emptyDroplets)
      next
    test <- matrix(clusterMeans[c(emptyDroplets,minimum,i),], ncol=2)
    if (abs(det(cbind(a, test/100))) < (dimensions[1]/(NumberOfSinglePos*20)) && clusterMeans[i,2] < (1-NumberOfSinglePos^2/100)*dimensions[2])
      collinear1 <- c(collinear1, i)
  }
  selection <- which(clusterMeans[,1] <= max(clusterMeans[collinear1,1]))
  selection <- selection[!selection %in% c(emptyDroplets, remove)]
  selection <- selection[(dataTable[selection] > max(dataTable[selection])/4)]
  realFirstCluster1 <- selection[which.max(clusterMeans[selection, 2])]
  
  # collinear1 <- realFirstCluster1
  # for (i in 1:nrow(clusterMeans)) {
  #   if (i == realFirstCluster1 || i == emptyDroplets)
  #     next
  #   test <- matrix(clusterMeans[c(emptyDroplets,realFirstCluster1,i),], ncol=2)
  #   if (abs(det(cbind(a, test/100))) < (dimensions[1]/(NumberOfSinglePos*20)))
  #     collinear1 <- c(collinear1, i)
  # }
  
  clusterMeans2 <- clusterMeans[-c(emptyDroplets, collinear1, remove), , drop=F]
  if (length(clusterMeans2) == 0) return(realFirstCluster1)
  minimum <- match(min(clusterMeans2[,2]), clusterMeans[,2])
  a <- matrix(1, nrow=3)
  realFirstCluster2 <- vector()
  collinear2 <- minimum
  for (i in 1:nrow(clusterMeans)) {
    if (i == minimum || i == emptyDroplets)
      next
    test <- matrix(clusterMeans[c(emptyDroplets,minimum,i),], ncol=2)
    #   print(abs(det(cbind(a, test/100))))
    if (abs(det(cbind(a, test/100))) < (dimensions[2]/(NumberOfSinglePos*20)) && clusterMeans[i,1] < (1-NumberOfSinglePos^2/100)*dimensions[1])
      collinear2 <- c(collinear2, i)
  }
  selection <- which(clusterMeans[,2] <= max(clusterMeans[collinear2,2]))
  selection <- selection[!selection %in% c(emptyDroplets, collinear1, remove)]
  selection <- selection[(dataTable[selection] > max(dataTable[selection])/4)]
  realFirstCluster2 <- selection[which.max(clusterMeans[selection, 1])]
  
  if (dataTable[realFirstCluster2] < dataTable[realFirstCluster1]/10) {
    cat(paste("Something wrong with primary cluster detection (Ch2: only", dataTable[realFirstCluster2], "droplets). Trying again..."))
    remove <- c(remove, collinear2)
    return(findPrimaryClusters(data, clusterMeans, emptyDroplets, remove, 1.1*dimensions, File, f, NumberOfSinglePos))
  }

  if (dataTable[realFirstCluster1] < dataTable[realFirstCluster2]/10) {
    cat(paste("Something wrong with primary cluster detection (Ch1: only", dataTable[realFirstCluster1], "droplets). Trying again..."))
    remove <- c(remove, collinear1)
    return(findPrimaryClusters(data, clusterMeans, emptyDroplets, remove, 1.1*dimensions, File, f, NumberOfSinglePos))
  }
  
  x_leftPrim <- clusterMeans[realFirstCluster1,1]
  y_leftPrim <- clusterMeans[realFirstCluster1,2]
  x_rightPrim <- clusterMeans[realFirstCluster2,1]
  y_rightPrim <- clusterMeans[realFirstCluster2,2]
  
  mSlope <- (y_leftPrim - y_rightPrim)/(x_leftPrim - x_rightPrim)
  theta <- abs(atan(mSlope))
  R <- matrix( c(cos(theta), sin(theta), -sin(theta), cos(theta)) ,2 ,2)
  
  Rot_xy_leftPrim <- R %*% c(x_leftPrim, y_leftPrim) # coordinates of rotated left  primary cluster
  Rot_xy_rightPrim <- R %*% c(x_rightPrim, y_rightPrim) # coordinates of rotated right primary cluster
  
  f_findExtremes_temp@exprs[,c(1,2)] <- t(R %*% t(f_findExtremes_temp@exprs[,c(1,2)]))
  f_remNeg_temp@exprs    [,c(1,2)] <- t(R %*% t(f_remNeg_temp@exprs    [,c(1,2)]))
  upSlantmax <- deGate(f_findExtremes_temp, c(2), percentile=0.999, use.percentile=T)
  upSlantmin <- deGate(f_findExtremes_temp, c(2), percentile=0.001, use.percentile=T)
  # Chop off the very positive events
  ScaleChop <-  (upSlantmax - upSlantmin) / max(File)
  indices <- which(f_remNeg_temp@exprs[,2] <= (max(Rot_xy_leftPrim[2], Rot_xy_rightPrim[2]) + CutAbovePrimary*ScaleChop))
  f_remNeg_temp@exprs <- f_remNeg_temp@exprs[indices,]
  
  # find thresholds to divide the primary clusters
  highP <- 1
  lowP <- 0
  repeat { # newton iteration
    newP <- (highP+lowP)/2
    if (newP < 2*epsilon) break
    Cuts <- deGate(f_remNeg_temp, 1, all.cut=T, tinypeak.removal=newP, adjust=0.1)
    if (length(Cuts) < (NumberOfSinglePos - 1)){
      highP <- newP
    } else if (length(Cuts) > (NumberOfSinglePos - 1)) {
      lowP <- newP
    } else {
      break
    }
  }
  Xc <- Yc <- NULL
  Cuts <- c(min(f_remNeg_temp@exprs[,1]), Cuts, max(f_remNeg_temp@exprs[,1]))
  
  # find the coordinates of the primary clusters
  for ( r1 in 1:(length(Cuts)-1) ) {
    indices <- intersect(which(f_remNeg_temp@exprs[,1] >= Cuts[r1]), which(f_remNeg_temp@exprs[,1] < Cuts[r1+1]))
    f_Cuts_temp <- f_remNeg_temp; f_Cuts_temp@exprs <- f_Cuts_temp@exprs[indices,]
    f_Cuts_temp@exprs[,c(1,2)] <- t(t(R) %*% t(f_Cuts_temp@exprs[,c(1,2)]))
    Xx <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,1], width=1000), tinypeak.removal=newP)$Peaks
    Yy <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,2], width=1000), tinypeak.removal=newP)$Peaks
    if (length(Xx) > 1 || length(Yy) > 1) {
      highP <- 1
      lowP <- 0
      repeat { # newton iteration
        newP <- (highP+lowP)/2
        if (newP < epsilon) break
        Xx <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,1], width=1000), tinypeak.removal=newP)$Peaks
        Yy <- flowDensity:::.getPeaks(density(f_Cuts_temp@exprs[,2], width=1000), tinypeak.removal=newP)$Peaks
        if (length(Xx) > 1 || length(Yy) > 1){
          highP <- newP
        } else if (length(Xx) < 1 || length(Yy) < 1) {
          lowP <- newP
        } else {
          break
        }
      }
    }
    Xc <- c(Xc, Xx[1])
    Yc <- c(Yc, Yy[1])
  }
  
  densPrimaries <- cbind(Xc, Yc)
  firstClusters <- vector()
  distMatrix <- vector()
  for (i in 1:nrow(densPrimaries)) {
    temp <- lapply(1:nrow(clusterMeans), function(x) return(clusterMeans[x,] - densPrimaries[i,]))
    rowSums <- sapply(1:length(temp), function(x) return(sum(abs(unlist(temp[x])))))
    distMatrix <- rbind(distMatrix, rowSums)
  }
  firstClusters <- solve_LSAP(distMatrix)
}

# Find the primary clusters based on their position.
findPrimaryClusters_old <- function(data, clusterMeans, emptyDroplets, remove=0, dimensions) {
  
  dataTable <- table(data)
  clusterMeans2 <- clusterMeans[-c(emptyDroplets, remove), , drop=F]
  if (length(clusterMeans2) == 0) return(NULL)
  minimum <- match(min(clusterMeans2[,1]), clusterMeans[,1])
  a <- matrix(1, nrow=3)
  realFirstCluster1 <- vector()
  collinear1 <- minimum
  for (i in 1:nrow(clusterMeans)) {
    if (i == minimum || i == emptyDroplets)
      next
    test <- matrix(clusterMeans[c(emptyDroplets,minimum,i),], ncol=2)
    if (abs(det(cbind(a, test/100))) < (dimensions[1]/100))
      collinear1 <- c(collinear1, i)
  }
  selection <- which(clusterMeans[,1] <= max(clusterMeans[collinear1,1]))
  selection <- selection[!selection %in% c(emptyDroplets, remove)]
  selection <- selection[(dataTable[selection] > max(dataTable[selection])/4)]
  realFirstCluster1 <- selection[which.min(clusterMeans[selection, 1])]
  
  clusterMeans2 <- clusterMeans[-c(emptyDroplets, collinear1, remove), , drop=F]
  if (length(clusterMeans2) == 0) return(realFirstCluster1)
  minimum <- match(min(clusterMeans2[,2]), clusterMeans[,2])
  a <- matrix(1, nrow=3)
  realFirstCluster2 <- vector()
  collinear2 <- minimum
  for (i in 1:nrow(clusterMeans)) {
    if (i == minimum || i == emptyDroplets)
      next
    test <- matrix(clusterMeans[c(emptyDroplets,minimum,i),], ncol=2)
    #   print(abs(det(cbind(a, test/100))))
    if (abs(det(cbind(a, test/100))) < (dimensions[2]/100))
      collinear2 <- c(collinear2, i)
  }
  selection <- which(clusterMeans[,2] <= max(clusterMeans[collinear2,2]))
  selection <- selection[!selection %in% c(emptyDroplets, collinear1, remove)]
  selection <- selection[(dataTable[selection] > max(dataTable[selection])/4)]
  realFirstCluster2 <- selection[which.min(clusterMeans[selection, 2])]
  
  if (dataTable[realFirstCluster2] < dataTable[realFirstCluster1]/4) {
    warning("Something wrong with initial cluster detection, to few events for Ch1 amplitude, try again...")
    remove <- c(remove, collinear2)
    return(findPrimaryClusters(data, clusterMeans, emptyDroplets, remove, dimensions))
  }
  
  if (dataTable[realFirstCluster1] < dataTable[realFirstCluster2]/4) {
    warning("Something wrong with initial cluster detection, to few events for Ch1 amplitude, try again....")
    remove <- c(remove, collinear1)
    return(findPrimaryClusters(data, clusterMeans, emptyDroplets, remove, dimensions))
  }
  
  clusterMeans2 <- clusterMeans[-c(emptyDroplets, collinear1, collinear2), , drop=F]
  moreFirstClusters <- vector()
  counter <- clusterMeans[realFirstCluster2, 1]
  clusterMeans2 <- clusterMeans2[clusterMeans2[,1]  < counter, ,drop = FALSE]
  counter <- clusterMeans[realFirstCluster1, 2]
  temp <- clusterMeans2[clusterMeans2[,2]  < counter, ,drop = FALSE]
  repeat {
    if (length(temp) == 0) {
      break
    } else {
      minimum <- which.min(temp[,1])
      cluster <- match(temp[minimum], clusterMeans)
      if (dataTable[cluster] < (mean(dataTable[c(realFirstCluster1, realFirstCluster2)]))/5) {
        temp <- temp[-minimum, , drop=F]
        next
      } else {
        moreFirstClusters <- c(moreFirstClusters, cluster)
        counter <- clusterMeans[cluster, 2]
        temp <- clusterMeans2[clusterMeans2[,2]  < counter, ,drop = FALSE]
      }
    }
  }
  
  #   newTable <- dataTable[-c(emptyDroplets, realFirstCluster1, realFirstCluster2)]
  #   moreFirstClusters <- names(newTable[newTable > (mean(dataTable[c(realFirstCluster1, realFirstCluster2)]))/2])
  #
  realFirstClusters <- as.integer(c(realFirstCluster1, moreFirstClusters, realFirstCluster2))
}

# Find the secondary clusters based on the positions of the primary clusters.
findSecondaryClusters <- function(firstClusters, clusterMeans, emptyDroplets, remove, dimensions, counts) {
  threshold <- dimensions/5
  clusterMeans2 <- clusterMeans[-c(emptyDroplets, firstClusters, remove),]
  if (length(clusterMeans2) == 0 ) return(list(clusters=0, correctionFactor=0))
  secondClusters <- vector()
  correction <- list()
  correction[[length(firstClusters)+1]] <- 0
  
  if(nrow(clusterMeans2) >= sum(1:(length(firstClusters)-1))) {
    distMatrix <- vector()
    for (a in firstClusters) {
      for (b in firstClusters) {
        if (match(a, firstClusters) >= match(b, firstClusters)) next
        temp <- lapply(1:nrow(clusterMeans2), function(x) return(clusterMeans2[x,] - (clusterMeans[a,] + clusterMeans[b,] - clusterMeans[emptyDroplets,])))
        #       temp2 <- lapply(1:length(temp), function(x) return(c(temp[[x]][1]*3, temp[[x]][2]/3)))
        rowSums <- sapply(1:length(temp), function(x) return(sum(abs(unlist(temp[x])))))
        distMatrix <- rbind(distMatrix, rowSums)
      }
    }
    tmp <- solve_LSAP(distMatrix)
    secondClusters <- match(clusterMeans2[tmp], clusterMeans)
    index <- 0
    for (a in firstClusters) {
      for (b in firstClusters) {
        if (match(a, firstClusters) >= match(b, firstClusters)) next
        index <- index + 1
        distance <- clusterMeans[secondClusters[index],] - (clusterMeans[a,] + clusterMeans[b,] - clusterMeans[emptyDroplets,])
        if (sum(abs(distance)) > threshold) {
          secondClusters[index] = 0
        } else {
          correction[[match(a, firstClusters)]] <- c(correction[[match(a, firstClusters)]], list(distance))
          correction[[match(b, firstClusters)]] <- c(correction[[match(b, firstClusters)]], list(distance))
        }
      }
    }
    
  } else {
    
    for (a in firstClusters) {
      corTemp <- 0
      for (b in firstClusters) {
        if (match(a, firstClusters) >= match(b, firstClusters)) {
          next
        } else {
          if (length(clusterMeans2) == 0 ) {
            cluster <- 0
            secondClusters <- c(secondClusters, cluster)
            next
          }
          temp <- lapply(1:nrow(clusterMeans2), function(x) return(clusterMeans2[x,] - (clusterMeans[a,] + clusterMeans[b,] - clusterMeans[emptyDroplets,])))
          rowSums <- sapply(1:length(temp), function(x) return(sum(abs(unlist(temp[x])))))
          if (min(rowSums) < threshold) {
            minimum <- which.min(rowSums)
            cluster <- match(clusterMeans2[minimum], clusterMeans)
            clusterMeans2 <- clusterMeans2[-minimum, , drop=F]
            correction[[match(a, firstClusters)]] <- c(correction[[match(a, firstClusters)]], temp[minimum])
            correction[[match(b, firstClusters)]] <- c(correction[[match(b, firstClusters)]], temp[minimum])
            corTemp <- temp[[minimum]]
          } else {
            temp <- lapply(1:nrow(clusterMeans2), function(x) return(clusterMeans2[x,] - (clusterMeans[a,] + clusterMeans[b,] - clusterMeans[emptyDroplets,] + corTemp)))
            rowSums <- sapply(1:length(temp), function(x) return(sum(abs(unlist(temp[x])))))
            if (min(rowSums) < threshold) {
              minimum <- which.min(rowSums)
              cluster <- match(clusterMeans2[minimum], clusterMeans)
              clusterMeans2 <- clusterMeans2[-minimum, , drop=F]
              correction[[match(a, firstClusters)]] <- c(correction[[match(a, firstClusters)]], temp[minimum])
              correction[[match(b, firstClusters)]] <- c(correction[[match(b, firstClusters)]], temp[minimum])
              corTemp <- temp[[minimum]]
            } else {
              cluster <- 0
              #correctionFactor <- c(correctionFactor, 0)
            }
          }
          
          secondClusters <- c(secondClusters, cluster)
        }
      }
    }
  }
  correctionFactor <- list()
  for (i in 1:(length(correction)-1)) {
    if (length(correction[[i]]) == 0) {
      correctionFactor[[i]] <- 0
    } else {
      temp <- t(do.call(cbind, correction[[i]]))
      correctionFactor[[i]] <- colMeans(temp)
    }
  }
  
  if (max(secondClusters) == 0) return(list(clusters=secondClusters, correctionFactor=correctionFactor))
  
  q <- quantile(counts[secondClusters])
  
  if ((IQR(counts[secondClusters]) > q[3]) && (q[3] > 5) && (q[1] < 5)) {
    remove <- c(remove, as.integer(names(which.min(counts[secondClusters]))))
    print(counts[secondClusters])
    return(findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, remove, dimensions, counts))
  }
  if (min(counts[secondClusters]) < q[2]/5) {
    remove <- c(remove, as.integer(names(which.min(counts[secondClusters]))))
    print(counts[secondClusters])
    return(findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, remove, dimensions, counts))
  }
  if (max(counts[secondClusters]) > q[4]*5) {
    remove <- c(remove, as.integer(names(which.max(counts[secondClusters]))))
    print(counts[secondClusters])
    return(findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, remove, dimensions, counts))
  }
  
  ## the correction factor is the avarage distance between the estimated positions of secondClusters and the real positions
  list(clusters=secondClusters, correctionFactor=correctionFactor)
}

# Find the tertiary clusters based on the positions of the primary and secondary clusters.
findTertiaryClusters <- function(emptyDroplets, firstClusters, secondClusters, remove, clusterMeans, correctionFactor, dimensions, counts) {
  threshold <- dimensions/6
  clusterMeans2 <- clusterMeans[-c(emptyDroplets, firstClusters, secondClusters, remove), ,drop=F]
  if (length(clusterMeans2) == 0 ) return(list(clusters=rep(0,4)))
  thirdClusters1 <- vector()
  thirdClusters2 <- vector()
  
  if (nrow(clusterMeans2) >= 4) {
    
    for (i in 1:4) {
      if (i<3) {
        cluster <- tail(secondClusters, n=1)
        if (cluster == 0 || length(clusterMeans2) == 0) {
          thirdClusters1 <- rbind(thirdClusters1, rep(dimensions, nrow(clusterMeans2)))
          next
        }
        temp <- lapply(1:nrow(clusterMeans2), function(x) return(clusterMeans2[x,] - (clusterMeans[cluster,] + clusterMeans[firstClusters[i],] - clusterMeans[emptyDroplets,] + correctionFactor[[i]])))
        rowSums <- sapply(1:length(temp), function(x) return(sum(abs(unlist(temp[x])))))
        if (min(rowSums) > threshold) {
          thirdClusters1 <- rbind(thirdClusters1, rep(dimensions, nrow(clusterMeans2)))
        } else {
          thirdClusters1 <- rbind(thirdClusters1, rowSums)
        }
      } else {
        cluster <- secondClusters[1]
        if (cluster == 0 || length(clusterMeans2) == 0) {
          thirdClusters2 <- rbind(thirdClusters2, rep(dimensions, nrow(clusterMeans2)))
          next
        }
        temp <- lapply(1:nrow(clusterMeans2), function(x) return(clusterMeans2[x,] - (clusterMeans[cluster,] + clusterMeans[firstClusters[i],] - clusterMeans[emptyDroplets,] + correctionFactor[[i]])))
        rowSums <- sapply(1:length(temp), function(x) return(sum(abs(unlist(temp[x])))))
        if (min(rowSums) > threshold) {
          thirdClusters2 <- rbind(thirdClusters2, rep(dimensions, nrow(clusterMeans2)))
        } else {
          thirdClusters2 <- rbind(thirdClusters2, rowSums)
        }
      }
    }
    distMatrix <- rbind(thirdClusters2, thirdClusters1)
    tmp <- solve_LSAP(distMatrix)
    thirdClusters <- match(clusterMeans2[tmp], clusterMeans)
    for (i in 1:length(tmp)) {
      if (distMatrix[i, tmp[i]] > threshold) {
        thirdClusters[i] = 0
      }
    }
    
  } else {
    
    for (i in 1:4) {
      if (i<3) {
        cluster <- tail(secondClusters, n=1)
        if (cluster == 0 || length(clusterMeans2) == 0) {
          thirdClusters1 <- c(thirdClusters1, 0)
          next
        }
        temp <- lapply(1:nrow(clusterMeans2), function(x) return(clusterMeans2[x,] - (clusterMeans[cluster,] + clusterMeans[firstClusters[i],] - 2*clusterMeans[emptyDroplets,] + correctionFactor[[i]])))
        rowSums <- sapply(1:length(temp), function(x) return(sum(abs(unlist(temp[x])))))
        if (min(rowSums) < threshold) {
          minimum <- which.min(rowSums)
          cluster <- match(clusterMeans2[minimum], clusterMeans)
          clusterMeans2 <- clusterMeans2[-minimum, , drop=F]
        } else {
          cluster <- 0
        }
        thirdClusters1 <- c(thirdClusters1, cluster)
      } else {
        cluster <- secondClusters[1]
        if (cluster == 0 || length(clusterMeans2) == 0) {
          thirdClusters2 <- c(thirdClusters2, 0)
          next
        }
        temp <- lapply(1:nrow(clusterMeans2), function(x) return(clusterMeans2[x,] - (clusterMeans[cluster,] + clusterMeans[firstClusters[i],] - 2*clusterMeans[emptyDroplets,] + correctionFactor[[i]])))
        rowSums <- sapply(1:length(temp), function(x) return(sum(abs(unlist(temp[x])))))
        if (min(rowSums) < threshold) {
          minimum <- which.min(rowSums)
          cluster <- match(clusterMeans2[minimum], clusterMeans)
          clusterMeans2 <- clusterMeans2[-minimum, , drop=F]
        } else {
          cluster <- 0
        }
        thirdClusters2 <- c(thirdClusters2, cluster)
      }
    }
    
    thirdClusters <- c(thirdClusters2, thirdClusters1)
    
  }
  
  if(max(thirdClusters) == 0) return(list(clusters=thirdClusters))
  
  q <- quantile(counts[thirdClusters])
  
  if ((IQR(counts[thirdClusters]) > q[3]) && (q[3] > 5) && (q[1] < 5)) {
    remove <- c(remove, as.integer(names(which.min(counts[thirdClusters]))))
    print(counts[thirdClusters])
    return(findTertiaryClusters(emptyDroplets, firstClusters, secondClusters, remove, clusterMeans, correctionFactor, dimensions, counts))
  }
  if (min(counts[thirdClusters]) < q[2]/5) {
    remove <- c(remove, as.integer(names(which.min(counts[thirdClusters]))))
    print(counts[thirdClusters])
    return(findTertiaryClusters(emptyDroplets, firstClusters, secondClusters, remove, clusterMeans, correctionFactor, dimensions, counts))
  }
  if (max(counts[thirdClusters]) > q[4]*5) {
    remove <- c(remove, as.integer(names(which.max(counts[thirdClusters]))))
    print(counts[thirdClusters])
    return(findTertiaryClusters(emptyDroplets, firstClusters, secondClusters, remove, clusterMeans, correctionFactor, dimensions, counts))
  }
  
  list(clusters=thirdClusters)
}

# Find the quarternary cluster based on its position.
findQuaternaryCluster <- function(clusterMeans, emptyDroplets, remove, firstClusters, secondClusters=0, thirdClusters=0) {
  ## TODO really find the cluster?
  clusterMeans2 <- clusterMeans[-c(emptyDroplets, firstClusters, secondClusters, thirdClusters), ,drop=F]
  if (length(clusterMeans2) == 0 ) return(0)
  rowSums <- sapply(1:nrow(clusterMeans), function(x) return(sum(clusterMeans[x,])))
  rowSums2 <- sapply(1:nrow(clusterMeans2), function(x) return(sum(clusterMeans2[x,])))
  if (match(max(rowSums), rowSums) == match(max(rowSums2), rowSums)) {
    return(which.max(rowSums))
  } else {
    return(0)
  }
}

# Merge clusters that a close to each other, based on the parameter p.
mergeClusters <- function(result, clusterMeans, finalResult, remove, p = 12) {
  #threshold <-  dimensions/(length(finalResult)*2)
  clusterMeans2 <- clusterMeans[-c(finalResult, remove), , drop=F]
  if(nrow(clusterMeans2) > 0) {
    for (i in 1:nrow(clusterMeans2)) {
      for (j in 1:nrow(clusterMeans)) {
        #threshold <- clusterMeans[j,]/as.integer(sqrt(length(finalResult))*1.5)
        threshold <- sqrt(clusterMeans[j,])*p
        if (abs(clusterMeans2[i,1] - clusterMeans[j,1]) < threshold[1] && abs(clusterMeans2[i,2] - clusterMeans[j,2]) < threshold[2]) {
          #print(paste("merging", ColoursUsed[match(clusterMeans2[i], clusterMeans)], ColoursUsed[j]))
          result[result == match(clusterMeans2[i], clusterMeans)] = j
        }
      }
    }
  }
  return(result)
}

# Adjust the cluster centres based on the local density function
adjustClusterMeans <- function(data, clusterMeans, result, clusters) {
  data_dir <- system.file("extdata", package = "flowDensity")
  load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE)) # load f to copy over later so we have an FCS file to use flowDensity
  for (i in clusters) {
    if(i == 0) next
    f@exprs <- as.matrix(data[result == i,])
    Xc <- flowDensity:::.getPeaks(density(f@exprs[,1], width=1000), tinypeak.removal=0.2)$Peaks[1]
    Yc <- flowDensity:::.getPeaks(density(f@exprs[,2], width=1000), tinypeak.removal=0.2)$Peaks[1]
    if(!is.null(Xc) && !is.null(Yc)) {
      clusterMeans[i,] <- c(Xc, Yc)
    }
  }
  return(clusterMeans)
}

# Find the rain and assign it based on the distance to vector lines connecting the cluster centres.
assignRain <- function(clusterMeans, data, result, emptyDroplets, firstClusters, secondClusters, thirdClusters, fourthCluster, flowDensity) {
  remove <- vector()
  sdeviations <- list()
  rownames(data) <- 1:nrow(data)
  scalingHelper <- mean(scalingParam/4)
  ## empty to primary clusters:
  sdEmpties <- c(sd(data[result==emptyDroplets,1], na.rm = T) + scalingParam[1], sd(data[result==emptyDroplets,2], na.rm = T)+scalingParam[2])
  
  # Calculate standard deviations:
  newData <- subset(data, !result %in% c(secondClusters,thirdClusters,fourthCluster))
  newData <- subset(newData, !(newData[,1] < (clusterMeans[emptyDroplets,1]+sdEmpties[1]) & newData[,2] < (clusterMeans[emptyDroplets,2]+sdEmpties[2])))
  for (c in firstClusters) {
    if (c==0) next
    sdeviation <- 2*c(sd(data[which(result == c),1]), 
                      sd(data[which(result == c),2]))
    if (anyNA(sdeviation)) sdeviation <- scalingParam/2
    if (any(sdeviation > scalingParam)) sdeviation <- scalingParam
    sdeviations[[c]] <- sdeviation
  }
  for (c in secondClusters) {
    if (c==0) next
    sdeviation <- 2*c(sd(data[which(result == c),1]), 
                      sd(data[which(result == c),2]))
    if (anyNA(sdeviation)) sdeviation <- scalingParam/2
    if (any(sdeviation > scalingParam)) sdeviation <- scalingParam
    sdeviations[[c]] <- sdeviation
  }
  for (c in thirdClusters) {
    if (c==0) next
    sdeviation <- 2*c(sd(data[which(result == c),1]), 
                      sd(data[which(result == c),2]))
    if (anyNA(sdeviation)) sdeviation <- scalingParam/2
    if (any(sdeviation > scalingParam)) sdeviation <- scalingParam
    sdeviations[[c]] <- sdeviation
  }
  
  # Go through each cluster:
  for (c in firstClusters) {
    if (c==0) next
    sdC <- sdeviations[[c]]
    result[as.numeric(rownames(newData[(newData[,1] < clusterMeans[c,1]+sdC[1] & newData[,1] > clusterMeans[c,1]-sdC[1]
                                        & newData[,2] < clusterMeans[c,2]+sdC[2] & newData[,2] > clusterMeans[c,2]-sdC[2]),]))] <- c
    newData <- subset(newData, !(newData[,1] > clusterMeans[c,1]-sdC[1] & newData[,2] > clusterMeans[c,2]-sdC[2]))
    #    newData <- newData[-intersect(as.numeric(rownames(newData)), which(result==c)),]
  }
  if (nrow(newData) > 0) {
    m <- matrix(nrow = nrow(newData), ncol = length(firstClusters))
    for (c in 1:length(firstClusters)) {
      if (firstClusters[c]==0) next
      for (row in 1:nrow(newData)) {
        m[row, c] <- distToLineSegment(as.numeric(newData[row,]), as.numeric(clusterMeans[emptyDroplets,]), as.numeric(clusterMeans[firstClusters[c],]))
      }
    }
    minimalDistance <- apply(m, 1, function (x) which(x == min(x, na.rm=T))[1])
    for (r in 1:length(minimalDistance)) {
      distPoint <- euc.dist(newData[r,], clusterMeans[emptyDroplets,])
      distTotal <- distPoint + euc.dist(newData[r,], clusterMeans[firstClusters[minimalDistance[r]],])
      if (distPoint <= 0.2*distTotal) {
        result[as.numeric(rownames(newData)[r])] <- emptyDroplets
      } else {
        if(ncol(m) > 1) {
          minAnd2ndMin <- sort(m[r,], index.return=T )
          if (minAnd2ndMin$x[1]/minAnd2ndMin$x[2] >= 0.95
              || (minAnd2ndMin$x[1] <= scalingHelper && minAnd2ndMin$x[2] <= scalingHelper)) {
            remove <- c(remove, as.numeric(rownames(newData)[r]))
          }
        }
        result[as.numeric(rownames(newData)[r])] <- firstClusters[minimalDistance[r]]
      }
    }
  }
  
  ## primary to secondary clusters (3-plex and 4-plex only):
  if(!is.null(secondClusters) && sum(secondClusters) > 0) {
    if(flowDensity) {
      newData <- subset(data, !result %in% c(emptyDroplets,thirdClusters,fourthCluster))
      for (c in firstClusters) {
        if (c==0) next
        sdC <- sdeviations[[c]]
        newData <- subset(newData, !(newData[,1] < (clusterMeans[c,1]+sdC[1]) & newData[,2] < (clusterMeans[c,2]+sdC[2])))
        #   # if (clusterMeans[c,1] < clusterMeans[c,2]) {
        #   #   newData <- subset(newData, !(newData[,1] < (clusterMeans[c,1]+3*sdC[1]) & newData[,2] < (clusterMeans[c,2]-2*sdC[2]))) 
        #   # } else {
        #   #   newData <- subset(newData, !(newData[,1] < (clusterMeans[c,1]-2*sdC[1]) & newData[,2] < (clusterMeans[c,2]+3*sdC[2]))) 
        #   # }
      } 
    } else {
      newData <- subset(data, !result %in% c(emptyDroplets,firstClusters,thirdClusters,fourthCluster))
    }
    for (c in secondClusters) {
      if (c==0) next
      sdC <- sdeviations[[c]]
      result[as.numeric(rownames(newData[(newData[,1] < clusterMeans[c,1]+sdC[1] & newData[,1] > clusterMeans[c,1]-sdC[1]
                                          & newData[,2] < clusterMeans[c,2]+sdC[2] & newData[,2] > clusterMeans[c,2]-sdC[2]),]))] <- c
      newData <- subset(newData, !(newData[,1] > clusterMeans[c,1]-sdC[1] & newData[,2] > clusterMeans[c,2]-sdC[2]))
    }
    if (nrow(newData) > 0) {
      i <- 1
      j <- 2
      m <- matrix(nrow = nrow(newData), ncol = 2*length(secondClusters))
      for (c in 1:length(secondClusters)) {
        if (secondClusters[c]==0) next
        for (row in 1:nrow(newData)) {
          m[row, c] <- distToLineSegment(as.numeric(newData[row,]), as.numeric(clusterMeans[firstClusters[i],]), as.numeric(clusterMeans[secondClusters[c],]))
          m[row, c+length(secondClusters)] <- distToLineSegment(as.numeric(newData[row,]), as.numeric(clusterMeans[firstClusters[j],]), as.numeric(clusterMeans[secondClusters[c],]))
        }
        j <- j+1
        if(j>length(firstClusters)) {
          i <- i + 1
          j <- i + 1
        }
      }
      minimalDistance <- apply(m, 1, function (x) which(x == min(x, na.rm=T))[1])
      
      m2 <- matrix(nrow = nrow(newData), ncol = length(firstClusters) + length(secondClusters))
      for (row in 1:nrow(newData)) {
        for (cl1 in 1:length(firstClusters)) {
          if (firstClusters[cl1]==0) next
          m2[row, cl1] <- distToRect(as.numeric(clusterMeans[firstClusters[cl1],])-2*sdeviations[[firstClusters[cl1]]], as.numeric(clusterMeans[firstClusters[cl1],])+2*sdeviations[[firstClusters[cl1]]], as.numeric(newData[row,]))
        }
        for (cl2 in 1:length(secondClusters)) {
          if (secondClusters[cl2]==0) next
          col <- cl2+length(firstClusters)
          m2[row, col] <- distToRect(as.numeric(clusterMeans[secondClusters[cl2],])-2*sdeviations[[secondClusters[cl2]]], as.numeric(clusterMeans[secondClusters[cl2],])+2*sdeviations[[secondClusters[cl2]]], as.numeric(newData[row,]))
        }
      }
      minimalClDistance <- apply(m2, 1, function (x) which(x == min(x, na.rm=T))[1])
      
      
      helper <- c(1,1,1,2,2,3,2,3,4,3,4,4)
      if (length(secondClusters) == 3) helper <- c(1,1,2,2,3,3)
      helper2 <- rep(1:length(secondClusters), 2)
      helper5 <- c(firstClusters, secondClusters)
      
      for (r in 1:length(minimalDistance)) {
        
        if(m2[r,minimalClDistance[r]] <= m[r,minimalDistance[r]]) {
          result[as.numeric(rownames(newData)[r])] <- helper5[minimalClDistance[r]]
          next
        }
        
        distPoint <- euc.dist(newData[r,], clusterMeans[firstClusters[helper[minimalDistance[r]]],]+sdeviations[[firstClusters[helper[minimalDistance[r]]]]])
        distTotal <- distPoint + euc.dist(newData[r,], clusterMeans[secondClusters[helper2[minimalDistance[r]]],])
        
        if (distPoint <= 0.2*distTotal) {
          result[as.numeric(rownames(newData)[r])] <- firstClusters[helper[minimalDistance[r]]]
        } else {
          if(ncol(m) > 1) {
            minAnd2ndMin <- sort(m[r,], index.return=T )
            if ((minAnd2ndMin$x[1]/minAnd2ndMin$x[2] >= 0.95 && helper2[minAnd2ndMin$ix[1]] != helper2[minAnd2ndMin$ix[2]])
                || ((minAnd2ndMin$x[1] <= scalingHelper && minAnd2ndMin$x[2] <= scalingHelper) && helper2[minAnd2ndMin$ix[1]] != helper2[minAnd2ndMin$ix[2]])) {
              remove <- c(remove, as.numeric(rownames(newData)[r]))
            }
          }
          result[as.numeric(rownames(newData)[r])] <- secondClusters[helper2[minimalDistance[r]]]
        }
      }
    }
  }
  
  ## secondary to tertiary clusters (4-plex only):
  if(!is.null(thirdClusters) && sum(thirdClusters) > 0) {
    if(flowDensity) {
      newData <- subset(data, !result %in% c(emptyDroplets,firstClusters,fourthCluster))
      for (c in secondClusters) {
        if (c==0) next
        sdC <- sdeviations[[c]]
        newData <- subset(newData, !(newData[,1] < (clusterMeans[c,1]+sdC[1]) & newData[,2] < (clusterMeans[c,2]+sdC[2])))
        #   # if (clusterMeans[c,1] < clusterMeans[c,2]) {
        #   #   newData <- subset(newData, !(newData[,1] < (clusterMeans[c,1]+2*sdC[1]) & newData[,2] < (clusterMeans[c,2]-sdC[2]))) 
        #   # } else {
        #   #   newData <- subset(newData, !(newData[,1] < (clusterMeans[c,1]-sdC[1]) & newData[,2] < (clusterMeans[c,2]+2*sdC[2]))) 
        #   # }
      }
    } else {
      newData <- subset(data, !result %in% c(emptyDroplets,firstClusters,secondClusters,fourthCluster))
    }
    for (c in thirdClusters) {
      if (c==0) next
      sdC <- sdeviations[[c]]
      result[as.numeric(rownames(newData[(newData[,1] < clusterMeans[c,1]+sdC[1] & newData[,1] > clusterMeans[c,1]-sdC[1]
                                          & newData[,2] < clusterMeans[c,2]+sdC[2] & newData[,2] > clusterMeans[c,2]-sdC[2]),]))] <- c
      newData <- subset(newData, !(newData[,1] > clusterMeans[c,1]-sdC[1] & newData[,2] > clusterMeans[c,2]-sdC[2]))
    }
    if (nrow(newData) > 0) {
      m <- matrix(nrow = nrow(newData), ncol = 3*length(thirdClusters))
      helper3 <- c(1,1,2,4,2,3,3,5,4,5,6,6)
      helper4 <- rep(1:length(thirdClusters), 3)
      for (c in 1:length(thirdClusters)) {
        if (thirdClusters[c]==0) next
        for (row in 1:nrow(newData)) {
          m[row, c] <- distToLineSegment(as.numeric(newData[row,]), as.numeric(clusterMeans[secondClusters[helper3[c]],]), as.numeric(clusterMeans[thirdClusters[c],]))
          m[row, c+length(thirdClusters)] <- distToLineSegment(as.numeric(newData[row,]), as.numeric(clusterMeans[secondClusters[helper3[c+length(thirdClusters)]],]), as.numeric(clusterMeans[thirdClusters[c],]))
          m[row, c+2*length(thirdClusters)] <- distToLineSegment(as.numeric(newData[row,]), as.numeric(clusterMeans[secondClusters[helper3[c+2*length(thirdClusters)]],]), as.numeric(clusterMeans[thirdClusters[c],]))
        }
      }
      minimalDistance <- apply(m, 1, function (x) which(x == min(x, na.rm=T))[1])
      m2 <- matrix(nrow = nrow(newData), ncol = length(secondClusters) + length(thirdClusters))
      for (row in 1:nrow(newData)) {
        for (cl1 in 1:length(secondClusters)) {
          if (secondClusters[cl1]==0) next
          m2[row, cl1] <- distToRect(as.numeric(clusterMeans[secondClusters[cl1],])-2*sdeviations[[secondClusters[cl1]]], as.numeric(clusterMeans[secondClusters[cl1],])+2*sdeviations[[secondClusters[cl1]]], as.numeric(newData[row,]))
        }
        for (cl2 in 1:length(thirdClusters)) {
          if (thirdClusters[cl2]==0) next
          col <- cl2+length(secondClusters)
          m2[row, col] <- distToRect(as.numeric(clusterMeans[thirdClusters[cl2],])-2*sdeviations[[thirdClusters[cl2]]], as.numeric(clusterMeans[thirdClusters[cl2],])+2*sdeviations[[thirdClusters[cl2]]], as.numeric(newData[row,]))
        }
      }
      minimalClDistance <- apply(m2, 1, function (x) which(x == min(x, na.rm=T))[1])
      
      helper5 <- c(secondClusters, thirdClusters)
      for (r in 1:length(minimalDistance)) {
        if(m2[r,minimalClDistance[r]] <= m[r,minimalDistance[r]]) {
          result[as.numeric(rownames(newData)[r])] <- helper5[minimalClDistance[r]]
          next
        }
        distPoint <- euc.dist(newData[r,], clusterMeans[secondClusters[helper3[minimalDistance[r]]],]+sdeviations[[secondClusters[helper3[minimalDistance[r]]]]])
        distTotal <- distPoint + euc.dist(newData[r,], clusterMeans[thirdClusters[helper4[minimalDistance[r]]],])
        if (distPoint <= 0.2*distTotal) {
          result[as.numeric(rownames(newData)[r])] <- secondClusters[helper3[minimalDistance[r]]]
        } else {
          if(ncol(m) > 1) {
            minAnd2ndMin <- sort(m[r,], index.return=T )
            if ((minAnd2ndMin$x[1]/minAnd2ndMin$x[2] >= 0.95 && helper4[minAnd2ndMin$ix[1]] != helper4[minAnd2ndMin$ix[2]])
                || ((minAnd2ndMin$x[1] <= scalingHelper && minAnd2ndMin$x[2] <= scalingHelper) && helper4[minAnd2ndMin$ix[1]] != helper4[minAnd2ndMin$ix[2]])) {
              remove <- c(remove, as.numeric(rownames(newData)[r]))
            }
          }
          result[as.numeric(rownames(newData)[r])] <- thirdClusters[helper4[minimalDistance[r]]]
        }
      }
    }
  }
  
  ## tertiary to quaternary cluster (4-plex), secondary to quarternary (3-plex), primary to quarternary (2-plex) :
  if(!is.null(fourthCluster) && sum(fourthCluster) > 0) {
    if(!is.null(thirdClusters) && sum(thirdClusters) > 0) {
      if(flowDensity) {
        newData <- subset(data, !result %in% c(emptyDroplets, firstClusters, secondClusters))
      } else {
        newData <- subset(data, !result %in% c(emptyDroplets, firstClusters, secondClusters, thirdClusters))
      }
      prevClusters <- thirdClusters
    } else if(!is.null(secondClusters) && sum(secondClusters) > 0) {
      if(flowDensity) {
        newData <- subset(data, !result %in% c(emptyDroplets, firstClusters))
      } else {
        newData <- subset(data, !result %in% c(emptyDroplets, firstClusters, secondClusters))
      }
      prevClusters <- secondClusters
    } else if(!is.null(firstClusters) && sum(firstClusters) > 0){
      if(flowDensity) {
        newData <- subset(data, !result %in% c(emptyDroplets))
      } else {
        newData <- subset(data, !result %in% c(emptyDroplets, firstClusters))
      }
      prevClusters <- firstClusters
    } else {
      return(list(result=result, remove=unique(remove)))
    }
    for (c in prevClusters) {
      if (c==0) next
      sdC <- sdeviations[[c]]
      newData <- subset(newData, !(newData[,1] < (clusterMeans[c,1]+sdC[1]) & newData[,2] < (clusterMeans[c,2]+sdC[2])))
    }
    for (c in fourthCluster) {
      sdC <- 2*c(sd(data[result==c,1], na.rm = T), sd(data[result==c,2], na.rm = T))
      if (is.na(sum(sdC))) next
      result[as.numeric(rownames(newData[(newData[,1] < clusterMeans[c,1]+sdC[1] & newData[,1] > clusterMeans[c,1]-sdC[1]
                                          & newData[,2] < clusterMeans[c,2]+sdC[2] & newData[,2] > clusterMeans[c,2]-sdC[2]),]))] <- c
      newData <- subset(newData, !(newData[,1] < clusterMeans[c,1]+sdC[1] & newData[,1] > clusterMeans[c,1]-sdC[1]
                                   & newData[,2] < clusterMeans[c,2]+sdC[2] & newData[,2] > clusterMeans[c,2]-sdC[2]))
    }
    if (nrow(newData) > 0) {
      m <- matrix(nrow = nrow(newData), ncol = length(prevClusters))
      for (c in 1:length(prevClusters)) {
        # if (prevClusters[c]==0) next
        for (row in 1:nrow(newData)) {
          m[row, c] <- distToLineSegment(as.numeric(newData[row,]), as.numeric(clusterMeans[fourthCluster,]), as.numeric(clusterMeans[prevClusters[c],]))
        }
      }
      minimalDistance <- apply(m, 1, function (x) which(x == min(x, na.rm=T))[1])
      for (r in 1:length(minimalDistance)) {
        distPoint <- euc.dist(newData[r,], clusterMeans[prevClusters[minimalDistance[r]],]+sdeviations[[prevClusters[minimalDistance[r]]]])
        distTotal <- distPoint + euc.dist(newData[r,], clusterMeans[fourthCluster,])
        if (distPoint <= 0.2*distTotal) {
          result[as.numeric(rownames(newData)[r])] <- prevClusters[minimalDistance[r]]
        } else {
          if(ncol(m) > 1) {
            minAnd2ndMin <- sort(m[r,], index.return=T )
            if (minAnd2ndMin$x[1]/minAnd2ndMin$x[2] >= 0.95
                || (minAnd2ndMin$x[1] <= scalingHelper && minAnd2ndMin$x[2] <= scalingHelper)) {
              remove <- c(remove, as.numeric(rownames(newData)[r]))
            }
          }
          result[as.numeric(rownames(newData)[r])] <- fourthCluster
        }
      }
    }
  }
  list(result=result, remove=unique(remove))
}


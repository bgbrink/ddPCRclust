## Part of the ddPCRclust algorithm
## Author: Benedikt G Brink, Bielefeld University
## Contributor: Justin Meskas
## November 2017

# Find the primary clusters based on their density peaks found by flowDensity.
findPrimaryClustersDensity <- function(f, File, f_remNeg, NumOfMarkers, scalingParam, 
    epsilon) {
    threshold <- 5 * epsilon
    CutAbovePrimary <- mean(scalingParam) * 2
    
    # Calculate the threshold to remove events that are too positive. This makes it
    # easier to find single positives.
    f_findExtremes_temp <- f
    f_remNeg_temp <- f_remNeg
    
    
    # find left and right primary
    tinyP <- 0.6
    repeat {
        tinyP <- tinyP/2
        if (tinyP < epsilon) 
            break
        leftPrim <- deGate(f_remNeg_temp, 1, all.cuts = TRUE, tinypeak.removal = tinyP, 
            adjust.dens = 0.1)[1]
        indices <- which(flowCore::exprs(f_remNeg_temp)[, 1] <= leftPrim)
        f_leftPrim <- f_remNeg_temp
        flowCore::exprs(f_leftPrim) <- flowCore::exprs(f_leftPrim)[indices, ]
        x_leftPrim <- max(flowDensity::getPeaks(stats::density(flowCore::exprs(f_leftPrim)[, 
            1], width = 1000), tinypeak.removal = tinyP)$Peaks)
        y_leftPrim <- max(flowDensity::getPeaks(stats::density(flowCore::exprs(f_leftPrim)[, 
            2], width = 1000), tinypeak.removal = tinyP)$Peaks)
        if (x_leftPrim < (min(File[, 1]) + (max(File[, 1]) - min(File[, 1]))/6)) {
            threshold <- tinyP
            break
        }
    }
    
    tinyP <- 0.6
    repeat {
        tinyP <- tinyP/2
        if (tinyP < epsilon) 
            break
        rightPrim <- deGate(f_remNeg_temp, 2, all.cuts = TRUE, tinypeak.removal = tinyP, 
            adjust.dens = 0.1)[1]
        indices <- which(flowCore::exprs(f_remNeg_temp)[, 2] <= rightPrim)
        f_rightPrim <- f_remNeg_temp
        flowCore::exprs(f_rightPrim) <- flowCore::exprs(f_rightPrim)[indices, ]
        x_rightPrim <- max(flowDensity::getPeaks(stats::density(flowCore::exprs(f_rightPrim)[, 
            1], width = 1000), tinypeak.removal = tinyP)$Peaks)
        y_rightPrim <- max(flowDensity::getPeaks(stats::density(flowCore::exprs(f_rightPrim)[, 
            2], width = 1000), tinypeak.removal = tinyP)$Peaks)
        if (y_rightPrim < (min(File[, 2]) + (max(File[, 2]) - min(File[, 2]))/6)) {
            threshold <- min(threshold, tinyP)
            break
        }
    }
    
    mSlope <- (y_leftPrim - y_rightPrim)/(x_leftPrim - x_rightPrim)
    theta <- abs(atan(mSlope))
    rotate <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
    
    Rot_xy_leftPrim <- rotate %*% c(x_leftPrim, y_leftPrim)  # coordinates of rotated left  primary cluster
    Rot_xy_rightPrim <- rotate %*% c(x_rightPrim, y_rightPrim)  # coordinates of rotated right primary cluster
    
    flowCore::exprs(f_findExtremes_temp)[, c(1, 2)] <- t(rotate %*% t(flowCore::exprs(f_findExtremes_temp)[, 
        c(1, 2)]))
    flowCore::exprs(f_remNeg_temp)[, c(1, 2)] <- t(rotate %*% t(flowCore::exprs(f_remNeg_temp)[, 
        c(1, 2)]))
    upSlantmax <- deGate(f_findExtremes_temp, c(2), percentile = 0.999, use.percentile = TRUE)
    upSlantmin <- deGate(f_findExtremes_temp, c(2), percentile = 0.001, use.percentile = TRUE)
    ScaleChop <- (upSlantmax - upSlantmin)/max(File)
    # Chop off the very positive events
    indices <- which(flowCore::exprs(f_remNeg_temp)[, 2] <= (max(Rot_xy_leftPrim[2], 
        Rot_xy_rightPrim[2]) + CutAbovePrimary * ScaleChop))
    flowCore::exprs(f_remNeg_temp) <- flowCore::exprs(f_remNeg_temp)[indices, ]
    
    # find thresholds to divide the primary clusters
    highP <- 1
    lowP <- 0
    repeat {
        # newton iteration
        newP <- (highP + lowP)/2
        if (newP < 2 * epsilon) 
            break
        Cuts <- deGate(f_remNeg_temp, 1, all.cuts = TRUE, tinypeak.removal = newP, 
            adjust.dens = 0.1)
        if (length(Cuts) < (NumOfMarkers - 1)) {
            highP <- newP
        } else if (length(Cuts) > (NumOfMarkers - 1)) {
            lowP <- newP
        } else {
            break
        }
    }
    Xc <- Yc <- NULL
    Cuts <- c(min(flowCore::exprs(f_remNeg_temp)[, 1]), Cuts, max(flowCore::exprs(f_remNeg_temp)[, 
        1]))
    deviationX <- vector()
    deviationY <- vector()
    
    # find the coordinates of the primary clusters
    for (r1 in seq_len(length(Cuts) - 1)) {
        indices <- intersect(which(flowCore::exprs(f_remNeg_temp)[, 1] >= Cuts[r1]), 
            which(flowCore::exprs(f_remNeg_temp)[, 1] < Cuts[r1 + 1]))
        f_Cuts_temp <- f_remNeg_temp
        flowCore::exprs(f_Cuts_temp) <- flowCore::exprs(f_Cuts_temp)[indices, ]
        flowCore::exprs(f_Cuts_temp)[, c(1, 2)] <- t(t(rotate) %*% t(flowCore::exprs(f_Cuts_temp)[, 
            c(1, 2)]))
        Xx <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
            1], width = 1000), tinypeak.removal = newP)$Peaks
        Yy <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
            2], width = 1000), tinypeak.removal = newP)$Peaks
        if (length(Xx) > 1 || length(Yy) > 1) {
            highP <- 1
            lowP <- 0
            repeat {
                # newton iteration
                newP <- (highP + lowP)/2
                if (newP < epsilon) 
                  break
                Xx <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
                  1], width = 1000), tinypeak.removal = newP)$Peaks
                Yy <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
                  2], width = 1000), tinypeak.removal = newP)$Peaks
                if (length(Xx) > 1 || length(Yy) > 1) {
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
        deviationX <- c(deviationX, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
            1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), intersect(which(flowCore::exprs(f)[, 
            2] <= Yy + 1000), which(flowCore::exprs(f)[, 2] >= Yy - 1000))), 1]))
        deviationY <- c(deviationY, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
            1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), intersect(which(flowCore::exprs(f)[, 
            2] <= Yy + 1000), which(flowCore::exprs(f)[, 2] >= Yy - 1000))), 2]))
        Xc <- c(Xc, Xx)
        Yc <- c(Yc, Yy)
    }
    if (length(Xc) == 4 && length(Yc) == 4) {
        Rot_xy_midLeftPrim <- rotate %*% c(Xc[2], Yc[2])  # coordinates of rotated middle left  primary cluster
        Rot_xy_midRightPrim <- rotate %*% c(Xc[3], Yc[3])  # coordinates of rotated middle right primary cluster
        
        test_left <- abs(det(rbind(cbind(1, 1, 1), cbind(Rot_xy_leftPrim, Rot_xy_midLeftPrim, 
            rbind(0, 0)))))/100
        test_right <- abs(det(rbind(cbind(1, 1, 1), cbind(Rot_xy_rightPrim, Rot_xy_midRightPrim, 
            rbind(0, 0)))))/100
        
        if (!is.na(test_left) && test_left < max(File)/2 && test_left > max(File)/20) {
            Xc[1] <- mean(Xc[seq_len(2)])
            Yc[1] <- mean(Yc[seq_len(2)])
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
    return(list(clusters = cbind(Xc, Yc), deviation = cbind(deviationX, deviationY)))
}

# Find the secondary clusters based on their density peaks found by flowDensity.
findSecondaryClustersDensity <- function(f, file, f_remNegPrim, emptyDroplets, firstClusters, 
    scalingParam, epsilon) {
    CutAbovePrimary <- mean(scalingParam) * 2
    f_remNegPrim_chopTertQuat <- f_remNegPrim
    NumberOfSinglePos <- nrow(firstClusters$clusters)
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
    
    flowCore::exprs(f_remNegPrim_chopTertQuat)[, c(1, 2)] <- t(rotate %*% t(flowCore::exprs(f_remNegPrim_chopTertQuat)[, 
        c(1, 2)]))
    indices <- which(flowCore::exprs(f_remNegPrim_chopTertQuat)[, 2] <= (upSlantmin + 
        2 * ((max(Rot_xy_leftPrim[2], Rot_xy_rightPrim[2]) + CutAbovePrimary * ScaleChop) - 
            upSlantmin)))
    flowCore::exprs(f_remNegPrim_chopTertQuat) <- flowCore::exprs(f_remNegPrim_chopTertQuat)[indices, 
        ]
    flowCore::exprs(f_remNegPrim_chopTertQuat)[, c(1, 2)] <- t(t(rotate) %*% t(flowCore::exprs(f_remNegPrim_chopTertQuat)[, 
        c(1, 2)]))
    
    deviation2X <- deviation2Y <- vector()
    Xc2g <- Yc2g <- NULL
    
    if (nrow(f_remNegPrim) >= NumberOfSinglePos^2 * 6) {
        # If less than 96 (4-plex) points for doub, trip and quad pops together, then use
        # vector addition
        
        adjustDens <- 0.3
        repeat {
            # find thresholds to divide the secondary clusters
            highP <- 0.4
            lowP <- 0
            repeat {
                # tinyP from 0.4 down to 0.04
                tinyP <- (highP + lowP)/2
                if ((tinyP + epsilon) > highP) 
                  break
                theta <- 0
                repeat {
                  # theta from 0 to pi/2
                  f_rotate_tinyP_temp <- f_remNegPrim_chopTertQuat
                  flowCore::exprs(f_rotate_tinyP_temp)[, c(1, 2)] <- t(rotate %*% 
                    t(flowCore::exprs(f_rotate_tinyP_temp)[, c(1, 2)]))
                  Cuts <- deGate(f_rotate_tinyP_temp, 1, all.cuts = TRUE, tinypeak.removal = tinyP, 
                    adjust.dens = adjustDens)
                  if (length(Cuts) == (choose(NumberOfSinglePos, 2) - 1)) 
                    break
                  
                  if (theta >= pi/2) 
                    break
                  
                  theta <- theta + pi/16
                }
                if (length(Cuts) < (choose(NumberOfSinglePos, 2) - 1)) {
                  highP <- tinyP
                } else if (length(Cuts) > (choose(NumberOfSinglePos, 2) - 1)) {
                  lowP <- tinyP
                } else {
                  break
                }
            }
            if (length(Cuts) == (choose(NumberOfSinglePos, 2) - 1)) 
                break
            
            if (adjustDens > 0.5) {
                message("double positive detection failed.")
                break
                
            }
            adjustDens <- adjustDens + 0.05
        }
        if (length(Cuts) == (choose(NumberOfSinglePos, 2) - 1)) {
            # find coordiates of secondary populations
            Cuts <- c(min(flowCore::exprs(f_rotate_tinyP_temp)[, 1]), Cuts, max(flowCore::exprs(f_rotate_tinyP_temp)[, 
                1]))
            
            for (r1 in seq_len(length(Cuts) - 1)) {
                indices <- intersect(which(flowCore::exprs(f_rotate_tinyP_temp)[, 
                  1] >= Cuts[r1]), which(flowCore::exprs(f_rotate_tinyP_temp)[, 1] < 
                  Cuts[r1 + 1]))
                switch(r1, {
                  XxA <- (firstClusters$clusters[1, 1] + firstClusters$clusters[2, 
                    1] - emptyDroplets[1])
                  YyA <- (firstClusters$clusters[1, 2] + firstClusters$clusters[2, 
                    2] - emptyDroplets[2])
                }, {
                  XxA <- (firstClusters$clusters[1, 1] + firstClusters$clusters[3, 
                    1] - emptyDroplets[1])
                  YyA <- (firstClusters$clusters[1, 2] + firstClusters$clusters[3, 
                    2] - emptyDroplets[2])
                }, {
                  XxA <- (firstClusters$clusters[2, 1] + firstClusters$clusters[3, 
                    1] - emptyDroplets[1])
                  YyA <- (firstClusters$clusters[2, 2] + firstClusters$clusters[3, 
                    2] - emptyDroplets[2])
                }, {
                  XxA <- (firstClusters$clusters[1, 1] + firstClusters$clusters[4, 
                    1] - emptyDroplets[1])
                  YyA <- (firstClusters$clusters[1, 2] + firstClusters$clusters[4, 
                    2] - emptyDroplets[2])
                }, {
                  XxA <- (firstClusters$clusters[2, 1] + firstClusters$clusters[4, 
                    1] - emptyDroplets[1])
                  YyA <- (firstClusters$clusters[2, 2] + firstClusters$clusters[4, 
                    2] - emptyDroplets[2])
                }, {
                  XxA <- (firstClusters$clusters[3, 1] + firstClusters$clusters[4, 
                    1] - emptyDroplets[1])
                  YyA <- (firstClusters$clusters[3, 2] + firstClusters$clusters[4, 
                    2] - emptyDroplets[2])
                })
                if (length(indices) < 2) {
                  deviation2X <- c(deviation2X, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                    1] <= XxA + 1000), which(flowCore::exprs(f)[, 1] >= XxA - 1000)), 
                    intersect(which(flowCore::exprs(f)[, 2] <= YyA + 1000), which(flowCore::exprs(f)[, 
                      2] >= YyA - 1000))), 1]))
                  deviation2Y <- c(deviation2Y, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                    1] <= XxA + 1000), which(flowCore::exprs(f)[, 1] >= XxA - 1000)), 
                    intersect(which(flowCore::exprs(f)[, 2] <= YyA + 1000), which(flowCore::exprs(f)[, 
                      2] >= YyA - 1000))), 2]))
                  Xc2g <- c(Xc2g, XxA)
                  Yc2g <- c(Yc2g, YyA)
                  next
                }
                f_Cuts_temp <- f_rotate_tinyP_temp
                flowCore::exprs(f_Cuts_temp) <- flowCore::exprs(f_Cuts_temp)[indices, 
                  ]
                flowCore::exprs(f_Cuts_temp)[, c(1, 2)] <- t(t(rotate) %*% t(flowCore::exprs(f_Cuts_temp)[, 
                  c(1, 2)]))
                
                highPXx <- 1
                lowPXx <- 0
                repeat {
                  tinyPXx <- (highPXx + lowPXx)/2
                  if ((tinyPXx + epsilon) > highPXx) {
                    Xx <- Xx[1]
                    message(paste0("Error finding X value for peak ", r1, " of ", 
                      choose(NumberOfSinglePos, 2)))
                    break
                  }
                  Xx <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
                    1], width = 1000), tinypeak.removal = tinyPXx)$Peaks
                  if (length(Xx) > 1) {
                    lowPXx <- tinyPXx
                  } else if (length(Xx) < 1) {
                    highPXx <- tinyPXx
                  } else {
                    break
                  }
                }
                highPYy <- 1
                lowPYy <- 0
                repeat {
                  tinyPYy <- (highPYy + lowPYy)/2
                  if ((tinyPYy + epsilon) > highPXx) {
                    Yy <- Yy[1]
                    message(paste0("Error finding Y value for peak ", r1, " of ", 
                      choose(NumberOfSinglePos, 2)))
                    break
                  }
                  Yy <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
                    2], width = 1000), tinypeak.removal = tinyPYy)$Peaks
                  if (length(Yy) > 1) {
                    lowPYy <- tinyPYy
                  } else if (length(Yy) < 1) {
                    highPYy <- tinyPYy
                  } else {
                    break
                  }
                }
                if (abs(Xx - XxA) > 3 * scalingParam[1] || abs(Yy - YyA) > 3 * scalingParam[2]) {
                  deviation2X <- c(deviation2X, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                    1] <= XxA + 1000), which(flowCore::exprs(f)[, 1] >= XxA - 1000)), 
                    intersect(which(flowCore::exprs(f)[, 2] <= YyA + 1000), which(flowCore::exprs(f)[, 
                      2] >= YyA - 1000))), 1]))
                  deviation2Y <- c(deviation2Y, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                    1] <= XxA + 1000), which(flowCore::exprs(f)[, 1] >= XxA - 1000)), 
                    intersect(which(flowCore::exprs(f)[, 2] <= YyA + 1000), which(flowCore::exprs(f)[, 
                      2] >= YyA - 1000))), 2]))
                  Xc2g <- c(Xc2g, XxA)
                  Yc2g <- c(Yc2g, YyA)
                } else {
                  deviation2X <- c(deviation2X, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                    1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), 
                    intersect(which(flowCore::exprs(f)[, 2] <= Yy + 1000), which(flowCore::exprs(f)[, 
                      2] >= Yy - 1000))), 1]))
                  deviation2Y <- c(deviation2Y, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                    1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), 
                    intersect(which(flowCore::exprs(f)[, 2] <= Yy + 1000), which(flowCore::exprs(f)[, 
                      2] >= Yy - 1000))), 2]))
                  Xc2g <- c(Xc2g, Xx)
                  Yc2g <- c(Yc2g, Yy)
                }
            }
        } else {
            # vector additon
            for (r1 in seq_len(NumberOfSinglePos)) {
                for (r2 in seq_len(NumberOfSinglePos)) {
                  if (r1 >= r2) 
                    next
                  Xx <- firstClusters$clusters[r1, 1] + firstClusters$clusters[r2, 
                    1] - emptyDroplets[1]
                  Yy <- firstClusters$clusters[r1, 2] + firstClusters$clusters[r2, 
                    2] - emptyDroplets[2]
                  deviation2X <- c(deviation2X, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                    1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), 
                    intersect(which(flowCore::exprs(f)[, 2] <= Yy + 1000), which(flowCore::exprs(f)[, 
                      2] >= Yy - 1000))), 1]))
                  deviation2Y <- c(deviation2Y, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                    1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), 
                    intersect(which(flowCore::exprs(f)[, 2] <= Yy + 1000), which(flowCore::exprs(f)[, 
                      2] >= Yy - 1000))), 2]))
                  Xc2g <- c(Xc2g, Xx)
                  Yc2g <- c(Yc2g, Yy)
                }
            }
        }
    } else {
        # vector additon
        for (r1 in seq_len(NumberOfSinglePos)) {
            for (r2 in seq_len(NumberOfSinglePos)) {
                if (r1 >= r2) 
                  next
                Xx <- firstClusters$clusters[r1, 1] + firstClusters$clusters[r2, 
                  1] - emptyDroplets[1]
                Yy <- firstClusters$clusters[r1, 2] + firstClusters$clusters[r2, 
                  2] - emptyDroplets[2]
                deviation2X <- c(deviation2X, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                  1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), 
                  intersect(which(flowCore::exprs(f)[, 2] <= Yy + 1000), which(flowCore::exprs(f)[, 
                    2] >= Yy - 1000))), 1]))
                deviation2Y <- c(deviation2Y, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                  1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), 
                  intersect(which(flowCore::exprs(f)[, 2] <= Yy + 1000), which(flowCore::exprs(f)[, 
                    2] >= Yy - 1000))), 2]))
                Xc2g <- c(Xc2g, Xx)
                Yc2g <- c(Yc2g, Yy)
            }
        }
    }
    # switch ordering of secondary populations if if staement is true
    if (length(Xc2g) > 3 && (Xc2g[3]^2 + Yc2g[3]^2) >= 1.2 * (Xc2g[4]^2 + Yc2g[4]^2)) {
        message("Switching the middle two double populations.")
        Xc2g[c(3, 4)] <- Xc2g[c(4, 3)]
        Yc2g[c(3, 4)] <- Yc2g[c(4, 3)]
    }
    if (length(Xc2g) > 1 && Xc2g[1] > Xc2g[2]) {
        message("Switching the top left two double populations.")
        Xc2g[c(1, 2)] <- Xc2g[c(2, 1)]
        Yc2g[c(1, 2)] <- Yc2g[c(2, 1)]
    }
    if (length(Yc2g) > 5 && Yc2g[5] < Yc2g[6]) {
        message("Switching the bottom right two double populations.")
        Xc2g[c(5, 6)] <- Xc2g[c(6, 5)]
        Yc2g[c(5, 6)] <- Yc2g[c(6, 5)]
    }
    
    return(list(clusters = cbind(Xc2g, Yc2g), deviation = cbind(deviation2X, deviation2Y)))
}

# Find the tertiary clusters based on their density peaks found by flowDensity.
findTertiaryClustersDensity <- function(f, f_remNegPrimSec, emptyDroplets, firstClusters, 
    secondaryClusters, scalingParam, epsilon) {
    NumberOfSinglePos <- nrow(firstClusters$clusters)
    XcTert <- YcTert <- NULL
    deviation2X <- deviation2Y <- vector()
    
    x_leftPrim <- firstClusters$clusters[1, 1]
    y_leftPrim <- firstClusters$clusters[1, 2]
    x_rightPrim <- firstClusters$clusters[NumberOfSinglePos, 1]
    y_rightPrim <- firstClusters$clusters[NumberOfSinglePos, 2]
    
    mSlope <- (y_leftPrim - y_rightPrim)/(x_leftPrim - x_rightPrim)
    theta <- abs(atan(mSlope))
    rotate <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
    
    if (nrow(f_remNegPrimSec) >= NumberOfSinglePos^2 * 3) {
        # If less than 50 points for trip and quad pops together, then use vector
        # addition
        
        # find thresholds to divide the tertiary clusters
        highP <- 1
        lowP <- 0
        repeat {
            # newton iteration
            newP <- (highP + lowP)/2
            Cuts <- deGate(f_remNegPrimSec, 1, all.cuts = TRUE, tinypeak.removal = newP, 
                adjust.dens = 0.1)
            if (length(Cuts) < (NumberOfSinglePos - 1)) {
                highP <- newP
            } else if (length(Cuts) > (NumberOfSinglePos - 1)) {
                lowP <- newP
            } else {
                break
            }
        }
        
        Cuts <- c(min(flowCore::exprs(f_remNegPrimSec)[, 1]), Cuts, max(flowCore::exprs(f_remNegPrimSec)[, 
            1]))
        
        # find coordiates of tertiary populations
        for (r1 in seq_len(length(Cuts) - 1)) {
            indices <- intersect(which(flowCore::exprs(f_remNegPrimSec)[, 1] >= Cuts[r1]), 
                which(flowCore::exprs(f_remNegPrimSec)[, 1] < Cuts[r1 + 1]))
            
            switch(r1, {
                XxA <- mean(secondaryClusters$clusters[1, 1] + firstClusters$clusters[3, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[2, 1] + firstClusters$clusters[2, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[3, 1] + firstClusters$clusters[1, 
                  1] - emptyDroplets[1])
                YyA <- mean(secondaryClusters$clusters[1, 2] + firstClusters$clusters[3, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[2, 2] + firstClusters$clusters[2, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[3, 2] + firstClusters$clusters[1, 
                  2] - emptyDroplets[2])
            }, {
                XxA <- mean(secondaryClusters$clusters[1, 1] + firstClusters$clusters[4, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[4, 1] + firstClusters$clusters[2, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[5, 1] + firstClusters$clusters[1, 
                  1] - emptyDroplets[1])
                YyA <- mean(secondaryClusters$clusters[1, 2] + firstClusters$clusters[4, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[4, 2] + firstClusters$clusters[2, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[5, 2] + firstClusters$clusters[1, 
                  2] - emptyDroplets[2])
            }, {
                XxA <- mean(secondaryClusters$clusters[2, 1] + firstClusters$clusters[4, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[4, 1] + firstClusters$clusters[3, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[6, 1] + firstClusters$clusters[1, 
                  1] - emptyDroplets[1])
                YyA <- mean(secondaryClusters$clusters[2, 2] + firstClusters$clusters[4, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[4, 2] + firstClusters$clusters[3, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[6, 2] + firstClusters$clusters[1, 
                  2] - emptyDroplets[2])
            }, {
                XxA <- mean(secondaryClusters$clusters[3, 1] + firstClusters$clusters[4, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[5, 1] + firstClusters$clusters[3, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[6, 1] + firstClusters$clusters[2, 
                  1] - emptyDroplets[1])
                YyA <- mean(secondaryClusters$clusters[3, 2] + firstClusters$clusters[4, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[5, 2] + firstClusters$clusters[3, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[6, 2] + firstClusters$clusters[2, 
                  2] - emptyDroplets[2])
            })
            if (length(indices) < 2) {
                deviation2X <- c(deviation2X, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                  1] <= XxA + 1000), which(flowCore::exprs(f)[, 1] >= XxA - 1000)), 
                  intersect(which(flowCore::exprs(f)[, 2] <= YyA + 1000), which(flowCore::exprs(f)[, 
                    2] >= YyA - 1000))), 1]))
                deviation2Y <- c(deviation2Y, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                  1] <= XxA + 1000), which(flowCore::exprs(f)[, 1] >= XxA - 1000)), 
                  intersect(which(flowCore::exprs(f)[, 2] <= YyA + 1000), which(flowCore::exprs(f)[, 
                    2] >= YyA - 1000))), 2]))
                XcTert <- c(XcTert, XxA)
                YcTert <- c(YcTert, YyA)
                next
            }
            f_remNegPrimSec_temp <- f_remNegPrimSec
            flowCore::exprs(f_remNegPrimSec_temp) <- flowCore::exprs(f_remNegPrimSec_temp)[indices, 
                ]
            flowCore::exprs(f_remNegPrimSec_temp)[, c(1, 2)] <- t(t(rotate) %*% t(flowCore::exprs(f_remNegPrimSec_temp)[, 
                c(1, 2)]))
            
            highPXx <- 1
            lowPXx <- 0
            repeat {
                tinyPXx <- (highPXx + lowPXx)/2
                if ((tinyPXx + epsilon) > highPXx) {
                  Xx <- Xx[1]
                  message(paste0("Default X value for tertiary peak ", r1, " of 4"))
                  break
                }
                Xx <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_remNegPrimSec_temp)[, 
                  1], width = 1000), tinypeak.removal = tinyPXx)$Peaks
                if (length(Xx) > 1) {
                  lowPXx <- tinyPXx
                } else if (length(Xx) < 1) {
                  highPXx <- tinyPXx
                } else {
                  break
                }
            }
            highPYy <- 1
            lowPYy <- 0
            repeat {
                tinyPYy <- (highPYy + lowPYy)/2
                if ((tinyPYy + epsilon) > highPYy) {
                  Yy <- Yy[1]
                  message(paste0("Default Y value for tertiary peak ", r1, " of 4"))
                  break
                }
                Yy <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_remNegPrimSec_temp)[, 
                  2], width = 1000), tinypeak.removal = tinyPYy)$Peaks
                if (length(Yy) > 1) {
                  lowPYy <- tinyPYy
                } else if (length(Yy) < 1) {
                  highPYy <- tinyPYy
                } else {
                  break
                }
            }
            if (abs(Xx - XxA) > 4 * scalingParam[1] || abs(Yy - YyA) > 4 * scalingParam[2]) {
                deviation2X <- c(deviation2X, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                  1] <= XxA + 1000), which(flowCore::exprs(f)[, 1] >= XxA - 1000)), 
                  intersect(which(flowCore::exprs(f)[, 2] <= YyA + 1000), which(flowCore::exprs(f)[, 
                    2] >= YyA - 1000))), 1]))
                deviation2Y <- c(deviation2Y, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                  1] <= XxA + 1000), which(flowCore::exprs(f)[, 1] >= XxA - 1000)), 
                  intersect(which(flowCore::exprs(f)[, 2] <= YyA + 1000), which(flowCore::exprs(f)[, 
                    2] >= YyA - 1000))), 2]))
                XcTert <- c(XcTert, XxA)
                YcTert <- c(YcTert, YyA)
            } else {
                deviation2X <- c(deviation2X, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                  1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), 
                  intersect(which(flowCore::exprs(f)[, 2] <= Yy + 1000), which(flowCore::exprs(f)[, 
                    2] >= Yy - 1000))), 1]))
                deviation2Y <- c(deviation2Y, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                  1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), 
                  intersect(which(flowCore::exprs(f)[, 2] <= Yy + 1000), which(flowCore::exprs(f)[, 
                    2] >= Yy - 1000))), 2]))
                XcTert <- c(XcTert, Xx)
                YcTert <- c(YcTert, Yy)
            }
        }
        
    } else {
        for (r1 in seq_len(4)) {
            switch(r1, {
                Xx <- mean(secondaryClusters$clusters[1, 1] + firstClusters$clusters[3, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[2, 1] + firstClusters$clusters[2, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[3, 1] + firstClusters$clusters[1, 
                  1] - emptyDroplets[1])
                Yy <- mean(secondaryClusters$clusters[1, 2] + firstClusters$clusters[3, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[2, 2] + firstClusters$clusters[2, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[3, 2] + firstClusters$clusters[1, 
                  2] - emptyDroplets[2])
            }, {
                Xx <- mean(secondaryClusters$clusters[1, 1] + firstClusters$clusters[4, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[4, 1] + firstClusters$clusters[2, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[5, 1] + firstClusters$clusters[1, 
                  1] - emptyDroplets[1])
                Yy <- mean(secondaryClusters$clusters[1, 2] + firstClusters$clusters[4, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[4, 2] + firstClusters$clusters[2, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[5, 2] + firstClusters$clusters[1, 
                  2] - emptyDroplets[2])
            }, {
                Xx <- mean(secondaryClusters$clusters[2, 1] + firstClusters$clusters[4, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[4, 1] + firstClusters$clusters[3, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[6, 1] + firstClusters$clusters[1, 
                  1] - emptyDroplets[1])
                Yy <- mean(secondaryClusters$clusters[2, 2] + firstClusters$clusters[4, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[4, 2] + firstClusters$clusters[3, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[6, 2] + firstClusters$clusters[1, 
                  2] - emptyDroplets[2])
            }, {
                Xx <- mean(secondaryClusters$clusters[3, 1] + firstClusters$clusters[4, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[5, 1] + firstClusters$clusters[3, 
                  1] - emptyDroplets[1], secondaryClusters$clusters[6, 1] + firstClusters$clusters[2, 
                  1] - emptyDroplets[1])
                Yy <- mean(secondaryClusters$clusters[3, 2] + firstClusters$clusters[4, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[5, 2] + firstClusters$clusters[3, 
                  2] - emptyDroplets[2], secondaryClusters$clusters[6, 2] + firstClusters$clusters[2, 
                  2] - emptyDroplets[2])
            })
            deviation2X <- c(deviation2X, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), intersect(which(flowCore::exprs(f)[, 
                2] <= Yy + 1000), which(flowCore::exprs(f)[, 2] >= Yy - 1000))), 
                1]))
            deviation2Y <- c(deviation2Y, stats::sd(flowCore::exprs(f)[intersect(intersect(which(flowCore::exprs(f)[, 
                1] <= Xx + 1000), which(flowCore::exprs(f)[, 1] >= Xx - 1000)), intersect(which(flowCore::exprs(f)[, 
                2] <= Yy + 1000), which(flowCore::exprs(f)[, 2] >= Yy - 1000))), 
                2]))
            XcTert <- c(XcTert, Xx)
            YcTert <- c(YcTert, Yy)
        }
    }
    return(list(clusters = cbind(XcTert, YcTert), deviation = cbind(deviation2X, 
        deviation2Y)))
}

# Find the quaternary clusters based on their density peaks found by flowDensity.
findQuaternaryClusterDensity <- function(f, f_onlyQuad, emptyDroplets, firstClusters, 
    secondaryClusters, tertiaryClusters) {
    NumberOfSinglePos <- nrow(firstClusters$clusters)
    if (nrow(f_onlyQuad) >= 4) {
        highP <- 1
        lowP <- 0
        repeat {
            tinyP <- (highP + lowP)/2
            XcQuad <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_onlyQuad)[, 
                1], width = 1000), tinypeak.removal = tinyP)$Peaks
            if (length(XcQuad) < 1) {
                highP <- tinyP
            } else if (length(XcQuad) > 1) {
                lowP <- tinyP
            } else {
                break
            }
        }
        highP <- 1
        lowP <- 0
        repeat {
            tinyP <- (highP + lowP)/2
            YcQuad <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_onlyQuad)[, 
                2], width = 1000), tinypeak.removal = tinyP)$Peaks
            if (length(YcQuad) < 1) {
                highP <- tinyP
            } else if (length(YcQuad) > 1) {
                lowP <- tinyP
            } else {
                break
            }
        }
    } else {
        if (is.null(tertiaryClusters)) {
            XcQuad <- mean(secondaryClusters$clusters[1, 1] + firstClusters$clusters[3, 
                1] - emptyDroplets[1], secondaryClusters$clusters[2, 1] + firstClusters$clusters[2, 
                1] - emptyDroplets[1], secondaryClusters$clusters[3, 1] + firstClusters$clusters[1, 
                1] - emptyDroplets[1])
            YcQuad <- mean(secondaryClusters$clusters[1, 2] + firstClusters$clusters[3, 
                2] - emptyDroplets[2], secondaryClusters$clusters[2, 2] + firstClusters$clusters[2, 
                2] - emptyDroplets[2], secondaryClusters$clusters[3, 2] + firstClusters$clusters[1, 
                2] - emptyDroplets[2])
        } else if (is.null(secondaryClusters)) {
            XcQuad <- firstClusters$clusters[1, 1] + firstClusters[2, 1] - emptyDroplets[1]
            YcQuad <- firstClusters$clusters[1, 2] + firstClusters[2, 2] - emptyDroplets[2]
        } else {
            XcQuad <- mean(tertiaryClusters$clusters[1, 1] + firstClusters$clusters[4, 
                1] - emptyDroplets[1], tertiaryClusters$clusters[2, 1] + firstClusters$clusters[3, 
                1] - emptyDroplets[1], tertiaryClusters$clusters[3, 1] + firstClusters$clusters[2, 
                1] - emptyDroplets[1], tertiaryClusters$clusters[4, 1] + firstClusters$clusters[1, 
                1] - emptyDroplets[1])
            YcQuad <- mean(tertiaryClusters$clusters[1, 2] + firstClusters$clusters[4, 
                2] - emptyDroplets[2], tertiaryClusters$clusters[2, 2] + firstClusters$clusters[3, 
                2] - emptyDroplets[2], tertiaryClusters$clusters[3, 2] + firstClusters$clusters[2, 
                2] - emptyDroplets[2], tertiaryClusters$clusters[4, 2] + firstClusters$clusters[1, 
                2] - emptyDroplets[2])
        }
    }
    return(list(clusters = cbind(XcQuad, YcQuad)))
}

# Find the primary clusters based on their density peaks found by flowDensity.
findPrimaryClusters <- function(data, clusterMeans, emptyDroplets, remove = 0, dimensions, 
    File, f, NumberOfSinglePos = 4, scalingParam, epsilon) {
    CutAbovePrimary <- mean(scalingParam) * 2
    NumOfClusters <- NumberOfSinglePos^2
    DataRemoved <- FinalResults <- NULL
    
    up1max <- deGate(f, c(1), percentile = 0.999, use.percentile = TRUE)
    up1min <- deGate(f, c(1), percentile = 0.001, use.percentile = TRUE)
    up2max <- deGate(f, c(2), percentile = 0.999, use.percentile = TRUE)
    up2min <- deGate(f, c(2), percentile = 0.001, use.percentile = TRUE)
    
    indices <- unique(c(which(flowCore::exprs(f)[, 1] >= 0.15 * (up1max - up1min) + 
        up1min), which(flowCore::exprs(f)[, 2] >= 0.15 * (up2max - up2min) + up2min)))
    
    # keep the non 15% bottom left corner (rn for removed neg)
    f_remNeg <- f
    flowCore::exprs(f_remNeg) <- flowCore::exprs(f_remNeg)[indices, ]  
    
    # keep the 15% bottom left corner
    f_onlyNeg <- f
    flowCore::exprs(f_onlyNeg) <- flowCore::exprs(f_onlyNeg)[-indices, ]  
    
    ClusterCentres <- vector()
    
    #---- find 1st gen clusters-------------------#
    
    # Calculate the threshold to remove events that are too positive. This makes it
    # easier to find single positives.
    f_findExtremes_temp <- f
    f_remNeg_temp <- f_remNeg
    threshold <- 0.1
    
    # find left and right primary
    dataTable <- table(data)
    clusterMeans2 <- clusterMeans[-c(emptyDroplets, remove), , drop = FALSE]
    if (!length(clusterMeans2)) 
        return(NULL)
    minimum <- match(min(clusterMeans2[, 1]), clusterMeans[, 1])
    a <- matrix(1, nrow = 3)
    realFirstCluster1 <- vector()
    collinear1 <- minimum
    for (i in seq_len(nrow(clusterMeans))) {
        if (i == minimum || i == emptyDroplets) 
            next
        test <- matrix(clusterMeans[c(emptyDroplets, minimum, i), ], ncol = 2)
        if (abs(det(cbind(a, test/100))) < (dimensions[1]/(NumberOfSinglePos * 20)) && 
            clusterMeans[i, 2] < (1 - NumberOfSinglePos^2/100) * dimensions[2]) 
            collinear1 <- c(collinear1, i)
    }
    selection <- which(clusterMeans[, 1] <= max(clusterMeans[collinear1, 1]))
    selection <- selection[!selection %in% c(emptyDroplets, remove)]
    selection <- selection[(dataTable[selection] > max(dataTable[selection])/4)]
    realFirstCluster1 <- selection[which.max(clusterMeans[selection, 2])]
    
    # collinear1 <- realFirstCluster1 for (i in seq_len(nrow(clusterMeans))) { if (i
    # == realFirstCluster1 || i == emptyDroplets) next test <-
    # matrix(clusterMeans[c(emptyDroplets,realFirstCluster1,i),], ncol=2) if
    # (abs(det(cbind(a, test/100))) < (dimensions[1]/(NumberOfSinglePos*20)))
    # collinear1 <- c(collinear1, i) }
    
    clusterMeans2 <- clusterMeans[-c(emptyDroplets, collinear1, remove), , drop = FALSE]
    if (!length(clusterMeans2)) 
        return(realFirstCluster1)
    minimum <- match(min(clusterMeans2[, 2]), clusterMeans[, 2])
    a <- matrix(1, nrow = 3)
    realFirstCluster2 <- vector()
    collinear2 <- minimum
    for (i in seq_len(nrow(clusterMeans))) {
        if (i == minimum || i == emptyDroplets) 
            next
        test <- matrix(clusterMeans[c(emptyDroplets, minimum, i), ], ncol = 2)
        if (abs(det(cbind(a, test/100))) < (dimensions[2]/(NumberOfSinglePos * 20)) && 
            clusterMeans[i, 1] < (1 - NumberOfSinglePos^2/100) * dimensions[1]) 
            collinear2 <- c(collinear2, i)
    }
    selection <- which(clusterMeans[, 2] <= max(clusterMeans[collinear2, 2]))
    selection <- selection[!selection %in% c(emptyDroplets, collinear1, remove)]
    selection <- selection[(dataTable[selection] > max(dataTable[selection])/4)]
    realFirstCluster2 <- selection[which.max(clusterMeans[selection, 1])]
    
    if (dataTable[realFirstCluster2] < dataTable[realFirstCluster1]/10) {
        remove <- c(remove, collinear2)
        return(findPrimaryClusters(data, clusterMeans, emptyDroplets, remove, 1.1 * 
            dimensions, File, f, NumberOfSinglePos, scalingParam, epsilon))
    }
    
    if (dataTable[realFirstCluster1] < dataTable[realFirstCluster2]/10) {
        remove <- c(remove, collinear1)
        return(findPrimaryClusters(data, clusterMeans, emptyDroplets, remove, 1.1 * 
            dimensions, File, f, NumberOfSinglePos, scalingParam, epsilon))
    }
    
    x_leftPrim <- clusterMeans[realFirstCluster1, 1]
    y_leftPrim <- clusterMeans[realFirstCluster1, 2]
    x_rightPrim <- clusterMeans[realFirstCluster2, 1]
    y_rightPrim <- clusterMeans[realFirstCluster2, 2]
    
    mSlope <- (y_leftPrim - y_rightPrim)/(x_leftPrim - x_rightPrim)
    theta <- abs(atan(mSlope))
    rotate <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
    
    Rot_xy_leftPrim <- rotate %*% c(x_leftPrim, y_leftPrim)  # coordinates of rotated left  primary cluster
    Rot_xy_rightPrim <- rotate %*% c(x_rightPrim, y_rightPrim)  # coordinates of rotated right primary cluster
    
    flowCore::exprs(f_findExtremes_temp)[, c(1, 2)] <- t(rotate %*% t(flowCore::exprs(f_findExtremes_temp)[, 
        c(1, 2)]))
    flowCore::exprs(f_remNeg_temp)[, c(1, 2)] <- t(rotate %*% t(flowCore::exprs(f_remNeg_temp)[, 
        c(1, 2)]))
    upSlantmax <- deGate(f_findExtremes_temp, c(2), percentile = 0.999, use.percentile = TRUE)
    upSlantmin <- deGate(f_findExtremes_temp, c(2), percentile = 0.001, use.percentile = TRUE)
    # Chop off the very positive events
    ScaleChop <- (upSlantmax - upSlantmin)/max(File)
    indices <- which(flowCore::exprs(f_remNeg_temp)[, 2] <= (max(Rot_xy_leftPrim[2], 
        Rot_xy_rightPrim[2]) + CutAbovePrimary * ScaleChop))
    flowCore::exprs(f_remNeg_temp) <- flowCore::exprs(f_remNeg_temp)[indices, ]
    
    # find thresholds to divide the primary clusters
    highP <- 1
    lowP <- 0
    repeat {
        # newton iteration
        newP <- (highP + lowP)/2
        if (newP < 2 * epsilon) 
            break
        Cuts <- deGate(f_remNeg_temp, 1, all.cuts = TRUE, tinypeak.removal = newP, 
            adjust.dens = 0.1)
        if (length(Cuts) < (NumberOfSinglePos - 1)) {
            highP <- newP
        } else if (length(Cuts) > (NumberOfSinglePos - 1)) {
            lowP <- newP
        } else {
            break
        }
    }
    Xc <- Yc <- NULL
    Cuts <- c(min(flowCore::exprs(f_remNeg_temp)[, 1]), Cuts, max(flowCore::exprs(f_remNeg_temp)[, 
        1]))
    
    # find the coordinates of the primary clusters
    for (r1 in seq_len(length(Cuts) - 1)) {
        indices <- intersect(which(flowCore::exprs(f_remNeg_temp)[, 1] >= Cuts[r1]), 
            which(flowCore::exprs(f_remNeg_temp)[, 1] < Cuts[r1 + 1]))
        f_Cuts_temp <- f_remNeg_temp
        flowCore::exprs(f_Cuts_temp) <- flowCore::exprs(f_Cuts_temp)[indices, ]
        flowCore::exprs(f_Cuts_temp)[, c(1, 2)] <- t(t(rotate) %*% t(flowCore::exprs(f_Cuts_temp)[, 
            c(1, 2)]))
        Xx <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
            1], width = 1000), tinypeak.removal = newP)$Peaks
        Yy <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
            2], width = 1000), tinypeak.removal = newP)$Peaks
        if (length(Xx) > 1 || length(Yy) > 1) {
            highP <- 1
            lowP <- 0
            repeat {
                # newton iteration
                newP <- (highP + lowP)/2
                if (newP < epsilon) 
                  break
                Xx <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
                  1], width = 1000), tinypeak.removal = newP)$Peaks
                Yy <- flowDensity::getPeaks(stats::density(flowCore::exprs(f_Cuts_temp)[, 
                  2], width = 1000), tinypeak.removal = newP)$Peaks
                if (length(Xx) > 1 || length(Yy) > 1) {
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
    for (i in seq_len(nrow(densPrimaries))) {
        temp <- lapply(seq_len(nrow(clusterMeans)), function(x) return(clusterMeans[x, 
            ] - densPrimaries[i, ]))
        rowSums <- vapply(seq_along(temp), function(x) return(sum(abs(unlist(temp[x])))), 
            double(length = 1))
        distMatrix <- rbind(distMatrix, rowSums)
    }
    firstClusters <- solve_LSAP(distMatrix)
}

# Find the primary clusters based on their position.
findPrimaryClusters_old <- function(data, clusterMeans, emptyDroplets, remove = 0, 
    dimensions) {
    dataTable <- table(data)
    clusterMeans2 <- clusterMeans[-c(emptyDroplets, remove), , drop = FALSE]
    if (!length(clusterMeans2)) 
        return(NULL)
    minimum <- match(min(clusterMeans2[, 1]), clusterMeans[, 1])
    a <- matrix(1, nrow = 3)
    realFirstCluster1 <- vector()
    collinear1 <- minimum
    for (i in seq_len(nrow(clusterMeans))) {
        if (i == minimum || i == emptyDroplets) 
            next
        test <- matrix(clusterMeans[c(emptyDroplets, minimum, i), ], ncol = 2)
        if (abs(det(cbind(a, test/100))) < (dimensions[1]/100)) 
            collinear1 <- c(collinear1, i)
    }
    selection <- which(clusterMeans[, 1] <= max(clusterMeans[collinear1, 1]))
    selection <- selection[!selection %in% c(emptyDroplets, remove)]
    selection <- selection[(dataTable[selection] > max(dataTable[selection])/4)]
    realFirstCluster1 <- selection[which.min(clusterMeans[selection, 1])]
    
    clusterMeans2 <- clusterMeans[-c(emptyDroplets, collinear1, remove), , drop = FALSE]
    if (!length(clusterMeans2)) 
        return(realFirstCluster1)
    minimum <- match(min(clusterMeans2[, 2]), clusterMeans[, 2])
    a <- matrix(1, nrow = 3)
    realFirstCluster2 <- vector()
    collinear2 <- minimum
    for (i in seq_len(nrow(clusterMeans))) {
        if (i == minimum || i == emptyDroplets) 
            next
        test <- matrix(clusterMeans[c(emptyDroplets, minimum, i), ], ncol = 2)
        if (abs(det(cbind(a, test/100))) < (dimensions[2]/100)) 
            collinear2 <- c(collinear2, i)
    }
    selection <- which(clusterMeans[, 2] <= max(clusterMeans[collinear2, 2]))
    selection <- selection[!selection %in% c(emptyDroplets, collinear1, remove)]
    selection <- selection[(dataTable[selection] > max(dataTable[selection])/4)]
    realFirstCluster2 <- selection[which.min(clusterMeans[selection, 2])]
    
    if (dataTable[realFirstCluster2] < dataTable[realFirstCluster1]/5) {
        remove <- c(remove, collinear2)
        return(findPrimaryClusters(data, clusterMeans, emptyDroplets, remove, dimensions))
    }
    
    if (dataTable[realFirstCluster1] < dataTable[realFirstCluster2]/5) {
        remove <- c(remove, collinear1)
        return(findPrimaryClusters(data, clusterMeans, emptyDroplets, remove, dimensions))
    }
    
    clusterMeans2 <- clusterMeans[-c(emptyDroplets, collinear1, collinear2), , drop = FALSE]
    moreFirstClusters <- vector()
    counter <- clusterMeans[realFirstCluster2, 1]
    clusterMeans2 <- clusterMeans2[clusterMeans2[, 1] < counter, , drop = FALSE]
    counter <- clusterMeans[realFirstCluster1, 2]
    temp <- clusterMeans2[clusterMeans2[, 2] < counter, , drop = FALSE]
    repeat {
        if (!length(temp)) {
            break
        } else {
            minimum <- which.min(temp[, 1])
            cluster <- match(temp[minimum], clusterMeans)
            if (dataTable[cluster] < (mean(dataTable[c(realFirstCluster1, realFirstCluster2)]))/5) {
                temp <- temp[-minimum, , drop = FALSE]
                next
            } else {
                moreFirstClusters <- c(moreFirstClusters, cluster)
                counter <- clusterMeans[cluster, 2]
                temp <- clusterMeans2[clusterMeans2[, 2] < counter, , drop = FALSE]
            }
        }
    }
    
    # newTable <- dataTable[-c(emptyDroplets, realFirstCluster1, realFirstCluster2)]
    # moreFirstClusters <- names(newTable[newTable >
    # (mean(dataTable[c(realFirstCluster1, realFirstCluster2)]))/2])
    realFirstClusters <- as.integer(c(realFirstCluster1, moreFirstClusters, realFirstCluster2))
}

# Find the secondary clusters based on the positions of the primary clusters.
findSecondaryClusters <- function(firstClusters, clusterMeans, emptyDroplets, remove, 
    dimensions, counts) {
    threshold <- dimensions/5
    clusterMeans2 <- clusterMeans[-c(emptyDroplets, firstClusters, remove), ]
    if (!length(clusterMeans2)) 
        return(list(clusters = 0, correctionFactor = 0))
    secondClusters <- vector()
    correction <- list()
    correction[[length(firstClusters) + 1]] <- 0
    
    if (nrow(clusterMeans2) >= sum(seq_len(length(firstClusters) - 1))) {
        distMatrix <- vector()
        for (a in firstClusters) {
            for (b in firstClusters) {
                if (match(a, firstClusters) >= match(b, firstClusters)) 
                  next
                temp <- lapply(seq_len(nrow(clusterMeans2)), function(x) return(clusterMeans2[x, 
                  ] - (clusterMeans[a, ] + clusterMeans[b, ] - clusterMeans[emptyDroplets, 
                  ])))
                # temp2 <- lapply(seq_along(temp), function(x) return(c(temp[[x]][1]*3,
                # temp[[x]][2]/3)))
                rowSums <- vapply(seq_along(temp), function(x) return(sum(abs(unlist(temp[x])))), 
                  double(length = 1))
                distMatrix <- rbind(distMatrix, rowSums)
            }
        }
        tmp <- solve_LSAP(distMatrix)
        secondClusters <- match(clusterMeans2[tmp], clusterMeans)
        index <- 0
        for (a in firstClusters) {
            for (b in firstClusters) {
                if (match(a, firstClusters) >= match(b, firstClusters)) 
                  next
                index <- index + 1
                distance <- clusterMeans[secondClusters[index], ] - (clusterMeans[a, 
                  ] + clusterMeans[b, ] - clusterMeans[emptyDroplets, ])
                if (sum(abs(distance)) > threshold) {
                  secondClusters[index] = 0
                } else {
                  correction[[match(a, firstClusters)]] <- c(correction[[match(a, 
                    firstClusters)]], list(distance))
                  correction[[match(b, firstClusters)]] <- c(correction[[match(b, 
                    firstClusters)]], list(distance))
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
                  if (!length(clusterMeans2)) {
                    cluster <- 0
                    secondClusters <- c(secondClusters, cluster)
                    next
                  }
                  temp <- lapply(seq_len(nrow(clusterMeans2)), function(x) return(clusterMeans2[x, 
                    ] - (clusterMeans[a, ] + clusterMeans[b, ] - clusterMeans[emptyDroplets, 
                    ])))
                  rowSums <- vapply(seq_along(temp), function(x) return(sum(abs(unlist(temp[x])))), 
                    double(length = 1))
                  if (min(rowSums) < threshold) {
                    minimum <- which.min(rowSums)
                    cluster <- match(clusterMeans2[minimum], clusterMeans)
                    clusterMeans2 <- clusterMeans2[-minimum, , drop = FALSE]
                    correction[[match(a, firstClusters)]] <- c(correction[[match(a, 
                      firstClusters)]], temp[minimum])
                    correction[[match(b, firstClusters)]] <- c(correction[[match(b, 
                      firstClusters)]], temp[minimum])
                    corTemp <- temp[[minimum]]
                  } else {
                    temp <- lapply(seq_len(nrow(clusterMeans2)), function(x) return(clusterMeans2[x, 
                      ] - (clusterMeans[a, ] + clusterMeans[b, ] - clusterMeans[emptyDroplets, 
                      ] + corTemp)))
                    rowSums <- vapply(seq_along(temp), function(x) return(sum(abs(unlist(temp[x])))), 
                      double(length = 1))
                    if (min(rowSums) < threshold) {
                      minimum <- which.min(rowSums)
                      cluster <- match(clusterMeans2[minimum], clusterMeans)
                      clusterMeans2 <- clusterMeans2[-minimum, , drop = FALSE]
                      correction[[match(a, firstClusters)]] <- c(correction[[match(a, 
                        firstClusters)]], temp[minimum])
                      correction[[match(b, firstClusters)]] <- c(correction[[match(b, 
                        firstClusters)]], temp[minimum])
                      corTemp <- temp[[minimum]]
                    } else {
                      cluster <- 0
                      # correctionFactor <- c(correctionFactor, 0)
                    }
                  }
                  
                  secondClusters <- c(secondClusters, cluster)
                }
            }
        }
    }
    correctionFactor <- list()
    for (i in seq_len(length(correction) - 1)) {
        if (!length(correction[[i]])) {
            correctionFactor[[i]] <- 0
        } else {
            temp <- t(do.call(cbind, correction[[i]]))
            correctionFactor[[i]] <- colMeans(temp)
        }
    }
    
    if (max(secondClusters) == 0) 
        return(list(clusters = secondClusters, correctionFactor = correctionFactor))
    
    q <- stats::quantile(counts[secondClusters])
    
    if ((stats::IQR(counts[secondClusters]) > q[3]) && (q[3] > 5) && (q[1] < 5)) {
        remove <- c(remove, as.integer(names(which.min(counts[secondClusters]))))
        return(findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, 
            remove, dimensions, counts))
    }
    if (min(counts[secondClusters]) < q[2]/5) {
        remove <- c(remove, as.integer(names(which.min(counts[secondClusters]))))
        return(findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, 
            remove, dimensions, counts))
    }
    if (max(counts[secondClusters]) > q[4] * 5) {
        remove <- c(remove, as.integer(names(which.max(counts[secondClusters]))))
        return(findSecondaryClusters(firstClusters, clusterMeans, emptyDroplets, 
            remove, dimensions, counts))
    }
    
    ## the correction factor is the avarage distance between the estimated positions
    ## of secondClusters and the real positions
    list(clusters = secondClusters, correctionFactor = correctionFactor)
}

# Find the tertiary clusters based on the positions of the primary and secondary
# clusters.
findTertiaryClusters <- function(emptyDroplets, firstClusters, secondClusters, remove, 
    clusterMeans, correctionFactor, dimensions, counts) {
    threshold <- dimensions/6
    clusterMeans2 <- clusterMeans[-c(emptyDroplets, firstClusters, secondClusters, 
        remove), , drop = FALSE]
    if (!length(clusterMeans2)) 
        return(list(clusters = rep(0, 4)))
    thirdClusters1 <- vector()
    thirdClusters2 <- vector()
    
    if (nrow(clusterMeans2) >= 4) {
        for (i in seq_len(4)) {
            if (i < 3) {
                cluster <- utils::tail(secondClusters, n = 1)
                if (cluster == 0 || length(clusterMeans2) == 0) {
                  thirdClusters1 <- rbind(thirdClusters1, rep(dimensions, nrow(clusterMeans2)))
                  next
                }
                temp <- lapply(seq_len(nrow(clusterMeans2)), function(x) return(clusterMeans2[x, 
                  ] - (clusterMeans[cluster, ] + clusterMeans[firstClusters[i], ] - 
                  clusterMeans[emptyDroplets, ] + correctionFactor[[i]])))
                rowSums <- vapply(seq_along(temp), function(x) return(sum(abs(unlist(temp[x])))), 
                  double(length = 1))
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
                temp <- lapply(seq_len(nrow(clusterMeans2)), function(x) return(clusterMeans2[x, 
                  ] - (clusterMeans[cluster, ] + clusterMeans[firstClusters[i], ] - 
                  clusterMeans[emptyDroplets, ] + correctionFactor[[i]])))
                rowSums <- vapply(seq_along(temp), function(x) return(sum(abs(unlist(temp[x])))), 
                  double(length = 1))
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
        for (i in seq_along(tmp)) {
            if (distMatrix[i, tmp[i]] > threshold) {
                thirdClusters[i] = 0
            }
        }
        
    } else {
        for (i in seq_len(4)) {
            if (i < 3) {
                cluster <- utils::tail(secondClusters, n = 1)
                if (cluster == 0 || length(clusterMeans2) == 0) {
                  thirdClusters1 <- c(thirdClusters1, 0)
                  next
                }
                temp <- lapply(seq_len(nrow(clusterMeans2)), function(x) return(clusterMeans2[x, 
                  ] - (clusterMeans[cluster, ] + clusterMeans[firstClusters[i], ] - 
                  2 * clusterMeans[emptyDroplets, ] + correctionFactor[[i]])))
                rowSums <- vapply(seq_along(temp), function(x) return(sum(abs(unlist(temp[x])))), 
                  double(length = 1))
                if (min(rowSums) < threshold) {
                  minimum <- which.min(rowSums)
                  cluster <- match(clusterMeans2[minimum], clusterMeans)
                  clusterMeans2 <- clusterMeans2[-minimum, , drop = FALSE]
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
                temp <- lapply(seq_len(nrow(clusterMeans2)), function(x) return(clusterMeans2[x, 
                  ] - (clusterMeans[cluster, ] + clusterMeans[firstClusters[i], ] - 
                  2 * clusterMeans[emptyDroplets, ] + correctionFactor[[i]])))
                rowSums <- vapply(seq_along(temp), function(x) return(sum(abs(unlist(temp[x])))), 
                  double(length = 1))
                if (min(rowSums) < threshold) {
                  minimum <- which.min(rowSums)
                  cluster <- match(clusterMeans2[minimum], clusterMeans)
                  clusterMeans2 <- clusterMeans2[-minimum, , drop = FALSE]
                } else {
                  cluster <- 0
                }
                thirdClusters2 <- c(thirdClusters2, cluster)
            }
        }
        
        thirdClusters <- c(thirdClusters2, thirdClusters1)
        
    }
    
    if (max(thirdClusters) == 0) 
        return(list(clusters = thirdClusters))
    
    q <- stats::quantile(counts[thirdClusters])
    
    if ((stats::IQR(counts[thirdClusters]) > q[3]) && (q[3] > 5) && (q[1] < 5)) {
        remove <- c(remove, as.integer(names(which.min(counts[thirdClusters]))))
        return(findTertiaryClusters(emptyDroplets, firstClusters, secondClusters, 
            remove, clusterMeans, correctionFactor, dimensions, counts))
    }
    if (min(counts[thirdClusters]) < q[2]/5) {
        remove <- c(remove, as.integer(names(which.min(counts[thirdClusters]))))
        return(findTertiaryClusters(emptyDroplets, firstClusters, secondClusters, 
            remove, clusterMeans, correctionFactor, dimensions, counts))
    }
    if (max(counts[thirdClusters]) > q[4] * 5) {
        remove <- c(remove, as.integer(names(which.max(counts[thirdClusters]))))
        return(findTertiaryClusters(emptyDroplets, firstClusters, secondClusters, 
            remove, clusterMeans, correctionFactor, dimensions, counts))
    }
    
    list(clusters = thirdClusters)
}

# Find the quarternary cluster based on its position.
findQuaternaryCluster <- function(clusterMeans, emptyDroplets, remove, firstClusters, 
    secondClusters = 0, thirdClusters = 0) {
    ## TODO really find the cluster?
    clusterMeans2 <- clusterMeans[-c(emptyDroplets, firstClusters, secondClusters, 
        thirdClusters), , drop = FALSE]
    if (!length(clusterMeans2)) 
        return(0)
    rowSums <- vapply(seq_len(nrow(clusterMeans)), function(x) return(sum(clusterMeans[x, 
        ])), double(length = 1))
    rowSums2 <- vapply(seq_len(nrow(clusterMeans2)), function(x) return(sum(clusterMeans2[x, 
        ])), double(length = 1))
    if (match(max(rowSums), rowSums) == match(max(rowSums2), rowSums)) {
        return(which.max(rowSums))
    } else {
        return(0)
    }
}

# Merge clusters that a close to each other, based on the parameter p.
mergeClusters <- function(result, clusterMeans, finalResult, remove, p = 12) {
    # threshold <- dimensions/(length(finalResult)*2)
    clusterMeans2 <- clusterMeans[-c(finalResult, remove), , drop = FALSE]
    if (nrow(clusterMeans2) > 0) {
        for (i in seq_len(nrow(clusterMeans2))) {
            for (j in seq_len(nrow(clusterMeans))) {
                # threshold <- clusterMeans[j,]/as.integer(sqrt(length(finalResult))*1.5)
                threshold <- sqrt(clusterMeans[j, ]) * p
                if (abs(clusterMeans2[i, 1] - clusterMeans[j, 1]) < threshold[1] && 
                  abs(clusterMeans2[i, 2] - clusterMeans[j, 2]) < threshold[2]) {
                  result[result == match(clusterMeans2[i], clusterMeans)] = j
                }
            }
        }
    }
    return(result)
}

# Adjust the cluster centres based on the local density function
adjustClusterMeans <- function(data, clusterMeans, result, clusters) {
    for (i in clusters) {
        if (i == 0) 
            next
        f <- flowCore::flowFrame(exprs = as.matrix(data[result == i, ]))
        Xc <- flowDensity::getPeaks(stats::density(flowCore::exprs(f)[, 1], width = 1000), 
            tinypeak.removal = 0.2)$Peaks[1]
        Yc <- flowDensity::getPeaks(stats::density(flowCore::exprs(f)[, 2], width = 1000), 
            tinypeak.removal = 0.2)$Peaks[1]
        if (!is.null(Xc) && !is.null(Yc)) {
            clusterMeans[i, ] <- c(Xc, Yc)
        }
    }
    return(clusterMeans)
}

# Find the rain and assign it based on the distance to vector lines connecting
# the cluster centres.
assignRain <- function(clusterMeans, data, result, emptyDroplets, firstClusters, 
    secondClusters, thirdClusters, fourthCluster, scalingParam, similarityParam = 0.95, 
    distanceParam = 0.2) {
    remove <- vector()
    sdeviations <- list()
    rownames(data) <- seq_len(nrow(data))
    scalingHelper <- mean(scalingParam/4)
    
    # Calculate standard deviations:
    for (cl in firstClusters) {
        if (cl == 0) 
            next
        sdeviation <- 2 * c(stats::sd(data[which(result == cl), 1]), stats::sd(data[which(result == 
            cl), 2]))
        if (anyNA(sdeviation)) 
            sdeviation <- scalingParam/2
        if (any(sdeviation > scalingParam)) 
            sdeviation <- scalingParam
        sdeviations[[cl]] <- sdeviation
    }
    for (cl in secondClusters) {
        if (cl == 0) 
            next
        sdeviation <- 2 * c(stats::sd(data[which(result == cl), 1]), stats::sd(data[which(result == 
            cl), 2]))
        if (anyNA(sdeviation)) 
            sdeviation <- scalingParam/2
        if (any(sdeviation > scalingParam)) 
            sdeviation <- scalingParam
        sdeviations[[cl]] <- sdeviation
    }
    for (cl in thirdClusters) {
        if (cl == 0) 
            next
        sdeviation <- 2 * c(stats::sd(data[which(result == cl), 1]), stats::sd(data[which(result == 
            cl), 2]))
        if (anyNA(sdeviation)) 
            sdeviation <- scalingParam/2
        if (any(sdeviation > scalingParam)) 
            sdeviation <- scalingParam
        sdeviations[[cl]] <- sdeviation
    }
    
    ## Empties to primary clusters:
    newData <- subset(data, !result %in% c(secondClusters, thirdClusters, fourthCluster))
    thisCluster <- which(rownames(newData) %in% which(result == 1))
    S <- stats::var(newData[thisCluster, , drop = FALSE])
    mahaDist <- stats::mahalanobis(newData, clusterMeans[1, ], S)
    posEv <- which(mahaDist < mean(mahaDist[thisCluster]))
    temp1 <- which(newData[, 1] <= clusterMeans[1, 1])
    temp2 <- which(newData[, 2] <= clusterMeans[1, 2])
    posEv <- c(posEv, intersect(temp1, temp2))
    result[as.numeric(rownames(newData[posEv, ]))] <- 1
    newData <- newData[-posEv, , drop = FALSE]
    # Go through each cluster:
    for (cl in firstClusters) {
        if (cl == 0) 
            next
        thisCluster <- which(rownames(newData) %in% which(result == cl))
        S <- stats::var(newData[thisCluster, , drop = FALSE])
        mahaDist <- stats::mahalanobis(newData, clusterMeans[cl, ], S)
        posEv <- which(mahaDist < mean(mahaDist[thisCluster]))
        result[as.numeric(rownames(newData[posEv, ]))] <- cl
        newData <- newData[-posEv, , drop = FALSE]
        sdC <- sdeviations[[cl]]
        newData <- subset(newData, !(newData[, 1] > clusterMeans[cl, 1] - sdC[1] & 
            newData[, 2] > clusterMeans[cl, 2] - sdC[2]))
        # newData <- newData[-intersect(as.numeric(rownames(newData)),
        # which(result==cl)),]
    }
    if (nrow(newData) > 0) {
        m <- matrix(nrow = nrow(newData), ncol = length(firstClusters))
        for (cl in seq_along(firstClusters)) {
            if (firstClusters[cl] == 0) 
                next
            for (row in seq_len(nrow(newData))) {
                m[row, cl] <- distToLineSegment(as.numeric(newData[row, ]), as.numeric(clusterMeans[emptyDroplets, 
                  ]), as.numeric(clusterMeans[firstClusters[cl], ]))
            }
        }
        minimalDistance <- apply(m, 1, function(x) which(x == min(x, na.rm = TRUE))[1])
        for (r in seq_along(minimalDistance)) {
            distPoint <- euc.dist(newData[r, ], clusterMeans[emptyDroplets, ])
            distTotal <- distPoint + euc.dist(newData[r, ], clusterMeans[firstClusters[minimalDistance[r]], 
                ])
            if (distPoint <= distanceParam * distTotal) {
                result[as.numeric(rownames(newData)[r])] <- emptyDroplets
            } else {
                if (ncol(m) > 1) {
                  minAnd2ndMin <- sort(m[r, ], index.return = TRUE)
                  if (minAnd2ndMin$x[1]/minAnd2ndMin$x[2] >= similarityParam || (minAnd2ndMin$x[1] <= 
                    scalingHelper && minAnd2ndMin$x[2] <= scalingHelper)) {
                    remove <- c(remove, as.numeric(rownames(newData)[r]))
                  }
                }
                result[as.numeric(rownames(newData)[r])] <- firstClusters[minimalDistance[r]]
            }
        }
    }
    
    ## primary to secondary clusters (3-plex and 4-plex only):
    if (!is.null(secondClusters) && sum(secondClusters) > 0) {
        newData <- subset(data, !result %in% c(emptyDroplets, thirdClusters, fourthCluster))
        for (cl in firstClusters) {
            if (cl == 0) 
                next
            sdC <- sdeviations[[cl]]
            newData <- subset(newData, !(newData[, 1] < (clusterMeans[cl, 1] + sdC[1]) & 
                newData[, 2] < (clusterMeans[cl, 2] + sdC[2])))
            # # if (clusterMeans[cl,1] < clusterMeans[cl,2]) { # newData <- subset(newData,
            # !(newData[,1] < (clusterMeans[cl,1]+3*sdC[1]) & newData[,2] <
            # (clusterMeans[cl,2]-2*sdC[2]))) # } else { # newData <- subset(newData,
            # !(newData[,1] < (clusterMeans[cl,1]-2*sdC[1]) & newData[,2] <
            # (clusterMeans[cl,2]+3*sdC[2]))) # }
        }
        for (cl in secondClusters) {
            if (cl == 0) 
                next
            thisCluster <- which(rownames(newData) %in% which(result == cl))
            S <- stats::var(newData[thisCluster, , drop = FALSE])
            mahaDist <- stats::mahalanobis(newData, clusterMeans[cl, ], S)
            posEv <- which(mahaDist < mean(mahaDist[thisCluster]))
            result[as.numeric(rownames(newData[posEv, ]))] <- cl
            newData <- newData[-posEv, , drop = FALSE]
            sdC <- sdeviations[[cl]]
            newData <- subset(newData, !(newData[, 1] > clusterMeans[cl, 1] - sdC[1] & 
                newData[, 2] > clusterMeans[cl, 2] - sdC[2]))
        }
        if (nrow(newData) > 0) {
            i <- 1
            j <- 2
            m <- matrix(nrow = nrow(newData), ncol = 2 * length(secondClusters))
            for (cl in seq_along(secondClusters)) {
                if (secondClusters[cl] == 0) 
                  next
                for (row in seq_len(nrow(newData))) {
                  m[row, cl] <- distToLineSegment(as.numeric(newData[row, ]), as.numeric(clusterMeans[firstClusters[i], 
                    ]), as.numeric(clusterMeans[secondClusters[cl], ]))
                  m[row, cl + length(secondClusters)] <- distToLineSegment(as.numeric(newData[row, 
                    ]), as.numeric(clusterMeans[firstClusters[j], ]), as.numeric(clusterMeans[secondClusters[cl], 
                    ]))
                }
                j <- j + 1
                if (j > length(firstClusters)) {
                  i <- i + 1
                  j <- i + 1
                }
            }
            minimalDistance <- apply(m, 1, function(x) which(x == min(x, na.rm = TRUE))[1])
            
            m2 <- matrix(nrow = nrow(newData), ncol = length(firstClusters) + length(secondClusters))
            for (row in seq_len(nrow(newData))) {
                for (cl1 in seq_along(firstClusters)) {
                  if (firstClusters[cl1] == 0) 
                    next
                  m2[row, cl1] <- distToRect(as.numeric(clusterMeans[firstClusters[cl1], 
                    ]) - sdeviations[[firstClusters[cl1]]], as.numeric(clusterMeans[firstClusters[cl1], 
                    ]) + sdeviations[[firstClusters[cl1]]], as.numeric(newData[row, 
                    ]))
                }
                for (cl2 in seq_along(secondClusters)) {
                  if (secondClusters[cl2] == 0) 
                    next
                  col <- cl2 + length(firstClusters)
                  m2[row, col] <- distToRect(as.numeric(clusterMeans[secondClusters[cl2], 
                    ]) - sdeviations[[secondClusters[cl2]]], as.numeric(clusterMeans[secondClusters[cl2], 
                    ]) + sdeviations[[secondClusters[cl2]]], as.numeric(newData[row, 
                    ]))
                }
            }
            minimalClDistance <- apply(m2, 1, function(x) which(x == min(x, na.rm = TRUE))[1])
            
            
            helper <- c(1, 1, 1, 2, 2, 3, 2, 3, 4, 3, 4, 4)
            if (length(secondClusters) == 3) 
                helper <- c(1, 1, 2, 2, 3, 3)
            helper2 <- rep(seq_along(secondClusters), 2)
            helper5 <- c(firstClusters, secondClusters)
            
            for (r in seq_along(minimalDistance)) {
                if (m2[r, minimalClDistance[r]] <= m[r, minimalDistance[r]]) {
                  result[as.numeric(rownames(newData)[r])] <- helper5[minimalClDistance[r]]
                  next
                }
                
                distPoint <- euc.dist(newData[r, ], clusterMeans[firstClusters[helper[minimalDistance[r]]], 
                  ] + sdeviations[[firstClusters[helper[minimalDistance[r]]]]])
                distTotal <- distPoint + euc.dist(newData[r, ], clusterMeans[secondClusters[helper2[minimalDistance[r]]], 
                  ])
                
                if (distPoint <= distanceParam * distTotal) {
                  result[as.numeric(rownames(newData)[r])] <- firstClusters[helper[minimalDistance[r]]]
                } else {
                  if (ncol(m) > 1) {
                    minAnd2ndMin <- sort(m[r, ], index.return = TRUE)
                    if ((minAnd2ndMin$x[1]/minAnd2ndMin$x[2] >= similarityParam && 
                      helper2[minAnd2ndMin$ix[1]] != helper2[minAnd2ndMin$ix[2]]) || 
                      ((minAnd2ndMin$x[1] <= scalingHelper && minAnd2ndMin$x[2] <= 
                        scalingHelper) && helper2[minAnd2ndMin$ix[1]] != helper2[minAnd2ndMin$ix[2]])) {
                      remove <- c(remove, as.numeric(rownames(newData)[r]))
                    }
                  }
                  result[as.numeric(rownames(newData)[r])] <- secondClusters[helper2[minimalDistance[r]]]
                }
            }
        }
    }
    
    ## secondary to tertiary clusters (4-plex only):
    if (!is.null(thirdClusters) && sum(thirdClusters) > 0) {
        newData <- subset(data, !result %in% c(emptyDroplets, firstClusters, fourthCluster))
        for (cl in secondClusters) {
            if (cl == 0) 
                next
            sdC <- sdeviations[[cl]]
            newData <- subset(newData, !(newData[, 1] < (clusterMeans[cl, 1] + sdC[1]) & 
                newData[, 2] < (clusterMeans[cl, 2] + sdC[2])))
            # # if (clusterMeans[cl,1] < clusterMeans[cl,2]) { # newData <- subset(newData,
            # !(newData[,1] < (clusterMeans[cl,1]+2*sdC[1]) & newData[,2] <
            # (clusterMeans[cl,2]-sdC[2]))) # } else { # newData <- subset(newData,
            # !(newData[,1] < (clusterMeans[cl,1]-sdC[1]) & newData[,2] <
            # (clusterMeans[cl,2]+2*sdC[2]))) # }
        }
        for (cl in thirdClusters) {
            if (cl == 0) 
                next
            thisCluster <- which(rownames(newData) %in% which(result == cl))
            S <- stats::var(newData[thisCluster, , drop = FALSE])
            mahaDist <- stats::mahalanobis(newData, clusterMeans[cl, ], S)
            posEv <- which(mahaDist < mean(mahaDist[thisCluster]))
            result[as.numeric(rownames(newData[posEv, ]))] <- cl
            newData <- newData[-posEv, , drop = FALSE]
            sdC <- sdeviations[[cl]]
            newData <- subset(newData, !(newData[, 1] > clusterMeans[cl, 1] - sdC[1] & 
                newData[, 2] > clusterMeans[cl, 2] - sdC[2]))
        }
        if (nrow(newData) > 0) {
            m <- matrix(nrow = nrow(newData), ncol = 3 * length(thirdClusters))
            helper3 <- c(1, 1, 2, 4, 2, 3, 3, 5, 4, 5, 6, 6)
            helper4 <- rep(seq_along(thirdClusters), 3)
            for (cl in seq_along(thirdClusters)) {
                if (thirdClusters[cl] == 0) 
                  next
                for (row in seq_len(nrow(newData))) {
                  m[row, cl] <- 
                    distToLineSegment(as.numeric(newData[row, ]), 
                                      as.numeric(clusterMeans[secondClusters[helper3[cl]],]), 
                                      as.numeric(clusterMeans[thirdClusters[cl], ]))
                  m[row, cl + length(thirdClusters)] <- 
                    distToLineSegment(as.numeric(newData[row,]), 
                                      as.numeric(clusterMeans[secondClusters[helper3[cl + length(thirdClusters)]],]), 
                                      as.numeric(clusterMeans[thirdClusters[cl], ]))
                  m[row, cl + 2 * length(thirdClusters)] <- 
                    distToLineSegment(as.numeric(newData[row,]), 
                                      as.numeric(clusterMeans[secondClusters[helper3[cl + 2 * length(thirdClusters)]],]), 
                                      as.numeric(clusterMeans[thirdClusters[cl], ]))
                }
            }
            minimalDistance <- apply(m, 1, function(x) which(x == min(x, na.rm = TRUE))[1])
            m2 <- matrix(nrow = nrow(newData), ncol = length(secondClusters) + length(thirdClusters))
            for (row in seq_len(nrow(newData))) {
                for (cl1 in seq_along(secondClusters)) {
                  if (secondClusters[cl1] == 0) 
                    next
                  m2[row, cl1] <- distToRect(as.numeric(clusterMeans[secondClusters[cl1], 
                    ]) - sdeviations[[secondClusters[cl1]]], as.numeric(clusterMeans[secondClusters[cl1], 
                    ]) + sdeviations[[secondClusters[cl1]]], as.numeric(newData[row, 
                    ]))
                }
                for (cl2 in seq_along(thirdClusters)) {
                  if (thirdClusters[cl2] == 0) 
                    next
                  col <- cl2 + length(secondClusters)
                  m2[row, col] <- distToRect(as.numeric(clusterMeans[thirdClusters[cl2], 
                    ]) - sdeviations[[thirdClusters[cl2]]], as.numeric(clusterMeans[thirdClusters[cl2], 
                    ]) + sdeviations[[thirdClusters[cl2]]], as.numeric(newData[row, 
                    ]))
                }
            }
            minimalClDistance <- apply(m2, 1, function(x) which(x == min(x, na.rm = TRUE))[1])
            
            helper5 <- c(secondClusters, thirdClusters)
            for (r in seq_along(minimalDistance)) {
                if (m2[r, minimalClDistance[r]] <= m[r, minimalDistance[r]]) {
                  result[as.numeric(rownames(newData)[r])] <- helper5[minimalClDistance[r]]
                  next
                }
                distPoint <- euc.dist(newData[r, ], clusterMeans[secondClusters[helper3[minimalDistance[r]]], 
                  ] + sdeviations[[secondClusters[helper3[minimalDistance[r]]]]])
                distTotal <- distPoint + euc.dist(newData[r, ], clusterMeans[thirdClusters[helper4[minimalDistance[r]]], 
                  ])
                if (distPoint <= distanceParam * distTotal) {
                  result[as.numeric(rownames(newData)[r])] <- secondClusters[helper3[minimalDistance[r]]]
                } else {
                  if (ncol(m) > 1) {
                    minAnd2ndMin <- sort(m[r, ], index.return = TRUE)
                    if ((minAnd2ndMin$x[1]/minAnd2ndMin$x[2] >= similarityParam && 
                      helper4[minAnd2ndMin$ix[1]] != helper4[minAnd2ndMin$ix[2]]) || 
                      ((minAnd2ndMin$x[1] <= scalingHelper && minAnd2ndMin$x[2] <= 
                        scalingHelper) && helper4[minAnd2ndMin$ix[1]] != helper4[minAnd2ndMin$ix[2]])) {
                      remove <- c(remove, as.numeric(rownames(newData)[r]))
                    }
                  }
                  result[as.numeric(rownames(newData)[r])] <- thirdClusters[helper4[minimalDistance[r]]]
                }
            }
        }
    }
    
    ## tertiary to quaternary cluster (4-plex), secondary to quarternary (3-plex),
    ## primary to quarternary (2-plex) :
    if (!is.null(fourthCluster) && sum(fourthCluster) > 0) {
        if (!is.null(thirdClusters) && sum(thirdClusters) > 0) {
            newData <- subset(data, !result %in% c(emptyDroplets, firstClusters, 
                secondClusters))
            prevClusters <- thirdClusters
        } else if (!is.null(secondClusters) && sum(secondClusters) > 0) {
            newData <- subset(data, !result %in% c(emptyDroplets, firstClusters))
            prevClusters <- secondClusters
        } else if (!is.null(firstClusters) && sum(firstClusters) > 0) {
            newData <- subset(data, !result %in% c(emptyDroplets))
            prevClusters <- firstClusters
        } else {
            return(list(result = result, remove = unique(remove)))
        }
        for (cl in prevClusters) {
            if (cl == 0) 
                next
            sdC <- sdeviations[[cl]]
            newData <- subset(newData, !(newData[, 1] < (clusterMeans[cl, 1] + sdC[1]) & 
                newData[, 2] < (clusterMeans[cl, 2] + sdC[2])))
        }
        for (cl in fourthCluster) {
            sdC <- 2 * c(stats::sd(data[result == cl, 1], na.rm = TRUE), stats::sd(data[result == 
                cl, 2], na.rm = TRUE))
            if (is.na(sum(sdC))) 
                next
            thisCluster <- which(rownames(newData) %in% which(result == cl))
            S <- stats::var(newData[thisCluster, , drop = FALSE])
            mahaDist <- stats::mahalanobis(newData, clusterMeans[cl, ], S)
            posEv <- which(mahaDist < mean(mahaDist[thisCluster]))
            result[as.numeric(rownames(newData[posEv, ]))] <- cl
            newData <- newData[-posEv, , drop = FALSE]
            newData <- subset(newData, !(newData[, 1] < clusterMeans[cl, 1] + sdC[1] & 
                newData[, 1] > clusterMeans[cl, 1] - sdC[1] & newData[, 2] < clusterMeans[cl, 
                2] + sdC[2] & newData[, 2] > clusterMeans[cl, 2] - sdC[2]))
        }
        if (nrow(newData) > 0) {
            m <- matrix(nrow = nrow(newData), ncol = length(prevClusters))
            for (cl in seq_along(prevClusters)) {
                # if (prevClusters[cl]==0) next
                for (row in seq_len(nrow(newData))) {
                  m[row, cl] <- distToLineSegment(as.numeric(newData[row, ]), as.numeric(clusterMeans[fourthCluster, 
                    ]), as.numeric(clusterMeans[prevClusters[cl], ]))
                }
            }
            minimalDistance <- apply(m, 1, function(x) which(x == min(x, na.rm = TRUE))[1])
            for (r in seq_along(minimalDistance)) {
                distPoint <- euc.dist(newData[r, ], clusterMeans[prevClusters[minimalDistance[r]], 
                  ] + sdeviations[[prevClusters[minimalDistance[r]]]])
                distTotal <- distPoint + euc.dist(newData[r, ], clusterMeans[fourthCluster, 
                  ])
                if (distPoint <= distanceParam * distTotal) {
                  result[as.numeric(rownames(newData)[r])] <- prevClusters[minimalDistance[r]]
                } else {
                  if (ncol(m) > 1) {
                    minAnd2ndMin <- sort(m[r, ], index.return = TRUE)
                    if (minAnd2ndMin$x[1]/minAnd2ndMin$x[2] >= similarityParam || 
                      (minAnd2ndMin$x[1] <= scalingHelper && minAnd2ndMin$x[2] <= 
                        scalingHelper)) {
                      remove <- c(remove, as.numeric(rownames(newData)[r]))
                    }
                  }
                  result[as.numeric(rownames(newData)[r])] <- fourthCluster
                }
            }
        }
    }
    list(result = result, remove = unique(remove))
}

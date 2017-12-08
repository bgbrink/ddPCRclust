## Part of the ddPCRclust algorithm
## Author: Benedikt G Brink, Bielefeld University
## Contributor: Justin Meskas
## November 2017

# A collection of useful functions

#------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
  start_time <- as.POSIXct(start_time)
  dt <- difftime(Sys.time(), start_time, units = "secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt, tz = "GMT"), "%H:%M:%S")
}
TimeOutput(Sys.Date())
#------------------------------------------------------------------------------------------------------------------------------------------------------------#
distToLineSegment <- function(x, v, w) {
  # Return minimum distance between line segment vw and point x
  l2 <- distSquared(v, w)  # i.e. |w-v|^2 -  avoid a sqrt
  if (is.na(l2)) {
    return (Inf)
  } else if (l2 == 0.0) {
    return (euc.dist(x, v))   # v == w case
  }
  # Consider the line extending the segment, parameterized as v + t (w - v).
  # We find projection of point x onto the line.
  # It falls where t = [(x-v) . (w-v)] / |w-v|^2
  #   t = dotprod(x - v, w - v) / l2
  t = ((x[2] - v[2]) * (w[2] - v[2]) + (x[1] - v[1]) * (w[1] - v[1])) / l2
  if (t < 0.0) {
    return (euc.dist(x, v))       # Beyond the 'v' end of the segment
  } else if (t > 1.0) {
    return (euc.dist(x, w))  # Beyond the 'w' end of the segment
  }
  projection <- v + t * (w - v)  # Projection falls on the segment
  #   projection <- v*(1-t) + t * (w)  # Projection falls on the segment
  return (euc.dist(x, projection))
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------#
euc.dist <- function(x1, x2) {
  return(sqrt(sum((x1 - x2) ^ 2)))
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------#
distSquared <- function(a, b) {
  return((a[1] - b[1]) ^ 2 + (a[2] - b[2]) ^ 2)
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------#
dotprod <- function(a, b) {
  return(t(a) %*% b)
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------#
unitLen <- function (a, b) {
  if (euc.dist(a, b) == 0) {
    return (0)
  }
  return(c((b[1] - a[1]), (b[2] - a[2])) / euc.dist(a, b))
}
#------------------------------------------------------------------------------------------------------------------------------------------------------------#
insertRow <- function(existingDF, newrow, r) {
  existingDF <- rbind(existingDF, newrow)
  existingDF <- existingDF[order(c(1:(nrow(existingDF) - 1), r - 0.5)), ]
  row.names(existingDF) <- 1:nrow(existingDF)
  return(existingDF)
}

distToRect <- function(rectMin, rectMax, point) {
  dx <- max(rectMin[1] - point[1], 0, point[1] - rectMax[1])
  dy <- max(rectMin[2] - point[2], 0, point[2] - rectMax[2])
  return(sqrt(dx ^ 2 + dy ^ 2))
}

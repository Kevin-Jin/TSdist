median.series <- function(sequences) {
  # first detrend each series
  detrended <- lapply(sequences, function(series) {
    # average change in elevation for each sample in our time series
    step <- (tail(series, 1) - head(series, 1)) / length(series)
    # subtract that change from each sample, starting with the second sample
    series - (step * (0:(length(series) - 1)))
  })
  
  # indices of the global maximum of each series
  globalMaxLocations <- unlist(lapply(detrended, which.max))
  globalMinLocations <- unlist(lapply(detrended, which.min))
  
  # pick the series whose extremes are the closest to the median extremes
  return(sequences[[which.min(apply(cbind(min = globalMinLocations, max = globalMaxLocations), 1, function(seriesExtremes) {
    (seriesExtremes["min"] - median(globalMinLocations)) ^ 2 + (seriesExtremes["max"] - median(globalMaxLocations)) ^ 2
  }))]])
}

barycenter.average <- function(sequences, initial = median.series(sequences), distance = "dtw", maxIterations = 100, ...) {
  average <- initial
  totalDist <- -1
  iteration <- 0
  
  repeat {
    paths <- lapply(sequences, function(series) list(series = series, path = tsDistances(average, series, distance = distance, lead.lag.info = TRUE, ...)[["full.path"]]))
    # a sample on the new average is simply the arithmetic mean of the samples,
    # from sequences, matched to the point on the old average
    sumPoints <- rep(0, length(average))
    nPoints <- rep(0, length(average))
    for (seriesAndPath in paths) {
      path <- seriesAndPath[["path"]]
      path <- path[, path[1, ] != 0 & path[2, ] != 0]
      series <- seriesAndPath[["series"]]
      for (coordinate in split(path, col(path))) {
        sumPoints[coordinate[1]] <- sumPoints[coordinate[1]] + series[coordinate[2]]
        nPoints[coordinate[1]] <- nPoints[coordinate[1]] + 1
      }
    }
    # calculate our centroid (barycenter)
    average <- sumPoints / nPoints
    
    # if Euclidean distance reaches equilibrium, then stop
    prevTotalDist <- totalDist
    totalDist <- sum(unlist(lapply(sequences, function(series) sum(unlist(lapply(1:length(series), function(i) (series[i] - average[i]) ^ 2))))))
    # if we've hit the maximum number of iterations, then stop
    iteration <- iteration + 1
    if (totalDist == prevTotalDist || iteration >= maxIterations)
      break
  }
  
  return(average)
}

erp.barycenter.average <- function(sequences) {
  barycenter.average(sequences, distance = "erp", g = 0)
}

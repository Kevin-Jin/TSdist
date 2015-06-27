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
    
    # each sample on the new average is simply the arithmetic mean of all samples
    # that are "matched" (by time warping) to the sample on the old average with
    # the corresponding timestamp
    individualSeriesMatches <- lapply(paths, function(seriesAndPath) {
      thisSeries <- seriesAndPath[["series"]]
      path <- seriesAndPath[["path"]]
      # name the rows so code below looks a little more sane
      rownames(path) <- c("Timestamp.On.Avg", "Timestamp.On.ThisSeries")
      
      # 0 indicates that a point wasn't matched from average series or thisSeries
      # to the other series, so just ignore matchings with timestamp index of 0
      path <- path[, path["Timestamp.On.Avg", ] != 0 & path["Timestamp.On.ThisSeries", ] != 0]
      
      # creates a numeric matrix: first column is some timestamp index on average series to find matches on,
      # second column is sum of sample values matched on thisSeries to the corresponding average series sample,
      # third column is number of samples on thisSeries matched to the corresponding average series sample
      thisSeriesMatches <- as.matrix(do.call(data.frame, aggregate(x = path["Timestamp.On.ThisSeries", ], by = list(path["Timestamp.On.Avg", ]), function(Matched.Timestamps.On.ThisSeries)
        c(sum(thisSeries[Matched.Timestamps.On.ThisSeries]), length(Matched.Timestamps.On.ThisSeries))
      )))
      colnames(thisSeriesMatches) <- c("Timestamp.On.Avg", "Sum.Matched.Samples", "Num.Matched.Samples")
      thisSeriesMatches
    })
    initialSumAndNum <- list(matrix(0, nrow = length(average), ncol = 2))
    combinedSamples <- Reduce(function(accSeries, addSeries) {
      # ensure addSeries has the correct length and timestamps are in the correct order
      addSeries <- addSeries[match(1:nrow(accSeries), addSeries[, "Timestamp.On.Avg"]), c("Sum.Matched.Samples", "Num.Matched.Samples")]
      
      # use 0s where addSeries is missing matches to average series
      addSeries[is.na(addSeries)] <- 0
      accSeries + addSeries
    }, c(initialSumAndNum, individualSeriesMatches))
    # once sum(x) and num(x) is found, we can calculate our new average series!
    average <- combinedSamples[, "Sum.Matched.Samples"] / combinedSamples[, "Num.Matched.Samples"]
    
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

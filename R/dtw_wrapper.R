# implement this in path_info.c if R performance is too poor
append.path.info <- function(x, y, dtw) {
  # for consistency with x.leads and y.leads of ERP distance,
  # we leave x.leads[1] == y.leads[1] == 0. Unlike ERP distance,
  # DTW distance must start on index 1 of both time series and
  # cannot leave either time series out in the beginning for the
  # sake of minimizing distances.
  i <- 2
  j <- 2
  x.leads <- numeric(length(x) + 1)
  y.leads <- numeric(length(y) + 1)
  # impossible to leave any time series out, so index on either
  # time series can never be 0
  x.leads[1] <- y.leads[1] <- 0
  x.gaps <- 0
  y.gaps <- 0
  
  for (step.type in dtw$stepsTaken) {
    if (step.type == 1) {
      # no.gap
      i <- i + 1
      j <- j + 1
      x.leads[i] <- y.leads[j] <- 0
    } else if (step.type == 2) {
      x.gaps <- x.gaps + 1
      j <- j + 1
      x.leads[i] <- x.leads[i] + 1
      y.leads[j] <- 0
    } else if (step.type == 3) {
      y.gaps <- y.gaps + 1
      i <- i + 1
      x.leads[i] <- 0
      y.leads[j] <- y.leads[j] + 1
    }
  }
  full.path <- cbind(x = dtw$index1, y = dtw$index2)
  augment.results(dtw$distance, list(full.path = full.path, x.leads = x.leads, y.leads = y.leads), x, y)
}

#Dynamic time warping distance is calculated by using dtw package
dtwDistance <- function(x, y, lead.lag.info = FALSE, ...){
  tryCatch({
    if (lead.lag.info)
      append.path.info(x, y, dtw(x, y, distance.only = FALSE, ...))
    else
      dtw(x, y, distance.only = TRUE, ...)$distance
  },
  error=function(e){print(e);NA})   
}
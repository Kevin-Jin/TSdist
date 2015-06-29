
#Dynamic time warping distance is calculated by using dtw package

dtwDistance <- function(x, y, lead.lag.info = FALSE, ...){
  tryCatch({
  if (lead.lag.info) {
    dtw <- dtw(x, y, ...)
    list(distance = dtw$distance, full.path = rbind(dtw$index1, dtw$index2))
  } else {
    dtw(x, y, ...)$distance
  }},
  error=function(e){print(e);NA})   
}
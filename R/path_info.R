match.last <- function(x, table, ...) {
  if (any(is.na(x))) {
    # since we replace duplicate values with NA, there will be ambiguity
    stop("Cannot match NA in match.last")
  }
  table[duplicated(table, fromLast = TRUE)] <- NA
  match(x, table, ...)
}

cumulative.leads <- function(full.path, tamx, tamy) {
  xCumLeads <- full.path[match.last(1:tamx, full.path[, "x"]), "y"] - 1:tamx
  yCumLeads <- full.path[match.last(1:tamy, full.path[, "y"]), "x"] - 1:tamy
  
  list(x.cum.leads = xCumLeads, y.cum.leads = yCumLeads)
}

lead.scores <- function(cumulative.leads) {
  x.sum.leads <- sum(cumulative.leads$x.cum.leads, na.rm = TRUE)
  y.sum.leads <- sum(cumulative.leads$y.cum.leads, na.rm = TRUE)
  if (x.sum.leads == 0 && y.sum.leads == 0)
    list(x.lead.score = 0, y.lead.score = 0)
  else
    list(x.lead.score = x.sum.leads / abs(y.sum.leads), y.lead.score = y.sum.leads / abs(x.sum.leads))
}

# agg.fun is applied to the vector of all values on one times series matched to
# a single value on the other time series.
get.matched <- function(full.path, x, y, agg.fun = mean) {
  # ERP distance can use index 0 when not stepping into the beginning of a time series
  to.use.x.start <- min(which(full.path[, "x"] > 0))
  to.use.y.start <- min(which(full.path[, "y"] > 0))
  expanded <- cbind(
    x = c(rep(NA, to.use.x.start - 1), x[full.path[to.use.x.start:nrow(full.path), "x"]]),
    y = c(rep(NA, to.use.y.start - 1), y[full.path[to.use.y.start:nrow(full.path), "y"]])
  )
  
  to.use.x.start <- max(1, full.path[to.use.y.start, "x"])
  to.use.y.start <- max(1, full.path[to.use.x.start, "y"])
  y.by.x <- c(rep(NA, to.use.x.start - 1), unlist(lapply(to.use.x.start:length(x), function(i) agg.fun(y[full.path[full.path[, "x"] == i, "y"]]))))
  x.by.y <- c(rep(NA, to.use.y.start - 1), unlist(lapply(to.use.y.start:length(y), function(j) agg.fun(x[full.path[full.path[, "y"] == j, "x"]]))))
  
  list(matched.y.by.x = cbind(x = x, y = y.by.x), matched.x.by.y = cbind(x = x.by.y, y = y), matched.expanded = expanded)
}

goodness <- function(matched) {
  par.ts <- matched$matched.y.by.x
  rows <- apply(par.ts, 1, function(z) all(!is.na(z)))
  y.by.x.cor <- cor(par.ts[rows, "x"], par.ts[rows, "y"])
  
  par.ts <- matched$matched.x.by.y
  rows <- apply(par.ts, 1, function(z) all(!is.na(z)))
  x.by.y.cor <- cor(par.ts[rows, "x"], par.ts[rows, "y"])
  
  par.ts <- matched$matched.expanded
  rows <- apply(par.ts, 1, function(z) all(!is.na(z)))
  expanded.cor <- cor(par.ts[rows, "x"], par.ts[rows, "y"])
  
  list(y.by.x.cor = y.by.x.cor, x.by.y.cor = x.by.y.cor, expanded.cor = expanded.cor)
}

augment.results <- function(distance, result, x, y) {
  cum.leads <- cumulative.leads(result$full.path, length(x), length(y))
  matched <- get.matched(result$full.path, x, y)
  extras <- c(result, cum.leads, lead.scores(cum.leads), matched, goodness(matched))
  attributes(distance) <- extras
  distance
}

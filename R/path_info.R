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

augment.results <- function(result, tamx, tamy) {
  cum.leads <- cumulative.leads(result$full.path, tamx, tamy)
  c(result, cum.leads, lead.scores(cum.leads))
}


erpDistance <- function(x, y, g, sigma, lead.lag.info = FALSE){
  
  if (class(try(erpInitialCheck(x, y, g, sigma)))=="try-error"){
    return(NA)
  }else{

  #The length of the series are defined
  tamx <- length(x)
  tamy <- length(y)
  
  #The local distance matrix is defined by using the Euclidean distance.
  distMatrix <- as.vector(t(proxy::dist(x, y, method="euclidean")))

  #The cost matrix is initialized and converted into a vector
  costMatrix <- c(1:((tamx+1) * (tamy+1))) * 0 + (max(distMatrix) * 
                                                length(distMatrix))
  pathMatrix <- rep(0, (tamx+1) * (tamy+1))
  #Benefit of keeping track of both gap measures is that we can determine the
  #length of the optimized path by using only xGaps(length(xGaps)),
  #yGaps(length(yGaps)), and tamy or tamx
  xGaps <- rep(0, (tamx+1) * (tamy+1))
  yGaps <- rep(0, (tamx+1) * (tamy+1))

  
  #The case with no temporal constraint
  if (missing(sigma)){
    #The cost matrix is computed using dynammic programming.
    resultList<-.C("erpnw", as.double(x), as.double(y), as.integer(tamx),
                   as.integer(tamy), as.double(costMatrix), 
                   as.double(distMatrix), as.integer(pathMatrix), as.integer(xGaps), as.integer(yGaps), as.double(g))
    costMatrix<-resultList[[5]]
    pathMatrix<-resultList[[7]]
    xGaps<-resultList[[8]]
    yGaps <- resultList[[9]]
    
    #The case with a temporal constraint
  } else {
    #The cost matrix is computed using dynammic programming.
    resultList<-.C("erp", as.double(x), as.double(y), as.integer(tamx),
                   as.integer(tamy), as.integer(sigma),as.double(costMatrix), 
                   as.double(distMatrix), as.integer(pathMatrix), as.integer(xGaps), as.integer(yGaps), as.double(g))
    costMatrix <- resultList[[6]]
    pathMatrix <- resultList[[8]]
    xGaps <- resultList[[9]]
    yGaps <- resultList[[10]]
  }

  #The last position of the cost matrix is returned as the distance between 
  #the series.
  d<-costMatrix[length(costMatrix)]
  if (!lead.lag.info) {
    return(d)
  } else {
    return(augment.results(c(distance = d, .Call("ts_xLeadOverY", as.integer(tamx), as.integer(tamy), as.integer(pathMatrix), as.integer(xGaps), as.integer(yGaps), package = "TSdist")), tamx, tamy))
  }
  }
}


# This function checks for possible initial errors: 
erpInitialCheck <- function(x, y, g, sigma){
  
  if (!is.numeric(x) | !is.numeric(y)){
    stop('The series must be numeric', call.=FALSE)
  }
  if (!is.vector(x) | !is.vector(y)){
    stop('The series must be univariate vectors', call.=FALSE)
  }
  if (!is.numeric(g)){
    stop('g must be numeric', call.=FALSE)
  }
  if (length(x) < 1 | length(y) < 1){
    stop('The series must have at least one point', call.=FALSE)
  }
  if (any(is.na(x)) | any(is.na(y))){
    stop('There are missing values in the series', call.=FALSE)
  } 
  if(!missing(sigma)){
    if((sigma)<=0){
      stop('The window size must be positive', call.=FALSE)
    }
    if((sigma + 1) > length(x)){
      stop('The window size exceeds the the length of the first series', call.=FALSE)
    }
    if((sigma + 1) > length(y)){
      stop('The window size exceeds the the length of the second series', call.=FALSE)
    }
    if(sigma < abs(length(x) - length(y))){
      stop('The window size can not be lower than the difference between the series lengths', call.=FALSE)
    }
  }
}

#' Find precipitation events
#'
#' @description This function identifies consecutive rainy days based on adily precipitation series.
#' @param data (numeric) data series of daily precipitation values.
#' @param thresholdA minimum amount of precipitation considered as rainy day (default 0).
#' @param thresholdB same as thresholdA, used to set the ending of an event.
#' @return matrix with indices corresponding to beginning, endings and length of the events.
#'
#' @export


findevents <- function(data, thresholdA=0, thresholdB=0) {
  n <- length(data)
  A <- thresholdA
  B <- thresholdB
  isevent <- F
  starts <- vector(mode='numeric', length=0)
  ends <- vector(mode='numeric', length=0)
  norain <- 0
  israin <- data>A
  start <- 1
  for (i in c(1:n)) {
    if (!is.na(israin[i])) {
      if (israin[i] && !isevent) {
        # start new event
        isevent <- T
        start <- i
      }
      if (!is.na(sum(!israin[start:i]))) {
        if (!israin[i] && isevent &&
            #sum(!is.na(data[i:{i+A}]))==0){
            sum(!israin[start:i])>B){
          # end event
          isevent <- F
          fin <- i-1
          starts <- append(starts,start)
          ends <- append(ends,fin)
        }			
      } else {
        # end event
        isevent <- F				
      }
    }
  }
  if (is.na(sum(israin))) warning('The data contained NAs')
  return(cbind(start=starts,end=ends,lenght=ends-starts+1))
}

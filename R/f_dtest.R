#' Determination of significance between models
#'
#' @description D.test determines the significance of non-stationary over stationary POT model based on the likelihood test, at significant level sign. 
#' @param nsModel non-stationary model (gpd object).
#' @param sModel stationary model (gpd object).
#' @param sign (optional) significance level (0.05 as default).
#' @return D.test result.
#' @details Hypothesis: the NSPOT model is better -> yields TRUE; H0: the stationary model is better -> yields FALSE.
#'
#' @export

D.test <- function(nsModel,sModel,sign=0.05)
{
  # likelihood distance between both models
  dist <- -2*(nsModel$nllh-sModel$nllh)
  # critical distance at the desired significance level
  #dcrit <- qchisq(1-alpha,length(nsModel$mle)-length(sModel$mle))
  dcrit <- qchisq({1-sign},length(nsModel$mle)-length(sModel$mle))
  # likelihood test: dist is significantly different from 0 if dist > dcrit
  return(dist>=dcrit)
}
#' Produce optimal contrast
#'
#' This function is a help function to automatrically produce mean vector for linear.
#' @param mu numeric values
#' @param Sin inverse of covariance matrix
#' @return optimal contrast
#'


optC<-function(mu, Sin){
  optcont<-Sin%*%(mu-sum(mu*rowSums(Sin))/sum(Sin))
  contMat<-(optcont-sum(optcont))/sqrt(sum(optcont^2))
}

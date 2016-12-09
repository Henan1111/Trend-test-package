#' Produce mean vector for linear
#'
#' This function is a help function to automatrically produce mean vector for linear.
#' @param k a numeric value specifies the number of active doses.
#' @param doses a numeric vector of specified doses.
#' @return mean vector for linear
#' @export
#' @examples k<-4
#' doses<-seq(0,k)
#' mu.linear(doses)
#'




mu.linear<-function(doses){
  #transDose<-doses
  #transDose<-c(1:k)
  #c(0, 1/k*seq(1,k, by=1))
  k<-length(doses[-1])
  mu<-matrix(doses/doses[k+1], nrow=1,ncol=k+1)
  colnames(mu)<-doses
  mu
}




#' Produce mean matrix for CPL
#'
#' This function is a help function to automatrically produce mean matrix for CPL.
#' @param k a numeric value specifies the number of active doses.
#' @param doses a numeric vector of specified doses.
#' @return mean matrix for CPL
#' @export
#' @examples k<-4
#' doses<-seq(0,k)
#' mu.piecewise(doses)
#'
mu.piecewise<-function(doses){
  #transDose<-doses
  #transDose<-c( 1:k)
  doses<-doses[-1]
  k<-length(doses)
  mu<-NULL
  for (i in 1:k){
    for (j in (i+1):(k+1)){
      #mu1<-c(rep(0, i), rep(1/(j-i), j-i-1), rep(1, k-j+2))
      #mu1<-c(rep(0, i), transDose[max(j-i-2, 0):max(j-i-1, 0)]/transDose[j-i], rep(1, k-j+2))
      mu1<-c(rep(0, i), doses/doses[j-i])
      mu1<-mu1[1:(k+1)]
      mu1[j:(k+1)]<-1
      #			print(mu1)
      mu<-rbind(mu,mu1)
      rownames(mu)<-NULL
    }
  }
  colnames(mu)<-c(0, doses )
  mu
}

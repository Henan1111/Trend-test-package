#' Calculate power for multiple contrast tests
#'
#' Calculate power for a multiple contrast test for a specified trend type.
#' @param trueMu mean values in all doses groups including placebo. It could be a vector of arbitrary values or an object from CPL.
#' @param trueModel an object of class "Mods", see Mods in DoseFinding package for detail. It is used to define mean values. Specify trueMu if you want to an arbitrary mu or an object from CPL,
#'                  or specify trueModel if you want to use an object from MCP Mods. If trueMu and trueModel both are missing, then trueMu will be an object from MCP Mod if trend is specified as
#'                  MCP Mod or an object from CPL for all other trends.
#' @param n a numeric vector of sample sizes in all doses groups including placebo.
#' @param sigma a single numeric value of residual standard deviation. It is assumed computation are made for a normal homoscedastic model with group sample sizes given by n and residual standard deviation sigma,
#'              i.e. the covariance matrix used for estiamtes is thus sigma^2*diag(1/n)
#' @param S convariance matrix for the estimates, either n and sigma or S need to be specified.
#' @param df specify degrees of freedom for multivariate t. If this argument is missing, df=0 corresponding multivariate normal.
#' @param evalDoses a numeric vector of specified doses. When this argument is missing, doses are 0 to k by 1. You can assign less levels of doses (subset of evalDoses). True mu will be adjusted correspondingly.
#' @param trend a single character string, determining the trend for the multiple contrast trend test, one of "CPL", "MCP Mod","Aug Williams", "Williams", "Dunnett", "Linear".
#' Note that when "MCP_Mod" is specified for trend, model needs to be specified.
#' @param evalModel an object of calss "Mods" if you choose trend to be MCP_Mod trend test, see Mods for details.
#' @param alternative a single character string, specifying the alternative for the multiple contrast trend test. If this argument is missing, alternative="one.sided".
#' @param alpha significance level to use. If this argument is missing, alpha=0.05.
#' @return true mean, optimal contrast, mean matrix, power.
#' @export
#' @examples N<-150
#' n0<-30
#' n1<-30
#' k<-5
#' n<-c(n0, rep(n1,k))
#' #doses
#' ds<-seq(0,5)
#'
#'
#' models<-Mods(emax=c(log(1+2*0.01),log(1+2*0.25),log(1+2*4)),
#'           linear=NULL,
#'           exponential=c(log(1+2*0.25), log(1+2*0.75), log(1+2*3.5)),
#'           sigEmax=rbind(c(log(1+2*1),25), c(log(1+2*5),25), c(log(1+2*15),25),c(log(1+2*1),5), c(log(1+2*5),5), c(log(1+2*15),5) ),
#'           doses=ds, placEff=0, maxEff=0.7)
#'
#'
#'MCP_models<-Mods(emax=c(log(1+2*0.01),log(1+2*0.25),log(1+2*4)),
#'                 linear=NULL,
#'                 exponential=c(log(1+2*0.25), log(1+2*0.75), log(1+2*3.5)),
#'                 sigEmax=rbind( c(log(1+2*1),5), c(log(1+2*5),5), c(log(1+2*15),5) ),
#'                 doses=ds, placEff=0, maxEff=0.7)
#'
#'
#'trendpow<-trendpow(trueModel=models,n=rep(30,6), sigma=1,evalDoses=ds, alpha=0.05, trend="MCP Mod", evalModel=MCP_models)
#'
#'mufrmCPL<-mu.piecewise(ds)
#'
#'trendpow<-trendpow(trueMu=mufrmCPL,n=rep(30,6), sigma=1,evalDoses=ds, alpha=0.05, trend="MCP Mod", evalModel=MCP_models)

trendpow<-function( trueMu=NULL, trueModel=NULL,  n=NULL, sigma=NULL, S=NULL, df=NULL, evalDoses=NULL,trend=c("CPL","MCP Mod","Aug Williams", "Williams", "Dunnett", "Linear"), evalModel=NULL, alternative="one.sided",alpha=0.05){

  #trueMu: used for arbitrary trueMu or object from MCP Mod
  #trueModel: an object from MCP Mod
  if(trend[1]=="MCP Mod"){
    if (is.null(evalModel))
      stop("evalModel needs to specified")
  } else {
    if(!is.null(evalModel))
      stop ("when evalModel is specified, trend needs to be MCP Mod")
  }

  ## extract covariance matrix
  if(is.null(S)){
    if(is.null(n) | is.null(sigma))
      stop("Either S or n and sigma need to be specified")
    #    if(length(n) != nD)
    #      stop("n needs to be of length nrow(contMat)")
    S <- sigma^2*diag(1/n)
    k<-length(diag(S))-1
    #df<-as.integer(sum(n)-ncol(S))
  } else { k<-length(n)-1
  if(!is.null(n)|!is.null(sigma))
    stop("Need to specify exactly one of \"S\" or \"n\" and \"sigma\"")
  if(nrow(S) != ncol(S))
    stop("S needs to be a square matrix")
  }


  #### when use an object from CPL or MCP Mod with less doses, need doses levels

  if(is.null(trueMu)){
    if(is.null(trueModel)){
      #when true mu is not given then true mu will be object from MCP Mod for MCP Mod trend or object from CPL for other trend tests
      #stop("Either trueMu or trueModel needs to be specified")
      if(trend[1]=="MCP Mod") {
        trueMu<-as.matrix((getResp(evalModel)))
      }else{
        if(is.null(evalDoses))
          evalDoses<-seq(0,k)
        trueMu<-t(mu.piecewise(evalDoses))
      }
    }else{
      #if trueModel is given
      if(!is.null(evalDoses)){
        #doses is given
        if(length(which(attributes(trueModel)$doses%in%as.character(evalDoses)))!=length(attributes(trueModel)$doses))
          clnm<-rownames(getResp(trueModel))
        idx<-which(clnm%in%as.character(evalDoses))
        trueMu<-getResp(trueModel)[idx,]
        S<-S[idx, idx]
        nD<-nrow(trueMu)
        k<-nD-1
      } else{
        evalDoses<-attributes(trueModel)$doses
        #print("Use doses in the specified model")
        trueMu<-getResp(trueModel)
        nD<-nrow(trueMu)
        k<-nD-1
      }
    }
  }else{
    # trueMu is given
    trueMu<-t(as.matrix(trueMu))
    if (!is.null(trueModel)){
      stop("Need to specify only one of \"trueMu\" or \"trueModel\"")
    }else{
      # print("An arbitrary trueMu or object from CPL is given")
      # print("trueMu need to be a matrix with nrow= number of doses")
      if(is.null(evalDoses)){
        ### dose is not given
        nD<-nrow(trueMu)  # of doses
        k<-nD-1
        evalDoses<-seq(0,k)
        print("Use default doses from 0 to k by 1")
      }else{
        nD<-length(evalDoses)
        k<-nD-1
        clnm<-colnames(trueMu)
        idx<-which(clnm%in%as.character(evalDoses))
        trueMu<-matrix(trueMu[,idx], nrow=length(idx))
        S<-S[idx, idx]
        # if(is.null(n) & is.null(sigma)){
        #   df<-df
        # }else{
        #   df<-as.integer(sum(n[idx])-nD)
        # }
      }
    }
  }







  ####prepare muMat and contMat for williams
  cont.Mat.Will<-matrix(0, nrow=k, ncol=k+1)
  mu.Mat.Will<-matrix(0, nrow=k, ncol=k+1)
  for(j in 1:k){

    Sin<-solve(S)[c(1,((k-j+2):(k+1))),c(1,((k-j+2):(k+1)))]

    #Sin<-diag(1/(diag(S)[c(1,((k-j+2):(k+1)))]))
    mu.Mat.Will[j,]<-c(0, rep(NA, k-j), rep(1,j))
    mu1<-c(0, rep(1,j))
    cont.Mat.Will[j,c(1,((k-j+2):(k+1)))]<-optC(mu1, Sin)
  }
  #####prepare muMat and contMat for Aug-Williams
  cont.Mat1<-matrix(0, nrow=k, ncol=k+1)
  mu.Mat1<-matrix(0, nrow=k, ncol=k+1)
  for(j in 1:k){
    Sin<-solve(S)[c(1:j,k+1),c(1:j,k+1)]
    #	Sin<-diag(1/(diag(S)[c(1:j,k+1)]))
    mu1<-c(rep(0,j),1)
    mu.Mat1[j,]<-c(rep(0,j), rep(NA, k-j), 1)
    cont.Mat1[j,c(1:j,k+1)]<-optC(mu1, Sin)
  }
  cont.Mat2<-matrix(0, nrow=k-1, ncol=k+1)
  mu.Mat2<-matrix(0, nrow=k-1, ncol=k+1)
  for(j in 1:(k-1)){

    Sin<-solve(S)[c(1,((k-j+1):(k+1))), c(1,((k-j+1):(k+1)))]
    #Sin<-diag(1/(diag(S)[c(1,((k-j+1):(k+1)))]))

    mu1<-c(0, rep(1, j+1))
    mu.Mat2[j,]<-c(0, rep(NA, k-j-1),rep(1, j+1))
    cont.Mat2[j,c(1,((k-j+1):(k+1)))]<-optC(mu1, Sin)
  }
  mu.Mat.AugWill<-rbind(mu.Mat1, mu.Mat2)
  cont.Mat.AugWill<-rbind(cont.Mat1, cont.Mat2)
  ####prepare muMat and contMat for Dunnett
  cont.Mat.Dun<-matrix(0, nrow=k, ncol=k+1)
  mu.Mat.Dun<-matrix(NA, nrow=k, ncol=k+1 )
  for (j in 1:k){
    #Sin<-solve(S)[c(1, j+1),c(1, j+1)]
    Sin<-diag(1/(diag(S)[c(1, j+1)]))
    mu1<-c(0,1)
    mu.Mat.Dun[j, c(1, (j+1))]<-c(0,1)
    cont.Mat.Dun[j,c(1, j+1)]<-optC(mu1, Sin)
  }


  muMat <- switch(trend[1],
                  "MCP Mod"=as.matrix(t(getResp(evalModel))),
                  "Linear"=mu.linear(evalDoses),
                  "Williams"=mu.Mat.Will,
                  "Aug Williams"=mu.Mat.AugWill,
                  "Dunnett"=mu.Mat.Dun,
                  "CPL"=mu.piecewise(evalDoses)
  )

  contMat <- switch(trend[1],
                    "MCP Mod"=optContr(evalModel, doses=evalDoses, S=S)$contMat,
                    "Linear"=apply(muMat,1, optC,Sin=solve(S)),
                    "Williams"=t(cont.Mat.Will),
                    "Aug Williams"=t(cont.Mat.AugWill),
                    "Dunnett"=t(cont.Mat.Dun),
                    "CPL"=apply(muMat,1, optC, Sin=solve(S))
  )





  cont<-contMat
  nC<-ncol(cont)  # of contrasts


  #calculate non-centrality parameter
  deltaMat<-t(cont)%*%trueMu
  covMat<-t(cont)%*%S%*%cont
  den<-sqrt(diag(covMat))
  deltaMat<-deltaMat/den
  if(alternative == "two.sided"){
    deltaMat <- abs(deltaMat)
  }

  diag(covMat)<-diag(covMat)+min(diag(covMat))*0.0001
  #corMat<-cov2cor(covMat+min(diag(covMat))*0.000001)
  corMat<-cov2cor(covMat)



  if(is.null(df))
    df<-0

  ###calculate critical values
  tail<-ifelse(alternative == "two.sided",
               "both.tails", "lower.tail")
  critV<-qmvt(1-alpha, tail=tail, df=df,delta=0, corr=corMat,algorithm=GenzBretz() )$quantile
  if(alternative == "two.sided"){
    lower <- rep(-critV, nC)
  } else {
    lower <- rep(-Inf, nC)
  }
  upper<-rep(critV, nC)


  #calculate power
  pow<-rep(NA,ncol(deltaMat))
  for (i in 1: ncol(deltaMat)){
    pow[i]<-1-pmvt(lower=lower, upper=upper,df=df,corr=corMat, delta=deltaMat[,i], algorithm=GenzBretz())
  }
  mylist<-list(t(trueMu), contMat, muMat, pow)
  names(mylist)<-c("true_mu", "contrast_matrix", "mu_matrix", "power")
  mylist
}


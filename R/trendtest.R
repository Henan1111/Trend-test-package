#'Performs multiple contrast test
#'
#'This function performs a multiple contrast test.  Optimal contrasts are derived based on the trend you choose.
#'@param muHat a numeric vector of expected values in all doses groups including placebo.
#'@param S convariance matrix.
#'@param df specify degree freedom for multivariate t. If this argument is missing, df=0 corresponding multivariate normal.
#'@param doses a numeric vector of specified doses. When this argument is missing, doses are the object in the specified model if model is given for MCP Mod trend test, or doses are 0 to k by 1 for all other trend test.
#'@param alternative a single character string, specifying the alternative for the multiple contrast trend test. If this argument is missing, alternative="one.sided".
#'@param trend a single character string, determining the trend for the multiple contrast trend test, one of "CPL", "MCP Mod","Aug Williams", "Williams", "Dunnett", "Linear","Arbitrary".
#' Note that when "MCP_Mod" is specified for trend, model needs to be specified. When "Arbitrary" is specified, muMat needs to be specified.
#'@param model specify models if you choose trend to be MCP Mod trend test.
#'@param muMat a numeric vector or matrix specified for mu when trend is specified to be "Arbitrary".
#'@return doses, optimal contrast, correlation contrast, test statisticis, p values
#'@export
#'@examples Dose<-c(0,2.5,5,10,20,50,100,200)
#'
#'dat<-data.frame(dose=as.factor(Dose),n=c(133,32,44,63,63,65,59,58),
#'                respn=c(13,4,5,16,12,14,14,21))
#'dat$resprate<-dat$respn/dat$n
#'
#'dat$txeff<-dat$resprate-dat$resprate[1]
#'
#'logfit<-glm(cbind(respn,n-respn)~dose+0,family=binomial,data=dat)
#'muHat<-coef(logfit)
#'S<-vcov(logfit)
#'
#'model<-Mods(
#'  logistic=c(log(1+2*10),0.5),
#'  sigEmax=c(log(1+2*10),4),
#'  emax=c(log(1+2*0.1),log(1+2*1.0)),
#'  linear=NULL,
#'  exponential=log(1+2*1.0),
#'  quadratic= -1/log(1+2*550),
#'  doses=Dose)
#'
#'test_MCP<-trendtest(muHat=muHat, S=S, doses=Dose, trend="MCP Mod", model=model)
#'
#'#with arbitrary trend
#'#arbitrary mu
#'mu<-mu.linear(doses=seq(0,7))
#'test_MCP<-trendtest(muHat=muHat, S=S, doses=transDose, trend="Arbitrary", muMat=mu)


trendtest<-function(muHat, S, df=NULL, doses=NULL, alternative=c("one.sided", "two.sided"),trend=c("CPL","MCP Mod","Aug Williams", "Williams", "Dunnett", "Linear", "Arbitrary"), model=NULL, muMat=NULL){

  if(!is.null(doses)){
    if(!is.null(model) & sum(doses[-1]==attributes(model)$doses[-1])==0)
      stop("doses in the model specified for MCP Mod are different from the specified doses")
  }else{
    if(!is.null(model)){
      doses<-attributes(model)$doses
      print("doses are not given, use doses in the model specified for MCP Mod")
    }else{
      doses<-seq(0, length(muHat)-1)
      print("doses are not given, use default values from 0 to k by 1")
    }
  }





  if(trend[1]=="MCP Mod"){
    if (is.null(model))
      stop("when trend is MCP Mod, model needs to be specified")
  } else {
    if(!is.null(model))
      stop ("when model is specified, trend needs to be MCP Mod")
  }

  if(trend[1]=="Arbitrary"){
    if(is.null(muMat))
      stop("when trend is Arbitrary, muMat needs to be specified")
  }else{
    if(!is.null(muMat))
      stop("when muMat is given, trend needs to be Arbitrary")
  }



  #if (is.null(model)) {  #create a fake model, won't use it anyway
  #   	model<-Mods(
  ##          logistic=c(log(1+2*10),0.5),
  #         doses=doses)
  # } else{
  # 	model<-model
  # }





  #alternative<-c("one.sided", "two.sided")
  #trend<-c("linear", "Williams", "Aug Williams", "Dunnett","CPL")

  #number of active doses
  k<-length(doses[-1])

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
  ####prepare muMat and contMat for arbitrary

  if(!is.null(muMat)) {
    trend<-"Arbitrary"
    mu.Mat.Arb<-matrix(muMat, ncol=length(doses))
    cont.Mat.Arb<-apply(mu.Mat.Arb, 1, optC, Sin=solve(S))
  }






  muMat <- switch(trend[1],
                  "MCP Mod"=as.matrix(t(getResp(model))),
                  "Linear"=mu.linear(doses),
                  "Williams"=mu.Mat.Will,
                  "Aug Williams"=mu.Mat.AugWill,
                  "Dunnett"=mu.Mat.Dun,
                  "CPL"=mu.piecewise(doses),
                  "Arbitrary"=mu.Mat.Arb
  )

  contMat <- switch(trend[1],
                    "MCP Mod"=optContr(model, doses=doses, S=S)$contMat,
                    "Linear"=apply(muMat, 1, optC, Sin=solve(S)),
                    "Williams"=t(cont.Mat.Will),
                    "Aug Williams"=t(cont.Mat.AugWill),
                    "Dunnett"=t(cont.Mat.Dun),
                    "CPL"=apply(muMat,1,optC, Sin=solve(S)),
                    "Arbitrary"=cont.Mat.Arb )




  nC<-ncol(contMat)


  corMat<-cov2cor(t(contMat) %*% S %*% contMat)
  tStat<-as.vector(muHat%*%contMat)/sqrt(diag(t(contMat) %*% S %*% contMat))

  if (is.null(df)) df<-0 #normal case
  lower <- switch(alternative[1],
                  one.sided = matrix(rep(-Inf, nC^2), nrow = nC),
                  two.sided = matrix(rep(-tStat, each = nC), nrow = nC))
  upper <- switch(alternative[1],
                  one.sided = matrix(rep(tStat, each = nC), nrow = nC),
                  two.sided = matrix(rep(tStat, each = nC), nrow = nC))
  pVals <- numeric(nC)
  for(i in 1:nC){
    pVals[i] <-   1 - pmvt(lower[,i], upper[,i], df = df,
                           corr = corMat)
  }
  res<-list(doses=doses, Contrast=contMat, Correlation=corMat, Statistics=tStat, Pvalues=pVals, min_pvalue=min(pVals))

  res
}


###################################












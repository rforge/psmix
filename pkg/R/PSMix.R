PSMix = function(K=2, Geno, itMax=1e4, eps=1e-6, seed=NULL, MarkerVar=FALSE, verbose=FALSE){
  if(!is.null(seed)) set.seed(seed)

  Geno = apply( Geno, 2, function(x) as.integer(factor(x)) )
  Geno[is.na(Geno)] = -9
  NM = dim(Geno)
  N = NM[1]/2; M = NM[2]

  Mn = as.vector( apply( Geno, 2, function(xx) length(unique(xx)) ) )
  MnCS = c( 0, cumsum(Mn) )
  Gkmj = matrix(0, K, sum(Mn))
  for( m in 1:M ){
    from = MnCS[m]+1; to = MnCS[m+1]
    a = matrix( 1+abs(rcauchy(K*Mn[m], scale=1000)), K, Mn[m] )
    Gkmj[, from:to] = a/rowSums(a)
  }

  if(MarkerVar){
    PIimk = array( 0, dim=c(N,M,K) )
    for(m in 1:M){
      a = matrix(1+abs(rcauchy(N*K,scale=1000)), N,K)
      PIimk[,m,] = a/rowSums(a)
    }
  
    obj = .C("ps_admix2",
            PIimk = as.double(PIimk),
            Zimak = double(2*N*M*K),
            Gkmj = as.double(Gkmj),
            as.integer(K),
            as.integer(N),
            as.integer(M),
            as.integer(Mn),
            as.integer(MnCS),
            as.integer(Geno),
            as.integer(itMax),
            as.double(eps),
            DUP=TRUE,
	    PACKAGE="PSMix")[1:3]

    a0 = array( obj$PIimk, dim=c(N,M,K) )
    a1 = array( obj$Zimak, dim=c(2*N,M,K) )
    a2 = vector("list", M)

    a3 = obj$Gkmj
    LogLik = 0
    for(m in 1:M){
      from = MnCS[m]*K+1; to = MnCS[m+1]*K
      a2[[m]] = matrix( a3[from:to], K, Mn[m] )
      j = Geno[,m]
      id1 = which(j>0)
      id2 = cbind( rep(1:N,each=2), j )
      id = id2[id1,]
      LogLik = LogLik + sum(log( (a0[,m,]%*%a2[[m]])[id] ))
    }
    if(verbose) cat("LogLikelihood = ", LogLik, "\n")
    PIimk = a0
    Gimk = apply(PIimk, 1:2, which.max)
    ## attributes(PIimk, "Zimak") = a1
    ## attributes(PIimk, "Gkmj") = a2
    return( list(AmPr=PIimk, AmId=Gimk) )
  } else{
    PIik = matrix( 1+abs(rcauchy(N*K, scale=1000)), N, K )
    PIik = PIik/rowSums(PIik)
  
    obj = .C("ps_admix1",
            PIik = as.double(PIik),
            Zimak = double(2*N*M*K),
            Gkmj = as.double(Gkmj),
            as.integer(K),
            as.integer(N),
            as.integer(M),
            as.integer(Mn),
            as.integer(MnCS),
            as.integer(Geno),
            as.integer(itMax),
            as.double(eps),
            DUP=TRUE,
	    PACKAGE="PSMix")[1:3]

    a0 = matrix(obj$PIik, N, K)
    a1 = array( obj$Zimak, dim=c(2*N,M,K) )
    a2 = vector("list", M)

    a3 = obj$Gkmj
    LogLik = 0
    for(m in 1:M){
      from = MnCS[m]*K+1; to = MnCS[m+1]*K
      a2[[m]] = matrix( a3[from:to], K, Mn[m] )
      j = Geno[,m]
      id1 = which(j>0)
      id2 = cbind( rep(1:N,each=2), j )
      id = id2[id1,]
      LogLik = LogLik + sum(log( (a0%*%a2[[m]])[id] ))
    }
    if(verbose) cat("LogLikelihood = ", LogLik, "\n")
    PIik = a0
    Gik = apply(PIik, 1, which.max)
    ## 
    return( list(AmPr=PIik, AmId=Gik) )
  }
}

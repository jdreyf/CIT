#' CIT with continuous outcome, option for one-sided test
#'
#' CIT w/ permutation results, continuous outcome, continuous L and G, option for one-sided test of TassocL
#' 
#' @param alternative indicates the alternative hypothesis for the association with outcome, T.
#' Must be one of "two.sided", "greater" or "less". You can specify just the initial letter. 
#' "greater" corresponds to positive association, "less" to negative association.
#' @export

cit_cp_alt = function(L, G, T, C=NULL, n.resampl=50, n.perm=0, rseed=NULL, 
                       alternative = c("two.sided", "less", "greater")){
  alternative = match.arg(arg=alternative, choices=c("two.sided", "less", "greater"))
  stopifnot(length(L)==length(G), length(G)==length(T))
  
  if( n.resampl < n.perm ) n.resampl = n.perm
  if( !is.null(C) ){
    mydat = as.data.frame(cbind( L, G, T, C ))
  } else mydat = as.data.frame(cbind( L, G, T ))
  for( i in 1:ncol(mydat) ) mydat[, i ] = as.numeric( mydat[, i ]  )
  
  if( !is.null(C) ){
    if(is.vector(C)) {
      C = as.data.frame( matrix( C, ncol=1) )
    } else {
      C = as.data.frame( as.matrix(C) )
    }
  }
  
  L.nms = "L" 
  C.nms=NULL
  if( !is.null(C) ) C.nms = paste("C", 1:ncol(C), sep="") 
  names(mydat) = c( L.nms,"G","T",C.nms )
  
  pvec = rep(NA,4)
  
  # pval for T ~ L
  #nm.y = "T"
  #nms.full = c(L.nms, C.nms)
  #pvec[1] = linreg( nms.full, nms.redu=C.nms, nm.y, mydat )
  pvec[1] = assoc(x=L,y=T,z=C,alternative = alternative)

  # pval for T ~ G|L
  # T ~ G is two.sided test
  nm.y = "T"
  nms.full = c("G", L.nms, C.nms)
  nms.redu = c(L.nms, C.nms)
  pvec[2] = linreg( nms.full, nms.redu, nm.y, mydat )
  
  # pval for G ~ L|T
  nm.y = "G"
  nms.full = c("T", L.nms )
  nms.redu = "T"
  pvec[3] = linreg( nms.full, nms.redu, nm.y, mydat )
  
  mydat1 = na.exclude(mydat)
  tmp = c( "T ~ G", C.nms )
  formula = paste( tmp, collapse="+" )
  fit3 = lm( formula, data=mydat1)
  tmp = c( "T ~ G", L.nms, C.nms )
  formula = paste( tmp, collapse="+" )
  fit5 = lm( formula, data=mydat1 )
  f.ind = anova(fit3,fit5)$F[2]
  
  vrs.1 = paste( L.nms, collapse="+" )
  formula1 = paste( "G ~ ", vrs.1, sep="")
  fitG = lm( formula1, data=mydat, na.action=na.exclude)
  
  coef.g = rep(NA, length(L.nms) + 1)
  coef.g[ 1 ] = summary(fitG)$coefficients["(Intercept)",1]
  #for( i in 1:length(L.nms) ) coef.g[ i + 1 ] = summary(fitG)$coefficients[ L.nms[ i ],1]
  
  for( i in 1:length(L.nms) ) {
    tmp = try( summary(fitG)$coefficients[ L.nms[ i ],1], silent = TRUE )
    tmp = strsplit( as.character( tmp ), " ", fixed=TRUE )[[ 1 ]]
    coef.g[ i + 1 ] = ifelse( length( tmp ) == 1, as.numeric(tmp), 0 )
  } # End L.nms loop
  
  mydat[, "G.r"] = resid(fitG)   
  
  fvecr = rep(NA,n.resampl)
  
  set.seed(rseed)
  
  #resample for 4th p-value
  for(rep in 1:n.resampl){
    
    nni  = sample( 1:nrow(mydat) ) 
    
    tmp = rep(0, nrow(mydat) )
    for( i in 1:length(L.nms) ) tmp = tmp + coef.g[ i + 1 ] * mydat[, L.nms[ i ] ]
    mydat[, "G.n"] = coef.g[ 1 ] + tmp + mydat[ nni, "G.r"] 
    
    # F for T ~ L|G.n
    mydat1 = na.exclude(mydat)
    tmp = c( "T ~ G.n", C.nms )
    formula = paste( tmp, collapse="+" )
    fit_0 = lm( formula, data=mydat1 )
    
    tmp = c( "T ~ G.n", L.nms, C.nms )
    formula = paste( tmp, collapse="+" )
    fit_1 = lm( formula, data=mydat1 )
    fvecr[ rep ] = anova(fit_0,fit_1)$F[2]
    
  } # End rep loop
  
  #####F Method
  fvecr = fvecr[!is.na(fvecr)]
  df1 = anova(fit3,fit5)$Df[2]
  df2 = anova(fit3,fit5)$Res.Df[2]
  fncp = mean(fvecr,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
  if(fncp < 0) fncp = 0
  
  ######### Transform F to normal
  npvals = pf(fvecr,df1,df2,ncp=fncp,lower.tail=TRUE)
  nfvecr = qnorm(npvals)
  
  npf = pf(f.ind,df1,df2,ncp=fncp,lower.tail=TRUE) #Transform observed F
  zf = qnorm(npf)
  pvec[4] = pnorm(zf,mean=mean(nfvecr),sd=sd(nfvecr))
  
  pvalc = max(pvec)  ###Causal p-value
  
  pvals = c( pvalc, pvec )
  names(pvals) = c( "p_cit", "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG")
  
  if( n.perm > 1 ){
    p.perm.ind = NA
    rep = n.resampl + 1
    
    nni  = sample( 1:nrow(mydat) )
    tmp = rep(0, nrow(mydat) )
    for( i in 1:length(L.nms) ) tmp = tmp + coef.g[ i + 1 ] * mydat[, L.nms[ i ] ]
    mydat[, "G.n"] = coef.g[ 1 ] + tmp + mydat[ nni, "G.r"] 
    
    # F for T ~ L|G.n
    mydat1 = na.exclude(mydat)
    tmp = c( "T ~ G.n", C.nms )
    formula = paste( tmp, collapse="+" )
    fit_0 = lm( formula, data=mydat1 )
    tmp = c( "T ~ G.n", L.nms, C.nms )
    formula = paste( tmp, collapse=" + " )
    fit_1 = lm( formula, data=mydat1 )
    fvecr[ rep ] = anova(fit_0,fit_1)$F[2]
    
    rand.v = sample( 1:length(fvecr) )
    
    for( perm in 1:n.perm){
      
      p.ind = rand.v[ perm ] 
      
      f.ind = fvecr[ p.ind ]
      fvecr.p = fvecr[ -p.ind ]
      fncp = mean(fvecr.p,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
      if(fncp < 0) fncp = 0
      
      ######### Transform F to normal
      npvals = pf(fvecr,df1,df2,ncp=fncp,lower.tail=TRUE)
      nfvecr = qnorm(npvals)
      
      npf = pf(f.ind,df1,df2,ncp=fncp,lower.tail=TRUE) #Transform perm stat F
      zf = qnorm(npf)
      p.perm.ind[ perm ] = pnorm(zf,mean=mean(nfvecr),sd=sd(nfvecr))
    } # End perm loop
    
    ########## permutation pvals for T ~ L, T ~ G|L, and G ~ L|T
    # compute residuals and coefficients from fit
    p.perm.TasscL = NA
    p.perm.TasscGgvnL = NA
    p.perm.GasscLgvnT = NA
    
    nm.y.1 = "T"
    nms.full.1 = c( L.nms, C.nms)
    nms.redu.1 = C.nms
    
    nm.y.2 = "T"
    nms.full.2 = c("G", L.nms, C.nms)
    nms.redu.2 = c(L.nms, C.nms)
    
    nm.y.3 = "G"
    nms.full.3 = c("T", L.nms)
    nms.redu.3 = "T"
    
    for( perm in 1:n.perm){ 
      
      nni  = sample( 1:nrow(mydat) ) 
      mydat.p = mydat
      
      mydat.p[ , L.nms ] = mydat[ nni , L.nms ]
      #p.perm.TasscL[perm] = linreg( nms.full.1, nms.redu.1, nm.y.1, mydat.p )
      if (!is.null(C)){
        p.perm.TasscL[perm] = assoc(mydat.p$L, mydat.p$T, alternative = alternative)
      } else {
        p.perm.TasscL[perm] = assoc(mydat.p$L, mydat.p$T, mydat.p[,C.nms], alternative = alternative)
      }
      mydat.p[ , L.nms ] = mydat[ , L.nms ]
      
      tmp.nms = nms.full.2[ !is.element( nms.full.2, nms.redu.2 ) ]
      mydat.p[ , tmp.nms ] = mydat[ nni , tmp.nms ]
      p.perm.TasscGgvnL[perm] = linreg( nms.full.2, nms.redu.2, nm.y.2, mydat.p )
      mydat.p[ , tmp.nms ] = mydat[ , tmp.nms ]
      
      tmp.nms = nms.full.2[ !is.element( nms.full.3, nms.redu.3 ) ]
      mydat.p[ , tmp.nms ] = mydat[ nni , tmp.nms ]
      p.perm.GasscLgvnT[perm] = linreg( nms.full.3, nms.redu.3, nm.y.3, mydat.p )
      
    } # End perm loop
    
    rslts = as.data.frame( matrix(NA, ncol=(length(pvals) + 1) ) )
    names(rslts) = c( "perm", names(pvals) )
    rslts[ 1, "perm" ] = 0
    rslts[ 1, names(pvals) ] = pvals
    rslts[ 2:(n.perm + 1), "perm" ] = 1:n.perm
    rslts[ 2:(n.perm + 1), "p_TassocL" ] = p.perm.TasscL
    rslts[ 2:(n.perm + 1), "p_TassocGgvnL" ] = p.perm.TasscGgvnL
    rslts[ 2:(n.perm + 1), "p_GassocLgvnT" ] = p.perm.GasscLgvnT
    rslts[ 2:(n.perm + 1), "p_LindTgvnG" ] = p.perm.ind
    for(i in 2:(n.perm+1)) rslts[ i, "p_cit" ] = max( rslts[ i, c( "p_TassocL", "p_TassocGgvnL", "p_GassocLgvnT", "p_LindTgvnG") ] )
    pvals = rslts
    
  } # End if n.perm
  
  return(pvals)
  
} # End function cit.cp
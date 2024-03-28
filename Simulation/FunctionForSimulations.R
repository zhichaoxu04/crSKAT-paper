

genData <- function(seed=NULL, n, alpha1, alpha2, beta1, beta2) {
  if (!is.null(seed)) { 
    set.seed(seed)
  }
  
  # assume as in Hudgens et al. paper that \sum Fk(\infty) = 1, so one of events will always happen,
  # of course alphas and betas are the determinants of whether this is true
  pi1 <- 1 - exp(beta1 / alpha1)
  pi2 <- 1 - exp(beta2 / alpha2)
  tempType <- rbinom(n=n, size=1, prob=pi2) + 1
  tempUnif <- runif(n=n)
  tempTime <- ifelse(tempType == 1, log( 1 - log((1 - tempUnif * pi1)) * alpha1 / beta1 ) / alpha1,
                     log( 1 - log((1 - tempUnif * pi2)) * alpha2 / beta2 ) / alpha2)
  
  return(list(tempTime = tempTime, tempType = tempType))
}



# I have to get rid of the 999, for right-censored observations the new right time is now the old left time 
# and the new left time is 0
createIntNew <- function(obsTimes, eventTime) {
  # order the times in case the random portion causes them to go out of order
  orderedTimes <- sort(obsTimes)
  nzero <- sum(orderedTimes == 0)

  # left end of interval
  minIdx <- which(orderedTimes < eventTime)
  if (length(minIdx)-nzero == 0) {
    minTime <- 0
    tpos <- 0
  } else {
    minTime <- orderedTimes[max(minIdx)]
    tpos <- 1
  }
  # right end of interval
  maxIdx <- which(orderedTimes >= eventTime)
  if (length(maxIdx) == 0) {   
    maxTime <- max(orderedTimes)  # Right Censored
    minTime <- 0
    obsind <- 0
  }else{
    maxTime <- orderedTimes[min(maxIdx)]
    obsind <- 1
  }
  
  return(c(minTime, maxTime, tpos, obsind))
}

# I have to get rid of the 999, for right-censored observations the old left time is now the new right time 
makeICdmat <- function(xMat, lt, rt, obsInd, quant_r = NULL, nKnots = 1) {
  
  # place the knots at equally spaced quantiles
  if (is.null(quant_r)) {
    quant_r <- stats::quantile(log(rt[obsInd == 1]), probs=seq(from=0, to=1, length.out=nKnots+2))
  }
  
  # a0 and a1 are always there
  right_a0 <- 1
  #right_a1 <- ifelse(obsInd == 0, 999, log(rt))
  right_a1 <- log(rt)
  if (is.null(xMat)) {
    right_dmat <- cbind(right_a0, right_a1)
  } else {
    right_dmat <- cbind(xMat, right_a0, right_a1)
  }
  
  # if lt = 0, then the cumulative hazard is necessarily 0 so set it all to 0
  left_a0 <- ifelse(lt == 0, 0, 1)
  left_a1 <- ifelse(lt == 0, 0, log(lt))
  if (is.null(xMat)) {
    left_dmat <- cbind(left_a0, left_a1)
  } else {
    left_dmat <- cbind(xMat, left_a0, left_a1)
  }
  
  kmax <- max(quant_r)
  kmin <- min(quant_r)
  # place the knots
  for (j in 1:nKnots) {
    ej <- (kmax - quant_r[j+1]) / (kmax - kmin)
    
    #right_aj <-  ifelse(obsInd == 0, 999,
    #                    pmax(0, (right_a1 - quant_r[j+1])**3) - ej * pmax(0, (right_a1 - kmin)**3) -
    #                      (1 - ej) * pmax(0, (right_a1 - kmax)**3))
    right_aj <- pmax(0, (right_a1 - quant_r[j+1])**3) - ej * pmax(0, (right_a1 - kmin)**3) -
      (1 - ej) * pmax(0, (right_a1 - kmax)**3)
    right_dmat <- cbind(right_dmat, right_aj)
    
    left_aj <- ifelse(lt == 0, 0, pmax(0, (left_a1 - quant_r[j+1])**3) - ej * pmax(0, (left_a1 - kmin)**3) -
                        (1 - ej) * pmax(0, (left_a1 - kmax)**3))
    left_dmat <- cbind(left_dmat, left_aj)
  }
  
  return(list(right_dmat=right_dmat, left_dmat=left_dmat, quant_r = quant_r))
}


# score equation using spline model
# if putting the genetic fixed effect in for the second outcome, specify the gMat
# deltaVec should have 0, 1, 2
scoreEqSpline <- function(x, leftDmat, rightDmat, gSummed=NULL, gMat=NULL, estG=FALSE, leftTimes, deltaVec) {
  
  # different designs for different outcomes
  # here there should be no "infinity" type terms in the rightDmat
  
  leftDmat2 <- cbind(leftDmat, gSummed)
  rightDmat2 <- cbind(rightDmat, gSummed)
  if (!is.null(gMat) & estG) {
    leftDmat1 <- cbind(leftDmat, gMat)
    rightDmat1 <- cbind(rightDmat, gMat)
  } else {
    leftDmat1 <- leftDmat
    rightDmat1 <- rightDmat
  }
  
  theta1 <- x[1:ncol(leftDmat1)]
  theta2 <- x[(ncol(leftDmat1)+1):length(x)]
  
  # separate deltas
  delta0 <- ifelse(deltaVec == 0, 1, 0)
  delta1 <- ifelse(deltaVec == 1, 1, 0)
  delta2 <- ifelse(deltaVec == 2, 1, 0)
  
  # H and S terms
  H1L <- ifelse(leftTimes == 0, 0, exp(leftDmat1 %*% theta1))
  #H1L <- exp(leftDmat %*% theta1)
  # have to use 999 because Inf * 0  = Inf
  #H1R <- ifelse(deltaVec == 0, 999, exp(rightDmat %*% theta1))
  H1R <- exp(rightDmat1 %*% theta1)
  S1L <- ifelse(leftTimes == 0, 1, exp(-H1L))
  S1R <- exp(-H1R)
  H2L <- ifelse(leftTimes == 0, 0, exp(leftDmat2 %*% theta2))
  #H2L <- exp(leftDmat %*% theta2)
  #H2R <- ifelse(deltaVec == 0, 999, exp(rightDmat %*% theta2))
  H2R <- exp(rightDmat2 %*% theta2)
  S2L <- ifelse(leftTimes == 0, 1, exp(-H2L))
  S2R <- exp(-H2R)
  
  # to be swept
  L1term <- ifelse(leftTimes == 0, 0, -delta1 * H1L * S1L / (S1L - S1R))
  R1term <- delta1 * H1R * S1R / (S1L - S1R)
  E1term <- -delta0 * H1R * S1R / (S1R + S2R - 1)
  L2term <- ifelse(leftTimes == 0, 0, -delta2 * H2L * S2L / (S2L - S2R))
  R2term <- delta2 * H2R * S2R / (S2L - S2R)
  E2term <- -delta0 * H2R * S2R / (S1R + S2R - 1)
  # sweep it
  U1swept <- sweep(leftDmat1, MARGIN=1, STATS=L1term, FUN="*") + 
    sweep(rightDmat1, MARGIN=1, STATS=R1term, FUN="*") + 
    sweep(rightDmat1, MARGIN=1, STATS=E1term, FUN="*")
  U1 <- colSums(U1swept)
  U2swept <- sweep(leftDmat2, MARGIN=1, STATS=L2term, FUN="*") + 
    sweep(rightDmat2, MARGIN=1, STATS=R2term, FUN="*") +
    sweep(rightDmat2, MARGIN=1, STATS=E2term, FUN="*")
  U2 <- colSums(U2swept)
  
  # this is just to get the score equations for \gamma under the null
  if (!is.null(gMat) & !estG) {
    UgSwept <- sweep(gMat, MARGIN=1, STATS=L1term + R1term + E1term, FUN="*") 
    uVec <- c(U1, colSums(UgSwept), U2)
  } else {
    uVec <- c(U1, U2)
  }
  return(uVec)
  
}


#---------------------------------------------------------------------#
# function to calculate the information
calcInfo <- function(leftDmat, rightDmat, gSummed, gMat = NULL, leftTimes, theta1, theta2, deltaVec) {
  
  # design matrices
  leftDmat2 <- cbind(leftDmat, gSummed)
  rightDmat2 <- cbind(rightDmat, gSummed)
  if (!is.null(gMat)) {
    leftDmat1 <- cbind(leftDmat, gMat)
    rightDmat1 <- cbind(rightDmat, gMat)
  } else {
    leftDmat1 <- leftDmat
    rightDmat1 <- rightDmat
  }
  
  # separate deltas
  delta0 <- ifelse(deltaVec == 0, 1, 0)
  delta1 <- ifelse(deltaVec == 1, 1, 0)
  delta2 <- ifelse(deltaVec == 2, 1, 0)
  
  # H and S terms
  H1L <- ifelse(leftTimes == 0, 0, exp(leftDmat1 %*% theta1))
  #H1L <- exp(leftDmat %*% theta1)
  # have to use 999 because Inf * 0  = Inf
  #H1R <- ifelse(deltaVec == 0, 999, exp(rightDmat %*% theta1))
  H1R <- exp(rightDmat1 %*% theta1)
  S1L <- ifelse(leftTimes == 0, 1, exp(-H1L))
  S1R <- exp(-H1R)
  H2L <- ifelse(leftTimes == 0, 0, exp(leftDmat2 %*% theta2))
  #H2L <- exp(leftDmat %*% theta2)
  #H2R <- ifelse(deltaVec == 0, 999, exp(rightDmat %*% theta2))
  H2R <- exp(rightDmat2 %*% theta2)
  S2L <- ifelse(leftTimes == 0, 1, exp(-H2L))
  S2R <- exp(-H2R)
  
  # dt1dt1 part
  d11ToSweep1 <- ifelse(leftTimes == 0, 0, delta1 * H1L * S1L * (1 - H1L)) / (S1L - S1R)
  d11ToSweep2 <- delta1 * H1R * S1R * (1 - H1R) / (S1L - S1R)
  d11ToSweep3 <- ifelse(leftTimes == 0, 0, delta1 * H1L * S1L) / (S1L - S1R)^1
  d11ToSweep4 <- delta1 * H1R * S1R / (S1L - S1R)^1
  d11ToSweep5 <- delta0 * H1R * S1R * (1 - H1R) / (S1R + S2R - 1)
  d11ToSweep6 <- delta0 * H1R * S1R / (S1R + S2R - 1)^1
  
  dt1dt1 <- -t(leftDmat1) %*% sweep(leftDmat1, MARGIN=1, STATS=d11ToSweep1, FUN="*") + 
    t(rightDmat1) %*% sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep2, FUN="*") + 
    -t(sweep(-leftDmat1, MARGIN=1, STATS=d11ToSweep3, FUN="*") + sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep4, FUN="*")) %*%
    (sweep(-leftDmat1, MARGIN=1, STATS=d11ToSweep3, FUN="*") + sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep4, FUN="*")) + 
    -t(rightDmat1) %*% sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep5, FUN="*") + 
    -t(sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep6, FUN="*")) %*% sweep(rightDmat1, MARGIN=1, STATS=d11ToSweep6, FUN="*")
  
  # dt2dt2 part
  d22ToSweep1 <- ifelse(leftTimes == 0, 0, delta2 * H2L * S2L * (1 - H2L)) / (S2L - S2R)
  d22ToSweep2 <- delta2 * H2R * S2R * (1 - H2R) / (S2L - S2R)
  d22ToSweep3 <- ifelse(leftTimes == 0, 0, delta2 * H2L * S2L) / (S2L - S2R)^1
  d22ToSweep4 <- delta2 * H2R * S2R / (S2L - S2R)^1
  d22ToSweep5 <- delta0 * H2R * S2R * (1 - H2R) / (S1R + S2R - 1)
  d22ToSweep6 <- delta0 * H2R * S2R / (S1R + S2R - 1)^1
  dt2dt2 <- -t(leftDmat2) %*% sweep(leftDmat2, MARGIN=1, STATS=d22ToSweep1, FUN="*") + 
    t(rightDmat2) %*% sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep2, FUN="*") + 
    -t(sweep(-leftDmat2, MARGIN=1, STATS=d22ToSweep3, FUN="*") + sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep4, FUN="*")) %*%
    (sweep(-leftDmat2, MARGIN=1, STATS=d22ToSweep3, FUN="*") + sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep4, FUN="*")) + 
    -t(rightDmat2) %*% sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep5, FUN="*") + 
    -t(sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep6, FUN="*")) %*% sweep(rightDmat2, MARGIN=1, STATS=d22ToSweep6, FUN="*")
  
  # off diagonal part
  dt1dt2 <- -t(sweep(rightDmat1, MARGIN=1, STATS=delta0 * H1R * S1R / (S1R + S2R - 1)^1, FUN="*")) %*% 
    sweep(rightDmat2, MARGIN=1, STATS=delta0 * H2R * S2R / (S1R + S2R - 1)^1, FUN="*")
  
  # put it together
  iMat <- rbind(cbind(dt1dt1, dt1dt2), cbind(t(dt1dt2), dt2dt2))
  
  return(iMat)
}


# ---- Create a class for result storage
multiResultClass <- function(thetaMatFull=NULL, thetaMatPartial=NULL, scoreMat=NULL, infoArrayFull=NULL, infoArrayPartial=NULL,
                             IgtArray=NULL, IggArray=NULL, skatDF=NULL){
  me <- list(
    thetaMatFull = thetaMatFull,
    thetaMatPartial = thetaMatPartial,
    scoreMat = scoreMat,
    infoArrayFull = infoArrayFull,
    infoArrayPartial = infoArrayPartial,
    IgtArray = IgtArray,
    IggArray = IggArray,
    skatDF = skatDF
  )
  
  ## Set the name for the class
  class(me) <- append(class(me), "multiResultClass")
  return(me)
}

multiResultClass_1 <- function(thetaMatPartial=NULL, skatDF=NULL){
  me <- list(
    # thetaMatFull = thetaMatFull,
    thetaMatPartial = thetaMatPartial,
    # scoreMat = scoreMat,
    # infoArrayFull = infoArrayFull,
    # infoArrayPartial = infoArrayPartial,
    # IgtArray = IgtArray,
    # IggArray = IggArray,
    skatDF = skatDF
  )
  
  ## Set the name for the class
  class(me) <- append(class(me), "multiResultClass_1")
  return(me)
}


sim_gmat <- function(n, q, rho, maxMAF=0.05){
  # Construct a binary correlation matrix
  cmat <- stats::toeplitz(c(1, rep(rho, q - 1))) # q * q
  meanparam1 <- runif(q, .01, maxMAF)
  meanparam2 <- runif(q, .01, maxMAF)
  x <- suppressWarnings(bindata::rmvbin(n, margprob = meanparam2, bincorr = cmat)  + 
                          bindata::rmvbin(n, margprob = meanparam2, bincorr = cmat))
  
  return(x)
}

sim_gmat_ld <- function(n, q, rho) {
  # Generate a correlation matrix based on the specified rho
  cmat <- matrix(rho, nrow=q, ncol=q)
  diag(cmat) <- 1
  
  # Simulate multivariate normal data based on the correlation matrix
  library(MASS) # For mvrnorm
  norm_data <- MASS::mvrnorm(n = n, mu = rep(0, q), Sigma = cmat)
  
  # Convert normal data to binary genotype data using a threshold
  # Assuming a biallelic system, alleles above the median are set to 1, below to 0
  geno_data <- ifelse(norm_data > median(norm_data), 1, 0)
  
  return(geno_data)
}


genData_pwr <- function(seed=NULL, n, alpha1, alpha2, beta1, beta2, gMatCausal, effectSizes) {
  if (!is.null(seed)) { 
    set.seed(seed)
  }
  
  # MAF <- apply(gMatCausal, 2, function(x) mean(x)/2)
  # Bk <- effectSizes * abs(log10(MAF))
  if(length(effectSizes) != ncol(gMatCausal)){stop("Length of effect Sizes is not equal to the gMatCausal!")}
  geneticVec <- gMatCausal %*% effectSizes
  
  # assume as in Hudgens et al. paper that \sum Fk(\infty) = 1, so one of events will always happen,
  # of course alphas and betas are the determinants of whether this is true
  pi1 <- 1 - exp(beta1 / alpha1)
  pi2 <- 1 - exp(beta2 / alpha2)
  tempType <- rbinom(n=n, size=1, prob=pi2) + 1
  tempUnif <- runif(n=n)
  tempTime <- ifelse(tempType == 1, log( 1 - alpha1 / beta1 * log((1 - tempUnif * pi1)) / exp(geneticVec) ) / alpha1,
                     log( 1 - alpha2 / beta2 * log((1 - tempUnif * pi2)) / exp(geneticVec) ) / alpha2)
  
  return(list(tempTime = tempTime, tempType = tempType))
}







# to flip SNPs
flipSNPs <- function(x) {2 - x}


# ----- WVIC -----

#----------------------- G-function -----------------------#
Gfunc=function(x,r){
  if(r==0){return(x)}
  else{return(log(1+r*x)/r)}
}

#----------------------- Derivative of G-function -----------------------#
Gdfunc=function(x,r){
  if(r==0){return(1)}
  else{return(1/(1+r*x))}
}

#----------------------- WVIC -----------------------#
s_cal_adj=function(x,jumppts,cumhaz,regcoef,r){
  Ltime=x[1]; Rtime=x[2]
  expterm=exp(t(x[-c(1,2)])%*%regcoef)
  cumhaz_L=ifelse(Ltime==0,0,cumhaz[tail(which(Ltime>=jumppts),1)])
  cumhaz_R=ifelse(Rtime==Inf,Inf,cumhaz[tail(which(Rtime>=jumppts),1)])
  chLexp=cumhaz_L*expterm
  chRexp=cumhaz_R*expterm
  if(Ltime==0)
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp)/(1-exp(-Gfunc(chRexp,r))))
  else if(Rtime==Inf)
    return(-Gdfunc(chLexp,r)*chLexp)
  else
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp-exp(-Gfunc(chLexp,r))*Gdfunc(chLexp,r)*chLexp)/(exp(-Gfunc(chLexp,r))-exp(-Gfunc(chRexp,r))))
}

s_cal_noadj=function(x,jumppts,cumhaz,r){
  Ltime=x[1]; Rtime=x[2]
  expterm=1
  cumhaz_L=ifelse(Ltime==0,0,cumhaz[tail(which(Ltime>=jumppts),1)])
  cumhaz_R=ifelse(Rtime==Inf,Inf,cumhaz[tail(which(Rtime>=jumppts),1)])
  chLexp=cumhaz_L*expterm
  chRexp=cumhaz_R*expterm
  if(Ltime==0)
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp)/(1-exp(-Gfunc(chRexp,r))))
  else if(Rtime==Inf)
    return(-Gdfunc(chLexp,r)*chLexp)
  else
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp-exp(-Gfunc(chLexp,r))*Gdfunc(chLexp,r)*chLexp)/(exp(-Gfunc(chLexp,r))-exp(-Gfunc(chRexp,r))))
}


#----------------------- WVICLT -----------------------#
s_cal_adj_trunc=function(x,jumppts,cumhaz,regcoef,r){
  V0time=x[1]; Ltime=x[2]; Rtime=x[3]
  expterm=exp(t(x[-c(1,2,3)])%*%regcoef)
  cumhaz_V0toL=ifelse(Ltime==V0time,0,cumhaz[tail(which(Ltime>=jumppts),1)]-cumhaz[tail(which(V0time>=jumppts),1)])
  cumhaz_V0toR=ifelse(Rtime==Inf,Inf,cumhaz[tail(which(Rtime>=jumppts),1)]-cumhaz[tail(which(V0time>=jumppts),1)])
  chLexp=cumhaz_V0toL*expterm
  chRexp=cumhaz_V0toR*expterm
  if(Ltime==V0time)
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp)/(1-exp(-Gfunc(chRexp,r))))
  else if(Rtime==Inf)
    return(-Gdfunc(chLexp,r)*chLexp)
  else
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp-exp(-Gfunc(chLexp,r))*Gdfunc(chLexp,r)*chLexp)/(exp(-Gfunc(chLexp,r))-exp(-Gfunc(chRexp,r))))
}

s_cal_noadj_trunc=function(x,jumppts,cumhaz,r){
  V0time=x[1]; Ltime=x[2]; Rtime=x[3]
  expterm=1
  cumhaz_V0toL=ifelse(Ltime==V0time,0,cumhaz[tail(which(Ltime>=jumppts),1)]-cumhaz[tail(which(V0time>=jumppts),1)])
  cumhaz_V0toR=ifelse(Rtime==Inf,Inf,cumhaz[tail(which(Rtime>=jumppts),1)]-cumhaz[tail(which(V0time>=jumppts),1)])
  chLexp=cumhaz_V0toL*expterm
  chRexp=cumhaz_V0toR*expterm
  if(Ltime==V0time)
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp)/(1-exp(-Gfunc(chRexp,r))))
  else if(Rtime==Inf)
    return(-Gdfunc(chLexp,r)*chLexp)
  else
    return((exp(-Gfunc(chRexp,r))*Gdfunc(chRexp,r)*chRexp-exp(-Gfunc(chLexp,r))*Gdfunc(chLexp,r)*chLexp)/(exp(-Gfunc(chLexp,r))-exp(-Gfunc(chRexp,r))))
}


Gfunc = function(x, r){
  if(r==0){ x }
  else if(r>0){ log(1+r*x)/r }
}

###################### Function for semi-param transformation model fit ######################
DentalEM_r0_adj = function(X, bound_time, order_bound_times, r_alpha, maxiter){
  #-----------------------------------------------------------#
  # X: covariates
  # bound_time: left, right bound times
  # order_bound_times: ordered baseline hazard jump time pts
  # r_alpha: semi-parametric transformation parameter
  # maxiter: maximum numer of iter
  #-----------------------------------------------------------#
  n = nrow(X)
  p = ncol(X)
  # Censor rate
  left.censor.rate = sum(bound_time[,1]==0)/n
  right.censor.rate = sum(bound_time[,2]==Inf)/n
  #---------------------- Initial Values -----------------------#
  L = length(order_bound_times)
  lambda.ini = rep(1/L, L)
  
  gamma.ini = rep(0, p)
  
  Params_old = rep(0, p)
  Params_new = gamma.ini
  
  lambda_old = rep(0, length(order_bound_times))
  lambda_new = lambda.ini
  
  #-------------------- Apply EM Algorithm --------------------#
  ### Initial values
  iter = 0
  
  ### collect all loglik values
  loglik_all = NULL
  
  message(" Now: Obtaining the unpenalized nonparametric MLE")
  ### Start iteration
  while (iter < maxiter){
    
    iter = iter + 1
    
    Params_old = Params_new
    lambda_old = lambda_new
    
    # plot(order_bound_times, cumsum(lambda_new))
    # lines(order_bound_times, order_bound_times/2)
    
    gamma = as.matrix(Params_old)
    
    Rij.star = apply(bound_time,1,function(x){ max(x[x<Inf]) })
    
    exp.term = c(exp(X %*% gamma))       # length n array
    
    ### index of jump time points between LR bounds
    jumptidx_in_LR = (outer(order_bound_times,bound_time[,1],'>') & outer(order_bound_times,bound_time[,2],'<='))+0    # len(order_bound_time) by n matrix
    
    sum.lambda.exp.left = rowSums(exp.term * t(outer(order_bound_times,bound_time[,1],'<=')*lambda_old))       #length n array
    sum.lambda.exp.right = rowSums(exp.term * t(outer(order_bound_times,bound_time[,2],'<=')*lambda_old))      #length n array
    sum.lambda.exp.between = matrix(colSums(lambda_old*jumptidx_in_LR)*exp.term, ncol=1)    #n by 1 matrix
    
    ### logLikelihood
    lik_i = rep(NA,n)
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) - exp(-Gfunc(sum.lambda.exp.right[i], r_alpha)) }
      else if(bound_time[i,2] == Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) }
    }
    loglik_i = log(lik_i)
    loglik = sum(loglik_i)
    loglik_all = c(loglik_all, loglik)
    
    ### E[C_ijt]
    E.Cijt = matrix(nrow=n, ncol=length(order_bound_times))      # n by len(order_bound_time) matrix
    for(i in 1:n){
      if(bound_time[i,2]==Inf){ E.Cijt[i,] = rep(0, length(order_bound_times)) }
      if(bound_time[i,2]!=Inf){
        E.Cijt[i,] = ((1/(1-exp(-sum.lambda.exp.between[i,])))*exp.term[i]*lambda_old*jumptidx_in_LR[,i])
      }
    }
    
    ### Score
    sum.exp.term = c(outer(order_bound_times, Rij.star, '<=') %*% matrix(exp.term, ncol=1))
    
    Score_vec = NULL
    for(j in 1:p){
      score.step1 = sapply(1:length(order_bound_times), function(y){ sum((order_bound_times[y]<=Rij.star) * X[,j] * exp.term) })
      score = sum(sapply(1:n, function(x){
        sum((order_bound_times<=Rij.star[x]) * X[x,j] * E.Cijt[x,]) - sum((order_bound_times<=Rij.star[x]) * E.Cijt[x,] * score.step1/sum.exp.term)
      }))
      Score_vec = rbind(Score_vec, score)
    }
    
    ### Hessian
    sum.X.exp.term = NULL
    for(j in 1:p){
      sum.X.exp.term = rbind(sum.X.exp.term, sapply(1:length(order_bound_times), function(x){ sum(((order_bound_times[x]<=Rij.star)+0)*X[,j]*exp.term) }))
    }
    
    sum.XX.exp.term = NULL
    for(j in 1:p){
      for(k in j:p){
        sum.XX.exp.term = rbind(sum.XX.exp.term, sapply(1:length(order_bound_times), function(x){ sum(((order_bound_times[x]<=Rij.star)+0)*X[,j]*X[,k]*exp.term) }))
      }
    }
    
    Hessian = matrix(0, nrow=p, ncol=p)
    # X single columns
    idx = 0
    for(i in 1:p){
      for(j in i:p){
        idx = idx + 1
        Hessian[i,j] = sum(sapply(1:n, function(x){
          sum(E.Cijt[x,]*((sum.XX.exp.term[idx,]*sum.exp.term - sum.X.exp.term[i,]*sum.X.exp.term[j,]) / (sum.exp.term)^2)*(order_bound_times<=Rij.star[x]))
        }))
      }
    }
    
    #Hessian_mat = Hessian
    Hessian_mat = Hessian + t(Hessian)
    diag(Hessian_mat) = 0.5*diag(Hessian_mat)
    
    ### Update parameters
    #Params_new = Params_old + c(solve(Hessian_mat) %*% Score_vec)
    Params_new = Params_old + c(ginv(Hessian_mat) %*% Score_vec)
    
    ####################################
    ###         Lambda Update
    ####################################
    gamma.l = as.matrix(Params_new)
    exp.term.l = c(exp(X %*% gamma.l))       # length n array
    
    ### Lambda (jump size) Update
    lambda.new.num = sapply(1:length(order_bound_times), function(x){ sum(E.Cijt[,x] * (order_bound_times[x] <= Rij.star)) })
    lambda.new.denom = sapply(1:length(order_bound_times), function(x){ sum(exp.term.l * (order_bound_times[x] <= Rij.star)) })
    lambda_new = lambda.new.num / lambda.new.denom
    
    ################################################
    ###            Stopping Criteria
    ################################################
    error.param = max(abs(c(Params_new - Params_old)))
    error.lambda = max(abs(c(lambda_new - lambda_old)))
    if(max(c(error.param, error.lambda))<=0.0001){ break }
    #if(max(error.param) <= 0.0005){ break }
    
    # print(sprintf("iter# %d", iter))
    # print(sprintf("Param change: %0.6f",error.param))
    # print(sprintf("lambda change: %0.6f", error.lambda))
    # print(sprintf("loglik: %0.6f", loglik))
    # newpar = as.data.frame(matrix(Params_new, nrow=1))
    # colnames(newpar) = c(paste0('gamma',c(1:p)))
    # print(newpar)
    
  } # end of while loop #
  
  ###################### Output results ######################
  result_all = list()
  result_all[['iter']] = iter
  result_all[['param']] = Params_new
  result_all[['lambda_est']] = lambda_new
  result_all[['order_bt']] = order_bound_times
  result_all[['loglik_vec']] = loglik_i
  result_all[['lcr']] = left.censor.rate
  result_all[['rcr']] = right.censor.rate
  return(result_all)
  
} #end of DentalEM_r0_adj function#

###################### Function for semi-param transformation model fit (w/ fixed gamma) ######################
dentalEM_fixgamma_r0_adj = function(X, gamma_est, bound_time, order_bound_times, r_alpha, maxiter){
  #-----------------------------------------------------------
  # X: covariates
  # bound_time: left, right bound times
  # order_bound_times: ordered baseline hazard jump time pts
  # r_alpha: semi-parametric transformation parameter
  # maxiter: maximum numer of iter
  # gamma_est: gamma estimates
  #-----------------------------------------------------------
  n = nrow(X)
  p = ncol(X)
  ###################### Initial Values ######################
  L = length(order_bound_times)
  lambda.ini = rep(1/L, L)
  
  lambda_old = rep(0, length(order_bound_times))
  lambda_new = lambda.ini
  
  ###################### Apply EM Algorithm ######################
  iter = 0
  while (iter < maxiter){
    
    iter = iter + 1
    lambda_old = lambda_new
    gamma = as.matrix(gamma_est)
    
    Rij.star = apply(bound_time,1,function(x){ max(x[x<Inf]) })
    
    exp.term = c(exp(X %*% gamma))       # length n array
    
    ### index of jump time points between LR bounds
    jumptidx_in_LR = (outer(order_bound_times,bound_time[,1],'>') & outer(order_bound_times,bound_time[,2],'<='))+0    # len(order_bound_time) by n matrix
    
    sum.lambda.exp.left = rowSums(exp.term * t(outer(order_bound_times,bound_time[,1],'<=')*lambda_old))       #length n array
    sum.lambda.exp.right = rowSums(exp.term * t(outer(order_bound_times,bound_time[,2],'<=')*lambda_old))      #length n array
    sum.lambda.exp.between = matrix(colSums(lambda_old*jumptidx_in_LR)*exp.term, ncol=1)     #n by n.l matrix
    
    ### logLikelihood
    lik_i = rep(NA,n)
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) - exp(-Gfunc(sum.lambda.exp.right[i], r_alpha)) }
      else if(bound_time[i,2] == Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) }
    }
    loglik_i = log(lik_i)
    
    ### E[C_ijt]
    E.Cijt = matrix(nrow=n, ncol=length(order_bound_times))      # n by len(order_bound_time) matrix
    for(i in 1:n){
      if(bound_time[i,2]==Inf){ E.Cijt[i,] = rep(0, length(order_bound_times)) }
      if(bound_time[i,2]!=Inf){
        E.Cijt[i,] = ((1/(1-exp(-sum.lambda.exp.between[i,])))*exp.term[i]*lambda_old*jumptidx_in_LR[,i])
      }
    }
    
    ###################### Lambda Update ######################
    ### Lambda (jump size) Update
    lambda.new.num = sapply(1:length(order_bound_times), function(x){ sum(E.Cijt[,x] * (order_bound_times[x] <= Rij.star)) })
    lambda.new.denom = sapply(1:length(order_bound_times), function(x){ sum(exp.term * (order_bound_times[x] <= Rij.star)) })
    lambda_new = lambda.new.num / lambda.new.denom
    
    ###################### Stopping Criteria ######################
    error.lambda = max(abs(c(lambda_new - lambda_old)))
    if(max(error.lambda)<=0.0001){ break }
    
  } # end of while loop #
  
  return(loglik_i)
  
} #end of DentalEM_r0_adj (fixed gamma) function#

###################### Main function ######################
semipar_trans_fit_r0_adj = function(X, bound_time, perturb, r_alpha, maxiter, cov_method){
  #-----------------------------------------------------------
  # X: covariate matrix
  # bound_time: left and right bound time
  # perturb: integer 1, 5, 10, for perturbation
  # r_alpha: semi-parametrix transformation parameter
  # maxiter: maximum EM iteration allowed
  # cov_method = 0: do not estimate covariance_mat
  #            = integer: boostrap method for covariance_mat
  #            = 'ZL1': 1st-order (Zeng 2017)
  #            = 'ZL2': 2nd-order (Zeng 2016)
  #-----------------------------------------------------------
  order_bound_times = sort(unique(c(bound_time[c(bound_time)!=0 & c(bound_time)!=Inf])))
  n = nrow(X)
  p = ncol(X)
  ###################### Fit semi-param transformation model ######################
  fit_bl = DentalEM_r0_adj(X, bound_time, order_bound_times, r_alpha, maxiter)
  
  param_bl = fit_bl$param          #parameter estimate
  loglik_bl = fit_bl$loglik_vec    #logliklihood vector (each represent one subject)
  lambda_est = fit_bl$lambda_est   #baseline hazard estimate
  order_bt = fit_bl$order_bt       #ordered baseline hazard jump time
  #lcr = fit_bl$lcr                 #left censor rate
  #rcr = fit_bl$rcr                 #right censor rate
  #iternum = fit_bl$iter            #iteration before algorithm stops
  
  ###################### Zeng 2017 method --> covariance matrix estimate ######################
  if(cov_method == 'ZL1'){
    #Add perturbation to "param_bl"
    delta = perturb/sqrt(n)
    param_delta = matrix(rep(param_bl,p),byrow=T,nrow=p) + diag(delta,p,p)
    
    #Gradient of logliklihood
    loglik_grad = NULL
    loglik_0 = dentalEM_fixgamma_r0_adj(X, param_bl, bound_time, order_bound_times, r_alpha, maxiter)
    for(i in 1:nrow(param_delta)){
      loglik_delta = dentalEM_fixgamma_r0_adj(X, param_delta[i,], bound_time, order_bound_times, r_alpha, maxiter)
      loglik_grad = rbind(loglik_grad, (loglik_delta-loglik_0)/delta)
    }
    
    # V matrix
    V_mat = 0
    for(j in 1:n){ V_mat = V_mat + tcrossprod(as.matrix(loglik_grad[,j])) }
    
    # covariance matrix estimate
    covest = solve(V_mat)
    
  } #end of if(cov_method=='ZL1')#
  ###################### Zeng 2016 method --> covariance matrix estimate ######################
  if(cov_method == 'ZL2'){
    #Add perturbation to "param_bl"
    delta = perturb/sqrt(n)
    
    #Gradient of logliklihood
    V_mat = matrix(NA, nrow=p, nrow=p)
    for(i in 1:p){
      for(j in i:p){
        param_delta_i = param_bl + delta*(c(1:p)==i)
        param_delta_j = param_bl + delta*(c(1:p)==j)
        param_delta_ij = param_bl + delta*(c(1:p)==i) + delta*(c(1:p)==j)
        V_mat[i,j] = sum(dentalEM_fixgamma_r0_adj(X, param_bl, bound_time, order_bound_times, n.l, r_alpha, maxiter)) - 
          sum(dentalEM_fixgamma_r0_adj(X, param_delta_i, bound_time, order_bound_times, n.l, r_alpha, maxiter)) - 
          sum(dentalEM_fixgamma_r0_adj(X, param_delta_j, bound_time, order_bound_times, n.l, r_alpha, maxiter)) +
          sum(dentalEM_fixgamma_r0_adj(X, param_delta_ij, bound_time, order_bound_times, n.l, r_alpha, maxiter))
      }
    }
    V_mat = V_mat + t(V_mat)
    diag(V_mat) = 0.5*diag(V_mat)
    
    # covariance matrix estimate
    covest = -1*solve(V_mat)
    
  } #end of if(cov_method=='ZL2')#
  ###################### Do not need covariance matrix estimate ######################
  if(cov_method == 0){ covest = NULL }
  
  ###################### Output results ######################
  output = list()
  output[['b']] = c(param_bl)
  output[['covest']] = covest
  output[['order_bt']] = c(order_bt)
  output[['lambda_est']] = c(lambda_est)
  output[['loglik']] = sum(loglik_bl)
  return(output)
  
} #end of output function#

Gfunc = function(x, r){
  if(r==0){ x }
  else if(r>0){ log(1+r*x)/r }
}

###################### Function for semi-param transformation model fit ######################
DentalEM_r0_noadj = function(bound_time, order_bound_times, r_alpha, maxiter){
  #-----------------------------------------------------------#
  # bound_time: left, right bound times
  # order_bound_times: ordered baseline hazard jump time pts
  # r_alpha: semi-parametric transformation parameter
  # maxiter: maximum numer of iter
  #-----------------------------------------------------------#
  n = nrow(bound_time)
  # Censor rate
  left.censor.rate = sum(bound_time[,1]==0)/n
  right.censor.rate = sum(bound_time[,2]==Inf)/n
  #---------------------- Initial Values -----------------------#
  L = length(order_bound_times)
  lambda.ini = rep(1/L, L)
  
  lambda_old = rep(0, length(order_bound_times))
  lambda_new = lambda.ini
  
  #-------------------- Apply EM Algorithm --------------------#
  ### Initial values
  iter = 0
  
  ### collect all loglik values
  loglik_all = NULL
  
  message(" Now: Obtaining the unpenalized nonparametric MLE")
  ### Start iteration
  while (iter < maxiter){
    
    iter = iter + 1
    
    lambda_old = lambda_new
    
    # plot(order_bound_times, cumsum(lambda_new))
    # lines(order_bound_times, order_bound_times/2)
    
    Rij.star = apply(bound_time,1,function(x){ max(x[x<Inf]) })
    
    exp.term = rep(1,n)       # length n array
    
    ### index of jump time points between LR bounds
    jumptidx_in_LR = (outer(order_bound_times,bound_time[,1],'>') & outer(order_bound_times,bound_time[,2],'<='))+0    # len(order_bound_time) by n matrix
    
    sum.lambda.exp.left = rowSums(exp.term * t(outer(order_bound_times,bound_time[,1],'<=')*lambda_old))       #length n array
    sum.lambda.exp.right = rowSums(exp.term * t(outer(order_bound_times,bound_time[,2],'<=')*lambda_old))      #length n array
    sum.lambda.exp.between = matrix(colSums(lambda_old*jumptidx_in_LR)*exp.term, ncol=1)    #n by 1 matrix
    
    ### logLikelihood
    lik_i = rep(NA,n)
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) - exp(-Gfunc(sum.lambda.exp.right[i], r_alpha)) }
      else if(bound_time[i,2] == Inf){ lik_i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) }
    }
    loglik_i = log(lik_i)
    loglik = sum(loglik_i)
    loglik_all = c(loglik_all, loglik)
    
    ### E[C_ijt]
    E.Cijt = matrix(nrow=n, ncol=length(order_bound_times))      # n by len(order_bound_time) matrix
    for(i in 1:n){
      if(bound_time[i,2]==Inf){ E.Cijt[i,] = rep(0, length(order_bound_times)) }
      if(bound_time[i,2]!=Inf){
        E.Cijt[i,] = ((1/(1-exp(-sum.lambda.exp.between[i,])))*exp.term[i]*lambda_old*jumptidx_in_LR[,i])
      }
    }
    
    ####################################
    ###         Lambda Update
    ####################################
    ### Lambda (jump size) Update
    lambda.new.num = sapply(1:length(order_bound_times), function(x){ sum(E.Cijt[,x] * (order_bound_times[x] <= Rij.star)) })
    lambda.new.denom = sapply(1:length(order_bound_times), function(x){ sum(exp.term * (order_bound_times[x] <= Rij.star)) })
    lambda_new = lambda.new.num / lambda.new.denom
    
    ################################################
    ###            Stopping Criteria
    ################################################
    error.lambda = max(abs(c(lambda_new - lambda_old)))
    if(max(error.lambda)<=0.0001){ break }
    
    # print(sprintf("iter# %d", iter))
    # print(sprintf("lambda change: %0.6f", error.lambda))
    # print(sprintf("loglik: %0.6f", loglik))
    
  } # end of while loop #
  
  ###################### Output results ######################
  result_all = list()
  result_all[['iter']] = iter
  result_all[['lambda_est']] = lambda_new
  result_all[['order_bt']] = order_bound_times
  result_all[['loglik_vec']] = loglik_i
  result_all[['lcr']] = left.censor.rate
  result_all[['rcr']] = right.censor.rate
  return(result_all)
  
} #end of DentalEM_r0 function#

###################### Main function ######################
semipar_trans_fit_r0_noadj = function(bound_time, r_alpha, maxiter){
  #-----------------------------------------------------------
  # bound_time: left and right bound time
  # r_alpha: semi-parametrix transformation parameter
  # maxiter: maximum EM iteration allowed
  #-----------------------------------------------------------
  order_bound_times = sort(unique(c(bound_time[c(bound_time)!=0 & c(bound_time)!=Inf])))
  ###################### Fit semi-param transformation model ######################
  fit_bl = DentalEM_r0_noadj(bound_time, order_bound_times, r_alpha, maxiter)
  
  loglik_bl = fit_bl$loglik_vec    #logliklihood vector (each represent one subject)
  lambda_est = fit_bl$lambda_est   #baseline hazard estimate
  order_bt = fit_bl$order_bt       #ordered baseline hazard jump time
  #lcr = fit_bl$lcr                 #left censor rate
  #rcr = fit_bl$rcr                 #right censor rate
  #iternum = fit_bl$iter            #iteration before algorithm stops
  
  ###################### Output results ######################
  output = list()
  output[['order_bt']] = c(order_bt)
  output[['lambda_est']] = c(lambda_est)
  output[['loglik']] = sum(loglik_bl)
  return(output)
} #end of output function#

Gfunc = function(x, r){
  if(r==0){ x }
  else if(r>0){ log(1+r*x)/r }
}

###################### Function for semi-param transformation model fit ######################
DentalEM_rp_adj = function(X, bound_time, order_bound_times, n.l, r_alpha, maxiter){
  #-----------------------------------------------------------#
  # X: covariates
  # bound_time: left, right bound times
  # order_bound_times: ordered baseline hazard jump time pts
  # n.l: number of Gauss-Laguerre quadrature points
  # r_alpha: semi-parametric transformation parameter
  # maxiter: maximum numer of iter
  #-----------------------------------------------------------#
  n = nrow(X)
  p = ncol(X)
  # Censor rate
  left.censor.rate = sum(bound_time[,1]==0)/n
  right.censor.rate = sum(bound_time[,2]==Inf)/n
  #---------------------- Initial Values -----------------------#
  L = length(order_bound_times)
  lambda.ini = rep(1/L, L)
  
  gamma.ini = rep(0, p)
  
  Params_old = rep(0, p)
  Params_new = gamma.ini
  
  lambda_old = rep(0, length(order_bound_times))
  lambda_new = lambda.ini
  
  #-------------------- Gaussian Quadrature --------------------#
  #Gauss-Laguerre: 1-D --> Alpha
  GLQ = gauss.quad(n.l, kind="laguerre")
  node_l = as.matrix(GLQ$nodes)
  weight_l = as.matrix(GLQ$weights)
  
  #----------------------- Density Alpha ----------------------#
  #Density for alpha; not including exp{-x} part
  density_Alpha = function(x, r){ (1/r)^(1/r) * x^(1/r-1) / factorial(1/r-1) }
  
  #-------------------- Apply EM Algorithm --------------------#
  ### Initial values
  iter = 0
  
  ### collect all loglik values
  loglik_all = NULL
  
  message(" Now: Obtaining the unpenalized nonparametric MLE")
  ### Start iteration
  while (iter < maxiter){
    
    iter = iter + 1
    
    Params_old = Params_new
    lambda_old = lambda_new
    
    # plot(order_bound_times, cumsum(lambda_new))
    # lines(order_bound_times, order_bound_times/2)
    
    gamma = as.matrix(Params_old)
    
    Rij.star = apply(bound_time,1,function(x){ max(x[x<Inf]) })
    
    gamma.X = c(X %*% gamma)
    alpha.all = r_alpha*node_l[,1]
    density.alpha.all = density_Alpha(alpha.all,r_alpha)
    
    exp.term = c(exp(gamma.X))       # length n array
    
    ### index of jump time points between LR bounds
    jumptidx_in_LR = (outer(order_bound_times,bound_time[,1],'>') & outer(order_bound_times,bound_time[,2],'<='))+0    # len(order_bound_time) by n matrix
    
    sum.lambda.exp.left = rowSums(exp.term * t(outer(order_bound_times,bound_time[,1],'<=')*lambda_old))       #length n array
    sum.lambda.exp.right = rowSums(exp.term * t(outer(order_bound_times,bound_time[,2],'<=')*lambda_old))      #length n array
    sum.lambda.alpha.exp.between = matrix(colSums(lambda_old*jumptidx_in_LR)*exp.term, ncol=1) %*% matrix(alpha.all,nrow=1)     #n by n.l matrix
    
    density.AB.rfs.i = matrix(nrow=n, ncol=n.l)    # n x n.l matrix
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ density.AB.rfs.i[i,] = exp(-sum.lambda.exp.left[i]*alpha.all) - exp(-sum.lambda.exp.right[i]*alpha.all) }
      else if(bound_time[i,2] == Inf){ density.AB.rfs.i[i,] = exp(-sum.lambda.exp.left[i]*alpha.all) }
    }
    
    density.all = t(density.alpha.all * t(density.AB.rfs.i))
    
    int.all.i = rep(NA, n)  # length n array
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ int.all.i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) - exp(-Gfunc(sum.lambda.exp.right[i], r_alpha)) }
      else if(bound_time[i,2] == Inf){ int.all.i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) }
    }
    
    ### logLikelihood
    loglik_i = log(int.all.i)
    loglik = sum(loglik_i)
    loglik_all = c(loglik_all, loglik)
    
    ### E[alpha * exp{}]
    E.alpha.exp = c(((matrix(exp.term,ncol=1) %*% matrix(alpha.all,nrow=1)) * density.all) %*% weight_l) / int.all.i  # length n array
    
    ### E[C_ijt]
    E.Cijt = matrix(nrow=n, ncol=length(order_bound_times))      # n by len(order_bound_time) matrix
    for(i in 1:n){
      if(bound_time[i,2]==Inf){ E.Cijt[i,] = rep(0, length(order_bound_times)) }
      if(bound_time[i,2]!=Inf){
        E.Cijt[i,] = c(t(weight_l*density.all[i,]) %*% ((1/(1-exp(-sum.lambda.alpha.exp.between[i,]))) * (matrix(alpha.all,ncol=1) %*% matrix(exp.term[i]*lambda_old*jumptidx_in_LR[,i],nrow=1)))) / int.all.i[i]
      }
    }
    
    ### Score
    sum.E.alpha.exp = c(outer(order_bound_times, Rij.star, '<=') %*% matrix(E.alpha.exp, ncol=1))
    
    Score_vec = NULL
    for(j in 1:p){
      score.step1 = sapply(1:length(order_bound_times), function(y){ sum((order_bound_times[y]<=Rij.star) * X[,j] * E.alpha.exp) })
      score = sum(sapply(1:n, function(x){
        sum((order_bound_times<=Rij.star[x]) * X[x,j] * E.Cijt[x,]) - sum((order_bound_times<=Rij.star[x]) * E.Cijt[x,] * score.step1/sum.E.alpha.exp)
      }))
      Score_vec = rbind(Score_vec, score)
    }
    
    ### Hessian
    sum.X.E.alpha.exp = NULL
    for(j in 1:p){
      sum.X.E.alpha.exp = rbind(sum.X.E.alpha.exp, sapply(1:length(order_bound_times), function(x){ sum(((order_bound_times[x]<=Rij.star)+0)*X[,j]*E.alpha.exp) }))
    }
    
    sum.XX.E.alpha.exp = NULL
    for(j in 1:p){
      for(k in j:p){
        sum.XX.E.alpha.exp = rbind(sum.XX.E.alpha.exp, sapply(1:length(order_bound_times), function(x){ sum(((order_bound_times[x]<=Rij.star)+0)*X[,j]*X[,k]*E.alpha.exp) }))
      }
    }
    
    Hessian = matrix(0, nrow=p, ncol=p)
    # X single columns
    idx = 0
    for(i in 1:p){
      for(j in i:p){
        idx = idx + 1
        Hessian[i,j] = sum(sapply(1:n, function(x){
          sum(E.Cijt[x,]*((sum.XX.E.alpha.exp[idx,]*sum.E.alpha.exp - sum.X.E.alpha.exp[i,]*sum.X.E.alpha.exp[j,]) / (sum.E.alpha.exp)^2)*(order_bound_times<=Rij.star[x]))
        }))
      }
    }
    
    #Hessian_mat = Hessian
    Hessian_mat = Hessian + t(Hessian)
    diag(Hessian_mat) = 0.5*diag(Hessian_mat)
    
    ### Update parameters
    #Params_new = Params_old + c(solve(Hessian_mat) %*% Score_vec)
    Params_new = Params_old + c(ginv(Hessian_mat) %*% Score_vec)
    
    ####################################
    ###         Lambda Update
    ####################################
    gamma.l = as.matrix(Params_new)
    gammaX.l = c(X %*% gamma.l)
    exp.term.l = c(exp(gammaX.l))       # length n array
    
    ### E[alpha * exp{}]
    E.alpha.exp.l = c(((matrix(exp.term.l,ncol=1) %*% matrix(alpha.all,nrow=1)) * density.all) %*% weight_l) / int.all.i  # length n array
    
    ### Lambda (jump size) Update
    lambda.new.num = sapply(1:length(order_bound_times), function(x){ sum(E.Cijt[,x] * (order_bound_times[x] <= Rij.star)) })
    lambda.new.denom = sapply(1:length(order_bound_times), function(x){ sum(E.alpha.exp.l * (order_bound_times[x] <= Rij.star)) })
    lambda_new = lambda.new.num / lambda.new.denom
    
    ################################################
    ###            Stopping Criteria
    ################################################
    error.param = max(abs(c(Params_new - Params_old)))
    error.lambda = max(abs(c(lambda_new - lambda_old)))
    if(max(c(error.param, error.lambda))<=0.0001){ break }
    #if(max(error.param) <= 0.0005){ break }
    
    # print(sprintf("iter# %d", iter))
    # print(sprintf("Param change: %0.6f",error.param))
    # print(sprintf("lambda change: %0.6f", error.lambda))
    # print(sprintf("loglik: %0.6f", loglik))
    # newpar = as.data.frame(matrix(Params_new, nrow=1))
    # colnames(newpar) = c(paste0('gamma',c(1:p)))
    # print(newpar)
    
  } # end of while loop #
  
  ###################### Output results ######################
  result_all = list()
  result_all[['iter']] = iter
  result_all[['param']] = Params_new
  result_all[['lambda_est']] = lambda_new
  result_all[['order_bt']] = order_bound_times
  result_all[['loglik_vec']] = loglik_i 
  result_all[['lcr']] = left.censor.rate
  result_all[['rcr']] = right.censor.rate
  return(result_all)
  
} #end of DentalEM_rp_adj function#

###################### Function for semi-param transformation model fit (w/ fixed gamma) ######################
DentalEM_fixgamma_rp_adj = function(X, gamma_est, bound_time, order_bound_times, n.l, r_alpha, maxiter){
  #-----------------------------------------------------------
  # X: covariates
  # bound_time: left, right bound times
  # order_bound_times: ordered baseline hazard jump time pts
  # n.l: number of Gauss-Laguerre quadrature points
  # r_alpha: semi-parametric transformation parameter
  # maxiter: maximum numer of iter
  # gamma_est: gamma estimates
  #-----------------------------------------------------------
  Gfunc = function(x, r){
    if(r==0){ x }
    else if(r>0){ log(1+r*x)/r }
  }
  n = nrow(X)
  p = ncol(X)
  ###################### Initial Values ######################
  L = length(order_bound_times)
  lambda.ini = rep(1/L, L)
  
  lambda_old = rep(0, length(order_bound_times))
  lambda_new = lambda.ini
  
  ###################### Gaussian Quadrature ######################
  #Gauss-Laguerre: 1-D --> Alpha
  GLQ = gauss.quad(n.l, kind="laguerre")
  node_l = as.matrix(GLQ$nodes)
  weight_l = as.matrix(GLQ$weights)
  
  ###################### Density Alpha ######################
  #Density for alpha; not including exp{-x} part
  density_Alpha = function(x, r){ (1/r)^(1/r) * x^(1/r-1) / factorial(1/r-1) }
  
  ###################### Apply EM Algorithm ######################
  iter = 0
  while (iter < maxiter){
    
    iter = iter + 1
    lambda_old = lambda_new
    gamma = as.matrix(gamma_est)
    
    Rij.star = apply(bound_time,1,function(x){ max(x[x<Inf]) })
    
    gamma.X = c(X %*% gamma)
    alpha.all = r_alpha*node_l[,1]
    density.alpha.all = density_Alpha(alpha.all,r_alpha)
    
    exp.term = c(exp(gamma.X))       # length n array
    
    ### index of jump time points between LR bounds
    jumptidx_in_LR = (outer(order_bound_times,bound_time[,1],'>') & outer(order_bound_times,bound_time[,2],'<='))+0    # len(order_bound_time) by n matrix
    
    sum.lambda.exp.left = rowSums(exp.term * t(outer(order_bound_times,bound_time[,1],'<=')*lambda_old))       #length n array
    sum.lambda.exp.right = rowSums(exp.term * t(outer(order_bound_times,bound_time[,2],'<=')*lambda_old))      #length n array
    sum.lambda.alpha.exp.between = matrix(colSums(lambda_old*jumptidx_in_LR)*exp.term, ncol=1) %*% matrix(alpha.all,nrow=1)     #n by n.l matrix
    
    density.AB.rfs.i = matrix(nrow=n, ncol=n.l)    # n x n.l matrix
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ density.AB.rfs.i[i,] = exp(-sum.lambda.exp.left[i]*alpha.all) - exp(-sum.lambda.exp.right[i]*alpha.all) }
      else if(bound_time[i,2] == Inf){ density.AB.rfs.i[i,] = exp(-sum.lambda.exp.left[i]*alpha.all) }
    }
    
    density.all = t(density.alpha.all * t(density.AB.rfs.i))
    
    int.all.i = rep(NA, n)  # length n array
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ int.all.i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) - exp(-Gfunc(sum.lambda.exp.right[i], r_alpha)) }
      else if(bound_time[i,2] == Inf){ int.all.i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) }
    }
    
    ### logLikelihood
    loglik_i = log(int.all.i)
    #loglik = sum(loglik_i)
    
    ### E[alpha * exp{}]
    E.alpha.exp = c(((matrix(exp.term,ncol=1) %*% matrix(alpha.all,nrow=1)) * density.all) %*% weight_l) / int.all.i  # length n array
    
    ### E[C_ijt]
    E.Cijt = matrix(nrow=n, ncol=length(order_bound_times))      # n by len(order_bound_time) matrix
    for(i in 1:n){
      if(bound_time[i,2]==Inf){ E.Cijt[i,] = rep(0, length(order_bound_times)) }
      if(bound_time[i,2]!=Inf){
        E.Cijt[i,] = c(t(weight_l*density.all[i,]) %*% ((1/(1-exp(-sum.lambda.alpha.exp.between[i,]))) * (matrix(alpha.all,ncol=1) %*% matrix(exp.term[i]*lambda_old*jumptidx_in_LR[,i],nrow=1)))) / int.all.i[i]
      }
    }
    
    ###################### Lambda Update ######################
    ### Lambda (jump size) Update
    lambda.new.num = sapply(1:length(order_bound_times), function(x){ sum(E.Cijt[,x] * (order_bound_times[x] <= Rij.star)) })
    lambda.new.denom = sapply(1:length(order_bound_times), function(x){ sum(E.alpha.exp * (order_bound_times[x] <= Rij.star)) })
    lambda_new = lambda.new.num / lambda.new.denom
    
    ###################### Stopping Criteria ######################
    error.lambda = max(abs(c(lambda_new - lambda_old)))
    if(max(error.lambda)<=0.0001){ break }
    
  } # end of while loop #
  
  return(loglik_i)
  
} #end of DentalEM_rp_adj (fixed gamma) function#

###################### Main function ######################
semipar_trans_fit_rp_adj = function(X, bound_time, perturb, n.l, r_alpha, maxiter, cov_method){
  #-----------------------------------------------------------
  # X: covariate matrix
  # bound_time: left and right bound time
  # perturb: integer 1, 5, 10, for perturbation
  # n.l: Gauss Laguerre quadrature point number
  # r_alpha: semi-parametrix transformation parameter
  # maxiter: maximum EM iteration allowed
  # cov_method = 0: do not estimate covariance_mat
  #            = integer: boostrap method for covariance_mat
  #            = 'ZL1': 1st-order (Zeng 2017)
  #            = 'ZL2': 2nd-order (Zeng 2016)
  #-----------------------------------------------------------
  order_bound_times = sort(unique(c(bound_time[c(bound_time)!=0 & c(bound_time)!=Inf])))
  n = nrow(X)
  p = ncol(X)
  ###################### Fit semi-param transformation model ######################
  fit_bl = DentalEM_rp_adj(X, bound_time, order_bound_times, n.l, r_alpha, maxiter)
  
  param_bl = fit_bl$param          #parameter estimate
  loglik_bl = fit_bl$loglik_vec    #logliklihood vector (each represent one subject)
  lambda_est = fit_bl$lambda_est   #baseline hazard estimate
  order_bt = fit_bl$order_bt       #ordered baseline hazard jump time
  #lcr = fit_bl$lcr                 #left censor rate
  #rcr = fit_bl$rcr                 #right censor rate
  #iternum = fit_bl$iter            #iteration before algorithm stops
  
  ###################### Zeng 2017 method --> covariance matrix estimate ######################
  if(cov_method == 'ZL1'){
    message(" Now: calculating the covariance matrix")
    #Add perturbation to "param_bl"
    delta = perturb/sqrt(n)
    param_delta = matrix(rep(param_bl,p),byrow=T,nrow=p) + diag(delta,p,p)
    
    #Gradient of logliklihood
    loglik_grad = NULL
    loglik_0 = DentalEM_fixgamma_rp_adj(X, param_bl, bound_time, order_bound_times, n.l, r_alpha, maxiter)
    for(i in 1:nrow(param_delta)){
      loglik_delta = DentalEM_fixgamma_rp_adj(X, param_delta[i,], bound_time, order_bound_times, n.l, r_alpha, maxiter)
      loglik_grad = rbind(loglik_grad, (loglik_delta-loglik_0)/delta)
    }
    
    # V matrix
    V_mat = 0
    for(j in 1:n){ V_mat = V_mat + tcrossprod(as.matrix(loglik_grad[,j])) }
    
    # covariance matrix estimate
    covest = solve(V_mat)
    
  } #end of if(cov_method=='ZL1')#
  ###################### Zeng 2016 method --> covariance matrix estimate ######################
  if(cov_method == 'ZL2'){
    #Add perturbation to "param_bl"
    delta = perturb/sqrt(n)
    
    #Gradient of logliklihood
    V_mat = matrix(NA, nrow=p, nrow=p)
    for(i in 1:p){
      for(j in i:p){
        param_delta_i = param_bl + delta*(c(1:p)==i)
        param_delta_j = param_bl + delta*(c(1:p)==j)
        param_delta_ij = param_bl + delta*(c(1:p)==i) + delta*(c(1:p)==j)
        V_mat[i,j] = sum(DentalEM_fixgamma_rp_adj(X, param_bl, bound_time, order_bound_times, n.l, r_alpha, maxiter)) - 
          sum(DentalEM_fixgamma_rp_adj(X, param_delta_i, bound_time, order_bound_times, n.l, r_alpha, maxiter)) - 
          sum(DentalEM_fixgamma_rp_adj(X, param_delta_j, bound_time, order_bound_times, n.l, r_alpha, maxiter)) +
          sum(DentalEM_fixgamma_rp_adj(X, param_delta_ij, bound_time, order_bound_times, n.l, r_alpha, maxiter))
      }
    }
    V_mat = V_mat + t(V_mat)
    diag(V_mat) = 0.5*diag(V_mat)
    
    # covariance matrix estimate
    covest = -1*solve(V_mat)
    
  } #end of if(cov_method=='ZL2')#
  ###################### Do not need covariance matrix estimate ######################
  if(cov_method == 0){ covest = NULL }
  
  ###################### Output results ######################
  output = list()
  output[['b']] = c(param_bl)
  output[['covest']] = covest
  output[['order_bt']] = c(order_bt)
  output[['lambda_est']] = c(lambda_est)
  output[['loglik']] = sum(loglik_bl)
  return(output)
} #end of output function#

Gfunc = function(x, r){
  if(r==0){ x }
  else if(r>0){ log(1+r*x)/r }
}

###################### Function for semi-param transformation model fit ######################
DentalEM_rp_noadj = function(bound_time, order_bound_times, n.l, r_alpha, maxiter){
  #-----------------------------------------------------------#
  # bound_time: left, right bound times
  # order_bound_times: ordered baseline hazard jump time pts
  # n.l: number of Gauss-Laguerre quadrature points
  # r_alpha: semi-parametric transformation parameter
  # maxiter: maximum numer of iter
  #-----------------------------------------------------------#
  n = nrow(bound_time)
  # Censor rate
  left.censor.rate = sum(bound_time[,1]==0)/n
  right.censor.rate = sum(bound_time[,2]==Inf)/n
  #---------------------- Initial Values -----------------------#
  L = length(order_bound_times)
  lambda.ini = rep(1/L, L)
  
  lambda_old = rep(0, length(order_bound_times))
  lambda_new = lambda.ini
  
  #-------------------- Gaussian Quadrature --------------------#
  #Gauss-Laguerre: 1-D --> Alpha
  GLQ = gauss.quad(n.l, kind="laguerre")
  node_l = as.matrix(GLQ$nodes)
  weight_l = as.matrix(GLQ$weights)
  
  #----------------------- Density Alpha ----------------------#
  #Density for alpha; not including exp{-x} part
  density_Alpha = function(x, r){ (1/r)^(1/r) * x^(1/r-1) / factorial(1/r-1) }
  
  #-------------------- Apply EM Algorithm --------------------#
  ### Initial values
  iter = 0
  
  ### collect all loglik values
  loglik_all = NULL
  
  message(" Now: Obtaining the unpenalized nonparametric MLE")
  ### Start iteration
  while (iter < maxiter){
    
    iter = iter + 1
    
    lambda_old = lambda_new
    
    # plot(order_bound_times, cumsum(lambda_new))
    # lines(order_bound_times, order_bound_times/2)
    
    Rij.star = apply(bound_time,1,function(x){ max(x[x<Inf]) })
    
    alpha.all = r_alpha*node_l[,1]
    density.alpha.all = density_Alpha(alpha.all,r_alpha)
    
    exp.term = rep(1,n)       # length n array
    
    ### index of jump time points between LR bounds
    jumptidx_in_LR = (outer(order_bound_times,bound_time[,1],'>') & outer(order_bound_times,bound_time[,2],'<='))+0    # len(order_bound_time) by n matrix
    
    sum.lambda.exp.left = rowSums(exp.term * t(outer(order_bound_times,bound_time[,1],'<=')*lambda_old))       #length n array
    sum.lambda.exp.right = rowSums(exp.term * t(outer(order_bound_times,bound_time[,2],'<=')*lambda_old))      #length n array
    sum.lambda.alpha.exp.between = matrix(colSums(lambda_old*jumptidx_in_LR)*exp.term, ncol=1) %*% matrix(alpha.all,nrow=1)     #n by n.l matrix
    
    density.AB.rfs.i = matrix(nrow=n, ncol=n.l)    # n x n.l matrix
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ density.AB.rfs.i[i,] = exp(-sum.lambda.exp.left[i]*alpha.all) - exp(-sum.lambda.exp.right[i]*alpha.all) }
      else if(bound_time[i,2] == Inf){ density.AB.rfs.i[i,] = exp(-sum.lambda.exp.left[i]*alpha.all) }
    }
    
    density.all = t(density.alpha.all * t(density.AB.rfs.i))
    
    int.all.i = rep(NA, n)  # length n array
    for(i in 1:n){
      if(bound_time[i,2] != Inf){ int.all.i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) - exp(-Gfunc(sum.lambda.exp.right[i], r_alpha)) }
      else if(bound_time[i,2] == Inf){ int.all.i[i] = exp(-Gfunc(sum.lambda.exp.left[i], r_alpha)) }
    }
    
    ### logLikelihood
    loglik_i = log(int.all.i)
    loglik = sum(loglik_i)
    loglik_all = c(loglik_all, loglik)
    
    ### E[alpha * exp{}]
    E.alpha.exp = c(((matrix(exp.term,ncol=1) %*% matrix(alpha.all,nrow=1)) * density.all) %*% weight_l) / int.all.i  # length n array
    
    ### E[C_ijt]
    E.Cijt = matrix(nrow=n, ncol=length(order_bound_times))      # n by len(order_bound_time) matrix
    for(i in 1:n){
      if(bound_time[i,2]==Inf){ E.Cijt[i,] = rep(0, length(order_bound_times)) }
      if(bound_time[i,2]!=Inf){
        E.Cijt[i,] = c(t(weight_l*density.all[i,]) %*% ((1/(1-exp(-sum.lambda.alpha.exp.between[i,]))) * (matrix(alpha.all,ncol=1) %*% matrix(exp.term[i]*lambda_old*jumptidx_in_LR[,i],nrow=1)))) / int.all.i[i]
      }
    }
    
    ####################################
    ###         Lambda Update
    ####################################
    ### Lambda (jump size) Update
    lambda.new.num = sapply(1:length(order_bound_times), function(x){ sum(E.Cijt[,x] * (order_bound_times[x] <= Rij.star)) })
    lambda.new.denom = sapply(1:length(order_bound_times), function(x){ sum(E.alpha.exp * (order_bound_times[x] <= Rij.star)) })
    lambda_new = lambda.new.num / lambda.new.denom
    
    ################################################
    ###            Stopping Criteria
    ################################################
    error.lambda = max(abs(c(lambda_new - lambda_old)))
    if(max(error.lambda)<=0.0001){ break }
    
    # print(sprintf("iter# %d", iter))
    # print(sprintf("lambda change: %0.6f", error.lambda))
    # print(sprintf("loglik: %0.6f", loglik))
    
  } # end of while loop #
  
  ###################### Output results ######################
  result_all = list()
  result_all[['iter']] = iter
  result_all[['lambda_est']] = lambda_new
  result_all[['order_bt']] = order_bound_times
  result_all[['loglik_vec']] = loglik_i 
  result_all[['lcr']] = left.censor.rate
  result_all[['rcr']] = right.censor.rate
  return(result_all)
  
} #end of DentalEM_rp_noadj function#

###################### Main function ######################
semipar_trans_fit_rp_noadj = function(bound_time, n.l, r_alpha, maxiter){
  #-----------------------------------------------------------
  # bound_time: left and right bound time
  # n.l: Gauss Laguerre quadrature point number
  # r_alpha: semi-parametrix transformation parameter
  # maxiter: maximum EM iteration allowed
  #-----------------------------------------------------------
  order_bound_times = sort(unique(c(bound_time[c(bound_time)!=0 & c(bound_time)!=Inf])))
  ###################### Fit semi-param transformation model ######################
  fit_bl = DentalEM_rp_noadj(bound_time, order_bound_times, n.l, r_alpha, maxiter)
  
  loglik_bl = fit_bl$loglik_vec    #logliklihood vector (each represent one subject)
  lambda_est = fit_bl$lambda_est   #baseline hazard estimate
  order_bt = fit_bl$order_bt       #ordered baseline hazard jump time
  #lcr = fit_bl$lcr                 #left censor rate
  #rcr = fit_bl$rcr                 #right censor rate
  #iternum = fit_bl$iter            #iteration before algorithm stops
  
  ###################### Output results ######################
  output = list()
  output[['order_bt']] = c(order_bt)
  output[['lambda_est']] = c(lambda_est)
  output[['loglik']] = sum(loglik_bl)
  return(output)
} #end of output function#

unpencoxIC.default <- function(lowerIC, upperIC, X, trunc = NULL, normalize.X = TRUE, covmat = TRUE, cl = NULL, tol = 1e-3, niter = 1e5, string.cen = Inf, string.missing = NA, ...) {
  
  match.call()
  
  if(missing(trunc)) {
    trunc <- NULL
    ind.trunc <- FALSE
    smallest.trunc <- 0
  } else {
    ind.trunc <- TRUE
    smallest.trunc <- min(trunc)
  }
  
  if (!is.null(cl)) {
    if (.Platform$OS.type == "windows") {
      if (!inherits(cl, "cluster"))
        cl <- NULL
    } else {
      if (inherits(cl, "cluster")) {
        if (length(cl) < 2L)
          cl <- NULL
      } else {
        if (cl < 2)
          cl <- NULL
      }
    }
  }
  
  xnames <- colnames(X)
  
  arglist <- fun_arglist(lowerIC, upperIC, X, trunc, normalize.X, tol, niter)
  arglist$initial_lambda <- rep(1/nrow(arglist$set), nrow(arglist$set))
  
  message(" Now: Obtaining the unpenalized nonparametric MLE")
  unpen <- fun_unpenSurvIC(rep(0, ncol(arglist$z)), arglist)
  final.b0 <- unpen$b
  final.lambda <- unpen$lambda
  log_pen <- log_penlikelihood(final.b0, arglist)
  
  arglist$initial_lambda <- final.lambda
  
  if (covmat == TRUE) {
    message(" Now: calculating the covariance matrix")
    cov <- fun_cov_parallel(b = final.b0, theta = 0, var.h = 5, arglist, cl)
  } else {
    cov <- rep(NA, ncol(arglist$z))
  }
  
  message(" Done.")
  
  if (!is.null(cl)) stopCluster(cl)
  
  if (normalize.X == TRUE) {
    atrue_sd <- (arglist$true_sd)
    final.b <- final.b0/atrue_sd
    final.cov <- cov / (atrue_sd %*%t(atrue_sd))
  } else {
    final.b <- final.b0
    final.cov <- cov
  }
  
  
  results <- list()
  results$xnames <- xnames
  results$n <- nrow(X)
  results$b <- final.b
  results$se <- sqrt(diag(final.cov))
  results$cov <- final.cov
  results$lambda <- final.lambda
  results$lambda.set <- arglist$set
  results$convergence <- unpen$convergence
  results$iteration <- unpen$iteration
  results$ind.trunc <- ind.trunc
  results$smallest.trunc <- ifelse(ind.trunc, min(trunc), 0)
  results$normalize.X <- normalize.X
  results$log_likelihood <- log_pen
  
  class(results) <- "unpencoxIC"
  
  return(results)
}


WVIC_test = function(X, Z=NULL, bound_times, trunct=NULL, Gsim, Fsim=NULL, r, lim=5*10^5, acc=10^(-5), covmat=FALSE){
  #-------------------------------------------------------------------------------#
  # X:           genetic variable to be tested
  # Z:           adjustment covariates (NULL if no adjustment covariates)
  # bound_times: 2-d matrix of (L,R) bound times (col2=Inf if right censored)
  # trunct:      left truncation times
  # Gsim:        Genetic similarity
  # Fsim:        Background similarity (NULL if no heterogeneity)
  # r:           semi-param transformation model parameter (r>=0)
  # lim, acc:    parameters used in Davies method
  # covmat:      TRUE or FALSE (calculate covariance matrix or not)
  #-------------------------------------------------------------------------------#
  
  if(r>0 & !is.null(trunct)){stop('truncation times not supported for semi-parametric transformation model with r>0')}
  
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  
  L_end = bound_times[,1]
  R_end = bound_times[,2]
  
  # left-, right-censor rates
  lc_rate = sum(L_end==0)/nrow(bound_times)
  rc_rate = sum(R_end==Inf)/nrow(bound_times)
  
  # Consider adjustment covariates
  if(!is.null(Z)){
    Z_name = colnames(as.data.frame(Z))
    Z = as.matrix(Z)
    if(r == 0){
      if(covmat){
        fit0 <- suppressMessages(unpencoxIC.default(L_end, R_end, Z, normalize.X=FALSE, covmat=TRUE))
        covar_matrix = fit0$cov
      } else {
        fit0 <- suppressMessages(unpencoxIC.default(L_end, R_end, Z, normalize.X=FALSE, covmat=FALSE))
        covar_matrix = NULL
      }
      jumppts = c(0,baseline(fit0)$upper.set)
      cumhaz = c(0,baseline(fit0)$clambda)
    }
    if(r > 0){
      if(covmat){
        fit0 <- suppressMessages(semipar_trans_fit_rp_adj(Z, bound_time=cbind(L_end,R_end), perturb=5, n.l=3, r_alpha=r, maxiter=100, cov_method='ZL1'))
        covar_matrix = fit0$covest
      } else {
        fit0 <- suppressMessages(semipar_trans_fit_rp_adj(Z, bound_time=cbind(L_end,R_end), perturb=5, n.l=3, r_alpha=r, maxiter=100, cov_method=0))
        covar_matrix = NULL
      }
      jumppts = c(0,fit0$order_bt)
      cumhaz = c(0,cumsum(fit0$lambda_est))
    }
    regcoef = fit0$b
  }
  # Do not consider adjustment covariates
  if(is.null(Z)){
    if(r == 0){
      fit0 <- suppressMessages(semipar_trans_fit_r0_noadj(bound_time=cbind(L_end,R_end), r_alpha=r, maxiter=100))
      covar_matrix = NULL
      jumppts = c(0,fit0$order_bt)
      cumhaz = c(0,cumsum(fit0$lambda_est))
    }
    if(r > 0){
      fit0 <- suppressMessages(semipar_trans_fit_rp_noadj(bound_time=cbind(L_end,R_end), n.l=3, r_alpha=r, maxiter=100))
      covar_matrix = NULL
      jumppts = c(0,fit0$order_bt)
      cumhaz = c(0,cumsum(fit0$lambda_est))
    }
    regcoef = NULL
  }
  
  #message(" Now: calculating the p-value")
  
  if(is.null(Fsim)){
    R = Gsim
  } else {
    R = Gsim*(1+Fsim)
  }
  
  if(is.null(Z)){
    I = diag(n)
    J = matrix(1/n,n,n)
    R.res = (I-J) %*% R %*% (I-J)
  } else {
    I = diag(n)
    Z.bar = cbind(rep(1,n),Z)
    hatmatrix = Z.bar %*% solve(t(Z.bar)%*%Z.bar) %*% t(Z.bar)
    R.res = (I-hatmatrix) %*% R %*% (I-hatmatrix)
  }
  
  if(is.null(Z)){
    if(is.null(trunct)){
      s = apply(cbind(L_end,R_end), 1, FUN=s_cal_noadj, jumppts=jumppts, cumhaz=cumhaz, r=r)
    } else {
      s = apply(cbind(trunct,L_end,R_end), 1, FUN=s_cal_noadj_trunc, jumppts=jumppts, cumhaz=cumhaz, r=r)
    }
  } else {
    if(is.null(trunct)){
      s = apply(cbind(L_end,R_end,Z), 1, FUN=s_cal_adj, jumppts=jumppts, cumhaz=cumhaz, regcoef=regcoef, r=r)
    } else {
      s = apply(cbind(trunct,L_end,R_end,Z), 1, FUN=s_cal_adj_trunc, jumppts=jumppts, cumhaz=cumhaz, regcoef=regcoef, r=r)
    }
  }
  eta = eigen(R.res, symmetric=TRUE, only.values=TRUE)$values
  
  lambda1 = as.numeric(t(s) %*% s)
  nV = c(t(s) %*% R.res %*% s / n)
  
  if(is.null(Z)){
    Davies = CompQuadForm::davies(nV, lambda1*eta/(n^2), lim=lim, acc=acc)
  } else {
    Davies = CompQuadForm::davies(nV, lambda1*eta/(n*(n-ncol(Z)-1)), lim=lim, acc=acc)
  }
  
  output <- list()
  output$coef_est = regcoef         #regression coef estimates
  output$CovEst = covar_matrix      #covariance matrix estimate for regression coef
  output$pval = Davies$Qq           #p-value
  output$ifault = Davies$ifault     #p-value calculation error
  output$lcr = lc_rate              #left-censor rate
  output$rcr = rc_rate              #right-censor rate
  return(output)
} #end of function#




cpICSCAT_pwr <- function(alpha1, alpha2, beta1, beta2, numCausal, effectSizes, gMat=NULL, WVIC=FALSE, maxMAF=0.05,
                         n, q=NULL, rho=NULL, init_1 = NULL, quants = NULL, numcores = numcores, B1, B){
  library(tidyverse)
  library(cmprsk)
  library(ggplot2)
  library(nleqslv)
  library(CompQuadForm)
  library(ICSKAT)
  library(bindata)
  library(data.table)
  library(purrr)
  library(rlist)
  library(foreach)
  require(doSNOW)
  library(parallel)
  library(iterators)
  library(snow)
  
  # Start parallel computing
  cl <- makeCluster(numcores, type = "SOCK", outfile = "clusterout.txt")
  registerDoSNOW(cl)
  
  oper <- foreach(sim_it = B1:B, .inorder = F, .errorhandling = "pass", .export=ls(.GlobalEnv)) %dopar% {
    
    # Set seed for covariates, genotypes
    set.seed(sim_it)
    
    idx <- 1:numCausal
    
    xMat <- cbind(rnorm(2*n), rbinom(n=2*n, size=1, prob=0.5))
    if(is.null(gMat)){
      gMat <- sim_gmat(2*n, q, rho, maxMAF=maxMAF)
    }else{
      gMat <- rbind(gMat, gMat)
    }
    
    q <- ncol(gMat)
    
    gSummed <- matrix(data=apply(gMat, 1, sum), ncol=1)
    gMatCausal <- gMat[, idx]
    
    # Generate the data and then select the initial value automatically
    tempDat_Large <- genData_pwr(seed=sim_it, n=2*n, alpha1, alpha2, beta1, beta2, 
                                 gMatCausal=gMatCausal, effectSizes=effectSizes)
    
    nonNAN <- which(!is.nan(tempDat_Large$tempTime))[1:n]
    tempDat <- list(tempTime = tempDat_Large$tempTime[nonNAN],
                    tempType = tempDat_Large$tempType[nonNAN])
    xMat <- xMat[nonNAN, ]
    gMat <- gMat[nonNAN, ]
    gSummed <- matrix(data=apply(gMat, 1, sum), ncol=1)
    gMatCausal <- gMat[, idx]
    
    
    
    if(is.null(quants)){
      quants <- stats::quantile(log(tempDat$tempTime), probs=seq(from=0, to=1, length.out=3))
    }else{
      quants <- quants
    }
    
    # ---- Get the estimate of initial value ---
    ej <- (max(quants) - quants[2]) / (max(quants) - min(quants))
    logH1 <- rep(log(-beta1 / alpha1), length(tempDat$tempTime)) + log(1 - exp(alpha1 * tempDat$tempTime))
    logH2 <- rep(log(-beta2 / alpha2), length(tempDat$tempTime)) + log(1 - exp(alpha2 * tempDat$tempTime))
    regDmat <- cbind(1, log(tempDat$tempTime), 
                     pmax(0, (log(tempDat$tempTime) - quants[2])**3) - ej * pmax(0, (log(tempDat$tempTime) - quants[1])**3) -
                       (1 - ej) * pmax(0, (log(tempDat$tempTime) - quants[3])**3))
    if(is.null(init_1)){
      init_1 <- c(as.numeric(summary(lm(logH1 ~ cbind(xMat, regDmat) - 1))$coef[, 1]),
                  as.numeric(summary(lm(logH2 ~ cbind(xMat, regDmat, gSummed) - 1))$coef[, 1]))
    } else{
      init_1 <- init_1
    }
    
    init_2 <- c(0,0,0,0,0.1)
    # -------------------------
    
    
    madeVisit <- matrix(data=rbinom(n=n*7, size=1, prob=0.9), nrow=n, ncol=7)
    visitTime <- sweep(matrix(data=runif(n=n*7, min=-1, max=1), nrow=n, ncol=7),
                       MARGIN=2, STATS=seq(from=4, to=28, by=4), FUN="+")
    # get all visits for each subject
    allVisits <- madeVisit * visitTime
    # make the interval for each subject - USING creatIntNew!
    allInts <- t(mapply(FUN=createIntNew, obsTimes = data.frame(t(allVisits)), eventTime=tempDat$tempTime))
    leftTimes <- allInts[, 1]
    rightTimes <- allInts[, 2]
    # tpos <- allInts[, 3]
    # obsind <- allInts[, 4]
    
    # new indicator of whether event was observed
    deltaVecSimple <- ifelse(rightTimes > tempDat$tempTime, 1, 0)
    deltaVec <- deltaVecSimple * tempDat$tempType
    
    # make spline terms
    # for interval censoring, our earliest visit time is 4 and our latest is 28, windows of size 1
    
    dmats <- makeICdmat(xMat=xMat, lt = leftTimes, rt = rightTimes, obsInd = deltaVecSimple, 
                        quant_r = quants, nKnots=1) 
    leftDmat <- dmats$left_dmat
    rightDmat <- dmats$right_dmat
    
    # solve score equations under null with no random effect
    solPartial <- nleqslv::nleqslv(x = init_1, 
                                   fn=scoreEqSpline, leftDmat = leftDmat,  rightDmat = rightDmat, leftTimes = leftTimes,
                                   deltaVec = deltaVec, gSummed = gSummed, gMat=NULL, estG = FALSE)
    
    if (solPartial$termcd > 2) {
      next
    }
    
    # Ugamma with null coefficients
    forUgamma <- scoreEqSpline(x=solPartial$x, leftDmat = leftDmat, 
                               rightDmat = rightDmat, leftTimes = leftTimes, deltaVec = deltaVec,
                               gSummed = gSummed, gMat = gMat, estG = FALSE)
    
    # information
    # partial information, under null
    iMatPartial <- calcInfo(leftDmat = leftDmat, rightDmat = rightDmat, leftTimes = leftTimes, 
                            theta1 = solPartial$x[1:5], theta2 = solPartial$x[6:11], deltaVec = deltaVec,
                            gSummed=gSummed, gMat = NULL)
    
    # Igt
    forIgt <- calcInfo(leftDmat = leftDmat, rightDmat = rightDmat, leftTimes = leftTimes, 
                       theta1 = c(solPartial$x[1:5], rep(0, q)), theta2 = solPartial$x[6:11], deltaVec = deltaVec,
                       gSummed = gSummed, gMat = gMat)
    
    # try SKAT
    skatQ <- t(forUgamma[6:(5+q)]) %*% forUgamma[6:(5+q)]
    burdenQ <- (sum(forUgamma[6:(5+q)]))^2
    Itt <- -iMatPartial
    Igg <- -forIgt[6:(5+q), 6:(5+q)]
    Igt <- -forIgt[c(1:5, (6+q):nrow(forIgt)), 6:(5+q)]
    sig_mat <- Igg - t(Igt) %*% solve(Itt) %*% (Igt)
    lambdaQ <- eigen(sig_mat)$values
    p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=1e-7)$Qq
    B_burden <- burdenQ / sum(sig_mat)
    p_burden <- 1 - stats::pchisq(B_burden, df = 1)
    
    if (!is.na(p_SKAT)) {
      if (p_SKAT > 1) {
        paramDF <- data.frame(expand.grid(lim = c(10000, 20000, 50000), acc=c(1e-7, 1e-6, 1e-5, 1e-4)))
        paramCounter <- 1
        while(p_SKAT > 1) {
          tempLim <- paramDF$lim[paramCounter]
          tempAcc <- paramDF$acc[paramCounter]
          p_SKAT <- CompQuadForm::davies(q=skatQ, lambda=lambdaQ, delta=rep(0,length(lambdaQ)), acc=tempAcc, lim=tempLim)$Qq
          paramCounter <- paramCounter + 1
          if (paramCounter > nrow(paramDF)) {break}
        }
        errCode <- 22
        errMsg <- "Had to adjust parameters on CompQuadForm"
      }
    }
    
    
    # ----- ICSKAT
    # Here there is some modification
    # The leftTime and RightTime are different because the createIntNew function
    # use left time as right time for the right-censored subjects.
    # also use 0 as left time for the right-censored subjects.
    obsind <- ifelse(deltaVec == 1, 1, 0)    # Treat cause2 as censoring
    tpos <- ifelse(leftTimes == 0, 0, 1)
    # tpos <- allInts[, 3]
    # obsind <- allInts[, 4]
    leftTimes1 <- allInts[, 1]
    rightTimes1 <- allInts[, 2]
    leftTimes1[which(obsind==0)] <- rightTimes1[which(obsind==0)]
    rightTimes1[which(obsind==0)] <- 999
    
    dmatsICSKAT <- ICSKAT::make_IC_dmat(xMat=xMat, lt=leftTimes1, rt=rightTimes1, obs_ind=obsind, tpos_ind=tpos, nKnots=1)
    nullFit <- ICSKAT::ICSKAT_fit_null(init_beta=init_2, 
                                       left_dmat=dmatsICSKAT$left_dmat, right_dmat=dmatsICSKAT$right_dmat,
                                       obs_ind=obsind, tpos_ind=tpos, lt=leftTimes1, rt=rightTimes1)
    
    if (is.na(as.numeric(nullFit$beta_fit[1]))){
      p_SKAT2 <- NA
      p_burden2 <- NA
      skatQ2 <- NA
      burdenQ2 <- NA
      err <- nullFit$errMsg
    } else{
      # perform the ICSKAT and Burden tests
      icskatOut <- ICSKAT::ICskat(left_dmat=dmatsICSKAT$left_dmat, right_dmat=dmatsICSKAT$right_dmat, 
                                  lt=leftTimes1, rt=rightTimes1,
                                  obs_ind=obsind, tpos_ind=tpos, gMat = gMat, 
                                  null_beta = as.numeric(nullFit$beta_fit), Itt = nullFit$Itt)
      p_SKAT2 <- as.numeric(icskatOut$p_SKAT)
      p_burden2 <- as.numeric(icskatOut$p_burden)
      skatQ2 <- as.numeric(icskatOut$skatQ)
      burdenQ2 <- as.numeric(icskatOut$burdenQ)
      err <- icskatOut$errMsg
    }
    
    
    ### ---- WVIC
    if(WVIC == TRUE){
      Gsim <- as.matrix(gMat) %*% t(as.matrix(gMat))
      leftTimes2 <- allInts[, 1]
      rightTimes2 <- allInts[, 2]
      ltWVIC <- ifelse(obsind == 1, leftTimes2, rightTimes2)
      rtWVIC <- ifelse(obsind == 1, rightTimes2, Inf)
      bound_times <- cbind(ltWVIC, rtWVIC)
      test_results <- WVIC_test(X=as.matrix(gMat),
                                Z=NULL,
                                bound_times=bound_times,
                                trunct=NULL,
                                Gsim=Gsim, Fsim=NULL,
                                r=0, lim=5*10^5, acc=10^(-5),
                                covmat=FALSE)
      pWVIC <- as.numeric(test_results$pval)
      rm(Gsim)
      
      if(!is.numeric(pWVIC)){
        pWVIC <- -99
      }
    }else{
      pWVIC <- -99
    }
    
    
    # # Record
    # cp_skat_pval[i] <- as.numeric(p_SKAT)
    # cp_burden_pval[i] <- as.numeric(p_burden)
    # cp_skatQ[i] <- as.numeric(skatQ)
    # cp_burdenQ[i] <- as.numeric(burdenQ)
    # ic_skat_pval[i] <- p_SKAT2
    # ic_burden_pval[i] <- p_burden2
    # ic_skatQ[i] <- skatQ2
    # ic_burdenQ[i] <- burdenQ2
    # errMsg[i] <- err
    
    res_parallel <- as.numeric(c(p_SKAT, p_burden, skatQ, burdenQ, p_SKAT2, p_burden2, skatQ2, burdenQ2, pWVIC, sim_it))
    return(res_parallel)
    
    # if(i %in% seq(0, runs, 20)){cat(i/runs*100, "% has been processed", "\n")}
    
  }
  
  # Stop parallel computing
  stopCluster(cl)
  filtered <- Filter(function(x) length(x) == 10, oper)
  skatDF <- as.data.frame(do.call(rbind, filtered))
  # oper <- oper[which(sapply(1:(B-B1+1), function(i) is.vector(oper[[i]])))]
  # oper <- t(as.data.frame(oper))
  colnames(skatDF) <- c("crICSKAT_pval", "crBurden_pval", "crICSKAT_Q", "crBurden_Q",
                        "ICSKAT_pval", "ICBurden_pval", "ICSKAT_Q", "ICBurden_Q", "pWVIC", "sim_it")
  rownames(skatDF) <- NULL
  skatDF <- as.data.frame(skatDF)
  
  return(skatDF)
  
}






































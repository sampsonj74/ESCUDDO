###################################################################################################
###################################################################################################

#' Estimate the Probability that a participant had both (i) no baseline (0/6 month) infection (ii) an incident/persistent HPV-16 (or HPV-18) Infection
#'
#' Estimate the Probability that a participant had both (i) no baseline (0/6 month) infection (ii) an incident/persistent HPV-16 (or HPV-18) Infection
#'
#' @export
#' @param dMat16x A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are HPV 16 (or HPV 18) infection status (w/ -9 for missing)
#' @param sMatx:  A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are of sexual activity status (w/ -9 for missing)
#' @param ATPx:   A vector of length number.of.subjects. Entries show whether the subject was part of the ATP group.
#' @return    infecOUT. A matrix (nrow=length number.of.subjects, ncol = 2). The first column indicates if an event was observed. The second column shows the probability that
#' an event occurred. An event is two baseline negative results (i.e. 0 and 6 months) and an incident prevalent infection (two consecutive positive visits -- no intervening
#' negative or missing visits). non-ATP individuals are set to missing. The second column is based on the methods described in Sampson, Gail, CCT.
#'
#' @examples
#' trialRes <- simTrialArm(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",yr=5,nsub=5000,pdis=0.01,pdo6=0.05,pdo=0.02,pmv=0.17)
#' dMat16x   <- trialRes$dMat
#' sMatx     <- trialRes$sMat
#' ATPx      <- trialRes$ATPx
#' probInf2   <- probInf(dMat16x,sMatx,ATPx)

probInf     <- function(dMat16x,sMatx,ATPx){

  ### Any women with an infection must be sexually active (and sexually active at all future visits)
  sMatx     <- ifelse(dMat16x==1,1,sMatx)
  for (i in 2:ncol(sMatx)) sMatx[,i] <- ifelse(sMatx[,i-1]==1,1,sMatx[,i])

  ### The output -- observed events and probability of an event
  infecOut <- matrix(nrow=nrow(dMat16x),ncol=2)
  colnames(infecOut) <- c("ObservedEvent","ProbEvent")

  ###ATP population
  ATP      <- c(1:nrow(dMat16x))[ATPx==1]
  dMat16y  <- dMat16x[ATP,-c(1:2)]
  sMaty    <- sMatx[ATP,-c(1:2)]

  ###Defaults so missing baseline sexual activity status is set to 0
  s2       <- ifelse(is.na(sMatx[ATP,2]),0,sMatx[ATP,2])
  p        <- ncol(dMat16y)
  nsub     <- nrow(dMat16y)
  dMat16y  <- ifelse(dMat16y==-9,NA,dMat16y)

  ###We want to get baseline infection status
  dMat16x2 <- ifelse(is.na(dMat16x[,c(1:2)]),0,dMat16x[,c(1:2)])
  blinfy   <- ifelse(rowSums(dMat16x2[ATP,c(1,2)])>0,1,0)

  ###identify the starting position of each new gap  and the position of when that gaps ends
  endGap     <- newGap <- matrix(nrow=nsub,ncol=p)
  newGap[,1] <- ifelse(is.na(dMat16y[,1]),1,0)
  for (i in 2:p)     newGap[,i]   <- ifelse(is.na(dMat16y[,i]) & !is.na(dMat16y[,i-1]),1,0)
  for (i in (p-1):1) endGap[,i]   <- ifelse(is.na(dMat16y[,i+1]),endGap[,i+1],i+1)


  ###dMat2 is designed to fill in all values that do not matter (i.e. that couldn't define an infection) w/ zeros
  dMat2     <- dMat16y
  dMat2[,1] <- ifelse(is.na(dMat16y[,1]) & dMat16y[,2]  ==0 & !is.na(dMat16y[,2]),  0,dMat16y[,1])
  dMat2[,p] <- ifelse(is.na(dMat16y[,p]) & dMat16y[,p-1]==0 & !is.na(dMat16y[,p-1]),0,dMat16y[,p])
  for (i in 2:(p-1)) dMat2[,i] <- ifelse(is.na(dMat16y[,i]) & dMat16y[,i-1]==0 & !is.na(dMat16y[,i-1]) &
                                                              dMat16y[,i+1]==0 & !is.na(dMat16y[,i+1]),0,dMat16y[,i])
  dM1       <- dMat2[,2:p]
  dM2       <- dMat2[,1:(p-1)]

  ###infec identifes the infection status of subjects that are clearly infected or uninfected
  noInfecObs   <- ifelse(rowSums( (1-dM1)*(1-dM2)) == ncol(dM1) & !is.na(rowSums( (1-dM1)*(1-dM2))),1,0)
  infecObs     <- rowSums(dM1*dM2,na.rm=T)>0 & blinfy == 0
  infec        <- ifelse(blinfy==1,0,ifelse(infecObs==1,1,ifelse(noInfecObs==1,0,NA)))
  missInf      <- c(1:nsub)[is.na(infec)]
  infecOut[ATP,1] <- infecObs

  ####for each column, among all individuals who are 1, calculate the probability that they were infected at the previous visit
  m1 <- matrix(nrow=2,ncol=p)
  for (s in 1:2) for (i in 2:p)     m1[s,i] <- mean(dMat16y[dMat16y[,i]==1 & !is.na(dMat16y[,i]) & sMaty[,i] == (s-1) & blinfy==0 ,i-1] ,na.rm=T)
  m1  <- ifelse(is.na(m1),1,m1)


  ####for each column, among all individuals who are 1, calculate the probability that they were infected at the next  visit
  p1  <- matrix(nrow=2,ncol=p)
  for (s in 1:2) for (i in 1:(p-1)) p1[s,i] <- mean(dMat16y[dMat16y[,i]==1 & !is.na(dMat16y[,i]) & sMaty[,i-1] == (s-1) & blinfy==0,i+1],na.rm=T)
  p1  <- ifelse(is.na(p1),1,p1)


  ####for each column, among all individuals who are 0, calculate the probability that they were infected at some later time
  l1 <- matrix(nrow=2,ncol=p)
  for (s in 1:2){
  for (i in 1:(p-2)) {ok   <- c(1:nsub)[rowSums(is.na(dMat16y[,i:p]))==0 & dMat16y[,i]==0 & sMaty[,i]==(s-1) & blinfy==0]
                      dMl1       <- dMat16y[,(i+1):p]
                      dMl2       <- dMat16y[,i:(p-1)]
                      if (i < (p-1))  l1[s,i]      <- mean(rowSums(dMl1[ok,]*dMl2[ok,])>0) }
  }
  l1 <- ifelse(is.na(l1),0,l1)



  ####for each subject w/ missing status, we walk through all gaps, identify the appropriate group of comparison individuals,
  ####calculate the probability that they are infected during the gap
  ###Note: puninf2 is only used to estimate the variance of imputation
  noMatch <- type1 <- type2 <- type3 <- type4 <- c()
  for (i in missInf) {         puninf <- 1
                               for (k in c(1:p)[newGap[i,]==1]){
                                     if (k==1 & !is.na(endGap[i,k]) )              {sub    <- c(1:nsub)[rowSums(is.na(dMat16y[,1:endGap[i,1]]))   == 0 & blinfy==0 &
                                                                                                        dMat16y[,endGap[i,1]]                     == dMat16y[i,endGap[i,1]] &
                                                                                                        sMaty[,endGap[i,1]]                       == sMaty[i,endGap[i,1]] &
                                                                                                        s2                                        == s2[i] ]
                                                                                    dMat.T <- dMat16y[sub,1:endGap[i,1]]
                                                                                    type1 <- c(type1,i)}
                                     if (k >1  & !is.na(endGap[i,k]) )             {sub    <- c(1:nsub)[rowSums(is.na(dMat16y[,(k-1):endGap[i,k]]))  == 0   & blinfy==0 &
                                                                                                     dMat16y[,(k-1)]                                == dMat16y[i,k-1] &
                                                                                                     dMat16y[,endGap[i,k]]                          == dMat16y[i,endGap[i,k]]  &
                                                                                                     sMaty[,(k-1)]                                  == sMaty[i,k-1] &
                                                                                                     sMaty[,endGap[i,k]]                            == sMaty[i,endGap[i,k]]]
                                                                                    dMat.T <- dMat16y[sub,(k-1):endGap[i,k]]
                                                                                    type2 <- c(type2,i)}
                                     if (k >1 &   is.na(endGap[i,k]) )             { sub    <- c(1:nsub)[rowSums(is.na(dMat16y[,(k-1):p]))==0 & blinfy==0 &
                                                                                                         dMat16y[,(k-1)] == dMat16y[i,k-1] &
                                                                                                         sMaty[,(k-1)]   == sMaty[i,k-1]   ]
                                                                                    dMat.T <- dMat16y[sub,(k-1):p]
                                                                                    type3 <- c(type3,i)}
                                     if (k==1 &   is.na(endGap[i,k]) )             {sub    <- c(1:nsub)[rowSums(is.na(dMat16y[,1:p]))==0 & blinfy==0 &
                                                                                                        s2                           == s2[i]]
                                                                                    dMat.T <- dMat16y[sub,1:p]
                                                                                    type4 <- c(type4,i)}
                                     if (length(sub)==1) dMat.T <- matrix(dMat.T,nrow=1)
                                     if (length(sub)>0)  dM1.T  <- matrix(dMat.T[,2:ncol(dMat.T)],nrow=length(sub))
                                     if (length(sub)>0)  dM2.T  <- matrix(dMat.T[,1:(ncol(dMat.T)-1)],nrow=length(sub))
                                     if (length(sub)>0)  puninf <- puninf*(1-mean(rowSums(dM1.T*dM2.T)>0))
                                     if (length(sub)==0){
                                      noMatch <- c(noMatch,i)
                                      sseg <- sspg <- NA
                                      if (!is.na(endGap[i,k])) sseg <- sMaty[i,endGap[i,k]] + 1
                                      if (k > 1)               sspg <- sMaty[i,k-1] + 1
                                      if (k == 1 & !is.na(endGap[i,k])) if (dMat16y[i,endGap[i,k]]==0)   puninf <- puninf
                                      if (k == 1 & !is.na(endGap[i,k])) if (dMat16y[i,endGap[i,k]]==1)   puninf <- 1-m1[sseg,endGap[i,k]]
                                      if (k >  1 & !is.na(endGap[i,k])) if (dMat16y[i,endGap[i,k]]==0 & dMat16y[i,k-1]==0)   puninf <- puninf
                                      if (k >  1 & !is.na(endGap[i,k])) if (dMat16y[i,endGap[i,k]]==1 & dMat16y[i,k-1]==1)   puninf <- 0
                                      if (k >  1 & !is.na(endGap[i,k])) if (dMat16y[i,endGap[i,k]]==1 & dMat16y[i,k-1]==0)   puninf <- puninf*(1-m1[sseg,endGap[i,k]])
                                      if (k >  1 & !is.na(endGap[i,k])) if (dMat16y[i,endGap[i,k]]==0 & dMat16y[i,k-1]==1)   puninf <- puninf*(1-p1[sspg,k-1])
                                      if (k == 1 &  is.na(endGap[i,k])) puninf <- NA
                                      if (k >  1 &  is.na(endGap[i,k])) if (dMat16y[i,k-1]==1 )                           puninf <- puninf*(1-p1[sspg,k-1])
                                      if (k >  1 &  is.na(endGap[i,k])) if (dMat16y[i,k-1]==0 )                           puninf <- puninf*(1-l1[sspg,k-1])
                                     }
                               }
                               infec[i]  <- 1-puninf
  }


infecOut[ATP,2] <- ifelse(is.na(infec),mean(infec,na.rm=TRUE),infec)
return(infecOut)
}

###################################################################################################
###################################################################################################

#' Obtain the observed and expected number of events in one arm of the trial
#'
#' Obtain the observed and expected number of events in one arm of the trial
#'
#' @export
#'
#' @param dMat16x A matrix (nrow=number.of.subjects, ncol=number.of.visits.including.baseline). Entries are 0 or 1, indicating HPV 16 infection status (w/ -9 for missing)
#' @param dMat18x A matrix (nrow=number.of.subjects, ncol=number.of.visits.including.baseline). Entries are 0 or 1, indicating HPV 18 infection status (w/ -9 for missing)
#' @param sMatx:  A matrix (nrow=number.of.subjects, ncol=number.of.visits.including.baseline). Entries are 0 or 1, indicating if a woman has had her sexual debut (w/ -9 for missing)
#' @param ATPx:   A vector of length number.of.subjects. Entries are 0 or 1, indicating whether the subject was part of the ATP group.
#' @return    N:      A single value. The number of women in the ATP cohort.
#' @return    nobs:   A single value. The number of ATP women OBSERVED to have an event  (i.e. no 0/6 month infection and a incident persistent infection)
#' @return    nexp:   A single value. The number of ATP women EXPECTED to have an event  (i.e. no 0/6 month infection and a incident persistent infection)
#' @return    infFact:   A single value. The inflation factor = nexp/nobs
#' @return    Nequiv:    A single value. The equivalent sample size = N/infFact
#' @examples
#' trialRes <- simTrialArm(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",yr=5,nsub=5000,pdis=0.01,pdo6=0.05,pdo=0.02,pmv=0.17)
#' dMat16x   <- trialRes$dMat16
#' dMat18x   <- trialRes$dMat18
#' sMatx     <- trialRes$sMat
#' ATPx      <- trialRes$ATPx
#' ecount    <- eventCount(dMat16x,dMat18x,sMatx,ATPx)

eventCount <- function(dMat16x,dMat18x,sMatx,ATPx){

  ###The observed and expected infection status for each woman for each individual type
  p16                <- probInf(dMat16x,sMatx,ATPx)
  p18                <- probInf(dMat18x,sMatx,ATPx)

  ###The observed and expected infection status for each woman for both types combined
  infecOut           <- matrix(nrow=nrow(dMat16x),ncol=2)
  colnames(infecOut) <- c("ObservedEvent","ProbEvent")


  infecOut[,1]       <- ifelse( (p16[,1] == 1 & !is.na(p16[,1])) |
                                (p18[,1] == 1 & !is.na(p18[,1])), 1,
                                 ifelse(is.na(p16[,1]) & is.na(p18[,1]),NA,0))

  infecOut[,2]       <- ifelse(is.na(p16[,2]) & is.na(p18[,2]),NA,
                              ifelse(is.na(p16[,2]) & !is.na(p18[,2]),p18[,2],
                                   ifelse(!is.na(p16[,2]) & is.na(p18[,2]),p16[,2],1-(1-p16[,2])*(1-p18[,2]))))

  N       <- sum(ATPx)
  nobs    <- sum(infecOut[,1],na.rm=TRUE)
  nexp    <- sum(infecOut[,2],na.rm=TRUE)
  infFact <- nexp/nobs
  Nequiv  <- N/infFact

  list("N"=N,"nobs"=round(nobs),"nexp"=round(nexp),"infFact"=infFact,"Nequiv"=round(Nequiv))
}


###################################################################################################
###################################################################################################

#' Generate a set of trial participants
#'
#' Generate a set of trial participants
#' @export
#' @param ofile  A string. The file name for the file containing the key parameters (generated by calcParam()).
#' @param yr     A single value. The number of years of the trial (either  4 or 5)
#' @param nsub   A single value. The number of subjects the arm.
#' @param pdis   A single value. The risk of being baseline uninfected and acquiring a new 16 or 18 persistent infection.
#' @param pdisL2 A single value. The risk of being baseline uninfected and having a HPV 16 or 18 persistent infection at the last two visits
#' (NOTE: only pdis or pdisL2 should be set to a non-missing value)
#' @param pdo6   A single value. The probability of dropping out before the 6 month visit
#' @param pdo    A single value. The probability of dropping out before any other follow-up visit
#' @param pmv    A single value. The probability of missing a visit (or specimen failure) at all follow-up vists after 6 months
#' @return dMat:      A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are infection HPV16/18 status (w/ -9 for missing)
#' @return dMat16:    A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are HPV 16 infection status (w/ -9 for missing)
#' @return dMat18:    A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are HPV 18 infection status (w/ -9 for missing)
#' @return sMat:      A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are of sexual activity status (no missing).
#' @return infecTrue: A vecor of length number.of.subjects. Entries show whether the subject was both baseline negative and had a persistent infection.
#' @return ATPx:      A vecor of length number.of.subjects. Entries show whether the subject received both vaccinations.
#' @examples
#' trialRes <- simTrialArm(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",yr=5,nsub=5000,pdis=0.01,pdo6=0.05,pdo=0.02,pmv=0.17)

simTrialArm <- function(ofile=NA,yr=5,nsub=5000,pdis=NA,pdisL2=NA,pdo6=0.05,pdo=0.02,pmv=0.17,seed=NA){

        if (!is.na(seed)) set.seed(seed)
        if (!file.exists(ofile)) stop('Error: ofile does not exist; have you run calcParam?')
        load(ofile)

        ###Recalibrate infection probabilities
        ###Initially they were calculated based on the best assumptions using calcParam
        ###Now, we ensure that pdis or pdisL2 is true

        if (yr == 4) {prob    <- infProb
                      poss    <- infPoss
                      nnew    <- nnew.4yr
                      infPers <- infPers.4yr
                      svisit  <- sexVisit
                      nc      <- 9}

        if (yr == 5) {prob    <- infProb11
                      poss    <- infPoss11
                      nnew    <- nnew.5yr
                      infPers <- infPers.5yr
                      svisit  <- sexVisit11
                      nc      <- 11}

        if (!is.na(pdis) & !is.na(pdisL2)) stop('Error: Either pdis and pdisL2 must be set to NA')

        ### The definitiion of infPers -- these are states (i.e rows of poss) that qualify
        ### as diseased (differs for survey or trial)
        if (!is.na(pdisL2)) infPers <- ifelse(poss[,nc-1] == 1 & poss[,nc] == 1 & poss[,1] == 0,1,0)
        pNewInf <- sum(prob[infPers==1])
        if (pdis   > pNewInf & !is.na(pdis))   stop(paste0("Disease Rate is too high: please lower pdis below ",round(pNewInf,2)))
        if (pdisL2 > pNewInf & !is.na(pdisL2)) stop(paste0("Disease Rate is too high: please lower pdisL2 below ",round(pNewInf,2)))

        ###This is the caliration step (i.e. we change prob to prob3 -- so that we achieve pdis or pdisL2)
        pTerm     <- 1
        pdisx     <- ifelse(!is.na(pdis),pdis,pdisL2)
        while (sum(pTerm^nnew[infPers==1]*prob[infPers==1])>pdisx) pTerm=pTerm-0.001
        prob3 <- prob2 <- prob*pTerm^nnew
        dif   <-  1-sum(prob2)
        prob3[rowSums(poss)==0] <- prob2[rowSums(poss)==0]   + dif

        ###Check to make sure infection probability matches pdis
        ###c(sum(prob3[infPers==1]),pdis)

        ###Identify the true disease mat of these subjects
        dType <- sample(c(1:nrow(poss)),nsub,replace=TRUE,prob=prob3)
        dMat <- poss[dType,]
        infecTrue <- ifelse(rowSums(dMat[,-c(1:3)]*dMat[,-c(1:2,ncol(dMat))])>0 & rowSums(dMat[,c(1:2)])==0,1,0)

        ###Identify the sexual activity status of each woman
        sMat <- matrix(0,nrow=nsub,ncol=ncol(poss))
        svisit <- ifelse(!is.na(svisit),svisit,1/ncol(poss))
        for (i in 1:nsub) { firstVal <- sample(c(1:(ncol(poss)+1)),1,prob=svisit[dType[i],])
                            if (firstVal <= ncol(poss)) sMat[i,c(firstVal:ncol(poss))] <- 1 }

        ###Matrix showing missed visits
        npv     <- ncol(poss)
        mMat    <- matrix(nrow=nsub,ncol=npv)
        mMat[,1] <- 0
        ATPx    <- rep(1,nsub)
        for (i in 2:npv) {  pdox <- ifelse(i==2,pdo6,pdo)
                            mMat[,i] <- ifelse(mMat[,i-1]==1,1,rbinom(nsub,1,pdox))
                            if (i==2) ATPx <- 1-mMat[,i]}
        for (i in 3:npv)    mMat[,i] <- ifelse(mMat[,i]==1,  1,rbinom(nsub,1,pmv))
        dMat[mMat==1] <- -9
        sMat[mMat==1] <- -9
        for (i in 2:npv) sMat[,i] <- ifelse(sMat[,i-1]==1,1,sMat[,i])
        dMat16 <- dMat18 <- dMat


        ###Want separate matrices for HPV16 and HPV18
        infWomen <- c(1:nrow(dMat))[rowSums(dMat==1) > 0]
        inf16    <- sample(infWomen,round(0.67*length(infWomen)))
        inf18    <- infWomen[!is.element(infWomen,inf16)]
        wMat     <- matrix(c(1:nsub),nrow=nsub,ncol=ncol(dMat),byrow=FALSE)
        dMat16[!is.element(wMat,inf16) & dMat==1] <- 0
        dMat18[!is.element(wMat,inf18) & dMat==1] <- 0
        return(list("dMat"=dMat,"sMat"=sMat,"dMat16"=dMat16,"dMat18"=dMat18,"ATPx"=ATPx,"infecTrue"=infecTrue))
}


###################################################################################################
###################################################################################################

#' Obtain the confidence interval for the difference in two rates (method 3: C.P. Farrington and G. Manning
#' Stat. Med., 9 (1990), pp. 1447-1454).
#'
#' Obtain the confidence interval for the difference in two rates
#' @export
#' @param n1     A single value. Number of events in the first group
#' @param N1     A single value. Number of individuals in the first group
#' @param n2     A single value. Number of events in the second group
#' @param N2     A single value. Number of individuals in the second group
#' @param  c     A single value. A confidence level (default c = 0.95)
#' @return L     A single value. Lower Bound
#' @return U     A single value. Upper Bound
#' @return D     A single value. The difference in rates (n1/N1-n2/N2)
#' @examples
#' confInt(5,100,4,100)
#'
confInt <- function(n1,N1,n2,N2,c=0.95){
  dObs <- L <- U <- NA
  dObs <- calc.dif(n1,N1,n2,N2)
  try(if (!is.na(dObs)) L <- calc.bound(n1,N1,n2,N2,bound="lb",c=c))
  try(if (!is.na(dObs)) U <- calc.bound(n1,N1,n2,N2,bound="ub",c=c))
  list("D"=dObs,"L"=L,"U"=U)
}



###Calculates the difference in event probabilities

calc.dif <- function(n1,N1,n2,N2){
  p1 <- n1/N1
  p2 <- n2/N2
  p1-p2
}

###Calculates each bound of the 95% confidence interval for difference

calc.bound <- function(n1,N1,n2,N2,bound="ub",c=0.95){
  alpha <- (1-c)/2
  N  <- N1+N2
  p1 <- n1/N1
  p2 <- n2/N2
  m1 <- n1+n2
  if (bound=="ub") d0 <- seq(p1-p2,  0.1,length.out=10000)
  if (bound=="lb") d0 <- seq(p1-p2, -0.1,length.out=10000)
  x21<- n2
  L3 <- N
  L2 <- (N+n2)*d0-N-m1
  L1 <- (n2*d0-N-2*x21)*d0+m1
  L0 <- x21*d0*(1-d0)
  C  <- L2^3/(27*L3^3)-L1*L2/(6*L3^2)+L0/(2*L3)
  B  <- sign(C)*sqrt(L2^2/(9*L3^2)-L1/(3*L3))
  A  <- (1/3)*(pi+acos(C/B^3))
  p2.t <- 2*B*cos(A)-L2/(3*L3)
  p1.t <- p2.t+d0
  z = (p1-p2-d0)/sqrt(p1.t*(1-p1.t)/N1+p2.t*(1-p2.t)/N2)
  if (bound=="ub") b <- min(d0[1-pnorm(abs(z))< alpha])
  if (bound=="lb") b <- max(d0[1-pnorm(abs(z))< alpha])
  b
}


###################################################################################################
###################################################################################################

#' Generate a set of survey participants
#'
#' Generate a set of survey participants
#' @export
#' @param ofile  A string. The file name for the file containing the key parameters.
#' @param yr     A single value. The number of years of the comparison trial group (either  4 or 5)
#' @param nsub   A single value. The number of in the survey.
#' @param pdisL2 A single value. The risk of having a recently acquired HPV 16 or 18 persistent infection at the last two visits
#' (usually this the value pdisL2.trial/(1-VE) where pdisL2.trial is the value input the trial generation)
#' @param pmv2   A single value. The probability of missing visit 2 or a failed sample at visit 2.
#' @param pmv1   A single value. The probability of a failed sample at visit 1.
#' @return dMat:      A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are infection HPV16/18 status (w/ -9 for missing)
#' @return dMat16:    A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are HPV 16 infection status (w/ -9 for missing)
#' @return dMat18:    A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are HPV 18 infection status (w/ -9 for missing)
#' @return sMat:      A matrix (nrow=number.of.subjects, ncol=number.of.visits). Entries are of sexual activity status (no missing).
#' @examples
#' survRes <- simSurvey(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",yr=5,nsub=5000,pdisL2=0.03,pmv2=0.17,seed=NA)

simSurvey <- function(ofile=NA,yr=5,nsub=5000,pdisL2=0.03,pmv2=0.17,pmv1=2/3*pmv2,seed=NA){

  if (!is.na(seed)) set.seed(seed)
  if (!file.exists(ofile)) stop('Error: ofile does not exist; have you run calcParam?')
  load(ofile)

  ###Recalibrate infection probabilities
  ###Initially they were calculated based on the best assumptions using calcParam
  ###Now, we ensure that pdis is true

  if (yr == 4) {prob    <- infProb
  poss    <- infPoss
  nnew    <- nnew.4yr
  infPers <- infPers.4yr
  svisit  <- sexVisit
  nc      <- 9}

  if (yr == 5) {prob    <- infProb11
  poss    <- infPoss11
  nnew    <- nnew.5yr
  infPers <- infPers.5yr
  svisit  <- sexVisit11
  nc      <- 11}

  ###Recalibrate the initial disease estimates to correspond with pdisL2
  infPers <- ifelse(poss[,nc-1] == 1 & poss[,nc] == 1 & poss[,1] == 0,1,0)
  pNewInf <- sum(prob[infPers==1])

  if (pdisL2 < pNewInf){
    pTerm     <- 1
    while (sum(pTerm^nnew[infPers==1]*prob[infPers==1])>pdisL2) pTerm=pTerm-0.001
    prob3 <- prob2 <- prob*pTerm^nnew
    dif   <-  1-sum(prob2)
    prob3[rowSums(poss)==0] <- prob2[rowSums(poss)==0]   + dif
  }

  if (pdisL2 > pNewInf){
    pTerm     <- 1
    while (sum(pTerm^nnew[infPers==1]*prob[infPers==1])<pdisL2) pTerm=pTerm+0.001
    prob3 <- prob2 <- prob*pTerm^nnew
    dif   <-  sum(prob2)-1
    prob3[rowSums(poss)==0] <- prob2[rowSums(poss)==0]   - dif
  }

  ###Check to make sure infection probability matches pdis
  ###c(sum(prob3[infPers==1]),pdisL2)

  ###Identify the true disease mat of these subjects
  dType <- sample(c(1:nrow(poss)),nsub,replace=TRUE,prob=prob3)
  dMat <- poss[dType,(nc-1):nc]

  ###Identify the sexual activity status
  sMat2 <- matrix(0,nrow=nsub,ncol=ncol(poss))
  svisit <- ifelse(!is.na(svisit),svisit,1/ncol(poss))
  for (i in 1:nsub) { firstVal <- sample(c(1:(ncol(poss)+1)),1,prob=svisit[dType[i],])
                      if (firstVal <= ncol(poss)) sMat2[i,c(firstVal:ncol(poss))] <- 1 }
  sMat  <- sMat2[,(nc-1):nc]

  ###Matrix showing missed visits
  mMat    <- matrix(nrow=nsub,ncol=2)
  mMat[,1] <- rbinom(nsub,1,pmv1)
  mMat[,2] <- rbinom(nsub,1,pmv2)
  dMat[mMat==1] <- -9
  sMat[mMat==1] <- -9
  dMat16 <- dMat18 <- dMat

  ###Want separate matrices for HPV16 and HPV18
  infWomen <- c(1:nrow(dMat))[rowSums(dMat==1) > 0]
  inf16    <- sample(infWomen,round(0.67*length(infWomen)))
  inf18    <- infWomen[!is.element(infWomen,inf16)]
  wMat     <- matrix(c(1:nsub),nrow=nsub,ncol=ncol(dMat),byrow=FALSE)
  dMat16[!is.element(wMat,inf16) & dMat==1] <- 0
  dMat18[!is.element(wMat,inf18) & dMat==1] <- 0
  return(list("dMat"=dMat,"sMat"=sMat,"dMat16"=dMat16,"dMat18"=dMat18))
}


###################################################################################################
###################################################################################################

#' Combine survey and trial participants
#'
#' Combine survey and trial participants
#' @export
#' @param trialRes   An object. The output from simTrialArm.
#' @param surveyRes  An object. The output from simSurvey
#' @return data      A matrix. The columns are GROUP (0=survey, 1=trial); DIS1.16 and DIS2.16 (HPV 16 status at 2 key visits);
#' DIS1.18 and DIS2.18 (HPV 18 status at 2 key visits); DIS0.16 (HPV 16 status at baseline; missing for survey);
#' DIS0.18 (HPV18 status at baseline; missing for survey); SEX1 and SEX2 (sexual activity status at two key visits);
#' X1 and X2 (two covariates); ATP (ATP status for trial participants; All survey participants have a value of 1);
#' @examples
#' trial1dose <- simTrialArm(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",yr=5,nsub=5000,pdisL2=0.04*(1-0.8),pdo6=0.05,pdo=0.02,pmv=0.17)
#' survey     <- simSurvey(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",  yr=5,nsub=5000,pdisL2=0.04,pmv2=0.17,pmv1=0.10)
#' study      <- simStudy(trialRes=trial1dose,surveyRes=survey)
simStudy <- function(trialRes,surveyRes){
     nc    <-   ncol(trialRes$dMat)
     GROUP <-   c(rep(1,nrow(trialRes$dMat)),rep(0,nrow(surveyRes$dMat)))
     SEX1  <-   c(trialRes$sMat[,nc-1],surveyRes$sMat[,1])
     SEX2  <-   c(trialRes$sMat[,nc],surveyRes$sMat[,2])
     DIS1.16  <-   c(trialRes$dMat16[,nc-1],surveyRes$dMat16[,1])
     DIS2.16  <-   c(trialRes$dMat16[,nc],surveyRes$dMat16[,2])
     DIS0.16  <-   c(trialRes$dMat16[,1],rep(-9,nrow(surveyRes$dMat16)))
     DIS1.18  <-   c(trialRes$dMat18[,nc-1],surveyRes$dMat18[,1])
     DIS2.18  <-   c(trialRes$dMat18[,nc],surveyRes$dMat18[,2])
     DIS0.18  <-   c(trialRes$dMat18[,1],rep(-9,nrow(surveyRes$dMat18)))
     X1    <-   rbinom(nrow(trialRes$dMat) + nrow(surveyRes$dMat),1,0.5)
     X2    <-   rbinom(nrow(trialRes$dMat) + nrow(surveyRes$dMat),1,0.5)
     ATP      <- c(trialRes$ATPx,rep(1,length=nrow(surveyRes$dMat)))
     veDATA <- as.matrix(data.frame(GROUP,SEX1,SEX2,DIS1.16,DIS2.16,DIS0.16,DIS1.18,DIS2.18,DIS0.18,ATP,X1,X2))
veDATA
}

###################################################################################################
###################################################################################################

#' Estimate the VE
#'
#' Estimate the VE
#' @export
#' @param study              A matrix. The columns GROUP (0=survey, 1=trial); DIS1.16 and DIS2.16 (HPV 16 status at 2 visits);
#' DIS1.18 and DIS2.18 (HPV 18 status at 2 visits); DIS0.16 (HPV 16 status at baseline; missing for survey);
#' DIS0.18 (HPV18 status at baseline; missing for survey); SEX1 and SEX2 (sexual activity status at two visits);
#' ATP (indicator for ATP status, all survey individuals must equal 1) covariates for the propensity score model;
#' missing values should be set to -9
#' @param psCov               A vector. A vector containing the variable names for the covariates to be used in the propensity score model.
#' @return VE                 A single value. Estimated VE
#' @return ntrial.obs.L2      A single value. Number of trial particpants with the same type of infection at the last vists
#' @return ntrial.obs.L2pB    A single value. Number of trial particpants with the same type of infection at the last vists and at baseline
#' @return nsurvey.obs.L2     A single value. Number of survey particpants with the same type of infection at both vists
#' @return Ntrial             A single value. Number of individuals in the ATP cohort who completed at least one of the last two visits
#' @return Nsurvey            A single value. Number of individuals in the survey who completed at least one of the two visits
#' @return Ntrial.both        A single value. Number of individuals in the ATP cohort who completed both of the last two visits
#' @return Nsurvey.both       A single value. Number of individuals in the survey who completed both of the two visits
#' @return ntrial.est.L2      A single value. Estimated number of trial particpants with the same type of infection at the last vists
#' @return ntrial.est.L2pB    A single value. Estimated number of trial particpants with the same type of infection at the last vists and at baseline
#' @return nsurvey.west.L2    A single value. Weighted estimated number of survey particpants with the same type of infection at the last vists
#' @examples
#' trial1dose <- simTrialArm(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",yr=5,nsub=5000,pdisL2=0.04*(1-0.8),pdo6=0.05,pdo=0.02,pmv=0.17)
#' survey     <- simSurvey(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",  yr=5,nsub=5000,pdisL2=0.04,pmv2=0.17,pmv1=0.10)
#' study      <- simStudy(trialRes=trial1dose,surveyRes=survey)
#' estVE(study,psCov=c("X1","X2"))

estVE <- function(veDATA,psCov=NA){
  require(tidyverse)
  veDATA2           <- as.data.frame(veDATA)
  colnames(veDATA2) <- colnames(veDATA)
  if (sum(is.na(veDATA2$GROUP) | !is.element(veDATA2$GROUP,c(0,1)))>0) stop("ERROR: All GROUP values are not 0 or 1")

  ###We will first define "naive events" (aka observed events)
  ###Missing unless we observed both visits for at last one HPV type

  veDATA2                <-   veDATA2 %>% mutate(nevent16=NA,nevent18=NA)
  veDATA2[,"nevent18"]   <-   ifelse(veDATA2[,"DIS1.18"]== 1 & veDATA2[,"DIS2.18"]==1,1,
                                   ifelse(veDATA2[,"DIS1.18"]!=-9 & veDATA2[,"DIS2.18"]!=-9,0,NA))
  veDATA2[,"nevent16"]   <-   ifelse(veDATA2[,"DIS1.16"]== 1 & veDATA2[,"DIS2.16"]==1,1,
                                ifelse(veDATA2[,"DIS1.16"]!=-9 & veDATA2[,"DIS2.16"]!=-9,0,NA))
  veDATA2[,"nevent"]     <-   ifelse( (veDATA2[,"nevent16"]== 1 & !is.na(veDATA2[,"nevent16"])) |
                                      (veDATA2[,"nevent18"]== 1 & !is.na(veDATA2[,"nevent18"])),1,
                                 ifelse( (veDATA2[,"nevent16"]== 0 & !is.na(veDATA2[,"nevent16"])) |
                                         (veDATA2[,"nevent18"]== 0 & !is.na(veDATA2[,"nevent18"])),0,NA))

  ###Change -9 to Missing Values
  for (i in 1:ncol(veDATA)) veDATA2[,i] <- ifelse(veDATA2[,i]==-9,NA,veDATA2[,i])



  ###For a given HPV type, create a matrix with two columns
  ###First column is the expected probability disease at the last two visits
  ###First column is the expected probability disease at the last two visits AND disease at the first visit

  expTrial <- function(veDATA2,type=16){

    ###Return the results in the same orderr
    veDATA2[,"ordInternal"] <- c(1:nrow(veDATA2))

    ###Get disease of appropriate HPV type
    if (type == 16) { veDATA2[,"DIS0"] <- veDATA2[,"DIS0.16"]
                      veDATA2[,"DIS1"] <- veDATA2[,"DIS1.16"]
                      veDATA2[,"DIS2"] <- veDATA2[,"DIS2.16"] }
    if (type == 18) { veDATA2[,"DIS0"] <- veDATA2[,"DIS0.18"]
                      veDATA2[,"DIS1"] <- veDATA2[,"DIS1.18"]
                      veDATA2[,"DIS2"] <- veDATA2[,"DIS2.18"] }

    ###Set all non-ATP individuals to have missing values
    veDATA2 <- veDATA2 %>% mutate(DIS0 = ifelse(ATP==0,NA,DIS0),DIS1 = ifelse(ATP==0,NA,DIS1),DIS2 = ifelse(ATP==0,NA,DIS2))

    ###DIVIDE RESULTS
    trial   <- veDATA2 %>% filter(GROUP==1)
    survey  <- veDATA2 %>% filter(GROUP==0)

    ###TRIAL RESULTS
    ###GET THE AVERAGES (FOR LATER IMPUTATION for PEOPLE WITH MISSING VALUES)
    ###(b0 = basline uninfected; m1 = missing the first of the relevant visits; m2 = missing 2nd of relevant visits)
    b0.m1   <- trial %>% filter(DIS0==0 & !is.na(DIS0)) %>% filter(DIS2==1 & !is.na(DIS2)) %>% summarize(mean=mean(DIS1,na.rm=TRUE))
    b0.m2   <- trial %>% filter(DIS0==0 & !is.na(DIS0)) %>% filter(DIS1==1 & !is.na(DIS1)) %>% summarize(mean=mean(DIS2,na.rm=TRUE))
    if (is.na(b0.m1)) b0.m1 <- b0.m2
    if (is.na(b0.m2)) b0.m2 <- b0.m1
    if(b0.m1==1 & !is.na(b0.m1)) b0.m1 <- 0.99
    if(b0.m2==1 & !is.na(b0.m2)) b0.m2 <- 0.99

    b1.m1   <- trial %>% filter(DIS0==1 & !is.na(DIS0)) %>% filter(DIS2==1 & !is.na(DIS2)) %>% summarize(mean=mean(DIS1,na.rm=TRUE))
    b1.m2   <- trial %>% filter(DIS0==1 & !is.na(DIS0)) %>% filter(DIS1==1 & !is.na(DIS1)) %>% summarize(mean=mean(DIS2,na.rm=TRUE))
    if (is.na(b1.m1)) b1.m1 = b1.m1
    if (is.na(b1.m2)) b1.m2 = b1.m2
    if(b1.m1==1 & !is.na(b1.m1)) b1.m1 <- 0.99
    if(b1.m2==1 & !is.na(b1.m2)) b1.m2 <- 0.99

    ###FILL in expected values for each trial participant
    anul <- function(t) as.numeric(unlist(t))
    event <- matrix(nrow=nrow(trial),ncol=1)
    event <- ifelse( (!is.na(trial$DIS1) & trial$DIS1==1) &
                       (!is.na(trial$DIS2) & trial$DIS2==1),1,event)
    event <- ifelse( (!is.na(trial$DIS1) & trial$DIS1==0) |
                       (!is.na(trial$DIS2) & trial$DIS2==0),0,event)
    event <- ifelse( (!is.na(trial$DIS1) & trial$DIS1==1) &
                       ( is.na(trial$DIS2) & trial$DIS0==0),anul(b0.m2),event)
    event <- ifelse( (!is.na(trial$DIS1) & trial$DIS1==1) &
                       ( is.na(trial$DIS2) & trial$DIS0==1),anul(b1.m2),event)
    event <- ifelse( (!is.na(trial$DIS2) & trial$DIS2==1) &
                       ( is.na(trial$DIS1) & trial$DIS0==0),anul(b0.m1),event)
    event <- ifelse( (!is.na(trial$DIS2) & trial$DIS2==1) &
                       ( is.na(trial$DIS1) & trial$DIS0==1),anul(b1.m1),event)

    ###survey RESULTS


    ###GET THE AVERAGES (FOR LATER IMPUTATION for PEOPLE WITH MISSING VALUES)
    m1   <- survey %>% filter(DIS2==1 & !is.na(DIS2)) %>%  summarize(mean=mean(DIS1,na.rm=TRUE))
    m2   <- survey %>% filter(DIS1==1 & !is.na(DIS1)) %>%  summarize(mean=mean(DIS2,na.rm=TRUE))
    if (is.na(m1)) m1 <- m2
    if (is.na(m2)) m2 <- m1
    if(m1==1 & !is.na(m1)) m1 <- 0.99
    if(m2==1 & !is.na(m2)) m2 <- 0.99

    event.s <- matrix(nrow=nrow(survey),ncol=1)
    event.s <- ifelse( (!is.na(survey$DIS1) & survey$DIS1==1) &
                         (!is.na(survey$DIS2) & survey$DIS2==1),1,event.s)
    event.s <- ifelse( (!is.na(survey$DIS1) & survey$DIS1==0) |
                         (!is.na(survey$DIS2) & survey$DIS2==0),0,event.s)
    event.s <- ifelse( (!is.na(survey$DIS1) & survey$DIS1==1) &
                         is.na(survey$DIS2),anul(m2),event.s)
    event.s <- ifelse( (!is.na(survey$DIS2) & survey$DIS2==1) &
                         is.na(survey$DIS1),anul(m1),event.s)

    trial2   <- trial  %>% mutate(event=event,   eventx=event*DIS0)
    survey2  <- survey %>% mutate(event=event.s, eventx=event.s*DIS0)
    veDATA3  <- trial2 %>% bind_rows(survey2) %>% arrange(ordInternal) %>% select(event,eventx)
    veDATA3
  }
  event16 <- expTrial(veDATA2,16)
  event18 <- expTrial(veDATA2,18)

  ###Clean-up results
  event16[,"event"] <- ifelse(is.na(event16[,"event"]) & !event18[,"event"],0,event16[,"event"])
  event18[,"event"] <- ifelse(is.na(event18[,"event"]) & !event16[,"event"],0,event18[,"event"])
  event16[,"eventx"] <- ifelse(is.na(event16[,"eventx"]) & !event18[,"eventx"],0,event16[,"eventx"])
  event18[,"eventx"] <- ifelse(is.na(event18[,"eventx"]) & !event16[,"eventx"],0,event18[,"eventx"])

  event16v  <- as.numeric(unlist(event16[,"event"]))
  event18v  <- as.numeric(unlist(event18[,"event"]))
  event16xv <- as.numeric(unlist(event16[,"eventx"]))
  event18xv <- as.numeric(unlist(event18[,"eventx"]))

  ###Get the group variable
  group     <- as.numeric(unlist(veDATA2$GROUP))

  ###Combine results (i.e. HPV16 or HPV18 event)
  event     <- 1 - (1-event16v)*(1-event18v)
  eventx    <- 1 - (1-event16xv)*(1-event18xv)

  ###Add variables
  veDATA3 <- veDATA2 %>% mutate(event=event,eventx=eventx) %>% filter(!is.na(event))

  ###Get the propensities for being in the trial
  fv       <- glm(GROUP~.,family=binomial(link="logit"),veDATA3[,c("GROUP",psCov)])$fitted.values
  veDATA4  <- veDATA3 %>% mutate(wevent = event*fv/(1-fv), weight = fv/(1-fv))

  ###REPORT EXPECTED VALUES
  nev1     <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 1) %>% summarize(nev1=sum(event))))
  nev1x    <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 1) %>% summarize(nev1x=sum(eventx))))
  nev0     <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 0) %>% summarize(nev0=sum(wevent))))
  pop0     <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 0) %>% summarize(pop0=sum(weight))))
  VE      <- 1-(nev1-nev1x)/(nev0-nev1x)


  ###REPORT OBSERVED VALUES
  ###(note: only true events have event == 1)
  nev1.obs     <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 1) %>% summarize(nev1.obs=sum(event==1))))
  nev1x.obs    <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 1) %>% summarize(nev1x.obs=sum(eventx==1))))
  nev0.obs     <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 0) %>% summarize(nev0.obs=sum(event==1))))

  nATP        <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 1) %>% summarize(nATP=sum(ATP==1))))
  nSURVEY     <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 0) %>% summarize(nSURVEY = length(GROUP) )))

  ###Only individuals with disease status at both visits have non-missing values
  nATP.both   <- as.numeric(unlist(veDATA4 %>% filter(GROUP == 1 & !is.na(nevent)) %>% summarize(nATP=sum(ATP==1))))
  nSURVEY.both<- as.numeric(unlist(veDATA4 %>% filter(GROUP == 0 & !is.na(nevent)) %>% summarize(nSURVEY = length(GROUP) )))
  return(list("VE"=VE,       "ntrial.obs.L2"=nev1.obs,"ntrial.obs.L2pB"=nev1x.obs,"nsurvey.obs.L2"=nev0.obs,
              "Ntrial"=nATP, "Nsurvey"=nSURVEY, "Ntrial.both"=nATP.both, "Nsurvey.both"=nSURVEY.both,
              "ntrial.est.L2"=nev1,"ntrial.est.L2pB"=nev1x,"nsurvey.west.L2"=nev0))
}



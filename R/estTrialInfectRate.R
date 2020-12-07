### All notation is based on (https://pubmed.ncbi.nlm.nih.gov/29474934/)

################################################################################################

#' Calculate the probability of acquiring a new infection in 6 months
#'
#' Calculate the probability of acquiring a new infection in 6 months
#' @returns The estimated probability of acquiring a new infection 6 months after the last visit


calcProbNewInf <- function(){

  ### First, we define the rate of clearance following an infection.
  lambda=-2*log(.666)

  ### Second, for this value of lambda and an infection incidence rate k,
  ### we can determine the prevalence of new infections after a given time T

  prev<-function(k,lambda,T)
  {
      work=k*exp(-T*lambda)
      work=work/(lambda-k)
      work=work*(exp(T*(lambda-k))-1)
      return(work)
  }

  ###Third, we want to find the value of k that gives us a prevalence of 0.073
  ###after one year

  fsolve<-function(k,prevrate,lambda,T)
  {
      work=prevrate-prev(k,lambda,T)
      return(work)
  }
  k=uniroot(fsolve,c(0,1),p=0.073,lambda,1)$root

  ###Fourth, we use this value of lambda, k, and T = 0.5, to determine the prevalence
  ###of new infections after one year

  p.5 <- prev(k,lambda,0.5)

return(p.5)
}

################################################################################################

#' Calculate the key parameters for ESCUDDO
#'
#' Calculate the key parameters for ESCUDDO
#' @export
#' @param pSA   A vector of 11 values. Probability of starting sexual activity by ages 12:22
#' @param pEA   A vector of 5 values.  Proportion of subjects who enroll at ages 12:16
#' @param hI1   A vector of 9 values.  Herd Immunity for individuals in 2018:2026 (<= 18 yo )
#' @param hI2   A vector of 9 values.  Herd Immunity for individuals in 2018:2026 (>  18 yo )
#' @param npop  A single value. Total number of women in the simulated population
#' @param ofile A string. The file name of where to store key parameters for simulation.
#' @return p0   A value. The proportion of women who have an infection at baseline
#' @return pi0.amongBaselineUnifected.4yr A value. Among baseline (0/6 month) uninfected women, this value  is the proportion who have an incident persistent infection in four years
#' @return pi0.amongBaselineUnifected.5yr A value. Among baseline (0/6 month) uninfected women, this value  is the proportion who have an incident persistent infection in five years
#' @return delta0.4yr A value. Among baseline-0 (only 0 month) uninfected women, this value  is the proportion who have a persistent infection spanning 3.5 and 4.0  years
#' @return delta0.5yr A value. Among baseline-0 (only 0 month) uninfected women, this value  is the proportion who have a persistent infection spanning 4.5 and 5.0  years
#' @examples
#' myRes <- calcParam(
#' hI1     = rep(0,9),
#' hI2     = rep(0,9),
#' npop    = 5000000,
#' ofile   = "/volumes/data/Projects/HPV/out/Rparams.Rdat")

calcParam <-  function(pSA     = c(0,    0,   0.10,0.26,0.39,0.53,0.66,0.77,0.84,0.91,0.95),
                       pEA     = c(0.25,0.20,0.19,0.19,0.17),
                       hI1     = c(0,0,0.03,0.04,0.08,0.10,0.10,0.10,0.10),
                       hI2     = c(0,0,0.03,0.03,0.05,0.07,0.07,0.07,0.07),
                       npop    = 500000,
                       ofile   = NA){

  ###Probability that an uninfected sexually active woman
  ###acquires an HPV 16 or 18 infection six months later
  pAcquire  <- calcProbNewInf()

  ###Age range for consideration
  age       <- c(12:22)

  ###Probability of starting sexual activity in a given year
  pStart    <- pSA-c(0,pSA[-length(pSA)])

  ###percentage of each starting age
  enroll     <- matrix(nrow=5,ncol=2)
  enroll[,1] <- c(12:16)
  enroll[,2] <- pEA

  ###herd immunity (calendar year, <= 18 y0, > 18 yo)
  herd       <- matrix(nrow=9,ncol=3)
  herd[,1]   <- c(2018:2026)
  herd[,2]   <- hI1
  herd[,3]   <- hI2

  ###Probability that an infection persists for 6 months
  pPersist   <- 0.66

  ###We will have 9 or 11 visits (i.e. 2 visits x 4 years + baseline or 2 visits/year x 5 years + baseline)
  ###We consider all positive infection states at those 9 or 11 visits

  ###all possible infection scenarios over 9 visits
  infPoss   <- rep(rep(c(0,1),each=2^(9-1)),2^(1-1))
  for (i in 2:9) infPoss <- cbind(infPoss,rep(rep(c(0,1),each=2^(9-i)),2^(i-1) ))

  ###all possible infection scenarios over 11 visits
  infPoss11 <- rep(rep(c(0,1),each=2^(11-1)),2^(1-1))
  for (i in 2:11) infPoss11 <- cbind(infPoss11,rep(rep(c(0,1),each=2^(11-i)),2^(i-1) ))

  ###get the probability of each of the infection patterns
  infProbx       <- rep(0,512)
  infProb11x     <- rep(0,2048)

  ###get the distribution of first sexual visit.
  ###(i.e. probability that a woman has had her sexual debut prior to each visit
  ###given the infection status)
  sexVisitx       <- matrix(0,nrow=512, ncol=10)
  sexVisit11x     <- matrix(0,nrow=2048,ncol=12)

   ###We will ideally want a simulation of 50 million women
   ###However, we will divide 10 simulations of 5 million womne
   nruns <- 10
   runMat <- matrix(nrow=nruns,ncol=5)
   ###p0:        Proportion of women who have an infection at visit 0.
   ###pi0.4y:    Proportion of baseline (i.e. visits 0/6) uninfected women who will have a persistent infect during visits 3-9
   ###pi0.5y:    Proportion of baseline (i.e. visits 0/6) uninfected women who will have a persistent infect during visits 3-11
   ###delta0.4y: Proportion of visit-0 uninfected women who will have a persistent infection at visit 8/9
   ###delta0.5y: Proportion of visit-0 uninfected women who will have a persistent infection at visit 10/11
   colnames(runMat) <- c("p0","pi0.4yr","pi0.5yr","delta0.4yr","delta0.5yr")

      for (tr in 1:nruns){

            ###We run simulations with 5 million women
            nwomen      <- floor(npop/10)

            ###age at study start
            startingAge <- sample(c(12:16),nwomen,replace=TRUE,prob=pEA)

            ###age at sexual debut
            sexAge      <- sample(c(age,23),nwomen,prob=c(pStart,1-sum(pStart)),replace=T)

            ###check to make sure distribution of sexual age looks right
            #prop <- table(sexAge)/nwomen
            #want <- c(pStart[-c(1:2)],1-sum(pStart))
            #cbind(prop,want)

            ###calendar year of study start
            calYear    <- 2018+rbinom(nwomen,1,0.5)

            ###dMat shows the disease status of each women (row) at each age
            ###(columns;i.e. ages 12, 12.5, 13, 13.5, ..., 19.5,20,20.5,21)
            dMat      <- matrix(nrow=nwomen,ncol=19)
            agePoss   <- seq(12,21,by=0.5)
            dMat[,1]  <- rbinom(nwomen,1,(12 >= sexAge)*pAcquire)
            pCheck    <- c()
            for (i in 2:19) {
                  curYear  <- floor(calYear + (agePoss[i]-startingAge))
                  ###Note i == 13 corresponds to age 18
                  if (i <= 13) pmult2   <- 1-herd[match(curYear,herd[,1]),2]
                  if (i  > 13) pmult2   <- 1-herd[match(curYear,herd[,1]),3]
                  pmult    <- ifelse(is.na(pmult2),1,pmult2)
                  dMat[,i] <- rbinom(nwomen,1,
                                   dMat[,i-1]*pPersist +
                                   (1-dMat[,i-1])*(12+(i-1)*0.5>=sexAge)*pAcquire*pmult)
            }

            ###sMat shows the sexual activity status (1 indicates sexually active) of each girl
            ###(columns;i.e. ages 12, 12.5, 13, 13.5, ..., 19.5,20,20.5,21)
            sMat      <- matrix(nrow=nwomen,ncol=ncol(dMat))
            for (i in 1:ncol(dMat)) sMat[,i] <- ifelse(sexAge > 12 + (i-1)/2,0,1)

            ###number of women with a disease while not sexually active
            ###sum(dMat[sMat==0])

            ###dMat2 shows the disease status of each women (row) at each study visit
            ###sMat2 shows the sexual status of each women (row) at each study visit
            ###(columns; visits 1, 2, ..., 8, 9, 10, 11)
            dMat2 <- sMat2 <- matrix(nrow=nwomen,ncol=11)

            for (i in 1:nwomen) {
                        fc        <- 2*(startingAge[i] - 12)+1
                        dMat2[i,] <- dMat[i,fc:(fc+10)]
                        sMat2[i,] <- sMat[i,fc:(fc+10)]
                        }

            for (i in 1:nwomen) {
                     poi <- c(1:512)[infPoss%*%matrix(dMat2[i,1:9],ncol=1)+
                                    (1-infPoss)%*%(1-matrix(dMat2[i,1:9],ncol=1))==9]
                     infProbx[poi] <- infProbx[poi] +1
                     sexVisitx[poi,min(c(1:9)[sMat2[i,1:9]==1],10)] = sexVisitx[poi,min(c(1:9)[sMat2[i,1:9]==1],10)] + 1

                     poi11 <- c(1:2048)[infPoss11%*%matrix(dMat2[i,],ncol=1)+
                                        (1-infPoss11)%*%(1-matrix(dMat2[i,],ncol=1))==11]
                     infProb11x[poi11] <- infProb11x[poi11] +1
                     sexVisit11x[poi,min(c(1:11)[sMat2[i,1:11]==1],12)] = sexVisit11x[poi,min(c(1:11)[sMat2[i,1:11]==1],12)] + 1
                     }

            ###indicator of whether a woman has a baseline infection
            baselineInfection  <- ifelse(dMat2[,1]+dMat2[,2]>0,1,0)
            initInf            <- dMat2[,1]
            ###Note: Baseline infections are about to disappear for young women...careful


            ###for purposes pcr-dectable infection, a woman cannot have an infection prior to age 15
            for (i in 1:nwomen) if (startingAge[i] <= 14) {
                                               lastUT <- 6 - (startingAge[i]-12)*2
                                               dMat2[i,1:lastUT] <- 0
            }

          ###Infection status at next (non-missing) study visit
          nextNM <- matrix(0,nrow=nrow(dMat2),ncol=ncol(dMat2))
          for (i in 10:1) nextNM[,i] <- ifelse(!is.na(dMat2[,i+1]),dMat2[,i+1],nextNM[,i+1])

          ###Indicator for a persistent infection
          pMat <- ifelse(dMat2==1 & nextNM==1 & !is.na(dMat2) & !is.na(nextNM),1,0)

          ###check to make sure the probability of a persistent infection is 0.666
          ###mean(pMat[,1:10][dMat2[,1:10]==1 & !is.na(dMat2[,1:10])],na.rm=T)

          ###proportion of uninfected women who will have an incident persistent infection during the trial
          nipa.4yr      <- 0
          nipa.5yr      <- 0

          for (i in 1:nrow(pMat)) nipa.4yr <- nipa.4yr + max(pMat[i,3:8])*(1-baselineInfection[i])
          runMat[tr,1] <- nipa.4yr/sum(1-baselineInfection)

          for (i in 1:nrow(pMat)) nipa.5yr <- nipa.5yr  + max(pMat[i,3:10])*(1-baselineInfection[i])
          runMat[tr,4] <- nipa.5yr/sum(1-baselineInfection)

          runMat[tr,2] <- mean(initInf)

          runMat[tr,3] <- mean(pMat[initInf==0,8] ==1)
          runMat[tr,5] <- mean(pMat[initInf==0,10]==1)

          ###print(tr)
          ####for (tr in 1:nruns){
          }

        etir.p0           <- mean(runMat[,2])
        etir.pi0.4yr      <- mean(runMat[,1])
        etir.pi0.5yr      <- mean(runMat[,4])
        etir.d0.4yr       <- mean(runMat[,3])
        etir.d0.5yr       <- mean(runMat[,5])

        ####NUMBER OF NEW INFECTIONS
        nnew.4yr          <- infPoss[,1]
        nnew.5yr          <- infPoss11[,1]
        for (i in 2:9)   nnew.4yr    <- nnew.4yr + ifelse(infPoss[,i-1]  ==0 & infPoss[,i]  ==1,1,0)
        for (i in 2:11)  nnew.5yr    <- nnew.5yr + ifelse(infPoss11[,i-1]==0 & infPoss11[,i]==1,1,0)

        ####A set with no baseline infection and new persistent infection
        infPers.4yr       <- ifelse(rowSums(infPoss[,3:8]*infPoss[,4:9])       > 0 & rowSums(infPoss[,1:2])  ==0,1,0)
        infPers.5yr       <- ifelse(rowSums(infPoss11[,3:10]*infPoss11[,4:11]) > 0 & rowSums(infPoss11[,1:2])==0,1,0)

        ###The probability of each infection combination
        infProb         <- infProbx/sum(infProbx)
        infProb11       <- infProb11x/sum(infProb11x)

        ###Distribution of First Sexual Visit
        sexVisit   <- sexVisitx/  matrix(rowSums(sexVisitx),  nrow=nrow(sexVisitx),  ncol=ncol(sexVisitx),  byrow=FALSE)
        sexVisit11 <- sexVisit11x/matrix(rowSums(sexVisit11x),nrow=nrow(sexVisit11x),ncol=ncol(sexVisit11x),byrow=FALSE)



    if (!is.na(ofile)) save(infPoss,infProb,infPoss11,infProb11,
                  infPers.4yr,infPers.5yr,
                  runMat,nnew.4yr,nnew.5yr,
                  etir.p0,
                  etir.pi0.4yr,etir.d0.4yr,
                  etir.pi0.5yr,etir.d0.5yr,
                  sexVisit,sexVisit11,
                  file=ofile)


        list("p0"=etir.p0,
          "pi0.amongBaselineUnifected.4yr"=etir.pi0.4yr,
          "pi0.amongBaselineUnifected.5yr"=etir.pi0.5yr,
          "delta0.4yr"=etir.d0.4yr,
          "delta0.5yr"=etir.d0.5yr)
}

#######################################################################################################################################
#######################################################################################################################################




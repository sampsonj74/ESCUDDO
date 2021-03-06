% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estVE.R
\name{estVE}
\alias{estVE}
\title{Estimate the VE}
\usage{
estVE(veDATA, psCov = NA)
}
\arguments{
\item{psCov}{A vector. A vector containing the variable names for the covariates to be used in the propensity score model.}

\item{study}{A matrix. The columns GROUP (0=survey, 1=trial); DIS1.16 and DIS2.16 (HPV 16 status at 2 visits);
DIS1.18 and DIS2.18 (HPV 18 status at 2 visits); DIS0.16 (HPV 16 status at baseline; missing for survey);
DIS0.18 (HPV18 status at baseline; missing for survey); SEX1 and SEX2 (sexual activity status at two visits);
ATP (indicator for ATP status, all survey individuals must equal 1) covariates for the propensity score model;
missing values should be set to -9}
}
\value{
VE                 A single value. Estimated VE

ntrial.obs.L2      A single value. Number of trial particpants with the same type of infection at the last vists

ntrial.obs.L2pB    A single value. Number of trial particpants with the same type of infection at the last vists and at baseline

nsurvey.obs.L2     A single value. Number of survey particpants with the same type of infection at both vists

Ntrial             A single value. Number of individuals in the ATP cohort who completed at least one of the last two visits

Nsurvey            A single value. Number of individuals in the survey who completed at least one of the two visits

Ntrial.both        A single value. Number of individuals in the ATP cohort who completed both of the last two visits

Nsurvey.both       A single value. Number of individuals in the survey who completed both of the two visits

ntrial.est.L2      A single value. Estimated number of trial particpants with the same type of infection at the last vists

ntrial.est.L2pB    A single value. Estimated number of trial particpants with the same type of infection at the last vists and at baseline

nsurvey.west.L2    A single value. Weighted estimated number of survey particpants with the same type of infection at the last vists
}
\description{
Estimate the VE
}
\examples{
trial1dose <- simTrialArm(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",yr=5,nsub=5000,pdisL2=0.04*(1-0.8),pdo6=0.05,pdo=0.02,pmv=0.17)
survey     <- simSurvey(ofile="/volumes/data/Projects/HPV/out/Rparams.Rdat",  yr=5,nsub=5000,pdisL2=0.04,pmv2=0.17,pmv1=0.10)
study      <- simStudy(trialRes=trial1dose,surveyRes=survey)
estVE(study,psCov=c("X1","X2"))
}

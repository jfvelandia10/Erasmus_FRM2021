#' ---
#' title: "Example Code Estimation EVPI and EVPPI"
#' author: "Juan Velandia"
#' date: "August 18, 2021"
#' ---

##Load required libraries ------------------------------------------------------
library(BCEA)
install.packages("INLA",repos=c(getOption("repos"),
      INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA) 
library(splancs)
library(GrassmannOptim)
#Download ldr package from here as tar.gz
#https://cran.r-project.org/src/contrib/Archive/ldr/ and install
library(ldr) 
library(devtools)
install_github("voi") #Install VOI package
library("voi")

##============================Implementation====================================
#NB Correspond to the net benefits or utility. 
#In this case is the difference between benefits and costs
#Each row is an iteration and each column is each one of the alternatives.
NB<-totalBenefitsAll-totalCostsAll
nSims<-nrow(NB)

#Estimate EVPI and EVPPI
#EVPI - Using the results of the PA ============================================
#Selected option under current knowledge - maximum utility 
maxUtilityCK<-max(colMeans(NB))
selectOption<-which(colMeans(NB)==maxUtilityCK)

#Calculate max utility, selected option, opportunity loss and VOI per iteration
VOIResults<-matrix(0, nrow = nSims, ncol = 4)

for(i in 1:nSims){
  VOIResults[i,1]<-max(NB[i,])
  VOIResults[i,2]<-which(NB[i,]==VOIResults[i,1])
  
  VOIResults[i,3]<-VOIResults[i,1]-NB[i,selectOption]
  if(VOIResults[i,3]<0){
    VOIResults[i,3]<-0
  }
  
  VOIResults[i,4]<-VOIResults[i,1]-maxUtilityCK
}
colnames(VOIResults)<-c('Max_Util_Scenario','Sel_Option','OL','VOI')

#Calculate EVPI
EVPIMeth1<-mean(VOIResults[,1])-maxUtilityCK

#EVPPI - It depends on the combination of sources of uncertainty ===============
#Using the VOI package 
#(Read more on https://chjackson.github.io/voi/articles/voi.html)

#Corresponds to the parameters of the PA. 
#Each row is an iteration and each column a parameter
UncParameterInput

#Matrix of EVPPIs - 
#Estimate hazard, exposure, vulnerability, cost and combination between them
matrixEVPPIS<-matrix(0,nrow=4,ncol=4)
listParMat<-list(c(1:7),c(8:12),c(13:16),c(17:18))

#Function evppi receives;
#-NB
#-Parameter uncertainty matrix
#-Name of the columns of parameters to be analyzed
#-Regression method (reg method)

#In this case reg method: GAM and SO for single parameters, INLA for more than 2

#The returned object allows to check the performance of the fit.
matrixEVPPIS[1,1]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[1]])],verbose=T)$evppi

matrixEVPPIS[1,2]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[1]],listParMat[[2]])],
          verbose=T)$evppi

matrixEVPPIS[1,3]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[1]],listParMat[[3]])],
          verbose=T)$evppi

matrixEVPPIS[1,4]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[1]],listParMat[[4]])],
          verbose=T)$evppi


matrixEVPPIS[2,2]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[2]])],verbose=T)$evppi

matrixEVPPIS[2,3]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[2]],listParMat[[3]])],
          verbose=T)$evppi

matrixEVPPIS[2,4]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[2]],listParMat[[4]])],
          verbose=T)$evppi

matrixEVPPIS[3,3]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[3]])],verbose=T)$evppi

matrixEVPPIS[3,4]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[3]],listParMat[[4]])]
          ,verbose=T)$evppi

matrixEVPPIS[4,4]<-evppi(NB,UncParameterInput,
          pars=colnames(UncParameterInput)[c(listParMat[[4]])],verbose=T)$evppi

colnames(matrixEVPPIS)<-c('Hazard','Exposure','Vulnerability','Costs')
rownames(matrixEVPPIS)<-c('Hazard','Exposure','Vulnerability','Costs')


#Single parameter estimation (SO and GAM) - In this case 18 parameters

singleParamEVPPI<-rep(0,18)
for (i in 1:7){
  singleParamEVPPI[i]<-evppi(NB,UncParameterInput,
                pars=colnames(UncParameterInput)[i],verbose=T)$evppi
}

for (i in 8:18){
  singleParamEVPPI[i]<-evppi(NB,UncParameterInput,
    pars=colnames(UncParameterInput)[i],verbose=T,method="so",n.blocks=5)$evppi
}

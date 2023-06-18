#Load libraries
library(extraDistr)
library(BCEA)
library(INLA) #install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(splancs)
library(GrassmannOptim)
library(ldr) #https://cran.r-project.org/src/contrib/Archive/ldr/
library(devtools)
install_github("voi")
library("voi")
library('immunarch')


##Estimate costs per measure-------------------------------------------------------
#Unitary construction costs are considered a triangular distribution
parUnitCosts<-c(15,22,30)
#Operation and Maintenance Cost O&M
parOMCosts<-c(0.01,0.03,0.05)
#uncertainty Parameter costs


#height dike varies between 0.5 and 5 ---- CHANGE TO NUMBER OF MEASURES
heightDike<-seq(0.5,3.5,0.5)
#Volume dike following La Mojana design and a length of 4.5km
volumeDike<-(3.5*heightDike+1*heightDike^2/2+1*heightDike^2/2)*2500

#Generate simulation prices
costsMatrix<-matrix(, nrow = nSims, ncol = (numMeasures+2))
costsMatrix[,1]<-rtriang(nSims,a=parUnitCosts[1],b=parUnitCosts[3],c=parUnitCosts[2])
costsMatrix[,2]<-rtriang(nSims,a=parOMCosts[1],b=parOMCosts[3],c=parOMCosts[2])
costsMatrix[,3:(numMeasures+2)]<-costsMatrix[,1]*t(replicate(nSims,volumeDike))
colnames(costsMatrix)<-c('Unit Cost Cons','OM Per Cost',as.character(heightDike))

##Estimate net benefits - lifecycle 100 years (NPV) --------------------------------------------------------
rDis<-0.09
allEAD<-cbind(EAD,EADMeasures)
totalCostsAll<-cbind(rep(0,nSims),costsMatrix[,3:(numMeasures+2)]*(1+costsMatrix[,2]/rDis))
benefitsAll<-cbind(rep(0,nSims),EAD-EADMeasures)
totalBenefitsAll<-benefitsAll/rDis
NB<-totalBenefitsAll-totalCostsAll






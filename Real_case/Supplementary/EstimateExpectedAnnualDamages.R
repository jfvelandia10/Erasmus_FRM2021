#Load libraries

#Estimate risk of BAU and measures
##Risk analysis received hazard result, exposure parameters and vulnerability parameters uncertainty

#Function that calculates damage based on water level,exposed value and vulnerability curve
#House Low - use 1, material 5,7,8
#House Med - use 1, material 1
#House High - use 1, material 2
#COmmerce Med - use 2, material 1
#Commerce High - use other, material 1 or 2

estimateDamage<-function(WL,expValuesPar,vulValues){
  if (is.na(WL)){return(c(NA,NA))
    
  }else{
  #Calculate exposed value of simulation
  expValues<-matrix(0,nrow = nrow(tabExp),ncol=2)
  expValues[which(tabExp$Use ==1 & (tabExp$Material==5 | tabExp$Material==7 | tabExp$Material==8)),1]<-expValuesPar[1]
  expValues[which(tabExp$Use ==1 & tabExp$Material==1),1]<-expValuesPar[2]
  expValues[which(tabExp$Use ==1 & tabExp$Material==2),1]<-expValuesPar[3]
  expValues[which(tabExp$Use ==2 & tabExp$Material==1),1]<-expValuesPar[4]
  expValues[which(tabExp$Use>2),1]<-expValuesPar[5]
  expValues[,2]<-expValues[,1]*tabExp$TotalArea
  
  #Calculate water depths of event-simulation
  WDValues<-matrix(0,nrow=nrow(tabExp),ncol=1)
  WDValues<-WL-tabExp$ZWSTerrain
  WDValues[WDValues<0]<-0
  
  #Apply vulnerability functions
  damFactor<-matrix(0,nrow=nrow(tabExp),ncol=1)
  damValues<-matrix(0,nrow=nrow(tabExp),ncol=1)
  if (vulValues[1]>0.5) {
    #Brick 1 floor
    damFactor[which(tabExp$Material==1 & tabExp$NoFloors==1)]<-brick1floor(WDValues[which(tabExp$Material==1 & tabExp$NoFloors==1)])*(1+vulValues[2])
    #Brick 2-3 floors
    damFactor[which(tabExp$Material==1 & tabExp$NoFloors>1)]<-brick23floor(WDValues[which(tabExp$Material==1 & tabExp$NoFloors>1)])*(1+vulValues[2])
    #Concrete 1 floor
    damFactor[which(tabExp$Material==2 & tabExp$NoFloors==1)]<-concrete1floor(WDValues[which(tabExp$Material==2 & tabExp$NoFloors==1)])*(1+vulValues[2])
    #Concrete 2-3 floors
    damFactor[which(tabExp$Material==2 & tabExp$NoFloors>1)]<-concrete2floor(WDValues[which(tabExp$Material==2 & tabExp$NoFloors>1)])*(1+vulValues[2])
    #Mud
    damFactor[which(tabExp$Material>2)]<-mud1floor(WDValues[which(tabExp$Material>2)])*(1+vulValues[2])
  } else {
    damFactor[which(tabExp$Use==1)]<-resJRC(WDValues[which(tabExp$Use==1)])*(1+vulValues[3])
    damFactor[which(tabExp$Use!=1)]<-commJRC(WDValues[which(tabExp$Use!=1)])*(1+vulValues[3])
  }
  
  #Calculate damages
  damFactor[damFactor>1]<-1
  damFactor[damFactor<0]<-0
  damValues<-damFactor*expValues[,2]
  totalDamage<-sum(damValues)
  
  #Return total damage and total exposed value 
  return(c(totalDamage,sum(expValues[,2])))}
}


##Estimate uncertainty of direct-total losses damages
#Uncertainty of direct-total damages
parIndDam<-c(2,3)
indDamMatrix<-runif(nSims,min=parIndDam[1],max=parIndDam[2])

#Estimate damages for the entire simulations and events
allInputs<-cbind(wdepthWetlandBAU,exposureParMatrix,vulParMatrix)

#Iterate
resMatrix<-apply(allInputs,1,function(x) 
    apply(as.matrix(x[1:nEvents]),1,estimateDamage,expValues=x[(nEvents+1):(nEvents+5)],vulValues=tail(x, n=3)))

#Rearrange results
totalExposure<-t(resMatrix[seq(2,nEvents*2,2),])
dirDamMatrix<-t(resMatrix[seq(1,nEvents*2-1,2),])

##Estimate total damages by multliplying the direct-total factor
totalDamMatrix<-dirDamMatrix*indDamMatrix

##Estimate EAD per column, omit NAs and assume that number of years is valid years*50/44
EAD<- rowSums(totalDamMatrix,na.rm = T)/(rowSums(!is.na(totalDamMatrix))*50/44)

##Estimate convergence of EAD and SD
cumMeanEAD<-cumsum(EAD)/seq(1:nSims)
cumSDEAD<-sqrt(cumsum(EAD^2)/seq(1:nSims)-(cumsum(EAD)/seq(1:nSims))^2)

#####RISk for measures - dike at 10 different levels (0.5,1,1.5,2,2.5,3,3.5,4,4.5,5)----------------------

lapply(list.files('Results',full.names = T),load,.GlobalEnv)
numMeasures<-7
wdepthMeas<-replicate(7,wdepthWetlandBAU)
wdepthMeas[,,1]<-wdepthDike05
wdepthMeas[,,2]<-wdepthDike1
wdepthMeas[,,3]<-wdepthDike15
wdepthMeas[,,4]<-wdepthDike2
wdepthMeas[,,5]<-wdepthDike25
wdepthMeas[,,6]<-wdepthDike3
wdepthMeas[,,7]<-wdepthDike3

wdepthMeas[wdepthMeas<7 |wdepthMeas>20]<-NA

dirDamMeasMatrix<-array(,dim=c(nSims,nEvents,numMeasures))
totalDamMeasMatrix<-array(,dim=c(nSims,nEvents,numMeasures))
EADMeasures<-matrix(0,nrow = nSims,ncol=numMeasures)

for (i in 1:numMeasures){
  allInputsIter<-cbind(wdepthMeas[,,i],exposureParMatrix,vulParMatrix)
  resMatrixIterations<-apply(allInputsIter,1,function(x) 
    apply(as.matrix(x[1:nEvents]),1,estimateDamage,expValues=x[(nEvents+1):(nEvents+5)],vulValues=tail(x, n=3)))
  
  dirDamMeasMatrix[,,i]<-t(resMatrixIterations[seq(1,nEvents*2-1,2),])
  totalDamMeasMatrix[,,i]<-dirDamMeasMatrix[,,i]*indDamMatrix
  
  ##Estimate EAD per column, omit NAs and assume that number of years is valid years*50/44
  EADMeasures[,i]<- rowSums(totalDamMeasMatrix[,,i],na.rm = T)/(rowSums(!is.na(totalDamMeasMatrix[,,i]))*50/44)
}




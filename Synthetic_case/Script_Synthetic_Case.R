#Load libraries, scripts and fonts
library(distr)
library(extraDistr)
library(pracma)
library(RColorBrewer)
library(plot3D)
library(plotly)
library(ggplot2)
library(Hmisc)
library(fanplot)
library(zoo)
library(reshape2)
library(metR)
library(cowplot)
library(gridExtra)
source("FilledContourFunc.r")
source("FilledLegend.r")
source("fanJF.r")
library(BCEA)
windowsFonts(Times=windowsFont("Times New Roman"))

##Define inputs
#Uncertainties:
  #Flood level
  #Depth-damage curve
  #Cost of the house

#Measures:
  #BAU
  #Increase dyke
  #Raise houses

numIt=1000
set.seed(5)

###Define Hazard function for different return periods (10,25,100)-----
TR<-c(2,5,10,25,50,100,200)
 
#Function parameters of a GEV treated as normal variables
locationDist<-Norm(mean=562,sd=23)
scaleDist<-Norm(mean=160,sd=13)
shapeDist<-Norm(mean=0.02423,sd=0.07265)

#Simulation of parameters
parGEV<-cbind(
locMCr<-r(locationDist)(numIt),
scaMCr<-r(scaleDist)(numIt),
shaMCr<-r(shapeDist)(numIt))

qFreqCurve<-matrix(0, nrow=numIt, ncol=length(TR))
colnames(qFreqCurve)<-TR

#Estimation of discharge values for different return periods
for(i in 1:numIt){
  for(j in 1:length(TR)){
    qFreqCurve[i,j]<-qgev(1-1/TR[j],mu=parGEV[i,1],sigma=parGEV[i,2],xi=parGEV[i,3])
  }
}  

#River rating curve
hFreqCurve<-(qFreqCurve/56)^(1/1.7)

#Simple volume model assuming a peak duration of 1 hour and a pond
#The discharge entering to the urban area is estimated by a weir model
dikeLevelBAU<-5
lengthArea<-500
widthArea<-1000
urbanArea<-lengthArea*widthArea
timePeak<-60*60
lengthOverflow<-100

inflowBAU<-1.1*(hFreqCurve-dikeLevelBAU)^1.5
inflowBAU[is.nan(inflowBAU)] <- 0
floodLevelBAU<-inflowBAU*lengthOverflow*timePeak/(urbanArea)

#Plot results
#Discharge
dev.off()
par(family="Times")
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
matplot(TR,t(qFreqCurve),log="x",type='l',xlab='Return Period (years)',
        ylab=expression('Discharge (m'^3* '/s)'),col = brewer.pal(n = 9, name = "Greys"),
        main = "Discharge Frequency Curve", ylim = c(500,2500),yaxt='n',xaxt='n',
        panel.first=abline(v = TR,h =seq(500,2500,by=250), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(TR,rev(TR)),y=c(apply(qFreqCurve,2,min),rev(apply(qFreqCurve,2,max))),
                                               border = NA, col=rangecolor))
        )
axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.8)
axis(2, at = seq(500,2500,by=250), labels = prettyNum(seq(500,2500,by=250),big.mark = ' '),cex.axis=0.8)

#Stage
dev.off()
windowsFonts(Times=windowsFont("Times New Roman"))
par(family="Times")
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
rangecolor <- rgb(209,229,240,alpha=175,maxColorValue = 255)
matplot(TR,t(hFreqCurve),log="x",type='l',xlab='Return Period (years)',
        ylab='Stage (m)',col = brewer.pal(n = 9, name = "Greys"),
        main = "Stage Frequency Curve", ylim = c(3,9),yaxt='n',xaxt='n',
        panel.first=abline(v = TR,h =seq(3,9,by=1), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(TR,rev(TR)),y=c(apply(hFreqCurve,2,min),rev(apply(hFreqCurve,2,max))),
                                               border = NA, col=rangecolor))
)
axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.8)
axis(2, at = seq(3,9,by=1), labels = prettyNum(seq(3,9,by=1),big.mark = ' '),cex.axis=0.8)


##Discharge - stage - Simulation lines
dev.off()
png("Discharge - Stage Frequency Curves.png",width = 150, height = 130, units='mm', res=500)
par(family="Times",mar=c(0.35, 2.5, 2, 0.5),mgp=c(1.25,0.25,0),mfrow=c(2,1),tck=-0.02)
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
matplot(TR,t(qFreqCurve),log="x",type='l',xlab='Return Period (years)',
        ylab=expression('Discharge (m'^3* '/s)'),cex.lab=0.8,lwd = 0.25,
        ylim = c(500,2500),yaxt='n',xaxt='n',col = "lightgray",
        panel.first=abline(v = TR,h =seq(500,2500,by=250), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(TR,rev(TR)),y=c(apply(qFreqCurve,2,min),rev(apply(qFreqCurve,2,max))),
                                               border = NA, col=rangecolor))
)
title("Discharge", line = -1,cex.main = 0.9)
axis(2, at = seq(500,2500,by=500), labels = prettyNum(seq(500,2500,by=500),big.mark = ' '),cex.axis=0.7)
rug(seq(750,2250,by=500), ticksize = -0.02, side = 2)
rug(TR, ticksize = -0.02, side = 1)
mtext(expression(bold("Frequency Curves of the River")), side = 3,line = -1.75, outer = TRUE,cex = 1.1,adj =0.5)

par(mar=c(2.1,2.5,0.35, 0.5))
rangecolor <- rgb(255,238,153,alpha=175,maxColorValue = 255)
matplot(TR,t(hFreqCurve),log="x",type='l',xlab='Return Period (years)',
        ylab='Stage (m)',cex.lab=0.8,lwd = 0.25,
        ylim = c(3,11),yaxt='n',xaxt='n',
        panel.first=abline(v = TR,h =seq(3,11,by=1), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(TR,rev(TR)),y=c(apply(hFreqCurve,2,min),rev(apply(hFreqCurve,2,max))),
                                               border = NA, col=rangecolor))
)
title("Stage", line = -1,cex.main =0.9)
#axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.5)
mgp.axis.labels(c(0.25,0.2,0), type='x')
mgp.axis(1,at=TR,labels=prettyNum(TR,big.mark = ' '),cex.axis =0.7)
axis(2, at = seq(3,11,by=1), labels = prettyNum(seq(3,11,by=1),big.mark = ' '),cex.axis=0.7)
rug(TR, ticksize = -0.02, side = 3)
dev.off()

##Discharge - stage - Final as FAN plot
dev.off()
png("Discharge - Stage Frequency Curves - Thesis.png",width = 150, height = 130, units='mm', res=500)
par(family="Times",mar=c(0.35, 2.5, 2, 0.5),mgp=c(1.25,0.25,0),mfrow=c(2,1),tck=-0.02)
plot(NULL,xlim=c(min(TR),max(TR)),ylim=c(500,2250),log='x',xlab='Return Period (years)',
     ylab=expression('Discharge (m'^3* '/s)'),cex.lab=0.8,yaxt='n',xaxt='n')
fanJF(data = qFreqCurve, TR=TR, type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = c(5,25,75,95),med.ln=TRUE,medlab=bquote(~bar('X')),roffset=0.1,rcex=0.5, med.col='white',
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Reds"))),alpha=0.55)
abline(v = TR,h =seq(500,2250,by=250), lty = "dotted", col = "lightgray",lwd=0.8)
title("Discharge", line = -1,cex.main = 0.9)
axis(2, at = seq(500,2250,by=500), labels = prettyNum(seq(500,2250,by=500),big.mark = ' '),cex.axis=0.7)
rug(seq(750,2250,by=500), ticksize = -0.02, side = 2)
rug(TR, ticksize = -0.02, side = 1)
mtext(expression(bold("Frequency Curves of the River")), side = 3,line = -1.75, outer = TRUE,cex = 1.1,adj =0.5)

par(mar=c(2.1,2.5,0.35, 0.5))
plot(NULL,xlim=c(min(TR),max(TR)),ylim=c(3.5,8.5),log='x',xlab='Return Period (years)',
     ylab='Water Stage (m)',cex.lab=0.8,yaxt='n',xaxt='n')
fanJF(data = hFreqCurve, TR=TR, type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = c(5,25,75,95),med.ln=TRUE,medlab=expression(bar(' X')),roffset=0.1,rcex=0.5, med.col='white',
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Blues"))),alpha=0.55)
abline(v = TR,h =seq(3.5,8.5,by=0.5), lty = "dotted", col = "lightgray",lwd=0.8)
title("Water Stage", line = -1,cex.main =0.9)
mgp.axis.labels(c(0.25,0.2,0), type='x')
mgp.axis(1,at=TR,labels=prettyNum(TR,big.mark = ' '),cex.axis =0.7)
axis(2, at = seq(3.5,8.5,by=1), labels = prettyNum(seq(3.5,8.5,by=1),big.mark = ' '),cex.axis=0.7)
rug(seq(4,8,by=1), ticksize = -0.02, side = 2)
rug(TR, ticksize = -0.02, side = 3)
dev.off()



#Flood level BAU
dev.off()
png("Flood Depth Frequency Curves.png",width = 110, height = 100, units='mm', res=300)
par(family="Times")
par(mgp=c(1.3,0.4,0))
par(mar=c(2.5, 2.5,2, 0.5),tck=-0.02)
rangecolor <- rgb(33,102,172,alpha=100,maxColorValue = 255)
matplot(TR,t(floodLevelBAU),log="x",type='l',xlab='Return Period (years)',
        ylab='Water Depth (m)',cex.lab=0.7,cex.main=0.9,
        main = "Flood Depth Frequency Curve - BAU", ylim = c(0,7),yaxt='n',xaxt='n',lwd=0.25,
        panel.first=abline(v = TR,h =seq(0,7,by=0.5), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(TR,rev(TR)),y=c(apply(floodLevelBAU,2,min),rev(apply(floodLevelBAU,2,max))),
                                               border = NA, col=rangecolor))
)
#axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.6)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1,at=TR,labels=prettyNum(TR,big.mark = ' '),cex.axis =0.7)
axis(2, at = seq(0,7,by=1), labels = prettyNum(seq(0,7,by=1),big.mark = ' '),cex.axis=0.6)
rug(seq(0.5,6.5,by=1), ticksize = -0.015, side = 2)
dev.off()

###Define vulnerability function ----------------------------------
#Function is a lognormal distribution following the GLobal Flood Depth-Damage
#function from the JRC
#Parameters mean and sd follow an uniform distribution - Bayesian approach
meanDistV<-Unif(-0.69-0.4*0.69,-0.69+0.4*0.69)
sdDistV<-Unif(0.95-0.4*0.95,0.95+0.4*0.95)

#Iterate to obtain a set of vulnerability functions
vulValues<-c(0,0.5,1,1.5,2,3,4,5,6)
vulResults<-matrix(vulValues, nrow=length(vulValues), ncol=numIt)


#Simulation of parameters
parLN<-cbind(
  meanMCr<-r(meanDistV)(numIt),
  sdMCr<-r(sdDistV)(numIt))

  
for(i in 1:numIt){
  L <- Lnorm(meanlog=parLN[i,1],sdlog=parLN[i,2])
  for(j in 1:length(vulValues)){
    vulResults[j,i]=p(L)(vulValues[j])
  }
}  

vulResults<-cbind(vulValues,vulResults)

#Plot vulnerability function simulations - lines
dev.off()
png("Vulnerability Curves - BAU.png",width = 110, height = 100, units='mm', res=300)
par(family="Times")
par(mgp=c(1.3,0.4,0))
par(mar=c(2.5, 2.5,2, 0.5),tck=-0.02)
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
matplot(vulResults[,1],vulResults[,-1],type='l',lwd = 0.01,xlab='Water Depth (m)',
        ylab='Damage Factor (-)',cex.lab=0.7,cex.main=0.9,
        main = "Vulnerability Function - BAU", ylim = c(0,1),yaxt='n',xaxt='n',
        panel.first=abline(v = vulResults[,1],h =seq(0,1,by=0.1), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(vulResults[,1],rev(vulResults[,1])),y=c(apply(vulResults[,-1],1,min),rev(apply(vulResults[,-1],1,max))),
                                               border = NA, col=rangecolor))
)
#axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.6)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = vulResults[,1], labels = prettyNum(vulResults[,1],big.mark = ' '),cex.axis =0.6)
axis(2, at = seq(0,1,by=0.1), labels = prettyNum(seq(0,1,by=0.1),big.mark = ' '),cex.axis=0.6)
dev.off()

#Plot vulnerability function simulations - FINAL FAN PLOT
dev.off()
png("Vulnerability Curves.png",width = 110, height = 100, units='mm', res=300)
par(family="Times",mgp=c(1.3,0.4,0),mar=c(2.5, 2.5,2, 0.5),tck=-0.02)
plot(NULL,xlim=c(0,6),ylim=c(0,1),xlab='Water Depth (m)',
     ylab='Damage Factor (-)',cex.lab=0.7,yaxt='n',xaxt='n',
     main = "Vulnerability Function",cex.main=0.9)
fanJF(data = t(vulResults[,-1]), TR=vulResults[,1], type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = c(5,95),med.ln=TRUE, med.col='white',rcex=0.0001,rcol='white',
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Reds"))),alpha=0.55)
text(x=1.5,y=0.76,labels='5%',cex= 0.6, pos=1,offset=-0.1,srt=40,family="Times")
text(x=1.5,y=0.9,labels=bquote(bar('X')),cex= 0.6, pos=3,offset=-0.1,srt=35,family="Times",col='white')
text(x=1.5,y=0.99,labels='95%',cex= 0.6, pos=3,offset=-0.1,srt=25,family="Times")
abline(v=vulResults[,1],h =seq(0,1,by=0.1), lty = "dotted", col = "lightgray",lwd=0.8)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = vulResults[,1], labels = prettyNum(vulResults[,1],big.mark = ' '),cex.axis =0.6)
axis(2, at = seq(0,1,by=0.1), labels = prettyNum(seq(0,1,by=0.1),big.mark = ' '),cex.axis=0.6)
dev.off()


###Define number of houses and prices------------------------------------
numHouses<-400
pricepersqDist<-rtriang(numIt,a=250,b=1700,c=950)
areaHouse<-150 

#Simulation of prices
priceMC<-pricepersqDist*numHouses*areaHouse

#Plot histogram of prices
dev.off()
png("Exposed Value.png",width = 110, height = 100, units='mm', res=300)
par(mgp=c(1.3,0.4,0),mar=c(2.5, 2.5,2, 0.5),tck=-0.02,family="Times")
rangecolor <- rgb(235,219,199,alpha=255,maxColorValue = 255)
hist(priceMC/10^6, xlab = paste("Value (","\u20AC"," Million)"),
     main = 'Total Exposed Value',col =rangecolor,ylim = c(0,150),xlim=c(10,110),breaks = 20,
     xaxt='n',yaxt='n',cex.lab=0.7,cex.main=0.9)
abline(v = seq(0,110,by=10) ,h =seq(0,150,by=12.5), lty = "dotted", col = "lightgray")
axis(1,at = seq(10,110,by=10), labels = prettyNum(seq(10,110,by=10),big.mark = ' '),cex.axis=0.6)
axis(2,at = seq(0,150,by=25), labels = prettyNum(seq(0,150,by=25),big.mark = ' '),cex.axis=0.6)
rug(seq(12.5,137.5,by=25), ticksize = -0.015, side = 2)
#grid()
box()
hist(priceMC/10^6, xlab = paste("Value (","\u20AC"," Million)"),
     main = 'Exposed Value',col =rangecolor,ylim = c(0,150),xlim=c(10,110),breaks = 20,
     xaxt='n',yaxt='n',cex.lab=0.7,cex.main=0.9,add = TRUE)

dev.off()

#Plot cummulative histogram - ECDF



###Estimate Risk level - LEC curve and EAD------------------------------

ELCResultsBAU<-matrix(TR, nrow=length(TR), ncol=numIt)
vulInondVal<-matrix(TR, nrow=length(TR), ncol=numIt)

for(i in 1:numIt){
  for (j in 1:length(TR)){
    
    inondVal<-floodLevelBAU[i,j]
    L <- Lnorm(meanlog=parLN[i,1],sdlog=parLN[i,2])
    vulInondVal[j,i]<-p(L)(inondVal)
    ELCResultsBAU[j,i]=vulInondVal[j,i]*priceMC[i]
    
  }
}

ELCResultsBAU<-cbind(TR,ELCResultsBAU)

  
#Plot results LEC and Risk Curve
#Risk curve -- NOT PNG
dev.off()
windowsFonts(Times=windowsFont("Times New Roman"))
par(family="Times")
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
rangecolor <- rgb(214,96,77,alpha=200,maxColorValue = 255)
matplot(ELCResultsBAU[,1],ELCResultsBAU[,-1]/10^6,log="x",type='l',lwd = 0.1,xlab='Return Period (years)',
        ylab=paste("Losses (","\u20AC"," Million)"),col = brewer.pal(n = 9, name = "Greys"),
        main = "Risk Curve - BAU", ylim = c(0,100),yaxt='n',xaxt='n',
        panel.first=abline(v = ELCResultsBAU[,1],h =seq(0,100,by=10), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(ELCResultsBAU[,1],rev(ELCResultsBAU[,1])),y=c(apply(ELCResultsBAU[,-1]/10^6,1,min),rev(apply(ELCResultsBAU[,-1]/10^6,1,max))),
                                               border = NA, col=rangecolor))
)
axis(1, at = ELCResultsBAU[,1], labels = prettyNum(ELCResultsBAU[,1],big.mark = ' '),cex.axis=0.8)
axis(2, at = seq(0,100,by=10), labels = prettyNum(seq(0,100,by=10),big.mark = ' '),cex.axis=0.8)

#Loss exceedance and Risk Curve curve
dev.off()
png("Loss exceedance curves - BAU.png",width = 110, height = 120, units='mm', res=300)
par(family="Times")
par(mgp=c(1.3,0.4,0))
par(mar=c(6, 2.5,2, 0.5),tck=-0.02)
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
matplot(1/ELCResultsBAU[,1],ELCResultsBAU[,-1]/10^6,log="x",type='l',xlab='Exceedance rate (1/year)',
        ylab=paste("Losses (","\u20AC"," Million)"),cex.lab=0.7,cex.main=0.9,
        main = "Loss Exceedance Curve - BAU", ylim = c(0,100),yaxt='n',xaxt='n',lwd=0.25,
        panel.first=abline(v = 1/ELCResultsBAU[,1],h =seq(0,100,by=10),lty = "dotted",  col = "lightgray",
                           anel.first=polygon(x=c(1/ELCResultsBAU[,1],rev(1/ELCResultsBAU[,1])),y=c(apply(ELCResultsBAU[,-1]/10^6,1,min),rev(apply(ELCResultsBAU[,-1]/10^6,1,max))),
                                              border = NA, col=rangecolor))
)
#axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.6)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = 1/ELCResultsBAU[,1], labels = prettyNum(1/ELCResultsBAU[,1],big.mark = ' '),cex.axis =0.6)
axis(2, at = seq(0,100,by=10), labels = prettyNum(seq(0,100,by=10),big.mark = ' '),cex.axis=0.6)

axis(1,line=3,at=1/ELCResultsBAU[,1],labels=TR,tick=TRUE,mgp=c(1.3,0.4,0),cex.axis=0.6)
mtext("Return Period (years)",1,line=4.3,adj=0.5,cex=0.7)
dev.off()


#Statistics
summary(t(ELCResultsBAU[,-1])/10^6)

#Estimate EAD
EADResultsBAU<-rep(0, numIt)
for(i in 1:numIt){
  EADResultsBAU[i]<-trapz(rev(1/ELCResultsBAU[,1]),rev(ELCResultsBAU[,i+1]))/10^6
}

#Plot EAD
dev.off()
png("EAD - BAU.png",width = 55, height = 100, units='mm', res=300)
par(mgp=c(1.3,0.4,0),mar=c(0.5, 2.5,2, 1.4),tck=-0.02,family="Times")
rangecolor <- rgb(244,165,130,alpha=255,maxColorValue = 255)
boxplot(EADResultsBAU,main = "Expected Annual Damage - BAU" ,
        ylab=paste("Losses/year (","\u20AC"," Million)"),border = NA, 
        xaxt='n', yaxt = "n", frame = FALSE,col=NA,ylim=c(0,9.2),cex.lab=0.7,cex.main=0.8)
abline(h =seq(0,9,by=0.5), lty = "dotted", col = "lightgray")
axis(2, at = seq(0,9,by=1), labels = prettyNum(seq(0,9,by=1),big.mark = ' '),cex.axis=0.6)
rug(seq(0.5,9,by=1), ticksize = -0.015, side = 2)
boxplot(EADResultsBAU,yaxt = "n",ann = FALSE, outcex=0.6,
        border=rangecolor,boxlwd = 1,col=NA,add = TRUE,ylim=c(0,9.2))
text(mean(EADResultsBAU),labels=bquote(bar('EAD')~'='~.(round(mean(EADResultsBAU),digits = 1))),cex= 0.6, pos=3,offset=-0.1,family="Times")
dev.off()

summary(EADResultsBAU)

####Measures risk calculation---------------------------------------------------------
##Measure 1 - Dike 6 meters (25-50 years) 
#Hazard 
#Simple volume model assuming a peak duration of 1 hour and a pond
#The discharge entering to the urban area is estimated by a weir model
dikeLevelMea1<-6
inflowMea1<-1.1*(hFreqCurve-dikeLevelMea1)^1.5
inflowMea1[is.nan(inflowMea1)] <- 0
floodLevelMea1<-inflowMea1*lengthOverflow*timePeak/(urbanArea)

#Flood level 
dev.off()
windowsFonts(Times=windowsFont("Times New Roman"))
par(family="Times")
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
rangecolor <- rgb(33,102,172,alpha=100,maxColorValue = 255)
matplot(TR,t(floodLevelMea1),log="x",type='l',xlab='Return Period (years)',
        ylab='Water Depth (m)',col = brewer.pal(n = 9, name = "Greys"),
        main = "Flood Depth Frequency Curve - Measure 1", ylim = c(0,5),yaxt='n',xaxt='n',
        panel.first=abline(v = TR,h =seq(0,5,by=0.5), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(TR,rev(TR)),y=c(apply(floodLevelMea1,2,min),rev(apply(floodLevelMea1,2,max))),
                                               border = NA, col=rangecolor))
)
axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.8)
axis(2, at = seq(0,5,by=0.5), labels = prettyNum(seq(0,5,by=0.5),big.mark = ' '),cex.axis=0.8)


#Same vulnerability function as the original

#Estimate risk when implementing measure 1
ELCResultsMea1<-matrix(TR, nrow=length(TR), ncol=numIt)
vulInondValMea1<-matrix(TR, nrow=length(TR), ncol=numIt)

for(i in 1:numIt){
  for (j in 1:length(TR)){
    
    inondVal<-floodLevelMea1[i,j]
    L <- Lnorm(meanlog=parLN[i,1],sdlog=parLN[i,2])
    vulInondValMea1[j,i]<-p(L)(inondVal)
    ELCResultsMea1[j,i]=vulInondValMea1[j,i]*priceMC[i]
    
  }
}

ELCResultsMea1<-cbind(TR,ELCResultsMea1)

#Estimate EAD - Measure 1
EADResultsMea1<-rep(0, numIt)
for(i in 1:numIt){
  EADResultsMea1[i]<-trapz(rev(1/ELCResultsMea1[,1]),rev(ELCResultsMea1[,i+1]))/10^6
}

#Plot EAD
dev.off()
windowsFonts(Times=windowsFont("Times New Roman"))
par(family="Times")
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
rangecolor <- rgb(244,165,130,alpha=255,maxColorValue = 255)
boxplot(EADResultsMea1,main = "Expected Annual Damage - Measure 1" ,
        ylab=paste("Losses/year (","\u20AC"," Million)"),border = NA, 
        xaxt='n', yaxt = "n", frame = FALSE,col=NA)
abline(h =seq(0,8,by=1), lty = "dotted", col = "lightgray")
axis(2, at = seq(0,8,by=1), labels = prettyNum(seq(0,8,by=1),big.mark = ' '),cex.axis=0.8)
boxplot(EADResultsMea1,main = "Expected Annual Damage - Measure 1" ,yaxt = "n",
        ylab=paste("Losses per year (","\u20AC"," Million)"),ann = FALSE,
        border=rangecolor,boxlwd = 1,col=NA,add = TRUE)
text(mean(EADResultsMea1),labels=round(mean(EADResultsMea1),digits = 2), cex= 0.8, pos=3,offset=-0.1,family="Times")

summary(EADResultsMea1)

##Measure 2 - Raise houses 50cm
#Same hazard as the original one

#Vulnerability function change - same as the original but displaced to the right to start at 0.5m
houseRaisMea2<-0.5
vulBAU<-vulResults

vulResultsMea2<-cbind(vulBAU[,1]+houseRaisMea2,vulBAU[,-1])
vulResultsMea2<-rbind(rep(0,numIt+1),vulResultsMea2)

#Plot vulnerability function simulations
matplot(vulResultsMea2[,1],vulResultsMea2[,-1],type='l')
dev.off()
windowsFonts(Times=windowsFont("Times New Roman"))
par(family="Times")
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
matplot(vulResultsMea2[,1],vulResultsMea2[,-1],type='l',lwd = 0.1,xlab='Water Depth (m)',
        ylab='Damage Factor (-)',col = brewer.pal(n = 9, name = "Greys"),
        main = "Vulnerability Function - BAU", ylim = c(0,1),yaxt='n',xaxt='n',
        panel.first=abline(v = vulResultsMea2[,1],h =seq(0,1,by=0.1), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(vulResultsMea2[,1],rev(vulResultsMea2[,1])),y=c(apply(vulResultsMea2[,-1],1,min),rev(apply(vulResultsMea2[,-1],1,max))),
                                               border = NA, col=rangecolor))
)
axis(1, at = vulResultsMea2[,1], labels = prettyNum(vulResultsMea2[,1],big.mark = ' '),cex.axis=0.8)
axis(2, at = seq(0,1,by=0.1), labels = prettyNum(seq(0,1,by=0.1),big.mark = ' '),cex.axis=0.8)

#Estimate risk when implementing measure 2
ELCResultsMea2<-matrix(TR, nrow=length(TR), ncol=numIt)
vulInondValMea2<-matrix(TR, nrow=length(TR), ncol=numIt)

for(i in 1:numIt){
  for (j in 1:length(TR)){
    
    inondVal<-floodLevelBAU[i,j]
    
    if (inondVal<0.5){
      inondVal=0
    }
    else{
      inondVal-0.5
    }
    
    L <- Lnorm(meanlog=parLN[i,1],sdlog=parLN[i,2])
    vulInondValMea2[j,i]<-p(L)(inondVal)
    ELCResultsMea2[j,i]=vulInondValMea2[j,i]*priceMC[i]
    
  }
}

ELCResultsMea2<-cbind(TR,ELCResultsMea2)

#Estimate EAD - Measure 1
EADResultsMea2<-rep(0, numIt)
for(i in 1:numIt){
  EADResultsMea2[i]<-trapz(rev(1/ELCResultsMea2[,1]),rev(ELCResultsMea2[,i+1]))/10^6
}

#Plot EAD
dev.off()
windowsFonts(Times=windowsFont("Times New Roman"))
par(family="Times")
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
rangecolor <- rgb(244,165,130,alpha=255,maxColorValue = 255)
boxplot(EADResultsMea2,main = "Expected Annual Damage - Measure 2" ,
        ylab=paste("Losses/year (","\u20AC"," Million)"),border = NA, 
        xaxt='n', yaxt = "n", frame = FALSE,col=NA)
abline(h =seq(0,8,by=1), lty = "dotted", col = "lightgray")
axis(2, at = seq(0,8,by=1), labels = prettyNum(seq(0,8,by=1),big.mark = ' '),cex.axis=0.8)
boxplot(EADResultsMea2,main = "Expected Annual Damage - Measure 2" ,yaxt = "n",
        ylab=paste("Losses per year (","\u20AC"," Million)"),ann = FALSE,
        border=rangecolor,boxlwd = 1,col=NA,add = TRUE)
text(mean(EADResultsMea2),labels=round(mean(EADResultsMea2),digits = 2), cex= 0.8, pos=3,offset=-0.1,family="Times")

summary(EADResultsMea2)

###Plot results risk for BAU and Measures-------------------------------------------
##Hazard
##Original all lines-
dev.off()
png("Flood Depth Frequency Curves - All measures.png",width = 160, height = 65, units='mm', res=300)
par(family="Times")
par(mgp=c(1.5,0.4,0))
par(mar=c(2.5, 2.5,1.25, 0.2),tck=-0.02)
par(mfcol=c(1,3))
#BAU
rangecolor <- rgb(33,102,172,alpha=100,maxColorValue = 255)
matplot(TR,t(floodLevelBAU),log="x",type='l',xlab='',
        ylab=expression(bold('Water Depth (m)')),cex.lab=0.9,cex.main=1.1,
        main = "BAU", ylim = c(0,7),yaxt='n',xaxt='n',lwd=0.5,
        panel.first=abline(v = TR,h =seq(0,7,by=0.5), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(TR,rev(TR)),y=c(apply(floodLevelBAU,2,min),rev(apply(floodLevelBAU,2,max))),
                                               border = NA, col=rangecolor))
)
#axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.6)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1,at=TR,labels=prettyNum(TR,big.mark = ' '),cex.axis =0.9)
axis(2, at = seq(0,7,by=1), labels = prettyNum(seq(0,7,by=1),big.mark = ' '),cex.axis=0.9)
rug(seq(0.5,6.5,by=1), ticksize = -0.015, side = 2)

#measure 1
par(mar=c(2.5, 0.2,1.25, 0.2),tck=-0.02)
rangecolor <- rgb(33,102,172,alpha=100,maxColorValue = 255)
matplot(TR,t(floodLevelMea1),log="x",type='l',xlab=expression(bold('Return Period (years)')),
        ylab='NA',cex.lab=0.9,cex.main=1.1,
        main = "Measure 1", ylim = c(0,7),yaxt='n',xaxt='n',lwd=0.5,
        panel.first=abline(v = TR,h =seq(0,7,by=0.5), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(TR,rev(TR)),y=c(apply(floodLevelMea1,2,min),rev(apply(floodLevelMea1,2,max))),
                                               border = NA, col=rangecolor))
)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1,at=TR,labels=prettyNum(TR,big.mark = ' '),cex.axis =0.9)
#axis(2, at = seq(0,7,by=1), labels = prettyNum(seq(0,7,by=1),big.mark = ' '),cex.axis=0.6)
rug(seq(0.5,6.5,by=1), ticksize = -0.015, side = 2)

#Measure 2
par(mar=c(2.5, 0.2,1.25, 0.5),tck=-0.02)
rangecolor <- rgb(33,102,172,alpha=100,maxColorValue = 255)
matplot(TR,t(floodLevelBAU),log="x",type='l',xlab='',
        ylab='NA',cex.lab=0.9,cex.main=1.1,
        main = "Measure 2", ylim = c(0,7),yaxt='n',xaxt='n',lwd=0.5,
        panel.first=abline(v = TR,h =seq(0,7,by=0.5), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(TR,rev(TR)),y=c(apply(floodLevelBAU,2,min),rev(apply(floodLevelBAU,2,max))),
                                               border = NA, col=rangecolor))
)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1,at=TR,labels=prettyNum(TR,big.mark = ' '),cex.axis =0.9)
#axis(2, at = seq(0,7,by=1), labels = prettyNum(seq(0,7,by=1),big.mark = ' '),cex.axis=0.6)
rug(seq(0.5,6.5,by=1), ticksize = -0.015, side = 2)
dev.off()

##FINAL FAN plot-
dev.off()
png("Flood Depth Frequency Curves - All measures-Thesis.png",width = 160, height = 65, units='mm', res=300)
par(family="Times",mgp=c(1.3,0.4,0),mar=c(2.5, 2.5,1.25, 0.2),tck=-0.02,mfcol=c(1,3),las=1)
#BAU
plot(NULL,xlim=c(min(TR),max(TR)+50),ylim=c(0,5),log='x',xlab='',
     ylab='Water Depth (m)',cex.lab=0.9,yaxt='n',xaxt='n',cex.main=1.1,main="BAU")
fanJF(data = floodLevelBAU, TR=TR, type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = c(5,25,75,95),med.ln=TRUE,medlab=bquote(~bar('X')),roffset=0.15,rcex=0.7, med.col='white',rpos=4,
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Blues"))),alpha=0.55)
abline(v = TR,h =seq(0,7,by=0.5), lty = "dotted", col = "lightgray",lwd=0.8)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1,at=TR,labels=prettyNum(TR,big.mark = ' '),cex.axis =0.9)
axis(2, at = seq(0,7,by=1), labels = prettyNum(seq(0,7,by=1),big.mark = ' '),cex.axis=0.9)
rug(seq(0.5,6.5,by=1), ticksize = -0.015, side = 2)

#Measure 1
par(mar=c(2.5, 0.2,1.25, 0.2),tck=-0.02)
plot(NULL,xlim=c(min(TR),max(TR)+50),ylim=c(0,5),log='x',xlab=expression(bold('Return Period (years)')),
     ylab='NA',cex.lab=0.9,yaxt='n',xaxt='n',cex.main=1.1,main="Raise dike")
fanJF(data = floodLevelMea1, TR=TR, type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = c(5,25,75,95),med.ln=TRUE,medlab=bquote(~bar('X')),roffset=0.15,rcex=0.7, med.col='white',rpos=4,
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Blues"))),alpha=0.55)
abline(v = TR,h =seq(0,7,by=0.5), lty = "dotted", col = "lightgray",lwd=0.8)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1,at=TR,labels=prettyNum(TR,big.mark = ' '),cex.axis =0.9)
rug(seq(0,6,by=0.5), ticksize = -0.015, side = 2)

#Measure 2
par(mar=c(2.5, 0.2,1.25,1),tck=-0.02)
plot(NULL,xlim=c(min(TR),max(TR)+50),ylim=c(0,5),log='x',xlab='',
     ylab='NA',cex.lab=0.9,yaxt='n',xaxt='n',cex.main=1.1,main="Elevate houses")
fanJF(data = floodLevelBAU, TR=TR, type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = c(5,25,75,95),med.ln=TRUE,medlab=bquote(~bar('X')),roffset=0.15,rcex=0.7, med.col='white',rpos=4,
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Blues"))),alpha=0.55)
abline(v = TR,h =seq(0,7,by=0.5), lty = "dotted", col = "lightgray",lwd=0.8)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1,at=TR,labels=prettyNum(TR,big.mark = ' '),cex.axis =0.9)
axis(4, at = seq(0,7,by=1), labels = prettyNum(seq(0,7,by=1),big.mark = ' '),cex.axis=0.9)
rug(seq(0,6,by=0.5), ticksize = -0.015, side = 2)

dev.off()

##Vulnerability and Exposure
dev.off()
png("Vulnerability-Exposure.png",width = 160, height = 80, units='mm', res=300)
par(family="Times")
par(mfcol=c(1,2))
par(mgp=c(1.5,0.4,0))
#Exposure
par(mar=c(2.5, 2.5,1.5, 0.2),tck=-0.02)
rangecolor <- rgb(235,219,199,alpha=255,maxColorValue = 255)
hist(priceMC/10^6, xlab = expression(bold(paste("Value (","\u20AC"," Million)"))),
     main = 'Total Maximum Damage Value',ylim = c(0,0.03),xlim=c(10,110),breaks = 20,
     xaxt='n',yaxt='n',cex.lab=0.7,cex.main=0.8,prob=T,ylab = expression(bold(Density)))
abline(v = seq(0,110,by=10) ,h =seq(0,0.03,by=0.0025), lty = "dotted", col = "lightgray")
axis(1,at = seq(10,110,by=10), labels = prettyNum(seq(10,110,by=10),big.mark = ' '),cex.axis=0.6)
axis(2,at = seq(0,0.03,by=0.005), labels = prettyNum(seq(0,0.03,by=0.005),big.mark = ' '),cex.axis=0.6)
rug(seq(0.0025,0.03-0.005,by=0.005), ticksize = -0.015, side = 2)
#grid()
box()
hist(priceMC/10^6, xlab = paste("Value (","\u20AC"," Million)"),
     main = 'Exposed Value',ylim = c(0,0.03),xlim=c(10,110),breaks = 20,
     xaxt='n',yaxt='n',cex.lab=0.7,cex.main=0.9,add = TRUE, prob=T,col = 'white')
lines(density(priceMC/10^6), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")
abline(v=mean(priceMC/10^6),col="black",lwd=1.5,lty=5)
abline(v = median(priceMC/10^6),col = "red",lwd=1.5,lty=5)
legend(x = "topright", # location of legend within plot area
       c("Density plot", "Mean", "Median"),
       col = c("chocolate3", "black", "red"),
       lwd = c(1, 1, 1), cex = 0.6)

#Vulnerability - fan plot
par(mar=c(4.5, 2,1.5, 1),tck=-0.02,mgp=c(1.1,0.4,0))
plot(NULL,xlim=c(min(vulResults[,1]),max(vulResults[,1])+0.2),ylim=c(0,100),xlab=expression(bold('Water depth (m) - BAU and Raise dike')),
     ylab=expression(bold(paste("Damage factor (%)"))),cex.lab=0.6,yaxt='n',xaxt='n',cex.main=0.8,main="Depth-Damage Function")
fanJF(data = t(vulResults[,-1]*100), TR=vulResults[,1], type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = 50,med.ln=TRUE,medlab=bquote(~bar('X')),roffset=-10,rcex=0.7, med.col='white',rpos=4,rcol = 'white',
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Reds"))),alpha=0.55)
abline(v = vulResults[,1],h =seq(0,100,by=10), lty = "dotted", col = "lightgray",lwd=0.8)
mtext(expression(bar('X')),1,line=-9.8,adj=0.98,cex=0.5)
mgp.axis.labels(c(0.2,0.2,0), type='x')
mgp.axis(1, at = vulResults[,1], labels = prettyNum((vulResults[,1]),big.mark = ' '),cex.axis =0.6)
axis(1,line=2.3,at=vulResults[,1],labels=prettyNum((vulResults[,1]+0.5),big.mark = ' '),tick=TRUE,mgp=c(1.2,0.4,0),cex.axis=0.6)
axis(2, at = seq(0,100,by=10), labels = prettyNum(seq(0,100,by=10),big.mark = ' '),cex.axis=0.6)
rug(seq(0,100,by=10), ticksize = -0.015, side = 2)
mtext(expression(bold("Water depth (m) - Elevate houses")),1,line=3.5,adj=0.5,cex=0.6)
dev.off()


##Risk - LEC -----
#Lines
dev.off()
png("LEC - All measures.png",width = 160, height = 80, units='mm', res=300)
par(family="Times")
par(mfcol=c(1,3))
par(mgp=c(1.5,0.4,0))
#BAU
par(mar=c(6, 2.5,1.5, 0.2),tck=-0.02)
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
matplot(1/ELCResultsBAU[,1],ELCResultsBAU[,-1]/10^6,log="x",type='l',xlab='',
        ylab=expression(bold(paste("Losses (","\u20AC"," Million)"))),cex.lab=0.9,cex.main=1.1,
        main = "BAU", ylim = c(0,100),yaxt='n',xaxt='n',lwd=0.5,
        panel.first=abline(v = 1/ELCResultsBAU[,1],h =seq(0,100,by=10),lty = "dotted",  col = "lightgray",
                           anel.first=polygon(x=c(1/ELCResultsBAU[,1],rev(1/ELCResultsBAU[,1])),y=c(apply(ELCResultsBAU[,-1]/10^6,1,min),rev(apply(ELCResultsBAU[,-1]/10^6,1,max))),
                                              border = NA, col=rangecolor))
)
#axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.6)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = 1/ELCResultsBAU[,1], labels = prettyNum(1/ELCResultsBAU[,1],big.mark = ' '),cex.axis =0.9)
axis(2, at = seq(0,100,by=10), labels = prettyNum(seq(0,100,by=10),big.mark = ' '),cex.axis=0.9)
axis(1,line=3,at=1/ELCResultsBAU[,1],labels=TR,tick=TRUE,mgp=c(1.3,0.4,0),cex.axis=0.9)

#Measure 1
par(mar=c(6, 0.2,1.5, 0.2),tck=-0.02)
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
matplot(1/ELCResultsMea1[,1],ELCResultsMea1[,-1]/10^6,log="x",type='l',xlab=expression(bold('Exceedance rate (1/year)')),
        ylab='',cex.lab=0.9,cex.main=1.1,
        main = "Measure 1", ylim = c(0,100),yaxt='n',xaxt='n',lwd=0.5,
        panel.first=abline(v = 1/ELCResultsBAU[,1],h =seq(0,100,by=10),lty = "dotted",  col = "lightgray",
                           anel.first=polygon(x=c(1/ELCResultsBAU[,1],rev(1/ELCResultsBAU[,1])),y=c(apply(ELCResultsMea1[,-1]/10^6,1,min),rev(apply(ELCResultsMea1[,-1]/10^6,1,max))),
                                              border = NA, col=rangecolor))
)
#axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.6)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = 1/ELCResultsBAU[,1], labels = prettyNum(1/ELCResultsBAU[,1],big.mark = ' '),cex.axis =0.9)
axis(1,line=3,at=1/ELCResultsBAU[,1],labels=TR,tick=TRUE,mgp=c(1.3,0.4,0),cex.axis=0.9)
mtext(expression(bold("Return Period (years)")),1,line=4.5,adj=0.5,cex=0.6)

#Measure 2
par(mar=c(6, 0.2,1.5, 0.5),tck=-0.02)
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
matplot(1/ELCResultsMea2[,1],ELCResultsMea2[,-1]/10^6,log="x",type='l',xlab='',
        ylab='',cex.lab=0.9,cex.main=1.1,
        main = "Measure 2", ylim = c(0,100),yaxt='n',xaxt='n',lwd=0.5,
        panel.first=abline(v = 1/ELCResultsBAU[,1],h =seq(0,100,by=10),lty = "dotted",  col = "lightgray",
                           anel.first=polygon(x=c(1/ELCResultsMea2[,1],rev(1/ELCResultsMea2[,1])),y=c(apply(ELCResultsMea2[,-1]/10^6,1,min),rev(apply(ELCResultsMea2[,-1]/10^6,1,max))),
                                              border = NA, col=rangecolor))
)
#axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.6)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = 1/ELCResultsBAU[,1], labels = prettyNum(1/ELCResultsBAU[,1],big.mark = ' '),cex.axis =0.9)
axis(1,line=3,at=1/ELCResultsBAU[,1],labels=TR,tick=TRUE,mgp=c(1.3,0.4,0),cex.axis=0.9)
dev.off()

#FINAL FAN PLOT --- ##FINAL FAN plot-
dev.off()
png("LEC - All measures - Final.png",width = 160, height = 80, units='mm', res=300)
par(family="Times",mfcol=c(1,3),mgp=c(1.5,0.4,0))
#BAU
par(mar=c(6, 2.5,1.5, 0.2),tck=-0.02)
plot(NULL,xlim=c(min(1/ELCResultsBAU[,1])-0.001,max(1/ELCResultsBAU[,1])),ylim=c(0,90),log='x',xlab='',
     ylab=expression(bold(paste("Damages (","\u20AC"," Million)"))),cex.lab=0.9,yaxt='n',xaxt='n',cex.main=1.1,main="BAU")
fanJF(data = t(ELCResultsBAU[,-1]/10^6), TR=1/ELCResultsBAU[,1], type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = c(5,25,75,95),med.ln=TRUE,medlab=bquote(~bar('X')),roffset=0.2,rcex=0.7, med.col='white',rpos=2,
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Greys"))),alpha=0.55)
abline(v = 1/ELCResultsBAU[,1],h =seq(0,100,by=10), lty = "dotted", col = "lightgray",lwd=0.8)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = 1/ELCResultsBAU[,1], labels = prettyNum(1/ELCResultsBAU[,1],big.mark = ' '),cex.axis =0.7)
axis(2, at = seq(0,100,by=10), labels = prettyNum(seq(0,100,by=10),big.mark = ' '),cex.axis=0.7)
axis(1,line=3,at=1/ELCResultsBAU[,1],labels=TR,tick=TRUE,mgp=c(1.3,0.4,0),cex.axis=0.7)

#Measure 1
par(mar=c(6, 0.2,1.5, 0.2),tck=-0.02)
plot(NULL,xlim=c(min(1/ELCResultsMea1[,1])-0.001,max(1/ELCResultsMea1[,1])),ylim=c(0,90),log='x',xlab=expression(bold('Exceedance probability')),
     ylab='',cex.lab=0.9,yaxt='n',xaxt='n',cex.main=1.1,main="Raise dike")
fanJF(data = t(ELCResultsMea1[,-1]/10^6), TR=1/ELCResultsMea1[,1], type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = c(5,25,75,95),med.ln=TRUE,medlab=bquote(~bar('X')),roffset=0.2,rcex=0.7, med.col='white',rpos=2,
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Greys"))),alpha=0.55)
abline(v = 1/ELCResultsMea1[,1],h =seq(0,100,by=10), lty = "dotted", col = "lightgray",lwd=0.8)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = 1/ELCResultsBAU[,1], labels = prettyNum(1/ELCResultsBAU[,1],big.mark = ' '),cex.axis =0.7)
axis(1,line=3,at=1/ELCResultsBAU[,1],labels=TR,tick=TRUE,mgp=c(1.3,0.4,0),cex.axis=0.7)
rug(seq(0,100,by=10), ticksize = -0.015, side = 2)
mtext(expression(bold("Return Period (years)")),1,line=4.5,adj=0.5,cex=0.6)

#Measure 2
par(mar=c(6, 0.2,1.5,2),tck=-0.02)
plot(NULL,xlim=c(min(1/ELCResultsMea2[,1])-0.001,max(1/ELCResultsMea2[,1])),ylim=c(0,90),log='x',xlab='',
     ylab='',cex.lab=0.9,yaxt='n',xaxt='n',cex.main=1.1,main="Elevate houses")
fanJF(data = t(ELCResultsMea2[,-1]/10^6), TR=1/ELCResultsMea2[,1], type = "percentile", probs = seq(0.01,0.99,0.01),
      ln = c(5,25,75,95),med.ln=TRUE,medlab=bquote(~bar('X')),roffset=0.2,rcex=0.7, med.col='white',rpos=2,
      fan.col = colorRampPalette(colors = rev(brewer.pal(9, "Greys"))),alpha=0.55)
abline(v = 1/ELCResultsMea2[,1],h =seq(0,100,by=10), lty = "dotted", col = "lightgray",lwd=0.8)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = 1/ELCResultsBAU[,1], labels = prettyNum(1/ELCResultsBAU[,1],big.mark = ' '),cex.axis =0.7)
axis(4, at = seq(0,100,by=10), labels = prettyNum(seq(0,100,by=10),big.mark = ' '),cex.axis=0.7)
axis(1,line=3,at=1/ELCResultsBAU[,1],labels=TR,tick=TRUE,mgp=c(1.3,0.4,0),cex.axis=0.7)
rug(seq(0,100,by=10), ticksize = -0.015, side = 2)
dev.off()

##Risk - EAD

dev.off()
png("EAD - BAU and Measures.png",width = 120, height = 60, units='mm', res=300)
par(mgp=c(1.3,0.4,0),mar=c(1.5, 2.5,1.5, 0.5),tck=-0.02,family="Times")
rangecolor <- rgb(244,165,130,alpha=255,maxColorValue = 255)
boxplot(EADAllResults,main =NA ,
        ylab=paste("EAD (","\u20AC"," Million/year)"), border = NA,
        xaxt='n', yaxt = "n", frame = FALSE,col=NA,ylim=c(0,9.2),cex.lab=0.7,cex.main=0.8)
abline(h =seq(0,9,by=0.5), lty = "dotted", col = "lightgray")
axis(2, at = seq(0,9,by=1), labels = prettyNum(seq(0,9,by=1),big.mark = ' '),cex.axis=0.6)
rug(seq(0.5,9,by=1), ticksize = -0.015, side = 2)
boxplot(EADAllResults,yaxt = "n",ann = FALSE, outcex=0.6,names=c("BAU","Raise dike","Elevate houses"),cex.axis=0.7,
        border=c(1,"#92C5DE","#F4A582"),boxlwd = 1,col=NA,add = TRUE,ylim=c(0,9.2))

text(x=1,y=mean(EADAllResults[,1]),labels=bquote(bar('EAD')~'='~.(round(mean(EADAllResults[,1]),digits = 1))),cex= 0.6, pos=3,offset=0.05,family="Times")
text(x=2,y=mean(EADAllResults[,2]),labels=bquote(bar('EAD')~'='~.(round(mean(EADAllResults[,2]),digits = 1))),cex= 0.6, pos=3,offset=0.15,family="Times")
text(x=3,y=mean(EADAllResults[,3]),labels=bquote(bar('EAD')~'='~.(round(mean(EADAllResults[,3]),digits = 1))),cex= 0.6, pos=3,offset=0.05,family="Times")
dev.off()


###Utility

dev.off()
png("Utility - BAU and Measures.png",width = 120, height = 60, units='mm', res=300)
par(mgp=c(1.3,0.4,0),mar=c(1.5, 2.5,1.5, 0.5),tck=-0.02,family="Times")
rangecolor <- rgb(244,165,130,alpha=255,maxColorValue = 255)
boxplot(utilityAll,main = NA ,
        ylab=paste("Utility (","\u20AC"," Million)"), border = NA,
        xaxt='n', yaxt = "n", frame = FALSE,col=NA,ylim=c(-10,50),cex.lab=0.7,cex.main=0.8)
abline(h =seq(-10,50,by=5), lty = "dotted", col = "lightgray")
axis(2, at = seq(-10,50,by=10), labels = prettyNum(seq(-10,50,by=10),big.mark = ' '),cex.axis=0.6)
rug(seq(-5,45,by=10), ticksize = -0.015, side = 2)
boxplot(utilityAll,yaxt = "n",ann = FALSE, outcex=0.6,names=c("BAU","Raise dike","Elevate houses"),cex.axis=0.7,
        border=c(1,"#92C5DE","#F4A582"),boxlwd = 1,col=NA,add = TRUE,ylim=c(0,9.2))

text(x=1,y=0,labels=bquote(bar('u')~'='~.(round(mean(utilityAll[,1]),digits = 1))),cex= 0.6, pos=3,offset=0.1,family="Times")
text(x=2,y=mean(utilityAll[,2]),labels=bquote(bar('u')~'='~.(round(mean(utilityAll[,2]),digits = 1))),cex= 0.6, pos=3,offset=0.05,family="Times")
text(x=3,y=mean(utilityAll[,3]),labels=bquote(bar('u')~'='~.(round(mean(utilityAll[,3]),digits = 1))),cex= 0.6, pos=3,offset=0.05,family="Times")
dev.off()

###Costs and discount rate analysis------------------------------------------------------------------------------------------------------
#Define costs per measure (constant per simulation but it should vary)
cMea1<-10
cMea2<-1
discRate<-0.12

###VOI analysis------------------------------------------------------------------------------------------------------
##Results without uncertainty---------------------------------------------------------------------------------------
#Example without uncertainty
##Hazard
#BAU
qFreq0Unc<-qgev(1-1/TR,mu=562,sigma=160,xi=0.02423)
hFreq0Unc<-(qFreq0Unc/56)^(1/1.7)
inflowBAU0Unc<-1.1*(hFreq0Unc-dikeLevelBAU)^1.5
inflowBAU0Unc[is.nan(inflowBAU0Unc)] <- 0
floodLevelBAU0Unc<-inflowBAU0Unc*lengthOverflow*timePeak/(urbanArea)
#Measure 1
inflowMea10Unc<-1.1*(hFreq0Unc-dikeLevelMea1)^1.5
inflowMea10Unc[is.nan(inflowMea10Unc)] <- 0
floodLevelMea10Unc<-inflowMea10Unc*lengthOverflow*timePeak/(urbanArea)
#Measure 2
floodLevelMea20Unc<-floodLevelBAU0Unc
floodLevelMea20Unc[floodLevelMea20Unc<0.5]<-0
floodLevelMea20Unc[floodLevelMea20Unc>=0.5]<-floodLevelMea20Unc[floodLevelMea20Unc>=0.5]-0.5

##Vulnerability
L <- Lnorm(meanlog=-0.69,sdlog=0.95)

##Exposure
priceMC0Unc<-(250+1700+950)/3*numHouses*areaHouse

##Risk
#BAU
ELCResultsBAU0Unc=p(L)(floodLevelBAU0Unc)*priceMC0Unc
EADResultsBAU0Unc<-trapz(rev(1/TR),rev(ELCResultsBAU0Unc))/10^6
#Measure 1
ELCResultsMea10Unc=p(L)(floodLevelMea10Unc)*priceMC0Unc
EADResultsMea10Unc<-trapz(rev(1/TR),rev(ELCResultsMea10Unc))/10^6
#Measure 2
ELCResultsMea20Unc=p(L)(floodLevelMea20Unc)*priceMC0Unc
EADResultsMea20Unc<-trapz(rev(1/TR),rev(ELCResultsMea20Unc))/10^6

##Utility
EADAllResults0Unc<-cbind(EADResultsBAU0Unc,EADResultsMea10Unc,EADResultsMea20Unc)
utilityAll0Unc<-(EADAllResults0Unc[1]-EADAllResults0Unc)/discRate-cbind(0,cMea1,cMea2)
maxUtility0Unc<-max(utilityAll0Unc)
selectOption0Unc<-which(utilityAll0Unc==maxUtility0Unc)-1

####Selected option
VOIhist <- hist(VOIResults[,2]) # or hist(x,plot=FALSE) to avoid the plot of the histogram
VOIhist$density <- VOIhist$counts/sum(VOIhist$counts)*100
#vALUES OBTAINED FROM voiHIST 
perOptions<-replicate(2,c(1.3,64,34.7))

dev.off()
png("Selectoption.png",width = 110, height = 90, units='mm', res=300)
barplot(perOptions,horiz=TRUE,width=c(0.01,0),col=c("#F7F7F7","#92C5DE","#F4A582"))
dev.off()

dev.off()
png("EVPI and EVPPI plots.png",width = 110, height = 90, units='mm', res=300)
par(family="Times")
par(mgp=c(1.3,0.4,0))
par(mar=c(2.5,2,0.5,0.5),tck=-0.02)
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
barplot(c(EVPPIredExpUnMeth1,EVPPIredVulUnMeth1,EVPPIredHazUnMeth1,EVPIMeth1),horiz=TRUE,width=c(0.01,0.01,0.01,0.01),
        names.arg=rev(c("EVPI", expression("EVPPI\nHazard"), expression("EVPPI\nVulnerability"),expression("EVPPI\nExposure"))),
        cex.lab=0.6,cex.axis=0.6,cex.names=0.6,xlab=paste("Value (","\u20AC"," Million)"),xlim = c(0,1.5),
        xaxt='n',col=c("#92C5DE","#F4A582","#2166AC","#D6604D"),border=NA)
axis(1,at=seq(0,1.4,by=0.1),labels = prettyNum(seq(0,1.4,by=0.1),big.mark = ' '),cex.axis=0.6)
dev.off()

##EVPI----------------------------------------------------------------------------

#Organize all results to have EAD per simulation and costs
EADAllResults<-cbind(EADResultsBAU,EADResultsMea1,EADResultsMea2)
colnames(EADAllResults)<-c('EAD_BAU','EAD_Measure1','EAD_Measure2')
costAll<-cbind(rep(0.0001,numIt),rep(cMea1,numIt),rep(cMea2,numIt))
colnames(costAll)<-c('Cost_BAU','Cost_Measure1','Cost_Measure2')

#Calculate utility 
utilityAll<-(EADAllResults[,1]-EADAllResults)/discRate-costAll
colnames(utilityAll)<-c('Util_BAU','Util_Measure1','util_Measure2')

#Selected option - utility current knowledge
maxUtilityCK<-max(colMeans(utilityAll))
selectOption<-which(colMeans(utilityAll)==maxUtilityCK)-1

#Calculate max utility per simulation, selected option, opportunity loss and VOI
VOIResults<-matrix(0, nrow = numIt, ncol = 4)

for(i in 1:numIt){
  VOIResults[i,1]<-max(utilityAll[i,])
  VOIResults[i,2]<-which(utilityAll[i,]==VOIResults[i,1])-1
  
  VOIResults[i,3]<-VOIResults[i,1]-utilityAll[i,selectOption+1]
  if(VOIResults[i,3]<0){
    VOIResults[i,3]<-0
  }
  
  VOIResults[i,4]<-VOIResults[i,1]-maxUtilityCK
}
colnames(VOIResults)<-c('Max_Util_Scenario','Sel_Option','OL','VOI')

#Calculate EVPI
EVPIMeth1<-mean(VOIResults[,1])-maxUtilityCK
EVPIMeth2<-mean(VOIResults[,3])

##EVPPI----------------------------------------------------------------------------
##Groups definition - hazard, vulnerability and exposed value
#Obtain nIntxnInt simulations

##Group 1 - Hazard ---
start_time <- Sys.time()

numIt<-100
ELCResultsBAUredHazUn<-array(0, dim=c(length(TR),numIt,numIt)) #[TR,TargetParam,ComplementParam]
EADResultsBAUredHazUn<-matrix(0, nrow=numIt,ncol=numIt) #[TargetParam,ComplementParam]

ELCResultsMea1redHazUn<-array(0, dim=c(length(TR),numIt,numIt)) #[TR,TargetParam,ComplementParam]
EADResultsMea1redHazUn<-matrix(0, nrow=numIt,ncol=numIt) #[TargetParam,ComplementParam]

ELCResultsMea2redHazUn<-array(0, dim=c(length(TR),numIt,numIt)) #[TR,TargetParam,ComplementParam]
EADResultsMea2redHazUn<-matrix(0, nrow=numIt,ncol=numIt) #[TargetParam,ComplementParam]

#Iterating target parameters
for (k in 1:numIt){
  #Iterating complementary parameters
  for(i in 1:numIt){
    #Iterating Return periods
    for (j in 1:length(TR)){
     
      #Estimate BAU risk
      inondVal<-floodLevelBAU[k,j]
      L <- Lnorm(meanlog=parLN[i,1],sdlog=parLN[i,2])
      ELCResultsBAUredHazUn[j,i,k]=p(L)(inondVal)*priceMC[i]
      
      #Estimate Measure 1 Risk
      inondVal<-floodLevelMea1[k,j]
      ELCResultsMea1redHazUn[j,i,k]=p(L)(inondVal)*priceMC[i]
      
      #Estimate Measure 2 Risk
      inondVal<-floodLevelBAU[k,j]
      
      if (inondVal<0.5){
        inondVal=0
      }
      else{
        inondVal-0.5
      }
      
      ELCResultsMea2redHazUn[j,i,k]=p(L)(inondVal)*priceMC[i]
     
      
    }
    
    #Estimate EAD 
    EADResultsBAUredHazUn[k,i]<-trapz(rev(1/TR),rev(ELCResultsBAUredHazUn[,i,k]))/10^6
    EADResultsMea1redHazUn[k,i]<-trapz(rev(1/TR),rev(ELCResultsMea1redHazUn[,i,k]))/10^6
    EADResultsMea2redHazUn[k,i]<-trapz(rev(1/TR),rev(ELCResultsMea2redHazUn[,i,k]))/10^6
  }
}

EADAllredHazUn<-cbind(rowMeans(EADResultsBAUredHazUn),rowMeans(EADResultsMea1redHazUn),rowMeans(EADResultsMea2redHazUn))

#Calculate utility
utilityAllredHazUn<-sweep((EADAllredHazUn[,1]-EADAllredHazUn)/discRate,2,c(0,cMea1,cMea2))
colnames(utilityAllredHazUn)<-c('Util_BAU','Util_Measure1','util_Measure2')


#Selected option - utility current knowledge
maxUtilityCKredHazUn<-max(colMeans(utilityAllredHazUn))
selectOptionredHazUn<-which(colMeans(utilityAllredHazUn)==maxUtilityCKredHazUn)-1

#Calculate max utility per simulation, selected option, opportunity loss and VOI
VOIResultsredHazUn<-matrix(0, nrow = numIt, ncol = 4)

for(i in 1:numIt){
  VOIResultsredHazUn[i,1]<-max(utilityAllredHazUn[i,])
  VOIResultsredHazUn[i,2]<-(which(utilityAllredHazUn[i,]==VOIResultsredHazUn[i,1])-1)[1]
  
  VOIResultsredHazUn[i,3]<-VOIResultsredHazUn[i,1]-utilityAllredHazUn[i,selectOptionredHazUn+1]
  if(VOIResultsredHazUn[i,3]<0){
    VOIResultsredHazUn[i,3]<-0
  }
  
    VOIResultsredHazUn[i,4]<-VOIResultsredHazUn[i,1]-maxUtilityCKredHazUn
}
colnames(VOIResultsredHazUn)<-c('Max_Util_Scenario','Sel_Option','OL','VOI')

#Calculate EVPPI
EVPPIredHazUnMeth1<-mean(VOIResultsredHazUn[,1])-maxUtilityCKredHazUn
EVPPIredHazUnMeth2<-mean(VOIResultsredHazUn[,3])

end_time <- Sys.time()
end_time - start_time

##Group 2 - Vulnerability ------------
start_time <- Sys.time()

numIt<-100
ELCResultsBAUredVulUn<-array(0, dim=c(length(TR),numIt,numIt)) #[TR,TargetParam,ComplementParam]
EADResultsBAUredVulUn<-matrix(0, nrow=numIt,ncol=numIt) #[TargetParam,ComplementParam]

ELCResultsMea1redVulUn<-array(0, dim=c(length(TR),numIt,numIt)) #[TR,TargetParam,ComplementParam]
EADResultsMea1redVulUn<-matrix(0, nrow=numIt,ncol=numIt) #[TargetParam,ComplementParam]

ELCResultsMea2redVulUn<-array(0, dim=c(length(TR),numIt,numIt)) #[TR,TargetParam,ComplementParam]
EADResultsMea2redVulUn<-matrix(0, nrow=numIt,ncol=numIt) #[TargetParam,ComplementParam]

#Iterating target parameters
for (k in 1:numIt){
  #Iterating complementary parameters
  for(i in 1:numIt){
    #Iterating Return periods
    for (j in 1:length(TR)){
      
      #Estimate BAU risk
      inondVal<-floodLevelBAU[i,j]
      L <- Lnorm(meanlog=parLN[k,1],sdlog=parLN[k,2])
      ELCResultsBAUredVulUn[j,i,k]=p(L)(inondVal)*priceMC[i]
      
      #Estimate Measure 1 Risk
      inondVal<-floodLevelMea1[i,j]
      ELCResultsMea1redVulUn[j,i,k]=p(L)(inondVal)*priceMC[i]
      
      #Estimate Measure 2 Risk
      inondVal<-floodLevelBAU[i,j]
      
      if (inondVal<0.5){
        inondVal=0
      }
      else{
        inondVal-0.5
      }
      
      ELCResultsMea2redVulUn[j,i,k]=p(L)(inondVal)*priceMC[i]
      
      
    }
    
    #Estimate EAD 
    EADResultsBAUredVulUn[k,i]<-trapz(rev(1/TR),rev(ELCResultsBAUredVulUn[,i,k]))/10^6
    EADResultsMea1redVulUn[k,i]<-trapz(rev(1/TR),rev(ELCResultsMea1redVulUn[,i,k]))/10^6
    EADResultsMea2redVulUn[k,i]<-trapz(rev(1/TR),rev(ELCResultsMea2redVulUn[,i,k]))/10^6
  }
}

EADAllredVulUn<-cbind(rowMeans(EADResultsBAUredVulUn),rowMeans(EADResultsMea1redVulUn),rowMeans(EADResultsMea2redVulUn))

#Calculate utility
utilityAllredVulUn<-sweep((EADAllredVulUn[,1]-EADAllredVulUn)/discRate,2,c(0,cMea1,cMea2))
colnames(utilityAllredVulUn)<-c('Util_BAU','Util_Measure1','util_Measure2')


#Selected option - utility current knowledge
maxUtilityCKredVulUn<-max(colMeans(utilityAllredVulUn))
selectOptionredVulUn<-which(colMeans(utilityAllredVulUn)==maxUtilityCKredVulUn)-1

#Calculate max utility per simulation, selected option, opportunity loss and VOI
VOIResultsredVulUn<-matrix(0, nrow = numIt, ncol = 4)

for(i in 1:numIt){
  VOIResultsredVulUn[i,1]<-max(utilityAllredVulUn[i,])
  VOIResultsredVulUn[i,2]<-which(utilityAllredVulUn[i,]==VOIResultsredVulUn[i,1])-1
  
  VOIResultsredVulUn[i,3]<-VOIResultsredVulUn[i,1]-utilityAllredVulUn[i,selectOptionredVulUn+1]
  if(VOIResultsredVulUn[i,3]<0){
    VOIResultsredVulUn[i,3]<-0
  }
  
  VOIResultsredVulUn[i,4]<-VOIResultsredVulUn[i,1]-maxUtilityCKredVulUn
}
colnames(VOIResultsredVulUn)<-c('Max_Util_Scenario','Sel_Option','OL','VOI')

#Calculate EVPPI
EVPPIredVulUnMeth1<-mean(VOIResultsredVulUn[,1])-maxUtilityCKredVulUn
EVPPIredVulUnMeth2<-mean(VOIResultsredVulUn[,3])

end_time <- Sys.time()
end_time - start_time

##Group 3 - Exposed value ------------
start_time <- Sys.time()

numIt<-100
ELCResultsBAUredExpUn<-array(0, dim=c(length(TR),numIt,numIt)) #[TR,TargetParam,ComplementParam]
EADResultsBAUredExpUn<-matrix(0, nrow=numIt,ncol=numIt) #[TargetParam,ComplementParam]

ELCResultsMea1redExpUn<-array(0, dim=c(length(TR),numIt,numIt)) #[TR,TargetParam,ComplementParam]
EADResultsMea1redExpUn<-matrix(0, nrow=numIt,ncol=numIt) #[TargetParam,ComplementParam]

ELCResultsMea2redExpUn<-array(0, dim=c(length(TR),numIt,numIt)) #[TR,TargetParam,ComplementParam]
EADResultsMea2redExpUn<-matrix(0, nrow=numIt,ncol=numIt) #[TargetParam,ComplementParam]

#Iterating target parameters
for (k in 1:numIt){
  #Iterating complementary parameters
  for(i in 1:numIt){
    #Iterating Return periods
    for (j in 1:length(TR)){
      
      #Estimate BAU risk
      inondVal<-floodLevelBAU[i,j]
      L <- Lnorm(meanlog=parLN[i,1],sdlog=parLN[i,2])
      ELCResultsBAUredExpUn[j,i,k]=p(L)(inondVal)*priceMC[k]
      
      #Estimate Measure 1 Risk
      inondVal<-floodLevelMea1[i,j]
      ELCResultsMea1redExpUn[j,i,k]=p(L)(inondVal)*priceMC[k]
      
      #Estimate Measure 2 Risk
      inondVal<-floodLevelBAU[i,j]
      
      if (inondVal<0.5){
        inondVal=0
      }
      else{
        inondVal-0.5
      }
      
      ELCResultsMea2redExpUn[j,i,k]=p(L)(inondVal)*priceMC[k]
      
      
    }
    
    #Estimate EAD 
    EADResultsBAUredExpUn[k,i]<-trapz(rev(1/TR),rev(ELCResultsBAUredExpUn[,i,k]))/10^6
    EADResultsMea1redExpUn[k,i]<-trapz(rev(1/TR),rev(ELCResultsMea1redExpUn[,i,k]))/10^6
    EADResultsMea2redExpUn[k,i]<-trapz(rev(1/TR),rev(ELCResultsMea2redExpUn[,i,k]))/10^6
  }
}

EADAllredExpUn<-cbind(rowMeans(EADResultsBAUredExpUn),rowMeans(EADResultsMea1redExpUn),rowMeans(EADResultsMea2redExpUn))

#Calculate utility
utilityAllredExpUn<-sweep((EADAllredExpUn[,1]-EADAllredExpUn)/discRate,2,c(0,cMea1,cMea2))
colnames(utilityAllredExpUn)<-c('Util_BAU','Util_Measure1','util_Measure2')


#Selected option - utility current knowledge
maxUtilityCKredExpUn<-max(colMeans(utilityAllredExpUn))
selectOptionredExpUn<-which(colMeans(utilityAllredExpUn)==maxUtilityCKredExpUn)-1

#Calculate max utility per simulation, selected option, opportunity loss and VOI
VOIResultsredExpUn<-matrix(0, nrow = numIt, ncol = 4)

for(i in 1:numIt){
  VOIResultsredExpUn[i,1]<-max(utilityAllredExpUn[i,])
  VOIResultsredExpUn[i,2]<-which(utilityAllredExpUn[i,]==VOIResultsredExpUn[i,1])-1
  
  VOIResultsredExpUn[i,3]<-VOIResultsredExpUn[i,1]-utilityAllredExpUn[i,selectOptionredExpUn+1]
  if(VOIResultsredExpUn[i,3]<0){
    VOIResultsredExpUn[i,3]<-0
  }
  
  VOIResultsredExpUn[i,4]<-VOIResultsredExpUn[i,1]-maxUtilityCKredExpUn
}
colnames(VOIResultsredExpUn)<-c('Max_Util_Scenario','Sel_Option','OL','VOI')

#Calculate EVPPI
EVPPIredExpUnMeth1<-mean(VOIResultsredExpUn[,1])-maxUtilityCKredExpUn
EVPPIredExpUnMeth2<-mean(VOIResultsredExpUn[,3])

end_time <- Sys.time()
end_time - start_time


##Plot EVPPI


dev.off()
png("EVPI and EVPPI plots.png",width = 110, height = 90, units='mm', res=300)
par(family="Times")
par(mgp=c(1.3,0.4,0))
par(mar=c(2.5,2,0.5,0.5),tck=-0.02)
rangecolor <- rgb(235,219,199,alpha=175,maxColorValue = 255)
barplot(c(EVPPIredExpUnMeth1,EVPPIredVulUnMeth1,EVPPIredHazUnMeth1,EVPIMeth1),horiz=TRUE,width=c(0.01,0.01,0.01,0.01),
        names.arg=rev(c("EVPI", expression("EVPPI\nHazard"), expression("EVPPI\nVulnerability"),expression("EVPPI\nExposure"))),
        cex.lab=0.6,cex.axis=0.6,cex.names=0.6,xlab=paste("Value (","\u20AC"," Million)"),xlim = c(0,1.5),
        xaxt='n',col=c("#92C5DE","#F4A582","#2166AC","#D6604D"),border=NA)
axis(1,at=seq(0,1.4,by=0.1),labels = prettyNum(seq(0,1.4,by=0.1),big.mark = ' '),cex.axis=0.6)
dev.off()





matplot(vulResults[,1],vulResults[,-1],type='l',lwd = 0.01,xlab='Water Depth (m)',
        ylab='Damage Factor (-)',cex.lab=0.7,cex.main=0.9,
        main = "Vulnerability Function - BAU", ylim = c(0,1),yaxt='n',xaxt='n',
        panel.first=abline(v = vulResults[,1],h =seq(0,1,by=0.1), lty = "dotted", col = "lightgray",
                           panel.first=polygon(x=c(vulResults[,1],rev(vulResults[,1])),y=c(apply(vulResults[,-1],1,min),rev(apply(vulResults[,-1],1,max))),
                                               border = NA, col=rangecolor))
)
#axis(1, at = TR, labels = prettyNum(TR,big.mark = ' '),cex.axis=0.6)
mgp.axis.labels(c(0.4,0.4,0), type='x')
mgp.axis(1, at = vulResults[,1], labels = prettyNum(vulResults[,1],big.mark = ' '),cex.axis =0.6)
axis(2, at = seq(0,1,by=0.1), labels = prettyNum(seq(0,1,by=0.1),big.mark = ' '),cex.axis=0.6)
dev.off()

###Sensitivity analysis ------------------------------------------------------------
##EVPI
#Define costs per measure (constant per simulation but it should vary)
cMea1Sens<-seq(0.1,30,by=0.5)
cMea2Sens<-seq(0.1,30,by=0.5)
cComb<-expand.grid(cMea1Sens,cMea2Sens)
cComb<-data.matrix(cComb)
numComb<-nrow(cComb)
discRate<-0.12

#Function to estimate EVPI
estimateEVPI<- function(EAD, costs, rRate) {
  #Calculate utility
  utility<-sweep((EAD[,1]-EAD)/rRate,2,costs)

  #Selected option - utility current knowledge
  maxUtility<-max(colMeans(utility))
  selectOption<-which(colMeans(utility)==maxUtility)-1
  
  
  #Calculate max utility per simulation, selected option, opportunity loss and VOI
  
  VOIResults<-matrix(0, nrow = nrow(utility), ncol = 4)
  
  for(i in 1:nrow(utility)){
    

    VOIResults[i,1]<-max(utility[i,])
    
    
    VOIResults[i,2]<-(which(utility[i,]==VOIResults[i,1])-1)[1]
    
    VOIResults[i,3]<-VOIResults[i,1]-utility[i,selectOption+1]
    if(VOIResults[i,3]<0){
      VOIResults[i,3]<-0
    }
    #browser()
    VOIResults[i,4]<-VOIResults[i,1]-maxUtility
    
  }

  
  #Calculate EVPI
  EVPI<-mean(VOIResults[,1])-maxUtility
  
  return(c(maxUtility,selectOption,EVPI))
 
}

#Iterate per combination of costs to calculate EVPI
EVPISens<-matrix(0,nrow = numComb,ncol = 3)

for (i in 1:numComb){
 
  EVPISens[i,]<-estimateEVPI(EADAllResults,c(0,cComb[i,]),discRate)
  
}


#Plot EVPI vs costs
dev.off() 
png("EVPI-Cost sensitivity analysis.png",width = 150, height = 100, units='mm', res=300)
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
par(family="Times")
filled.contour(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1),matrix(EVPISens[,3],length(cMea1Sens), byrow = TRUE),
              plot.title = {par(cex.main=1.2,cex.lab=0.8);title(main ='EVPI for different Cost Combinations',xlab = paste("Measure 2 - Costs (","\u20AC"," Million)"), ylab = paste("Measure 1 - Costs (","\u20AC"," Million)"))},
              nlevels = 9,col =hcl.colors(10, palette = "Reds",rev = TRUE),
              key.title = {par(cex.main=0.8,mar=c(5.1,0.1,2.5,1),cex.axis=0.8);title(main=paste('EVPI(',"\u20AC",'Million )'))})
dev.off()

#Plot Expected utility vs costs
dev.off() 
png("Expected utility - Cost sensitivity analysis.png",width = 150, height = 100, units='mm', res=300)
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
par(family="Times")
filled.contour(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1),matrix(EVPISens[,1],length(cMea1Sens), byrow = TRUE),
               plot.title = {par(cex.main=1,cex.lab=0.8);title(main ='Expected Utility for different Cost Combinations',xlab = paste("Measure 2 - Costs (","\u20AC"," Million)"), ylab = paste("Measure 1 - Costs (","\u20AC"," Million)"))},
               nlevels = 9,col = hcl.colors(10, palette = "Blues",rev = TRUE),
               key.title = {par(cex.main=0.6,mar=c(5.1,0.1,2.5,1.4),cex.axis=0.8);title(main=paste('Utility (',"\u20AC",'Million )'))})
dev.off()


#Plot Selected option vs costs
dev.off() 
png("Selected option - Cost sensitivity analysis.png",width = 150, height = 100, units='mm', res=300)
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
par(family="Times")
filled.contour(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1),matrix(EVPISens[,1],length(cMea1Sens), byrow = TRUE),
               plot.title = {par(cex.main=1,cex.lab=0.8);title(main ='Expected Utility for different Cost Combinations',xlab = paste("Measure 2 - Costs (","\u20AC"," Million)"), ylab = paste("Measure 1 - Costs (","\u20AC"," Million)"))},
               nlevels = 9,col = hcl.colors(10, palette = "Blues",rev = TRUE),
               key.title = {par(cex.main=0.6,mar=c(5.1,0.1,2.5,1.4),cex.axis=0.8);title(main=paste('Utility (',"\u20AC",'Million )'))})
dev.off()


#Plot ratio EVPI
dev.off() 

EVPIMaxUtilSens<-expUtilSens
EVPIMaxUtilSens$z<-rep(0,nrow(EVPISens))
EVPIMaxUtilSens$z[EVPISens[,3]==0]<-0
EVPIMaxUtilSens$z[((EVPISens[,3]-EVPISens[,1])<0) & EVPISens[,3]!=0 ]<-1
EVPIMaxUtilSens$z[(EVPISens[,3]-EVPISens[,1])>0]<-2
  

p5<-ggplot(EVPIMaxUtilSens, aes(Var1, Var2, z = z)) + 
  geom_raster(aes(fill = z),interpolate = TRUE) + 
  labs(x = paste("Raise dike - Costs (","\u20AC"," Million)"), 
       y = paste("Elevate houses - Costs (","\u20AC"," Million)"), 
       title="Benefits of aquiring information based on EVPI") + 
  theme(plot.title = element_text(size=7,face='bold',hjust = 0.5),
        legend.title = element_blank(), 
        legend.position = "top", panel.background = element_blank(), 
        axis.text = element_text(colour = "black", size = 6), 
        axis.title = element_text(size = 7),
        axis.ticks = element_line(size=0.5),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,-5,-5,-5),
        legend.text = element_text(size = 4), 
        legend.key.width = unit(0.1, "cm") ,
        legend.key.height = unit(0.1, "cm"),
        plot.margin=margin(t = 0, r = 2, b = 2, l = 2, unit = "pt"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  scale_fill_gradientn(colours=c("#F7F7F7","#92C5DE","#F4A582"),guide = "legend",breaks=c(0,1,2),labels=c('Not beneficial','Potentially beneficial','Highly potentially beneficial')) +
  scale_y_continuous(breaks=seq(0,30,by=5),expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,30,by=5),expand=c(0,0))+
  geom_point(data = data.frame(Var1=10,Var2=1,z=9.57), col = 'black',shape=4,size=1.5)


ggsave(filename = "Ratio EVPI - MAx utility- Cost sensitivity analysis.png", p5,
       width = 60, height = 70, dpi = 300, units = "mm", device='png')

###Plot all together--------------
#Plot Expected utility vs costs
dev.off() 
png("Expected utility - Cost sensitivity analysis.png",width = 150, height = 100, units='mm', res=300)
par(mgp=c(2,0.8,0))
par(mar=c(4.1, 4.1,2, 2.1))
par(family="Times")
contour(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1),t(matrix(EVPISens[,1],length(cMea1Sens), byrow = TRUE)),
               plot.title = {par(cex.main=1,cex.lab=0.8);title(main ='Maximum Expected Utility for different Cost Combinations',ylab = paste("Elevate houses - Costs (","\u20AC"," Million)"), xlab = paste("Raise dike - Costs (","\u20AC"," Million)"))},
               nlevels = 9,col = hcl.colors(10, palette = "Blues",rev = TRUE),
               key.title = {par(cex.main=0.6,mar=c(5.1,0.1,2.5,1.4),cex.axis=0.8);title(main=paste('Utility (',"\u20AC",'Million )'))})
dev.off()

expUtilSens<-expand.grid(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1))
expUtilSens$z<-EVPISens[,1]
voiUtilSens<-expUtilSens
voiUtilSens$z<-EVPISens[,3]
selUtilSens<-expUtilSens
selUtilSens$z<-EVPISens[,2]

windowsFonts(Times=windowsFont("Times New Roman"))
theme_set(theme_minimal(base_family = "Times"))
p1<-ggplot(expUtilSens, aes(Var1, Var2, z = z)) + 
  geom_raster(aes(fill = z),interpolate = TRUE) + 
  geom_contour(aes(z = z), colour = "grey50", size = 0.3, alpha = 0.5) + 
  geom_text_contour(aes(z = z),  colour = "black",size =2,label.placement=label_placement_fraction(0.4)) +
  labs(x = paste("Raise dike - Costs (","\u20AC"," Million)"), 
       y = paste("Elevate houses - Costs (","\u20AC"," Million)"), 
       fill = paste('Utility (',"\u20AC",'Million )'),
       title="Maximum Expected Utility") + 
  theme(plot.title = element_text(size=9,face='bold',hjust = 0.5),
    legend.title = element_text(size = 8, face = "bold"), 
        legend.position = "bottom", panel.background = element_blank(), 
        axis.text = element_text(colour = "black", size = 7), 
        axis.title = element_text(size = 8),
        axis.ticks = element_line(size=0.5),
    axis.title.x=element_text(margin = margin(t=5)),
    axis.title.y=element_text(margin = margin(r=5)),
    legend.margin=margin(t = 0,b=0),
    legend.box.margin=margin(0,0,0,0),
        legend.text = element_text(size = 7), legend.key = element_blank(),
    legend.key.width = unit(0.5, "cm") ,
    legend.key.height = unit(0.3, "cm"),
    plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  guides(fill = guide_colorbar(title.position = "bottom",title.hjust = 0.5))+
  scale_fill_gradientn(colours=hcl.colors(10, palette = "Blues",rev = TRUE)) + 
  scale_y_continuous(breaks=seq(0,30,by=5),expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0,30,by=5),expand=c(0,0))+
  geom_point(data = data.frame(Var1=10,Var2=1,z=9.57), col = 'black',shape=4,size=1.5)
  

p2<-ggplot(selUtilSens, aes(Var1, Var2, z = z)) + 
  geom_raster(aes(fill = z),interpolate = TRUE) + 
  labs(x = paste("Raise dike - Costs (","\u20AC"," Million)"), 
       y = NULL, 
       fill = paste('Selected Option'),
       title="Current Optimal Option") + 
  theme(plot.title = element_text(size=9,face='bold',hjust = 0.5),
        legend.title = element_text(size = 8, face = "bold"), 
        legend.position = "bottom", panel.background = element_blank(), 
        axis.text = element_text(colour = "black", size = 7), 
        axis.title = element_text(size = 8),
        axis.ticks = element_line(size=0.5),
        axis.ticks.y.right = element_line(size=0.5),
        axis.text.y = element_blank(),
        axis.title.x=element_text(margin = margin(t=8)),
        legend.margin=margin(t = 2.5,b=4),
        legend.text = element_text(size = 7), 
        legend.key.width = unit(0.3, "cm") ,
        legend.key.height = unit(0.3, "cm"),
        plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + 
  scale_fill_gradientn(colours=c("#F7F7F7","#92C5DE","#F4A582"),guide = "legend",breaks=c(0,1,2),labels=c('BAU','Raise dike','Elevate houses')) +
  guides(fill = guide_legend(title.position = "bottom",title.hjust = 0.5))+
  scale_y_continuous(breaks=seq(0,30,by=5),expand=c(0,0),sec.axis =  dup_axis()) +
  scale_x_continuous(breaks=seq(0,30,by=5),expand=c(0,0))+
  geom_point(data = data.frame(Var1=10,Var2=1,z=9.57), col = 'black',shape=4,size=1.5)

p3<-ggplot(voiUtilSens, aes(Var1, Var2, z = z)) + 
  geom_raster(aes(fill = z),interpolate = TRUE) + 
  geom_contour(aes(z = z), colour = "grey50", size = 0.3, alpha = 0.5) + 
  geom_text_contour(aes(z = z),  colour = "black",size =2,label.placement=label_placement_fraction(0.4)) +
  labs(x = paste("Raise dike - Costs (","\u20AC"," Million)"), 
       y = paste("Elevate houses - Costs (","\u20AC"," Million)"), 
       fill = paste('EVPI (',"\u20AC",'Million )'),
       title="EVPI") + 
  theme(plot.title = element_text(size=9,face='bold',hjust = 0.5),
        legend.title = element_text(size = 8, face = "bold"), 
        legend.position = "bottom", panel.background = element_blank(), 
        axis.text = element_text(colour = "black", size = 7), 
        axis.title = element_text(size = 8),
        axis.title.x=element_text(margin = margin(t=5)),
        axis.title.y.right=element_text(margin = margin(l=5)),
        axis.ticks = element_line(size=0.5),
        legend.margin=margin(t = 0,b=0),
        legend.box.margin=margin(0,0,0,0),
        legend.text = element_text(size = 7), legend.key = element_blank(),
        legend.key.width = unit(0.5, "cm") ,
        legend.key.height = unit(0.3, "cm"),
        plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  guides(fill = guide_colorbar(title.position = "bottom",title.hjust = 0.5))+
  scale_fill_gradientn(colours=hcl.colors(10, palette = "Reds",rev = TRUE)) + 
  scale_y_continuous(breaks=seq(0,30,by=5),expand=c(0,0),position = "right") +
  scale_x_continuous(breaks=seq(0,30,by=5),expand=c(0,0))+
  geom_point(data = data.frame(Var1=10,Var2=1,z=9.57), col = 'black',shape=4,size=1.5)

#finalP<-grid.arrange(p1, p2,p3,nrow = 1)

finalP<-plot_grid(p1, p2, p3, nrow = 1, rel_widths  = c(0.35,0.3, 0.35))

ggsave(filename = "Sensitivity-Utility and EVPIr.png", finalP,
       width = 150, height = 100, dpi = 300, units = "mm", device='png')
##EVPPI--------------------------------------------------------------------------------------------------------------------------------------------------
#Iterate per combination of costs to calculate EVPPI per parameters
EVPPISensredHaz<-matrix(0,nrow = numComb,ncol = 3)
EVPPISensredVul<-matrix(0,nrow = numComb,ncol = 3)
EVPPISensredExp<-matrix(0,nrow = numComb,ncol = 3)

for (i in 1:numComb){
  
  EVPPISensredHaz[i,]<-estimateEVPI(EADAllredHazUn,c(0,cComb[i,]),discRate)
  EVPPISensredVul[i,]<-estimateEVPI(EADAllredVulUn,c(0,cComb[i,]),discRate)
  EVPPISensredExp[i,]<-estimateEVPI(EADAllredExpUn,c(0,cComb[i,]),discRate)
  
}

#Plot results for each EVPPI
dev.off()
png("EVPPI - Cost sensitivity analysis.png",width = 160, height = 65, units='mm', res=300)
plot.new()
#Hazard
par(new = "TRUE",mgp=c(0.8,0.15,0),mar=c(0.5, 0.5,0.7,0.5),plt=c(0.06,0.32,0.15,0.95),cex.axis = 0.6,family="Times",tck=-0.02)
filled.contour3(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1),t(matrix(EVPPISensredHaz[,3],length(cMea1Sens), byrow = TRUE)),
               plot.title = {par(cex.main=0.8,cex.lab=0.7);title(main ='Hazard Parameters',xlab = NA, ylab = expression(bold(paste("Elevate houses - Costs (","\u20AC"," Million)"))))},
               nlevels = 7,col =hcl.colors(8, palette = "Reds",rev = TRUE),zlim = c(0,max(EVPPISensredHaz[,3])),
               key.title = {par(cex.main=0.8,mar=c(5.1,0.1,2.5,1),cex.axis=0.8);title(main=paste('EVPI(',"\u20AC",'Million )'))})
grid(col='gray48',lwd=0.1)
points(10,1,pch =4,cex=0.8)

#Vulnerability
par(new = "TRUE",mgp=c(0.8,0.15,0),mar=c(0.5,0.5,0.7,0.5),plt=c(0.33,0.60,0.15,0.95),cex.axis = 0.6,yaxt='n')
filled.contour3(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1),t(matrix(EVPPISensredVul[,3],length(cMea1Sens), byrow = TRUE)),
                plot.title = {par(cex.main=0.8,cex.lab=0.7);title(main ='Vulnerability Parameters',xlab = expression(bold(paste("Raise dike - Costs (","\u20AC"," Million)"))), yaxt = "n",ylab = NA)},
                nlevels = 7,col =hcl.colors(8, palette = "Reds",rev = TRUE),zlim = c(0,max(EVPPISensredHaz[,3])),
                key.title = {par(cex.main=0.8,mar=c(5.1,0.1,2.5,1),cex.axis=0.8);title(main=paste('EVPI(',"\u20AC",'Million )'))})
grid(col='gray48',lwd=0.1)
points(10,1,pch =4,cex=0.8)

#Exposure
par(new = "TRUE",mgp=c(0.8,0.15,0),mar=c(0.5,0.5,0.7,0.5),plt=c(0.61,0.88,0.15,0.95),cex.axis = 0.6)
filled.contour3(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1),t(matrix(EVPPISensredExp[,3],length(cMea1Sens), byrow = TRUE)),
                plot.title = {par(cex.main=0.8,cex.lab=0.7);title(main ='Exposure Parameters',xlab = NA, ylab = NA)},
                nlevels = 7,col =hcl.colors(8, palette = "Reds",rev = TRUE),zlim = c(0,max(EVPPISensredHaz[,3])),
                key.title = {par(cex.main=0.8,mar=c(5.1,0.1,2.5,1),cex.axis=0.8);title(main=paste('EVPI(',"\u20AC",'Million )'))})
grid(col='gray48',lwd=0.1)
points(10,1,pch =4,cex=0.8)

#Legend
par(new = "TRUE",mar=c(0,0,0,0),plt=c(0.9,0.95,0.2,0.8),cex.axis = 0.5,yaxt='s',las=1)

filled.legend(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1),t(matrix(EVPPISensredHaz[,3],length(cMea1Sens), byrow = TRUE)),
                plot.title = {par(cex.main=0.5,cex.lab=0.5);title(main ='EVPI for different Cost Combinations',xlab = paste("Measure 2 - Costs (","\u20AC"," Million)"), ylab = paste("Measure 1 - Costs (","\u20AC"," Million)"))},
                nlevels = 9,col =hcl.colors(10, palette = "Reds",rev = TRUE),zlim = c(0,max(EVPPISensredHaz[,3])),
                key.title = {par(cex.main=0.5,mar=c(5.1,0.1,2.5,1),cex.axis=0.4);title(main=paste('EVPI(',"\u20AC",'Million )'))})
filled.legend(matrix(cMea1Sens,nrow =length(cMea1Sens),ncol=1),matrix(cMea2Sens,nrow =length(cMea2Sens),ncol=1),matrix(EVPPISensredHaz[,3],length(cMea1Sens), byrow = TRUE),
              plot.title = {par(cex.main=0.5,cex.lab=0.5);title(main ='EVPI for different Cost Combinations',xlab = paste("Measure 2 - Costs (","\u20AC"," Million)"), ylab = paste("Measure 1 - Costs (","\u20AC"," Million)"))},
              nlevels = 7,col =hcl.colors(8, palette = "Reds",rev = TRUE),zlim = c(0,max(EVPPISensredHaz[,3])),
              key.title = {par(cex.main=0.5,mar=c(5.1,0.1,2.5,1),cex.axis=0.4);title(main=paste('EVPI(',"\u20AC",'Million )'))})
par(xpd=NA)
text(x = 0.6,y = 3.7,expression(bold(paste('EVPPI\n(\u20AC Million)'))),cex = 0.5)
dev.off()


###########TEST EVPPI additive ---------------------------------------------------
NBAnalysis<-bcea(e=(EADAllResults[,1]-EADAllResults)/0.12,c=costAll,interventions='Measures',wtp=1)
UncParameterInput<-cbind(parGEV,parLN,pricepersqDist)
colnames(UncParameterInput)<-c('GEVLoc','GEVSca','GEVSha','LNMean','LnSd','ExpUnPr')

HazardEVPPI<-evppi(parameter=c(1:3),input=UncParameterInput,he=NBAnalysis)
VulEVPPI<-evppi(parameter=c(4:5),input=UncParameterInput,he=NBAnalysis)
ExpEVPPI<-evppi(parameter=6,input=UncParameterInput,he=NBAnalysis)

TotalEVPPI<-evppi(parameter=c(1:6),input=UncParameterInput,he=NBAnalysis,method='GP')


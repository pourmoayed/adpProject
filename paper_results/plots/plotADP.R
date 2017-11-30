require(tools)
library(reshape2)
library(ggplot2)
library(tikzDevice)
library(grid)
library(plyr)
library(dplyr)
library(tidyr)


#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------Plots for Experiment 2 (plots 5,6,7,8, and 9)-----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------

# Reading data for Experiment 2:
#read the data in the pen levels
numPigsPens <- read.csv2("numPigsPens.csv")
numPigsPens$weight <- as.numeric(as.character(numPigsPens$weight))
numPigsPens$pigWeight <- as.numeric(as.character(numPigsPens$pigWeight))
numPigsPens$growth <- as.numeric(as.character(numPigsPens$growth))
numPigsPens$sd <- as.numeric(as.character(numPigsPens$sd))

numPigsPens$secLabel <- numPigsPens$section
numPigsPens$secLabel[numPigsPens$secLabel==0]<-"Section 1"
numPigsPens$secLabel[numPigsPens$secLabel==1]<-"Section 2"
numPigsPens$secLabel[numPigsPens$secLabel==2]<-"Section 3"

#read the data in section levels
numPigsSecs <- read.csv2("numPigsSecs.csv")
numPigsSecs <- numPigsSecs[139:285,]
#numPigsSecs <- numPigsSecs[178:333,]
#numPigsSecs <- numPigsSecs[166:258,]
#numPigsSecs <- numPigsSecs[0:150,]
numPigsSecs$weight <- as.numeric(as.character(numPigsSecs$weight))
numPigsSecs$growth <- as.numeric(as.character(numPigsSecs$growth))
numPigsSecs$sd <- as.numeric(as.character(numPigsSecs$sd))

numPigsSecs$optimal <- numPigsSecs$opt
numPigsSecs$optimal[numPigsSecs$optimal==0] <- "C"
numPigsSecs$optimal[numPigsSecs$optimal==numPigsSecs$numPigs]<-"T"

numPigsSecs$secLabel <- numPigsSecs$section
numPigsSecs$secLabel[numPigsSecs$secLabel==0]<-"Section 1"
numPigsSecs$secLabel[numPigsSecs$secLabel==1]<-"Section 2"
numPigsSecs$secLabel[numPigsSecs$secLabel==2]<-"Section 3"

# data frame for transportation cost 
# numPigsSecsTr<-numPigsSecs 

# High transportation cost
# numPigsSecsNew <-read.csv2("high_transport_cost/numPigsSecs.csv")
# numPigsSecsNew <-numPigsSecsNew[139:285,]
# numPigsSecsTr$totalCullHigh<-numPigsSecsNew$totalCull
# rm(numPigsSecsNew)

# Low transportation cost
# numPigsSecsNew <-read.csv2("numPigsSecs.csv")
# numPigsSecsNew <-numPigsSecsNew[139:285,]
# numPigsSecsTr$totalCullLow<-numPigsSecsNew$totalCull
# rm(numPigsSecsNew)

#save(file="numPigsSecsTr", numPigsSecsTr)
#write.csv2(numPigsSecsTr,"numPigsSecsTr.csv", row.names = FALSE)


# plot 6
tikz("optADP.tex", width = 10, height = 5.5, standAlone = F)
plot<-ggplot(data=numPigsSecs, aes(x = epoch, y= numPigs, group=factor(secLabel)  ) ) + 
  scale_y_continuous(breaks=c(-20,seq(0,360,40)), label = c("$a^*$", seq(0,360,40)) ) + scale_x_continuous(breaks=seq(46,94,1), labels=46:94 ) +
  geom_bar(aes(fill=numPigs), stat="identity", fill = "dark grey", colour = "black", alpha = 1/3, position=position_dodge(), width=.5) +  
  geom_text(mapping=aes(x=epoch, y=-20, label=optimal), size=2.3) +
  facet_grid( secLabel ~ .) + 
  xlab("Week number") + ylab(" ") 
#g <- guide_legend("",nrow=1,byrow=TRUE, override.aes = list(fill=NA))
#g <- guide_legend("",nrow=1,byrow=TRUE)
plot + theme_bw() + 
  theme(legend.position="none", panel.background = element_blank(), 
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        axis.title.x= element_text(vjust = -0.7),
        axis.text=element_text(size=6),
        strip.background=element_rect(fill = NA)) 
dev.off()
tools::texi2pdf(file = "optADP.tex", clean = T)


#plot 5
tikz("weightInfo.tex", width = 10, height = 6, standAlone = F)
plot<-ggplot(data=numPigsPens, aes(x = t, y= pigWeight, group=factor(t)  ) ) + 
  scale_y_continuous(breaks=c(seq(0,160,10)), label = c(seq(0,160,10)) ) + scale_x_continuous(breaks=seq(0,14,1), labels=1:15 ) +
  geom_boxplot(aes(t) ) + 
  facet_wrap(~ secLabel) + # , scales="free"
  xlab("Week number") + ylab("Weight \\small(kg)") 
plot + theme_bw() + 
  theme(legend.position="none", panel.background = element_blank(), 
          panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
          axis.title.x= element_text(vjust = -0.7),
          strip.background=element_rect(fill = NA)) 
dev.off()
tools::texi2pdf(file = "weightInfo.tex", clean = T)


#plot 7
tikz("totalCull.tex", width = 10, height = 4, standAlone = F)
plot<-ggplot(data=numPigsSecs, aes(x = epoch, y= totalCull ) ) + 
  scale_y_continuous(breaks=c(-20,seq(0,400,50)), label = c("$z^*$", seq(0,400,50)) ) + scale_x_continuous(breaks=seq(46,94,1), labels=46:94 ) +
  geom_bar(aes(fill=totalCull), stat="identity", fill = "dark grey", colour = "black", alpha = 1/3, position=position_dodge(), width=.5) +  
  #geom_line() + #geom_point() +
  geom_hline(yintercept=205, linetype = 2) + 
  geom_text(mapping=aes(x=epoch, y=-20, label= truckNum), size=2) +
  xlab("Week number") + ylab(" ") 
plot + theme_bw() + 
  theme(legend.position="none", panel.background = element_blank(), 
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        axis.title.x= element_text(vjust = -0.7),
        axis.text=element_text(size=6),
        strip.background=element_rect(fill = NA)) 
dev.off()
tools::texi2pdf(file = "totalCull.tex", clean = T)

#plot 8
datTruck<-melt(data.frame(numPigsSecsTr),
          # ID variables - all the variables to keep but not split apart on
          id.vars="epoch",
          # The source columns
          measure.vars=c("totalCull","totalCullHigh"),
          # Name of the destination column that will identify the original
          # column that the measurement came from
          variable.name="name",
          value.name="y"
)
datTruck$name<-mapvalues(datTruck$name, from = c("totalCull","totalCullHigh"), to = c("$c^{\\mathtt{\\small truck}}=400$ DKK","$c^{\\mathtt{\\small truck}}=2000$ Dkk") )

tikz("truckCompare.tex", width = 10, height = 4, standAlone = F)
plot<-ggplot(data=datTruck, aes(x = epoch, y = y, group = name, fill =name ) ) + 
  scale_y_continuous(breaks=c(-20,seq(0,400,50)), label = c("$z^*$", seq(0,400,50)) ) + scale_x_continuous(breaks=seq(46,94,1), labels=46:94 ) +
  geom_bar(stat="identity", alpha = 1/3, position=position_dodge(), width=0.5) + 
  scale_fill_manual(values=c("gray70","black") ) + 
  #scale_y_continuous(breaks=c(-20,seq(0,400,50)), label = c("$z^*$", seq(0,400,50)) ) + scale_x_continuous(breaks=seq(22,73,1), labels=22:73 ) +
  #geom_histogram(stat="identity", data=datTruck, alpha = 1/4, colour=NA, width=0.25)+
  #geom_bar(aes(fill=totalCullLow), stat="identity", fill = "dark grey", colour = "black", alpha = 1/3, position=position_dodge(), width=0.5) +  
  #geom_bar(aes(fill=totalCullHigh), stat="identity", fill = "gray", colour = "black", alpha = 1/3, position=position_dodge(), width=.5) +  
  #geom_line() + #geom_point() +
  geom_hline(yintercept=205, linetype = 2) + 
  #geom_text(mapping=aes(x=epoch, y=-20, label= truckNum), size=2) +
  xlab("Week number") + ylab(" ")
  g <- guide_legend("", override.aes = list(fill=NA))
plot  + guides(shape = g, linetype=g) + theme_bw() + 
  theme( panel.background = element_blank(), 
        panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        legend.key = element_rect(fill = NA, colour = NA),
        #legend.key.width = unit(0.5, "cm"),
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_blank(),
        axis.title.x= element_text(vjust = -0.7),
        axis.text=element_text(size=6),
        strip.background=element_rect(fill = NA)) 
dev.off()
tools::texi2pdf(file = "truckCompare.tex", clean = T)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------Plots for Experiment 1 (plots 2,3, and 4)----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------

#read data for experiment 1

# Find optimal policy at pen level using R package in the first paper
library(hmdpFeedPigIT)
set.seed(234567)
# Set HMDP parameters: 
param<-setParameters(pigs=18, 
                     rations=1,
                     phases=1,
                     tMax=15,
                     tStartMarketing=9,
                     iniFeedMix=1,
                     minPhaseT=c(1),
                     disWeight=c(10,2),
                     disSD=c(4,1),
                     disGrowth=c(2,0.3),
                     priorGrowth=c(6),
                     cullActions = T
)
# Create GSSM
dlms<-list()
for (i in 1:param$rations){
  dlmMod<-iniDLM()
  Y<-matrix(c(28.91122,11.17050, 34.73722, 11.57909, 40.55647, 11.91425, 46.43588, 12.31422, 52.60592, 12.76750, 58.90283, 13.18425, 64.73567, 13.57245, 70.88710, 13.94777, 76.19884, 14.36087, 82.30708, 14.76617, 88.41322, 15.10727, 94.68169,  15.56838, 100.68169, 15.96838, 106.68169, 16.57838, 112.68169, 17.06838), nrow=2)
  dlms[[i]]<-buildDLM(dlmMod, Y)
}
# Get the nGSSM
dglmParam<-iniDGLM()
# Build the HMDP model
prefix<-"hmdp_"
BuildHMDP(prefix, param, dlms, dglmParam)
# Solve the HMDP model
wLbl<-"Reward"
durLbl<-"Time"
mdp<-loadMDP(prefix)
rate<-0.05              # discount rate
rateBase<-1             # rate base
valueIte(mdp, wLbl, durLbl, rate, rateBase, maxIte = 80000, eps = 10^-40)
#g<-policyIteAve(mdp,wLbl,durLbl)      # Finds the optimal policy using the avage reward per week (g) criterion 
policy<-getPolicy(mdp) # optimal policy for each sId
do.call(file.remove,list(list.files(pattern = prefix)))
rm(mdp)

# Calculate slope values from optimal policy at pen level: 
require(stringr)
policyExternalM<-policy[1:131545,]
aaa<-do.call(rbind, str_split(gsub( "\\(|\\)" , "" , unlist(policyExternalM$stateLabel) ),","))
policyExternalM$n<-as.integer(aaa[,1])
policyExternalM$iSW<-as.integer(aaa[,2])
policyExternalM$iSG<-as.integer(aaa[,3])
policyExternalM$iSSd<-as.integer(aaa[,4])
policyExternalM$t<-as.integer(aaa[,5])
# policyExternalM$RS<-as.integer(aaa[,5])
# policyExternalM$phase<-as.integer(aaa[,6])
# policyExternalM$iRation<-as.integer(aaa[,7])
# policyExternalM$t<-as.integer(aaa[,8])
slopeHMDP<-matrix(nrow = 6804,ncol = 6)
colnames(slopeHMDP) <- c("t","iL","iG","iSd", "slope", "intercept")
i<-1
sumError<-0
counter<-0
for(tt in 10:15){
  for(iSWw in 0:20){
    for(iSSdd in 0:8){
      for(iSGg in 0:5){
        #policyFilter = subset(policyExternalM,  iSW==iSWw &  iSG==iSGg & iSSd==iSSdd  & RS==1 & phase==0 & iRation==0 & t==tt )
        policyFilter = subset(policyExternalM,  iSW==iSWw &  iSG==iSGg & iSSd==iSSdd  & t==tt )
        #lines(policyFilter$n, policyFilter$weight, col=sample(rainbow(10)) )
        x<-policyFilter$n
        y<-policyFilter$weight
        dd<-data.frame(
          x,
          y
        )
        fit <- lm(formula=y~poly(x=x,degree=1,raw=TRUE),data=dd)
        #fitL <- lm(formula=y~0+x,data=dd)
        #h<-predict(fit, newdata=data.frame(x=x))
        #hL<-predict(fitL, newdata=data.frame(x=x))
        #sumError <- sqrt( sum( (hL-h)^2 ) )/ length(x) + sumError
        sumError <- summary(fit)$adj.r.squared + sumError
        counter <- counter + 1
        #lines(x, predict(fit, newdata=data.frame(x=x)))
        slope <- summary(fit)$coefficients[2, 1]
        intercept <- summary(fit)$coefficients[1, 1]
        slopeHMDP[i,1] = tt #-1
        slopeHMDP[i,2] = iSWw
        slopeHMDP[i,3] = iSGg
        slopeHMDP[i,4] = iSSdd
        slopeHMDP[i,5] = round(as.numeric(slope),3)
        slopeHMDP[i,6] = round(intercept,3)
        i<-i+1
      }
    }
  }
}
#sumError/counter
dataADP <- read.csv2("dataADP.csv")
slopeHMDP <- as.data.frame(slopeHMDP)
write.csv2(slopeHMDP,"slopeHMDPNewT.csv", row.names = FALSE)
dataADP$t<-as.numeric(as.character(dataADP$t))
dataADP$iL<-as.numeric(as.character(dataADP$iL))
dataADP$iG<-as.numeric(as.character(dataADP$iG))
dataADP$iSd<-as.numeric(as.character(dataADP$iSd))
dataADP$slope<-as.numeric(as.character(dataADP$slope))
sHMDP<-c()
for( z in 1:dim(dataADP)[1])
  sHMDP[z] <- unlist( subset(slopeHMDP, t==dataADP[z,]$t & iL== dataADP[z,]$iL & iG== dataADP[z,]$iG & iSd== dataADP[z,]$iSd)$slope )
dataADP$sHMDP <- sHMDP


#Plot 4 (slopes AVI vs VI )
library(tikzDevice)
tikz("slopesAVI.tex", width = 8.5, height = 3.5, standAlone = F)
plot(x = c(1,60), y= c(400,1000), xlab='', xaxt='n',  ylab='$\\hat{b}(t,w)$', pch=NA)
points(dataADP[35:90,]$slope,  col = "blue")
points(dataADP[35:90,]$sHMDP    , col = "green")
opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
            mar=c(0, 0, 0, 0), new=TRUE)
on.exit(par(opar))
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend("bottom", legend=c("Slopes values resulted by AVI algorithm","Slopes values resulted by VI algorithm"), pch =c(1,1), col = c("blue","green"), box.lty=0)
dev.off()
tools::texi2pdf(file = "slopesAVI.tex", clean = T)


# Plot 3 (converging trend of slopes)
dataValADP <- read.csv2("dataValADP.csv")
dataValADP$t<-as.numeric(as.character(dataValADP$t))
dataValADP$iL<-as.numeric(as.character(dataValADP$iL))
dataValADP$iG<-as.numeric(as.character(dataValADP$iG))
dataValADP$iSd<-as.numeric(as.character(dataValADP$iSd))

tikz("slopesUpdate.tex", width = 8.5, height = 5.5, standAlone = F)
plot(x = c(1,30), c(400, 1000), xlab='Updating frequency',  ylab='$\\hat{b}(t,w)$',  pch=NA)
for(i in 1:dim(dataValADP)[1] ){
  if(!is.na(dataValADP[i,]$iL) ){
    a<-i+1
    dataValADP[a,]$t<-450
    j<-i+1
    b<-0
    while( is.na(dataValADP[j,]$iL)   ){
      b<-j
      j<-j+1
      if(b>dim(dataValADP)[1]) break;
    }
  }
  if( (b-a>3) & (b-a<30)  )
    lines(y = dataValADP[a:b,]$t, x = 1:(b-a+1)) #  col=sample(rainbow(10))
  
}    
dev.off()
tools::texi2pdf(file = "slopesUpdate.tex", clean = T)

# Plot 2(slopes of VI) 
library(tikzDevice)
tikz("slopesHMDP.tex", width = 8.5, height = 5.5, standAlone = F)
plot(x = c(0,18), c(0,20000), xlab='$q$' , ylab='$\\nu(s)$', pch=NA, xaxt="n")
axis(1, las=1, at=c(0:18), labels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18") )
for(z in 36:dim(dataADP)[1] ){
  policyFilter = subset(policyExternalM,  iSW==dataADP[z,]$iL &  iSG==dataADP[z,]$iG & iSSd==dataADP[z,]$iSd  & RS==1 & phase==0 & iRation==0 & t==dataADP[z,]$t+1)
  
  #   x<-policyFilter$n
  #   y<-policyFilter$weight
  #   dd<-data.frame(
  #     x,
  #     y
  #   )
  #fit <- lm(formula=y~poly(x=x,degree=1,raw=TRUE),data=dd)
  #fit <- lm(formula=y~0+x,data=dd)
  #lines(x, predict(fit, newdata=data.frame(x=x)) )
  lines(policyFilter$n, policyFilter$weight )
}
dev.off()
tools::texi2pdf(file = "Slopes1.tex", clean = T)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------Plots for Experiment 4 (plot9)----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------


# aggregate data for all gropus of scenarios 
read_sec_dat<-function(file_path){
  numPigsSecs <- read.csv2(file_path)
  #numPigsSecs <- numPigsSecs[139:285,]
  #numPigsSecs <- numPigsSecs[178:333,]
  #numPigsSecs <- numPigsSecs[166:258,]
  #numPigsSecs <- numPigsSecs[0:150,]
  numPigsSecs$weight <- as.numeric(as.character(numPigsSecs$weight))
  numPigsSecs$growth <- as.numeric(as.character(numPigsSecs$growth))
  numPigsSecs$sd <- as.numeric(as.character(numPigsSecs$sd))
  numPigsSecs$optimal <- numPigsSecs$opt
  numPigsSecs$optimal[numPigsSecs$optimal==0] <- "C"
  numPigsSecs$optimal[numPigsSecs$optimal==numPigsSecs$numPigs]<-"T"
  numPigsSecs$secLabel <- numPigsSecs$section
  numPigsSecs$secLabel[numPigsSecs$secLabel==0]<-"Section 1"
  numPigsSecs$secLabel[numPigsSecs$secLabel==1]<-"Section 2"
  numPigsSecs$secLabel[numPigsSecs$secLabel==2]<-"Section 3"
  return(numPigsSecs)
}
basic<-read_sec_dat(file_path = "last_section/basic/numPigsSecs.csv")
feed_h <- read_sec_dat(file_path = "last_section/feed_20+/numPigsSecs.csv")
feed_l<-read_sec_dat(file_path = "last_section/feed_20-/numPigsSecs.csv")
pork_h<-read_sec_dat(file_path = "last_section/pork_20+/numPigsSecs.csv")
pork_l<-read_sec_dat(file_path = "last_section/pork_20-/numPigsSecs.csv")
transport_l<-read_sec_dat(file_path = "last_section/transport_102/numPigsSecs.csv")
transport_h<-read_sec_dat(file_path = "last_section/transport_305/numPigsSecs.csv")
main_dat_cull<-basic %>% dplyr::mutate(totalCull_feed_l=feed_l$totalCull,
                                  totalCull_feed_h=feed_h$totalCull,
                                  totalCull_pork_l=pork_l$totalCull,
                                  totalCull_pork_h=pork_h$totalCull,
                                  totalCull_transport_l=transport_l$totalCull,
                                  totalCull_transport_h=transport_h$totalCull
                                  ) %>% dplyr::select(epoch,totalCull,totalCull_feed_l,
                                                      totalCull_feed_h,
                                                      totalCull_pork_l,
                                                      totalCull_pork_h,
                                                      totalCull_transport_l,
                                                      totalCull_transport_h) %>% distinct()
main_dat_truck<-basic %>% dplyr::mutate(truckNum_feed_l=feed_l$truckNum,
                                  truckNum_feed_h=feed_h$truckNum,
                                  truckNum_pork_l=pork_l$truckNum,
                                  truckNum_pork_h=pork_h$truckNum,
                                  truckNum_transport_l=transport_l$truckNum,
                                  truckNum_transport_h=transport_h$truckNum
                    ) %>% dplyr::select(epoch,truckNum,truckNum_feed_l,
                                       truckNum_feed_h,
                                       truckNum_pork_l,
                                       truckNum_pork_h,
                                      truckNum_transport_l,
                                      truckNum_transport_h) %>% distinct()
main_dat_cull<-main_dat_cull[81:120,]
main_dat_truck<-main_dat_truck[81:120,]
dat_cull<-melt(data.frame(main_dat_cull),
               # ID variables - all the variables to keep but not split apart on
               id.vars="epoch",
               # The source columns
               measure.vars=c("totalCull","totalCull_feed_l","totalCull_feed_h","totalCull_pork_l","totalCull_pork_h","totalCull_transport_l","totalCull_transport_h"),
               # Name of the destination column that will identify the original
               # column that the measurement came from
               variable.name="name",
               value.name="y"
)
dat_truck<-melt(data.frame(main_dat_truck),
               # ID variables - all the variables to keep but not split apart on
               id.vars="epoch",
               # The source columns
               measure.vars=c("truckNum","truckNum_feed_l","truckNum_feed_h","truckNum_pork_l","truckNum_pork_h","truckNum_transport_l","truckNum_transport_h"),
               # Name of the destination column that will identify the original
               # column that the measurement came from
               variable.name="name",
               value.name="y"
)
dat<- dat_cull %>% dplyr::mutate(numTrucks=dat_truck$y)
dat$class<-NA
dat[dat$name %in%c('totalCull'),]$class<-'Basic'
dat[dat$name %in%c('totalCull_feed_l','totalCull_feed_h'),]$class<-'Group 2'
dat[dat$name %in%c('totalCull_pork_l','totalCull_pork_h'),]$class<-'Group 3'
dat[dat$name %in%c('totalCull_transport_l','totalCull_transport_h'),]$class<-'Group 1'
dat$name<-mapvalues(dat$name, from = c("totalCull","totalCull_feed_l","totalCull_feed_h","totalCull_pork_l","totalCull_pork_h","totalCull_transport_l","totalCull_transport_h"), 
                         to = c("Basic",
                                "$p^{\\mathtt{\\small feed}}\\downarrow$",
                                "$p^{\\mathtt{\\small feed}}\\uparrow$",
                                "$p^{\\mathtt{\\small pork}}\\downarrow$",
                                "$p^{\\mathtt{\\small prork}}\\uparrow$",
                                "$k^{\\mathtt{\\small truck}}\\downarrow$",
                                "$k^{\\mathtt{\\small truck}}\\uparrow$") )

#plot 9
tikz("sensitivity_cull.tex", width = 10, height = 5.5, standAlone = F)
plot<-ggplot(data=dat, aes(x = epoch, y = y, group = name, fill =name ) ) + 
  scale_y_continuous(breaks=seq(0,1100,100), label = seq(0,1100,100) ) +
  scale_x_continuous(breaks=seq(min(dat$epoch),max(dat$epoch),1), labels=(min(dat$epoch)+1):(max(dat$epoch)+1) ) +
  geom_bar(stat="identity", alpha = 1/3, position=position_dodge(), width=0.7, colour="black") + 
  #scale_fill_grey(start=0) +
  scale_fill_manual(
    values=c("gray0","khaki4", "linen","grey25","grey85", "gray95", "gray100" ) 
    ) + 
  xlab("Week number") + ylab(" ") +
  facet_grid(class ~ .)
plot  + guides(fill=guide_legend(nrow=1,byrow=TRUE)) + theme_bw() + 
  theme( panel.background = element_blank(), 
         panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
         legend.position="bottom",
         legend.direction="horizontal",
         legend.title = element_blank(),
         axis.title.x= element_text(vjust = -0.7, size = 7),
         axis.text=element_text(size=6),
         strip.background=element_rect(fill = NA)
         ) 
dev.off()
tools::texi2pdf(file = "sensitivity_cull.tex", clean = T)





























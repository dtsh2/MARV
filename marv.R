rm(list=ls())
## estimate parameters, inc R0
seroneg<-(1622-(40+250))/(1622)
R0<-1/seroneg
mu<-0.000510492
infprev<-40/1622
betaM<-mu/infprev*(R0-1)

## MARV SEIR model with density-dept transmission
## 25 years Sens Analysis

# instal r developer toolbox first (Rtools from Cran R)

## pomp test run lbv
getwd()
setwd("~/GitHub/MARV") # revise as necessary
# setwd("C:/Users/dtshayma/Dropbox")
library(pomp)

#Compiling C code and loading the dll
# dyn.unload("marburgmodel.dll")

#system("R CMD SHLIB marburgmodel.c")

dyn.load("marburgmodel.dll")

pomp(
  data = data.frame(
    time=seq(from=0,to=365*25,by=1),  # time for simulations to run
    X = NA # dummy variables
    ),
  times="time",
  t0=0,
  ## native routine for the process simulator:
  rprocess=euler.sim(
    step.fun="sir_euler_simulator",
    delta.t=1,
    #PACKAGE="pomp"  ## may be include if does not work - this is where to look for the c file 
  ## name of the shared-object library containing the PACKAGE="pomp",
  ),
  
  ## the order of the state variables assumed in the native routines:
  statenames=c("SUSJ","EXPJ","INFJ", "RECJ", "SUSA", "EXPA","INFA", "RECA","PREVA","PREVJ","SPA","SPJ"),
  ## the order of the parameters assumed in the native routines:
  paramnames=c("BETA","MU","DELTA","SIGMA","K","EPSILON","TAU","KAPPA","S","OMEGA","PHI",
               "SUSJ.0","EXPJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EXPA.0","INFA.0", "RECA.0","PREVA.0","PREVJ.0","SPA.0","SPJ.0"),
  initializer=function(params,t0,statenames,...){
    x0<-params[paste(statenames,".0",sep="")]
    names(x0)<-statenames
    return(x0)
  }
) -> seir

params <- c(
  BETA=0.004506855,
  MU=0.000510492,
  DELTA=0.002312247,
  SIGMA=1/21,
  K=40000,
  EPSILON=1/365,
  TAU=1/7,
  KAPPA=1.5/365,
  S=14.35,
  OMEGA=2/365,
  PHI=0.0,
  SUSJ.0=28000,EXPJ.0=1000,INFJ.0=1000,RECJ.0=10000,
  SUSA.0=28000, EXPA.0=1000,INFA.0=1000, RECA.0=10000,
  PREVA.0= 0.025, PREVJ.0=0.025 ,
  SPA.0=0.25,SPJ.0=0.25) # this adds to the initial conditions given the state variables

sim <- simulate(seir,params=c(params),#seed=3493885L,
                nsim=100,states=T,obs=F,as.data.frame=T) # 
class(seir) # pomp object
class(sim) # data frame - even if I remove "as.data.frame" in the above code (sim)
#sim <- simulate(sir,params=c(params),nsim=1,states=T,obs=F)#,as.data.frame=T) # saves as an array
# pf<-pfilter(sim,params=c(params),Np=1000) # won't work, because this is a data frame, not pomp object

inf.res.a<-matrix(sim$INFA,ncol=100)
inf.res.j<-matrix(sim$INFJ,ncol=100)
susj.res.j<-matrix(sim$SUSJ,ncol=100)

par(mfrow=c(1,1))
#matplot(inf.res.a,col="grey",pch=1,ylab="Numbers",xlab="Time",xaxt="n")
#lines(rowMeans(inf.res.a))
#matplot(inf.res.j,col="grey",pch=1,ylab="Numbers",xlab="Time",xaxt="n")
#lines(rowMeans(inf.res.j))

matplot(susj.res.j[8000:9000,],col="darkgrey",type="l",ylab="Numbers",xlab="Time",xaxt="n")
#lines(rowMeans(inf.res.j[8000:9000,]))
susj.res.j[susj.res.j == 0] <- NA
lines(rowMeans(susj.res.j[8000:9000,],na.rm=T),col="white")

matplot(inf.res.a[8000:9000,],col="lightgrey",type="l",#ylab="Numbers",xlab="Time",
        xaxt="n",add=T)
#lines(rowMeans(inf.res.a[8000:9000,]))
inf.res.a[inf.res.a == 0] <- NA
lines(rowMeans(inf.res.a[8000:9000,],na.rm=T))#,col="white")

matplot(inf.res.j[8000:9000,],col="black",type="l",#ylab="Numbers",xlab="Time",
        xaxt="n",add=T)
#lines(rowMeans(inf.res.a[8000:9000,]))
inf.res.j[inf.res.j == 0] <- NA
lines(rowMeans(inf.res.j[8000:9000,],na.rm=T),col="white")

legend("topright", c("Juvenile","Adult"), fill = c("grey","black"), col = c("black","white"),
       border = c("grey","black"), bty = "n")

#######################

inf.res.a[is.na(inf.res.a)] <- 0
inf.res.a[is.na(inf.res.a)] <- 0

resbp<-array(NA,dim=c(100))

for (i in 1:100){ # each stochastic run
  resbp[i]<-ifelse(sum(inf.res.a[9126,i],inf.res.j[9126,i])>0,1,0) # extinct for each run
}
mean(resbp)
#text(500, 200, "the text is CENTERED around (x,y) = (X,X) by default",cex = .8)

#########################################################
## carrying capacity

nonVarying = matrix(c(
  BETA=0.004506855,
  MU=0.000510492,
  DELTA=0.002312247,
  SIGMA=1/21,
  #K=40000,
  EPSILON=1/365,
  TAU=1/7,
  KAPPA=1.5/365,
  S=14.35,
  OMEGA=2/365,
  PHI=0.0,
  SUSJ.0=37000,
  EXPJ.0=1000,
  INFJ.0=1000,
  RECJ.0=1000,
  SUSA.0=37000,
  EXPA.0=1000,
  INFA.0=1000,
  RECA.0=1000,
  PREVA.0= 0.025,
  PREVJ.0=0.025 ,
  SPA.0=0.025,
  SPJ.0=0.025),
  ncol=22,
  nrow=40,
  byrow=T) #binded with non-varying parameters


dimnames(nonVarying)[[2]]=c("BETA","MU","DELTA","SIGMA",#"K",
                            "EPSILON","TAU","KAPPA","S","OMEGA","PHI",
                            "SUSJ.0","EXPJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EXPA.0","INFA.0", "RECA.0","PREVA.0","PREVJ.0","SPA.0","SPJ.0") # naming non-varying columns

## from other code
Kset=seq(from = 100, to=200000, by =5000)

fullParamSets = cbind(nonVarying,Kset) # full parameter set
#fullParamSets[,1] <- fullParamSets[,1]*fullParamSets[,29]
head(fullParamSets)
dim(fullParamSets)

dimnames(fullParamSets)[[2]]=c("BETA","MU","DELTA","SIGMA",#"K",
                               "EPSILON","TAU","KAPPA","S","OMEGA","PHI",
                               "SUSJ.0","EXPJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EXPA.0","INFA.0", "RECA.0","PREVA.0","PREVJ.0","SPA.0","SPJ.0","K")

# order for pomp/C model:  
BETA = fullParamSets[,1]
MU = fullParamSets[,2]
DELTA = fullParamSets[,3]
SIGMA = fullParamSets[,4]
EPSILON = fullParamSets[,5]
TAU = fullParamSets[,6]
KAPPA = fullParamSets[,7]
S = fullParamSets[,8]
OMEGA = fullParamSets[,9]
PHI = fullParamSets[,10]
SUSJ.0 = fullParamSets[,11]
EXPJ.0 = fullParamSets[,12]
INFJ.0 = fullParamSets[,13]
RECJ.0 = fullParamSets[,14]
SUSA.0 = fullParamSets[,15]
EXPA.0 = fullParamSets[,16]
INFA.0 = fullParamSets[,17]
RECA.0 = fullParamSets[,18]
PREVA.0 = fullParamSets[,19]
PREVJ.0 = fullParamSets[,20]
SPA.0 = fullParamSets[,21]
SPJ.0 = fullParamSets[,22]
K = fullParamSets[,23]

paramset<-cbind(BETA,MU,DELTA,SIGMA,K,EPSILON,TAU,KAPPA,S,OMEGA,PHI,
                SUSJ.0,EXPJ.0,INFJ.0, RECJ.0, SUSA.0,EXPA.0,INFA.0,RECA.0,PREVA.0,PREVJ.0,SPA.0,SPJ.0)
######################################################################################

## Calling requisite libraries for parallel computing

#library(foreach)
#library(doSNOW)

#Setting up "parallel backend"

#w<-makeCluster(3,type="SOCK") # makes the cluster, i.e. no of cores ABC = 8 cores, DEF = 12 see performance to see # of cores
#registerDoSNOW(w) # 

#Checks that the number of workers is set up correctly.

#getDoParWorkers()

#######################################################3
results<-array(NA,dim=c(40,1,5))

# for one parameter set....
out1 <-simulate(seir,params=c(paramset[1,]),
                seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #

outres1 <- out1[seq(from=9126,to=912600,by=9126),] # select last #s
N1 = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
for (i in 1:100){ # each stochastic run
  N1[i,1]<-sum(outres1[i,1:8]) # pop pers
  N1[i,2]<-(sum(outres1[i,3],outres1[i,7])/sum(outres1[i,1:8]))*100 # prevalence; total
  N1[i,3]<-((outres1[i,8])/(sum(outres1[i,5:8])))*100 # adult seroprevalence; total
  N1[i,4]<-ifelse(sum(outres1[i,1:8])>0,1,0) # population extinct for each run
  N1[i,5]<-ifelse(sum(outres1[i,3],outres1[i,7])>0,1,0) # pathogen extinction for each run
}
N1[is.na(N1)]<- 0
## now average
M1 = array(0,c(1,5))
M1[1] = mean(N1[1:100,1]) # population size
M1[2] = mean(N1[1:100,2]) # prevalence
M1[3] = mean(N1[1:100,3]) # adult seroprevalence
M1[4] = mean(N1[1:100,4]) # adult seroprevalence
M1[5] = mean(N1[1:100,5]) # adult seroprevalence
rm(out1)
M1
results[1,,]<-M1
results[1,,]

##########################################################33
# for all parameter sets....
results<-array(NA,dim=c(40,1,6))

for (j in 1:length(paramset[,1])){
  out <-simulate(seir,params=c(paramset[j,]),
                 seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #
  outres <- out[seq(from=9126,to=912600,by=9126),] # select last #s
  N = array(0,c(100,6)) # same dimensions as No. runs * outputs I want
  for (i in 1:100){ # each stochastic run
    N[i,1]<-sum(outres[i,1:8])
    N[i,2]<-(sum(outres[i,3],outres[i,7])/sum(outres[i,1:8]))*100 # prevalence; total
    N[i,3]<-((outres[i,8])/(sum(outres[i,5:8])))*100 # adult seroprevalence; total
    N[i,4]<-((outres[i,4])/(sum(outres[i,1:4])))*100 # juvenile seroprevalence; total
    N[i,5]<-ifelse(sum(outres[i,1:8])>0,1,0) # population extinct for each run
    N[i,6]<-ifelse(sum(outres[i,3],outres[i,7])>0,1,0) # pathogen extinction for each run
  }
  N[is.na(N)]<- 0
  ## now average
  M = array(0,c(1,5))
  M[1] = mean(N[1:100,1]) # population size
  M[2] = mean(N[1:100,2]) # prevalence
  M[3] = mean(N[1:100,3]) # adult seroprevalence
  M[4] = mean(N[1:100,4]) # juvenile seroprevalence
  M[5] = mean(N[1:100,5]) # mean pop extinction
  M[6] = mean(N[1:100,6]) # mean path extinction
  rm(out)
  results[j,,]<-M
}
#
#########################################################33

# w<-makeCluster(1,type="SOCK") # return to one core

############################################################3
## need matrix of results...
X<-aperm(results,c(1,2,3))
dim(X)<-c(40,6)
head(X)
tail(X)

################################################################################
# below for K
#
#plot(X[,1],X[,2])
# plot
#

par(omi=c(1,1.5,0.5,0.5))
par(mai=c(1.5,1.5,0.5,0.5))
par(mar=c(4, 4, 4, 4) + 0.1)
par(mfrow=c(1,3))
plot(X[1:20,1],X[1:20,2],pch=16,
     ylab="Mean prevalence (%)",xlab="",
     col="grey25", cex.lab=1)
abline(v=40000,lty=2)
mtext('a',font=2,side=3,line=1,at=-2,cex=1)
plot(X[1:20,1],X[1:20,3],pch=16,
     ylab="Mean seroprevalence (%)",xlab="",
     col="grey25", cex.lab=1)
points(X[1:20,1],X[1:20,4],pch=1,
     col="grey25")
abline(v=40000,lty=2)
#legend("right",c("adult","juvenile"),pch=c(16,1),bty="n")
mtext('b',font=2,side=3,line=1,at=-2,cex=1)
plot(X[1:20,1],X[1:20,6],pch=16,
     ylab="P[persist]",xlab="",
     col="grey25", cex.lab=1)
abline(v=40000,lty=2)
mtext("Population size",side=1,outer=T,cex=0.9)
mtext('c',font=2,side=3,line=1,at=-2,cex=1)

#plot(paramset[,5],X[,5],pch=16)
##########################################################################
## single birth pulse vs 2
## 7 vs 21 day inf period

nonVarying = matrix(c(
  BETA=0.004506855,
  MU=0.000510492,
  DELTA=0.002312247,
  #SIGMA=1/21,
  #K=40000,
  EPSILON=1/365,
  TAU=1/7,
  KAPPA=1.5/365,
  S=14.35,
  #OMEGA=1/365,
  PHI=0.0,
  SUSJ.0=37000,
  EXPJ.0=1000,
  INFJ.0=1000,
  RECJ.0=1000,
  SUSA.0=37000,
  EXPA.0=1000,
  INFA.0=1000,
  RECA.0=1000,
  PREVA.0= 0.025,
  PREVJ.0=0.025 ,
  SPA.0=0.025,
  SPJ.0=0.025),
  ncol=20,
  nrow=40,
  byrow=T) #binded with non-varying parameters


dimnames(nonVarying)[[2]]=c("BETA","MU","DELTA",#"SIGMA",#"K",
                            "EPSILON",
                            "TAU","KAPPA","S",#"OMEGA",
                            "PHI",
                            "SUSJ.0","EXPJ.0","INFJ.0", "RECJ.0", "SUSA.0", "EXPA.0","INFA.0", "RECA.0","PREVA.0","PREVJ.0","SPA.0","SPJ.0") # naming non-varying columns

## from other code

sigmaset4=rep(c(rep(1/7,40),rep(1/21,40)),2)
Kset=seq(from = 100, to=200000, by =5000)
omegaset4=c(rep(1/365,80),rep(2/365,80))
Kset4=rep(Kset,4)
nonVar4=rbind(nonVarying,nonVarying,nonVarying,nonVarying)

fullParamSets = cbind(nonVar4,Kset4,sigmaset4,omegaset4) # full parameter set
#fullParamSets[,1] <- fullParamSets[,1]*fullParamSets[,29]
head(fullParamSets)
dim(fullParamSets)

dimnames(fullParamSets)[[2]]=c("BETA", #1
                               "MU", #2
                               "DELTA", #3
                               "EPSILON",
                               "TAU", #
                               "KAPPA", #
                               "S", #
                               #"OMEGA", #
                               "PHI", #
                               "SUSJ.0", #
                               "EXPJ.0", #
                               "INFJ.0", #
                               "RECJ.0", #
                               "SUSA.0", #
                               "EXPA.0", #
                               "INFA.0", #
                               "RECA.0", #
                               "PREVA.0", #
                               "PREVJ.0", #
                               "SPA.0", #
                               "SPJ.0", #
                               "K", #
                               "SIGMA", #22
                               "OMEGA") #23

# order for pomp/C model:  
BETA = fullParamSets[,1]
MU = fullParamSets[,2]
DELTA = fullParamSets[,3]
SIGMA = fullParamSets[,22]
EPSILON = fullParamSets[,4]
TAU = fullParamSets[,5]
KAPPA = fullParamSets[,6]
S = fullParamSets[,7]
OMEGA = fullParamSets[,23]
PHI = fullParamSets[,8]
SUSJ.0 = fullParamSets[,9]
EXPJ.0 = fullParamSets[,10]
INFJ.0 = fullParamSets[,11]
RECJ.0 = fullParamSets[,12]
SUSA.0 = fullParamSets[,13]
EXPA.0 = fullParamSets[,14]
INFA.0 = fullParamSets[,15]
RECA.0 = fullParamSets[,16]
PREVA.0 = fullParamSets[,17]
PREVJ.0 = fullParamSets[,18]
SPA.0 = fullParamSets[,19]
SPJ.0 = fullParamSets[,20]
K = fullParamSets[,21]

paramset<-cbind(BETA,MU,DELTA,SIGMA,K,EPSILON,TAU,KAPPA,S,OMEGA,PHI,
                SUSJ.0,EXPJ.0,INFJ.0, RECJ.0, SUSA.0,EXPA.0,INFA.0,RECA.0,PREVA.0,PREVJ.0,SPA.0,SPJ.0)
######################################################################################

## Calling requisite libraries for parallel computing

#library(foreach)
#library(doSNOW)

#Setting up "parallel backend"

#w<-makeCluster(3,type="SOCK") # makes the cluster, i.e. no of cores ABC = 8 cores, DEF = 12 see performance to see # of cores
#registerDoSNOW(w) # 

#Checks that the number of workers is set up correctly.

#getDoParWorkers()

#######################################################3
results<-array(NA,dim=c(160,1,5))

# for one parameter set....
out1 <-simulate(seir,params=c(paramset[1,]),
                seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #

outres1 <- out1[seq(from=9126,to=912600,by=9126),] # select last #s
N1 = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
for (i in 1:100){ # each stochastic run
  N1[i,1]<-sum(outres1[i,1:8]) # pop pers
  N1[i,2]<-(sum(outres1[i,3],outres1[i,7])/sum(outres1[i,1:8]))*100 # prevalence; total
  N1[i,3]<-((outres1[i,8])/(sum(outres1[i,5:8])))*100 # adult seroprevalence; total
  N1[i,4]<-ifelse(sum(outres1[i,1:8])>0,1,0) # population extinct for each run
  N1[i,5]<-ifelse(sum(outres1[i,3],outres1[i,7])>0,1,0) # pathogen extinction for each run
}
N1[is.na(N1)]<- 0
## now average
M1 = array(0,c(1,5))
M1[1] = mean(N1[1:100,1]) # population size
M1[2] = mean(N1[1:100,2]) # prevalence
M1[3] = mean(N1[1:100,3]) # adult seroprevalence
M1[4] = mean(N1[1:100,4]) # adult seroprevalence
M1[5] = mean(N1[1:100,5]) # adult seroprevalence
rm(out1)
M1
results[1,,]<-M1
results[1,,]

##########################################################33
# for all parameter sets....
results<-array(NA,dim=c(160,1,5))

for (j in 1:length(paramset[,1])){
  out <-simulate(seir,params=c(paramset[j,]),
                 seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #
  outres <- out[seq(from=9126,to=912600,by=9126),] # select last #s
  N = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
  for (i in 1:100){ # each stochastic run
    N[i,1]<-sum(outres[i,1:8])
    N[i,2]<-(sum(outres[i,3],outres[i,7])/sum(outres[i,1:8]))*100 # prevalence; total
    N[i,3]<-((outres[i,7])/(sum(outres[i,5:8])))*100 # adult seroprevalence; total
    N[i,4]<-ifelse(sum(outres[i,1:8])>0,1,0) # population extinct for each run
    N[i,5]<-ifelse(sum(outres[i,3],outres[i,7])>0,1,0) # pathogen extinction for each run
  }
  N[is.na(N)]<- 0
  ## now average
  M = array(0,c(1,5))
  M[1] = mean(N[1:100,1]) # population size
  M[2] = mean(N[1:100,2]) # prevalence
  M[3] = mean(N[1:100,3]) # adult seroprevalence
  M[4] = mean(N[1:100,4]) # mean pop extinction
  M[5] = mean(N[1:100,5]) # mean path extinction
  rm(out)
  results[j,,]<-M
}
#
#########################################################33

#w<-makeCluster(1,type="SOCK") # return to one core

############################################################3
## need matrix of results...
X<-aperm(results,c(1,2,3))
dim(X)<-c(160,5)
head(X)
tail(X)

par(omi=c(1,1,0.5,0.5))
par(mai=c(0.8,0.8,0.8,0.8))

par(mfrow=c(2,2))
plot(X[1:40,1],X[1:40,5],pch=16,ylim=c(0,1),
     #ylab="P[persist]",xlab="Population size",
     col="grey25", cex.lab=1.2,ann=F)
text(150000,0.8,expression(paste("birth pulse/yr"==1)))
text(150000,0.6,expression(paste(1/sigma==7)))
mtext('a',font=2,side=3,line=1,at=-2,cex=1.2)

plot(X[41:80,1],X[41:80,5],pch=16,ylim=c(0,1),
     #ylab="P[persist]",xlab="Population size",
     col="grey25", cex.lab=1.2,ann=F)
text(150000,0.8,expression(paste("birth pulse/yr"==1)))
text(150000,0.6,expression(paste(1/sigma==21)))
mtext('b',font=2,side=3,line=1,at=-2,cex=1.2)

plot(X[81:120,1],X[81:120,5],pch=16,ylim=c(0,1),
     #ylab="P[persist]",xlab="Population size",
     col="grey25", cex.lab=1.2,ann=F)
text(150000,0.8,expression(paste("birth pulse/yr"==2)))
text(150000,0.6,expression(paste(1/sigma==7)))
mtext('c',font=2,side=3,line=1,at=-2,cex=1.2)

plot(X[121:160,1],X[121:160,5],pch=16,ylim=c(0,1),
     #ylab="P[persist]",xlab="Population size",
     col="grey25", cex.lab=1.2,ann=F)
text(150000,0.8,expression(paste("birth pulse/yr"==2)))
text(150000,0.6,expression(paste(1/sigma==21)))
mtext('d',font=2,side=3,line=1,at=-2,cex=1.2)

mtext("P (persist)",side=2,outer=T)
mtext("Population size",side=1,outer=T)

######################### rest of this code requires the MLE estimate for the rho and beta params
## for LHS parameter set
## Calling requisite libraries for parallel computing

#library(foreach)
#library(doSNOW)

#Setting up "parallel backend"

#w<-makeCluster(3,type="SOCK") # makes the cluster, i.e. no of cores ABC = 8 cores, DEF = 12 see performance to see # of cores
#registerDoSNOW(w) # 

#Checks that the number of workers is set up correctly.

#getDoParWorkers()

## ~~~~~~~~~~ LHS SAMPLING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#############
# NB not testing K - carrying cap, or rate of aging..

library(lhs)

nspaces=100 ## how many bins/intervals, 100 needed for PRCC

hypercube=randomLHS(n=nspaces, k=9) ## function with N columns
dimnames(hypercube)[[2]]=c("beta","mu","delta",
                           "sigma","K",#"epsilon",
                           "tau","k",
                           "s",#"omega",
                           "phi")  # named columns with parameter names


mins = c(   beta=0.0004506855,
            mu=0.0000510492,
            delta=0.0002312247,
            sigma=1/210,
            K=10000,
          #  epsilon=1/365,
            tau=1/70,
            kappa=0.15/365,
            s=1.435,
          #  omega=1/365,
            phi=0.0)

maxs = c(  beta=0.04506855,
           mu=0.00510492,
           delta=0.02312247,
           sigma=1/2.1,
           K=1000000,
          # epsilon=1/365,
           tau=1,
           kappa=15/365,
           s=140.35,
          # omega=2/365,
           phi=1)                      

diffs=maxs-mins ## range of each variable

hypercubeadj = hypercube # create copy of hypercube samples to modify, hypercube adjusted; i.e. new matrix
for (i in 1:ncol(hypercube)){
  hypercubeadj[,i]=hypercube[,i]*diffs[i]+mins[i] # scale samples to difference and add minimum
}

head(hypercubeadj)
dim(hypercubeadj)

dimnames(hypercubeadj)[[2]]=c("beta",
                              "mu",
                              "delta",
                              "sigma",
                              "K",
                              # epsilon=1/365,
                              "tau",
                              "kappa",
                              "s",
                              #"omega",
                              "phi")

nonVarying = matrix(c(
  #K= 1000000,			# K, however, population size fluctuates...up to >1*10^6
  omega=1/365, # # pulses
  epsilon= 2/365,	# rate of aging for those juveniles, should be ~ annual - per cap per year
  SUSJ.0=37000,EXPJ.0=1000,INFJ.0=1000,RECJ.0=1000,
  SUSA.0=37000, EXPA.0=1000,INFA.0=1000, RECA.0=1000,
  PREVA.0= 0.025, PREVJ.0=0.025 ,
  SPA.0=0.025,SPJ.0=0.025),
              ncol=14,
              nrow=100,
              byrow=T) #binded with non-varying parameters


dimnames(nonVarying)[[2]]=c("epsilon",  # rate of aging for those juveniles, should be ~ annual - per cap per year
                            "omega",
                            "SUSJ.0","EXPJ.0","INFJ.0","RECJ.0",
                            "SUSA.0", "EXPA.0","INFA.0", "RECA.0",
                            "PREVA.0", "PREVJ.0",
                            "SPA.0","SPJ.0") # naming non-varying columns

fullParamSets = cbind(hypercubeadj,nonVarying) # full parameter set

head(fullParamSets)
dim(fullParamSets)

# order for pomp/C model:  

BETA = fullParamSets[,1]
MU = fullParamSets[,2]
DELTA = fullParamSets[,3]
SIGMA = fullParamSets[,4]
K = fullParamSets[,5]
EPSILON = fullParamSets[,10]
TAU = fullParamSets[,6]
KAPPA = fullParamSets[,7]
S = fullParamSets[,8]
OMEGA = fullParamSets[,11]
PHI = fullParamSets[,9]
SUSJ.0 = fullParamSets[,12]
EXPJ.0 = fullParamSets[,13]
INFJ.0 = fullParamSets[,14]
RECJ.0 = fullParamSets[,15]
SUSA.0 = fullParamSets[,16]
EXPA.0 = fullParamSets[,17]
INFA.0 = fullParamSets[,18]
RECA.0 = fullParamSets[,19]
PREVA.0 = fullParamSets[,20]
PREVJ.0 = fullParamSets[,21]
SPA.0 =  fullParamSets[,22]
SPJ.0 =  fullParamSets[,23]

paramset<-cbind(BETA,
                MU,
                DELTA,
                SIGMA,
                K,
                EPSILON,
                TAU,
                KAPPA,
                S,
                OMEGA,
                PHI,
                SUSJ.0,
                EXPJ.0,
                INFJ.0,
                RECJ.0,
                SUSA.0,
                EXPA.0,
                INFA.0,
                RECA.0,
                PREVA.0,
                PREVJ.0,
                SPA.0,
                SPJ.0)

######################################################################################
results<-array(NA,dim=c(100,1,5))


# for one parameter set....
out1 <-simulate(seir,params=c(paramset[1,]),
                seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #

outres1 <- out1[seq(from=9126,to=912600,by=9126),] # select last #s
N1 = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
for (i in 1:100){ # each stochastic run
  N1[i,1]<-sum(outres1[i,1:8]) # pop pers
  N1[i,2]<-(sum(outres1[i,3],outres1[i,7])/sum(outres1[i,1:8]))*100 # prevalence; total
  N1[i,3]<-((outres1[i,8])/(sum(outres1[i,5:8])))*100 # adult seroprevalence; total
  N1[i,4]<-ifelse(sum(outres1[i,1:8])>0,1,0) # population extinct for each run
  N1[i,5]<-ifelse(sum(outres1[i,3],outres1[i,7])>0,1,0) # pathogen extinction for each run
}
N1[is.na(N1)]<- 0
## now average
M1 = array(0,c(1,5))
M1[1] = mean(N1[1:100,1]) # population size
M1[2] = mean(N1[1:100,2]) # prevalence
M1[3] = mean(N1[1:100,3]) # adult seroprevalence
M1[4] = mean(N1[1:100,4]) # adult seroprevalence
M1[5] = mean(N1[1:100,5]) # adult seroprevalence
rm(out1)
M1
results[1,,]<-M1
results[1,,]

##########################################################
##
##
##  check # of sims and dimensions
##
##
# for all parameter sets....
results<-array(NA,dim=c(100,1,5))

for (j in 1:length(paramset[,1])){
  out <-simulate(seir,params=c(paramset[j,]),
                 seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #
  outres <- out[seq(from=9126,to=912600,by=9126),] # select last #s
  N = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
  for (i in 1:100){ # each stochastic run
    N[i,1]<-sum(outres[i,1:8])
    N[i,2]<-(sum(outres[i,3],outres[i,7])/sum(outres[i,1:8]))*100 # prevalence; total
    N[i,3]<-((outres[i,7])/(sum(outres[i,5:8])))*100 # adult seroprevalence; total
    N[i,4]<-ifelse(sum(outres[i,1:8])>0,1,0) # population extinct for each run
    N[i,5]<-ifelse(sum(outres[i,3],outres[i,7])>0,1,0) # pathogen extinction for each run
  }
  N[is.na(N)]<- 0
  ## now average
  M = array(0,c(1,5))
  M[1] = mean(N[1:100,1]) # population size
  M[2] = mean(N[1:100,2]) # prevalence
  M[3] = mean(N[1:100,3]) # adult seroprevalence
  M[4] = mean(N[1:100,4]) # mean pop extinction
  M[5] = mean(N[1:100,5]) # mean path extinction
  rm(out)
  results[j,,]<-M
}
#
#########################################################33

w<-makeCluster(1,type="SOCK") # return to one core

  ############################################################3
## need matrix of results...
X<-aperm(results,c(1,2,3))
dim(X)<-c(100,5)
head(X)
tail(X)

################################################################################
#####Functions to calculate and plot partial-rank correlation coefficients
#####(PRCCs) between parameter values and model output.
#####
#####Written by: Michael Buhnerkempe
#####      Date: Oct. 7, 2011
#####
##### Functions:
#####          prcc - calculate the partial rank correlation coefficient
#####     plot.prcc - plot a prcc object
#####
##### A brief example is presented at the end
################################################################################

################################################################################
## prcc - function to calculate partial rank correlation coefficients between
##           each of p parameters and k model outputs using n different
##           observations (number of parameter sets)
##
##    Arguments:
##          par.mat = n x p matrix containing the parameter values obtained
##                        from Latin Hypercube Sampling
##      model.output = n x k matrix containing the model outputs
##           routine = how should the PRCCs be calculated? One of:
##                        "blower" - calculated according to Appendix A in
##                                   Blower and Dowlatabadi (1994). DEFAULT.
##                        "regression" - calculated using a regression approach.
##                                   Here, the partial correlation coefficient
##                                   is defined as the correlation between the
##                                   residuals after regressing a model output
##                                   on all of the parameters except for the
##                                   parameter of interest and the residuals
##                                   after regressing the parameter of interest
##                                   on all of the other parameters. This can
##                                   be interpreted as the correlation between
##                                   of a parameter and the model output when
##                                   the effects of all other parameters have
##                                   been removed.
##         par.names = names of parameters
##      output.names = names of model outputs
##               ... = additional arguments to be passed to functions called
##                     within this function
##
##
##    Output attributes:
##      $par.matrix = original matrix of parameter set
##      $model.output = original model output
##      $(model output name) = prcc results for the named model output

prcc = function( par.mat, model.output, routine = "blower",
                 par.names = NA,output.names = NA, ...){
  
  #Make sure par.mat and model.output are matrices
  par.mat = as.matrix(par.mat)
  model.output = as.matrix(model.output)
  
  #How many parameter sets are there?
  n = length(par.mat[,1])
  
  #How many parameters are there?
  p = length(par.mat[1,])
  
  #How many model outputs are we calculating PRCCs for?
  k = length(model.output[1,])
  
  #Find the ranks of the parameter values and the model output
  par.rank = apply(par.mat,2,rank,...)
  output.rank = apply(model.output,2,rank,...)
  
  #What is the average rank?
  ave.rank = (1 + n)/2
  
  #Create a list object to store the PRCCs
  results = list()
  
  results$num.sets = n
  
  #Try to automatically get parameter and output names if they are not
  #given
  if( sum(is.na(par.names)) > 0){par.names=dimnames(par.mat)[[2]]}
  if( sum(is.na(output.names)) > 0){output.names=dimnames(model.output)[[2]]}
  
  ########################################################################
  #Calculate the PRCCs using Appendix A from Blower and Dowlatabadi (1994)
  ########################################################################
  if( routine == "blower" ){
    
    #Do the calculation for each model output
    for( i in 1:k ){
      
      work.mat = cbind(par.rank,output.rank[,i])
      
      C.temp = matrix(0,nrow=p+1,ncol=p+1)
      
      #Calculate the C matrix
      for( j in 1:(p+1) ){
        for( l in 1:(p+1) ){
          
          C.temp[j,l]=sum((work.mat[,j]-ave.rank)*(work.mat[,l]-ave.rank))/
            sqrt(sum((work.mat[,j]-ave.rank)^2)*
            sum((work.mat[,l]-ave.rank)^2))
        }
      }
      
      #Calculate the B matrix (qr helps with inversion)
      B.temp = solve(qr(C.temp))
      
      coeff.val = rep(0,p)
      
      #Calculate the PRCC
      for( j in 1:p ){
        
        coeff.val[j] = -B.temp[j,p+1]/sqrt(B.temp[j,j]*B.temp[p+1,p+1])
        
      }
      
      #Calculate the t-test statistics and p-values
      t.val = coeff.val*sqrt((n-2)/1-coeff.val)
      p.val = 2*pt(abs(t.val),df=(n-2),lower.tail=F)
      
      #Output the results
      results[[output.names[i]]] = data.frame(
        prcc = coeff.val,
        t.value = t.val,
        p.value = p.val,
        row.names = par.names)
    }
    
    return(results)
  }
  
  ########################################################################
  #Calculate the PRCCs using regression methods
  ########################################################################
  else if( routine == "regression" ){
    
    #Do the calculation for each model output
    for( i in 1:k ){
      
      coeff.val = rep(0,p)
      
      #Calculate the PRCC
      for( j in 1:p ){
        
        #Variation in output that can not be explained by all other predictors
        #(except the predictor of interest)
        fit.y = lm(output.rank[,i] ~ par.rank[,-j])
        
        #Variation in the predictor of interest that can not be explained
        #by the other predictors
        fit.x = lm(par.rank[,j] ~ par.rank[,-j])
        
        #PRCC is the correlation between the residuals of the two
        #regressions above
        coeff.val[j] = cor(fit.y$residuals,fit.x$residuals)
        
      }
      
      #Calculate the t-test statistics and p-values
      t.val = coeff.val*sqrt((n-2)/1-coeff.val)
      p.val = 2*pt(abs(t.val),df=(n-2),lower.tail=F)
      
      #Output the results
      results[[output.names[i]]] = data.frame(
        prcc = coeff.val,
        t.value = t.val,
        p.value = p.val,
        row.names = par.names)
    }
    
    return(results)
  }
  
  else{ return("Error: Calculation type is invalid. Must be either 'blower' or 'regression'") }
  
}


################################################################################
## plot.prcc - function to plot a prcc object
##
##    Arguments:
##          prcc.obj = a prcc object from the 'prcc' function
##             alpha = level of significance desired for cut-off lines
##               ... = additional arguments to be passed to functions called
##                     within this function
##
##
##    Output:
##      A plot that has a bar graph for each output variable giving the PRCCs.
##      Dashed red lines give the cutoff values for significance. It parameter
##      names are specified correctly, axis labels will be smart.

plot.prcc = function(prcc.obj,alpha=0.05,...){
  
  x11()
  par(mfrow=c(ceiling((length(prcc.obj)-1)/3),min(c(length(prcc.obj)-1,3))))
  
  for( i in 2:length(prcc.obj) ){
    
    names.list=dimnames(results[[i]])[[1]]
    
    #Bar graph with the prcc values. The function in names.arg converts character
    #strings containing the names of greek letters to actual greek letters
    barplot(prcc.obj[[i]][,"prcc"],
            names.arg=sapply(names.list,
                             function(x) as.expression(substitute(list(a),list(a=as.symbol(x))))),
            main=names(results)[i],ylab="PRCC",cex.lab=1.1,...)
    
    #Plot lines to show the cutoff values for the alpha level of significance
    t.cutoff=qt(alpha/2,df=prcc.obj$num.sets-2)
    sig.cutoffs=c( (-(t.cutoff^2)-sqrt((t.cutoff^4) + 4*(prcc.obj$num.sets-2)*(t.cutoff^2)))/(2*(prcc.obj$num.sets-2)),
                   (-(t.cutoff^2)+sqrt((t.cutoff^4) + 4*(prcc.obj$num.sets-2)*(t.cutoff^2)))/(2*(prcc.obj$num.sets-2)))
    abline(h=sig.cutoffs,lty=2,col="red")
    abline(h=0,lty=1,col="black")
  }
}

# prcc with stoch simulation results...
res<-cbind(X[,1],
           #X[,2],
           #X[,3],
           # X[,4], # pop extinction
           X[,5])
#res<-as.list(res)
dimnames(res)[[2]]<-c("Population",
                      #"Prevalence",
                      #"Adult seroprevalence",
                      #"Population persistence",
                      "LBV persistence")

# write.csv(res, "results_25yr1000Sens2.csv", row.names=F, na="")

results=prcc(par.mat=hypercube,model.output=res ## results matrix here...
               ,routine="blower" # NB removed par names so uses symbols, add [par.names="",]
             ,output.names=c("Population size",
                           # "Prevalence",
                           # "Adult seroprevalence",
                            #"Population persistence",
                             "Pathogen persistence"
                            ))

plot.prcc(results,ylim=c(-1,1),cex.sub=1.2)

## all works if enough paramets sets.... [otherwise get singular matrix warning..]
#########################
##
## plotting birth pulse

library(deSolve)

## this is the BP function to add....

b <-function(t,s=50,omega=2,phi=0,k=1.2){ ## adding birth pulse
  kappa<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
}

plot(b)#,add=T) ## no "pulsed" over the 0-1

ss<-seq(20,70,by=5)
test<-for (i in 1:length(ss)){
  res<-integrate((adf<-function(t,s=ss[i],omega=2,phi=0,k=1.2){
    epsilon<-1/s
    x<-cos(pi*omega*t-phi)
    k*(1/sqrt(epsilon*pi)*exp(-x^2/epsilon))
  }),0,1)
  res
}

## NB s == controls synchrony - 
##    omega == # pulses
##    phi == position
##    k == heat of pulse
## t == time obviously

intBt=integrate((adf<-function(t,s=30,omega=2,phi=0,k=1.2){
  epsilon<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(epsilon*pi)*exp(-x^2/epsilon))
}),0,1)

intBt # to get value for integral of B(t) this should be annual birth rate

B.t = function(t,s,k,omega=2,phi=0) {
  k*sqrt(s/pi)*exp(-s*cos(pi*t*omega-phi)^2)
}

# with the birth rate set so that it matches expected values
setK = function(t,s,omega=2,phi=0,br){
  k=br/integrate(B.t,0,1,s,1,omega=2,phi=0)$value  
  # then calculate B(t) with this calculated value of k
  B.t(t,s,k,omega,phi)
}

res<-setK(t=1,s=30,omega=2,phi=0,br=0.48)

b <-function(t=1,s=30,omega=2,phi=0,k=res){ ## adding birth pulse
  kappa<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
}
plot(b)



# Birth pulse function
B.t = function(t,s,k,omega=2,phi=0) {
  k*sqrt(s/pi)*exp(-s*cos(pi*t*omega-phi)^2)
}

# Function to calculate the birth pulse function, with the birth rate set so that it balances the death rate, m
setB.t = function(t,s,m,omega=2,phi=0){
  k=m/integrate(B.t,0,1,s,1,omega,phi)$value
  k
}

ss<-seq(7,170,by=10)
res<-matrix(NA,ncol=1,nrow=length(ss))
for (i in 1:nrow(res)){
  res[i,]<-setB.t(t=1,s=ss[i],m=0.48)  
  res
}
res

par.plot<-cbind(res,ss,rep(2,length(ss)),rep(0,length(ss)))
dimnames(par.plot)[[2]]<-c("k","s","omega","phi")
par.plot

par(omi=c(1,0.5,0.5,0.5))
par(mai=c(0.8,0.8,0.8,0.8))
par(mar=c(1, 1, 1, 1) + 0.1)
par(mfrow=c(1,1))

b <-function(t=1,s=par.plot[17,2],omega=2,phi=0,k=par.plot[17,1]){ ## adding birth pulse
  kappa<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
}
plot(b,xlab="",ylab="",#xaxt="n",yaxt="n",
     ylim=c(0,15),lwd=2,bty="n")

#b <-function(t=1,s=par.plot[15,2],omega=2,phi=0,k=par.plot[15,1]){ ## adding birth pulse
#  kappa<-1/s
#  x<-cos(pi*omega*t-phi)
#  k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
#}
#plot(b,add=T,col="black",lty=2)

b <-function(t=1,s=par.plot[1,2],omega=2,phi=0,k=par.plot[1,1]){ ## adding birth pulse
  kappa<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
}
plot(b,add=T,col="black",lty=2,lwd=2)

b <-function(t=1,s=par.plot[17,2],omega=1,phi=0,k=par.plot[1,1]){ ## adding birth pulse
  kappa<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
}
plot(b,add=T,col="darkgrey",lty=1,lwd=2)


b <-function(t=1,s=par.plot[1,2],omega=1,phi=0,k=par.plot[1,1]){ ## adding birth pulse
  kappa<-1/s
  x<-cos(pi*omega*t-phi)
  k*(1/sqrt(kappa*pi)*exp(-x^2/kappa))
}

plot(b,add=T,col="darkgrey",lty=2,lwd=2)

legend(x=0.25,y=16,lty=rep(1:2,2),lwd=rep(2,4),
       legend=c(expression(paste(omega==1, ", k = 1.44, s = 7"),paste(omega==1, ", k = 1.50, s = 167"),
                           paste(omega==2, ", k = 1.44, s = 7"),paste(omega==2, ", k = 1.50, s = 167"))),
       bty="n",col=c("darkgrey","darkgrey",1,1),cex=0.8)

mtext("1 year",side=1,outer=F,line=3)


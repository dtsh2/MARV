## MARV SEIR model with density-dept transmission
## 25 years Sens Analysis
## clean up

rm(list=ls())

# instal r developer toolbox first (Rtools from Cran R)

## pomp test run lbv
getwd()
setwd("~/GitHub/MARV") # revise as necessary
setwd("C:/Users/dtshayma/Dropbox")
library(pomp)

#Compiling C code and loading the dll
# dyn.unload("marburgmodel.dll")

# system("R CMD SHLIB marburgmodel.c")

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
  BETA=18,
  MU=0.000510492,
  DELTA=0.002312247,
  SIGMA=1/21,
  K=40000,
  EPSILON=1/365,
  TAU=1/24,
  KAPPA=1.5/365,
  S=14.35,
  OMEGA=1/365,
  PHI=0.0,
  SUSJ.0=4000,EXPJ.0=0,INFJ.0=100,RECJ.0=10000,
  SUSA.0=50000, EEXPA.0=100,INFA.0=5000, RECA.0=50000,
  PREVA.0= , PREVJ.0= ,
  SPA.0=0.4994506,SPJ.0=0.5882353) # this adds to the initial conditions given the state variables

sim <- simulate(seir,params=c(params),seed=3493885L,nsim=1,states=T,obs=F,as.data.frame=T) # 
class(seir) # pomp object
class(sim) # data frame - even if I remove "as.data.frame" in the above code (sim)
#sim <- simulate(sir,params=c(params),nsim=1,states=T,obs=F)#,as.data.frame=T) # saves as an array
# pf<-pfilter(sim,params=c(params),Np=1000) # won't work, because this is a data frame, not pomp object

plot(sim$time,sim$SUSJ,type="l")
points(sim$time,sim$RECJ,col="green",type="l")
points(sim$time,sim$EXPJ,col="brown",type="l")
points(sim$time,sim$INFJ,col="red",type="l")

plot(sim$time,sim$SUSA,type="l")
points(sim$time,sim$RECA,col="green",type="l")
points(sim$time,sim$EXPA,col="brown",type="l")
points(sim$time,sim$INFA,col="red",type="l")

plot(sim$time,sim$SPA,type="l",col="green",ylim=c(0,1))
points(sim$time,sim$SPJ,type="l",col="red")
points(sim$time,sim$SPJ,type="l",col="red")
points(sim$time,sim$SPJ,type="l",col="red")

#########################################################
######################### rest of this code requires the MLE estimate for the rho and beta params
## for LHS parameter set
## Calling requisite libraries for parallel computing

library(foreach)
library(doSNOW)

#Setting up "parallel backend"

w<-makeCluster(3,type="SOCK") # makes the cluster, i.e. no of cores ABC = 8 cores, DEF = 12 see performance to see # of cores
registerDoSNOW(w) # 

#Checks that the number of workers is set up correctly.

getDoParWorkers()

## ~~~~~~~~~~ LHS SAMPLING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#############
# NB not testing K - carrying cap, or rate of aging..

library(lhs)

nspaces=100 ## how many bins/intervals, 100 needed for PRCC

hypercube=randomLHS(n=nspaces, k=10) ## function with N columns
dimnames(hypercube)[[2]]=c("beta","mu","delta","alpha","rho",
                           "sigma","K",#"epsilon",
                           "tau","psi",#"k",
                           "s")#,
                           #"phi")  # named columns with parameter names

mins = c( 				    ## set mins for each parameters-exclude those not to be varied if any
  beta= 1,		     # transmission 
  mu= 0.0000510492,  			# natural mortality from my CMR study
  delta= 0.0002312247, 			# juvenile mortality rate
  alpha= 0.02,			# dis induced mortality
  rho= 0.001,				# probability that exposure/infection will lead to infection & infectiousness (and dead)
  sigma= 0.002083333,			# incubation period 
  K= 100000,			      # K, however, population size fluctuates...up to >1*10^6
#  epsilon= 1/10,				# rate of aging for those juveniles, should be ~ annual
  tau= 0.004166667, 	      # rate of seroconversion
  psi= 10,			######### this will need to be 1/psi for analysis 
  #k=,        			# nb this is birth rate which halved
  s=77.82/10)#,   # very synchronous
#  phi=0.01)

maxs = c( 				    ## set mins for each parameters-exclude those not to be varied if any
  beta= 30,           # transmission
  mu= 0.00510492,  	          # natural mortality from my CMR study
  delta= 0.02312247,            # juvenile mortality rate
  alpha= 0.9,		    # dis induced mortality
  rho= 0.4,			    # probability 
  sigma= 0.2,          # incubation period 
  K= 100000*10,			    # K, however, population size fluctuates...up to >1*10^6
#  epsilon= 1*10,	                # rate of aging for those juveniles, should be ~ annual
  tau= 0.2,         # rate of seroconversion
  psi= 200,      ###### will need to be 1/psi for the model
  #k=1.5,                       # nb this is peak
  s=77.82*2)#,
#  phi=0.99)                      

diffs=maxs-mins ## range of each variable

hypercubeadj = hypercube # create copy of hypercube samples to modify, hypercube adjusted; i.e. new matrix
for (i in 1:ncol(hypercube)){
  hypercubeadj[,i]=hypercube[,i]*diffs[i]+mins[i] # scale samples to difference and add minimum
}

head(hypercubeadj)
dim(hypercubeadj)

dimnames(hypercubeadj)[[2]]=c("beta","mu","delta","alpha","rho",
                              "sigma","K",#"epsilon",
                              "tau","psi",#"k",
                              "s")#,"phi")

nonVarying = matrix(c(
  #K= 1000000,			# K, however, population size fluctuates...up to >1*10^6
  epsilon= 1/365,	# rate of aging for those juveniles, should be ~ annual - per cap per year
  omega= 1/365, # for pb per year - 1 per year
  phi= 0.0, # for timing - ...
  k=1.5,
  SUSJ.0=4000,MDAJ.0=4000, 
  SUSJM.0=1000,EIJ.0=1000,
  ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, 
  EIA.0=100, ERA.0=1000,
  INFA.0=5000, RECA.0=50000),
              ncol=16,
              nrow=100,
              byrow=T) #binded with non-varying parameters


dimnames(nonVarying)[[2]]=c("omega","epsilon","phi","k","SUSJ.0","MDAJ.0","SUSJM.0","EIJ.0",
                            "ERJ.0","INFJ.0","RECJ.0","SUSA.0",
                            "EIA.0","ERA.0","INFA.0","RECA.0") # naming non-varying columns

gamma.day <- (365-hypercubeadj[,9])
gamma <-1/gamma.day
hypercubeadj[,9] <-1/hypercubeadj[,9]
#hypercubeadj[,1] <-hypercubeadj[,1]*hypercubeadj[,7]
fullParamSets = cbind(nonVarying,hypercubeadj,gamma) # full parameter set

head(fullParamSets)
dim(fullParamSets)

dimnames(fullParamSets)[[2]]=c("OMEGA","EPSILON","PHI","KAPPA",
                               "SUSJ.0","MDAJ.0","SUSJM.0",
                               "EIJ.0","ERJ.0","INFJ.0",
                               "RECJ.0","SUSA.0","EIA.0",
                               "ERA.0","INFA.0","RECA.0",
                               "BETA","MU","DELTA",
                               "ALPHA","RHO","SIGMA",
                               "K","TAU","PSI",
                              # "KAPPA",
                               "S","GAMMA")

# order for pomp/C model:  
BETA = fullParamSets[,17]
MU = fullParamSets[,18]
DELTA = fullParamSets[,19]
ALPHA = fullParamSets[,20]
RHO = fullParamSets[,21]
SIGMA = fullParamSets[,22]
K = fullParamSets[,23]
EPSILON = fullParamSets[,2]
TAU = fullParamSets[,24]
PSI = fullParamSets[,25]
KAPPA = fullParamSets[,4]
S = fullParamSets[,26]
OMEGA = fullParamSets[,1]
PHI = fullParamSets[,3]
GAMMA = fullParamSets[,27]
SUSJ.0 = fullParamSets[,5]
MDAJ.0 = fullParamSets[,6]
SUSJM.0 = fullParamSets[,7]
EIJ.0 = fullParamSets[,8]
ERJ.0 = fullParamSets[,9]
INFJ.0 = fullParamSets[,10]
RECJ.0 = fullParamSets[,11]
SUSA.0 = fullParamSets[,12]
EIA.0 = fullParamSets[,13]
ERA.0 = fullParamSets[,14]
INFA.0 = fullParamSets[,15]
RECA.0 = fullParamSets[,16]

paramset<-cbind(BETA,MU,DELTA,ALPHA,RHO,SIGMA,K,EPSILON,TAU,PSI,KAPPA,S,OMEGA,PHI,GAMMA,
                SUSJ.0, MDAJ.0, SUSJM.0,EIJ.0,ERJ.0,INFJ.0, RECJ.0, SUSA.0, EIA.0,ERA.0,INFA.0, RECA.0)

######################################################################################
results<-array(NA,dim=c(100,1,5))

# for one parameter set....
out1 <-simulate(sir,params=c(paramset[1,]),
                seed=1493885L,nsim=1000,states=T,obs=F,as.data.frame=T) #

outres1 <- out1[seq(from=9126,to=912600,by=9126),] # select last #s
N1 = array(0,c(1000,5)) # same dimensions as No. runs * outputs I want
for (i in 1:1000){ # each stochastic run
  N1[i,1]<-sum(outres1[i,1:12])
  N1[i,2]<-(sum(outres1[i,6],outres1[i,11])/sum(outres1[i,1:12]))*100 # prevalence; total
  N1[i,3]<-((outres1[i,12])/(sum(outres1[i,8:12])))*100 # adult seroprevalence; total
  N1[i,4]<-ifelse(sum(outres1[i,1:12])>0,1,0) # population extinct for each run
  N1[i,5]<-ifelse(sum(outres1[i,6],outres1[i,11])>0,1,0) # pathogen extinction for each run
}
N1[is.na(N1)]<- 0
## now average
M1 = array(0,c(1,5))
M1[1] = mean(N1[1:1000,1]) # population size
M1[2] = mean(N1[1:1000,2]) # prevalence
M1[3] = mean(N1[1:1000,3]) # adult seroprevalence
M1[4] = mean(N1[1:1000,4]) # pop ext
M1[5] = mean(N1[1:1000,5]) # lbv ext
rm(out1)
M1
results[1,,]<-M1
results[1,,]

##########################################################33
# for all parameter sets....
results<-array(NA,dim=c(100,1,5))

for (j in 1:length(paramset[,1])){
out <-simulate(sir,params=c(paramset[j,]),
                seed=1493885L,nsim=1000,states=T,obs=F,as.data.frame=T) #
outres <- out[seq(from=9126,to=912600,by=9126),] # select last #s
N = array(0,c(1000,5)) # same dimensions as No. runs * outputs I want
for (i in 1:1000){ # each stochastic run
  N[i,1]<-sum(outres[i,1:12])
  N[i,2]<-(sum(outres[i,6],outres[i,11])/sum(outres[i,1:12]))*100 # prevalence; total
  N[i,3]<-((outres[i,12])/(sum(outres[i,8:12])))*100 # adult seroprevalence; total
  N[i,4]<-ifelse(sum(outres[i,1:12])>0,1,0) # population persistence for each run
  N[i,5]<-ifelse(sum(outres[i,6],outres[i,11])>0,1,0) # pathogen persistence for each run
}
N[is.na(N)]<- 0
## now average
M = array(0,c(1,5))
M[1] = mean(N[1:1000,1]) # population size
M[2] = mean(N[1:1000,2]) # prevalence
M[3] = mean(N[1:1000,3]) # adult seroprevalence
M[4] = mean(N[1:1000,4]) # mean pop persistence
M[5] = mean(N[1:1000,5]) # mean path persistence
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
res<-cbind(X[,1],X[,2],
           X[,3],
            X[,4], # pop extinction
           X[,5])
dimnames(res)[[2]]<-c("Population","Prevalence",
                      "Adult seroprevalence",
                      "Population persistence",
                      "LBV persistence")

 write.csv(res, "results_25yr1000Sens2.csv", row.names=F, na="")

res<-cbind(X[,1],#X[,2],
           #X[,3],
#            X[,4], # pop extinction
           X[,5])
## res<-read.csv("results_25yrSens.csv")#, row.names=F, na="")

results=prcc(par.mat=hypercube,model.output=res ## results matrix here...
               ,routine="blower" # NB removed par names so uses symbols, add [par.names="",]
             ,output.names=c("Population size",
                            #"Prevalence",
                             #"Adult seroprevalence",
                             #"Population persistence",
                             "LBV persistence"))

plot.prcc(results,ylim=c(-1,1),cex.sub=1.2)

## all works if enough paramets sets.... [otherwise get singular matrix warning..]

# NB not testing K - carrying cap, or rate of aging..
# no LHS needed


nonVarying = matrix(c(
  BETA=10,
  MU=0.000510492,
  DELTA=0.002312247,
  ALPHA=0.2,
  RHO=0.02,
  SIGMA=1/48,
  #K=1000000,
  EPSILON=1/365,
  TAU=1/24,
  KAPPA=1.5/365,
  S=77.82,
  OMEGA=1/365,
  PHI=0.5,
  SUSJ.0=4000,EIJ.0=1000,ERJ.0=1000,INFJ.0=1000,
  RECJ.0=10000,SUSA.0=50000, EIA.0=100,
  ERA.0=1000,INFA.0=5000, RECA.0=50000),
  ncol=22,
  nrow=40,
  byrow=T) #binded with non-varying parameters


dimnames(nonVarying)[[2]]=c("BETA",
                            "MU",
                            "DELTA",
                            "ALPHA",
                            "RHO",
                            "SIGMA",
                            "EPSILON",
                            "TAU",
                            "KAPPA",
                            "S",
                            "OMEGA",
                            "PHI",
                            "SUSJ.0",
                            "EIJ.0","ERJ.0","INFJ.0",
                            "RECJ.0","SUSA.0", "EIA.0",
                            "ERA.0","INFA.0","RECA.0") # naming non-varying columns

## from other code
Kset=seq(from = 100, to=200000, by =5000)

fullParamSets = cbind(nonVarying,Kset) # full parameter set
#fullParamSets[,1] <- fullParamSets[,1]*fullParamSets[,29]
head(fullParamSets)
dim(fullParamSets)

dimnames(fullParamSets)[[2]]=c("BETA","MU","DELTA","ALPHA","RHO",
                               "SIGMA","EPSILON","TAU",
                               "KAPPA","S","OMEGA","PHI",
                               "SUSJ.0","EIJ.0",
                               "ERJ.0","INFJ.0","RECJ.0","SUSA.0",
                               "EIA.0","ERA.0","INFA.0","RECA.0",
                               "K")

# order for pomp/C model:  
BETA = fullParamSets[,1]
MU = fullParamSets[,2]
DELTA = fullParamSets[,3]
ALPHA = fullParamSets[,4]
RHO = fullParamSets[,5]
SIGMA = fullParamSets[,6]
K = fullParamSets[,23]
EPSILON = fullParamSets[,7]
TAU = fullParamSets[,8]
KAPPA = fullParamSets[,9]
S = fullParamSets[,10]
OMEGA = fullParamSets[,11]
PHI = fullParamSets[,12]
SUSJ.0 = fullParamSets[,13]
EIJ.0 = fullParamSets[,14]
ERJ.0 = fullParamSets[,15]
INFJ.0 = fullParamSets[,16]
RECJ.0 = fullParamSets[,17]
SUSA.0 = fullParamSets[,18]
EIA.0 = fullParamSets[,19]
ERA.0 = fullParamSets[,20]
INFA.0 = fullParamSets[,21]
RECA.0 = fullParamSets[,22]

paramset<-cbind(BETA,MU,DELTA,ALPHA,RHO,SIGMA,K,EPSILON,TAU,KAPPA,S,OMEGA,PHI,
                SUSJ.0,EIJ.0,ERJ.0,INFJ.0, RECJ.0, SUSA.0, EIA.0,ERA.0,INFA.0, RECA.0)
######################################################################################

## Calling requisite libraries for parallel computing

library(foreach)
library(doSNOW)

#Setting up "parallel backend"

w<-makeCluster(3,type="SOCK") # makes the cluster, i.e. no of cores ABC = 8 cores, DEF = 12 see performance to see # of cores
registerDoSNOW(w) # 

#Checks that the number of workers is set up correctly.

getDoParWorkers()

#######################################################3
results<-array(NA,dim=c(40,1,5))

# for one parameter set....
out1 <-simulate(sir,params=c(paramset[1,]),
                seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #

outres1 <- out1[seq(from=10951,to=1095100,by=10951),] # select last #s
N1 = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
for (i in 1:100){ # each stochastic run
  N1[i,1]<-sum(outres1[i,1:10])
  N1[i,2]<-(sum(outres1[i,4],outres1[i,9])/sum(outres1[i,1:10]))*100 # prevalence; total
  N1[i,3]<-((outres1[i,10])/(sum(outres1[i,6:10])))*100 # adult seroprevalence; total
  N1[i,4]<-ifelse(sum(outres1[i,1:10])>0,1,0) # population extinct for each run
  N1[i,5]<-ifelse(sum(outres1[i,4],outres1[i,9])>0,1,0) # pathogen extinction for each run
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
results<-array(NA,dim=c(40,1,5))

for (j in 1:length(paramset[,1])){
  out <-simulate(sir,params=c(paramset[j,]),
                 seed=1493885L,nsim=100,states=T,obs=F,as.data.frame=T) #
  outres <- out[seq(from=10951,to=1095100,by=10951),] # select last #s
  N = array(0,c(100,5)) # same dimensions as No. runs * outputs I want
  for (i in 1:100){ # each stochastic run
    N[i,1]<-sum(outres[i,1:10])
    N[i,2]<-(sum(outres[i,4],outres[i,9])/sum(outres[i,1:10]))*100 # prevalence; total
    N[i,3]<-((outres[i,10])/(sum(outres[i,6:10])))*100 # adult seroprevalence; total
    N[i,4]<-ifelse(sum(outres[i,1:10])>0,1,0) # population extinct for each run
    N[i,5]<-ifelse(sum(outres[i,4],outres[i,9])>0,1,0) # pathogen extinction for each run
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
dim(X)<-c(40,5)
head(X)
tail(X)

################################################################################
# below for K

plot(X[,1],X[,2])
# plot

par(mar=c(5, 6, 4, 4) + 0.1)

plot(X[,1],X[,2],pch=16,
     ylab="Mean prevalence",xlab="Population size",
     col="grey25", cex.lab=1.2)

plot(X[,1],X[,3],pch=16,
     ylab="Mean seroprevalence",xlab="Population size",
     col="grey25", cex.lab=1.2)

plot(X[,1],X[,5],pch=16,
     ylab="P[persist]",xlab="Population size",
     col="grey25", cex.lab=1.2)

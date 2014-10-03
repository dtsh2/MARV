rm(list=ls())
getwd()

res<-read.csv("data.csv",header=T)

datres<-cbind(res$Positive,res$Negative)
model<-glm(datres~res$Pulses,binomial,na.action=na.omit)
model
summary(model)
plot(model)

model<-glm(datres~res$Births,binomial,na.action=na.omit)
model
summary(model)
plot(model)

###############################

sp<-function(x,y){
  sp<-x/y
}
sp<-res$Positive/res$Tested
res.new<-rbind(res$Positive,res$Tested,sp)
rownames(res.new)<-c("P","N","Sp")
res.new<-t(res.new)
class(res.new)

res.n<-res.new[c(-15),] 

res.u<-matrix(NA,ncol=70,nrow=2)
for (ii in 1:70){
  res.u[1,ii]<-binom.test(res.n[ii],res.n[ii+70])$conf.int[1]
  res.u[2,ii]<-binom.test(res.n[ii],res.n[ii+70])$conf.int[2]
}

res.u<-t(res.u)
res.all<-cbind(res.n,res.u)
dim(res.all)
colnames(res.all)<-c("P","N","Sp","L","U")

res.all<-as.data.frame(res.all)

## each sample get CIs

barplot(res.all$Sp,ylim=c(0,1))

####
library(gplots)
barplot2(res.all$Sp, plot.ci=TRUE, ci.l=res.all$L, ci.u=res.all$U,
         ylim=c(0,1))

## Add legends ##########

####
## select those with one birth pulse
single <- res[(!is.na(res$Births)) & res$Births==0,]
double <- res[(!is.na(res$Births)) & res$Births==1,]

single.sp<-single$Positive/single$Tested
single.res.new<-rbind(single$Positive,single$Tested,single.sp)
rownames(single.res.new)<-c("P","N","Sp")
single.res.new<-t(single.res.new)
class(single.res.new)
dim(single.res.new)
single.res.u<-matrix(NA,ncol=26,nrow=2)
for (ii in 1:26){
  single.res.u[1,ii]<-binom.test(single.res.new[ii],single.res.new[ii+26])$conf.int[1]
  single.res.u[2,ii]<-binom.test(single.res.new[ii],single.res.new[ii+26])$conf.int[2]
}

single.res.u<-t(single.res.u)
single.res.all<-cbind(single.res.new,single.res.u)
dim(single.res.all)
colnames(single.res.all)<-c("P","N","Sp","L","U")

single.res.all<-as.data.frame(single.res.all)

## each sample get CIs

barplot(single.res.all$Sp,ylim=c(0,1))

####
library(gplots)
barplot2(single.res.all$Sp, plot.ci=TRUE, ci.l=single.res.all$L, ci.u=single.res.all$U,
         ylim=c(0,1))

##
double.sp<-double$Positive/double$Tested
double.res.new<-rbind(double$Positive,double$Tested,double.sp)
rownames(double.res.new)<-c("P","N","Sp")
double.res.new<-t(double.res.new)
class(double.res.new)

## remove 2/0
double.res.n<-double.res.new[c(-9),] 
dim(double.res.n)
double.res.u<-matrix(NA,ncol=35,nrow=2)
for (ii in 1:35){
  double.res.u[1,ii]<-binom.test(double.res.n[ii],double.res.n[ii+35])$conf.int[1]
  double.res.u[2,ii]<-binom.test(double.res.n[ii],double.res.n[ii+35])$conf.int[2]
}

double.res.u<-t(double.res.u)
double.res.all<-cbind(double.res.n,double.res.u)
dim(double.res.all)
colnames(double.res.all)<-c("P","N","Sp","L","U")

double.res.all<-as.data.frame(double.res.all)

## each sample get CIs

barplot(double.res.all$Sp,ylim=c(0,1))

####
library(gplots)
barplot2(double.res.all$Sp, plot.ci=TRUE, ci.l=double.res.all$L, ci.u=double.res.all$U,
         ylim=c(0,1))

par(mfrow=c(1,2))
barplot2(double.res.all$Sp, plot.ci=TRUE, ci.l=double.res.all$L, ci.u=double.res.all$U,
         ylim=c(0,1),main="Double",ci.lty=2)

barplot2(single.res.all$Sp, plot.ci=TRUE, ci.l=single.res.all$L, ci.u=single.res.all$U,
         ylim=c(0,1),main="Single",ci.lty=2)

### check without repetition 
res<-read.csv("species.csv",header=T)
datres<-cbind(res$Positive,res$Negative)
#model<-glm(datres~res$Births,quasibinomial,na.action=na.omit)
#model
#summary(model)
#plot(model)
#model<-glm(datres~res$Births+res$Fruit.Bat,binomial,na.action=na.omit)
#summary(model)

model<-glm(datres~res$Births*res$Fruit.Bat,quasibinomial,na.action=na.omit)
summary(model)
model1<-update(model,~.-res$Births:res$Fruit.Bat)
anova(model,model1,test="F")
## interaction not significant
summary(model1)
model2<-update(model1,~.-res$Fruit.Bat)
anova(model1,model2,test="F")
summary(model2)
coef(model2)
1/(1+1/exp(coef(model2)[1]))
1/(1+1/exp(coef(model2)[1]+coef(model2)[2]))
res.t<-tapply(predict(model2,type="response"),na.omit(res$Births),mean)
res.t[2]/res.t[1]

## chi squared test
filo.positive<-c(587,12)
filo.total<-c(7857,673)
prop.test(filo.positive,filo.total)
filo.negative<-c(7270,661)
res.test<-cbind(filo.positive,filo.negative)
fisher.test(res.test)
chisq.test(res.test)

## test without R. aegypt
filo.positive<-c(215,12)
filo.total<-c(4559,673)
prop.test(filo.positive,filo.total)
filo.negative<-c(4344,661)
res.test<-cbind(filo.positive,filo.negative)
fisher.test(res.test)
chisq.test(res.test)


######

## fruit bat vs births
bres<-res[-c(8,18,48,51),]
birthtable<-table(bres$Births,bres$Fruit.Bat)
birthtable
fisher.test(birthtable)
###############################

sp<-function(x,y){
  sp<-x/y
}
sp<-res$Positive/res$Tested
res.new<-rbind(res$Positive,res$Tested,sp)
rownames(res.new)<-c("P","N","Sp")
res.new<-t(res.new)
class(res.new)

res.n<-res.new#[c(-2),] 
dim(res.n)
res.u<-matrix(NA,ncol=53,nrow=2)
for (ii in 1:53){
  res.u[1,ii]<-binom.test(res.n[ii],res.n[ii+53])$conf.int[1]
  res.u[2,ii]<-binom.test(res.n[ii],res.n[ii+53])$conf.int[2]
}

res.u<-t(res.u)
res.all<-cbind(res.n,res.u)
dim(res.all)
colnames(res.all)<-c("P","N","Sp","L","U")

res.all<-as.data.frame(res.all)

## each sample get CIs

barplot(res.all$Sp,ylim=c(0,1))

####
library(gplots)
barplot2(res.all$Sp, plot.ci=TRUE, ci.l=res.all$L, ci.u=res.all$U,
         ylim=c(0,1))

## Add legends ##########

####
## select those with one birth pulse
single <- res[(!is.na(res$Births)) & res$Births==0,]
double <- res[(!is.na(res$Births)) & res$Births==1,]

single.sp<-single$Positive/single$Tested
single.res.new<-rbind(single$Positive,single$Tested,single.sp)
rownames(single.res.new)<-c("P","N","Sp")
single.res.new<-t(single.res.new)
class(single.res.new)
dim(single.res.new)
single.res.u<-matrix(NA,ncol=24,nrow=2)
for (ii in 1:24){
  single.res.u[1,ii]<-binom.test(single.res.new[ii],single.res.new[ii+24])$conf.int[1]
  single.res.u[2,ii]<-binom.test(single.res.new[ii],single.res.new[ii+24])$conf.int[2]
}

single.res.u<-t(single.res.u)
single.res.all<-cbind(single.res.new,single.res.u)
dim(single.res.all)
colnames(single.res.all)<-c("P","N","Sp","L","U")

single.res.all<-as.data.frame(single.res.all)

## each sample get CIs

barplot(single.res.all$Sp,ylim=c(0,1))

####
library(gplots)
barplot2(single.res.all$Sp, plot.ci=TRUE, ci.l=single.res.all$L, ci.u=single.res.all$U,
         ylim=c(0,1))

##
double.sp<-double$Positive/double$Tested
double.res.new<-rbind(double$Positive,double$Tested,double.sp)
rownames(double.res.new)<-c("P","N","Sp")
double.res.new<-t(double.res.new)
class(double.res.new)
double.res.new
## remove 2/0
double.res.n<-double.res.new#[c(-2),] 
dim(double.res.n)
double.res.u<-matrix(NA,ncol=20,nrow=2)
for (ii in 1:20){
  double.res.u[1,ii]<-binom.test(double.res.n[ii],double.res.n[ii+20])$conf.int[1]
  double.res.u[2,ii]<-binom.test(double.res.n[ii],double.res.n[ii+20])$conf.int[2]
}

double.res.u<-t(double.res.u)
double.res.all<-cbind(double.res.n,double.res.u)
dim(double.res.all)
colnames(double.res.all)<-c("P","N","Sp","L","U")

double.res.all<-as.data.frame(double.res.all)

## each sample get CIs

barplot(double.res.all$Sp,ylim=c(0,1))

####
library(gplots)
barplot2(double.res.all$Sp, plot.ci=TRUE, ci.l=double.res.all$L, ci.u=double.res.all$U,
         ylim=c(0,1))

par(mfrow=c(1,2))
barplot2(double.res.all$Sp, plot.ci=TRUE, ci.l=double.res.all$L, ci.u=double.res.all$U,
         ylim=c(0,1),main="",ci.lty=2,names.arg=c(1:20))
mtext('a',font=2,side=3,line=2,at=-4,cex=1)
names.arg=c(double$Bat.species)

barplot2(single.res.all$Sp, plot.ci=TRUE, ci.l=single.res.all$L, ci.u=single.res.all$U,
         ylim=c(0,1),main="",ci.lty=2,names.arg=c(1:24))
mtext('b',font=2,side=3,line=2,at=-4,cex=1)

## test using species names
par(omi=c(1.5,0.2,0.2,0.2))
par(mai=c(1.5,1.5,0.5,0.5))
par(mar=c(4, 4, 4, 4) + 0.1)
par(mfrow=c(1,2))
par(las=2)
barplot2(double.res.all$Sp, plot.ci=TRUE, ci.l=double.res.all$L, ci.u=double.res.all$U,
         ylim=c(0,1),main="",ci.lty=2,names.arg=double$Bat.Species,cex.names=0.6)
par(las=1)
mtext('a',font=2,side=3,line=2,at=-4,cex=1)
par(las=2)
barplot2(single.res.all$Sp, plot.ci=TRUE, ci.l=single.res.all$L, ci.u=single.res.all$U,
         ylim=c(0,1),main="",ci.lty=2,names.arg=single$Bat.Species,cex.names=0.6)
par(las=1)
mtext('b',font=2,side=3,line=2,at=-4,cex=1)

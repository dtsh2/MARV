## just for revised prcc figure
test<-read.csv("prcc_25yr1000Sens2_v2.csv")
names.list=as.character(test$X)
names.list[7]<-"kappa"
#Bar graph with the prcc values. The function in names.arg converts character
#strings containing the names of greek letters to actual greek letters
barplot(test$prcc,
        names.arg=sapply(names.list,
                         function(x) as.expression(substitute(list(a),list(a=as.symbol(x))))),ylab="PRCC",cex.lab=1.1,
        ylim=c(-1,1))
abline(h=0.09,lty=2,col="red")
abline(h=0)
abline(h=-0.09,lty=2,col="red")

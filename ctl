# Premer, M.I. PhD Dissertation 
# Chapter 4 Analysis
# Examining CTL patterns according to site varaibles that covary with distance to CTL trail
# February 25, 2015

require(nlme)
require(MASS)
require(lattice)
require(dplyr)
require(ggplot2)

ctl.models <- read.csv("C:/Users/Mike/Desktop/Dissertation/Chapter4/raw/csv/final/WTH_Plot_CTL_REGEN.csv",head=TRUE)
ctl.models <- ctl.models[which(ctl.models$treat=="ctl"), ]
ctl.models["log.height"] <- log(ctl.models$max.ht)
ctl.models["log.density"] <- log(ctl.models$stock)
ctl.models["log.density.acres"] <- log(ctl.models$stock.ac)

plot(ctl.models$age,ctl.models$log.height,data=ctl.models)
plot(ctl.models$age,ctl.models$log.height,data=ctl.models,cex=2.5)
plot(ctl.models$age,ctl.models$ht.ft)

subsummary <- ctl.models %>%
  group_by(Subplot,age)%>%
  summarize(mean.ht = mean(ht.ft))

ctlsubsum1 <- filter(subsummary, Subplot == "1")
ctlsubsum2 <- filter(subsummary, Subplot == "2")
ctlsubsum3 <- filter(subsummary, Subplot =="3")

par(mar=c(5,5,5,5))
plot(ctlsubsum1$age,ctlsubsum1$mean.ht,cex=3,ylim=c(0,10),xlim=c(0,10), pch=1,
     ylab="Mean dominant height (ft)",xlab="Stand age (years)",cex.axis=1.5,cex.lab=1.5,col="black")
points(ctlsubsum2$age,ctlsubsum2$mean.ht,cex=3,pch=19,col="gray30")
points(ctlsubsum3$age,ctlsubsum3$mean.ht,cex=3,pch=19,col="gray60")
nd <- data.frame(age=seq(0.1,12,0.1))
CTL3mod <- 0+1.4793*nd$age
points(nd$age,CTL3mod,type="l",lwd=3,lty=3)
CTL1mod <- 0+1.1891*nd$age
points(nd$age,CTL1mod,type="l",lwd=3,lty=1)

legend(0,10,c("CTL center","CTL 3.3 ft / CTL 16.4 ft  "),lty=c(1,4),col=c("black","black"),lwd=c(3,3),cex=1.25)
text(9.25,9.5,labels=expression("(a)") ,pos=4,cex=2.5)

ctl.models <- groupedData(log.height~age|stand,data=ctl.models)
nullhtmod <- lme(log.height~age-1,random=~1|Stand/Plot,na.action="na.omit",method="REML",
                 weights=varIdent(form=~1|Stand),data=ctl.models)
fullhtmod1 <- lme(log.height~age*Subplot-1,random=~1|Stand/Plot,na.action="na.omit",method="REML",
                  weights=varIdent(form=~1|Stand),data=ctl.models)
fullhtmod2 <- lme(log.height~age*bulk.density*ddw.mass-1,random=~1|Stand/Plot,na.action="na.omit",method="REML",
                  weights=varIdent(form=~1|Stand),data=ctl.models)
fullhtmod3 <- lme(log.height~age*bulk.density-1,random=~1|Stand/Plot,na.action="na.omit",method="REML",
                  weights=varIdent(form=~1|Stand),data=ctl.models)
fullhtmod4 <- lme(log.height~age*ddw.mass-1,random=~1|Stand/Plot,na.action="na.omit",method="REML",
                  weights=varIdent(form=~1|Stand),data=ctl.models)

# fullhtmod better
anova(fullhtmod1)
summary(fullhtmod2)
anova(fullhtmod3)
anova(fullhtmod4)

AIC(fullhtmod2)
x <- c(AIC(fullhtmod3), AIC(fullhtmod4), AIC(fullhtmod2))
delta <- x-min(x)
L <- exp(-0.5*delta)
w <- L/sum(L)

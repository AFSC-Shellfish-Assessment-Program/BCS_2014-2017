

require(ggplot2)
require(GGally)
require(reshape2)
require(lme4)
require(compiler)
require(parallel)
require(boot)
require(mgcv)
library(ResourceSelection)
require(visreg)
require(glmulti)
#require(LogisticDx)


#getting data ready for analysis
aa1=read.csv("bcs14.csv")
str(aa1)
class(aa1)
aa1Classes=as.vector(lapply(aa1,class))
aa1Classes$Size="numeric"
aa1Classes

aa2=read.csv("bcs15.csv",colClasses=aa1Classes)
str(aa2)
colnames(aa2)=colnames(aa1)

aa3=read.csv("bcs16.csv",colClasses=aa1Classes)
colnames(aa3)=colnames(aa1)

cbind(colnames(aa1),colnames(aa2),colnames(aa3))

aa=rbind(aa1,aa2,aa3)
str(aa)


str(aa)
substr(aa$Specific_Location,10,10)
aa$index=as.numeric(substr(aa$Specific_Location,10,10))
#check for some Specific_Location not in index sites

str(aa)
#it appears there are some crabs from outside index areas?
unique(aa$index)
idxs=is.na(aa$index)
which(idxs==T)
aa=aa[!idxs,]
str(aa)
aa=aa[aa$PCR_result==0 | aa$PCR_result==1,]
str(aa)
unique(aa$Species_Name)
idxs=is.na(aa$Species_Name)
aa=aa[!idxs,]
str(aa)
#taking out clutches larger then 0.
unique(aa$Clutch)
idxs=which(aa$Clutch>0)
aa=aa[-idxs,]
str(aa)



#changingto appriopiate classes for glm analysis
aa$Sex=as.factor(aa$Sex)
aa$index=as.factor(aa$index)
aa$Year=as.factor(aa$Year)
aa$BOTTOM_DEPTH=as.numeric(aa$BOTTOM_DEPTH)
aa$PCR_result=as.factor(aa$PCR_result)


str(aa)
unique(aa$Species_Name)
bairdi=aa[aa$Species_Name=="Chionoecetes bairdi",]
bairdi=bairdi[which(bairdi$index==1|bairdi$index==2|bairdi$index==3),]
opilio=aa[aa$Species_Name=="Chionoecetes opilio",]
opilio=opilio[which(opilio$index==4|opilio$index==5|opilio$index==6),]
bairdi$PCR_result
str(bairdi)

bcpue=read.csv("cpue_bairdi.csv")
bcpue=bcpue[,c(1,2,4,9,10)]
str(bcpue)
bairdi=merge(bairdi,bcpue)
ocpue=read.csv("cpue_opilio.csv")
ocpue=ocpue[,c(1,2,4,9,10)]
str(ocpue)
opilio=merge(opilio,ocpue)

bairdi1=as.data.frame(cbind(RR=bairdi$PCR_result,Size=bairdi$Size,
	BT=bairdi$Bottom_Temp,ST=bairdi$SURFACE_TEMP,
	BD=bairdi$BOTTOM_DEPTH,II=bairdi$index,
	Sex=bairdi$Sex,YY=bairdi$Year,
	MM=bairdi$Maturity,DN=log(bairdi$CPUE_NOHA),DB=log(bairdi$CPUE_KGHA)))
str(bairdi1)
opilio1=as.data.frame(cbind(RR=opilio$PCR_result,Size=opilio$Size,
	BT=opilio$Bottom_Temp,ST=opilio$SURFACE_TEMP,
	BD=opilio$BOTTOM_DEPTH,II=opilio$index,
	Sex=opilio$Sex,YY=opilio$Year,
	MM=opilio$Maturity,DN=log(opilio$CPUE_NOHA),DB=log(opilio$CPUE_KGHA)))
str(opilio1)

bairdi1$RR=as.integer(bairdi1$RR-1)
bairdi1$Sex=as.factor(bairdi1$Sex)
bairdi1$II=as.factor(bairdi1$II)
bairdi1$YY=as.factor(bairdi$Year)
bairdi1$MM=as.factor(bairdi1$MM)
str(bairdi1)

opilio1$RR=as.integer(opilio1$RR-1)
opilio1$Sex=as.factor(opilio1$Sex)
opilio1$II=as.factor(opilio1$II)
opilio1$YY=as.factor(opilio$Year)
opilio1$MM=as.factor(opilio1$MM)
str(opilio1)



#full model glmulti
#runs very long can be omitted
#======================================================================
#======================================================================
#======================================================================
#DB and DN were ommited because of confounding with enioronmental data

bairdi_sele <-
    glmulti(RR ~ Size+BT+ST+BD+II+YY+MM+Sex, 
		data = bairdi1,
            level = 1,               # No interaction considered
            method = "h",            # Exhaustive approach
            crit = "aic",            # AIC as criteria
            confsetsize = 5,         # Keep 5 best models
            plotty = F, report = F,  # No plot or interim reports
            fitfunction = "glm",     # glm function
            family = binomial)       # binomial family for logistic regression

## Show 5 best models (Use @ instead of $ for an S4 object)
bairdi_sele@formulas


## Show result for the best model
summary(bairdi_sele)
summary(bairdi_sele@objects[[1]])

#full model glmulti AIC
bairdi_sele <-
    glmulti(RR ~ Size+BT+ST+BD+II+YY+MM+Sex, 
		data = bairdi1,
            level = 2,               # No interaction considered
            method = "g",            # Exhaustive approach
            crit = "aic",            # AIC as criteria
            confsetsize = 5,         # Keep 5 best models
            plotty = F, report = F,  # No plot or interim reports
            fitfunction = "glm",     # glm function
            family = binomial)       # binomial family for logistic regression

## Show 5 best models (Use @ instead of $ for an S4 object)
bairdi_sele@formulas
str(lbw)

## Show result for the best model
summary(bairdi_sele)
summary(bairdi_sele@objects[[1]])

#full model glmulti BIC
bairdi_sele.bic <-
    glmulti(RR ~ Size+BT+ST+BD+II+YY+MM+Sex, 
		data = bairdi1,
            level = 2,               # No interaction considered
            method = "g",            # Exhaustive approach
            crit = "bic",            # AIC as criteria
            confsetsize = 5,         # Keep 5 best models
            plotty = F, report = F,  # No plot or interim reports
            fitfunction = "glm",     # glm function
            family = binomial)       # binomial family for logistic regression

## Show 5 best models (Use @ instead of $ for an S4 object)
bairdi_sele.bic@formulas
str(lbw)

## Show result for the best model
summary(bairdi_sele.bic)
summary(bairdi_sele.bic@objects[[1]])

save(bairdi_sele,file="sele.aic")
save(bairdi_sele.bic,file="bairdi_sele.bic")

load("bairdi_sele.bic")
summary(bairdi_sele.bic)
#======================================================================
#======================================================================
#======================================================================

#OPILIO
#full model glmulti BIC
opilio_sele.bic <-
    glmulti(RR ~ Size+BT+ST+BD+DN+DB+II+YY+MM+Sex, 
		data = opilio1,
            level = 2,               # No interaction considered
            method = "g",            # Exhaustive approach
            crit = "bic",            # AIC as criteria
            confsetsize = 5,         # Keep 5 best models
            plotty = F, report = F,  # No plot or interim reports
            fitfunction = "glm",     # glm function
            family = binomial)       # binomial family for logistic regression

## Show 5 best models (Use @ instead of $ for an S4 object)
opilio_sele.bic@formulas
#str(lbw)

## Show result for the best model
summary(opilio_sele.bic)
summary(opilio_sele.bic@objects[[1]])

save(opilio_sele,file="opilio_sele.aic")
save(opilio_sele.bic,file="opilio_sele.bic")

load("opilio_sele.bic")
summary(opilio_sele.bic)
#======================================================================
#======================================================================
#======================================================================

#final glm model for bairdi
mod2 <- glm(RR ~ 1 + Size + ST + ST:BT + BD:ST + II:Size + YY:ST,
	 data = bairdi1, family = binomial)
summary(mod2)
plot(mod2)

#ordering variables from most to least important
anova(mod2)
summary(mod2)$dispersion

#checking for overdispersion
mod2int <- glm(RR ~ 1 + Size + ST + ST:BT + BD:ST + II:Size + YY:ST,
	data = bairdi1, family = quasibinomial)
summary(mod2int)
# no need to model overdispersion as the dispersion is = 1
# Dispersion parameter for quasibinomial family taken to be 1.007


aa1$STATIONID=as.factor(aa$STATIONID)
#mixed model final
#this one does not converge
mod2mixed=glmer(RR ~ 1 + Size + ST + ST:BT + BD:ST + II:Size + YY:ST
	+(1|STATIONID), data = bairdi1, family = binomial,
	control = glmerControl(optimizer="Nelder_Mead"))
summary(mod2mixed)
BIC(mod2mixed)
#is(mod2mixed, "merPredD")

# Instead we can also test for overdispersion
pchisq(summary(mod2)$dispersion * mod2$df.residual, 
       mod2$df.residual, lower = F)
#p-value above 0.05 indicate no overdispersion
resids=residuals(mod2)
preds=fitted(mod2)
plot(resids~preds)
lines(smooth.spline(preds, resids, df = 10), lty = 2, col = "red")
qqnorm(resids)


#ploting predictor data
pdf("predictor_data_bairdi.pdf",w=10,h=10)
ggpairs(bairdi1[, c("RR","YY","II","Size","BD","ST","BT")])
dev.off()
head(bairdi1)

#predictors contriubtion
str(bairdi)
IIsum=unique(bairdi[,c(7,17,32,34:36,38)])
str(IIsum)
avgs=aggregate(cbind(BD=BOTTOM_DEPTH,ST=SURFACE_TEMP,
	Size=Size,BT=Bottom_Temp)~Year+index,data=IIsum,mean)

summary(mod2)
anova(mod2)
print(mod2)
confint(mod2) 
predict(mod2)
 
termplot(mod2,se=T,ylim=c(-5,5))

#pdf("predictions_bairdi.pdf")
#visreg(mod2,print.cond=T,scale="response",
	ylab="Prevalence")
pdf("predictions_bairdi_yy_ii.pdf",width=8, height=6)
boxplot(predict(mod2,newdata=bairdi1,type="response")~mod2$data$II+mod2$data$YY)

dev.off()
pdf("predictions_bairdi_size.pdf",width=8, height=8)
par(mfrow=c(3,3))
for(ii in 1:9){
visreg(mod2,"Size",type="conditional",scale="response",
	print.cond=T,cond=list(II=avgs$index[ii],YY=avgs$Year[ii],
	BD=avgs$BD[ii],ST=avgs$ST[ii],BT=avgs$BT[ii]),
	ylab="Prevalence",main=c(avgs$index[ii],YY=
	as.character(avgs$Year[ii])))
}
dev.off()
pdf("predictions_bairdi_ST.pdf",width=8, height=3.5)
for(ii in 1:9){
visreg(mod2,"ST",by=c("Size"),type="conditional",scale="response",
	print.cond=T,cond=list(II=avgs$index[ii],YY=avgs$Year[ii],
	BD=avgs$BD[ii],BT=avgs$BT[ii]),
	ylab="Prevalence",main=c(avgs$index[ii],YY=
	as.character(avgs$Year[ii])))
}
dev.off()
pdf("predictions_bairdi_BD.pdf",width=8, height=3.5)
for(ii in 1:9){
visreg(mod2,"BD",by=c("Size"),type="conditional",scale="response",
	print.cond=T,cond=list(II=avgs$index[ii],YY=avgs$Year[ii],
	ST=avgs$ST[ii],BT=avgs$BT[ii]),
	ylab="Prevalence",main=c(avgs$index[ii],YY=
	as.character(avgs$Year[ii])))
}
dev.off()
pdf("predictions_bairdi_BT.pdf",width=8, height=3.5)
for(ii in 1:9){
visreg(mod2,"BT",by=c("Size"),type="conditional",scale="response",
	print.cond=T,cond=list(II=avgs$index[ii],YY=avgs$Year[ii],
	BD=avgs$BD[ii],ST=avgs$ST[ii]),
	ylab="Prevalence",main=c(avgs$index[ii],YY=
	as.character(avgs$Year[ii])))
}
dev.off()
pdf("predictions_bairdi_data.pdf",width=7, height=7)
plot(bairdi1$BD~bairdi1$Size,col=bairdi1$RR+1,
	xlab="Size", ylab="Depth")
dev.off()





visreg(mod2,"Size", by=c("II"),type="conditional",scale="response",
	print.cond=T,cond=list(),
	ylab="Prevalence")
visreg(mod2,"Size", by=c("YY"),type="conditional",scale="response",
	ylab="Prevalence")

visreg(mod2,"YY", by=c("II"),type="conditional",scale="response",
	print.cond=T,cond=list(),
	ylab="Prevalence",xlab="Year")

visreg(mod2,"II", by=c("YY"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Index site")
visreg(mod2,"BT", by=c("II"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom temperature")
visreg(mod2,"BT", by=c("YY"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom temperature")
visreg(mod2,"BD", by=c("II"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom depth")
visreg(mod2,"BD", by=c("YY"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom depth")
visreg(mod2,"ST", by=c("II"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Surface_temperature")
visreg(mod2,"ST", by=c("YY"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Surface_temperature")
visreg(mod2,"BD", by=c("Size"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom depth")
visreg(mod2,"Size", by=c("BD"),type="conditional",scale="response",
	ylab="Prevalence")
visreg(mod2,"BD", by=c("ST"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom depth")
visreg(mod2,"BT", by=c("ST"),type="conditional",scale="response",
	ylab="Prevalence")
visreg2d(mod2,"ST","BD",type="conditional",scale="response",
	ylab="Bottom depth",main="Prevalence")
visreg2d(mod2,"BT","ST",type="conditional",scale="response",
	xlab="Bottom temperature",main="Prevalence")
visreg2d(mod2,"ST","BD",type="conditional",scale="response",
	ylab="Bottom depth",xlab="Surface temperature",main="Prevalence")
visreg2d(mod2,"ST","BD",type="effect",scale="response",
	ylab="Bottom depth",xlab="Surface temperature",main="Prevalence")
visreg2d(mod2,"Size","BD",type="effect",scale="response",
	ylab="Bottom depth",xlab="Size",main="Prevalence")
plot(bairdi1$BD~bairdi1$Size,col=bairdi1$RR+1)
plot(bairdi1$BD~bairdi1$Size,pch=bairdi1$RR+1)
plot(bairdi1$ST~bairdi1$Size,col=bairdi1$RR+1)
plot(bairdi1$BT~bairdi1$Size,col=bairdi1$RR+1)
plot(bairdi1$BT~bairdi1$ST,col=bairdi1$RR+1)
plot(bairdi1$BT~bairdi1$BD,col=bairdi1$RR+1)
plot(bairdi1$ST~bairdi1$BD,col=bairdi1$RR+1)
#not very usefull
#visreg(mod2,"Size", by="index",type="contrast",scale="response")


#Goodness of fit assesments:
anova(mod2)
plot(mod2) #first plot here is not correct
resids=residuals(mod2)
fits=fitted(mod2)
plot(resids~fits)
lines(smooth.spline(fits, resids, df = 5), lty = 2, col = "red")

#plot(mod2,toPdf = T)
#dx(mod2)
#gof(mod2, g = 10, plotROC = TRUE)
#Hosmer–Lemeshow GOF test
hoslem.test(mod2$y, mod2$fitted)







#OPILIO
#=======================================================================
#=======================================================================
#final glm model for opilio
mod2 <- glm(RR ~ 1 + Size + BD + BT:Size + II:Size + YY:Size,
	 data = opilio1, family = binomial)
summary(mod2)
plot(mod2)

#ordering variables from most to least important
anova(mod2)

#checking for overdispersion
mod2int <- glm(RR ~ 1 + Size + BD + BT:Size + II:Size + YY:Size,
	 data = opilio1, family = quasibinomial)
summary(mod2int)
# no need to model overdispersion as the dispersion is = 1
# Dispersion parameter for quasibinomial family taken to be 0.93


# Instead we can also test for overdispersion
pchisq(summary(mod2)$dispersion * mod2$df.residual, 
       mod2$df.residual, lower = F)
#p-value above 0.05 indicate no overdispersion
resids=residuals(mod2)
preds=fitted(mod2)
plot(resids~preds)
lines(smooth.spline(preds, resids, df = 10), lty = 2, col = "red")
qqnorm(resids)

#ploting predictor data
pdf("predictor_data_opilio.pdf",w=10,h=10)
ggpairs(opilio1[, c("RR","YY","II","Size","BD","BT")])
dev.off()
head(opilio1)

#predictors contriubtion
str(opilio)
IIsum=unique(opilio[,c(7,17,32,34:36,38)])
str(IIsum)
avgs=aggregate(cbind(BD=BOTTOM_DEPTH,ST=SURFACE_TEMP,
	Size=Size,BT=Bottom_Temp)~Year+index,data=IIsum,mean)
avgs
summary(mod2)
anova(mod2)
print(mod2)
confint(mod2) 
predict(mod2)
 
termplot(mod2,se=T,ylim=c(-5,5))

#pdf("predictions_opilio.pdf")
#visreg(mod2,print.cond=T,scale="response",
	ylab="Prevalence")
pdf("predictions_opilio_yy_ii.pdf",width=8, height=6)
boxplot(predict(mod2,newdata=opilio1,type="response")~mod2$data$II+mod2$data$YY)

dev.off()
pdf("predictions_opilio_size.pdf",width=8, height=8)
par(mfrow=c(3,3))
for(ii in 1:9){
visreg(mod2,"Size",type="conditional",scale="response",
	print.cond=T,cond=list(II=avgs$index[ii],YY=avgs$Year[ii],
	BD=avgs$BD[ii],BT=avgs$BT[ii]),
	ylab="Prevalence",main=c(avgs$index[ii],YY=
	as.character(avgs$Year[ii])))
}
dev.off()
pdf("predictions_opilio_BD.pdf",width=8, height=3.5)
for(ii in 1:9){
visreg(mod2,"BD",by=c("Size"),type="conditional",scale="response",
	print.cond=T,cond=list(II=avgs$index[ii],YY=avgs$Year[ii],
	BT=avgs$BT[ii]),
	ylab="Prevalence",main=c(avgs$index[ii],YY=
	as.character(avgs$Year[ii])))
}
dev.off()
pdf("predictions_opilio_BT.pdf",width=8, height=3.5)
for(ii in 1:9){
visreg(mod2,"BT",by=c("Size"),type="conditional",scale="response",
	print.cond=T,cond=list(II=avgs$index[ii],YY=avgs$Year[ii],
	BD=avgs$BD[ii]),
	ylab="Prevalence",main=c(avgs$index[ii],YY=
	as.character(avgs$Year[ii])))
}
dev.off()
pdf("predictions_opilio_data.pdf",width=7, height=7)
plot(opilio1$BD~opilio1$Size,col=opilio1$RR+1,
	xlab="Size", ylab="Depth")
dev.off()




#visreg(mod2,print.cond=T,scale="response",
	ylab="Prevalence", cond=list(YY="2016"), rug=2, partial=F)
#visreg(mod2,print.cond=T,scale="response",type="contrast",
	ylab="Prevalence")
visreg(mod2,"BD",type="conditional",scale="response",
	print.cond=T,cond=list()
	ylab="Prevalence",xlab="Bottom Depth")
visreg(mod2,"BD", by=c("II"),type="conditional",scale="response",
	print.cond=T,
	ylab="Prevalence",xlab="Bottom Depth")
visreg(mod2,"BD", by=c("YY"),type="conditional",scale="response",
	print.cond=T,
	ylab="Prevalence",xlab="Bottom Depth")
visreg(mod2,"Size", by=c("YY"),type="conditional",scale="response",
	print.cond=T,
	ylab="Prevalence",xlab="Bottom Depth")
visreg(mod2,"Size", by=c("YY"),type="conditional",scale="response",
	print.cond=T,cond=list(II=4),
	ylab="Prevalence",xlab="Bottom Depth")
visreg(mod2,"Size", by=c("YY"),type="conditional",scale="response",
	print.cond=T,cond=list(II=6),
	ylab="Prevalence",xlab="Bottom Depth")



visreg(mod2,"YY", by=c("II"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Year")
visreg(mod2,"II", by=c("YY"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Index site")
visreg(mod2,"Size", by=c("II"),type="conditional",scale="response",
	ylab="Prevalence")
visreg(mod2,"Size", by=c("YY"),type="conditional",scale="response",
	ylab="Prevalence")
visreg(mod2,"BT", by=c("II"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom temperature")
visreg(mod2,"BT", by=c("YY"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom temperature")
visreg(mod2,"BD", by=c("II"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom depth")
visreg(mod2,"BD", by=c("YY"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom depth")
visreg(mod2,"BD", by=c("Size"),type="conditional",scale="response",
	ylab="Prevalence",xlab="Bottom depth")
visreg(mod2,"Size", by=c("BD"),type="conditional",scale="response",
	ylab="Prevalence")
plot(opilio1$BD~opilio1$Size,col=opilio1$RR+1)
plot(opilio1$BD~opilio1$Size,pch=opilio1$RR+46)
plot(opilio1$BT~opilio1$Size,col=opilio1$RR+1)
visreg2d(mod2,"YY","II",type="conditional",scale="response")
dev.off()

#not very usefull
#visreg(mod2,"Size", by="II",type="contrast",scale="response")


#Goodness of fit assesments:
anova(mod2)
plot(mod2) #first plot here is not correct
resids=residuals(mod2)
fits=fitted(mod2)
plot(resids~fits)
lines(smooth.spline(fits, resids, df = 5), lty = 2, col = "red")

#plot(mod2,toPdf = T)
#dx(mod2)
#gof(mod2, g = 10, plotROC = TRUE)
#Hosmer–Lemeshow GOF test
hoslem.test(mod2$y, mod2$fitted)






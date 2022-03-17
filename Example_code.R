###############################################################################
#
#
# Title: Fluctuating temperature modifies heat-mortality association in the globe 
# 
# Author: Yao Wu, Bo Wen, Yuming Guo, Shanshan Li & co-authors from the MCC network
#
# Affiliation: Monash University
#
# Example code
#
#
###############################################################################
library(dlnm);library(splines);library(gnm);library(tsModel);library(mgcv);library(mixmeta);
library(ggplot2);library(ggpubr);library(ggsci);library(dplyr)


## LOAD DATA
df.fr <- read.csv("Example_data.csv")

department_list <- unique(df.fr[,c("citycode","cityname")])
names(department_list) <- "city"

###############################################################################
##### FIRST STAGE #####
# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "ns"
vardegree <- NULL
varper <- c(50, 90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 10
lagnk <- 2

# DEGREE OF FREEDOM FOR SEASONALITY
dfseas <- 4

# SET MATRIX TO STORE RESULTS
coefm1<-matrix(data = NA, nrow = nrow(department_list), ncol = (length(varper)+1),
               dimnames = list(department_list$city,paste0("est",1:(length(varper)+1))))
coefm2<-matrix(data = NA, nrow = nrow(department_list), ncol = (length(varper)+1),
               dimnames = list(department_list$city,paste0("est",1:(length(varper)+1))))
coefm3<-matrix(data = NA, nrow = nrow(department_list), ncol = (length(varper)+1),
               dimnames = list(department_list$city,paste0("est",1:(length(varper)+1))))
coefm4<-matrix(data = NA, nrow = nrow(department_list), ncol = (length(varper)+1),
               dimnames = list(department_list$city,paste0("est",1:(length(varper)+1))))

vcovl1 <- vector("list",nrow(department_list))
names(vcovl1) <- department_list$city
vcovl2 <- vector("list",nrow(department_list))
names(vcovl2) <- department_list$city
vcovl3 <- vector("list",nrow(department_list))
names(vcovl3) <- department_list$city
vcovl4 <- vector("list",nrow(department_list))
names(vcovl4) <- department_list$city

# LOOP TO MODEL
for (i in seq(department_list$city)){
  cat(i,"\n ")
  data <- df.fr[df.fr$citycode==department_list$city[i],]
  
  ## TV group
  tv<-apply(with(data,cbind(Lag(tmin,0:1),Lag(tmax,0:1))),1,sd,na.rm=T)
  qk<-quantile(tv,na.rm=T)
  data$qgroup[tv<=qk[2]]<-"q1"
  data$qgroup[tv>qk[2]&tv<=qk[3]]<-"q2"
  data$qgroup[tv>qk[3]&tv<=qk[4]]<-"q3"
  data$qgroup[tv>qk[4]&tv<=qk[5]]<-"q4"
  
  data$date <- as.Date(data$date)

  data$year <- as.numeric(format(data$date,"%Y"))
  data$dow <- format(data$date,"%w")
  data$time<-1:nrow(data) # Create time
  df<-(max(data$year)-min(data$year)+1)*dfseas
  
  ## modelling
  bs.t<-crossbasis(data$tmean, lag=lag,
                   argvar=list(fun=varfun,
                               knots = quantile(data$tmean, varper / 100, na.rm = T),
                               Bound = range(data$tmean, na.rm = T)),
                   arglag=list(knots=logknots(lag, lagnk)))
  
  cen <- quantile(data$tmean,probs = 0.75,na.rm = T)
  
  fit<-glm(y~bs.t:qgroup+tv+ns(time,df=df)+as.factor(dow),family=quasipoisson,data=data)
  summary(fit)
  ## group 1
  indexcoef1<-grep("qgroupq1",names(coef(fit)))
  coef1<-coef(fit)[indexcoef1]
  
  indexvcov1<-grep("qgroupq1",rownames(vcov(fit)))
  vcov1<-vcov(fit)[indexvcov1,indexvcov1]
  
  pred1<-crossreduce(bs.t, coef=coef1, vcov=vcov1, type="overall",model.link="log",cen=cen)
  coefm1[i,]<-coef(pred1)
  vcovl1[[i]]<-vcov(pred1)
  
  
  ## group 2
  indexcoef2<-grep("qgroupq2",names(coef(fit)))
  coef2<-coef(fit)[indexcoef2]
  
  indexvcov2<-grep("qgroupq2",rownames(vcov(fit)))
  vcov2<-vcov(fit)[indexvcov2,indexvcov2]
  
  pred2<-crossreduce(bs.t, coef=coef2, vcov=vcov2, type="overall",model.link="log",cen=cen)
  coefm2[i,]<-coef(pred2)
  vcovl2[[i]]<-vcov(pred2)
  
  
  
  ## group 3
  indexcoef3<-grep("qgroupq3",names(coef(fit)))
  coef3<-coef(fit)[indexcoef3]
  
  indexvcov3<-grep("qgroupq3",rownames(vcov(fit)))
  vcov3<-vcov(fit)[indexvcov3,indexvcov3]
  
  pred3<-crossreduce(bs.t, coef=coef3, vcov=vcov3, type="overall",model.link="log",cen=cen)
  coefm3[i,]<-coef(pred3)
  vcovl3[[i]]<-vcov(pred3)
  
  
  
  ## group 4
  indexcoef4<-grep("qgroupq4",names(coef(fit)))
  coef4<-coef(fit)[indexcoef4]
  
  indexvcov4<-grep("qgroupq4",rownames(vcov(fit)))
  vcov4<-vcov(fit)[indexvcov4,indexvcov4]
  
  pred4<-crossreduce(bs.t, coef=coef4, vcov=vcov4, type="overall",model.link="log",cen=cen)
  coefm4[i,]<-coef(pred4)
  vcovl4[[i]]<-vcov(pred4)
}

###############################################################################
##### Second stage #####
country.name <- department_list$city

### META FOR NATIONAL RESULT
meta_nest1<-mixmeta(coefm1~1,vcovl1)
meta_nest2<-mixmeta(coefm2~1,vcovl2)
meta_nest3<-mixmeta(coefm3~1,vcovl3)
meta_nest4<-mixmeta(coefm4~1,vcovl4)



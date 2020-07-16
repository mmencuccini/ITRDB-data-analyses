setwd("D:/Dropbox/My Documents/science/projects/grants/China Tong Zheng/Tong Zheng folder/R script 20200621")

# start from R.Data in this folder to load Tong's quality controlled Workspace, then save under MM.RData
# after subsetting for Li data, if needed (i.e., start from line 1102 directly)


#some special packages
# install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
# install.packages("devtools")
# library(devtools)
# devtools::install_github("valentinitnelav/plotbiomes")
# library(plotbiomes)


####【Library】####
library(coin)
library(lmerTest) # this gives lsmeans; otherwise use emmeans package with more options
library(optimx)
library(tidyverse)
library('gridExtra')
library("mgcViz")
library(gamm4)
library(data.table) # loads data.table and set functions from library
library(performance)
library(SPEI)
library(raster)
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
library(ncdf.tools)
library(PCICt)
library(ncdf4.helpers)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggmap)
library(maps)
library(cowplot)
library(pointRes)
library(lme4)
library(nlme)
library(dplR)
library(MuMIn)
library(spam)
library(grid)
library(effects)
library(Taxonstand)
library(stringr)
library(stringi)
library(taxonlookup)
library(interplot)
library(visreg)
library(splines)
library(lattice)
library(latticeExtra)
library(car)
library(scales)
library(rlang)
library(Rmisc)
library(plotrix)
library(devtools)
library(plotbiomes)
library(sjPlot)
library(RColorBrewer)
library(emmeans)


meanNA <- function(x) mean(x, na.rm=TRUE)

medianNA <- function(x) median(x, na.rm=TRUE)


# is.NaN method for data frames
is.nan.data.frame <- function(x)
{
  do.call(cbind, lapply(x, is.nan))
} 

# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
lsos()





####【data backup 20200614】####
load(file = "INFOSET_T.r") #the largest dataset, including the information of the years from 1901 to 2010

####【Center the factors, preparing to set the model 】####
#If you want to set the model, do this first.
INFOSET_T$log_RWI4 <- with(INFOSET_T, ave(log_RWI4, site, FUN=function(x) x - mean(x,na.rm=TRUE))) #Log Mean_growth_4y should be group centered
INFOSET_T$year <- with(INFOSET_T, year - mean(year, na.rm=TRUE)) #year should be Global centered
INFOSET_T$TMP_Annual <- with(INFOSET_T, TMP_Annual - mean(TMP_Annual, na.rm=TRUE)) #TMP_Annual should be Global centered
INFOSET_T$TMP_Annual_fixed <- with(INFOSET_T, TMP_Annual_fixed - mean(TMP_Annual_fixed, na.rm=TRUE)) #TMP_Annual_fixed should be Global centered
INFOSET_T$PRE_Annual <- with(INFOSET_T, PRE_Annual - mean(PRE_Annual, na.rm=TRUE)) #PRE_Annual should be Global centered
INFOSET_T$PRE_Annual_fixed <- with(INFOSET_T, PRE_Annual_fixed - mean(PRE_Annual_fixed, na.rm=TRUE)) #PRE_Annual_fixed should be Global centered
INFOSET_T$Max_age <- with(INFOSET_T, Max_age - mean(Max_age, na.rm=TRUE)) #Max_age should be Global centered
INFOSET_T$continentality <- with(INFOSET_T, continentality - mean(continentality, na.rm=TRUE)) #continentality should be Global centered
INFOSET_T$continentality_fixed <- with(INFOSET_T, continentality_fixed - mean(continentality_fixed, na.rm=TRUE)) 
INFOSET_T$log_Max_age <- with(INFOSET_T, log_Max_age - mean(log_Max_age, na.rm=TRUE)) #log_Max_age should be Global centered
INFOSET_T$log_PRE_Annual_fixed <- with(INFOSET_T, log_PRE_Annual_fixed - mean(log_PRE_Annual_fixed, na.rm=TRUE)) 



# ####【r files backup】####
# load(file = paste(getwd(),"/r.save/name.infoset.r", sep = ""))
# load(file = paste(getwd(),"/r.save/siteinfo_Taxonomy.r", sep = ""))
# load(file = paste(getwd(),"/r.save/sumspei.r", sep = ""))
# load(file = paste(getwd(),"/r.save/sumspei3.r", sep = ""))
# load(file = paste(getwd(),"/r.save/sumpre.r", sep = ""))
# load(file = paste(getwd(),"/r.save/sumtmp.r", sep = ""))
# load(file = paste(getwd(),"/r.save/sumtmx.r", sep = ""))
# load(file = paste(getwd(),"/r.save/sumtmn.r", sep = ""))
# load(file = paste(getwd(),"/r.save/Species names.r", sep = ""))
# load(file = paste(getwd(),"/r.save/coord.r", sep = ""))












########################【【【【【Involving the SPEIprev SPEIpost and SPEId into Model II, 20200621】】】】】##########################
#this means we need to calculate SPEI 4 years before/after the drought as well
#But I'd like to calculate RWIprev again for double-check before doing it to SPEI
#First read the RWI file
seriesoriginal <- read.table(file = "series_without0_fixed.txt", header = T)

series <- data.frame()
i <- 1
for (i in 1:length(unique(INFOSET_T$site))) {
  value1 <- as.character(unique(INFOSET_T$site))[i] #site name
  value2 <- which(colnames(seriesoriginal) == value1) #which column is this site
  value3 <- seriesoriginal[1:118,value2] #the subset of this site, RWI
  value4 <- as.data.frame(cbind(rep(value1, 118), 1897:2014, value3)) #c(site, year, RWI)
  j <- 5
  value8 <- data.frame()
  value9 <- data.frame()
  for (j in 5:114) {
    value5 <- mean(as.numeric(value4$value3[(j-4):(j-1)]), na.action = na.omit) #RWIprev
    value6 <- value4$value3[j] #RWId
    value7 <- mean(as.numeric(value4$value3[(j+1):(j+4)]), na.action = na.omit) #RWIpost
    value8 <- rbind(value8, c(1896+j, value5, value6, value7)) #c(year, RWIprev, RWId, RWIpost)
  }
  value9 <- cbind(value4[-c(1:4,115:118),], value8)
  colnames(value9) <- c("site", "year", "RWI", "year2", "RWI4prev", "RWI2", "RWI4post")
  series <- rbind(series, value9)
  print(i)
}
series <- series[,-c(4,6)]
colnames(series) <- c("site", "year", "RWI", "RWI4prev", "RWI4post")
rownames(series) <- 1:nrow(series)



#Time to do it to SPEI, i.e. SPEIprev, SPEId and SPEIpost
#Read the SPEI list
SPEId_list <- sumspei[[3]]
colnames(SPEId_list) <- c("site", "Latitude", 1901:2015)

SPEIprevpost <- data.frame()
i <- 1
for (i in 1:length(unique(INFOSET_T$site))) {
  value1 <- as.character(unique(INFOSET_T$site))[i] #site name
  value2 <- which(SPEId_list$site == value1) #which row is this site
  value3 <- t(SPEId_list[value2, 3:117]) #the subset of this site, SPEI
  value4 <- as.data.frame(cbind(rep(value1, 115), 1901:2015, value3)) #c(site, year, SPEI)
  j <- 1
  value8 <- data.frame()
  value9 <- data.frame()
  for (j in 1:115) { #j is the sequence of year
    if (j >= 5) {value5 <- mean(as.numeric(value4[(j-4):(j-1),3]))} #SPEIprev
    if (j == 5) {value5 <- mean(as.numeric(value4[2:4,3]))}
    if (j == 4) {value5 <- mean(as.numeric(value4[2:3,3]))}
    if (j == 3) {value5 <- mean(as.numeric(value4[2,3]))}
    if (j == 2) {value5 <- NA}
    if (j == 1) {value5 <- NA}
    value6 <- value4[j,3] #SPEId
    value7 <- mean(as.numeric(value4[(j+1):(j+4),3]), na.action = na.omit) #SPEIpost
    value8 <- rbind(value8, c(1900+j, value5, value6, value7)) #c(year, SPEIprev, SPEId, SPEIpost)
  }
  value9 <- cbind(value4[-c(111:115),], value8[1:110,])
  colnames(value9) <- c("site", "year", "SPEI", "year2", "SPEIprev", "SPEI2", "SPEIpost")
  SPEIprevpost <- rbind(SPEIprevpost, value9)
  print(i)
}



SPEIprevpost <- SPEIprevpost[,-c(4,6)]
colnames(SPEIprevpost) <- c("site", "year", "SPEI", "SPEIprev", "SPEIpost")
rownames(SPEIprevpost) <- 1:nrow(SPEIprevpost)



# #some changes in INFOSET_T
INFOSET_T <- cbind(INFOSET_T[,1:5], series[,4:5], INFOSET_T[,7:39])
INFOSET_T <- INFOSET_T[,-(18:19)]
INFOSET_T$RWI4prev <- as.numeric(INFOSET_T$RWI4prev)
INFOSET_T$RWI4post <- as.numeric(INFOSET_T$RWI4post)



INFOSET_T <- cbind(INFOSET_T[,1:18], SPEIprevpost[,4:5], INFOSET_T[,19:38])
colnames(INFOSET_T)[18] <- "SPEId"
INFOSET_T$SPEIprev <- as.numeric(INFOSET_T$SPEIprev)
INFOSET_T$SPEIpost <- as.numeric(INFOSET_T$SPEIpost)


save(series, file = "series.r")
save(INFOSET_T, file = "INFOSET_T.r")


#Subsets of INFOSET_T
INFOSET <- subset(INFOSET_T, is.na(TRI) == FALSE) #all the years, with tree rings
db1.5_T <- subset(INFOSET_T, SPEId + 1.5 < 0) #SPEI < -1.5 drought years, with and without tree rings
db1.5 <- subset(INFOSET, SPEId + 1.5 < 0) #SPEI < -1.5 drought years, with tree rings








#####################【【【【【Main LME model, 20200621】】】】】#####################
########################【【【【【Build Model II again, 20200621】】】】】##########################
#### Updated on 20200603, we have RWI4 log-transformed in the model
#### Updated on 20200603, we have the interaction of Groups*log_RWI4 in the models mm4 to mm6 now
#### Updated on 20200614, we have correct these wrong data in SPEI3m_Min and SPEI3m_Min_Season, so let's run the main model again! 
load(file = "INFOSET_T.r") #the largest dataset, including the information of the years from 1901 to 2010

INFOSET_T$RWI4prev <- with(INFOSET_T, ave(RWI4prev, site, FUN=function(x) x - mean(x,na.rm=TRUE)))
INFOSET_T$RWI4post <- with(INFOSET_T, ave(RWI4post, site, FUN=function(x) x - mean(x,na.rm=TRUE)))
INFOSET_T$year <- with(INFOSET_T, year - mean(year, na.rm=TRUE))
INFOSET_T$TMP_Annual <- with(INFOSET_T, TMP_Annual - mean(TMP_Annual, na.rm=TRUE))
INFOSET_T$TMP_Annual_fixed <- with(INFOSET_T, TMP_Annual_fixed - mean(TMP_Annual_fixed, na.rm=TRUE))
INFOSET_T$PRE_Annual <- with(INFOSET_T, PRE_Annual - mean(PRE_Annual, na.rm=TRUE))
INFOSET_T$PRE_Annual_fixed <- with(INFOSET_T, PRE_Annual_fixed - mean(PRE_Annual_fixed, na.rm=TRUE))
INFOSET_T$Max_age <- with(INFOSET_T, Max_age - mean(Max_age, na.rm=TRUE))
INFOSET_T$continentality <- with(INFOSET_T, continentality - mean(continentality, na.rm=TRUE))
INFOSET_T$continentality_fixed <- with(INFOSET_T, continentality_fixed - mean(continentality_fixed, na.rm=TRUE)) 
INFOSET_T$log_Max_age <- with(INFOSET_T, log_Max_age - mean(log_Max_age, na.rm=TRUE))
INFOSET_T$log_PRE_Annual_fixed <- with(INFOSET_T, log_PRE_Annual_fixed - mean(log_PRE_Annual_fixed, na.rm=TRUE)) 

INFOSET <- subset(INFOSET_T, is.na(TRI) == FALSE)
db1.5_T <- subset(INFOSET_T, SPEId + 1.5 < 0)
db1.5 <- subset(INFOSET, SPEId + 1.5 < 0)
db1.0 <- subset(INFOSET, SPEId + 1.0 < 0)


#part 0: the model with only Year * Group
mm0 <- lme(log(Resistance) ~ year + Groups + year* Groups, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.0, control = list(opt = "optim"), method = "REML")
mm00 <- lme(log(Recovery) ~ year + Groups + year* Groups, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.0, control = list(opt = "optim"), method = "REML")
mm000 <- lme(log(Resilience) ~ year + Groups + year* Groups, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.0, control = list(opt = "optim"), method = "REML")

#part 1: the model without CUM4SPEI, SPEId and Season
mm1 <- lme(log(Resistance) ~ year + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + Groups +  
             year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + year* log_PRE_Annual_fixed + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.0, control = list(opt = "optim"), method = "REML")
mm2 <- lme(log(Recovery) ~ year + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + Groups + 
             year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + year* log_PRE_Annual_fixed + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.0, control = list(opt = "optim"), method = "REML")
mm3 <- lme(log(Resilience) ~ year + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + Groups +  
             year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + year* log_PRE_Annual_fixed + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.0, control = list(opt = "optim"), method = "REML")

#part 2: the model with CUM4SPEI, SPEId, Season and interactions
mm4 <- lme(log(Resistance) ~ year + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
             Groups + slope_TMPxYear + slope_PRExYear + year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + 
             year* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.0, control = list(opt = "optim"), method = "REML")
mm5 <- lme(log(Recovery) ~   year + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
             Groups + slope_TMPxYear + slope_PRExYear + year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + 
             year* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.0, control = list(opt = "optim"), method = "REML")
mm6 <- lme(log(Resilience) ~ year + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
             Groups + slope_TMPxYear + slope_PRExYear + year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + 
             year* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.0, control = list(opt = "optim"), method = "REML")

anova(mm0)
anova(mm00)
anova(mm000)
anova(mm1)
anova(mm2)
anova(mm3)
anova(mm4)
anova(mm5)
anova(mm6)

summary(mm0)
summary(mm00)
summary(mm000)
summary(mm1)
summary(mm2)
summary(mm3)
summary(mm4)
summary(mm5)
summary(mm6)

r.squaredGLMM(mm0)
r.squaredGLMM(mm00)
r.squaredGLMM(mm000)
r.squaredGLMM(mm1)
r.squaredGLMM(mm2)
r.squaredGLMM(mm3)
r.squaredGLMM(mm4)
r.squaredGLMM(mm5)
r.squaredGLMM(mm6)












#####################【【【【【Figure 1 Standardized effects sizes, 20200621】】】】】#####################
#Figure1
#version with Y-axis
Figure0 <- plot_model(mm0, show.data = T, type = "std", show.values = F, show.p = F, sort.est = F, title = "log(Rt)", axis.title = "Standardized effect size", 
                       vline.color = "grey", axis.labels = c("Year:Groups", "Groups", "Year")) + theme_bw() + ylim(-0.2,0.2)

Figure00 <- plot_model(mm00, show.data = T, type = "std", show.values = F, show.p = F, sort.est = F, title = "log(Rc)", axis.title = "Standardized effect size", 
                       vline.color = "grey", axis.labels = c("Year:Groups", "Groups", "Year")) + theme_bw() + ylim(-0.2,0.2)

Figure000 <- plot_model(mm000, show.data = T, type = "std", show.values = F, show.p = F, sort.est = F, title = "log(Rs)", axis.title = "Standardized effect size", 
                        vline.color = "grey", axis.labels = c("Year:Groups", "Groups", "Year")) + theme_bw() + ylim(-0.2,0.2)

Figure11 <- plot_model(mm1, show.data = T, type = "std", show.values = F, show.p = F, sort.est = F, title = "log(Rt)", axis.title = "Standardized effect size", 
                      vline.color = "grey", axis.labels = c("Groups:RWI4", "Year:TMP", "Year:PRE", "Year:CL", "Year:Groups", "Year:K", "Groups", "K", 
                                                            "PRE", "TMP", "CL", "RWI4", "Year")) + theme_bw() + ylim(-0.2,0.2)
Figure12 <- plot_model(mm2, show.data = T, type = "std", show.values = F, show.p = F, sort.est = F, title = "log(Rc)", axis.title = "Standardized effect size", 
                      vline.color = "grey", axis.labels = c("Groups:RWI4", "Year:TMP", "Year:PRE", "Year:CL", "Year:Groups", "Year:K", "Groups", "K", 
                                                            "PRE", "TMP", "CL", "RWI4", "Year")) + theme_bw() + ylim(-0.2,0.2)
Figure13 <- plot_model(mm3, show.data = T, type = "std", show.values = F, show.p = F, sort.est = F, title = "log(Rs)", axis.title = "Standardized effect size", 
                      vline.color = "grey", axis.labels = c("Groups:RWI4", "Year:TMP", "Year:PRE", "Year:CL", "Year:Groups", "Year:K", "Groups", "K", 
                                                            "PRE", "TMP", "CL", "RWI4", "Year")) + theme_bw() + ylim(-0.2,0.2)

Figure14 <- plot_model(mm4, show.data = T, type = "std", show.values = F, show.p = F, sort.est = F, title = "log(Rt)", axis.title = "Standardized effect size", 
                      vline.color = "grey", axis.labels = c("Groups:RWI4", "Groups:SPEIpost", "Groups:SPEIprev", "Groups:SPEId", "Year:PRE", "Year:TMP", "Year:CL", "Year:Groups", 
                                                            "Year:K", "SPRE", "STMP", "Groups", "K", "PRE", "TMP", "CL", "SeasonWinter", "SeasonSummer", "SeasonSpring", 
                                                            "SPEIpost", "SPEIprev", "SPEId","RWI4", "Year")) + theme_bw() + ylim(-0.2,0.2)

Figure15 <- plot_model(mm5, show.data = T, type = "std", show.values = F, show.p = F, sort.est = F, title = "log(Rc)", axis.title = "Standardized effect size", 
                       vline.color = "grey", axis.labels = c("Groups:RWI4", "Groups:SPEIpost", "Groups:SPEIprev", "Groups:SPEId", "Year:PRE", "Year:TMP", "Year:CL", "Year:Groups", 
                                                             "Year:K", "SPRE", "STMP", "Groups", "K", "PRE", "TMP", "CL", "SeasonWinter", "SeasonSummer", "SeasonSpring", 
                                                             "SPEIpost", "SPEIprev", "SPEId","RWI4", "Year")) + theme_bw() + ylim(-0.2,0.2)

Figure16 <- plot_model(mm6, show.data = T, type = "std", show.values = F, show.p = F, sort.est = F, title = "log(Rs)", axis.title = "Standardized effect size", 
                       vline.color = "grey", axis.labels = c("Groups:RWI4", "Groups:SPEIpost", "Groups:SPEIprev", "Groups:SPEId", "Year:PRE", "Year:TMP", "Year:CL", "Year:Groups", 
                                                             "Year:K", "SPRE", "STMP", "Groups", "K", "PRE", "TMP", "CL", "SeasonWinter", "SeasonSummer", "SeasonSpring", 
                                                             "SPEIpost", "SPEIprev", "SPEId", "RWI4", "Year")) + theme_bw() + ylim(-0.2,0.2)

#Export these 6 figures together
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,3))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(Figure0, vp = vplayout(1,1))   
print(Figure00, vp = vplayout(1,2))   
print(Figure000, vp = vplayout(1,3))   
# print(Figure11, vp = vplayout(1,1))   
# print(Figure12, vp = vplayout(1,2))   
# print(Figure13, vp = vplayout(1,3))  
print(Figure14, vp = vplayout(2,1))   
print(Figure15, vp = vplayout(2,2))   
print(Figure16, vp = vplayout(2,3)) 



















###########【【【【【Figure 2, Multi-slope plot for Resilience indices, 20200623 】】】】】############
###########【【【【【Calculate the Multi-slopes for Resilience indices, using the model with SPEI_August, CUM4SPEI and Season】】】】】############
#multislope for log(Resistance)
multislope_Rt <- data.frame()
i <- 1
for(i in 1:nrow(coef(mm1))){
  value0 <- rownames(coef(mm1)[i,]) #the site
  value1 <- coef(mm1)[["year"]][i] #slope of year
  value2 <- coef(mm1)[["year:continentality_fixed"]][i] #slope of year:K
  value3 <- coef(mm1)[["year:GroupsAngiosperms"]][i] #slope of year:Groups
  value4 <- coef(mm1)[["year:log_Max_age"]][i] #slope of year:MA
  value5 <- coef(mm1)[["year:log_PRE_Annual_fixed"]][i] #slope of year:PRE
  value6 <- coef(mm1)[["year:TMP_Annual_fixed"]][i] #slope of year:TMP
  value7 <- subset(db1.5, site == value0)[1,]$continentality
  value8 <- as.character(subset(db1.5, site == value0)[1,]$Groups)
  value9 <- subset(db1.5, site == value0)[1,]$log_Max_age
  value10 <- subset(db1.5, site == value0)[1,]$log_PRE_Annual_fixed
  value11 <- subset(db1.5, site == value0)[1,]$TMP_Annual_fixed
  if(value8 == "Angiosperms"){
    value12 <- value1 + value2*value7 + value3 + value4*value9 + value5*value10 + value6*value11
  }else{
    value12 <- value1 + value2*value7 + value4*value9 + value5*value10 + value6*value11
  }
  multislope_Rt <- rbind(multislope_Rt, c(value0, value1, value12)) 
  print(i)
}
rownames(multislope_Rt) <- 1:nrow(multislope_Rt)
colnames(multislope_Rt) <- c("site", "yearslope_Rt", "multislope_Rt")
View(multislope_Rt)


#multislope for log(Recovery)
multislope_Rc <- data.frame()
i <- 1
for(i in 1:nrow(coef(mm2))){
  value0 <- rownames(coef(mm2)[i,]) #the site
  value1 <- coef(mm2)[["year"]][i] #slope of year
  value2 <- coef(mm2)[["year:continentality_fixed"]][i] #slope of year:K
  value3 <- coef(mm2)[["year:GroupsAngiosperms"]][i] #slope of year:Groups
  value4 <- coef(mm2)[["year:log_Max_age"]][i] #slope of year:MA
  value5 <- coef(mm2)[["year:log_PRE_Annual_fixed"]][i] #slope of year:PRE
  value6 <- coef(mm2)[["year:TMP_Annual_fixed"]][i] #slope of year:TMP
  value7 <- subset(db1.5, site == value0)[1,]$continentality
  value8 <- as.character(subset(db1.5, site == value0)[1,]$Groups)
  value9 <- subset(db1.5, site == value0)[1,]$log_Max_age
  value10 <- subset(db1.5, site == value0)[1,]$log_PRE_Annual_fixed
  value11 <- subset(db1.5, site == value0)[1,]$TMP_Annual_fixed
  if(value8 == "Angiosperms"){
    value12 <- value1 + value2*value7 + value3 + value4*value9 + value5*value10 + value6*value11
  }else{
    value12 <- value1 + value2*value7 + value4*value9 + value5*value10 + value6*value11
  }
  multislope_Rc <- rbind(multislope_Rc, c(value0, value1, value12)) 
  print(i)
}
rownames(multislope_Rc) <- 1:nrow(multislope_Rc)
colnames(multislope_Rc) <- c("site", "yearslope_Rc", "multislope_Rc")
View(multislope_Rc)


#multislope for log(Resilience)
multislope_Rs <- data.frame()
i <- 1
for(i in 1:nrow(coef(mm3))){
  value0 <- rownames(coef(mm3)[i,]) #the site
  value1 <- coef(mm3)[["year"]][i] #slope of year
  value2 <- coef(mm3)[["year:continentality_fixed"]][i] #slope of year:K
  value3 <- coef(mm3)[["year:GroupsAngiosperms"]][i] #slope of year:Groups
  value4 <- coef(mm3)[["year:log_Max_age"]][i] #slope of year:MA
  value5 <- coef(mm3)[["year:log_PRE_Annual_fixed"]][i] #slope of year:PRE
  value6 <- coef(mm3)[["year:TMP_Annual_fixed"]][i] #slope of year:TMP
  value7 <- subset(db1.5, site == value0)[1,]$continentality
  value8 <- as.character(subset(db1.5, site == value0)[1,]$Groups)
  value9 <- subset(db1.5, site == value0)[1,]$log_Max_age
  value10 <- subset(db1.5, site == value0)[1,]$log_PRE_Annual_fixed
  value11 <- subset(db1.5, site == value0)[1,]$TMP_Annual_fixed
  if(value8 == "Angiosperms"){
    value12 <- value1 + value2*value7 + value3 + value4*value9 + value5*value10 + value6*value11
  }else{
    value12 <- value1 + value2*value7 + value4*value9 + value5*value10 + value6*value11
  }
  multislope_Rs <- rbind(multislope_Rs, c(value0, value1, value12)) 
  print(i)
}
rownames(multislope_Rs) <- 1:nrow(multislope_Rs)
colnames(multislope_Rs) <- c("site", "yearslope_Rs", "multislope_Rs")
View(multislope_Rs)

#convert into dataframe
multislope_Rt <- as.data.frame(multislope_Rt)
multislope_Rc <- as.data.frame(multislope_Rc)
multislope_Rs <- as.data.frame(multislope_Rs)

#【【Important update on 20200620, multislope_Rt, multislope_Rc, multislope_Rs used to be character, we need to convert it to numeric, 】】
#【【though it would not affect the result of plotting】】
multislope_Rt$yearslope_Rt <- as.numeric(multislope_Rt$yearslope_Rt)
multislope_Rt$multislope_Rt <- as.numeric(multislope_Rt$multislope_Rt)
multislope_Rc$yearslope_Rc <- as.numeric(multislope_Rc$yearslope_Rc)
multislope_Rc$multislope_Rc <- as.numeric(multislope_Rc$multislope_Rc)
multislope_Rs$yearslope_Rs <- as.numeric(multislope_Rs$yearslope_Rs)
multislope_Rs$multislope_Rs <- as.numeric(multislope_Rs$multislope_Rs)

#put them together by site, because they don't share the same sites
multislope <- data.frame()
multislope <- left_join(multislope_Rt, multislope_Rc, by = "site")
multislope <- left_join(multislope, multislope_Rs, by = "site")
View(multislope)

#put multislope into coord then plot
# load(file = "coord.r")
coord <- coord[,-c(13:18)]
coord <- left_join(coord, multislope, by = "site")
write.table(coord, "coord with Multislope.csv", sep = ",", row.names = FALSE)
save(coord, file = "coord.r")



#start plotting
#Figure2a
Figure2a <- NULL
mapWorld <- borders("world", colour="grey", fill="white") # create a layer of borders
Figure2a <- ggplot() + mapWorld
#Now Layer the cities on top
Figure2a <- Figure2a + geom_point(aes(x=coord$x.real, y=coord$y.real, color=coord$multislope_Rt), size=1) + theme_bw() + 
  labs(title = "", x = "Longitude (°)", y = "Latitude (°)") + 
  scale_colour_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0, name = "slope") + theme(legend.position = c(0.1, 0.3)) + 
  ggtitle("Model without SPEI_August and CUM4SPEI \nMulti slope of log(Resistance) for drought years")
Figure2a


#Figure2b
Figure2b <- NULL
mapWorld <- borders("world", colour="grey", fill="white") # create a layer of borders
Figure2b <- ggplot() + mapWorld
#Now Layer the cities on top
Figure2b <- Figure2b + geom_point(aes(x=coord$x.real, y=coord$y.real, color=coord$multislope_Rc), size=1) + theme_bw() + 
  labs(title = "", x = "Longitude (°)", y = "Latitude (°)") + 
  scale_colour_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0, name = "slope") + theme(legend.position = c(0.1, 0.3)) + 
  ggtitle("Model without SPEI_August and CUM4SPEI \nMulti slope of log(Recovery) for drought years")
Figure2b


#Figure3
Figure3 <- NULL
mapWorld <- borders("world", colour="grey", fill="white") # create a layer of borders
Figure3 <- ggplot() + mapWorld
#Now Layer the cities on top
Figure3 <- Figure3 + geom_point(aes(x=coord$x.real, y=coord$y.real, color=coord$multislope_Rs), size=1) + theme_bw() + 
  labs(title = "", x = "Longitude (°)", y = "Latitude (°)") + 
  scale_colour_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0, name = "slope")  + 
  theme(legend.position = c(0.1, 0.3)) + ggtitle("Model without SPEI_August and CUM4SPEI \nMulti slope of log(Resilience) for drought years")
Figure3





















################【【【【【Figure 4 extra infomation, 2020622, the significance test for the two groups of Figure 4, using emtrend】】】】】#############
## Learn experience from demo
# fiber.lm <- lm(strength ~ diameter*machine, data = fiber)
# emtrends(fiber.lm, pairwise ~ machine, var = "diameter")
# ef <- effect("diameter:machine", fiber.lm)
# x <- as.data.frame(ef)
# x$Groups <- as.factor(as.character(x$machine))
# ggplot(data = x, aes(x=diameter, y=fit, color=machine)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
#   theme_classic() + theme(legend.position="none") + ggtitle(" \n") 

value1 <- emtrends(mm1, pairwise ~ Groups, var = "year")
value2 <- emtrends(mm2, pairwise ~ Groups, var = "year")
value3 <- emtrends(mm3, pairwise ~ Groups, var = "year")
value4 <- emtrends(mm4, pairwise ~ Groups, var = "year")
value5 <- emtrends(mm5, pairwise ~ Groups, var = "year")
value6 <- emtrends(mm6, pairwise ~ Groups, var = "year")

################【【【【【Figure 4, 2020622, Effect plot for log(Rx) against year in difference of Angiosperm and Gymnosperm in two models】】】】】#############
ef <- effect("year: Groups", mm1)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
Figure41 <- ggplot(data = x, aes(x=year, y=fit, color=Groups)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  theme_classic() + theme(legend.position="none") + ggtitle("\nlog(Rt)") + annotate("text", x=42, y=-0.225, parse=FALSE, label= "p = 0.1557")

ef <- effect("year: Groups", mm2)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
Figure42 <- ggplot(data = x, aes(x=year, y=fit, color=Groups)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  theme_classic() + theme(legend.position="none") + ggtitle("\nlog(Rc)") + annotate("text", x=42, y=0.32, parse=FALSE, label= "p = 0.0018")

ef <- effect("year: Groups", mm3)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
Figure43 <- ggplot(data = x, aes(x=year, y=fit, color=Groups)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  theme_classic() + theme(legend.position="none") + ggtitle("\nlog(Rs)") + annotate("text", x=42, y=0.04, parse=FALSE, label= "p = 0.0001")

ef <- effect("year: Groups", mm4)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
Figure44 <- ggplot(data = x, aes(x=year, y=fit, color=Groups)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  theme_classic() + theme(legend.position="none") + ggtitle(" \n") + annotate("text", x=42, y=-0.225, parse=FALSE, label= "p = 0.5930")

ef <- effect("year: Groups", mm5)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
Figure45 <- ggplot(data = x, aes(x=year, y=fit, color=Groups)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  theme_classic() + theme(legend.position="none") + ggtitle(" \n") + annotate("text", x=42, y=0.3, parse=FALSE, label= "p = 0.0269")

ef <- effect("year: Groups", mm6)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
Figure46 <- ggplot(data = x, aes(x=year, y=fit, color=Groups)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  theme_classic() + theme(legend.position=c(0.2,0.2)) + ggtitle(" \n") + annotate("text", x=42, y=0.02, parse=FALSE, label= "p = 0.0002")


#Export these 6 figures together
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,3))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(Figure41, vp = vplayout(1,1))   
print(Figure42, vp = vplayout(1,2))   
print(Figure43, vp = vplayout(1,3)) 
print(Figure44, vp = vplayout(2,1))   
print(Figure45, vp = vplayout(2,2))   
print(Figure46, vp = vplayout(2,3))  






























######## Null model #########
null1 <- lme(log(Resistance) ~ year + Groups + year*Groups, random=~1|site, na.action = na.omit, data = db1.5, control = list(opt = "optim"), method = "ML")
null2 <- lme(log(Recovery) ~ year + Groups + year*Groups, random=~1|site, na.action = na.omit, data = db1.5, control = list(opt = "optim"), method = "ML")
null3 <- lme(log(Resilience) ~ year + Groups + year*Groups, random=~1|site, na.action = na.omit, data = db1.5, control = list(opt = "optim"), method = "ML")

anova(null1)
anova(null2)
anova(null3)

summary(null1)
summary(null2)
summary(null3)

r.squaredGLMM(null1)
r.squaredGLMM(null2)
r.squaredGLMM(null3)


# Effect plot of log(Rx)~year in two Groups for the null model, like Figure 4
ef <- effect("year: Groups", null1)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
Figure47 <- ggplot(data = x, aes(x=year, y=fit, color=Groups)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  theme_classic() + theme(legend.position="none") + ggtitle(" \n") 

ef <- effect("year: Groups", null2)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
Figure48 <- ggplot(data = x, aes(x=year, y=fit, color=Groups)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  theme_classic() + theme(legend.position="none") + ggtitle(" \n")

ef <- effect("year: Groups", null3)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
Figure49 <- ggplot(data = x, aes(x=year, y=fit, color=Groups)) + geom_line() + geom_ribbon(aes(ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  theme_classic() + theme(legend.position=c(0.3,0.85)) + ggtitle(" \n")

#Export these 3 figures together
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,3))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(Figure47, vp = vplayout(1,1))   
print(Figure48, vp = vplayout(1,2))   
print(Figure49, vp = vplayout(1,3))  



# Standardized effects sizes plot for the null model, like Figure 1
Figure17 <- plot_model(null1, show.data = T, show.values = F, show.p = F, sort.est = F, title = "log(Rt)", axis.title = "Standardized effect size", 
                       vline.color = "grey", axis.labels = c("Year:Groups", "Groups", "Year")) + theme_bw()
Figure18 <- plot_model(null2, show.data = T, show.values = F, show.p = F, sort.est = F, title = "log(Rc)", axis.title = "Standardized effect size", 
                       vline.color = "grey", axis.labels = c("Year:Groups", "Groups", "Year")) + theme_bw()
Figure19 <- plot_model(null3, show.data = T, show.values = F, show.p = F, sort.est = F, title = "log(Rs)", axis.title = "Standardized effect size", 
                       vline.color = "grey", axis.labels = c("Year:Groups", "Groups", "Year")) + theme_bw()
#Export these 3 figures together
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,3))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(Figure17, vp = vplayout(1,1))   
print(Figure18, vp = vplayout(1,2))   
print(Figure19, vp = vplayout(1,3))  










################【【【【【2020622, Explore the correlation between changes in all-years SPEI and Rx in Model II for the different models】】】】】#############
# So basically we need to calculate the site-specific time slope of SPEI
# slope for all-year SPEI
# To fill the INFOSET_T with its column
slope_SPEId <- data.frame()
i <- 1
for (i in 1:2383) {
  value1 <- as.character(unique(INFOSET_T$site)[i]) #site name 
  value2 <- subset(INFOSET_T, site == value1)[,c(2,18)] #subset of this site
  value3 <- lm(value2$SPEId ~ value2$year) #linear model
  value4 <- coef(value3)[[2]] #the slope of SPEI
  value5 <- as.data.frame(rep(value4, 110))
  slope_SPEId <- rbind(slope_SPEId, value5)
  print(i)
}
View(slope_SPEId)
INFOSET_T <- cbind(INFOSET_T, slope_SPEId)
colnames(INFOSET_T)[46] <- "slope_SPEId"



# Build an unique dataset for slope_SPEId
# slope_SPEId1
slope_SPEId <- data.frame()
i <- 1
for (i in 1:2383) {
  value1 <- as.character(unique(INFOSET_T$site)[i]) #site name 
  value2 <- subset(INFOSET_T, site == value1)[,c(2,18)] #subset of this site with useful information 
  value3 <- lm(value2$SPEId ~ value2$year) #linear model
  value4 <- coef(value3)[[2]] #the slope of SPEI
  value5 <- c(value1, value4)
  slope_SPEId <- rbind(slope_SPEId, value5)
  print(i)
}
colnames(slope_SPEId) <- c("site", "slope_SPEId")
View(slope_SPEId)
slope_SPEId$slope_SPEId <- as.numeric(slope_SPEId$slope_SPEId)
hist(slope_SPEId$slope_SPEId)



#slope_SPEId2
slope_SPEId2 <- data.frame()
mod_SPEId <- lme(SPEId ~ year + Groups + year*Groups, random=~year|site, na.action = na.omit, data=INFOSET_T, control = list(opt = "optim"), method = "ML")
slope_SPEId2 <- as.data.frame(cbind(as.character(unique(INFOSET_T$site)), coef(mod_SPEId)[[2]])) 
# use coef instead of ranef, where coef=ranef+fixef (already gives the total coefficient);
# ranef is centred around zero; coef is centred around the fixef value
View(slope_SPEId2)
colnames(slope_SPEId2) <- c("site", "slope_SPEId2")
slope_SPEId2$slope_SPEId2 <- as.numeric(slope_SPEId2$slope_SPEId2)
hist(slope_SPEId2$slope_SPEId2)






# slope： log(Rs) ~ year for null model of SDY (db1.5)
null1 <- lme(log(Resilience) ~ year + Groups + year*Groups, random=~year|site, na.action = na.omit, data = db1.5, control = list(opt = "optim"), method = "ML")
slope_Rs_null <- data.frame
slope_Rs_null <- as.data.frame(cbind(rownames(coef(null1)), coef(null1)[[2]]))
colnames(slope_Rs_null) <- c("site", "slope_Rs_null")
slope_Rs_null$site <- as.factor(slope_Rs_null$site)
slope_Rs_null$slope_Rs_null <- as.numeric(slope_Rs_null$slope_Rs_null)
hist(slope_Rs_null$slope_Rs_null)








# slope： log(Rs) ~ year for Model II of SDY (db1.5), plus multi-slope
slope_Rs_model2 <- data.frame()
i <- 1
for(i in 1:nrow(coef(mm6))){
  value0 <- rownames(coef(mm6)[i,]) #the site
  value1 <- coef(mm6)[["year"]][i] #slope of year
  value2 <- coef(mm6)[["year:continentality_fixed"]][i] #slope of year:K
  value3 <- coef(mm6)[["year:GroupsAngiosperms"]][i] #slope of year:Groups
  value4 <- coef(mm6)[["year:log_Max_age"]][i] #slope of year:MA
  value5 <- coef(mm6)[["year:log_PRE_Annual_fixed"]][i] #slope of year:PRE
  value6 <- coef(mm6)[["year:TMP_Annual_fixed"]][i] #slope of year:TMP
  value7 <- subset(db1.5, site == value0)[1,]$continentality
  value8 <- as.character(subset(db1.5, site == value0)[1,]$Groups)
  value9 <- subset(db1.5, site == value0)[1,]$log_Max_age
  value10 <- subset(db1.5, site == value0)[1,]$log_PRE_Annual_fixed
  value11 <- subset(db1.5, site == value0)[1,]$TMP_Annual_fixed
  if(value8 == "Angiosperms"){
    value12 <- value1 + value2*value7 + value3 + value4*value9 + value5*value10 + value6*value11
  }else{
    value12 <- value1 + value2*value7 + value4*value9 + value5*value10 + value6*value11
  }
  slope_Rs_model2 <- rbind(slope_Rs_model2, c(value0, value1, value12)) 
  print(i)
}
rownames(slope_Rs_model2) <- 1:nrow(slope_Rs_model2)
colnames(slope_Rs_model2) <- c("site", "slope_year_model2", "slope_Rs_model2")
View(slope_Rs_model2)
slope_Rs_model2$slope_year_model2 <- as.numeric(slope_Rs_model2$slope_year_model2)
slope_Rs_model2$slope_Rs_model2 <- as.numeric(slope_Rs_model2$slope_Rs_model2)
hist(slope_Rs_model2$slope_Rs_model2)


slope_compare <- data.frame()
slope_compare <- left_join(slope_SPEId, slope_SPEId2, by = "site")
slope_compare <- left_join(slope_compare, slope_Rs_null, by = "site")
slope_compare <- left_join(slope_compare, slope_Rs_model2, by = "site")
slope_compare <- slope_compare[,-5]
write.table(slope_compare, file = "slope_compare.csv", row.names = F, na = "NA", sep = ",")

coord <- left_join(coord, slope_compare, by = "site")







###########【【【【【Figure test1, Slope comparison plot for SPEId, log(Rs)_null and log(Rs)_model2, 20200623 】】】】】############
Figuret1 <- NULL
mapWorld <- borders("world", colour="grey", fill="white") # create a layer of borders
Figuret1 <- ggplot() + mapWorld
#Now Layer the cities on top
Figuret1 <- Figuret1 + geom_point(aes(x=coord$x.real, y=coord$y.real, color=coord$slope_SPEId), size=1) + theme_bw() + 
  labs(title = "", x = "Longitude (°)", y = "Latitude (°)") + 
  scale_colour_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0, name = "slope") + theme(legend.position = c(0.1, 0.3)) + 
  ggtitle("Site-specific slopes of SPEId ~ year, all-year data")
Figuret1



Figuret2 <- NULL
mapWorld <- borders("world", colour="grey", fill="white") # create a layer of borders
Figuret2 <- ggplot() + mapWorld
#Now Layer the cities on top
Figuret2 <- Figuret2 + geom_point(aes(x=coord$x.real, y=coord$y.real, color=coord$slope_SPEId2), size=1) + theme_bw() + 
  labs(title = "", x = "Longitude (°)", y = "Latitude (°)") + 
  scale_colour_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0, name = "slope") + theme(legend.position = c(0.1, 0.3)) + 
  ggtitle("Site-specific slopes of SPEId2 ~ year, all-year data")
Figuret2


Figuret3 <- NULL
mapWorld <- borders("world", colour="grey", fill="white") # create a layer of borders
Figuret3 <- ggplot() + mapWorld
#Now Layer the cities on top
Figuret3 <- Figuret3 + geom_point(aes(x=coord$x.real, y=coord$y.real, color=coord$slope_Rs_null), size=1) + theme_bw() + 
  labs(title = "", x = "Longitude (°)", y = "Latitude (°)") + 
  scale_colour_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0, name = "slope") + theme(legend.position = c(0.1, 0.3)) + 
  ggtitle("Slopes of log(Resilience) ~ year, null model from SDY data")
Figuret3



Figuret4 <- NULL
mapWorld <- borders("world", colour="grey", fill="white") # create a layer of borders
Figuret4 <- ggplot() + mapWorld
#Now Layer the cities on top
Figuret4 <- Figuret4 + geom_point(aes(x=coord$x.real, y=coord$y.real, color=coord$slope_Rs_model2), size=1) + theme_bw() + 
  labs(title = "", x = "Longitude (°)", y = "Latitude (°)") + 
  scale_colour_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 0, name = "slope") + theme(legend.position = c(0.1, 0.3)) + 
  ggtitle("Slopes of log(Resilience) ~ year, model II from SDY data")
Figuret4






###########【【【【【20200624, Spatial regression model for slope_Rx_model2 ~ slope_SPEI* Groups 】】】】】############

summary(lm(slope_Rs_null ~ slope_SPEId* Groups, data = coord, na.action = na.omit))
summary(lm(slope_Rs_model2 ~ slope_SPEId* Groups, data = coord, na.action = na.omit))

summary(lm(slope_Rs_null ~ slope_SPEId2* Groups, data = coord, na.action = na.omit))
summary(lm(slope_Rs_model2 ~ slope_SPEId2* Groups, data = coord, na.action = na.omit))









########################【【【【【Pre/post 1950 models】】】】】##########################
#### Updated on 20200603, we have RWI4 log-transformed in the model
#### Updated on 20200603, we have the interaction of Groups*log_RWI4 in the models mm4 to mm6 now
#### Updated on 20200614, we have correct these wrong data in SPEI3m_Min and SPEI3m_Min_Season, so let's run the main model again! 
load(file = "INFOSET_T.r") #the largest dataset, including the information of the years from 1901 to 2010

INFOSET_T$RWI4prev <- with(INFOSET_T, ave(RWI4prev, site, FUN=function(x) x - mean(x,na.rm=TRUE)))
INFOSET_T$RWI4post <- with(INFOSET_T, ave(RWI4post, site, FUN=function(x) x - mean(x,na.rm=TRUE)))
INFOSET_T$year <- with(INFOSET_T, year - mean(year, na.rm=TRUE))
INFOSET_T$TMP_Annual <- with(INFOSET_T, TMP_Annual - mean(TMP_Annual, na.rm=TRUE))
INFOSET_T$TMP_Annual_fixed <- with(INFOSET_T, TMP_Annual_fixed - mean(TMP_Annual_fixed, na.rm=TRUE))
INFOSET_T$PRE_Annual <- with(INFOSET_T, PRE_Annual - mean(PRE_Annual, na.rm=TRUE))
INFOSET_T$PRE_Annual_fixed <- with(INFOSET_T, PRE_Annual_fixed - mean(PRE_Annual_fixed, na.rm=TRUE))
INFOSET_T$Max_age <- with(INFOSET_T, Max_age - mean(Max_age, na.rm=TRUE))
INFOSET_T$continentality <- with(INFOSET_T, continentality - mean(continentality, na.rm=TRUE))
INFOSET_T$continentality_fixed <- with(INFOSET_T, continentality_fixed - mean(continentality_fixed, na.rm=TRUE)) 
INFOSET_T$log_Max_age <- with(INFOSET_T, log_Max_age - mean(log_Max_age, na.rm=TRUE))
INFOSET_T$log_PRE_Annual_fixed <- with(INFOSET_T, log_PRE_Annual_fixed - mean(log_PRE_Annual_fixed, na.rm=TRUE)) 

INFOSET <- subset(INFOSET_T, is.na(TRI) == FALSE)
db1.5_T <- subset(INFOSET_T, SPEId + 1.5 < 0)
db1.5 <- subset(INFOSET, SPEId + 1.5 < 0)

# subsets of time periods
db1.5_0020 <- subset(db1.5, (year + 35.5) <= 0)
db1.5_2040 <- subset(db1.5, (year + 15.5)  <= 0 & (year + 35.5) > 0)
db1.5_4060 <- subset(db1.5, (year - 5.5)  < 0 & (year + 15.5) > 0)
db1.5_6080 <- subset(db1.5, (year - 25.5)  < 0 & (year - 5.5) > 0)
db1.5_8000 <- subset(db1.5, (year - 45.5)  < 0 & (year - 25.5) > 0)
db1.5_over2000 <- subset(db1.5, (year - 45.5) >= 0)

db1.5_0020$pre_post   <- rep('00_20', length(db1.5_0020$year))
db1.5_2040$pre_post <-   rep('20_40',  length(db1.5_2040$year))
db1.5_4060$pre_post   <- rep('40_60', length(db1.5_4060$year))
db1.5_6080$pre_post   <- rep('60_80', length(db1.5_6080$year))
db1.5_8000$pre_post   <- rep('80_00', length(db1.5_8000$year))
db1.5_over2000$pre_post   <- rep('over2000', length(db1.5_over2000$year))

db1.5_pre_post <- rbind(db1.5_0020, db1.5_2040, db1.5_4060, db1.5_6080, db1.5_8000, db1.5_over2000)
db1.5_pre_post$pre_post <- as.factor(db1.5_pre_post$pre_post)
db1.5_pre_post$pre_post <- factor(db1.5_pre_post$pre_post, levels = c('00_20', '20_40', '40_60', 
                                                                      '60_80', '80_00', 'over2000'))

#part 1: the model without CUM4SPEI, SPEId and Season, with sTMP and sPRE and Species_code
mm1spp <- lme(log(Resistance) ~ year + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + Groups + slope_TMPxYear + slope_PRExYear + 
             year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + year* log_PRE_Annual_fixed + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.5_pre1950, control = list(opt = "optim"), method = "ML")
mm2spp <- lme(log(Recovery) ~ year + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + Groups + slope_TMPxYear + slope_PRExYear + 
             year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + year* log_PRE_Annual_fixed + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.5_pre1950, control = list(opt = "optim"), method = "ML")
mm3spp <- lme(log(Resilience) ~ year + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + Groups + slope_TMPxYear + slope_PRExYear + 
             year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + year* log_PRE_Annual_fixed + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.5_pre1950, control = list(opt = "optim"), method = "ML")

#part 2: the model with CUM4SPEI, SPEId, Season and interactions, with sTMP and sPRE and Species_code
mm4spp <- lme(log(Resistance) ~ year + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
             Groups + slope_TMPxYear + slope_PRExYear + year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + 
             year* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.5_pre1950, control = list(opt = "optim"), method = "ML")
mm5spp <- lme(log(Recovery) ~   year + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
             Groups + slope_TMPxYear + slope_PRExYear + year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + 
             year* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.5_pre1950, control = list(opt = "optim"), method = "ML")
mm6spp <- lme(log(Resilience) ~ year + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
             Groups + slope_TMPxYear + slope_PRExYear + year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + 
             year* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.5_pre1950, control = list(opt = "optim"), method = "ML")

anova(mm1spp)
anova(mm2spp)
anova(mm3spp)
anova(mm4spp)
anova(mm5spp)
anova(mm6spp)

summary(mm1spp)
summary(mm2spp)
summary(mm3spp)
summary(mm4spp)
summary(mm5spp)
summary(mm6spp)

r.squaredGLMM(mm1spp)
r.squaredGLMM(mm2spp)
r.squaredGLMM(mm3spp)
r.squaredGLMM(mm4spp)
r.squaredGLMM(mm5spp)
r.squaredGLMM(mm6spp)



#part 3: the model with CUM4SPEI, SPEId, Season and interactions, Species_code, 
# no Year, no sTMP no sPRE 
# but with factor pre/post1950
mm1_prepost <- lme(log(Resistance) ~ log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
                Groups + pre_post +  
                Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4 + Groups * pre_post, 
              random = ~1|Species_code/site, na.action = na.omit, data = db1.5_pre_post, control = list(opt = "optim"), method = "ML")
mm2_prepost <- lme(log(Recovery) ~   log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
                Groups + pre_post + 
                Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4 + Groups * pre_post, 
              random = ~1|Species_code/site, na.action = na.omit, data = db1.5_pre_post, control = list(opt = "optim"), method = "ML")
mm3_prepost <- lme(log(Resilience) ~ log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
                Groups + pre_post +  
                Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4 + Groups * pre_post, 
              random = ~1|Species_code/site, na.action = na.omit, data = db1.5_pre_post, control = list(opt = "optim"), method = "ML")


intervals(mm1_prepost)
intervals(mm2_prepost)
intervals(mm3_prepost)


anova(mm1_prepost)
anova(mm2_prepost)
anova(mm3_prepost)

summary(mm1_prepost)
summary(mm2_prepost)
summary(mm3_prepost)



ef <- effect("Groups : pre_post", mm1_prepost)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  ggtitle('Resistance') +
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)


ef <- effect("Groups : pre_post", mm2_prepost)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  ggtitle('Recovery') +
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)

ef <- effect("Groups : pre_post", mm3_prepost)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
    ggtitle('Resilience') +
  geom_point(position = position_dodge(width = 0.3)) + 
    geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)
  








########################【【【【【Testing Li et al calculations ��】】】】##########################



seriesoriginal <- read.table(file = "series_without0_fixed.txt", header = T)

series_Li <- data.frame()
i <- 1
for (i in 1:length(unique(INFOSET_T$site))) {
  value1 <- as.character(unique(INFOSET_T$site))[i] #site name
  value2 <- which(colnames(seriesoriginal) == value1) #which column is this site
  value3 <- seriesoriginal[5:114,value2] #the subset of this site, RWI
  value4 <- as.data.frame(cbind(rep(value1, 110), 1901:2010, value3)) #c(site, year, RWI)
  Lipost <- c(value3[-(1:2)], NA, NA) # takes Ye+2 value as in Li et al for each Ye
  value5 <- cbind(value4, Lipost)
  colnames(value5) <- c("site_Li", "year_Li", "RWI_Li", "Gpost_Li")
  series_Li <- rbind(series_Li, value5)
  print(i)
}


# employs Li et al (2020) definition of 10% quantile (quantile 1.28 gives 10.03% exact) for 'extreme' droughts
db1.28_nondrgt <- subset(INFOSET_T, SPEId + 1.28 > 0) # extracts non-drought years according to Li et al


# mean growth during non-drought years by site for use by Li et al
MeanYn <- with(db1.28_nondrgt, aggregate(TRI, by=list(site), FUN=meanNA))[,2]

# #some changes in INFOSET_T
INFOSET_T <- cbind(INFOSET_T, series_Li[,1:4])
save(series_Li, file = "series_Li")
INFOSET_T$year_Li <- as.numeric(INFOSET_T$year_Li)
INFOSET_T$RWI_Li <- as.numeric(INFOSET_T$RWI_Li)
INFOSET_T$MeanYn_Li <- rep(MeanYn, 1, each=110) # mean Ym by site according to Li et al
save(INFOSET_T, file = "INFOSET_T.r")


# subset INFOSET_T now after adding the Ym from Li et al
# but apparently they did not use this approach
db1.28_Li <- subset(INFOSET_T, SPEId + 1.28 <= 0)


# db1.28_Li2 <- subset(INFOSET_T, SPEId + 1.28 <= 0)

# subset INFOSET_T after changing old df to db1.28_Li2 and implemented their 10% percentile approach
# db1.28_Li <- as.data.frame(matrix(nrow = 0, ncol = ncol(db1.28_Li2)))
# colnames(db1.28_Li) <- colnames(db1.28_Li2)
# for (i in 1:length(unique(INFOSET_T$site))) {
#   value1 <- as.character(unique(INFOSET_T$site))[i] #site name
#   pippo <- INFOSET_T[which(INFOSET_T$site==value1),]
#   qle <- with(pippo, quantile(SPEId, 0.1, na.rm=TRUE))
#   pippo <- pippo[which(pippo$SPEId <= qle),]
#   db1.28_Li <- rbind(db1.28_Li, pippo)
#   print(i)
# }




# formulas as in Li et al (2020)
db1.28_Li$Li_Rt <- with(db1.28_Li, MeanYn_Li/abs(RWI_Li-MeanYn_Li))
db1.28_Li$Li_Rs <- with(db1.28_Li, abs((RWI_Li-MeanYn_Li)/(Gpost_Li-MeanYn_Li)))



# db1.28_Li$site <- as.character(db1.28_Li$site)
# 
# # ok, code below shows sites all match, both before and after cleaning
# yes=0
# for(i in 1:length(db1.28_Li$Li_Rt)) {
# 
#   if(db1.28_Li$site[i]==db1.28_Li$site_Li[i]) {yes=yes+1}
# print(yes)
#   }
# 
# db1.28_Li$site <- as.factor(db1.28_Li$site)



# some manual cleaning to do to do with NAs

# # eliminate +Inf in both Rt and Rs (uses library(data.table))
invisible(lapply(names(db1.28_Li),function(.name) set(db1.28_Li, which(is.infinite(db1.28_Li[[.name]])), 
                                                      j = .name,value =NA)))
# eliminate NAs in both Rt and Rs (deletes about 6,000 NAs)
# db1.28_Li <- db1.28_Li[complete.cases(db1.28_Li[ , 'Li_Rt']),]
db1.28_Li <- db1.28_Li[complete.cases(db1.28_Li[ , 'Li_Rs']),]

# # eliminate cases of Rs=0
# db1.28_Li$Li_Rt <- db1.28_Li$Li_Rt + 0.001
db1.28_Li$Li_Rs <- db1.28_Li$Li_Rs + 0.001



# subsets of time periods
db1.28_0229 <- subset(db1.28_Li, year_Li <= 1929)
db1.28_3049 <- subset(db1.28_Li, year_Li  <= 1949  & year_Li > 1929)
db1.28_5069 <- subset(db1.28_Li, year_Li  <= 1969  & year_Li > 1949)
db1.28_7089 <- subset(db1.28_Li, year_Li <= 1989 &  year_Li  > 1969)
db1.28_9008 <- subset(db1.28_Li,             year_Li  > 1989)

db1.28_0229$pre_post   <- rep('02_29', length(db1.28_0229$year)) #7229
db1.28_3049$pre_post <-   rep('30_49', length(db1.28_3049$year)) #5208
db1.28_5069$pre_post   <- rep('50_69', length(db1.28_5069$year)) #5354
db1.28_7089$pre_post   <- rep('70_89', length(db1.28_7089$year)) #3130
db1.28_9008$pre_post   <- rep('90_08', length(db1.28_9008$year)) #722


db1.28_pre_post <- rbind(db1.28_0229, db1.28_3049, db1.28_5069, db1.28_7089, db1.28_9008)
db1.28_pre_post$pre_post <- as.factor(db1.28_pre_post$pre_post)
db1.28_pre_post$pre_post <- factor(db1.28_pre_post$pre_post, levels = c('02_29', '30_49', '50_69', 
                                                                      '70_89', '90_08'))




####____________________ MODELS WITH ANGIOS and GYMNOS

#part 0: the model with only Year * Group
mm0 <- lme(log(Li_Rt) ~ year + Groups + year* Groups, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.28_Li, control = list(opt = "optim"), method = "REML")
mm000 <- lme(log(Li_Rs) ~ year + Groups + year* Groups, 
             random = ~1|Species_code/site, na.action = na.omit, data = db1.28_Li, control = list(opt = "optim"), method = "REML")

#part 1: the model without CUM4SPEI, SPEId and Season
mm1 <- lme(log(Li_Rt) ~ year + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + Groups +  
             year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + year* log_PRE_Annual_fixed + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.28_Li, control = list(opt = "optim"), method = "REML")
mm3 <- lme(log(Li_Rs) ~ year + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + Groups +  
             year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + year* log_PRE_Annual_fixed + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.28_Li, control = list(opt = "optim"), method = "REML")

#part 2: the model with CUM4SPEI, SPEId, Season and interactions
mm4 <- lme(log(Li_Rt) ~ year + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
             Groups + slope_TMPxYear + slope_PRExYear + year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + 
             year* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.28_Li, control = list(opt = "optim"), method = "REML")
mm6 <- lme(log(Li_Rs) ~ year + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed + 
             Groups + slope_TMPxYear + slope_PRExYear + year* continentality_fixed + year* Groups + year* log_Max_age + year* TMP_Annual_fixed + 
             year* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4, 
           random = ~1|Species_code/site, na.action = na.omit, data = db1.28_Li, control = list(opt = "optim"), method = "REML")

anova(mm0)
anova(mm000)
anova(mm1)
anova(mm3)
anova(mm4)
anova(mm6)

summary(mm0)
summary(mm000)
summary(mm1)
summary(mm3)
summary(mm4)
summary(mm6)

r.squaredGLMM(mm0)
r.squaredGLMM(mm000)
r.squaredGLMM(mm1)
r.squaredGLMM(mm3)
r.squaredGLMM(mm4)
r.squaredGLMM(mm6)






#part 3: the model with CUM4SPEI, SPEId, Season and interactions, Species_code, 
# no Year, no sTMP no sPRE 
# but with factor pre/post1950
# cannot have variable slopes at both levels
mm0_prepost_Li <- lmer(log(Li_Rt) ~ Groups + pre_post + Groups * pre_post +
                         (1 | Species_code/site), na.action = na.omit, data = db1.28_pre_post)

# cannot have variable slopes at both levels
mm1_prepost_Li <- lmer(log(Li_Rt) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                        Groups + slope_TMPxYear + slope_PRExYear + pre_post* Groups + pre_post* log_Max_age + pre_post* TMP_Annual_fixed + 
                        Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4 +
                      (1|Species_code/site), na.action = na.omit, data = db1.28_pre_post, REML=TRUE)


# pre_post : K,pre_post : log_PRE, pre_post: Groups, pre_post: SPEId highly collinear
mm3_prepost_Li <- lmer(log(Li_Rs) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                        Groups + slope_TMPxYear + slope_PRExYear + pre_post* log_Max_age + pre_post* TMP_Annual_fixed + 
                        Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4 +
                      (pre_post|Species_code), na.action = na.omit, data = db1.28_pre_post, REML = TRUE)

# mm0_prepost_Li <- lme(log(Li_Rt) ~ Groups + pre_post +Groups * pre_post, 
#                    random = ~1|Species_code/site, na.action = na.omit, data = db1.28_pre_post, control = list(opt = "optim"), method = "ML")
# 
# 
# 
# mm1_prepost_Li <- lme(log(Li_Rt) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
#                      Groups + slope_TMPxYear + slope_PRExYear + pre_post* continentality_fixed + pre_post* Groups + pre_post* log_Max_age + pre_post* TMP_Annual_fixed + 
#                pre_post* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4,
#              random = ~1|Species_code/site, na.action = na.omit, data = db1.28_pre_post, control = list(opt = "optim"), method = "REML")
# 

# mm00_prepost_Li <- lme(log(Li_Rs) ~ Groups + pre_post +Groups * pre_post, 
#                    random = ~1|Species_code, na.action = na.omit, data = db1.28_pre_post, control = list(opt = "optim"), method = "ML")
# mm3_prepost_Li <- lme(log(Li_Rs) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
#                      Groups + slope_TMPxYear + slope_PRExYear + pre_post* continentality_fixed + pre_post* Groups + pre_post* log_Max_age + pre_post* TMP_Annual_fixed + 
#                      pre_post* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4,
#                    random = ~1|Species_code, na.action = na.omit, data = db1.28_pre_post, control = list(opt = "optim"), method = "REML")


intervals(mm1_prepost_Li)
intervals(mm0_prepost_Li)
intervals(mm3_prepost_Li)
intervals(mm00_prepost_Li)


anova(mm1_prepost_Li)
anova(mm3_prepost_Li)

summary(mm1_prepost_Li)
summary(mm0_prepost_Li)
summary(mm3_prepost_Li)
summary(mm00_prepost_Li)

r.squaredGLMM(mm0_prepost_Li)
r.squaredGLMM(mm00_prepost_Li)
r.squaredGLMM(mm1_prepost_Li)
r.squaredGLMM(mm3_prepost_Li)

# checks model quality for mixed models using library(performance)
check_model(mm0_prepost_Li)
check_model(mm00_prepost_Li)
chks1 <- check_model(mm1_prepost_Li)
chks3 <- check_model(mm3_prepost_Li)

cbind(chks1$VIF$x, format(chks1$VIF$y, digits=2))
cbind(chks3$VIF$x, format(chks3$VIF$y, digits=2))

# model_performance(mm0_prepost_Li)
model_performance(mm00_prepost_Li)
model_performance(mm1_prepost_Li)
model_performance(mm3_prepost_Li)

ef <- effect("pre_post:Groups", mm1_prepost_Li)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  labs(title = "Isbell indices", x = "Period in 20th century", y = "Resistance") + 
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)

ef <- effect("Groups : pre_post", mm0_prepost_Li)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  labs(title = "Isbell indices", x = "Period in 20th century", y = "Resistance") + 
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)

ef <- effect("pre_post: Groups", mm3_prepost_Li)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  labs(title = "Isbell indices", x = "Period in 20th century", y = "Resilience") + 
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)

ef <- effect("Groups : pre_post", mm00_prepost_Li)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  labs(title = "Isbell indices", x = "Period in 20th century", y = "Resilience") + 
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)



#part 4: Tong's models

# mm0_prepost_Tong <- lmer(log(Resistance) ~ Groups + pre_post +Groups * pre_post +
#                         (1|Species_code/site), na.action = na.omit, data = db1.28_pre_post, REML = TRUE)
# 
# mm1_prepost_Tong <- lmer(log(Resistance) ~ pre_post + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
#                           Groups + slope_TMPxYear + slope_PRExYear + pre_post* continentality_fixed + pre_post* Groups + pre_post* log_Max_age + pre_post* TMP_Annual_fixed + 
#                           pre_post* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4 +
#                         (1|Species_code/site), na.action = na.omit, data = db1.28_pre_post, REML = TRUE)


mm00_prepost_Tong <- lmer(log(Resilience) ~ Groups + pre_post +Groups * pre_post +
                         (1|Species_code), na.action = na.omit, data = db1.28_pre_post, REML = TRUE)




# pre_post : K,pre_post : log_PRE highly collinear
mm3_prepost_Tong <- lmer(log(Resilience) ~ pre_post  + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                          Groups + slope_TMPxYear + slope_PRExYear + pre_post* Groups + pre_post* log_Max_age + pre_post* TMP_Annual_fixed + 
                          Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4 +
                        (1|Species_code/site), na.action = na.omit, data = db1.28_pre_post, REML = TRUE)



# mm0_prepost_Tong <- lme(log(Resistance) ~ Groups + pre_post +Groups * pre_post, 
#                    random = ~1|Species_code/site, na.action = na.omit, data = db1.28_pre_post, control = list(opt = "optim"), method = "ML")
# 
# mm1_prepost_Tong <- lme(log(Resistance) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
#                      Groups + slope_TMPxYear + slope_PRExYear + pre_post* continentality_fixed + pre_post* Groups + pre_post* log_Max_age + pre_post* TMP_Annual_fixed + 
#                      pre_post* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4,
#                    random = ~1|Species_code/site, na.action = na.omit, data = db1.28_pre_post, control = list(opt = "optim"), method = "REML")
# 
# 
# mm00_prepost_Tong <- lme(log(Resilience) ~ Groups + pre_post +Groups * pre_post, 
#                     random = ~1|Species_code, na.action = na.omit, data = db1.28_pre_post, control = list(opt = "optim"), method = "ML")
# mm3_prepost_Tong <- lme(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
#                      Groups + slope_TMPxYear + slope_PRExYear + pre_post* continentality_fixed + pre_post* Groups + pre_post* log_Max_age + pre_post* TMP_Annual_fixed + 
#                      pre_post* log_PRE_Annual_fixed + Groups* SPEId + Groups* SPEIprev + Groups* SPEIpost + Groups* log_RWI4,
#                    random = ~1|Species_code/site, na.action = na.omit, data = db1.28_pre_post, control = list(opt = "optim"), method = "REML")
# 

intervals(mm1_prepost_Tong)
intervals(mm0_prepost_Tong)
intervals(mm3_prepost_Tong)
intervals(mm00_prepost_Tong)

anova(mm1_prepost_Tong)
anova(mm3_prepost_Tong)

summary(mm1_prepost_Tong)
summary(mm0_prepost_Tong)
summary(mm3_prepost_Tong)
summary(mm00_prepost_Tong)

r.squaredGLMM(mm0_prepost_Tong)
r.squaredGLMM(mm00_prepost_Tong)
r.squaredGLMM(mm1_prepost_Tong)
r.squaredGLMM(mm3_prepost_Tong)

check_model(mm0_prepost_Tong)
check_model(mm00_prepost_Tong)
check_model(mm1_prepost_Tong)
chks3 <-check_model(mm3_prepost_Tong)
cbind(chks3$VIF$x, format(chks3$VIF$y, digits=2))

model_performance(mm0_prepost_Tong)
model_performance(mm00_prepost_Tong)
model_performance(mm1_prepost_Tong)
model_performance(mm3_prepost_Tong)

check_singularity(mm1_prepost_Tong)
check_singularity(mm3_prepost_Tong)

ef <- effect("pre_post:Groups", mm1_prepost_Tong)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  labs(title = "Lloret indices", x = "Period in 20th century", y = "Resistance") + 
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)

ef <- effect("Groups : pre_post", mm0_prepost_Tong)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  labs(title = "Lloret indices", x = "Period in 20th century", y = "Resistance") + 
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)

ef <- effect("pre_post:Groups", mm3_prepost_Tong)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  labs(title = "Lloret indices", x = "Period in 20th century", y = "Resilience") + 
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)

ef <- effect("Groups : pre_post", mm00_prepost_Tong)
x <- as.data.frame(ef)
x$Groups <- as.factor(as.character(x$Groups))
ggplot(data = x, aes(x=pre_post, y=fit, color=Groups)) + 
  labs(title = "Lloret indices", x = "Period in 20th century", y = "Resilience") + 
  geom_point(position = position_dodge(width = 0.3)) + 
  geom_errorbar(position = position_dodge(width = 0.3), aes(ymin=fit-se, ymax=fit+se), width=0.4) + theme_bw(base_size=12)

########################【【【【【Testing Li et al / only CONIFERS ��】】】】##########################

db1.28_Gymn <- subset(db1.28_pre_post, Groups == 'Gymnosperms')


# scale all numeric predictors before analysis
db1.28_Gymn$log_RWI4             <- as.vector(scale(db1.28_Gymn$log_RWI4))
db1.28_Gymn$log_Max_age          <- as.vector(scale(db1.28_Gymn$log_Max_age))
db1.28_Gymn$TMP_Annual_fixed     <- as.vector(scale(db1.28_Gymn$TMP_Annual_fixed))
db1.28_Gymn$log_PRE_Annual_fixed <- as.vector(scale(db1.28_Gymn$log_PRE_Annual_fixed))
db1.28_Gymn$continentality_fixed <- as.vector(scale(db1.28_Gymn$continentality_fixed))
db1.28_Gymn$slope_TMPxYear       <- as.vector(scale(db1.28_Gymn$slope_TMPxYear))
db1.28_Gymn$slope_PRExYear       <- as.vector(scale(db1.28_Gymn$slope_PRExYear))

# some obvious problems with outliers in these data; let's clean them
# > with(db1.28_Gymn, summary(Li_Rt))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
# 0.6965    2.9417    5.4670   20.2208   12.1771 2271.0000        27 
# > with(db1.28_Gymn, summary(Li_Rs))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0010    0.5341    1.1878    5.6955    2.7764 1959.0010 

with(db1.28_Gymn, hist(log(Li_Rt+1)))

#part 3: Li indices with factor pre/post1950
# cannot have variable slope at both levels; cannot even have intercept for site or slope for species
# can only have intercept for species; lmer converges for site intercept, but stdev=e-23!!
mm0_prepost_Li <- lmer(log(Li_Rs) ~ pre_post +
                         (1 | Species_code), 
                       control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
                       na.action = na.omit, data = db1.28_Gymn)

# cannot have variable slope at either levels; cannot even have intercept for site or slope for species
# can only have intercept for species
mm2_prepost_Li <- lmer(log(Li_Rs) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season +
                         log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                         # slope_TMPxYear + slope_PRExYear +
                         pre_post * log_Max_age + pre_post * TMP_Annual_fixed + pre_post * continentality_fixed +
                         (1|Species_code), na.action = na.omit, data = db1.28_Gymn, REML=TRUE)

# even with optimX no luck in having a variable slope or a random site
pip1 <- lmer(log(Li_Rs) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season +
               log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
               # slope_TMPxYear + slope_PRExYear +  
               pre_post * log_Max_age + pre_post * TMP_Annual_fixed + pre_post * continentality_fixed + 
               (1|Species_code/site),
             control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')),
             na.action = na.omit, data = db1.28_Gymn, REML=TRUE)


# type='I' is sequential SS, testing Time first
anova(mm0_prepost_Li, type='I') # F=1.6603 P=0.1562 for pre_post
anova(mm2_prepost_Li, type='I') # F=2.0431 P=0.08552 for pre_post

summary(mm0_prepost_Li)
summary(mm2_prepost_Li)

# r.squaredGLMM(mm0_prepost_Li)
r.squaredGLMM(mm0_prepost_Li) # R2m = 0.0005200253  , R2c = 0.01086297
r.squaredGLMM(mm2_prepost_Li) # R2m = 0.02200544,   R2c = 0.02543616

# checks model quality for mixed models using library(performance) (needs lme4)
# check_model(mm0_prepost_Li)
check_model(mm0_prepost_Li)
chks02 <- check_model(mm2_prepost_Li)
cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))

# model_performance(mm0_prepost_Li)
model_performance(mm0_prepost_Li) # AIC=64005.30 , R2c=7.00e-03, R2m=0, ICC=7.00e-03, RMSE=1.50
model_performance(mm2_prepost_Li) # AIC=62427.03 , R2c=0.02, R2m=0.02, ICC=4.00e-03, RMSE=1.49

check_singularity(mm2_prepost_Li) # FALSE
check_singularity(mm0_prepost_Li) # FALSE



#____________________________________________part 4: Tong's models


# variable slope much better than variable intercept / cannot have both variable
mm0_prepost_Tong <- lmer(log(Resilience) ~ pre_post +
                         (1|Species_code), 
                         control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
                         na.action = na.omit, data = db1.28_Gymn, REML = FALSE)

# pre_post : continentality_fixed and  pre_post : log_PRE_Annual_fixed highly multicollinear
mm2_prepost_Tong <- lmer(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
                            log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                           # slope_TMPxYear + slope_PRExYear + 
                           pre_post * log_Max_age + pre_post * continentality_fixed + pre_post * TMP_Annual_fixed +
                           (1|Species_code/site), na.action = na.omit, data = db1.28_Gymn, REML = FALSE)

# definitely the best !!!
pip4 <- lmer(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
               log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
               # slope_TMPxYear + slope_PRExYear + 
               pre_post * log_Max_age + pre_post * continentality_fixed + pre_post * TMP_Annual_fixed + 
               (pre_post|Species_code) + (1|Species_code:site), 
             control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
             na.action = na.omit, data = db1.28_Gymn, REML = FALSE)




anova(mm0_prepost_Tong, type='I') # F=24.2  P <.0001 for pre_post
anova(mm2_prepost_Tong, type='I') # F= 28.8 P <.0001 for pre_post

# interestingly, if I include randomly varying slopes (by Species), the pre_post effect disappears (P=0.45)
anova(pip4, type='I')  # I have loaded lmerTest; .
anova(pip5)              # F=1.716  P=0.1433 for pre_post; and P=0.1571 if I include interactions of pre_post with
# slope_TMPxYear and slope_PRExYear. VERY INTERESTING!! Hence, suggestion is that time effect disappears when 
# I include variable random slopes by species??

summary(mm0_prepost_Tong)
summary(mm2_prepost_Tong)
summary(pip4)

# r.squaredGLMM(mm0_prepost_Tong)
r.squaredGLMM(mm0_prepost_Tong) # R2m = 0.005643217   ;  R2c = 0.01453085
r.squaredGLMM(mm2_prepost_Tong) # R2m = 0.415607   ;   R2c = 0.4559308
r.squaredGLMM(pip4)     # R2m = 0.4270828     ;   R2c = 0.5214921 !!!
r.squaredGLMM(pip5)     # R2m = 0.4110719  ;   R2c = 0.5103781 !!!

# check_model(mm0_prepost_Tong)
check_model(mm0_prepost_Tong)
chks02 <- check_model(mm2_prepost_Tong)
cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))
check_model(pip4)

# model_performance(mm0_prepost_Tong)
model_performance(mm0_prepost_Tong) # AIC=-481.13  , R2c=0.01, R2m=6.00e-3, ICC=0.01, RMSE=0.24
model_performance(mm2_prepost_Tong) # AIC=-9563.00  , R2c=0.46, R2m=0.42,    ICC=0.07, RMSE=0.17
model_performance(pip4)             # AIC=-10097.14, R2c=0.51, R2m=0.42,    ICC=0.15, RMSE=0.17

check_singularity(mm2_prepost_Tong) # FALSE
check_singularity(mm0_prepost_Tong) # FALSE
check_singularity(pip4) # TRUE

pairs(emmeans(mm0_prepost_Li, 'pre_post', type = 'response'))
pairs(emmeans(mm0_prepost_Tong, 'pre_post', type = 'response'))

pairs(emmeans(mm2_prepost_Li, 'pre_post', type = 'response'))
pairs(emmeans(mm2_prepost_Tong, 'pre_post', type = 'response'))



# USE OF EMMEANS WITH LMER4
# There is an optional lmer.df argument that defaults to get_EMM_option("lmer.df") (which in turn defaults to
# "kenward-roger"). The possible values are "kenward-roger", "satterthwaite", and "asymptotic" (these are partially 
# matched and case-insensitive). With "kenward-roger", d.f. are obtained using code from the pbkrtest package, if 
# installed. With "satterthwaite", d.f. are obtained using code from the lmerTest package, if installed. With 
# "asymptotic", or if the needed package is not installed, d.f. are set to Inf. (For backward compatibility, the user 
# may specify mode in lieu of lmer.df.)
# A by-product of the Kenward-Roger method is that the covariance matrix is adjusted using pbkrtest::vcovAdj(). 
# This can require considerable computation; so to avoid that overhead, the user should opt for the Satterthwaite or 
# asymptotic method

pp_emm_Li <- emmeans(mm2_prepost_Li, 'pre_post', type = 'response')
pp_emm_Tong <- emmeans(mm2_prepost_Tong, 'pre_post', type = 'response')


o1 <- ggplot(summary(pp_emm_Tong), aes(pre_post, response)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.1, size = 0.5) +
  annotate("text", x=1, y=1.044, label= "n=7229") +
  annotate("text", x=2, y=1.044, label= "n=5208") +
  annotate("text", x=3, y=1.044, label= "n=5354") +
  annotate("text", x=4, y=1.044, label= "n=3130") +
  annotate("text", x=5, y=1.044, label= "n=722") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=1.04, label= "a") +
  annotate("text", x=2, y=1.04, label= "bd") +
  annotate("text", x=3, y=1.04, label= "c") +
  annotate("text", x=4, y=1.04, label= "ab") +
  annotate("text", x=5, y=1.04, label= "cd") +
  theme_classic() +
  labs(x = "Period",
       y = "Resilience",
       title = "Isbell index")

o2 <- ggplot(summary(pp_emm_Li), aes(pre_post, response)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.1, size = 0.5) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=1.53, label= "a") +
  annotate("text", x=2, y=1.53, label= "a") +
  annotate("text", x=3, y=1.53, label= "a") +
  annotate("text", x=4, y=1.53, label= "a") +
  annotate("text", x=5, y=1.53, label= "a") +
  theme_classic() +
  labs(x = "Period",
       y = "Resilience",
       title = "Isbell index")

# uses library('gridExtra')
grid.arrange(o1, o2, nrow = 2)
# 


pp_emm_Li   <- emmeans(mm0_prepost_Li, 'pre_post', type = 'response')
pp_emm_Tong <- emmeans(mm0_prepost_Tong, 'pre_post', type = 'response')

o3 <- ggplot(summary(pp_emm_Tong), aes(pre_post, response)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.1, size = 0.5) +
  annotate("text", x=1, y=1.044, label= "n=7229") +
  annotate("text", x=2, y=1.044, label= "n=5208") +
  annotate("text", x=3, y=1.044, label= "n=5354") +
  annotate("text", x=4, y=1.044, label= "n=3130") +
  annotate("text", x=5, y=1.044, label= "n=722") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=1.04, label= "a") +
  annotate("text", x=2, y=1.04, label= "b") +
  annotate("text", x=3, y=1.04, label= "c") +
  annotate("text", x=4, y=1.04, label= "ab") +
  annotate("text", x=5, y=1.04, label= "cd") +
  theme_classic() +
  labs(x = "Period",
       y = "Resilience",
       title = "Isbell index")

o4 <- ggplot(summary(pp_emm_Li), aes(pre_post, response)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, 
                    ymax = asymp.UCL), 
                width = 0.1, size = 0.5) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=1.53, label= "ab") +
  annotate("text", x=2, y=1.53, label= "a") +
  annotate("text", x=3, y=1.53, label= "b") +
  annotate("text", x=4, y=1.53, label= "ab") +
  annotate("text", x=5, y=1.53, label= "ab") +
  theme_classic() +
  labs(x = "Period",
       y = "Resilience",
       title = "Isbell index")

# uses library('gridExtra')
grid.arrange(o3, o4, nrow = 2)
# 

# trade-offs Rt versus Rs in Tong's models
# amazing; no trade-off found with Lloret's indices
# par(mfrow=c(3,2))
# plot(p1$data$predicted, p3$data$predicted)
# abline(lm(p3$data$predicted ~ p1$data$predicted))
# R2 <- summary(lm(p3$data$predicted ~ p1$data$predicted))$ r.squared
# text(labels = format(round(R2, 2), nsmall = 2), x = 0.7, y = 1.02)
# 
# plot(p2$data$predicted, p4$data$predicted)
# abline(lm(p4$data$predicted ~ p2$data$predicted))
# R2 <- summary(lm(p4$data$predicted ~ p2$data$predicted))$ r.squared
# text(labels = format(round(R2, 2), nsmall = 2), x = 9, y = 1.02)
# 
# plot(p5$data$predicted, p7$data$predicted)
# abline(lm(p7$data$predicted ~ p5$data$predicted))
# R2 <- summary(lm(p7$data$predicted ~ p5$data$predicted))$ r.squared
# text(labels = format(round(R2, 2), nsmall = 2), x = 0.7, y = 1.02)
# 
# plot(p6$data$predicted, p8$data$predicted)
# abline(lm(p8$data$predicted ~ p6$data$predicted))
# R2 <- summary(lm(p8$data$predicted ~ p6$data$predicted))$ r.squared
# text(labels = format(round(R2, 2), nsmall = 2), x = 9, y = 1.02)
# 
# plot(p9$data$predicted, p11$data$predicted)
# abline(lm(p11$data$predicted ~ p9$data$predicted))
# R2 <- summary(lm(p11$data$predicted ~ p9$data$predicted))$ r.squared
# text(labels = format(round(R2, 2), nsmall = 2), x = 0.7, y = 1.02)
# 
# plot(p10$data$predicted, p12$data$predicted)
# abline(lm(p12$data$predicted ~ p10$data$predicted))
# R2 <- summary(lm(p12$data$predicted ~ p10$data$predicted))$ r.squared
# text(labels = format(round(R2, 2), nsmall = 2), x = 9, y = 1.02)








########################【【【【【Testing Li et al / only CONIFERS w/ Rt<1 ��】】】】###############

db1.28_Gym_lowRt <- subset(db1.28_Gymn, Resistance < 1)

ggplot(db1.28_Gymn, aes(log(Li_Rs), fill = pre_post)) + geom_density(alpha = 0.2) + xlim(-8, 8)
ggplot(db1.28_Gym_lowRt, aes(log(Li_Rs), fill = pre_post)) + geom_density(alpha = 0.2) + xlim(-8, 8)

# systematic comparison of distributions, prior to and following censoring of Resistance <1
# it does not as if this changed the distribution of Li_Rt and Li_Rs at all
with(db1.28_Gym_lowRt, hist(Resistance))
with(db1.28_Gymn, hist(Resistance))
with(db1.28_Gym_lowRt, hist(log(Li_Rt)))
with(db1.28_Gymn, hist(log(Li_Rt)))
with(db1.28_Gym_lowRt, hist(log(Li_Rs)))
with(db1.28_Gymn, hist(log(Li_Rs)))



#part 3: Li indices with factor pre/post1950
# cannot have variable slope at both levels; cannot even have intercept for site or slope for species
# can only have intercept for species; lmer converges for site intercept, but stdev=e-23!!
mm0_prepost_Li <- lmer(log(Li_Rs) ~ pre_post +
                         (1 | Species_code), 
                       control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
                       na.action = na.omit, data = db1.28_Gym_lowRt)


# cannot have variable slope at either levels; 
# this one plots Rs without adding unity; much better normality of residuals
mm2_prepost_Li <- lmer(log(Li_Rs) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season +
                         log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                         slope_TMPxYear + slope_PRExYear +  pre_post * log_Max_age +
                         pre_post* TMP_Annual_fixed + 
                         (1|Species_code/site), na.action = na.omit, data = db1.28_Gym_lowRt, REML=TRUE)

# this one plots Rs+1 following Li et al; unclear why they do that
pip1 <- lmer(log(Li_Rs+1) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season +
               log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
               slope_TMPxYear + slope_PRExYear +  pre_post * log_Max_age +
               pre_post* TMP_Annual_fixed +
               (1|Species_code/site), 
             control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
             na.action = na.omit, data = db1.28_Gym_lowRt, REML=TRUE)


anova(mm0_prepost_Li, type='I') # P = 0.1562 for pre_post; , type='I' is sequential
anova(mm2_prepost_Li, type='I') # P = 0.037792   for pre_post
anova(pip1, type='I') # P = 0.0116420  for pre_post; , type='I' is sequential

summary(mm0_prepost_Li)
summary(mm2_prepost_Li)
summary(pip1)

r.squaredGLMM(mm0_prepost_Li) # R2m = 0.0005200253  , R2c = 0.01086297
r.squaredGLMM(mm2_prepost_Li) # R2m = 0.04500534 ,   R2c = 0.04936497
r.squaredGLMM(pip1) # R2m = 0.04530806 ,   R2c = 0.0588918

# checks model quality for mixed models using library(performance) (needs lme4)
# check_model(mm0_prepost_Li)
check_model(mm0_prepost_Li)
chks01 <- check_model(pip1)
cbind(chks01$VIF$x, format(chks01$VIF$y, digits=2))
chks02 <- check_model(mm2_prepost_Li)
cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))

# model_performance(mm0_prepost_Li)
model_performance(mm0_prepost_Li) # AIC=64005.30 , R2c=7.00e-03, R2m=0, ICC=7.00e-03, RMSE=1.50
model_performance(mm2_prepost_Li) # AIC=62427.03 , R2c=0.02, R2m=0.02, ICC=4.00e-03, RMSE=1.49
model_performance(pip1) # AIC=62427.03 , R2c=0.02, R2m=0.02, ICC=4.00e-03, RMSE=1.49

check_singularity(mm2_prepost_Li) # FALSE
check_singularity(mm0_prepost_Li) # FALSE
check_singularity(pip1) # FALSE

#____________________________________________part 4: Tong's models


# variable slope much better than variable intercept / cannot have both variable
mm0_prepost_Tong <- lmer(log(Resilience) ~ pre_post +
                           (1|Species_code), 
                         control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
                         na.action = na.omit, data = db1.28_Gym_lowRt, REML = FALSE)


# pre_post : continentality_fixed and  pre_post : log_PRE_Annual_fixed highly multicollinear
mm2_prepost_Tong <- lmer(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
                           log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                           slope_TMPxYear + slope_PRExYear + 
                           pre_post * log_Max_age + pre_post * continentality_fixed + pre_post * TMP_Annual_fixed + 
                           (1|Species_code/site), na.action = na.omit, data = db1.28_Gym_lowRt, REML = FALSE)

# definitely the best !!! However, difficult to interpret fixed pre_post trends with variable slope!!
# to get convergence I need to get rid of fixed interactions of pre_post
pip4 <- lmer(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
               log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
               slope_TMPxYear + slope_PRExYear + 
               # pre_post * log_Max_age + pre_post * continentality_fixed + pre_post * TMP_Annual_fixed + 
               (pre_post|Species_code) + (1|Species_code:site), 
             control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
             na.action = na.omit, data = db1.28_Gym_lowRt, REML = FALSE)


anova(mm0_prepost_Tong, type='I') # F=35.135   P <.0001 for pre_post
anova(mm2_prepost_Tong, type='I') # F= 39.1825  P <.0001 for pre_post
anova(pip4, type='I') # F= 1.4197   P = ns for pre_post !!! very nice

# VERY INTERESTING!! Hence, suggestion is that time effect disappears when 
# I include variable random slopes by species??

summary(mm00_prepost_Tong)
summary(mm02_prepost_Tong)
summary(pip4)

# r.squaredGLMM(mm0_prepost_Tong)
r.squaredGLMM(mm0_prepost_Tong) # R2m = 0.005691636  ;  R2c = 0.01458711
r.squaredGLMM(mm2_prepost_Tong) # R2m = 0.4043942  ;   R2c = 0.4438393
r.squaredGLMM(pip4)     # R2m = 0.4204657    ;   R2c = 0.5104287 !!!

# check_model(mm0_prepost_Tong)
check_model(mm0_prepost_Tong)
chks01 <- check_model(mm1_prepost_Tong)
chks02 <- check_model(mm2_prepost_Tong)

cbind(chks01$VIF$x, format(chks01$VIF$y, digits=2))
cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))
check_model(pip4)

# model_performance(mm0_prepost_Tong)
model_performance(mm0_prepost_Tong) # AIC=-481.13  , R2c=0.01, R2m=6.00e-3, ICC=0.01, RMSE=0.24
model_performance(mm2_prepost_Tong) # AIC=-9159.35 , R2c=0.44, R2m=0.40,    ICC=0.07, RMSE=0.18
model_performance(pip4)             # AIC=-10097.14, R2c=0.51, R2m=0.42,    ICC=0.15, RMSE=0.17

check_singularity(mm2_prepost_Tong) # FALSE
check_singularity(mm0_prepost_Tong) # FALSE
check_singularity(pip4) # TRUE
# 


pp_emm <- emmeans(mm0_prepost_Li, 'pre_post')
pairs(pp_emm, type='response')

pp_emm <- emmeans(mm2_prepost_Li, 'pre_post')
pairs(pp_emm, type='response')



pp_emm <- emmeans(mm0_prepost_Tong, 'pre_post')
pairs(pp_emm, type='response')

pp_emm <- emmeans(mm2_prepost_Tong, 'pre_post')
pairs(pp_emm, type='response')




# # type='pred' gives marginal effect of pre_post (i.e., the partial derivative, which includes interactions)

p33 <- plot_model(mm0_prepost_Tong, type = 'pred', terms = 'pre_post', title='Lloret Resilience',
                  axis.title = c('Period', 'Resilience'))
p44 <- plot_model(mm0_prepost_Li, type = 'pred', terms = 'pre_post', title='Isbell Resilience',
                  axis.title = c('Period', 'Resilience')) 
p3 <- p33 +
  annotate("text", x=1, y=1.044, label= "n=7229") +
  annotate("text", x=2, y=1.044, label= "n=5208") +
  annotate("text", x=3, y=1.044, label= "n=5354") +
  annotate("text", x=4, y=1.044, label= "n=3130") +
  annotate("text", x=5, y=1.044, label= "n=722") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=1.03, label= "a") +
  annotate("text", x=2, y=1.03, label= "be") +
  annotate("text", x=3, y=1.03, label= "c") +
  annotate("text", x=4, y=1.03, label= "a") +
  annotate("text", x=5, y=1.03, label= "bcd") +
  theme_classic()
p4 <- p44 +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=1.38, label= "a") +
  annotate("text", x=2, y=1.38, label= "a") +
  annotate("text", x=3, y=1.38, label= "a") +
  annotate("text", x=4, y=1.38, label= "a") +
  annotate("text", x=5, y=1.38, label= "a") +
  theme_classic() 

# uses library('gridExtra')
grid.arrange(p3, p4, nrow = 2)

with(db1.28_Gym_lowRt, length(Li_Rs[which(pre_post=='02_29')])) # 4480
with(db1.28_Gym_lowRt, length(Li_Rs[which(pre_post=='30_49')])) # 3009
with(db1.28_Gym_lowRt, length(Li_Rs[which(pre_post=='50_69')])) # 3135
with(db1.28_Gym_lowRt, length(Li_Rs[which(pre_post=='70_89')])) # 1769
with(db1.28_Gym_lowRt, length(Li_Rs[which(pre_post=='90_08')])) # 475



p11 <- plot_model(mm2_prepost_Tong, type = 'pred', terms = 'pre_post', title='Lloret Resilience',
                  axis.title = c('Period', 'Resilience'), show.values = T)
p12 <- plot_model(mm2_prepost_Li, type = 'pred', terms = 'pre_post', title='Isbell Resilience',
                  axis.title = c('Period', 'Resilience'))  
p1 <- p11 +
  annotate("text", x=1, y=1.044, label= "n=4480") +
  annotate("text", x=2, y=1.044, label= "n=3009") +
  annotate("text", x=3, y=1.044, label= "n=3135") +
  annotate("text", x=4, y=1.044, label= "n=1769") +
  annotate("text", x=5, y=1.044, label= "n=475") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=1.035, label= "a") +
  annotate("text", x=2, y=1.035, label= "b") +
  annotate("text", x=3, y=1.035, label= "c") +
  annotate("text", x=4, y=1.035, label= "bd") +
  annotate("text", x=5, y=1.035, label= "ce") +
  theme_classic()
p2 <- p12 +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=1.8, label= "a") +
  annotate("text", x=2, y=1.8, label= "b") +
  annotate("text", x=3, y=1.8, label= "ab") +
  annotate("text", x=4, y=1.8, label= "ab") +
  annotate("text", x=5, y=1.8, label= "ab") +
  theme_classic() 

# uses library('gridExtra')
grid.arrange(p1, p2, nrow = 2)


####################【【【【【Testing Li et al /  CONIFERS 5-droughts SUBSET ��】】】】 ####
# starts from db1.28_Gym_lowRt database and
# subset of sites with at least one drought for each of the periods:1950-1969, 1970-1989 and 1990-2009

# Part 1: site selection according to Li et al (2020)

site_lst <- unique(as.character(db1.28_Gymn$site))
site_lst <- as.data.frame(site_lst)
site_lst$select_1 <- rep(NA, length(site_lst[,1]))
site_lst$select_2 <- rep(NA, length(site_lst[,1]))


for(i in 1:length(unique(db1.28_Gymn$site))) {
  # if there was a drought in every one of the three periods, then OK, otherwise NO
  pippo <- subset(db1.28_Gymn, site == site_lst[i,1], select = c('site', 'year_Li', 'TRI'))
  if( any(data.table::between(lower=(1950), upper=(1969), pippo$year_Li) ) & 
      any(data.table::between(lower=(1970), upper=(1989), pippo$year_Li) ) & 
      any(data.table::between(lower=(1990), upper=(2009), pippo$year_Li) ) )
  { site_lst$select_1[i] <- 'OK' } else {site_lst$select_1[i] <- 'NO' }
  # if first criterion OK, then if n. droughts>4 in entire period, then OK, otherwise NO 
  if( site_lst$select_1[i]== 'OK'  & 
    sum(data.table::between(lower=(1950), upper=(2009), pippo$year_Li))>4)
    {site_lst$select_2[i] <- 'OK' } else {site_lst$select_2[i] <- 'NO' }
    
    print(i)
}

lrg_grp  <- site_lst[site_lst$select_1=='OK',]  # 271 sites (or 260 w/ 10 p.le) (186 with Rt<1)
smll_grp <- site_lst[site_lst$select_1=='OK' & site_lst$select_2=='OK',] # 222 sites (or 205  w/ 10 p.le) (132 with Rt<1)

db_lrg_grp <-  db1.28_Gymn[db1.28_Gymn$site %in% lrg_grp$site_lst,]
db_smll_grp <- db1.28_Gymn[db1.28_Gymn$site %in% smll_grp$site_lst,]



#part 2: One drought every 20 years for three periods____________________________________________

# cannot have a random factor here (singular fit)
mm00_prepost_Li <- lmer(log(Li_Rs) ~ pre_post +
                          (1|Species_code), na.action = na.omit, data = db_lrg_grp, REML = FALSE)


# cannot have a random factor here (singular fit)
mm02_prepost_Li <- lmer(log(Li_Rs) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
                          log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                          slope_TMPxYear + slope_PRExYear + pre_post * log_Max_age + 
                          (1|Species_code), na.action = na.omit, data = db_lrg_grp, REML = TRUE)


anova(mm00_prepost_Li) # F=0.6, P = 0.2886 for pre_post
anova(mm02_prepost_Li) # F=0.8, P = 0.2952 for pre_post

summary(mm00_prepost_Li)
summary(mm02_prepost_Li)

# r.squaredGLMM(mm0_prepost_Li)
r.squaredGLMM(mm00_prepost_Li) # R2m = 0.0008806325 ,  R2c = 0.01449188
r.squaredGLMM(mm02_prepost_Li) # R2m = 0.03883429 ,  R2c = 0.04353525

# checks model quality for mixed models using library(performance) (needs lme4)
check_model(mm00_prepost_Li)
chks01 <- check_model(mm01_prepost_Li)
cbind(chks01$VIF$x, format(chks01$VIF$y, digits=2))
chks02 <- check_model(mm02_prepost_Li)
cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))

# model_performance(mm0_prepost_Li)
model_performance(mm00_prepost_Li) # AIC=10655.10 , R2c=0.01 , R2m=1.00e-03, ICC=0.01 , RMSE=1.47
model_performance(mm02_prepost_Li) # AIC=10482.99 , R2c=0.04 , R2m=0.04 , ICC=5.00e-3, RMSE=1.46

check_singularity(mm01_prepost_Li) # FALSE
check_singularity(mm02_prepost_Li) # FALSE
check_singularity(mm00_prepost_Li) # FALSE



#part 3: Tong's models

# variable slope not possible; only variable intercept possible; however can have spp/site
mm00_prepost_Tong <- lmer(log(Resilience) ~ pre_post +
                            (1|Species_code/site), na.action = na.omit, data = db_lrg_grp, REML = FALSE)


# pre_post : continentality_fixed and  pre_post : log_PRE_Annual_fixed highly multicollinear
mm02_prepost_Tong <- lmer(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
                            log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                            slope_TMPxYear + slope_PRExYear + 
                            pre_post * log_Max_age + pre_post * TMP_Annual_fixed + 
                            (1|Species_code/site), na.action = na.omit, data = db_lrg_grp, REML = FALSE)

# does converge !!!
pip4 <- lmer(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
               log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
               slope_TMPxYear + slope_PRExYear + 
               pre_post * log_Max_age + pre_post * TMP_Annual_fixed + 
               (pre_post|Species_code:site)+(1|Species_code),
             control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
             na.action = na.omit, data = db_lrg_grp, REML = FALSE)


pip5 <- lme(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
              log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
              slope_TMPxYear + slope_PRExYear + 
              pre_post * log_Max_age + pre_post * TMP_Annual_fixed, 
            random = list(Species_code = ~ 1, site = ~ pre_post), control = list(opt = "optim"),
            na.action = na.omit, data = db_lrg_grp, method= 'ML')


anova(mm00_prepost_Tong) # F=6.4, P <.0001 for pre_post
anova(mm02_prepost_Tong) # F=7.9, P <.0001 for pre_post
anova(pip5)              # F=4.9, P <.0001 for pre_post

summary(mm00_prepost_Tong)
summary(mm02_prepost_Tong)
summary(pip5)

r.squaredGLMM(mm00_prepost_Tong) # R2m = 0.008912227 ; R2c = 0.03093105
r.squaredGLMM(mm02_prepost_Tong) # R2m = 0.4247691  ;  R2c = 0.469373
r.squaredGLMM(pip5)              # R2m = 0.4766428 ;   R2c = 0.6411904 !!!
AIC(pip5)

# check_model(mm0_prepost_Tong)
check_model(mm00_prepost_Tong)
chks01 <- check_model(mm01_prepost_Tong)
chks02 <- check_model(mm02_prepost_Tong)

cbind(chks01$VIF$x, format(chks01$VIF$y, digits=2))
cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))
check_model(pip4)

# model_performance(mm0_prepost_Tong)
model_performance(mm00_prepost_Tong) # AIC=-567.21,  R2c=0.03, R2m=9.00e-3, ICC=0.02, RMSE=0.21
model_performance(mm02_prepost_Tong) # AIC=-2039.13 , R2c=0.47, R2m=0.42,    ICC=0.08, RMSE=0.16
model_performance(pip4)             #  AIC=-2182.5 , R2c=0.64, R2m=0.48,    ICC=0.31, RMSE=0.13

check_singularity(mm02_prepost_Tong) # FALSE
check_singularity(mm00_prepost_Tong) # FALSE
check_singularity(pip4)              # FALSE
# 




p11 <- plot_model(mm02_prepost_Tong, type = 'pred', terms = 'pre_post', title='Lloret Resilience',
                  axis.title = c('Period', 'Resilience'), show.values = T)
p12 <- plot_model(mm02_prepost_Li, type = 'pred', terms = 'pre_post', title='Isbell Resilience',
                  axis.title = c('Period', 'Resilience'))  
p1 <- p11 +
  annotate("text", x=1, y=1.044, label= "n=710") +
  annotate("text", x=2, y=1.044, label= "n=486") +
  annotate("text", x=3, y=1.044, label= "n=675") +
  annotate("text", x=4, y=1.044, label= "n=568") +
  annotate("text", x=5, y=1.044, label= "n=498") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  theme_classic()
p2 <- p12 +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  theme_classic() 

# uses library('gridExtra')
grid.arrange(p1, p2, nrow = 2)




#part 4: 5 droughts in 50 years_____________________________________________________

# cannot have a random factor here (singular fit)
mm00_prepost_Li <- lmer(log(Li_Rs) ~ pre_post +
                          (1|Species_code), na.action = na.omit, data = db_smll_grp, REML = FALSE)

# cannot have a random factor here (singular fit)
mm01_prepost_Li <- lmer(log(Li_Rs) ~ pre_post + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + 
                          continentality_fixed + 
                          pre_post * log_Max_age + 
                          (1|Species_code), na.action = na.omit, data = db_smll_grp, REML = TRUE)

# cannot have a random factor here (singular fit)
mm02_prepost_Li <- lmer(log(Li_Rs) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
                          log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                          slope_TMPxYear + slope_PRExYear +  
                          (1|Species_code), 
                        na.action = na.omit, data = db_smll_grp, REML = TRUE)


anova(mm00_prepost_Li) # F=0.48, P = 0.2886 for pre_post
anova(mm01_prepost_Li) # F=0.61, P = 0.2966 for pre_post
anova(mm02_prepost_Li) # F=0.65, P = 0.2952 for pre_post

summary(mm00_prepost_Li) # no random variance at either site level
summary(mm01_prepost_Li) # no random variance at either site level
summary(mm02_prepost_Li) # no random variance at either site level

r.squaredGLMM(mm00_prepost_Li) # R2m = 0.0007878289 ,  R2c = 0.01128688
r.squaredGLMM(mm01_prepost_Li) # R2m = 0.01793089 ,  R2c = 0.01893975 
r.squaredGLMM(mm02_prepost_Li) # R2m = 0.03067418 ,  R2c = 0.03237119  

# checks model quality for mixed models using library(performance) (needs lme4)
check_model(mm00_prepost_Li)
chks01 <- check_model(mm01_prepost_Li)
cbind(chks01$VIF$x, format(chks01$VIF$y, digits=2))
chks02 <- check_model(mm02_prepost_Li)
cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))

# model_performance(mm0_prepost_Li)
model_performance(mm00_prepost_Li) # AIC=8821.88 , R2c=0.01, R2m=1.00e-03, ICC=0.01, RMSE=1.47
model_performance(mm01_prepost_Li) # AIC=8860.71 , R2c=0.02 , R2m=0.02 , ICC=1.00e-03, RMSE=1.47
model_performance(mm02_prepost_Li) # AIC=8707.93 , R2c=0.03 , R2m=0.03  , ICC=2.00e-03, RMSE=1.47

check_singularity(mm01_prepost_Li) # FALSE
check_singularity(mm02_prepost_Li) # FALSE
check_singularity(mm00_prepost_Li) # FALSE



#part 5: Tong's models

# variable slope not possible; only variable intercept possible and only at spp level
mm00_prepost_Tong <- lmer(log(Resilience) ~ pre_post +
                            (1|Species_code), 
                          control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
                          na.action = na.omit, data = db_smll_grp, REML = FALSE)

# not working
pip1 <- lmer(log(Resilience) ~ pre_post +
               (pre_post+0|Species_code/site), na.action = na.omit, data = db_smll_grp, REML = FALSE)

# not working
pip2 <- lmer(log(Resilience) ~ pre_post +
               (pre_post|Species_code), na.action = na.omit, data = db_smll_grp, REML = FALSE)

# not working
pip3 <- lmer(log(Resilience) ~ pre_post +
               (1|Species_code) + (pre_post|Species_code/site), na.action = na.omit, data = db_smll_grp, REML = FALSE)

# pre_post : continentality_fixed and  pre_post : log_PRE_Annual_fixed highly multicollinear
mm01_prepost_Tong <- lmer(log(Resilience) ~ pre_post + log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                            pre_post * log_Max_age + pre_post * TMP_Annual_fixed + 
                            (1|Species_code/site), na.action = na.omit, data = db_smll_grp, REML = FALSE)

# pre_post : continentality_fixed and  pre_post : log_PRE_Annual_fixed highly multicollinear
mm02_prepost_Tong <- lmer(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
                            log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                            slope_TMPxYear + slope_PRExYear + 
                            pre_post * log_Max_age + pre_post * TMP_Annual_fixed + 
                            (1|Species_code/site), 
                          control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')), 
                          na.action = na.omit, data = db_smll_grp, REML = FALSE)

# definitely the best !!!
pip4 <- lmer(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
               log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
               slope_TMPxYear + slope_PRExYear + 
               pre_post * log_Max_age + pre_post * TMP_Annual_fixed + 
               (pre_post|Species_code:site)+ (1|Species_code), 
             na.action = na.omit, data = db_smll_grp, REML = TRUE,
             control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb')))

# does converge !!!
pip5 <- lme(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
              log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
              slope_TMPxYear + slope_PRExYear + 
              pre_post * log_Max_age + pre_post * TMP_Annual_fixed, 
            random = list(Species_code = ~ 1, site = ~ pre_post), control = list(opt = "optim"),
            na.action = na.omit, data = db_smll_grp, method= 'ML')

anova(mm00_prepost_Tong) # F=6.4, P <.0001 for pre_post
anova(mm01_prepost_Tong) # F=10.9, P <.0001 for pre_post
anova(mm02_prepost_Tong) # F=7.5, P <.0001 for pre_post
anova(pip4)              # F=4-0, P <.0001 for pre_post
lsmeans(pip4)

summary(mm00_prepost_Tong)
summary(mm01_prepost_Tong)
summary(mm02_prepost_Tong)
summary(pip4)
summary(pip5)

intervals(pip5)

r.squaredGLMM(mm00_prepost_Tong) # R2m = 0.01082628  ; R2c = 0.02631512
r.squaredGLMM(mm01_prepost_Tong) # R2m = 0.3739672  ;   R2c = 0.4134006
r.squaredGLMM(mm02_prepost_Tong) # R2m = 0.4160626   ;  R2c = 0.4637758
r.squaredGLMM(pip4)              # R2m = 0.4696482  ;   R2c = 0.6664186 !!!
r.squaredGLMM(pip5)              # R2m = 0.4733085   ;   R2c = 0.6638215 !!!

# check_model(mm0_prepost_Tong)
check_model(mm00_prepost_Tong)
chks01 <- check_model(mm01_prepost_Tong)
chks02 <- check_model(mm02_prepost_Tong)

cbind(chks01$VIF$x, format(chks01$VIF$y, digits=2))
cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))
check_model(pip4)

# model_performance(mm0_prepost_Tong)
model_performance(mm00_prepost_Tong) # AIC=-467.61 ,  R2c=0.03, R2m=0.01, ICC=0.02, RMSE=0.22
model_performance(mm01_prepost_Tong) # AIC=-1518.74 , R2c=0.41, R2m=0.37,    ICC=0.06, RMSE=0.17
model_performance(mm02_prepost_Tong) # AIC=-1648.15 , R2c=0.46, R2m=0.42,    ICC=0.08, RMSE=0.16
model_performance(pip4)              # AIC=-1594.61 , R2c=0.67, R2m=0.47,    ICC=0.37, RMSE=0.12

check_singularity(mm01_prepost_Tong) # FALSE
check_singularity(mm02_prepost_Tong) # FALSE
check_singularity(mm00_prepost_Tong) # FALSE
check_singularity(pip4)              # FALSE
# 




p11 <- plot_model(mm02_prepost_Tong, type = 'pred', terms = 'pre_post', title='Lloret Resilience',
                  axis.title = c('Period', 'Resilience'), show.values = T)
p12 <- plot_model(mm02_prepost_Li, type = 'pred', terms = 'pre_post', title='Isbell Resilience',
                  axis.title = c('Period', 'Resilience'))  
p1 <- p11 +
  annotate("text", x=1, y=1.049, label= "n=517") +
  annotate("text", x=2, y=1.049, label= "n=356") +
  annotate("text", x=3, y=1.049, label= "n=610") +
  annotate("text", x=4, y=1.049, label= "n=512") +
  annotate("text", x=5, y=1.049, label= "n=439") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  theme_classic()
p2 <- p12 +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  theme_classic() 

# uses library('gridExtra')
grid.arrange(p1, p2, nrow = 2)


#____________________ testing distributions after 1950_______________####

db_past_1950 <- with(db_lrg_grp, db_lrg_grp[which(pre_post== '50_69' | pre_post== '70_89' | pre_post== '70_89' | pre_post== '90_08'),])
db_past_1950 <- droplevels(db_past_1950)



ggplot(db_past_1950, aes(log(Li_Rs+1), fill = pre_post)) + geom_density(alpha = 0.2) + xlim(0, 4)


# One-sample Kolmogorov-Smirnov test
# 
ks.test(log(db_past_1950$Li_Rs), "pnorm", mean=mean(log(db_past_1950$Li_Rs)), sd=sd(log(db_past_1950$Li_Rs)))


# data:  log(db_lrg_grp$Li_Rs)
# D = 0.045344, p-value = 4.957e-05
# alternative hypothesis: two-sided

shapiro.test(log(db_past_1950$Li_Rs))

# Shapiro-Wilk normality test
# 
# data:  log(db_lrg_grp$Li_Rs)
# W = 0.97986, p-value < 2.2e-16



# Using Kolmogorov-Sirnov which is not very powerful to detect shifts in means
# does no find significance
with(db_past_1950, ks.test(Li_Rs[which(pre_post=='50_69')],
                                 Li_Rs[which(pre_post=='70_89')], alternative = 'two.sided'))

with(db_past_1950, ks.test(Li_Rs[which(pre_post=='50_69')],
                                 Li_Rs[which(pre_post=='90_08')], alternative = 'two.sided'))

with(db_past_1950, ks.test(Li_Rs[which(pre_post=='70_89')],
                                 Li_Rs[which(pre_post=='90_08')], alternative = 'two.sided'))

# using more powerful Kruskal Wallis, either as asymptotic or as approximate (permutations in library(coin))
kruskal.test(Li_Rs ~ pre_post, data = db_past_1950)
kruskal.test(Resilience ~ pre_post, data = db_past_1950)

# data:  value by Periods
# Kruskal-Wallis chi-squared = 8.4312, df = 2, p-value = 0.01476

# Using permutation tests using library(coin)
it <- with(db_past_1950, independence_test(Li_Rs ~ pre_post, data = df1,
                                  ytrafo = function(data)
                                    trafo(data, numeric_trafo = rank_trafo),
                                  teststat = "quadratic"))
it
statistic(it, type = "standardized")
it <- with(db_past_1950, independence_test(Resilience ~ pre_post, data = df1,
                                  ytrafo = function(data)
                                    trafo(data, numeric_trafo = rank_trafo),
                                  teststat = "quadratic"))
it

statistic(it, type = "standardized")
# same result as Kruskal-Wallis: finds difference in Resilience, but not Li_Rs, with usual trend!

# Li_Rs
# 50_69  1.0434152
# 70_89 -0.8873709
# 90_08 -0.2550417


# Resilience
# 50_69 -0.0008066007
# 70_89  3.8325628576
# 90_08 -4.1902222572


# Asymptotic General Independence Test
# 
# data:  value by Periods (Rs_50_69, Rs_70_89, Rs_90_09)
# chi-squared = 8.4312, df = 2, p-value = 0.01476


# Exact null distribution of the Kruskal-Wallis test is approximated
# by Monte Carlo resampling using 10000 replicates via

kt <- kruskal_test(Li_Rs ~ pre_post, data = db_past_1950,
                   distribution = approximate(nresample = 10000))
kt
pvalue(kt)
kt <- kruskal_test(Resilience ~ pre_post, data = db_past_1950,
                   distribution = approximate(nresample = 10000))
kt
pvalue(kt)


ggplot(df1, aes(log(Li_Rs+1), fill = pre_post)) + geom_density(alpha = 0.2) + xlim(0, 4)




####################_____________Li et al import 332 sites data _____________________________####


#_____________________________________
# three elements: 1) distributions; 2) models with my data at those sites; 3) models with their data

# first issue: we looked at *.crn whereas they looked at *.rwl data (more datasets)
# when I run my models, I still get the same result as before on 139 of those 332 sites

#_____________________________________ Importing Shilong's data

setwd("D:/Dropbox/My Documents/science/projects/grants/China Tong Zheng/Response to Li")

Li_etal_Rssites <-data.frame(read.table("Li_etal_332_Rsvalues.csv", header=T, sep=","))
Li_etal_Rtsites <-data.frame(read.table("Li_etal_332_Rtvalues.csv", header=T, sep=","))


# ok, code below checks if sites match line by line

yes=0
for(i in 1:length(Li_etal_Rssites$sitename)) {

  if(Li_etal_Rssites$sitename[i]==Li_etal_Rtsites$sitename[i]) {yes=yes+1}
print(yes)
  }

Li_etal_sites <- merge(Li_etal_Rssites, Li_etal_Rtsites, 'sitename')
colnames(Li_etal_sites) <- c('site', 'Rs_50_69', 'Rs_70_89', 'Rs_90_09', 'Rt_50_69', 'Rt_70_89', 'Rt_90_09')

rm(Li_etal_Rssites)
rm(Li_etal_Rtsites)


setwd("D:/Dropbox/My Documents/science/projects/grants/China Tong Zheng/Tong Zheng folder/R script 20200621")


# uses library data.table to go from wide to long
Li_etal_sites_long <- melt(setDT(Li_etal_sites[,1:4]), id.vars = c("site"), variable.name = "Periods")


#____________________________________________________________________
# subsetting db1.28_Gymn to list in Li_etal_sites (we used *.crn instead of *.rwl files)
# end up with 138 sites instead of 332

Li_etal_sites$site <- as.character(Li_etal_sites$site)
Li_etal_sites$site <- paste(Li_etal_sites$site, '.crn', sep='')
db1.28_Gymn$site   <- as.character(db1.28_Gymn$site)

Li_etal_data <- as.data.frame(matrix(nrow = 0, ncol = ncol(db1.28_Gymn)))
for(i in 1:length(unique(Li_etal_sites$site))) {
  Site <- Li_etal_sites$site[i]
  if(length(db1.28_Gymn$site[db1.28_Gymn$site == Site])>0) {
    pippo <-   with(db1.28_Gymn,db1.28_Gymn[site == Site,])
    Li_etal_data <-   rbind(Li_etal_data, pippo)
  }
print(i)
  }

with(Li_etal_data, plot(Longitude, Latitude))


#____________________________________________________________________
# subsetting Li_etal_sites_long to after 1950; then 
# aggregating to one mean Rs value/site after 1950. Finally, 
# Merging data from Shilong to our data.


Lietal_past_1950 <- with(Li_etal_data, Li_etal_data[which(pre_post== '50_69' | pre_post== '70_89' | pre_post== '70_89' | pre_post== '90_08'),])
Lietal_past_1950 <- droplevels(Lietal_past_1950)

# Nesting -------------------------------------------------------------------
df1 <- tibble(Lietal_past_1950)

# uses dplyr to complete a data frame with missing combinations of data.
df2 <- complete(df1, site, pre_post)

#____________________________________________________________________
# Aggregating to one mean Rs value/site after 1950.


# uses dplyr
df3 <- df2 %>% group_by(site, pre_post) %>% summarise_if(is.numeric, meanNA)
df4 <- df2 %>% group_by(site) %>% summarise_if(is.numeric, meanNA) # gives 138 sites

spp_cds <- with(df1, aggregate(Species_code, by=list(site), FUN=unique))[,'x']

df3$Species_code <- rep(spp_cds, 1, each = 3)

df3[is.nan.data.frame(df3)] <- NA  # to eliminate all the NaNs in the df


#____________________________________________________________________
# finally merging the two (left_join)

# subsetting Li_etal_sites_long to our sites
Li_etal_sites_long2 <- Li_etal_sites_long %>% filter(site %in% df3$site)
colnames(Li_etal_sites_long2)[2] <- 'pre_post'
colnames(Li_etal_sites_long2)[3] <- 'Li_Rs_true'
levels(Li_etal_sites_long2$pre_post) <- c('50_69', '70_89', '90_08')

df3 <- df3[, !(names(df3) %in% 'Li_Rs_true')]

df3 <- left_join(df3, Li_etal_sites_long2)

# check two estimates of Li_Rs
with(df3, plot(log(Li_Rs), log(Li_Rs_true)))
abline(0,1, col='red')

with(df3, summary(lm(log(Li_Rs) ~ log(Li_Rs_true))))


#___________________________________________________________________________
# Using Kolmogorov-Sirnov which is not very powerful to detect shifts in means
# does no find significance
with(df3, ks.test(Li_Rs[which(pre_post=='50_69')],
                  Li_Rs[which(pre_post=='70_89')], alternative = 'two.sided'))

with(df3, ks.test(Li_Rs[which(pre_post=='50_69')],
                  Li_Rs[which(pre_post=='90_08')], alternative = 'two.sided'))

with(df3, ks.test(Li_Rs[which(pre_post=='70_89')],
                  Li_Rs[which(pre_post=='90_08')], alternative = 'two.sided'))

# using more powerful Kruskal Wallis, either as asymptotic or as approximate (permutations in library(coin))
kruskal.test(value ~ Periods, data = Li_etal_sites_long) # original 332 sites
kruskal.test(Li_Rs_true ~ pre_post, data = df3) # subset of 138 sites
kruskal.test(Li_Rs ~ pre_post, data = df3) # my calcs of Li_Rs
kruskal.test(Resilience ~ pre_post, data = df3) # sensu-Lloret Rs

# data:  value by Periods
# Kruskal-Wallis chi-squared = 8.4312, df = 2, p-value = 0.01476

# Using permutation tests using library(coin) on original data
it <- with(df3, independence_test(value ~ Periods, data = Li_etal_sites_long,
                                  ytrafo = function(data)
                                    trafo(data, numeric_trafo = rank_trafo),
                                  teststat = "quadratic"))
it
statistic(it, type = "standardized")

# Using permutation tests using library(coin) on my data
it <- with(df3, independence_test(Resilience ~ pre_post, data = df3,
                                  ytrafo = function(data)
                                    trafo(data, numeric_trafo = rank_trafo),
                                  teststat = "quadratic"))
it

statistic(it, type = "standardized")
# same result as Kruskal-Wallis: finds difference in Resilience, but not Li_Rs, with usual trend!

# Using permutation tests using library(coin) on subset of 138 sites
it <- with(df3, independence_test(Li_Rs_true ~ pre_post, data = df3,
                                  ytrafo = function(data)
                                    trafo(data, numeric_trafo = rank_trafo),
                                  teststat = "quadratic"))
it
statistic(it, type = "standardized")
# Li_Rs
# 50_69  1.846837
# 70_89 -1.572541
# 90_08 -0.329047


# Resilience
# 50_69 -1.082909
# 70_89  3.223093
# 90_08 -2.362185


# Asymptotic General Independence Test
# 
# data:  value by Periods (Rs_50_69, Rs_70_89, Rs_90_09)
# chi-squared = 8.4312, df = 2, p-value = 0.01476


# Exact null distribution of the Kruskal-Wallis test is approximated
# by Monte Carlo resampling using 10000 replicates via

kt <- with(df3, kruskal_test(Li_Rs_true ~ pre_post, data = df3,
                             distribution = approximate(nresample = 10000)))
kt
pvalue(kt)


kt <- with(df3, kruskal_test(Li_Rs ~ pre_post, data = df3,
                             distribution = approximate(nresample = 10000)))
kt
pvalue(kt)


kt <- with(df3, kruskal_test(Resilience ~ pre_post, data = df3,
                             distribution = approximate(nresample = 10000)))
kt
pvalue(kt)







#___________________________________________________________________________
#part 2: One drought every 20 years for three periods____________________________________________

# cannot have a random factor here (singular fit)
# mm00_prepost_Li <- lmer(log(Li_Rs) ~ pre_post +
#                           (1|Species_code), na.action = na.omit, data = df3, REML = FALSE)
# 
# 
# # cannot have a random factor here (singular fit)
# mm02_prepost_Li <- lmer(log(Li_Rs) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost +  
#                           log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
#                           pre_post * log_Max_age + pre_post * TMP_Annual_fixed + pre_post * continentality_fixed +
#                           (1|Species_code), na.action = na.omit, data = df3, REML = TRUE)

# cannot have a random factor here (singular fit)
mm00_prepost_Li <- lmer(log(Li_Rs_true) ~ pre_post +
                          (1|Species_code), na.action = na.omit, data = df3, REML = FALSE)

# mm02_prepost_Li <- lmer(log(Li_Rs) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
#                           log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
#                           pre_post * log_Max_age + pre_post * TMP_Annual_fixed + pre_post * continentality_fixed +
#                           (1|Species_code), na.action = na.omit, data = db_lrg_grp, REML = TRUE)
# 


# cannot have a random factor here (singular fit)
# note that I had to delete Season as it does not exist in this df
mm02_prepost_Li <- lmer(log(Li_Rs_true) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost +  
                          log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                          pre_post * log_Max_age + pre_post * TMP_Annual_fixed + pre_post * continentality_fixed +
                          (1|Species_code), na.action = na.omit, data = df3, REML = TRUE)



anova(mm00_prepost_Li, type = 'I') # F=0.4, P = 0.7555 for pre_post
anova(mm02_prepost_Li, type = 'I') # F=0.9, P = 0.4873355 for pre_post

summary(mm00_prepost_Li)
summary(mm02_prepost_Li)

# r.squaredGLMM(mm0_prepost_Li)
r.squaredGLMM(mm00_prepost_Li) # R2m = 0.001297039  ,  R2c = 0.01972575
r.squaredGLMM(mm02_prepost_Li) # R2m = 0.06957116 ,  R2c = 0.07489106

# checks model quality for mixed models using library(performance) (needs lme4)
check_model(mm00_prepost_Li)
chks02 <- check_model(mm02_prepost_Li)
cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))

# model_performance(mm0_prepost_Li)
model_performance(mm00_prepost_Li) # AIC=5378.22, R2c=0.02 , R2m=1.00e-03, ICC=0.02 , RMSE=1.49
model_performance(mm02_prepost_Li) # AIC=5299.14, R2c=0.07 , R2m=0.07 , ICC=6.00e-3, RMSE=1.47

check_singularity(mm02_prepost_Li) # FALSE
check_singularity(mm00_prepost_Li) # FALSE



#part 3: Tong's models

# variable slope not possible; only variable intercept possible; however can have spp/site
mm00_prepost_Tong <- lmer(log(Resilience) ~ pre_post +
                            (1|Species_code), na.action = na.omit, data = df3, REML = FALSE)


# pre_post : continentality_fixed and  pre_post : log_PRE_Annual_fixed highly multicollinear
mm02_prepost_Tong <- lmer(log(Resilience) ~ pre_post + log_RWI4 + SPEId + SPEIprev + SPEIpost  + 
                            log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                            pre_post * log_Max_age + pre_post * TMP_Annual_fixed + pre_post * continentality_fixed +
                            (1|Species_code/site), na.action = na.omit, data = df3, REML = FALSE)


anova(mm00_prepost_Tong, type = 'I') # F=7.22 , P <.0001 for pre_post
anova(mm02_prepost_Tong, type = 'I') # F=11.8061, P <.0001 for pre_post

summary(mm00_prepost_Tong)
summary(mm02_prepost_Tong)

r.squaredGLMM(mm00_prepost_Tong) # R2m = 0.01969529 ; R2c = 0.03971507
r.squaredGLMM(mm02_prepost_Tong) # R2m = 0.4545994  ; R2c = 0.5029409

check_model(mm00_prepost_Tong)
chks02 <- check_model(mm02_prepost_Tong)

cbind(chks02$VIF$x, format(chks02$VIF$y, digits=2))

model_performance(mm00_prepost_Tong) # AIC=-269.45,  R2c=0.04, R2m=0.02, ICC=0.02, RMSE=0.22
model_performance(mm02_prepost_Tong) # AIC=-1049.53, R2c=0.50, R2m=0.46,    ICC=0.09, RMSE=0.16

check_singularity(mm02_prepost_Tong) # FALSE
check_singularity(mm00_prepost_Tong) # FALSE
# 


pairs(emmeans(mm02_prepost_Li, 'pre_post', type = 'response'))
pairs(emmeans(mm02_prepost_Tong, 'pre_post', type = 'response'))

pp_emm_Li   <- emmeans(mm02_prepost_Li, 'pre_post', type = 'response')
pp_emm_Tong <- emmeans(mm02_prepost_Tong, 'pre_post', type = 'response')


p1 <- ggplot(summary(pp_emm_Tong), aes(pre_post, response)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, 
                    ymax = upper.CL), 
                width = 0.1, size = 0.5) +
  # coord_cartesian(ylim = c(0.85, 1.1)) +
  ylim(0.85, 1.1) +
  annotate("text", x=1, y=1.1, label= "n=430") +
  annotate("text", x=2, y=1.1, label= "n=222") +
  annotate("text", x=3, y=1.1, label= "n=355") +
  annotate("text", x=4, y=1.1, label= "n=239") +
  annotate("text", x=5, y=1.1, label= "n=221") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=1.07, label= "a") +
  annotate("text", x=2, y=1.07, label= "bc") +
  annotate("text", x=3, y=1.07, label= "c") +
  annotate("text", x=4, y=1.07, label= "ab") +
  annotate("text", x=5, y=1.07, label= "c") +
  theme_classic() +
  labs(x = "Period",
       y = "Resilience",
       title = "Lloret index")


p2 <- ggplot(summary(pp_emm_Li), aes(pre_post, response)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower.CL, 
                    ymax = upper.CL), 
                width = 0.1, size = 0.5) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate("text", x=1, y=2.2, label= "a") +
  annotate("text", x=2, y=2.2, label= "a") +
  annotate("text", x=3, y=2.2, label= "a") +
  annotate("text", x=4, y=2.2, label= "a") +
  annotate("text", x=5, y=2.2, label= "a") +
  theme_classic() +
  labs(x = "Period",
       y = "Resilience",
       title = "Isbell index")


# uses library('gridExtra')
grid.arrange(p1, p2, nrow = 2)





# testing distributions after 1950 in original data

ggplot(Li_etal_sites, aes(log(Rs_50_69+1))) + geom_density(alpha = 0.2) + xlim(0, 4) +
         geom_density(aes(log(Rs_70_89+1), col = 'red')) + 
         geom_density(aes(log(Rs_90_09+1), col = 'blue'))


# distributions in original dataset (Kolmogorov-Smirnov ONE-SAMPLE distribution test)
ks.test(log(Li_etal_sites$Rs_50_69), "pnorm", mean=mean(log(Li_etal_sites$Rs_50_69)), 
        sd=sd(log(Li_etal_sites$Rs_50_69)))
ks.test(log(Li_etal_sites$Rs_70_89), "pnorm", mean=mean(log(Li_etal_sites$Rs_70_89)), 
        sd=sd(log(Li_etal_sites$Rs_70_89)))
ks.test(log(Li_etal_sites$Rs_90_09), "pnorm", mean=mean(log(Li_etal_sites$Rs_90_09)), 
        sd=sd(log(Li_etal_sites$Rs_90_09)))
# all data
with(Li_etal_sites, ks.test(c(log(Rs_50_69), log(Rs_70_89), log(Rs_90_09)), 
         "pnorm", mean=mean(c(log(Rs_50_69), log(Rs_70_89), log(Rs_90_09))), 
       sd=sd(c(log(Rs_50_69), log(Rs_70_89), log(Rs_90_09)))))



# The Kruskal-Wallis test is the nonparametric rank-based equivalent to a one-way ANOVA. The Kruskal-Wallis is 
#  to the Wilcoxon-Mann-Whitney (or KS) two-sample test as one-way ANOVA is to a two-sample t-test. It is used 
#  when you want to test against the null that more than two groups have the same location.
# If you want to test specifically for a difference in means, neither the K-S test nor the W-M-W really does it 
#  (though with some additional assumptions the W-M-W is also a test for a difference in means).

#  The best way to test for a difference in means is probably to do a permutation test, as long as the distributions 
#  would be the same under the null (so if your alternative is a location-shift, you're basically assuming identical 
#  shapes apart from location). 

#  You could use a Smirnov test (a two sample K-S test) to test for any kind of difference between the two groups, 
#  but if your interest is a difference in means it's not a very powerful test.

# Yes, a Kruskal-Wallis applied to two samples would give the same result as a Wilcoxon-Mann-Whitney for a two-tailed
#  test (i.e. it doesn't allow a one-sided alternative, in the same way that a one-way ANOVA doesn't give you the 
#  one-sided alternative you can get with a two-sample t-test and a chi-square doesn't give the one-sided alternative
#  of a two-sample proportions test).


# Kruskal-Wallis: non-parametric equivalent of one-way ANOVA
Li_etal_sites_long$Periods <- as.character(Li_etal_sites_long$Periods)

# Li_etal_sites_long <- with(Li_etal_sites_long, Li_etal_sites_long[Periods=='Rs_70_89' | Periods=='Rs_90_09',])
# alternatively, tried this partitioning
# Li_etal_sites_long <- with(Li_etal_sites_long, Li_etal_sites_long[Periods=='Rs_50_69' | Periods=='Rs_70_89',])
# Li_etal_sites_long <- with(Li_etal_sites_long, Li_etal_sites_long[Periods=='Rs_50_69' | Periods=='Rs_90_09',])

kruskal.test(value ~ Periods, data = Li_etal_sites_long)

# data:  value by Periods
# Kruskal-Wallis chi-squared = 8.4312, df = 2, p-value = 0.01476

# Using permutation tests using library(coin)
it <- with(Li_etal_sites_long, independence_test(value ~ Periods, data = Li_etal_sites_long,
                        ytrafo = function(data)
                          trafo(data, numeric_trafo = rank_trafo),
                        teststat = "quadratic"))
it

statistic(it, type = "standardized")
# same result as Kruskal-Wallis

# Asymptotic General Independence Test
# 
# data:  value by Periods (Rs_50_69, Rs_70_89, Rs_90_09)
# chi-squared = 8.4312, df = 2, p-value = 0.01476


# Exact null distribution of the Kruskal-Wallis test is approximated
# by Monte Carlo resampling using 10000 replicates via

kt <- with(Li_etal_sites_long, kruskal_test(value ~ Periods, data = Li_etal_sites_long,
             distribution = approximate(nresample = 10000)))

pvalue(kt)
# same results as Kruskal-Wallis
# testing distribution after 1950 in merged dataset; don't understand why the use Rs+1 ???


# subset to after 1950 data
db_past_1950 <- with(Li_etal_data, Li_etal_data[which(pre_post== '50_69' | pre_post== '70_89' | pre_post== '70_89' | 
                                                        pre_post== '90_08'),])
db_past_1950 <- droplevels(db_past_1950)

# merged dataset
ks.test(log(db_past_1950$Li_Rs), "pnorm", mean=mean(log(db_past_1950$Li_Rs)), sd=sd(log(db_past_1950$Li_Rs)))


# One-sample Kolmogorov-Smirnov test
# 
# data:  log(Li_etal_data$Li_Rs)
# D = 0.045344, p-value = 4.957e-05
# alternative hypothesis: two-sided

shapiro.test(log(db_past_1950$Li_Rs))

# Shapiro-Wilk normality test
# 
# data:  log(Li_etal_data$Li_Rs)
# W = 0.97986, p-value < 2.2e-16









### GAMM models in CONIFERS ###  ########################
# gamms estimated as random effects as well as fixed effect 8from its own variance)
# hence, makes no sense to have year as fixed effect and s(year) as random effect
###
###




gamm0_Tong <- gamm4(log(Resistance) ~ s(year), 
                    random = ~ (1|Species_code), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)
viz <- getViz(gamm0_Tong$gam)
plot_Tong_m0 <- plot( sm(viz, 1))

gamm00_Tong <- gamm4(log(Resilience) ~ s(year), 
              random = ~ (1|Species_code), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)
# if random slope for year, then slope and intercept have r=1
# gives random var=0 for site within species
# hence only logical model is random = ~ (1|Species_code): AIC=-415.5766

# plot(gamm00_Tong$gam)
viz <- getViz(gamm00_Tong$gam)
# print(plot(viz, allTerms = T), pages = 1)

plot_Tong_m00 <- plot( sm(viz, 1))

with(plot_Tong_m00, plot(ggObj$data$x, exp(ggObj$data$y), ylim=c(0.8, 1.4), type='l', 
                        main="Lloret null model", xlab= "Year", ylab= 'Resilience',
                        cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, las = 1, lwd = 2, xaxt = 'n'))
axis(side = 1, at = c( -40, -20, 0, 20, 40),
     labels = c('1915', '1935', '1955', '1975', '1995'), cex.axis = 1.5)
with(plot_Tong_m00, lines(x=ggObj$data$x, y=exp(ggObj$data$y-5*ggObj$data$se), 
                         col='red', lwd = 2))
with(plot_Tong_m00, lines(x=ggObj$data$x, y=exp(ggObj$data$y+5*ggObj$data$se), 
                         col='red', lwd = 2))

anova(gamm00_Tong$gam)
anova(gamm00_Tong$mer)
summary(gamm00_Tong$gam)
summary(gamm00_Tong$mer)
model_performance(gamm00_Tong$mer) # R2c and R2m both = 0.39; ICC = 0; RMSE = 0.24
check_singularity(gamm00_Tong$mer)


gamm1_Tong <- gamm4(log(Resistance) ~ s(year) +
                      log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                      year : TMP_Annual_fixed + year : log_PRE_Annual_fixed + year : continentality_fixed + slope_TMPxYear + slope_PRExYear,
                    random = ~ (1|Species_code/site), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)
viz <- getViz(gamm1_Tong$gam)
plot_Tong_m1 <- plot( sm(viz, 1))


gamm01_Tong <- gamm4(log(Resilience) ~ s(year) +
                      log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                      year : TMP_Annual_fixed + year : log_PRE_Annual_fixed + year : continentality_fixed + slope_TMPxYear + slope_PRExYear,
                    random = ~ (1|Species_code/site), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)

anova(gamm01_Tong$gam)
anova(gamm01_Tong$mer)
summary(gamm01_Tong$gam)
summary(gamm01_Tong$mer)
model_performance(gamm01_Tong$mer) # R2c and R2m both = 0.47; ICC = 0; RMSE = 0.18
check_singularity(gamm01_Tong$mer)

# plot(gamm01_Tong$gam)

viz <- getViz(gamm01_Tong$gam)
# plot_Tong_m01 <- plot( viz, allTerms = T)
plot_Tong_m01 <- plot( sm(viz, 1))
with(plot_Tong_m01, plot(ggObj$data$x, exp(ggObj$data$y), ylim=c(0.7, 1.5), type='l', 
                        main="Lloret climate model", xlab= "Year", ylab= 'Resilience',
                        cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, las = 1, lwd = 2, xaxt = 'n'))
axis(side = 1, at = c( -40, -20, 0, 20, 40),
     labels = c('1915', '1935', '1955', '1975', '1995'), cex.axis = 1.5)
with(plot_Tong_m01, lines(x=ggObj$data$x, y=exp(ggObj$data$y-5*ggObj$data$se), 
                         col='red', lwd = 2))
with(plot_Tong_m01, lines(x=ggObj$data$x, y=exp(ggObj$data$y+5*ggObj$data$se), 
                         col='red', lwd = 2))



gamm2_Tong <- gamm4(log(Resistance) ~ s(year) +
                      log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
                      log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                      year : TMP_Annual_fixed + year : log_PRE_Annual_fixed + year : continentality_fixed + slope_TMPxYear + slope_PRExYear, 
                    random = ~ (1|Species_code/site), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)
viz <- getViz(gamm2_Tong$gam)
plot_Tong_m2 <- plot( sm(viz, 1))


gamm02_Tong <- gamm4(log(Resilience) ~ s(year) +
                      log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + 
                       log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                      year : TMP_Annual_fixed + year : slope_TMPxYear + year : continentality_fixed + slope_TMPxYear + slope_PRExYear, 
                    random = ~ (1|Species_code/site), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)


anova(gamm02_Tong$gam)
anova(gamm02_Tong$mer)
summary(gamm02_Tong$gam)
summary(gamm02_Tong$mer)
model_performance(gamm02_Tong$mer) # R2c and R2m both = 0.57; ICC = 0; RMSE = 0.17
check_singularity(gamm02_Tong$mer)

# plot(gamm02_Tong$gam)

viz <- getViz(gamm02_Tong$gam)
# plot_Tong_m2 <- plot( viz, allTerms = T)
plot_Tong_m02 <- plot( sm(viz, 1))
with(plot_Tong_m02, plot(ggObj$data$x, exp(ggObj$data$y), ylim=c(0.7, 1.2), type='l', 
                      main="Lloret SPEI model", xlab= "Year", ylab= 'Resilience',
                      cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, las = 1, lwd = 2, xaxt = 'n'))
axis(side = 1, at = c( -40, -20, 0, 20, 40),
     labels = c('1915', '1935', '1955', '1975', '1995'), cex.axis = 1.5)
with(plot_Tong_m02, lines(x=ggObj$data$x, y=exp(ggObj$data$y-5*ggObj$data$se), 
                       col='red', lwd = 2))
with(plot_Tong_m02, lines(x=ggObj$data$x, y=exp(ggObj$data$y+5*ggObj$data$se), 
                       col='red', lwd = 2))


#---

gamm0_Li <- gamm4(log(Li_Rt) ~ s(year), 
                  random = ~ (1|Species_code), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)
viz <- getViz(gamm0_Li$gam)
plot_Li_m0 <- plot( sm(viz, 1))


gamm00_Li <- gamm4(log(Li_Rs) ~ s(year), 
                    random = ~ (1|Species_code), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)

# gives random var of 0.000000 for for slope of year of site within species
# hence only logical model is random = ~ (1|Species_code): AIC=-415.5766

anova(gamm00_Li$gam) # p=0.822
anova(gamm00_Li$mer) # unable to fit a sensible random effect on s(year); gives variance =0.000000
summary(gamm00_Li$gam) # R2 adj= - 4.99e-05
summary(gamm00_Li$mer)
model_performance(gamm00_Li$mer) #  R2c = not av and R2m = 0; ICC = na; RMSE = 1.50
check_singularity(gamm00_Li$mer) # even without nested site gives true

# plot(gamm00_Li$gam)

viz <- getViz(gamm00_Li$gam)
# plot_Li_m00 <- plot( viz, allTerms = T)
plot_Li_m00 <- plot( sm(viz, 1))
with(plot_Li_m00, plot(ggObj$data$x, exp(ggObj$data$y), ylim=c(0.8, 1.4), type='l', 
                      main="Isbell null model", xlab= "Year", ylab= 'Resilience',
                      cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, las = 1, lwd = 2, xaxt = 'n'))
axis(side = 1, at = c( -40, -20, 0, 20, 40),
     labels = c('1915', '1935', '1955', '1975', '1995'), cex.axis = 1.5)
with(plot_Li_m00, lines(x=ggObj$data$x, y=exp(ggObj$data$y-5*ggObj$data$se), 
                       col='red', lwd = 2))
with(plot_Li_m00, lines(x=ggObj$data$x, y=exp(ggObj$data$y+5*ggObj$data$se), 
                       col='red', lwd = 2))



gamm1_Li <- gamm4(log(Li_Rt) ~ s(year) +
                    log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                    year : TMP_Annual_fixed + year : log_PRE_Annual_fixed + year : continentality_fixed + slope_TMPxYear + slope_PRExYear, 
                  random = ~ (1|Species_code), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)
viz <- getViz(gamm1_Li$gam)
plot_Li_m1 <- plot( sm(viz, 1))


gamm01_Li <- gamm4(log(Li_Rs) ~ s(year) +
                      log_RWI4 + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                      year : TMP_Annual_fixed + year : log_PRE_Annual_fixed + year : continentality_fixed + slope_TMPxYear + slope_PRExYear, 
                    random = ~ (1|Species_code), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)

# gives random var =0 for for intercept of year of site within species
# hence only logical model is random = ~ (1|Species_code): AIC=-415.5766



anova(gamm01_Li$gam) # p value for smooth terms p=0.951; only CL significant
anova(gamm01_Li$mer)
summary(gamm01_Li$gam) # R2=0.00165
summary(gamm01_Li$mer)
model_performance(gamm01_Li$mer)  # R2c = na; R2m = 2.00e-03; ICC = 0; RMSE = 1.50
check_singularity(gamm01_Li$mer) # still gives true

# plot(gamm01_Li$gam)

viz <- getViz(gamm01_Li$gam)
# plot_Li_m01 <- plot( viz, allTerms = T)
plot_Li_m01 <- plot( sm(viz, 1))
with(plot_Li_m1, plot(ggObj$data$x, exp(ggObj$data$y), ylim=c(0.7, 1.5), type='l', 
                      main="Isbell cilmate model", xlab= "Year", ylab= 'Resilience',
                      cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, las = 1, lwd = 2, xaxt = 'n'))
axis(side = 1, at = c( -40, -20, 0, 20, 40),
     labels = c('1915', '1935', '1955', '1975', '1995'), cex.axis = 1.5)
with(plot_Li_m1, lines(x=ggObj$data$x, y=exp(ggObj$data$y-5*ggObj$data$se), 
                       col='red', lwd = 2))
with(plot_Li_m1, lines(x=ggObj$data$x, y=exp(ggObj$data$y+5*ggObj$data$se), 
                       col='red', lwd = 2))


gamm2_Li <- gamm4(log(Li_Rs) ~ s(year) +
                    log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                    year : TMP_Annual_fixed + year : log_PRE_Annual_fixed + year : continentality_fixed + slope_TMPxYear + slope_PRExYear, 
                  random = ~ (1|Species_code), na.action = na.omit, data = db1.28_Gymn, REML = TRUE)
viz <- getViz(gamm2_Li$gam)
plot_Li_m2 <- plot( sm(viz, 1))


gamm02_Li <- gamm4(log(Li_Rs) ~ s(year) +
                    log_RWI4 + SPEId + SPEIprev + SPEIpost + SPEI3m_Min_Season + log_Max_age + TMP_Annual_fixed + log_PRE_Annual_fixed + continentality_fixed +
                    year : TMP_Annual_fixed + year : log_PRE_Annual_fixed + year : continentality_fixed + slope_TMPxYear + slope_PRExYear, 
                  random = ~ (1|Species_code), na.action = na.omit, data = db_lrg_grp, REML = TRUE)

# gives random var = 0 for for intercept of year of site within species
# hence only logical model is random = ~ (1|Species_code): AIC=39653.65



anova(gamm02_Li$gam) # p value for smooth terms p=0.969; only SPEId significant
anova(gamm02_Li$mer)
summary(gamm02_Li$gam) # R2= 0.00129
summary(gamm02_Li$mer)
model_performance(gamm02_Li$mer)  # R2c= na; R2m = 2.00e-03; ICC = 0; RMSE = 1.50
check_singularity(gamm02_Li$mer) # still gives TRUE!!

# plot(gamm02_Li$gam)

viz <- getViz(gamm02_Li$gam)
# plot_Li_m02 <- plot( viz, allTerms = T)
plot_Li_m02 <- plot( sm(viz, 1))
with(plot_Li_m02, plot(ggObj$data$x, exp(ggObj$data$y), ylim=c(0.7, 1.5), type='l', 
     main="Isbell SPEI model", xlab= "Year", ylab= 'Resilience',
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, las = 1, lwd = 2, xaxt = 'n'))
axis(side = 1, at = c( -40, -20, 0, 20, 40),
     labels = c('1915', '1935', '1955', '1975', '1995'), cex.axis = 1.5)
with(plot_Li_m02, lines(x=ggObj$data$x, y=exp(ggObj$data$y-5*ggObj$data$se), 
                       col='red', lwd = 2))
with(plot_Li_m02, lines(x=ggObj$data$x, y=exp(ggObj$data$y+5*ggObj$data$se), 
                       col='red', lwd = 2))



# trade-offs Rt versus Rs in Tong's models
# amazing; no trade-off found with Lloret's indices
par(mfrow=c(3,2))
with(plot_Li_m0, plot(exp(ggObj$data$y), exp(plot_Li_m00$ggObj$data$y), ylim=c(0.7, 1.5),  
                      main="Isbell null model", xlab= "Resistance", ylab= 'Resilience'))

with(plot_Tong_m0, plot(exp(ggObj$data$y), exp(plot_Tong_m00$ggObj$data$y), ylim=c(0.7, 1.5),  
                        main="Lloret null model", xlab= "Resistance", ylab= 'Resilience'))


with(plot_Li_m1, plot(exp(ggObj$data$y), exp(plot_Li_m01$ggObj$data$y), ylim=c(0.7, 1.5),  
                      main="Isbell climate model", xlab= "Resistance", ylab= 'Resilience'))


with(plot_Tong_m1, plot(exp(ggObj$data$y), exp(plot_Tong_m01$ggObj$data$y), ylim=c(0.7, 1.5),  
                        main="Lloret climate model", xlab= "Resistance", ylab= 'Resilience'))


with(plot_Li_m2, plot(exp(ggObj$data$y), exp(plot_Li_m02$ggObj$data$y), ylim=c(0.7, 1.5),  
                      main="Isbell SPEI model", xlab= "Resistance", ylab= 'Resilience'))


with(plot_Tong_m2, plot(exp(ggObj$data$y), exp(plot_Tong_m02$ggObj$data$y), ylim=c(0.7, 1.5),  
                      main="Lloret SPEI model", xlab= "Resistance", ylab= 'Resilience'))



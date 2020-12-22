
#### Libraries and functions ______________________________________________ ####
library(bestNormalize) # automatically selects best transformation
library(lmerTest) # this gives lsmeans; otherwise use emmeans package with more options
library(tidyverse)
library(data.table) # loads data.table and set functions from library
library(performance)
library(SPEI)
library(plyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(dplR)
library(MuMIn)
library(effects)
library(Taxonstand)
library(taxonlookup)
library(car)
library(scales)
library(Rmisc)
library(emmeans)


meanNA <- function(x) mean(x, na.rm=TRUE)
sterr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
medianNA <- function(x) median(x, na.rm=TRUE)


# is.NaN method for data frames
is.nan.data.frame <- function(x)
{
  do.call(cbind, lapply(x, is.nan))
} 

uniqueNA <- function(x) if(length((x))>1) {
  as.character(unique(x[!is.na(x)]))} else {unique(x)}









#### Drought definition ____________________________________________________####

# subset INFOSET_T to extract non-drought years
db1.28_nondrgt <- subset(INFOSET_T, SPEId + 1.28 > 0) # extracts non-drought years according to Li et al

# mean growth during non-drought years by site for use by Li et al
MeanYn <- with(db1.28_nondrgt, aggregate(TRI, by=list(site), FUN=meanNA))[,2]

# mean Ym by site according to Li et al
INFOSET_T$MeanYn_Li <- rep(MeanYn, 1, each=110)

# subset INFOSET_T to extract drought years
db1.28_Li <- subset(INFOSET_T, SPEId + 1.28 <= 0)

# formulas as in Li et al (2020)
db1.28_Li$Li_Rt <- with(db1.28_Li, MeanYn_Li/abs(RWI_Li-MeanYn_Li))
db1.28_Li$Li_Rs <- with(db1.28_Li, abs((RWI_Li-MeanYn_Li)/(Gpost_Li-MeanYn_Li)))


rm(MeanYn)


#### Pre_post dataset                                _______________________####

# # eliminate +Inf in both Rt and Rs (uses library(data.table))
invisible(lapply(names(db1.28_Li),function(.name) set(db1.28_Li, which(is.infinite(db1.28_Li[[.name]])), 
                                                      j = .name,value =NA)))
# eliminate NAs in both Rt and Rs
db1.28_Li  <- db1.28_Li[complete.cases(db1.28_Li[ , 'Li_Rs']),]


# subsets of time periods for both db1.28 and db1.28_2 df
db1.28_0229 <- subset(db1.28_Li, year_Li <= 1929)
db1.28_3049 <- subset(db1.28_Li, year_Li  <= 1949  & year_Li > 1929)
db1.28_5069 <- subset(db1.28_Li, year_Li  <= 1969  & year_Li > 1949)
db1.28_7089 <- subset(db1.28_Li, year_Li <= 1989 &  year_Li  > 1969)
db1.28_9009 <- subset(db1.28_Li,             year_Li  > 1989)


db1.28_0229$pre_post   <- rep('02_29', length(db1.28_0229$year))
db1.28_3049$pre_post <-   rep('30_49', length(db1.28_3049$year))
db1.28_5069$pre_post   <- rep('50_69', length(db1.28_5069$year))
db1.28_7089$pre_post   <- rep('70_89', length(db1.28_7089$year))
db1.28_9009$pre_post   <- rep('90_09', length(db1.28_9009$year))


db1.28_pre_post <- rbind(db1.28_0229, db1.28_3049, db1.28_5069, db1.28_7089, db1.28_9009)
db1.28_pre_post$pre_post <- as.factor(db1.28_pre_post$pre_post)
db1.28_pre_post$pre_post <- factor(db1.28_pre_post$pre_post, levels = c('02_29', '30_49', '50_69', 
                                                                      '70_89', '90_09'))


#### Subset to only CONIFERS ___________________________________________    ####

db1.28_Gymn <- subset(db1.28_pre_post, Groups == 'Gymnosperms')

# scale all numeric predictors before analysis
db1.28_Gymn$year                 <- as.vector(scale(db1.28_Gymn$year))
db1.28_Gymn$log_RWI4             <- as.vector(scale(db1.28_Gymn$log_RWI4))
db1.28_Gymn$log_RWI42            <- as.vector(scale(db1.28_Gymn$log_RWI4))
db1.28_Gymn$log_Max_age          <- as.vector(scale(db1.28_Gymn$log_Max_age))
db1.28_Gymn$TMP_Annual_fixed     <- as.vector(scale(db1.28_Gymn$TMP_Annual_fixed))
db1.28_Gymn$log_PRE_Annual_fixed <- as.vector(scale(db1.28_Gymn$log_PRE_Annual_fixed))
db1.28_Gymn$continentality_fixed <- as.vector(scale(db1.28_Gymn$continentality_fixed))
db1.28_Gymn$slope_TMPxYear       <- as.vector(scale(db1.28_Gymn$slope_TMPxYear))
db1.28_Gymn$slope_PRExYear       <- as.vector(scale(db1.28_Gymn$slope_PRExYear))

# normalise Lloret resilience and year
Rs_norm <- bestNormalize(db1.28_Gymn$Resilience, allow_lambert_s = TRUE, allow_lambert_h = TRUE) # eliminate original log scale
Rs_norm

db1.28_Gymn$Rs_norm <- Rs_norm $chosen_transform$x.t

#____________________________________________________________________________________________________________
Rs_norm <- bestNormalize(db1.28_Gymn$year, allow_lambert_s = TRUE, allow_lambert_h = TRUE) # eliminate original log scale
Rs_norm
db1.28_Gymn$year_norm <- Rs_norm $chosen_transform$x.t

rm(Rs_norm)


#### Subset to CONIFERS 5-droughts _________________________________________####
# starts from db1.28_Gym_lowRt database and
# subset of sites with at least one drought for each of the periods:1950-1969, 1970-1989 and 1990-2009

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

lrg_grp  <- site_lst[site_lst$select_1=='OK',]  # 368 sites
smll_grp <- site_lst[site_lst$select_1=='OK' & site_lst$select_2=='OK',] # 319 sites

db_lrg_grp <-  db1.28_Gymn[db1.28_Gymn$site %in% lrg_grp$site_lst,]
db_smll_grp <- db1.28_Gymn[db1.28_Gymn$site %in% smll_grp$site_lst,]

db_lrg_grp$site  <-  as.factor(db_lrg_grp$site)
db_smll_grp$site <-  as.factor(db_smll_grp$site)


db_past_1950 <- with(db_lrg_grp, db_lrg_grp[which(pre_post== '50_69' | pre_post== '70_89' | pre_post== '70_89' | pre_post== '90_09'),])
db_past_1950 <- droplevels(db_past_1950)

rm(pippo)


#### Import Li et al data from 332 sites ___________________________________####

# data are given as log(Rt) and log(Rs) already in csv

Li_etal_Rssites <-data.frame(read.table("Li_etal_332_Rsvalues.csv", header=T, sep=","))
Li_etal_Rtsites <-data.frame(read.table("Li_etal_332_Rtvalues.csv", header=T, sep=","))

Li_etal_sites <- merge(Li_etal_Rssites, Li_etal_Rtsites, 'sitename')
colnames(Li_etal_sites) <- c('site', 'Rs_50_69', 'Rs_70_89', 'Rs_90_09', 'Rt_50_69', 'Rt_70_89', 'Rt_90_09')

# uses library data.table to go from wide to long
Li_etal_sites_long <- melt(setDT(Li_etal_sites[,1:7]), id.vars = c("site"), variable.name = "Periods")

Li_etal_sites_long <- Li_etal_sites_long %>%
  mutate(Li_true = case_when(Periods == 'Rs_50_69' ~ 'Rs',
                             Periods == 'Rs_70_89' ~ 'Rs',
                             Periods == 'Rs_90_09' ~ 'Rs',
                             Periods == 'Rt_50_69' ~ 'Rt',
                             Periods == 'Rt_70_89' ~ 'Rt',
                             Periods == 'Rt_90_09' ~ 'Rt'))

Li_etal_sites_long2 <- Li_etal_sites_long[Li_etal_sites_long$Li_true=='Rs']
colnames(Li_etal_sites_long2)[3] <- 'Rs'
Li_etal_sites_long2$Rt <- Li_etal_sites_long$value[Li_etal_sites_long$Li_true=='Rt']
Li_etal_sites_long2$Li_true <- NULL
Li_etal_sites_long2 <- droplevels(Li_etal_sites_long2)


#____________________________________________________________________
# subsetting db1.28_Gymn to list in Li_etal_sites

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

Li_etal_data$site <- as.factor(Li_etal_data$site)
Li_etal_data <- droplevels(Li_etal_data)

Li_etal_data$pre_post <- as.character(Li_etal_data$pre_post)
Li_etal_data$pre_post[which(Li_etal_data$pre_post=='90_08')] <- '90_09'
Li_etal_data$pre_post <- factor(Li_etal_data$pre_post, levels = c('50_69', '70_89', '90_09'))

#____________________________________________________________________
# subsetting Li_etal_sites_long to after 1950; then 
# aggregating to one mean Rs value/site after 1950. Finally, 
# Merging data from Shilong to our data.


Lietal_past_1950 <- with(Li_etal_data, Li_etal_data[which(pre_post== '50_69' | pre_post== '70_89' | pre_post== '90_09'),])
Lietal_past_1950 <- droplevels(Lietal_past_1950)

# Switching to dplyr _______________________________________________________________
df1 <- tibble(Lietal_past_1950)

# make sure data frame is complete with missing combinations of data.
df2 <- complete(df1, site, pre_post)

# uses dplyr
df3 <- df2 %>% group_by(site, pre_post) %>% summarise_if(is.numeric, meanNA) # gives 298 * 3 periods = 894 obs.
df4 <- df2 %>% group_by(site) %>% summarise_if(is.numeric, meanNA) # gives 298 sites (313 if using SPEI<1.28)

# re-add Species_code to df3
spp_cds <- with(df1, aggregate(Species_code, by=list(site), FUN=unique))

df3$Species_code <- rep(spp_cds[,'x'], 1, each = 3)
df3$site2 <- rep(spp_cds[,'Group.1'], 1, each = 3)

df3[is.nan.data.frame(df3)] <- NA  # to eliminate all the NaNs in the df

#____________________________________________________________________
# finally merging the two (left_join)
Li_etal_sites_long2$site <- paste(Li_etal_sites_long2$site, '.crn', sep='')

# subsetting Li_etal_sites_long to our sites
Li_etal_sites_long2 <- Li_etal_sites_long2 %>% filter(site %in% df3$site)
colnames(Li_etal_sites_long2)[2] <- 'pre_post'
colnames(Li_etal_sites_long2)[3] <- 'logLi_Rs_true'
colnames(Li_etal_sites_long2)[4] <- 'logLi_Rt_true'

Li_etal_sites_long2$pre_post <- as.character(Li_etal_sites_long2$pre_post)
Li_etal_sites_long2$pre_post[which(Li_etal_sites_long2$pre_post=='Rs_50_69')] <- '50_69'
Li_etal_sites_long2$pre_post[which(Li_etal_sites_long2$pre_post=='Rs_70_89')] <- '70_89'
Li_etal_sites_long2$pre_post[which(Li_etal_sites_long2$pre_post=='Rs_90_09')] <- '90_09'
Li_etal_sites_long2$pre_post <- as.factor(Li_etal_sites_long2$pre_post)
levels(Li_etal_sites_long2$pre_post) <- c('50_69', '70_89', '90_09')

df3 <- df3[, !(names(df3) %in% 'logLi_Rs_true')]

df3 <- left_join(df3, Li_etal_sites_long2)


#_______________________________________________________________________
# Normalise Li et al data

Li_Rs_norm <- bestNormalize(exp(df3$logLi_Rs_true), allow_lambert_s = TRUE, allow_lambert_h = TRUE) # eliminate original log scale
Li_Rs_norm

df3$Li_Rs_norm <- Li_Rs_norm $chosen_transform$x.t

Li_Rt_norm <- bestNormalize(exp(df3$logLi_Rt_true), allow_lambert_s = TRUE, allow_lambert_h = TRUE) # eliminate original log scale
Li_Rt_norm

df3$Li_Rt_norm <- Li_Rt_norm $chosen_transform$x.t

rm(Li_Rs_norm, Li_Rt_norm)

#_______________________________________________________________________


#### Test co-variance of Rs and Rt against random time series ______________####

  
  library(RcppRoll)
  sd_G   <- with(db1.28_Gymn, sd(log(TRI)))
  mean_G <- with(db1.28_Gymn, mean(log(TRI), na.remove=T))
  
  # draw log-normal distribution centred at real values of mean and var for TRI
  randG <- (rlnorm(1000, meanlog = mean_G, sdlog = sd_G))
  
  rndm <- as.data.frame(randG)
  
  rndm <- rndm %>%
    rownames_to_column() %>%
    mutate(year = as.numeric(rowname)) %>%
    column_to_rownames()
  
  # creation of Gprev and Gpost variables
  rndm$Gpost <-  c(tail(roll_meanr(rndm$randG, n= 4), -4), rep(NA,4))
  rndm$Gprev <-  c(NA, head(roll_meanr(rndm$randG, n= 4), -1))
  
  
  # calculation of the two pairs of Rs and Rt indices
  rndm$Li_Rt   <- with(rndm, mean(randG)/abs(randG - mean(randG)))
  rndm$Li_Rs   <- with(rndm, abs((randG-mean(randG))/(Gpost-mean(randG))))

  rndm <- rndm %>%
    drop_na(Li_Rs, Li_Rt) %>%
    mutate(Li_Rt_norm = bestNormalize(Li_Rt, 
                                      allow_lambert_s = TRUE, 
                                      allow_lambert_h = TRUE)$chosen_transform$x.t) %>%
    mutate(Li_Rs_norm = bestNormalize(Li_Rs, 
                                      allow_lambert_s = TRUE, 
                                      allow_lambert_h = TRUE)$chosen_transform$x.t)
  
# plot Isbell against random draw of points

  library(smatr)
  with(rndm, summary(lm(Li_Rt_norm  ~  Li_Rs_norm, rndm)))
  with(rndm, summary(lm(log(Li_Rt)  ~  log(Li_Rs), rndm)))
  m <- with(rndm, sma(Li_Rt_norm  ~  Li_Rs_norm, rndm))
  m2 <- with(rndm, lm(Li_Rt_norm  ~  Li_Rs_norm, rndm))
  
  par(mfrow=c(1,1))
  par(mar = c(4,4,3,2) + 0.1, mgp = c(2, 1, 0), cex=1, cex.lab=1.6,
      cex.main=1.6, las = 1, cex.axis = 1.4)
  
  # plot Li et al
  with(db1.5_T_Gymn, plot(Li_Rs_norm ~ Li_Rt_norm, 
                          main='Random draw test', 
                          ylab='Resilience', xlab='Resistance'))
  text(x = 8, y = 5, labels = 'Li et al (2020)', cex = 2)
  text(x = 8, y = 3.5, labels = 'Random draw', col = 'red', cex = 2)
  points(x = rndm$Li_Rt_norm, y = rndm$Li_Rs_norm, col = 'red')
  text(x = 8, y = 2, labels = 'R2 = 0.50', col = 'red', cex = 1.6)
  
  
  # plot non-normalised Li et al
  with(db1.5_T_Gymn, plot(log(Li_Rs) ~ log(Li_Rt), 
                          main='Autocorrelation test', 
                          ylab='Resilience', xlab='Resistance'))
  text(x = 8, y = 5, labels = 'Li et al (2020)', cex = 2)
  text(x = 8, y = 3.5, labels = 'Random draw', col = 'red', cex = 2)
  points(x = log(rndm$Li_Rt), y = log(rndm$Li_Rs), col = 'red')
  text(x = 8, y = 2, labels = 'R2 = 0.50', col = 'red', cex = 1.6)
  
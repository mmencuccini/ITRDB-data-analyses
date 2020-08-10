## procedure to create an INFOSET_NEW df containing only the data from the Li et al sites
## and a new data frame (all series), containing all new and old sites (3889 series).

## libraries ________________________________####

library(rprojroot)
library(raster)
require(dplR)
library(stringr)
library(Hmisc)
library(ggplot2)
library(datasets)
library(splusTimeDate)
library(ggmap) # google country name to Lat/Long service (limited to max 2500 queries each day)
# West of Greenwich is -ve; E if +ve
# North equator is +ve; S is -ve

simpleCap <- function(x) { 
  s <- tolower(x) 
  s <- str_trim(str_replace_all(s, "[.]", "")) # drops all '.' and trims string
  s <- strsplit(s, " ")[[1]] 
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ") } 





## Reading Zhao et al 2019 rwl folder into crn format ___________________________#####

wd <- normalizePath(wd, winslash = "/", mustWork = NA)
setwd(wd)

list_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
list_dirs <- list_dirs[!list_dirs %in% './crn']

for(j in 1:length(list_dirs)) {
  if(j==1) {setwd(wd) } else {    setwd('..')   }                         #  reset to home directory
  setwd(list_dirs[j])                  # child directory
  list_files <- list.files(path = '.') # list of files
  # for each directory extracts and process all files and builds chronologies
  for(i in 1:length(list_files)) {
    file <- read.rwl(fname = list_files[i], format = 'tucson')
    if(!all(unlist(lapply(file, function(x) any(is.na(x) | any(x==0)))))) { # is not all NA and not all 0, then
          
      dtrdn_file <- detrend(file, make.plot = FALSE, method = 'ModNegExp', difference = FALSE)
      name <- strsplit(list_files[i], split = '.rwl')
      name <- paste(name, '.crn', sep='')
      df <- chron(x = dtrdn_file, biweight = FALSE, prewhiten = FALSE)
      df$years <- as.numeric(rownames(df))
      
      setwd('..')                          # resets home directory
      setwd('./crn')                     # moves to crn folder
      saveRDS(df, file = name)
      setwd('..')                            #  reset to home directory
      setwd(list_dirs[j])                  # reset to same child directory
      }
    }
}



## Reading new crn files from folder ______________________________####

setwd('..')

Zhao_metadt <- read.csv('rwl_metadata.csv', header = TRUE)
Zhao_metadt <- Zhao_metadt[Zhao_metadt$type == 'Ring Width',]
Zhao_metadt$id <- paste(Zhao_metadt$id, '.crn', sep = '')

setwd('./crn')

list_files <- list.files(path = '.')                                

length(Zhao_metadt$id)                                             
length(list_files)                                                
length(unique(site_lst$site_lst))                                   
length(unique(coord$file_name))                                     
length(unique(Li_etal_sites$site))                                 

length(which((Li_etal_sites$site %in% site_lst$site_lst)))          
Li_etal_sites$site[ !(Li_etal_sites$site %in% site_lst$site_lst)]   
Li_etal_sites$site[ !(Li_etal_sites$site %in% coord$file_name)]     

length(Li_etal_sites$site[Li_etal_sites$site %in% list_files])      
length(Li_etal_sites$site[Li_etal_sites$site %in% Zhao_metadt$id])  


sites_mssng <- setdiff(Li_etal_sites$site , site_lst$site_lst)      


availbl_sites <- intersect(sites_mssng, Zhao_metadt$id)            
unavailbl_sites <- sites_mssng[ !(sites_mssng %in% Zhao_metadt$id)]   
setdiff(sites_mssng, Zhao_metadt$id)                                


## SECOND ROUTINE TO CATCH MISSING FILES___________________####

setwd('..') 
wd <- normalizePath(wd, winslash = "/", mustWork = NA)
setwd(wd)

list_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
list_dirs <- list_dirs[!list_dirs %in% './crn']

availbl_sites_rwl <- paste(unlist(strsplit(availbl_sites, split = '.crn')), '.rwl', sep = '')

for(j in 1:length(list_dirs)) {
  if(j==1) {setwd(wd) } else {    setwd('..')   }             #  reset to home directory
  setwd(list_dirs[j])                                         # child directory
  list_files <- list.files(path = '.') # list of files
  # for each directory extracts and process all files and builds chronologies
  for(site in availbl_sites_rwl) {
    if(site %in% list_files) {
      file <- read.rwl(fname = site, format = 'tucson')
      
      # this option employed the Gazol et al 2016 approach of Neg Exp detrending
      dtrdn_file <- detrend(file, make.plot = FALSE, method = 'ModNegExp', difference = FALSE)
      name <- strsplit(site, split = '.rwl')
        name <- paste(name, '.crn', sep='')
        df <- chron(x = dtrdn_file, biweight = FALSE, prewhiten = FALSE)
        df$years <- as.numeric(rownames(df))
        
        setwd('..')                          # resets home directory
        setwd('./crn')                     # moves to crn folder
        saveRDS(df, file = name)
        setwd('..')                            #  reset to home directory
        setwd(list_dirs[j])                  # reset to same child directory
    }
  }
}
setwd('..')

# now repeat first step to check on number of available files____________####

getwd()
setwd('./crn')

list_files <- list.files(path = '.')                            
length(list_files)                                                
length(Zhao_metadt$id)                                              
length(unique(site_lst$site_lst))                                 
length(unique(coord$file_name))                                    
length(unique(Li_etal_sites$site))                                  

length(which((Li_etal_sites$site %in% site_lst$site_lst)))          
Li_etal_sites$site[ !(Li_etal_sites$site %in% site_lst$site_lst)]   
Li_etal_sites$site[ !(Li_etal_sites$site %in% coord$file_name)]   

length(Li_etal_sites$site[Li_etal_sites$site %in% list_files])      
availbl_sites <- Li_etal_sites$site[Li_etal_sites$site %in% list_files] 

length(Li_etal_sites$site[Li_etal_sites$site %in% Zhao_metadt$id])  

sites_mssng <- setdiff(Li_etal_sites$site , site_lst$site_lst)    


possibl_sites <- intersect(sites_mssng, Zhao_metadt$id)             
impossibl_sites <- sites_mssng[ !(sites_mssng %in% Zhao_metadt$id)]   
setdiff(sites_mssng, Zhao_metadt$id)                             

# crossing available sites from Zhao folder with missing sites in current list 
new_sites <- intersect(sites_mssng, availbl_sites)                  
new_sites


# create INFOSET_NEW folder with RWI data from only the new sites____________####

n_sites <- length(new_sites)
INFOSET_NEW <- as.data.frame(matrix(data = NA, nrow = (2010-1897+1), ncol = n_sites))

for(i in 1:length(new_sites)) {
  name  <- strsplit(new_sites[i], split = '.crn')
  file <- readRDS(new_sites[i])
  file <- with(file, file[which(years > 1896),'xxxstd'])        # first year read is 1897
  if(length(file)>114) {
    INFOSET_NEW[,i] <-  file[1:114]                             # cut to last year =2010
  } else {
    INFOSET_NEW[,i] <-  c(file, rep(NA, (114-length(file))))    # add NA until last year= 2010
  }
  colnames(INFOSET_NEW)[i] <- new_sites[i]

}

# create new complete text file of detrended files____________####

setwd('..')

n_sites <- length(list_files)
all_series <- as.data.frame(matrix(data = NA, nrow = (2010-1897+1), ncol = (1+n_sites)))
colnames(all_series)[1] <- 'year'
all_series$year <- seq(from = 1897, to = 2010, by =1)

for(i in 2:length(list_files)) {
  name  <- strsplit(list_files[i], split = '.crn')
  file <- readRDS(list_files[i])
  file <- with(file, file[which(years > 1896),'xxxstd'])        # first year read is 1897
  file[file==0] <- NA
  if(length(file)>114) {
    all_series[,i] <-  file[1:114]                             # cut to last year =2010
  } else {
    all_series[,i] <-  c(file, rep(NA, (114-length(file))))    # add NA until last year= 2010
  }
  colnames(all_series)[i] <- list_files[i]
  
}





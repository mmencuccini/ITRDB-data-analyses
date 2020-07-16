## procedure to extract missing files from Shilong Piao list of 332 sites in Li_etal_sites df
## libraries ________________________________####

library(rprojroot)
library(dplR)
library(raster)
require(dplR)
library(stringr)
library(Hmisc)
library(ggplot2)
library(maps)
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



wd <- setwd("D:/Dropbox/My Documents/science/projects/grants/China Tong Zheng/Zhao ITRDB dataset/AppendixS1/Cleaned datasets")
wd <- normalizePath(wd, winslash = "/", mustWork = NA)

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
      dtrdn_file <- detrend(file, make.plot = FALSE, method = 'Spline', f = 0.5, nyrs = 30, difference = FALSE)
      name <- strsplit(list_files[i], split = '.rwl')
      name <- paste(name, '.crn', sep='')
      df <- chron(x = dtrdn_file, biweight = FALSE, prewhiten = FALSE)
      df$years <- as.numeric(rownames(df))
      
      setwd('..')                          # resets home directory
      setwd('./crn')                     # moves to crn folder
      save(df, file = name)
      setwd('..')                            #  reset to home directory
      setwd(list_dirs[j])                  # reset to same child directory
      }
    }
}



## Reading new crn files from folder ______________________________####

setwd('..')
wd <- setwd("D:/Dropbox/My Documents/science/projects/grants/China Tong Zheng/Zhao ITRDB dataset/AppendixS1/Cleaned datasets")

Zhao_metadt <- read.csv('rwl_metadata.csv', header = TRUE)
Zhao_metadt <- Zhao_metadt[Zhao_metadt$type == 'Ring Width',]
Zhao_metadt$id <- paste(Zhao_metadt$id, '.crn', sep = '')

setwd('./crn')

list_files <- list.files(path = '.')                                # list of files in Zhao folder, *.crn

length(Zhao_metadt$id)                                              # 4488 files in Zhao csv metadata file, , *.crn
length(list_files)                                                  # 3792 sites in Zhao folder, *.crn
length(unique(site_lst$site_lst))                                   # 1946 currently selected files, *.crn
length(unique(coord$file_name))                                     # 2383 total number of downloaded files, *.crn
length(unique(Li_etal_sites$site))                                  # 332 sites employed by Li et al, *.crn

length(which((Li_etal_sites$site %in% site_lst$site_lst)))          # 139 sites         of Li etal found in list of 1946
Li_etal_sites$site[ !(Li_etal_sites$site %in% site_lst$site_lst)]   # 193 sites         of Li etal missing in list of 1946
Li_etal_sites$site[ !(Li_etal_sites$site %in% coord$file_name)]     # 189 sites         of Li etal missing in coord

length(Li_etal_sites$site[Li_etal_sites$site %in% list_files])      # 125 sites         of Li_etal found in Zhao folder
length(Li_etal_sites$site[Li_etal_sites$site %in% Zhao_metadt$id])  # 320 sites         of Li_etal found in Zhao metadat


sites_mssng <- setdiff(Li_etal_sites$site , site_lst$site_lst)      # 193 sites in Lietal, but not in list of 1946


availbl_sites <- intersect(sites_mssng, Zhao_metadt$id)             # 181 sites left are in Zhao metadt: 12 missing
unavailbl_sites <- sites_mssng[ !(sites_mssng %in% Zhao_metadt$id)] # remaining 12 missing sites   
setdiff(sites_mssng, Zhao_metadt$id)                                # same result


## SECOND ROUTINE TO CATCH MISSING FILES___________________####

wd <- setwd("D:/Dropbox/My Documents/science/projects/grants/China Tong Zheng/Zhao ITRDB dataset/AppendixS1/Cleaned datasets")
wd <- normalizePath(wd, winslash = "/", mustWork = NA)

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
      # if(!all(unlist(lapply(file, function(x) any(is.na(x) | any(x==0)))))) { # is not all NA and not all 0, then
        dtrdn_file <- detrend(file, make.plot = FALSE, method = 'Spline', f = 0.5, nyrs = 30, difference = FALSE)
        name <- strsplit(site, split = '.rwl')
        name <- paste(name, '.crn', sep='')
        df <- chron(x = dtrdn_file, biweight = FALSE, prewhiten = FALSE)
        df$years <- as.numeric(rownames(df))
        
        setwd('..')                          # resets home directory
        setwd('./crn')                     # moves to crn folder
        save(df, file = name)
        setwd('..')                            #  reset to home directory
        setwd(list_dirs[j])                  # reset to same child directory
        # }
    }
  }
}
setwd('..')

# now repeat first step to check on number of available files____________####

getwd()
setwd('./crn')

list_files <- list.files(path = '.')                                # list of files in Zhao folder, *.crn
length(list_files)                                                  # 3888 (3792) sites in Zhao folder, *.crn
length(Zhao_metadt$id)                                              # 4488 files in Zhao csv metadata file, , *.crn
length(unique(site_lst$site_lst))                                   # 1946 currently selected files, *.crn
length(unique(coord$file_name))                                     # 2383 total number of downloaded files, *.crn
length(unique(Li_etal_sites$site))                                  # 332 sites employed by Li et al, *.crn

length(which((Li_etal_sites$site %in% site_lst$site_lst)))          # 139 sites         of Li etal found in list of 1946
Li_etal_sites$site[ !(Li_etal_sites$site %in% site_lst$site_lst)]   # 193 sites         of Li etal missing in list of 1946
Li_etal_sites$site[ !(Li_etal_sites$site %in% coord$file_name)]     # 189 sites         of Li etal missing in coord

length(Li_etal_sites$site[Li_etal_sites$site %in% list_files])      # 221(125) sites    of Li_etal found in Zhao folder
availbl_sites <- Li_etal_sites$site[Li_etal_sites$site %in% list_files] # 221 sites     of Li_etal found in Zhao metadat

length(Li_etal_sites$site[Li_etal_sites$site %in% Zhao_metadt$id])  # 320 sites         of Li_etal found in Zhao metadat

sites_mssng <- setdiff(Li_etal_sites$site , site_lst$site_lst)      # 193 sites in Lietal, but not in list of 1946


possibl_sites <- intersect(sites_mssng, Zhao_metadt$id)             # 181 sites left are in Zhao metadt: 12 missing
impossibl_sites <- sites_mssng[ !(sites_mssng %in% Zhao_metadt$id)] # remaining 12 missing sites   
setdiff(sites_mssng, Zhao_metadt$id)                                # same result

# crossing available sites from Zhao folder with missing sites in current ilst of 1946
new_sites <- intersect(sites_mssng, availbl_sites)                  # 181 new sites found (181+139=320) cvd!!!
new_sites

# found all the new files


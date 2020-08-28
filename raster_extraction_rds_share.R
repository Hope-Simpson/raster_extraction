# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~ Extracting Pixel Values from Raster Stacks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ And other useful chunks of code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## Start Date: 26.08.2020
## Place: London, UK
## Project: LASER General

# clear workspace
rm(list=ls())
gc()

## 0. Set up a function to install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
# List of packages
packages <- c("sp", "raster","rgdal","proj4","ggplot2","gridExtra",
              "readxl", "magrittr", "exactextractr", "sf", "gdalUtils","hagc", "dplyr")
ipak(packages)

## 1. Set up the pathways to access to the source of data (spatial and others)
user_path <-  "C:/Users/Hope/Documents/GitHub"
path_shapefiles <- paste(user_path, "/raster-extraction/data/shapefiles", sep = "") 
path_vulnerability <- paste(user_path, "/raster-extraction/data/raster/vulnerability", sep = "") 
path_comorbidities <- paste(user_path, "/raster-extraction/data/raster/comorbidities", sep = "") 
path_pop <- paste(user_path, "/raster-extraction/data/raster/population", sep = "") 
path_outputs <- paste(user_path, "/raster-extraction/data/outputs", sep = "") 

# read shapefile of AFRO region
setwd(path_shapefiles)
AFRO <- st_read("AFRO_ADM0_201708.shp")
IU <- st_read("AFRO_IUs_201812.shp")

# create lookup field
IU$lookup <- paste0(IU$IU_ID2, IU$ADMIN0ISO3, sep = "")

# transform to pseudo mercator CRS
AFRO <- st_transform(AFRO, 3857)

## 2. Rasterize Framework for extraction of raster datasets
raster.AFRO <- raster(AFRO)
res(raster.AFRO) <- 5000 # 5000 m resolution
values(raster.AFRO) <- 1 
raster.AFRO <- rasterize(AFRO, raster.AFRO)

## 3. Working with the covariates 
#  Create a stack of total and U5 population estimates
#  nb U5 population estimates have been aggregated from 1km rasters 
setwd(path_pop)
list.files(pattern="*tif$")
pop.raster.files <- list.files(path_pop, 
                           pattern="*tif$", full.names=TRUE) 

pop <- raster(pop.raster.files[1])
u5p <- raster(pop.raster.files[2])

# 3.2 Iterate through the raster files, read, recode NAs, adjust extension and crop to the framework area, then pack them all as a raster Stack object.
pop.Covariates <- stack()
for (j in 1:length(pop.raster.files)){
  i <- raster(pop.raster.files[j])
  i[i < 0] <- NA
  i <- extend(i, raster.AFRO)
  i <- raster::resample(i, raster.AFRO, method = 'bilinear')
  i <- crop(i, raster.AFRO)
  pop.Covariates <- stack(pop.Covariates,i)
}

names(pop.Covariates)

# create population density layer
pop_dens <- pop.Covariates[[1]]

# estimate total pop within IUs
IU$total_pop <- exact_extract(pop_dens, IU, 'sum')

# create U5 population density layer
u5pop <- pop.Covariates[[2]]

# estimate total U5 pop within IUs
IU$U5_pop <- exact_extract(u5pop, IU, 'sum')

# 3.1 Load the vulnerability covariates 
setwd(path_vulnerability)
list.files(pattern="*tif$")
raster.files <- list.files(path_vulnerability, 
                           pattern="*tif$", full.names=TRUE) 

# 3.2 Iterate through the raster files, read, recode NAs, adjust extension and crop to the framework area, then pack them all as a raster Stack object.
vuln_Covariates <- stack()
for (j in 1:length(raster.files)){
  i <- raster(raster.files[j])
  i[i < 0] <- NA
  i <- extend(i, raster.AFRO)
  i <- raster::resample(i, raster.AFRO, method = 'bilinear')
  i <- crop(i, raster.AFRO)
  vuln_Covariates <- stack(vuln_Covariates,i)
}

names(vuln_Covariates)

# CALCULATE PROPORTION WITH IMPROVED WATER by summing piped and other improved categories
imp_water <-  vuln_Covariates[[5]] + vuln_Covariates[[6]]

# CALCULATE PROPORTION WITH IMPROVED SANITATION by summing piped and other improved categories
imp_san <-  vuln_Covariates[[2]] + vuln_Covariates[[3]]

# drop the original water and sanitation rasters from the stack
vuln_Covariates <- vuln_Covariates[[-c(2,3,5,6)]]

names(vuln_Covariates)

# HOUSING

# Calculate prev of unimproved housing per pixel
perc_unimp_hous <- 1 - vuln_Covariates[[1]] 

# estimate population with unimproved housing in each pixel
n_unimp_hous <- perc_unimp_hous * pop_dens

# sum total population with unimproved housing in each IU
IU$n_unimp_hous <- exact_extract(n_unimp_hous, IU, 'sum')

# Calculate mean prev of improved housing for each IU
IU$prev_unimp_hous <- IU$n_unimp_hous / IU$total_pop

# WASH 

# estimate population with improved water in each pixel
n_imp_water <- imp_water * pop_dens

# sum total population with improved water in each IU
IU$n_imp_water <- exact_extract(n_imp_water, IU, 'sum')

# Calculate mean prev of improved water for each IU
IU$prev_imp_water <- IU$n_imp_water / IU$total_pop

# estimate population with imp_san in each pixel
n_imp_san <- imp_san * pop_dens

# sum total population with imp_san in each IU
IU$n_imp_san <- exact_extract(n_imp_san, IU, 'sum')

# Calculate mean prev of imp_san for each IU
IU$prev_imp_san <- IU$n_imp_san / IU$total_pop

## NUTRITION

# estimate U5 population with stunting in each pixel
n_stunted <- vuln_Covariates[[2]] * u5pop

# sum U5 population with stunting in each IU
IU$n_stunted <- exact_extract(n_stunted, IU, 'sum')

# estimate  proportion U5  with stunting in each IU
IU$pr_stunted <- IU$n_stunted / IU$U5_pop

# estimate U5 population with wasting in each pixel
n_wasted <- vuln_Covariates[[3]] * u5pop

# sum U5 population with wasting in each IU
IU$n_wasted <- exact_extract(n_wasted, IU, 'sum')

# estimate proportion U5  with wasting in each IU
IU$pr_wasted <- IU$n_wasted / IU$U5_pop

# load the comorbidity rasters
setwd(path_comorbidities)
list.files(pattern="*tif$")
raster.files <- list.files(path_comorbidities, 
                           pattern="*tif$", full.names=TRUE) 

# 3.2 Iterate through the raster and adjust extension and crop for the framework area, and finally pack them all  as a raster Stack object.
cm.Covariates <- stack()
for (j in 1:length(raster.files)){
  i <- raster(raster.files[j])
  i[i < 0] <- NA
  i <- extend(i, raster.AFRO)
  i <- raster::resample(i, raster.AFRO, method = 'bilinear')
  i <- crop(i, raster.AFRO)
  cm.Covariates <- stack(cm.Covariates,i)
}

names(cm.Covariates)

# falciparum malaria

# estimate n of clinical pf infections in each pixel
# assume pf raster is number of cases / population per pixel  
n_pf_inf <- cm.Covariates[[1]] * pop_dens 

# sum total n of pf infections in each IU in 2017
IU$n_pf_inf <- exact_extract(n_pf_inf, IU, 'sum')

# estimate  incidence of pf infections in each IU in 2017 
IU$inc_pf_inf <- IU$n_pf_inf / IU$total_pop

# HIV

# estimate n with HIV in each pixel 2017
n_w_HIV <- cm.Covariates[[2]] * pop_dens / 100

# sum total n living with HIV per IU in 2017
IU$n_w_HIV <- exact_extract(n_w_HIV, IU, 'sum')

# estimate  prev pf HIV in each IU in 2017 
IU$pr_HIV <- IU$n_w_HIV / IU$total_pop

IU_table <- IU
st_geometry(IU_table) <- NULL

setwd(path_outputs)

write.csv(IU_table, "IU_table_contextual.csv", row.names = F)

## Extract rasters to points
raster.5km <- raster(as(AFRO, "Spatial"))
res(raster.5km) <- 5000 # 1000 m resolution
raster.5km <- rasterize(as(AFRO, "Spatial"), raster.5km)
coord.5km <- as.data.frame(rasterToPoints(raster.5km)[,c(1, 2)])
AFRO.map.5km <- rasterToPoints(raster.5km, spatial = TRUE)

# combine the 2 covariates stacks
all_covariates <- stack(vuln_Covariates, cm.Covariates)

# extract raster values to 5km grid
c <- data.frame(raster::extract(all_covariates, 
                                AFRO.map.5km,
                                method = "bilinear"))

c <- cbind(coord.5km, c)

cont.cov <- coord.5km

# rename fields
c <- c %>% rename(pr_imp_housing2015 = housing_proj,
                  pr_u5_stunting2015 = stunting_CLIP,
                  pr_u5_wasting2015 = wasting_CLIP,
                  inc_pf_clinical_2017 = X2019_pf_2017,
                  pr_HIV_2017 = IHME_AFRICA_HIV)

setwd(path_outputs)

saveRDS(c, file = "contextual_covariates.rds")
write.csv(c, "contextual_covariates_matrix.csv", row.names = F)

## Read the matrix with the predictors
setwd(path_outputs)
covariates <- readRDS("contextual_covariates.rds")
PCS <- raster.5km@crs

## Generate a raster stack object with the predictors
Covs <- stack()
for(i in 3:length(covariates)){
  raster1 <- rasterFromXYZ(covariates[,c(1,2,i)],
                           crs = PCS)
  Covs <- stack(Covs, raster1)
}

# check the rasters
names(Covs)
# plot(Covs)
is.factor(Covs)

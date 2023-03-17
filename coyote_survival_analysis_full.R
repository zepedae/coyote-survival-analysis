######### LOAD PACKAGES ###########

library(lubridate)
library(dplyr)
library(tidyr)
library(raster)
library(sp)
library(survival)



########## PREPARE MORTALITY DATA ############

mortality<- read.csv("mortality.csv")

colnames(mortality)[colnames(mortality) == "Coyote.ID."] <- "id"
colnames(mortality)[colnames(mortality) == "Easting"] <- "easting"
colnames(mortality)[colnames(mortality) == "Northing"] <- "northing"
colnames(mortality)[colnames(mortality) == "Date"] <- "mort_date"
colnames(mortality)[colnames(mortality) == "Mortality.cause"] <- "cause"

mortality$mort_date <- mdy(mortality$mort_date)
mort_data <- subset(mortality, select=c("id", "mort_date", "easting", "northing", 
                                        "cause"))

# remove any animals with more than one mort date
mort_data[duplicated(mort_data$id),] # so remove the "not dead" entry for 173 and 
 # both 453s b/c can't determine which is correct
mort_data <- mort_data[mort_data$cause != "Not dead" & mort_data$id != 453,]


########## PREPARE UTM DATA ############

utm1 <- read.csv("utm1.csv")
utm2 <- read.csv("utm2.csv")
utm <- rbind(utm1, utm2)

colnames(utm)[colnames(utm) == "Coyote.ID.."] <- "id"
colnames(utm)[colnames(utm) == "Date"] <- "date"
utm$date <- mdy(utm$date)
colnames(utm)[colnames(utm) == "Time"] <- "time"
colnames(utm)[colnames(utm) == "Easting"] <- "easting"
colnames(utm)[colnames(utm) == "Northing"] <- "northing"

# create a status column that indicates whether the animal was censored 
 # i.e. not recovered post-mortem (0) or confirmed dead (1)
utm <- left_join(utm[, c("id", "northing", "easting", "date")], 
             mort_data[, c("id", "mort_date")], by = "id", all.x = TRUE) 
utm$status <- ifelse(is.na(utm$mort_date), 0, 1)

# remove any points that were taken after mort_date
utm$false_points <- ifelse(is.na(utm$mort_date) | utm$date <= utm$mort_date, 0, 1)
utm_false <- utm[utm$false_points == 1,]
utm <- utm[utm$false_points == 0,]



######### CALCULATE SURVIVAL TIME ############

load("pup_parent_ids.rdata")
birth_ids <- pup_parent_ids$coyote_id

# subset utm data to ids with known birth year
utm <- utm[utm$id %in% birth_ids,]

# create an "end date" that's either the mort_date or last date of tracking, whichever comes last
utm <- utm %>% 
  group_by(id) %>% 
  mutate(max_date = max(date), 
         end_date = data.table::fifelse(is.na(mort_date), 
                               max(date), mort_date)) %>% ungroup                                                      
                                                                          
# calculate survival time 
utm <- pup_parent_ids %>% 
  dplyr::select(coyote_id, birth_year) %>%
  mutate(birth_year = mdy(paste("05-01-", birth_year))) %>% 
  left_join(., utm, by = c("coyote_id" = "id")) %>% 
  mutate(surv_time = as.numeric(difftime(end_date, birth_year))) %>% ungroup()




######### FILTER OUT INDIVIDUALS WITH TOO FEW LOCATIONS ############

# remove individuals whose last tracking day was more than 1 month from mortality
utm$gap <- as.numeric(difftime(utm$end_date-1, utm$max_date, units = "days"))
utm <- utm[utm$gap < 30,]

# global analysis includes data from an animal's entire life
length(unique(utm$id))
subLocations60 <- utm %>% group_by(id) %>% tally()
over60 <- ifelse(subLocations60$n > 60, subLocations60$id, NA)
over60 <- na.omit(over60)
utmGlobal <- utm[utm$id %in% over60,]

# local analysis includes individuals with 30 locations within the last six 
  # months of death or tracking if they weren't recovered post mortem
utmLocal <- utm[utm$date > utm$end_date - 180,]
length(unique(utmLocal$id))
subLocations30 <- utmLocal %>% group_by(id) %>% tally() %>% ungroup()
over30 <- ifelse(subLocations30$n > 30, subLocations30$id, NA)
over30 <- na.omit(over30)
utmLocal <- utmLocal[utmLocal$id %in% over30,]




########## EXTRACT SOCIAL AND ENVIRONMENTAL CHARACTERISTICS ###########

#load landscape raster with covariates
raster <- raster::stack("finalRaster524.grd")

# convert utmGlobal DF to SPDF 
xy <- utmGlobal[, c("easting", "northing")]
utmGlobalSPDF <- SpatialPointsDataFrame(coords = xy, data = utmGlobal, proj4string = CRS("+init=epsg:32616"))
utmGlobalSPDF <- spTransform(utmGlobalSPDF, raster::crs(raster))

# extract covariates
globalValues <- utmGlobalSPDF %>% raster::extract(raster, .) %>% as.data.frame() %>% cbind(., utmGlobalSPDF@data$id)
colnames(globalValues)[colnames(globalValues) == "utmGlobalSPDF@data$id"] <- "id"

# calculate mean and standardize variables
globalValsST <- globalValues %>% group_by(id) %>%
  summarise_all(list(mean))%>% ungroup() %>% 
  mutate_at(vars(-id), ~(scale(.) %>% as.vector))

# convert utmLocal DF to SPDF 
xy <- utmLocal[, c("easting", "northing")]
utmLocalSPDF <- SpatialPointsDataFrame(coords = xy, data = utmLocal, proj4string = CRS("+init=epsg:32616"))
utmLocalSPDF <- spTransform(utmLocalSPDF, raster::crs(raster))

# extract covariates
localValues <- utmLocalSPDF %>% raster::extract(raster, .) %>% as.data.frame() %>% cbind(., utmLocalSPDF@data$id)
colnames(localValues)[colnames(localValues) == "utmLocalSPDF@data$id"] <- "id"

# calculate mean and standardize variables
localValsST <- localValues %>% group_by(id) %>%
  summarise_all(list(mean))%>% ungroup() %>% 
  mutate_at(vars(-id), ~(scale(.) %>% as.vector))


###########  COMBINE SURVIVAL AND SOCIAL/ENVIRONMENTAL DATA  ##########

surv_data <- unique(utm[, c("id", "status", "surv_time")])
global_df <- merge(globalValsST, surv_data, by = "id")
local_df <- merge(localValsST, surv_data, by = "id")


############  COX PROPORTIONAL HAZARDS MODEL  #############

global_fit <- coxph(Surv(surv_time, status) ~ medInc + popDens + nat + dist + propW, data = global_df)
local_fit <- coxph(Surv(surv_time, status) ~ medInc + popDens + nat + dist, data = local_df)

c("id", "date", "time", "easting", "northing", 
  "mort_date", "status")

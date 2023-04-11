# Survival Responses of Coyotes to Urban Social and Environmental Characteristics

Understanding species responses to urbanization is key for gaining insights into 
the ecology and management of wildlife in these rapidly expanding environments. 
Survival is a critical process linking individual-level responses to broader-scale 
dynamics important for population persistence. Notably, the resources and risks
affecting survival in urban wildlife are in large part shaped by aspects of the 
human social system which influence the distribution of resources like green space
and risks like pollutants. I explored the effects of social 
and environmental characteristics on survival in urban coyotes - a species of 
particular management interest due to their conflicts with domestic animals and,
less frequently, humans.

Using location data from radiocollared coyotes captured as pups and tracked as adults, I estimated
the effect of natural and disturbed habitat availability, median income, and human population density
on their survival time. Each location's characteristics were extracted using a geospatial raster (https://github.com/zepedae/social-environmental-raster.git) and the *raster* package, averaged within individuals, and then used as predictors in a Cox proportiona hazards model.

Two models were constructed. The local model uses data from the coyote's last 6 months
of life and the global model uses the animal's entire set of locations. 



### Files
Data used in this study are the property of the Urban Coyote Research Project
and cannot be shared here. Below are descriptions of each file for those
interested in recreating the analysis with their own data.

- "mortality.csv": contains data on the approximate date of death for coyotes
recovered postmortem
- "parent_pup_ids": approximate birth date and RFID of individuals who were 
captured as pups
- "utm1.csv", "utm2.csv": contains the UTM coordinates, date, and time of locations
collected from VHF collared animals
- "finalRaster524.grd": a file containing geospatial social and environmental data
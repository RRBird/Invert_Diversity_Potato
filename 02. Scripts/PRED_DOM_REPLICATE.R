#Main Script for data analysis 

#LIBRARIES AND FUNCTIONS----
library("corrplot")
library("dplyr")
library('lme4')
library("AICcmodavg")
library('arm')
library("spatial")
library("spdep")
library("spatialreg")
library("nlme")
library('sp')
library('sf')

invisible(lapply(paste("Functions/", dir("Functions"), sep=""), function(x) source(x)))


#DATA----

getwd()

point <- read.csv("01. Data/Point_Data.csv")
head(point);dim(point)

field <- read.csv("01. Data/Survey_Data.csv")
head(field);dim(field)

morpho <- read.csv("01. Data/Morpho_Data.csv")
head(morpho);dim(morpho)

obs <- read.csv("01. Data/Observation_Data.csv")
head(obs);dim(obs)


#DATA EXPLORATION AND PREPERATION----


summary(obs) #removing an unneeded column
obs <- obs[,-which(names(obs) =="Notes")]

length(unique(obs$Obs_ID))

unique(obs$Stage)

length(unique(obs$Morphospecies))
length(unique(morpho$Morphospecies))

length(unique(obs$Site))
length(unique(point$Site)) #this looks right since there where a few sites that I didn't find any insects


summary(morpho) #removing some unneeded columns
morpho <- morpho[,-which(names(morpho) =="NOTES")]
morpho <- morpho[,-which(names(morpho) =="In_Database")]
morpho <- morpho[,-which(names(morpho) =="Duplicates")]


unique(morpho$Trophic)

length(morpho$Trophic[morpho$Trophic == 'Parasite'])
length(morpho$Trophic[morpho$Trophic == 'Predator'])

unique(morpho$Order)

summary(point)

length(unique(point$Survey_Field))
length(field$ID)

hist(point$Plant_Height)
hist(point$Ground_Cover)

summary(field) #removing some unneeded columns
field <- field[,-which(names(field) =="Notes")]
field <- field[,-which(names(field) =="NDVIsum_Field")] #decided this wasn't a good variable since the field size differed
field <- field[,-c(7:16)] #removing climate data which isn't needed for analysis 

#exploring some of the variables

hist(field$Field_Area_m2)
hist(field$Crop_Age_Days)
hist(field$Average_Height)
hist(field$Average_GC)

table(field$X500m_Prop_Build) #high proportion of zeros makes it unsuitable for modelling 
field <- field[,-which(names(field) =="X500m_Prop_Build")]

table(field$X500m_Prop_Nat_Graze) #high proportion of zeros makes it unsuitable for modelling
field <- field[,-which(names(field) =="X500m_Prop_Nat_Graze")]

table(field$X500m_Prop_Water)
table(field$X500m_Prop_Crops)

table(field$X1km_Prop_Build)#Not enough variation in data - unsuitable for modelling
field <- field[,-which(names(field) =="X1km_Prop_Build")]

table(field$X1km_Prop_Nat)#high proportion of zeros makes it unsuitable for modelling 
field <- field[,-which(names(field) =="X1km_Prop_Nat")]

table(field$X1km_Prop_Water)
table(field$X1km_Prop_Crops)
table(field$X1km_Prop_Nat_Graze) #high proportion of zeros makes it unsuitable for modelling 
field <- field[,-which(names(field) =="X1km_Prop_Nat_Graze")]

hist(point$Ground_Cover)
hist(point$Plant_Height)

##Variable Correlation----

cordata2 <- data.frame(ID = point$Survey_Field,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position)

cordata2 <- merge(cordata2,field, by = "ID")

cordata2 <- cordata2 %>% dplyr::select(Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X500m_Prop_Crops,X500m_Prop_Water,X1km_Prop_Crops,X1km_Prop_Water,NDVImean_Field,NDVIsum_500m,NDVImean_500m,NDVIsum_1km,NDVImean_1km)

cordata2$Position[cordata2$Position == 'Outer'] <- 1
cordata2$Position[cordata2$Position == 'Inner'] <- 2
cordata2$Position <- as.numeric(cordata2$Position)

str(cordata2)#checking for any other character variables
cor2 <- cor(cordata2,method = "spearman")
colnames(cor2) <- c("Height", "GC", "Position","Day","Crop Age","Field Area","Farming 500m","Water 500m","Farming 1000m","Water 1000m","Field NDVI (M)","NDVI 500m (S)","NDVI 500m (M)","NDVI 1000m (S)","NDVI 1000m (M)")
rownames(cor2) <- c("Height", "GC", "Position","Day","Crop Age","Field Area","Farming in 500m","Water in 500m","Farming in 1000m","Water in 1000m","Field NDVI (M)","NDVI 500m (S)","NDVI 500 (M)","NDVI 1000 (S)","NDVI 1000 (M)")

dev.new(height=8,width=11,dpi=80,pointsize=14,noRStudioGD = T)
corrplot::corrplot(cor2,method="color",  
                   type="upper",addCoef.col = 'black',number.cex = 0.6)
head(cordata2)

#Crops and water 500m are correlated with a design variable
#crops 1000m has a high correlation with a design variable also

#Uncorrelated variables listed below based on correlation plot

#Design variables: Crop age, day sampled and spatial position
#Point Variables: Ground cover, plant height
#Field variables: Field size, field NDVI (mean)
#Landscape variable: water 1km, NDVI (sum) (500m or 1km as these are correlated with each other)


##Prepping data for modelling----

#prep variables
variables <- data.frame(Site = point$Site,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position,ID = point$Survey_Field)

variables <- merge(variables,field, by = "ID")

head(variables);dim(variables)

variables <- variables %>% dplyr::select(ID,Site,Field,Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X1km_Prop_Crops,X1km_Prop_Water,NDVImean_Field,NDVIsum_500m,NDVIsum_1km)


#prep insect observations
head(obs);dim(obs)
head(morpho);dim(morpho)

invert <- data.frame(ID = obs$Site,Morphospecies = obs$Morphospecies)

invert <- merge(invert,morpho,by = "Morphospecies")
head(invert);dim(invert)


#Removing wasps, lepidoptera and Serpintine Leaf Miner to control for chemicals applied during some of the sampling period

length(invert$Trophic[invert$Trophic == 'Parasite'])
invert <- invert[invert$Trophic != 'Parasite',]

table(invert$Order)
length(invert$Trophic[invert$Order == 'Lepidoptera'])
invert <- invert[invert$Trophic != 'Lepidoptera',]

length(invert$Morphospecies[invert$Morphospecies == 'Serpintine_Leaf_Miner'])
invert <- invert[invert$Morphospecies != 'Serpintine_Leaf_Miner',]


invert <- invert %>% dplyr::select(ID, Trophic)
head(invert);dim(invert)

#remove tropic levels which aren't predator or herbivore
unique(invert$Trophic)

length(invert$Trophic[invert$Trophic == 'Unknown'])
invert <- invert[invert$Trophic != 'Unknown',]

length(invert$Trophic[invert$Trophic == 'Ominvore'])
invert <- invert[invert$Trophic != 'Ominvore',]

length(invert$Trophic[invert$Trophic == 'Pollinator'])
invert <- invert[invert$Trophic != 'Pollinator',]

length(invert$Trophic[invert$Trophic == 'Fungivore'])
invert <- invert[invert$Trophic != 'Fungivore',]

length(invert$Trophic[invert$Trophic == 'Hematophagous'])
invert <- invert[invert$Trophic != 'Hematophagous',]

head(invert);dim(invert)


#Counting herbivores and enemies for each site

herb_count <- invert %>%
  filter(Trophic == "Herbivore") %>%
  group_by(ID) %>%
  summarize(count = n()) 
head(herb_count);dim(herb_count)

herb_count <- data.frame(herb_count)
names(herb_count)[names(herb_count) == "ID"] <- "Site"
str(herb_count)

pred_count <- invert %>%
  filter(Trophic == "Predator") %>%
  group_by(ID) %>%
  summarize(count = n()) 
head(pred_count);dim(pred_count)
pred_count <- data.frame(pred_count)
names(pred_count)[names(pred_count) == "ID"] <- "Site"
str(pred_count)

#Add invert count to variables database
mod_data <- merge(variables,herb_count,by = 'Site',all.x = T)
head(mod_data);dim(mod_data)
mod_data$count[is.na(mod_data$count)] <- 0
names(mod_data)[names(mod_data) == "count"] <- "Herb_Count"
str(mod_data)

mod_data <- merge(mod_data,pred_count,by = 'Site',all.x = T)
head(mod_data);dim(mod_data)
mod_data$count[is.na(mod_data$count)] <- 0
names(mod_data)[names(mod_data) == "count"] <- "Pred_Count"
str(mod_data)

#Proportion of Preds

mod_data$Prop_Pred <- ifelse(mod_data$Pred_Count == 0, 0,mod_data$Pred_Count/(mod_data$Herb_Count+mod_data$Pred_Count))

head(mod_data);dim(mod_data)
summary(mod_data$Prop_Pred)
range(mod_data$Prop_Pred)

head(mod_data);dim(mod_data)

#SPATIAL AUTOCORRELATION CHECK----

#All data on specific locations has been removed for privacy reasons however the code is left below for how we undertook the analysis of spatial autocorrelation for point data 


head(point) 
coords <- data.frame(lat = mod_data$Lat, long = mod_data$Long)

head(coords);dim(coords)

##Convert to UTM----
coords_sf <- st_as_sf(coords, coords = c("long", "lat"), crs = 4326)
head(coords_sf);dim(coords_sf)
coords_utm <- st_transform(coords_sf, crs = 32755)

utmutmutm_coords <- st_coordinates(coords_utm)#Extract UTM coordinates (Easting and Northing)

head(utm_coords);dim(utm_coords)
str(utm_coords)
utm_coords <- data.frame(utm_coords)

UTM_Data <- data.frame(Site = mod_data$Site,X_UTM = utm_coords$X,Y_UTM = utm_coords$Y)

mod_data <- merge(mod_data,UTM_Data,by = 'Site',all.x = T)
head(mod_data);dim(mod_data)

#need to assess each field separately at different distances for neighbors

##Moran I----

results_50 <- list()

for (i in unique(mod_data$ID)) {
  field_data <- subset(mod_data, ID == i)
  neighbors <- dnearneigh(field_data[,c(20:21)], 0, 50)
  weights <- nb2listw(neighbors, style = "W")
  results_50[[i]]  <- moran.test(field_data[,which(names(field_data)=="Prop_Pred")], weights)
}
#error in calculations because 50m is too small and causes sub-groups


results_100 <- list()

for (i in unique(mod_data$ID)) {
  field_data <- subset(mod_data, ID == i)
  neighbors <- dnearneigh(field_data[,c(20:21)], 0, 100)
  weights <- nb2listw(neighbors, style = "W")
  results_100[[i]]  <- moran.test(field_data[,which(names(field_data)=="Prop_Pred")], weights)
}

results_100 #all fields have no spatial autocorrelation except S8_F13

results_150 <- list()

for (i in unique(mod_data$ID)) {
  field_data <- subset(mod_data, ID == i)
  neighbors <- dnearneigh(field_data[,c(20:21)], 0, 150)
  weights <- nb2listw(neighbors, style = "W")
  results_150[[i]]  <- moran.test(field_data[,which(names(field_data)=="Prop_Pred")], weights)
}

results_150 #all fields have no spatial autocorrelation

results_200 <- list()

for (i in unique(mod_data$ID)) {
  field_data <- subset(mod_data, ID == i)
  neighbors <- dnearneigh(field_data[,c(20:21)], 0, 200)
  weights <- nb2listw(neighbors, style = "W")
  results_200[[i]]  <- moran.test(field_data[,which(names(field_data)=="Prop_Pred")], weights)
}

results_200 #all fields have no spatial autocorrelation except S6_F11 (smallest field) which was unable to calculate p value so don't try any higher

#Since this analysis includes inverts with wings (i.e. higher dispersal capacity) I will choose the highest distance without introducing errors in calculations (avoid spatial averaging due to too many points included as neighbors) 
#therefore use distance to nearest neighbor as 150m

#No spatial autocorrelation in response variable


#checking height variable

results_Height <- list()

for (i in unique(mod_data$ID)) {
  field_data <- subset(mod_data, ID == i)
  neighbors <- dnearneigh(field_data[,c(20:21)], 0, 100)
  weights <- nb2listw(neighbors, style = "W")
  results_Height[[i]]  <- moran.test(field_data[,which(names(field_data)=="Height")], weights)
}

results_Height 
#No spatial autocorrelation in height


#checking GC variable

results_GC <- list()

for (i in unique(mod_data$ID)) {
  field_data <- subset(mod_data, ID == i)
  neighbors <- dnearneigh(field_data[,c(20:21)], 0, 100)
  weights <- nb2listw(neighbors, style = "W")
  results_GC[[i]]  <- moran.test(field_data[,which(names(field_data)=="GC")], weights)
}

results_GC 
#No spatial autocorrelation in GC


#checking field size variable 

Data_Field_Spatial <- mod_data[!duplicated(mod_data$Field), ]
neighbors <- dnearneigh(Data_Field_Spatial[,c(20:21)], 0, 3500)
weights <- nb2listw(neighbors, style = "W")
results_Field_Area  <- moran.test(Data_Field_Spatial[,which(names(Data_Field_Spatial)=="Field_Area_m2")], weights)



neighbors <- dnearneigh(Data_Field_Spatial[,c(20:21)], 0, 3500)
weights <- nb2listw(neighbors, style = "W")
results_Water  <- moran.test(Data_Field_Spatial[,which(names(Data_Field_Spatial)=="X1km_Prop_Water")], weights)



#MODELLING----

##Modelling design variables----

head(mod_data);dim(mod_data)
mod_data$Sample_Day_Scale <- scale(mod_data$Day_Sampled)
mod_data$Crop_Age_Scale <- scale(mod_data$Crop_Age_Days)

mod_data$Position <- as.factor(mod_data$Position)



Null_Dom <- glmer(Prop_Pred~ 1 + (1|Field),data = mod_data,family = binomial)

Position <- glmer(Prop_Pred~ Position + (1|Field),data = mod_data,family = binomial) 
Age <- glmer(Prop_Pred~ Crop_Age_Days + (1|Field),data = mod_data,family = binomial)
Day <- glmer(Prop_Pred~ Day_Sampled + (1|Field),data = mod_data,family = binomial)

Position_Age1 <- glmer(Prop_Pred~ Position + Crop_Age_Days + (1|Field),data = mod_data,family = binomial) 
Position_Age2 <- glmer(Prop_Pred~ Position * Crop_Age_Days + (1|Field),data = mod_data,family = binomial) 

Position_Day1 <- glmer(Prop_Pred~ Position + Day_Sampled + (1|Field),data = mod_data,family = binomial) 
Position_Day2 <- glmer(Prop_Pred~ Position * Day_Sampled + (1|Field),data = mod_data,family = binomial)

Age_Day1 <- glmer(Prop_Pred~ Crop_Age_Days + Day_Sampled + (1|Field),data = mod_data,family = binomial) 
Age_Day2 <- glmer(Prop_Pred~ Crop_Age_Scale * Sample_Day_Scale + (1|Field),data = mod_data,family = binomial) #needed re scaled variables

Position_Day_Age1 <- glmer(Prop_Pred~ Position + Day_Sampled + Crop_Age_Days + (1|Field),data = mod_data,family = binomial) 
Position_Day_Age2 <- glmer(Prop_Pred~ Position * Day_Sampled + Crop_Age_Days + (1|Field),data = mod_data,family = binomial) 
Position_Day_Age3 <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale+ (1|Field),data = mod_data,family = binomial) #needed re scaled variables
Position_Day_Age4 <- glmer(Prop_Pred~ Position + Crop_Age_Scale  *  Sample_Day_Scale + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence
Position_Day_Age5 <- glmer(Prop_Pred~ Position * Crop_Age_Scale *  Sample_Day_Scale + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence

aictab(list(null = Null_Dom,position = Position, age = Age, day = Day, position_Age = Position_Age1, positionxage = Position_Age2,position_day = Position_Day1, positionxday = Position_Day2, position_day_age = Position_Day_Age1, positionxday_age =Position_Day_Age2, positionxage_day = Position_Day_Age3, position_agexday = Position_Day_Age4,positionxagexday = Position_Day_Age5))

###predictions (Figure 2)----

summary(Position_Day_Age3)

Sample_Day_Pred <- seq(min(mod_data$Sample_Day_Scale),max(mod_data$Sample_Day_Scale),length.out=20) 
Crop_Age_Pred <- seq(min(mod_data$Crop_Age_Scale),max(mod_data$Crop_Age_Scale),length.out=20)
Position_Pred <- c("Outer","Inner")
Position_Pred <- as.factor(Position_Pred)

Design_Pred1 <- expand.grid(Sample_Day_Scale = Sample_Day_Pred, Crop_Age_Scale = Crop_Age_Pred, Position = Position_Pred)

Design_Pred2 <- predictSE(mod = Position_Day_Age3,newdata=Design_Pred1,se.fit = T, type = "link")

Design_Pred3<-data.frame(Design_Pred1,fit.link=Design_Pred2$fit,se.link=Design_Pred2$se.fit)
head(Design_Pred3);dim(Design_Pred3)

Design_Pred3$lci.link<-Design_Pred3$fit.link-(1.96*Design_Pred3$se.link)
Design_Pred3$uci.link<-Design_Pred3$fit.link+(1.96*Design_Pred3$se.link)

Design_Pred3$fit<-invlogit(Design_Pred3$fit.link)
Design_Pred3$se<-invlogit(Design_Pred3$se.link)
Design_Pred3$lci<-invlogit(Design_Pred3$lci.link)
Design_Pred3$uci<-invlogit(Design_Pred3$uci.link)

###Figure 2 ----

Crop_Age_Pred

xx <- Design_Pred3$Sample_Day_Scale == Sample_Day_Pred[10]
x <- Design_Pred3$Crop_Age_Scale == Crop_Age_Pred[10] & Design_Pred3$Position == "Inner"

head(Design_Pred3,3); dim(Design_Pred3)
head(mod_data,3); dim(mod_data)


dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = mod_data$Crop_Age_Scale,y = mod_data$Prop_Pred,xlab = "Crop Age (Days)",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,, xaxt = 'n')

polygon(x = c(Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Inner"],rev(Design_Pred3$Crop_Age_Scale[xx& Design_Pred3$Position == "Inner"])), y = c(Design_Pred3$lci[xx& Design_Pred3$Position == "Inner"],rev(Design_Pred3$uci[xx& Design_Pred3$Position == "Inner"])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Inner"],y = Design_Pred3$fit[xx & Design_Pred3$Position == "Inner"],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Outer"],rev(Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Outer"])), y = c(Design_Pred3$lci[xx& Design_Pred3$Position == "Outer"],rev(Design_Pred3$uci[xx& Design_Pred3$Position == "Outer"])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Outer"],y = Design_Pred3$fit[xx & Design_Pred3$Position == "Outer"],lwd = 2,lty = 2, col = 'grey30')
mtext(side=3,line=0,at = -2.5,'a)',cex=1)

axis(side=1, at=seq(from=min(Design_Pred3$Crop_Age_Scale),to=max(Design_Pred3$Crop_Age_Scale),length.out=6),labels=round(seq(from=min(mod_data$Crop_Age_Days),to=max(mod_data$Crop_Age_Days),length.out=6),-1))

legend('topleft',legend = c('Interior', "Edge"), lty = c(1,2), col = 'grey30',pt.cex = 1)


plot(x =mod_data$Sample_Day_Scale,y=mod_data$Prop_Pred,xlab = "Day Surveryed",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1,-0.3,1.25,-0.3, length =0.1)
mtext(side=3,line=0,at = -2.3,'b)',cex=1)

polygon(x = c(Design_Pred3$Sample_Day_Scale[x],rev(Design_Pred3$Sample_Day_Scale[x])), y = c(Design_Pred3$lci[x],rev(Design_Pred3$uci[x])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Design_Pred3$Sample_Day_Scale[x],y = Design_Pred3$fit[x],lwd = 2,col = 'grey30')

axis(side=1, at=seq(from=min(Design_Pred3$Sample_Day_Scale),to=max(Design_Pred3$Sample_Day_Scale),length.out=6),labels=round(seq(from=min(mod_data$Day_Sampled),to=max(mod_data$Day_Sampled),length.out=6),-1))

##Modelling environmental variables----

head(mod_data);dim(mod_data)
mod_data$Field_Area_Scaled <- scale(mod_data$Field_Area_m2)


summary(Position_Day_Age3)
#point-level variables
GC <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + GC + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence
Height <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + Height + (1|Field),data = mod_data,family = binomial)


#site-level variables
Field_Area <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + Field_Area_Scaled + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence
NDVI_Field  <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + NDVImean_Field + (1|Field),data = mod_data,family = binomial)

#landscape-level variables

Prop_Water  <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + X1km_Prop_Water + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence
NDVI_1000  <- glmer(Prop_Pred~ Position * Crop_Age_Days   +  Day_Sampled + NDVIsum_1km + (1|Field),data = mod_data,family = binomial)

aictab(list(null = Null_Dom, positionxage_day = Position_Day_Age3, GC = GC, height = Height, Field_Area = Field_Area, NDVI_Field = NDVI_Field, Proportion_Water = Prop_Water,NDVI_Landscape_1000 = NDVI_1000))


###Predictions (NDVI field) (Figure 3)----

summary(NDVI_Field)

NDVI_Field_Pred <- seq(min(mod_data$NDVImean_Field ),max(mod_data$NDVImean_Field ),length.out=20)


Enviro_Pred1 <- expand.grid(Sample_Day_Scale = Sample_Day_Pred, Crop_Age_Scale = Crop_Age_Pred, Position = Position_Pred, NDVImean_Field = NDVI_Field_Pred)

Enviro_Pred2 <- predictSE(mod = NDVI_Field,newdata=Enviro_Pred1,se.fit = T, type = "link")

Enviro_Pred3<-data.frame(Enviro_Pred1,fit.link=Enviro_Pred2$fit,se.link=Enviro_Pred2$se.fit)
head(Enviro_Pred3);dim(Enviro_Pred3)

Enviro_Pred3$lci.link<-Enviro_Pred3$fit.link-(1.96*Enviro_Pred3$se.link)
Enviro_Pred3$uci.link<-Enviro_Pred3$fit.link+(1.96*Enviro_Pred3$se.link)

Enviro_Pred3$fit<-invlogit(Enviro_Pred3$fit.link)
Enviro_Pred3$se<-invlogit(Enviro_Pred3$se.link)
Enviro_Pred3$lci<-invlogit(Enviro_Pred3$lci.link)
Enviro_Pred3$uci<-invlogit(Enviro_Pred3$uci.link)

###Figure 3 (NDVI field)----

yy <- Enviro_Pred3$Sample_Day_Scale == Sample_Day_Pred[10]
y <- Enviro_Pred3$Crop_Age_Scale == Crop_Age_Pred[10] & Enviro_Pred3$Position == "Inner" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]
yyy<- Enviro_Pred3$Crop_Age_Scale == Crop_Age_Pred[10] & Enviro_Pred3$Position == "Inner" & Enviro_Pred3$Sample_Day_Scale == Sample_Day_Pred[10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = mod_data$Crop_Age_Scale,y = mod_data$Prop_Pred,xlab = "Crop Age (Days)",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,, xaxt = 'n')

polygon(x = c(Enviro_Pred3$Crop_Age_Scale[yy & Enviro_Pred3$Position == "Inner" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],rev(Enviro_Pred3$Crop_Age_Scale[yy& Enviro_Pred3$Position == "Inner" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]])), y = c(Enviro_Pred3$lci[yy& Enviro_Pred3$Position == "Inner"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],rev(Enviro_Pred3$uci[yy& Enviro_Pred3$Position == "Inner"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Enviro_Pred3$Crop_Age_Scale[yy & Enviro_Pred3$Position == "Inner" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],y = Enviro_Pred3$fit[yy & Enviro_Pred3$Position == "Inner"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],lwd = 2,lty = 1, col = 'grey30')
axis(side=1, at=seq(from=min(Enviro_Pred3$Crop_Age_Scale),to=max(Enviro_Pred3$Crop_Age_Scale),length.out=6),labels=round(seq(from=min(mod_data$Crop_Age_Days),to=max(mod_data$Crop_Age_Days),length.out=6),-1),cex.axis = 1)
mtext("100", side =1, line =1, at =0.7, cex = 0.9)
mtext("a)", side =3, line =0.5, at =-2.75, cex = 0.9)

polygon(x = c(Enviro_Pred3$Crop_Age_Scale[yy & Enviro_Pred3$Position == "Outer" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],rev(Enviro_Pred3$Crop_Age_Scale[yy& Enviro_Pred3$Position == "Outer" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]])), y = c(Enviro_Pred3$lci[yy& Enviro_Pred3$Position == "Outer"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],rev(Enviro_Pred3$uci[yy& Enviro_Pred3$Position == "Outer"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Enviro_Pred3$Crop_Age_Scale[yy & Enviro_Pred3$Position == "Outer" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],y = Enviro_Pred3$fit[yy & Enviro_Pred3$Position == "Outer"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],lwd = 2,lty = 2, col = 'grey30')
legend('topleft',legend = c('Interior', "Edge"), lty = c(1,2), col = 'grey30',pt.cex = 1)


plot(x =mod_data$Sample_Day_Scale,y=mod_data$Prop_Pred,xlab = "Day Surveyed",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')

polygon(x = c(Enviro_Pred3$Sample_Day_Scale[y],rev(Enviro_Pred3$Sample_Day_Scale[y])), y = c(Enviro_Pred3$lci[y],rev(Enviro_Pred3$uci[y])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Enviro_Pred3$Sample_Day_Scale[y],y = Enviro_Pred3$fit[y],lwd = 2,col = 'grey30')
axis(side=1, at=seq(from=min(Enviro_Pred3$Sample_Day_Scale),to=max(Enviro_Pred3$Sample_Day_Scale),length.out=6),labels=round(seq(from=min(mod_data$Day_Sampled),to=max(mod_data$Day_Sampled),length.out=6),-1),cex.axis = 1)
mtext("120", side =1, line =1, at =0.9, cex = 0.9)
mtext("b)", side =3, line =0.5, at =-2.6, cex = 0.9)
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.7)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.7)
arrows(-0.70,-0.43,1.10,-0.43, length =0.1)



plot(x =mod_data$NDVImean_Field,y=mod_data$Prop_Pred,xlab = "Mean NDVI of Field",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')

polygon(x = c(Enviro_Pred3$NDVImean_Field[yyy],rev(Enviro_Pred3$NDVImean_Field[yyy])), y = c(Enviro_Pred3$lci[yyy],rev(Enviro_Pred3$uci[yyy])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Enviro_Pred3$NDVImean_Field[yyy],y = Enviro_Pred3$fit[yyy],lwd = 2,col = 'grey30')
axis(side=1, at=seq(from=min(Enviro_Pred3$NDVImean_Field),to=max(Enviro_Pred3$NDVImean_Field),length.out=6),labels=round(seq(from=min(Enviro_Pred3$NDVImean_Field),to=max(Enviro_Pred3$NDVImean_Field),length.out=6),2),cex.axis = 0.95)
mtext("0.22", side =1, line =1, at =0.22, cex = 0.8)
mtext("0.35", side =1, line =1, at =0.35, cex = 0.8)
mtext("0.49", side =1, line =1, at =0.49, cex = 0.8)
mtext("c)", side =3, line =0.5, at =0.04, cex = 0.9)










#END----


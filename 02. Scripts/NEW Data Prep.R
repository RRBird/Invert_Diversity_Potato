options(scipen = 999) #prevents r from automatically displaying large numbers with scientic notation


#This script contains the data prepping, exploration and manipulation before moving into actual analysis


#Libraries----

library("dplyr")
library('vegan')
library("corrplot")
library("tidyr")

#Data----
point <- read.csv("01. Data/Point_Data.csv")
head(point);dim(point)

field <- read.csv("01. Data/Survey_Data.csv")
head(field);dim(field)

morpho <- read.csv("01. Data/Morphospecies_Data.csv")
head(morpho);dim(morpho)

obs <- read.csv("01. Data/Observation_Data.csv")
head(obs);dim(obs)


#Prepping data----

##Removing unneeded columns----

summary(obs)
obs <- obs[,-which(names(obs) =="Notes")]
obs <- obs[,-which(names(obs) =="Wings")]
head(obs)

head(morpho)
summary(morpho)
morpho <- morpho[,-which(names(morpho)=='NOTES')]
morpho <- morpho[,-which(names(morpho)=='Duplicates')]
morpho <- morpho[,-which(names(morpho)=='Wings')]
morpho <- morpho[,-which(names(morpho)=='Intro_Native')]

summary(point)
hist(point$Plant_Height)
hist(point$Ground_Cover)

summary(field)
field <- field[,-which(names(field)=='Notes')]
field <- field[,-c(7:16)]
field <- field[,-which(names(field) =="NDVIsum_Field")] 
field <- field[,-which(names(field) =="Sample_Date")] 
field <- field[,-which(names(field) =="Potato_Varient_1")] 
field <- field[,-which(names(field) =="Potato_Varient_2")] 
 

table(field$X500m_Prop_Build) #high proportion of zeros makes it unsuitable for modelling 
field <- field[,-which(names(field) =="X500m_Prop_Build")]

table(field$X500m_Prop_Nat_Graze) #high proportion of zeros makes it unsuitable for modelling
field <- field[,-which(names(field) =="X500m_Prop_Nat_Graze")]

table(field$X500m_Rip_Prop)
table(field$X500m_Prop_Crops)
table(field$X500m_Prop_Water)


table(field$X1km_Prop_Build)#Not enough variation in data - unsuitable for modelling
field <- field[,-which(names(field) =="X1km_Prop_Build")]

table(field$X1km_Prop_Nat)#high proportion of zeros makes it unsuitable for modelling 
field <- field[,-which(names(field) =="X1km_Prop_Nat")]

table(field$X1km_Rip_Prop)
table(field$X1km_Prop_Crops)
table(field$X1km_Prop_Water)
table(field$X1km_Prop_Nat_Graze) #high proportion of zeros makes it unsuitable for modelling 
field <- field[,-which(names(field) =="X1km_Prop_Nat_Graze")]

##Some data fixes for functional groups----

table(morpho$Size)
morpho$Size[morpho$Size == "Unknown"] <- NA
morpho$Size[morpho$Size == "5-Oct"] <- "5-10"
morpho$Size[morpho$Size == "Oct-15"] <- "10-15"

#Since all levels greater then 10mm only have 1 morphospecies will combine together to  >10
morpho$Size[morpho$Size == ">30"] <- ">10"
morpho$Size[morpho$Size == "20-30"] <- ">10"
morpho$Size[morpho$Size == "15-20"] <- ">10"
morpho$Size[morpho$Size == "10-15"] <- ">10"

morpho$Size <- as.factor(morpho$Size)
levels(morpho$Size)
morpho$Size <- factor(morpho$Size, levels = c("0-2.5","2.5-5","5-10",">10"))

table(morpho$Trophic)
morpho$Trophic[morpho$Trophic == "Unknown"] <- NA
morpho$Trophic[morpho$Trophic == "Ominvore"] <- "Omnivore"


table(morpho$Hunting.Style)
morpho$Hunting.Style[morpho$Hunting.Style == "Unknown"] <- NA


#Correlation for environmental variables----

cordata <- data.frame(ID = point$Survey_Field,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position)

cordata <- merge(cordata,field, by = "ID")

head(cordata);dim(cordata)

cordata <- cordata %>% dplyr::select(Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X500m_Prop_Crops,X500m_Prop_Water,X1km_Prop_Crops,X1km_Prop_Water,NDVImean_Field,NDVIsum_500m,NDVIsum_1km,X500m_Rip_Prop,X1km_Rip_Prop)

cordata$Position[cordata$Position == 'Outer'] <- 1
cordata$Position[cordata$Position == 'Inner'] <- 2
cordata$Position <- as.numeric(cordata$Position)

str(cordata)#checking for any other character variables
cor <- cor(cordata,method = "spearman")
colnames(cor) <- c("Height", "GC", "Position","Day Sampled","Crop Age","Field Area","Cropping 500m","Water 500m","Cropping 1000m","Water 1000m","Field NDVI (M)","NDVI 500m (S)","NDVI 1000m (S)","Riparian 500m","Riparian 1000m")
rownames(cor) <- c("Height", "GC", "Position","Day Sampled","Crop Age","Field Area","Cropping 500m","Water 500m","Cropping 1000m","Water 1000m","Field NDVI (M)","NDVI 500m (S)","NDVI 1000m (S)","Riparian 500m","Riparian 1000m")

dev.new(height=8,width=8,dpi=80,pointsize=14,noRStudioGD = T)
corrplot::corrplot(cor,method="color",  
                   type="upper",addCoef.col = 'black',number.cex = 0.6)
head(cor)


#Day sampled is needed to account for different sampling days so the ones correlated with this variable need to be removed which are: Cropping 500m, Water 500m, cropping 1km, riparian 500m and riparian 1km

#Checking how many correlated variables are left when remove the ones above
head(cordata)
cordata2 <- cordata %>% dplyr::select(-X500m_Prop_Crops,-X500m_Prop_Water,-X1km_Prop_Crops,-X500m_Rip_Prop,-X1km_Rip_Prop)
cor2 <- cor(cordata2,method = "spearman")

dev.new(height=8,width=8,dpi=80,pointsize=14,noRStudioGD = T)
corrplot::corrplot(cor2,method="color",  
                   type="upper",addCoef.col = 'black',number.cex = 0.6)

#only correlated variable left is NDVI sum 500m and 1km (expected) - will proceed with 1km so it's on the same scale as other landscape variable (water 1km)

head(cor)

#Merging data bases----

variables <- data.frame(Site = point$Site,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position,ID = point$Survey_Field)

variables <- merge(variables,field, by = "ID")

head(variables);dim(variables)

variables <- variables %>% dplyr::select(ID,Site,Field,Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X1km_Prop_Water, NDVImean_Field,NDVIsum_1km)
head(variables);dim(variables)

invert <- data.frame(ID = obs$Site,Morphospecies = obs$Morphospecies)
head(invert)

invert <- merge(invert,morpho,by = "Morphospecies")
head(invert);dim(invert)

##removing wasps, Lepidoptera and serpintine leaf miner
dim(invert)

table(invert$Trophic)
length(which(invert$Trophic == 'Parasite'))
invert <- subset(invert, Trophic != 'Parasite' | is.na(Trophic))

table(invert$Order)
length(which(invert$Order == 'Lepidoptera'))
invert <- subset(invert, Order != 'Lepidoptera' | is.na(Order))

length(invert$Morphospecies[invert$Morphospecies == 'Serpintine_Leaf_Miner'])
invert <- subset(invert, Morphospecies != 'Serpintine_Leaf_Miner' | is.na(Morphospecies))

dim(morpho)

table(morpho$Trophic)
length(which(morpho$Trophic == 'Parasite'))
morpho <- subset(morpho, Trophic != 'Parasite' | is.na(Trophic))

table(morpho$Order)
length(which(morpho$Order == 'Lepidoptera'))
morpho <- subset(morpho, Order != 'Lepidoptera' | is.na(Order))

length(morpho$Morphospecies[morpho$Morphospecies == 'Serpintine_Leaf_Miner'])
morpho <- morpho[morpho$Morphospecies != 'Serpintine_Leaf_Miner',]

##Removing first survey observations
#diversity measures taken weren't reliable so excluded from analysis

head(variables);dim(variables)

length(variables$ID[variables$ID == 'S1_F3'])
variables <- variables[variables$ID != 'S1_F3',]

head(invert);dim(invert)

invert <- invert[!grepl("^S1", invert$ID), ]

#total observations = 2329

setdiff(unique(morpho$Morphospecies), unique(invert$Morphospecies))

#need to remove two from morpho to make sure that its got the same species are were observed
#double checking they aren't there
invert[invert$Morphospecies == 'Whitefly',]
invert[invert$Morphospecies == 'Brown_Weevil',] 
invert[invert$Morphospecies == 'Small_Ant',]

morpho <- morpho[morpho$Morphospecies != 'Brown_Weevil',]
morpho <- morpho[morpho$Morphospecies != 'Small_Ant',]
head(morpho);dim(morpho)

#Checking the groups----
head(invert)
length(point$Site) #total of 270 points

table(invert$Order) #must have more then number of points in order to be included in functional diversity analysis 

table(invert$Order,invert$ID)

table(invert$Order,invert$Trophic)
table(invert$Order,invert$Size)
table(invert$Order,invert$Hunting.Style)

#Prep for taxonomic Modelling----

names(invert)[names(invert) == "ID"] <- "Site"

TaxModel <- data.frame(Site = variables$Site)
head(TaxModel);dim(TaxModel)

##Speceis Richness----

richness <- aggregate(Morphospecies ~ Site, data = invert, FUN = function(x) length(unique(x)))
TaxModel <- merge(TaxModel,richness,by = "Site",all.x = T)
head(TaxModel);dim(TaxModel)
colnames(TaxModel)[2] <- "Species_Rich"
head(TaxModel);dim(TaxModel)

min(TaxModel$Species_Rich,na.rm=T)
max(TaxModel$Species_Rich,na.rm=T)

TaxModel$Species_Rich[is.na(TaxModel$Species_Rich)] <- 0

##Diversity (Inverse Simpson's diversity index)----

#diversity(table(x), index = "invsimpson")
#table (x) makes a table of species counts which is used to calculate diversity

diversity <- aggregate(Morphospecies ~ Site, data = invert, FUN = function(x) diversity(table(x), index = "invsimpson"))
TaxModel <- merge(TaxModel,diversity,by = "Site",all.x = T)
head(TaxModel);dim(TaxModel)
colnames(TaxModel)[3] <- "Diversity"
head(TaxModel);dim(TaxModel)

min(TaxModel$Diversity,na.rm=T)
max(TaxModel$Diversity,na.rm=T)

TaxModel$Diversity[is.na(TaxModel$Diversity)] <- 0

##Adding Variables and XY Coordinates----

TaxModel <- merge(TaxModel,variables,by = "Site")
head(TaxModel);dim(TaxModel)

#need XY in data to check for spatial autocorrelation in models later on in process

head(point);dim(point)

coords <- point %>% dplyr::select(Site, X_Cor,Y_Cor,Survey)
head(coords);dim(coords)

coords <- coords[coords$Survey != 1, ]
coords$Survey <- NULL

TaxModel <- merge(TaxModel,coords,by = "Site")
head(TaxModel);dim(TaxModel)


#Setting up for RLQ Analysis----
#TO DO----
#Libraries

library("dplyr")



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

summary(morpho)
morpho <- morpho[,-which(names(morpho)=='NOTES')]
morpho <- morpho[,-which(names(morpho)=='In_Database')]
morpho <- morpho[,-which(names(morpho)=='Duplicates')]

summary(point)
point <- point[,-which(names(point)=='X_Cor')]
point <- point[,-which(names(point)=='Y_Cor')]

summary(field)
field <- field[,-which(names(field)=='Notes')]
field <- field[,-c(7:16)]
field <- field[,-which(names(field) =="NDVIsum_Field")] 
field <- field[,-which(names(field) =="Potato_Varient_1")] 
field <- field[,-which(names(field) =="Potato_Varient_2")] 

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

#CORRELATION TO DO----

##Merging data bases----

variables <- data.frame(Site = point$Site,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position,ID = point$Survey_Field)

variables <- merge(variables,field, by = "ID")

head(variables);dim(variables)

variables <- variables %>% dplyr::select(ID,Site,Field,Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X1km_Prop_Water,NDVImean_Field,NDVIsum_1km)


invert <- data.frame(ID = obs$Site,Morphospecies = obs$Morphospecies)

invert <- merge(invert,morpho,by = "Morphospecies")
head(invert);dim(invert)


##removing wasps, Lepidoptera and serpintine leaf miner----
dim(invert)

length(invert$Trophic[invert$Trophic == 'Parasite'])
invert <- invert[invert$Trophic != 'Parasite',]

table(invert$Order)
length(invert$Trophic[invert$Order == 'Lepidoptera'])
invert <- invert[invert$Order != 'Lepidoptera',]

length(invert$Morphospecies[invert$Morphospecies == 'Serpintine_Leaf_Miner'])
invert <- invert[invert$Morphospecies != 'Serpintine_Leaf_Miner',]


##Removing first survey observations----
#diversity measures taken weren't reliable so excluded from analysis

head(variables)

length(variables$ID[variables$ID == 'S1_F3'])
variables <- variables[variables$ID != 'S1_F3',]

head(invert)

dim(invert)

invert <- invert[!grepl("^S1", invert$ID), ]

#Diversity measures for modelling----

##Richness----

###Araneae
richness <- aggregate(Morphospecies ~ ID, data = invert[invert$Order == "Araneae",], FUN = function(x) length(unique(x)))

##Add the richness number to the database

###Hemieptera

###Coleoptera

  #Diversity -> araneae, hemi, coleoptera
  #Community Comp -> araneae, hemi, coleoptera
  #Beta Diversity -> araneae, hemi, coleoptera
  #Functional abundance
    #araneae - hunting type
    #hemi - size and trophic
    #coleoptera - size and trophic

#Then need to check all these calculations for spatial autocorrelation
#Don't forget to remove all the cordinates etc. after this 

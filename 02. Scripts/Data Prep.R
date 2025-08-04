options(scipen = 999)

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

summary(morpho)
morpho <- morpho[,-which(names(morpho)=='NOTES')]
morpho <- morpho[,-which(names(morpho)=='Duplicates')]

summary(point)
hist(point$Plant_Height)
hist(point$Ground_Cover)

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

table(morpho$Hunting.Style)
morpho$Hunting.Style[morpho$Hunting.Style == "Unknown"] <- NA

table(morpho$Intro_Native) 
morpho$Intro_Native[morpho$Intro_Native == "Unknown"] <- NA
morpho$Intro_Native <- as.factor(morpho$Intro_Native)
levels(morpho$Intro_Native)
morpho$Intro_Native <- factor(morpho$Intro_Native, levels = c("Native","Introduced"))

table(morpho$Wings) 
morpho$Wings[morpho$Wings == "Unknown"] <- NA

#Correlation for environmental variables----

cordata <- data.frame(ID = point$Survey_Field,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position)

cordata <- merge(cordata,field, by = "ID")

head(cordata);dim(cordata)

cordata <- cordata %>% dplyr::select(Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X500m_Prop_Crops,X500m_Prop_Water,X1km_Prop_Crops,X1km_Prop_Water,NDVImean_Field,NDVIsum_500m,NDVImean_500m,NDVIsum_1km,NDVImean_1km)

cordata$Position[cordata$Position == 'Outer'] <- 1
cordata$Position[cordata$Position == 'Inner'] <- 2
cordata$Position <- as.numeric(cordata$Position)

str(cordata)#checking for any other character variables
cor <- cor(cordata,method = "spearman")
colnames(cor) <- c("Height", "GC", "Position","Day Sampled","Crop Age","Field Area","Cropping 500m","Water 500m","Cropping 1000m","Water 1000m","Field NDVI (M)","NDVI 500m (S)","NDVI 500m (M)","NDVI 1000m (S)","NDVI 1000m (M)")
rownames(cor) <- c("Height", "GC", "Position","Day Sampled","Crop Age","Field Area","Cropping 500m","Water 500m","Cropping 1000m","Water 1000m","Field NDVI (M)","NDVI 500m (S)","NDVI 500m (M)","NDVI 1000m (S)","NDVI 1000m (M)")

dev.new(height=8,width=8,dpi=80,pointsize=14,noRStudioGD = T)
corrplot::corrplot(cor,method="color",  
                   type="upper",addCoef.col = 'black',number.cex = 0.6)
head(cor)

#Crops and water 500m are correlated with a design variable
#crops 1000m has a high correlation with a design variable also

#Uncorrelated variables listed below based on correlation plot

#Design variables: Crop age, day sampled and spatial position
#Point Variables: Ground cover, plant height
#Field variables: Field size, field NDVI (mean)
#Landscape variable: water 1km, NDVI (sum) (500m or 1km as these are correlated with each other)


##Merging data bases----

variables <- data.frame(Site = point$Site,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position,ID = point$Survey_Field)

variables <- merge(variables,field, by = "ID")

head(variables);dim(variables)

variables <- variables %>% dplyr::select(ID,Site,Field,Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X1km_Prop_Water,NDVImean_Field,NDVIsum_1km)
head(variables);dim(variables)

invert <- data.frame(ID = obs$Site,Morphospecies = obs$Morphospecies)

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


#Diversity measures for modelling----

table(invert$Order)
dim(variables)


##Richness----

head(invert);dim(invert)
head(variables);dim(variables)

names(invert)[names(invert) == "ID"] <- "Site"

richness <- data.frame(Site = variables$Site)
head(richness);dim(richness)

#all species richness
temprich <- aggregate(Morphospecies ~ Site, data = invert, FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[2] <- "All"
head(richness);dim(richness)

#Trophic - Predator
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Trophic == "Predator",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[3] <- "Predator"

#Trophic - Herbivore
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Trophic == "Herbivore",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[4] <- "Herbivore"

#Trophic - Omnivore
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Trophic == "Ominvore",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[5] <- "Omnivore"

#Trophic - Fungivore
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Trophic == "Fungivore",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[6] <- "Fungivore"

#Trophic - Hematophagous
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Trophic == "Hematophagous",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[7] <- "Hematophagous"

#Hunting Style - Web
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Hunting.Style == "Web",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[8] <- "Web"

#Hunting Style - Active Hunting
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Hunting.Style == "Active_Hunting",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[9] <- "Active_Hunting"

#Hunting Style - Ambush Hunting
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Hunting.Style == "Ambush_Hunter",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[10] <- "Ambush_Hunting"

#Hunting Style - Hawking
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Hunting.Style == "Hawking",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[11] <- "Hawking"

#Size - 0-2.5mm
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Size == "0-2.5",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[12] <- "0-2.5"

#Size - 2.5-5mm
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Size == "2.5-5",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[13] <- "2.5-5"

#Size - 5-10mm
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Size == "5-10",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[14] <- "5-10"

#Size - >10mm
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Size == ">10",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[15] <- ">10"

#Native/intro - Native
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Intro_Native == "Native",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[16] <- "Native"

#Native/intro - Introduced 
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Intro_Native == "Introduced",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[17] <- "Introduced"

#Wings - Always Winged 
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Wings == "Always_Winged",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[18] <- "Always_Winged"

#Wings - Develops Wings 
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Wings == "Develops_Wings",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[19] <- "Develops_Wings"

#Wings - Wingless 
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Wings == "Wingless",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[20] <- "Wingless"

#Wings - Polymorphic 
temprich <- aggregate(Morphospecies ~ Site, data = invert[invert$Wings == "Polymorphic",], FUN = function(x) length(unique(x)))
richness <- merge(richness,temprich,by = "Site",all.x = T)
head(richness);dim(richness)
colnames(richness)[21] <- "Polymorphic"

#replacing all NAs with 0's 
richness[is.na(richness)] <- 0

###Creating modelling data----

head(variables);dim(variables)
head(richness);dim(richness)

ModelRich <- merge(richness,variables,by = "Site")
head(ModelRich);dim(ModelRich)

#workong out prop zero for each functional group
head(richness[,c(2:21)]);dim(richness)
sapply(richness[,c(2:21)], function(col) mean(col == 0, na.rm = TRUE))

table(morpho$Morphospecies)

sum(table(invert$Wings))


#Removing native and introduced - I think this could be very misleading results since only 42 of 215 inverts were able to be classified 
#None of the functional groups were able to be identified for all species but I think this is the one that's the most in disproportion
head(ModelRich);dim(ModelRich)
ModelRich <- ModelRich %>% dplyr::select(-Native,-Introduced)

##Abundance----
head(invert);dim(invert)

abundance <- data.frame(Site = variables$Site)
head(abundance);dim(abundance)

#all species 

tempabund <- aggregate(Morphospecies ~ Site, data = invert, FUN = function(x) length(x))
abundance <- merge(abundance,tempabund,by = "Site",all.x = T)
head(abundance);dim(abundance)
colnames(abundance)[2] <- "All"
head(abundance);dim(abundance)

#Trophic

trophic.list <- list("Predator", "Herbivore","Ominvore","Fungivore","Hematophagous")

p <- 3
for (i in trophic.list) {

  tempabund <- aggregate(Morphospecies ~ Site, 
                         data = invert[invert$Trophic == i,], 
                         FUN = function(x) length(x))
  
  abundance <- merge(abundance,tempabund,by = "Site",all.x = T)
  
  colnames(abundance)[p] <- i
  
  p <- p+1
}

head(abundance);dim(abundance)

colnames(abundance)[5] <- "Omnivore"

#Hunting Style

hunt.list <- list("Web","Active_Hunting","Ambush_Hunter","Hawking")

for (i in hunt.list) {
  
  tempabund <- aggregate(Morphospecies ~ Site, 
                         data = invert[invert$Hunting.Style == i,], 
                         FUN = function(x) length(x))
  
  abundance <- merge(abundance,tempabund,by = "Site",all.x = T)
  
  colnames(abundance)[p] <- i
  
  p <- p+1
}

head(abundance);dim(abundance)

#Size

size.list <- list("0-2.5","2.5-5","5-10",">10")

for (i in size.list) {
  
  tempabund <- aggregate(Morphospecies ~ Site, 
                         data = invert[invert$Size == i,], 
                         FUN = function(x) length(x))
  
  abundance <- merge(abundance,tempabund,by = "Site",all.x = T)
  
  colnames(abundance)[p] <- i
  
  p <- p+1
}

head(abundance);dim(abundance)

#Wings

wing.list <- list("Always_Winged","Develops_Wings","Wingless","Polymorphic")

for (i in wing.list) {
  
  tempabund <- aggregate(Morphospecies ~ Site, 
                         data = invert[invert$Wings == i,], 
                         FUN = function(x) length(x))
  
  abundance <- merge(abundance,tempabund,by = "Site",all.x = T)
  
  colnames(abundance)[p] <- i
  
  p <- p+1
}

head(abundance);dim(abundance)

###Creating modelling data----

abundance[is.na(abundance)] <- 0 #replacing all NAs with 0's 

head(variables);dim(variables)
head(abundance);dim(abundance)

ModelAbund <- merge(abundance,variables,by = "Site")
head(ModelAbund);dim(ModelAbund)

##Diversity (Inverse Simpson's diversity index)----

#diversity(table(x), index = "invsimpson")
#table (x) makes a table of species counts which is used to calulate diversity

head(invert);dim(invert)

diversity <- data.frame(Site = variables$Site)
head(diversity);dim(diversity)

#all species 

tempdiv <- aggregate(Morphospecies ~ Site, data = invert, FUN = function(x) diversity(table(x), index = "invsimpson"))
diversity <- merge(diversity,tempdiv,by = "Site",all.x = T)
head(diversity);dim(diversity)
colnames(diversity)[2] <- "All"
head(diversity);dim(diversity)

min(diversity$All,na.rm=T)
max(diversity$All,na.rm=T)

#Trophic

g <- 3
for (i in trophic.list) {
  
  tempdiv <- aggregate(Morphospecies ~ Site, 
                         data = invert[invert$Trophic == i,], 
                         FUN = function(x) diversity(
                           table(x),index = "invsimpson"))
  
  diversity <- merge(diversity,tempdiv,by = "Site",all.x = T)
  
  colnames(diversity)[g] <- i
  
  g <- g+1
}

head(diversity);dim(diversity)

colnames(diversity)[5] <- "Omnivore"

#Hunting style

for (i in hunt.list) {
  
  tempdiv <- aggregate(Morphospecies ~ Site, 
                       data = invert[invert$Hunting.Style == i,], 
                       FUN = function(x) diversity(
                         table(x),index = "invsimpson"))
  
  diversity <- merge(diversity,tempdiv,by = "Site",all.x = T)
  
  colnames(diversity)[g] <- i
  
  g <- g+1
}

head(diversity);dim(diversity)

#Size

for (i in size.list) {
  
  tempdiv <- aggregate(Morphospecies ~ Site, 
                       data = invert[invert$Size == i,], 
                       FUN = function(x) diversity(
                         table(x),index = "invsimpson"))
  
  diversity <- merge(diversity,tempdiv,by = "Site",all.x = T)
  
  colnames(diversity)[g] <- i
  
  g <- g+1
}

head(diversity);dim(diversity)


#Wings

for (i in wing.list) {
  
  tempdiv <- aggregate(Morphospecies ~ Site, 
                       data = invert[invert$Wings == i,], 
                       FUN = function(x) diversity(
                         table(x),index = "invsimpson"))
  
  diversity <- merge(diversity,tempdiv,by = "Site",all.x = T)
  
  colnames(diversity)[g] <- i
  
  g <- g+1
}

head(diversity);dim(diversity)


###Creating modelling data----


head(variables);dim(variables)
head(diversity);dim(diversity)

ModelDiv <- merge(diversity,variables,by = "Site")
head(ModelDiv);dim(ModelDiv)


#Binomial data----

head(invert);dim(invert)

occurance <- data.frame(Site = variables$Site)
head(occurance);dim(occurance)

#Not all functional groups are getting binomial models


#Trophic

trophic.list2 <- trophic.list[2:5]

#END----
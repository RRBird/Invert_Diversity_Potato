options(scipen = 999)

#Libraries----

library("dplyr")
library('vegan')



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

#Some data fixes for functional groups----

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


##UP TO HERE----
##Still need to clear the below to make sure it still fits into the analysis plans


#CORRELATION TO DO----

#transplant same correlation from pred dom since it's the same variables

##Merging data bases----

variables <- data.frame(Site = point$Site,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position,ID = point$Survey_Field)

variables <- merge(variables,field, by = "ID")

head(variables);dim(variables)

variables <- variables %>% dplyr::select(ID,Site,Field,Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X1km_Prop_Water,NDVImean_Field,NDVIsum_1km)


invert <- data.frame(ID = obs$Site,Morphospecies = obs$Morphospecies)

invert <- merge(invert,morpho,by = "Morphospecies")
head(invert);dim(invert)


##removing wasps, Lepidoptera and serpintine leaf miner
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

head(variables);dim(variables)

length(variables$ID[variables$ID == 'S1_F3'])
variables <- variables[variables$ID != 'S1_F3',]

head(invert)

dim(invert)

invert <- invert[!grepl("^S1", invert$ID), ]

#Diversity measures for modelling----

table(invert$Order)
dim(variables)

#We will only calculate diversity measures for taxon Orders with more observations than the number of sites

#Araneae, Coleoptera, Hemiptera

#remaining orders don't met the above rule even when combined together into one group (other)


##Richness----

head(invert);dim(invert)
head(variables);dim(variables)
head(richness);dim(richness)

###Araneae

names(invert)[names(invert) == "ID"] <- "Site"

richness <- aggregate(Morphospecies ~ Site, data = invert[invert$Order == "Araneae",], FUN = function(x) length(unique(x)))

variables <- merge(variables,richness, by = "Site",all.x = T)
names(variables)[names(variables) == "Morphospecies"] <- "Richness_A"
variables$Richness_A[is.na(variables$Richness_A)] <- 0


(table(variables$Richness_A)/sum(table(variables$Richness_A)))*100

dev.new(height=20,width=40,dpi=80,pointsize=14,noRStudioGD = T)
plot(as.factor(variables$Site), variables$Richness_A, type = "p", xlab = "Site", ylab = "Richness (spiders)",)

###Hemieptera

richness <- aggregate(Morphospecies ~ Site, data = invert[invert$Order == "Hemiptera",], FUN = function(x) length(unique(x)))

variables <- merge(variables,richness, by = "Site",all.x = T)
names(variables)[names(variables) == "Morphospecies"] <- "Richness_H"
variables$Richness_H[is.na(variables$Richness_H)] <- 0


(table(variables$Richness_H)/sum(table(variables$Richness_H)))*100

dev.new(height=20,width=40,dpi=80,pointsize=14,noRStudioGD = T)
plot(as.factor(variables$Site), variables$Richness_H, type = "p", xlab = "Site", ylab = "Richness (True Bugs)",)

###Coleoptera

richness <- aggregate(Morphospecies ~ Site, data = invert[invert$Order == "Coleoptera",], FUN = function(x) length(unique(x)))

variables <- merge(variables,richness, by = "Site",all.x = T)
names(variables)[names(variables) == "Morphospecies"] <- "Richness_C"
variables$Richness_C[is.na(variables$Richness_C)] <- 0



(table(variables$Richness_C)/sum(table(variables$Richness_C)))*100

dev.new(height=20,width=40,dpi=80,pointsize=14,noRStudioGD = T)
plot(as.factor(variables$Site), variables$Richness_C, type = "p", xlab = "Site", ylab = "Richness (Beetles)",)

##Diversity (Inverse Simpson's diversity index)----

#creating abundance matrixs
table(invert$Order)


diversity(x = invert[which(invert$Order=='Araneae'),] )


dim(invert[which(invert$Order=='Araneae'),])
###Araneae

Araneae <- invert[which(invert$Order=='Araneae'),]
head(Araneae);dim(Araneae)
Araneae_matrix <- table(Araneae$Site, Araneae$Morphospecies)
head(Araneae_matrix);dim(Araneae_matrix)

Araneae_diversity<-diversity(Araneae_matrix, index = "invsimpson")
Araneae_diversity<-data.frame(Site = names(Araneae_diversity), Diversity_A = as.numeric(Araneae_diversity))
head(Araneae_diversity);dim(Araneae_diversity)
variables <- merge(variables,Araneae_diversity, by = "Site",all.x = T)
head(variables);dim(variables)
variables$Diversity_A[is.na(variables$Diversity_A)] <- 0.000001


(table(variables$Diversity_A)/sum(table(variables$Diversity_A)))*100

dev.new(height=20,width=40,dpi=80,pointsize=14,noRStudioGD = T)
plot(as.factor(variables$Site), variables$Diversity_A, type = "p", xlab = "Site", ylab = "Diversity (spiders)",)

###Hemieptera

Hemiptera <- invert[which(invert$Order=='Hemiptera'),]
head(Hemiptera);dim(Hemiptera)
Hemiptera_matrix <- table(Hemiptera$Site, Hemiptera$Morphospecies)
head(Hemiptera_matrix);dim(Hemiptera_matrix)

Hemiptera_diversity<-diversity(Hemiptera_matrix, index = "invsimpson")
Hemiptera_diversity<-data.frame(Site = names(Hemiptera_diversity), Diversity_H = as.numeric(Hemiptera_diversity))
head(Hemiptera_diversity);dim(Hemiptera_diversity)
variables <- merge(variables,Hemiptera_diversity, by = "Site",all.x = T)
head(variables);dim(variables)
variables$Diversity_H[is.na(variables$Diversity_H)] <- 0.000001


(table(variables$Diversity_H)/sum(table(variables$Diversity_H)))*100

dev.new(height=20,width=40,dpi=80,pointsize=14,noRStudioGD = T)
plot(as.factor(variables$Site), variables$Diversity_H, type = "p", xlab = "Site", ylab = "Diversity (true bugs)",)


###Coleoptera

Coleoptera <- invert[which(invert$Order=='Coleoptera'),]
head(Coleoptera);dim(Coleoptera)
Coleoptera_matrix <- table(Coleoptera$Site, Coleoptera$Morphospecies)
head(Coleoptera_matrix);dim(Coleoptera_matrix)

Coleoptera_diversity<-diversity(Coleoptera_matrix, index = "invsimpson")
Coleoptera_diversity<-data.frame(Site = names(Coleoptera_diversity), Diversity_C = as.numeric(Coleoptera_diversity))
head(Coleoptera_diversity);dim(Coleoptera_diversity)
variables <- merge(variables,Coleoptera_diversity, by = "Site",all.x = T)
head(variables);dim(variables)
variables$Diversity_C[is.na(variables$Diversity_C)] <- 0.000001


(table(variables$Diversity_C)/sum(table(variables$Diversity_C)))*100

dev.new(height=20,width=40,dpi=80,pointsize=14,noRStudioGD = T)
plot(as.factor(variables$Site), variables$Diversity_C, type = "p", xlab = "Site", ylab = "Diversity (beetles)",)

##Community Composition----

###Araneae

Araneae_matrix2 <- as.matrix(Araneae_matrix)
str(Araneae_matrix2)
storage.mode(Araneae_matrix2) <- "numeric"  # Forces values to be numeric, not integer

pca_A<-prcomp(Araneae_matrix2,scale = T)
summary(pca_A)
pov_A <- summary(pca_A)$importance[2,]

sum(pov_A)

dev.new(height=20,width=20,dpi=80,pointsize=14,noRStudioGD = T)
plot(x=1:length(pov_A),y=pov_A,ylab="Propotion Varience Explained",xlab="Components",type="h",las=1)

###Hemieptera

Hemiptera_matrix2 <- as.matrix(Hemiptera_matrix)
str(Hemiptera_matrix2)
storage.mode(Hemiptera_matrix2) <- "numeric"  # Forces values to be numeric, not integer

pca_H<-prcomp(Hemiptera_matrix2,scale = T)
summary(pca_H)
pov_H <- summary(pca_H)$importance[2,]

sum(pov_H)

dev.new(height=20,width=20,dpi=80,pointsize=14,noRStudioGD = T)
plot(x=1:length(pov_H),y=pov_H,ylab="Propotion Varience Explained",xlab="Components",type="h",las=1)

###Coleoptera

Coleoptera_matrix2 <- as.matrix(Coleoptera_matrix)
str(Coleoptera_matrix2)
storage.mode(Coleoptera_matrix2) <- "numeric"  # Forces values to be numeric, not integer

pca_C<-prcomp(Coleoptera_matrix2,scale = T)
summary(pca_C)
pov_C <- summary(pca_C)$importance[2,]

sum(pov_C)

dev.new(height=20,width=20,dpi=80,pointsize=14,noRStudioGD = T)
plot(x=1:length(pov_C),y=pov_C,ylab="Propotion Varience Explained",xlab="Components",type="h",las=1)




#Beta Diversity -> araneae, hemi, coleoptera
  ###Araneae
  
  ###Hemieptera
  
  ###Coleoptera

#Functional abundance
    #araneae - hunting type
    #hemi - size and trophic
    #coleoptera - size and trophic

#Then need to check all these calculations for spatial autocorrelation
#Don't forget to remove all the coordinators etc. after this 

options(scipen = 999) #prevents r from automatically displaying large numbers with scientific notation

#Author: Rhiannon Bird
#Written under version R 4.5.1

#This script contains the data prepping, exploration and manipulation before moving into actual analysis


#Libraries----

library("dplyr")
library('vegan')
library("corrplot")
library("tidyr")
library("tibble")


#Data----
point <- read.csv("01. Data/Point_Data.csv")
head(point);dim(point)

field <- read.csv("01. Data/Survey_Data.csv")
head(field);dim(field)

morpho <- read.csv("01. Data/Morphospecies_Data.csv")
head(morpho);dim(morpho)

obs <- read.csv("01. Data/Observation_Data.csv")
head(obs);dim(obs)

head

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
morpho$Size[morpho$Size == "5-Oct"] <- "5-10"
morpho$Size[morpho$Size == "Oct-15"] <- "10-15"

#Since all levels greater then 10mm only have 1 morphospecies I'm combining together into  >10
morpho$Size[morpho$Size == ">30"] <- ">10"
morpho$Size[morpho$Size == "20-30"] <- ">10"
morpho$Size[morpho$Size == "15-20"] <- ">10"
morpho$Size[morpho$Size == "10-15"] <- ">10"

table(morpho$Size)
morpho$Size[is.na(morpho$Size)] <- "No_Size"

table(morpho$Trophic)
morpho$Trophic[morpho$Trophic == "Ominvore"] <- "Omnivore"


table(morpho$Hunting.Style)
morpho$Hunting.Style[is.na(morpho$Hunting.Style)] <- "Non-Predator"

#Correlation for environmental variables----

cordata <- data.frame(ID = point$Survey_Field,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position)

cordata <- merge(cordata,field, by = "ID")

head(cordata);dim(cordata)

cordata <- cordata %>% dplyr::select(Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X1km_Prop_Crops,NDVImean_Field,NDVIsum_1km,X1km_Rip_Prop)

cordata$Position[cordata$Position == 'Outer'] <- 1
cordata$Position[cordata$Position == 'Inner'] <- 2
cordata$Position <- as.numeric(cordata$Position)

str(cordata)#checking for any other character variables
cor <- cor(cordata,method = "spearman")
colnames(cor) <- c("Height", "Ground Cover", "Position","Day Sampled","Crop Age","Field Area","Crop 1000m","Field NDVI Mean","NDVI 1000m Sum","Riparian 1000m")
rownames(cor) <- c("Height", "Ground Cover", "Position","Day Sampled","Crop Age","Field Area","Crop 1000m","Field NDVI Mean","NDVI 1000m Sum","Riparian 1000m")

dev.new(height=8,width=8,dpi=80,pointsize=14,noRStudioGD = T)
corrplot::corrplot(cor,method="color",  
                   type="upper",addCoef.col = 'black',number.cex = 0.6)
head(cor)


#Merging data bases----

variables <- data.frame(Site = point$Site,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position,ID = point$Survey_Field)

variables <- merge(variables,field, by = "ID")

head(variables);dim(variables)

variables <- variables %>% dplyr::select(ID,Site,Field,Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X1km_Prop_Crops, NDVImean_Field,NDVIsum_1km,X1km_Rip_Prop)
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


#total observations = 2329

setdiff(unique(morpho$Morphospecies), unique(invert$Morphospecies))

#need to remove two from morpho to make sure that its got the same species are were observed
#double checking they aren't in the observations
invert[invert$Morphospecies == 'Brown_Weevil',] 
invert[invert$Morphospecies == 'Small_Ant',]
#They aren't so can remove from morphospecies list

morpho <- morpho[morpho$Morphospecies != 'Brown_Weevil',]
morpho <- morpho[morpho$Morphospecies != 'Small_Ant',]
head(morpho);dim(morpho)

#Checking the groups----
head(invert)
length(point$Site) #total of 220 points

table(invert$Order) 

table(invert$Order,invert$ID)

table(invert$Order,invert$Trophic)
table(invert$Order,invert$Size)
table(invert$Order,invert$Hunting.Style)

length(unique(invert$Morphospecies))
length(invert$Morphospecies)
table(invert$ID_Level)

length(unique(invert$Family))
length(unique(invert$Genus))
#Prep for taxonomic Modelling----

names(invert)[names(invert) == "ID"] <- "Site"

TaxModel <- data.frame(Site = variables$Site)
head(TaxModel);dim(TaxModel)

##Species Richness----

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

#filtered out morphospecies that are an unknown order

length(which(is.na(invert$Order)))

invert_filtered <- invert[-which(is.na(invert$Order)), ]

head(invert_filtered);dim(invert_filtered) 
dim(invert);dim(invert_filtered)

##L Data (Site x Species)----

head(invert_filtered);dim(invert_filtered)


#Count occurrences of each morphospecies at each site
L_table <- invert_filtered %>%
  group_by(Site, Morphospecies) %>%
  summarise(abundance = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Morphospecies, 
              values_from = abundance, 
              values_fill = 0)
head(L_table);dim(L_table)
L_table[1,1]

#Convert to data frame with sites as row names
L_matrix <- L_table %>%
  column_to_rownames("Site") %>%
  as.data.frame()
head(L_matrix);dim(L_matrix)

##Q Data (Species x Trait)----

head(invert_filtered)

#Get unique morphospecies with traits
Q_table <- invert_filtered %>%
  select(Morphospecies, Order, Size, Hunting.Style, Trophic) %>%
  distinct()
head(Q_table);dim(Q_table)

Q_matrix <- Q_table %>%
  column_to_rownames("Morphospecies") %>%
  as.data.frame()
head(Q_matrix);dim(Q_matrix)

table(Q_matrix$Size)


Q_matrix$Size[Q_matrix$Size == "Unknown"] <- "No_Size"


##R Data (Site x Environmental)----

head(variables);dim(variables)

R_table <- variables %>%
  select(-ID, -Field)
head(R_table);dim(R_table)

#to make the next part transform correctly
R_table <- as.data.frame(R_table)
rownames(R_table) <- NULL

R_matrix <- R_table %>%
  column_to_rownames("Site") %>%
  as.data.frame()

##Checking that it worked----

#first up is sites

sites_L <- rownames(L_matrix)
sites_R <- rownames(R_matrix)

if (!all(sites_L %in% sites_R)) {
  warning("Some sites in L table are not in R table")
  print(setdiff(sites_L, sites_R))
} #Needs fixing

if (!all(sites_R %in% sites_L)) {
  warning("Some sites in R table are not in L table")
  print(setdiff(sites_R, sites_L))
} #Needs fixing

#fix missing sites and make sure sites are in the same order for both
common_sites <- intersect(sites_L, sites_R)
L_matrix <- L_matrix[common_sites, ]
R_matrix <- R_matrix[common_sites, ]

sites_L.2 <- rownames(L_matrix)
sites_R.2 <- rownames(R_matrix)

if (!all(sites_R.2 %in% sites_L.2)) {
  warning("Some sites in R table are not in L table")
  print(setdiff(sites_R.2, sites_L.2))
} #fixed they match

#Now checking morphospecies 

species_L <- colnames(L_matrix)
species_Q <- rownames(Q_matrix)

if (!all(species_L %in% species_Q)) {
  warning("Some morphospecies in L table are not in Q table")
  print(setdiff(species_L, species_Q))
} #looks good

if (!all(species_Q %in% species_L)) {
  warning("Some morphospecies in Q table are not in L table")
  print(setdiff(species_Q, species_L))
} #looks good

#match Q order to L
Q_matrix <- Q_matrix[species_L, ]

##Fixing missing values----

any(is.na(L_matrix)) #No NA's
any(is.na(R_matrix)) #No NA's
any(is.na(Q_matrix)) #No NA's

head(L_matrix[,1:10]);dim(L_matrix)
head(R_matrix[,1:10]);dim(R_matrix)
head(Q_matrix);dim(Q_matrix)

#Traits need to be factors

Q_matrix[] <- lapply(Q_matrix, as.factor)
str(Q_matrix)

levels(Q_matrix$Size)
Q_matrix$Size <- factor(Q_matrix$Size, 
                        levels = c("0-2.5", "2.5-5", 
                                   "5-10", ">10", "No_Size"),
                        ordered = TRUE)
levels(Q_matrix$Size)

#then also for position
R_matrix$Position <- as.factor(R_matrix$Position)

#END


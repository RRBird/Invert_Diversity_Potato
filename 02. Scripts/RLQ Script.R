#Libraries----

library("ade4")
library("dplyr")
library("tidyverse")

#DOING A BASIC TEST TO SEE IF IT WILL WORK SO NEED TO GO THROUGH THROULY AND MAKE SURE IT'S CORRECT




#Data---
#Prepping RLQ Data----

##R Matrix (environment)----

head(variables);dim(variables)

environment <- variables %>% dplyr::select(-ID,-Field,-Position,-Day_Sampled,-Crop_Age_Days)
head(environment);dim(environment)

#not sure if my design variables should be included in my matrix so I've done two

environment2 <- variables %>% dplyr::select(-ID,-Field)
head(environment2);dim(environment2)

##L Matrix (Species abundance) ----

table(invert$Morphospecies, invert$Site)
table(invert$Order)


prespecies <- invert
head(invert);dim(invert)
prespecies <- invert %>% filter(Order %in% c("Araneae", "Coleoptera","Diptera","Hemiptera"))
head(prespecies);dim(prespecies)
table(prespecies$Order)


species <- prespecies %>%
  count(Morphospecies, Site) %>%
  pivot_wider(names_from = Site, values_from = n, values_fill = 0)

species <- table(prespecies$Morphospecies, prespecies$Site)
head(species)
str(species)
species <- data.frame(species)
colnames(species)[1] <- "Morphospecies"
colnames(species)[2] <- "Site"
colnames(species)[3] <- "Abundance"

head(species);dim(species)


## Q Matrix (traits)----

head(morpho);dim(morpho)

pretrait <- morpho %>% filter(Order %in% c("Araneae", "Coleoptera","Diptera","Hemiptera"))
head(pretrait);dim(pretrait)

pretrait$Size <- addNA(pretrait$Size) 
levels(pretrait$Size)[is.na(levels(pretrait$Size))] <- "No_Size"

pretrait$Trophic[is.na(pretrait$Trophic)] <- "Unknown"
pretrait$Hunting.Style[is.na(pretrait$Hunting.Style)] <- "No_Hunt"

traits <- pretrait %>% dplyr::select(Order,Trophic,Hunting.Style, Size)
head(traits);dim(traits)

str(traits) 



#correspondence analysis on matrix L
coa1 <- dudi.coa(species, scannf = F)
head(coa1)

pca.traits <- dudi.pca(traits,row.w = coa1$cw,scannf=F)







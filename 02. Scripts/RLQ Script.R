#Libraries----

library("ade4")

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

species <- invert %>%
  count(Morphospecies, Site) %>%
  pivot_wider(names_from = Site, values_from = n, values_fill = 0)

species <- table(invert$Morphospecies, invert$Site)
head(species)
str(species)
species <- data.frame(species)
colnames(species)[1] <- "Morphospecies"
colnames(species)[2] <- "Site"
colnames(species)[3] <- "Abundance"



## Q Matrix (traits)----

head(morpho);dim(morpho)

traits <- morpho %>% dplyr::select(Morphospecies,Trophic, Size, Wings)
head(traits);dim(traits)

str(traits) 
#might need to turn traits into factors I'm not sure

traits_clean <- na.omit(traits)
head(traits_clean);dim(traits_clean)


#Turns out this analysis doesn't allow for any NA's so need to figure that out. Would adding another category instead of NA which is Unknown fix that problem?



#correspondence analysis on matrix L
coa1 <- dudi.coa(species, scannf = F)
head(coa1)

pca.traits <- dudi.pca(traits,row.w = coa1$cw,scannf=F)







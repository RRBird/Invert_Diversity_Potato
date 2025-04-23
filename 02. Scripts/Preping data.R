
#Libraries


#Data----
point <- read.csv("01. Data/Point_Data.csv")
head(point);dim(point)

field <- read.csv("01. Data/Survey_Data.csv")
head(field);dim(field)

morpho <- read.csv("01. Data/Morphospecies_Data.csv")
head(morpho);dim(morpho)

obs <- read.csv("01. Data/Observation_Data.csv")
head(obs);dim(obs)


#prepping data
##remove unneeded columns
##Remove lepedoptera and wasps
##remove first survey


##Calculate for modelling
  #Richness -> araneae, hemi, coleoptera
  #Diversity -> araneae, hemi, coleoptera
  #Community Comp -> araneae, hemi, coleoptera
  #Beta Diversity -> araneae, hemi, coleoptera
  #Functional abundance
    #araneae - hunting type
    #hemi - size and trophic
    #coleoptera - size and trophic

#Then need to check all these calculations for spatial autocorrelation
#Don't forget to remove all the cordinates etc. after this 
  
#Libraries----


library("AICcmodavg")
library('lme4')


#modelling and predictions etc will go here


#Richness -> araneae, hemi, coleoptera
  #Diversity -> araneae, hemi, coleoptera
  #Community Comp -> araneae, hemi, coleoptera
  #Beta Diversity -> araneae, hemi, coleoptera
  #Functional abundance
    #araneae - hunting type
    #hemi - size and trophic
    #coleoptera - size and trophic




#Model Trial

variables$Sample_Day_Scale <- scale(variables$Day_Sampled)
variables$Crop_Age_Scale <- scale(variables$Crop_Age_Days)

Null_Rich <- glmer(Richness_A~ 1 + (1|Field),data = variables,family = poisson)

Position <- glmer(Richness_A~ Position + (1|Field),data = variables,family = poisson) 
Age <- glmer(Richness_A~ Crop_Age_Scale + (1|Field),data = variables,family = poisson)
Day <- glmer(Richness_A~ Sample_Day_Scale + (1|Field),data = variables,family = poisson)

Position_Age1 <- glmer(Richness_A~ Position + Crop_Age_Scale + (1|Field),data = variables,family = poisson) 
Position_Age2 <- glmer(Richness_A~ Position * Crop_Age_Scale + (1|Field),data = variables,family = poisson) 

Position_Day1 <- glmer(Richness_A~ Position + Day_Sampled + (1|Field),data = variables,family = poisson) 
Position_Day2 <- glmer(Richness_A~ Position * Sample_Day_Scale + (1|Field),data = variables,family = poisson)


aictab(list(null = Null_Rich,position = Position, age = Age, day = Day, position_Age = Position_Age1, positionxage = Position_Age2,position_day = Position_Day1, positionxday = Position_Day2))

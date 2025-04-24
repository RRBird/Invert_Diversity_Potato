#Libraries----

library(AICcmodavg)
library(lme4)
library(glmmADMB)


#modelling and predictions etc here


#Richness -> araneae, hemi, coleoptera
  #Diversity -> araneae, hemi, coleoptera
  #Community Comp -> araneae, hemi, coleoptera
  #Beta Diversity -> araneae, hemi, coleoptera
  #Functional abundance
    #araneae - hunting type
    #hemi - size and trophic
    #coleoptera - size and trophic





#Just some trial models




Null_Rich <- glmmadmb(Richness_A~ 1 + (1|Field),data = variables,family = nbinom)

Position <- glmer(Richness_A~ Position + (1|Field),data = variables,family = binomial) 
Age <- glmer(Richness_A~ Crop_Age_Days + (1|Field),data = variables,family = binomial)
Day <- glmer(Richness_A~ Day_Sampled + (1|Field),data = variables,family = binomial)

Position_Age1 <- glmer(Richness_A~ Position + Crop_Age_Days + (1|Field),data = variables,family = binomial) 
Position_Age2 <- glmer(Richness_A~ Position * Crop_Age_Days + (1|Field),data = variables,family = binomial) 

Position_Day1 <- glmer(Richness_A~ Position + Day_Sampled + (1|Field),data = variables,family = binomial) 
Position_Day2 <- glmer(Richness_A~ Position * Day_Sampled + (1|Field),data = variables,family = binomial)

Age_Day1 <- glmer(Richness_A~ Crop_Age_Days + Day_Sampled + (1|Field),data = variables,family = binomial) 
Age_Day2 <- glmer(Richness_A~ Crop_Age_Days * Day_Sampled + (1|Field),data = variables,family = binomial) 



aictab(list(null = Null_Dom,position = Position, age = Age, day = Day, position_Age = Position_Age1, positionxage = Position_Age2,position_day = Position_Day1, positionxday = Position_Day2, Age_Day = Age_Day1, Age*Day = Age_Day2 ))
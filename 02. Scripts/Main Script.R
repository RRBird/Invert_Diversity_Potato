#Libraries----


library("AICcmodavg")
library('lme4')
library("glmmTMB")

#Richness Modelling----

##Step 1 - Design Variables----

#position, age and day - all the various combinations of these 

#try to write a loop for this first part - can put the results into a list so then I can do list[1] and see the modelling results for just the first grouping
head(ModelRich2);dim(ModelRich2)
str(ModelRich2)
#loop is struggling with the size column names so will just update them
colnames(ModelRich2)[12] <- "Size_1"
colnames(ModelRich2)[13] <- "Size_2"
colnames(ModelRich2)[14] <- "Size_3"
colnames(ModelRich2)[15] <- "Size_4"


rich_names <- colnames(ModelRich2)[2:19]
length(rich_names)
Richlist1 <- list()
allsum <- list()

for (i in rich_names) {
#run models and collect summaries
    null <- glmmTMB(as.formula(paste(i, "~ 1 + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum [[1]] <- summary(null)
    P <- glmmTMB(as.formula(paste(i, "~ Position + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum [[2]] <- summary(P)
    A <- glmmTMB(as.formula(paste(i, "~ Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum [[3]] <- summary(A)
    D <- glmmTMB(as.formula(paste(i, "~ Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[4]] <- summary(D)
    PA <- glmmTMB(as.formula(paste(i, "~ Position + Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[5]] <- summary(PA)
    PD <- glmmTMB(as.formula(paste(i, "~ Position + Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[6]] <- summary(PD)
    DA <- glmmTMB(as.formula(paste(i, "~ Day_Sampled + Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[7]] <- summary(DA)
    PxA <- glmmTMB(as.formula(paste(i, "~ Position * Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[8]] <- summary(PxA)
    PxD <- glmmTMB(as.formula(paste(i, "~ Position * Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[9]] <- summary(PxD)
    DxA <- glmmTMB(as.formula(paste(i, "~ Day_Sampled * Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[10]] <- summary(DxA)
    PAD <- glmmTMB(as.formula(paste(i, "~ Position + Crop_Age_Days + Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[11]] <- summary(PAD)
    PxAD <- glmmTMB(as.formula(paste(i, "~ Position * Crop_Age_Days + Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[12]] <- summary(PxAD)
    PxDA <- glmmTMB(as.formula(paste(i, "~ Position * Day_Sampled + Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[13]] <- summary(PxDA)
    PAxD <- glmmTMB(as.formula(paste(i, "~ Position + Crop_Age_Days * Day_Sampled + (1 | Field)")), data = ModelRich2)
    allsum[[14]] <- summary(PAxD)
    PxAxD <- glmmTMB(as.formula(paste(i, "~ Position * Crop_Age_Days * Day_Sampled + (1 | Field)")), data = ModelRich2)
    allsum[[15]] <- summary(PxAxD)
    
  #collect models
  tempmodlist <- list("null" = null, "P" = P, "A" = A, "D" = D,
                      "PA" = PA, "PD" = PD, "DA" = DA, 
                      "PxA" = PxA, "PxD" = PxD, "DxA" = DxA, 
                      "PAD" = PAD, "PxAD" = PxAD, "PxDA" = PxDA, 
                      "PAxD" = PAxD, "PxAxD" = PxAxD)
  
  #check which models didn't coverge/had issue to exclude from AIC
  has_na <- sapply(allsum, function(model) is.na(model$AICtab[1]))
  na_models <- allsum[has_na]
  cleaned_models <- tempmodlist[!has_na] #remove models with NA as AIC
  
 Richlist1 [[i]] <- aictab(cleaned_models,modnames = NULL)
 
}

length(Richlist1)

Richlist1[[1]]

##Step 2 - Environmental Variables----

#can we do a similar thing but would need to be able to shift the variables already there for each 
#Will probably need to store the formula for each group that needs to be included in a list or something -> Group = Position * Age (or something)
#Possibly a for loop within a for loop???

#Can i have t = 1 (add to t each time)
#then inside model richformulas[[1]]
#e.g. glmmTMB(as.formula(paste(i, "~ richformulas[[1]] + Envrio + (1 | Field)")), family = nbinom2, data = ModelRich2)

#small trial
trial <- list()
trial [[1]] <- "Day_Sampled * Crop_Age_Days"
glmmTMB(All ~ trial [[1]] + GC + (1 | Field), family = nbinom2, data = ModelRich2)
#that didn't work so see if you can fiund a way around it

#Don't forget need to run the model from step 1 and all the environmental models 


##Step 3 - Check for spatial autocorrelation----


##Step 4 - Visualisation----

###Main Figures----

###Supporting Figures----


#Abundance Modelling----

#Diversity Modelling----


#Binomial Modelling----
















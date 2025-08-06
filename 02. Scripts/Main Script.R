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
    P+A <- glmmTMB(as.formula(paste(i, "~ Position + Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[5]] <- summary(P+A)
    P+D <- glmmTMB(as.formula(paste(i, "~ Position + Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[6]] <- summary(P+D)
    D+A <- glmmTMB(as.formula(paste(i, "~ Day_Sampled + Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[7]] <- summary(D+A)
    PxA <- glmmTMB(as.formula(paste(i, "~ Position * Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[8]] <- summary(PxA)
    PxD <- glmmTMB(as.formula(paste(i, "~ Position * Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[9]] <- summary(PxD)
    DxA <- glmmTMB(as.formula(paste(i, "~ Day_Sampled * Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[10]] <- summary(DxA)
    
    P+A+D <- glmmTMB(as.formula(paste(i, "~ Position + Crop_Age_Days + Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[11]] <- summary(DxA)
    PxA+D <- glmmTMB(as.formula(paste(i, "~ Position * Crop_Age_Days + Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[12]] <- summary(DxA)
    PxD+A <- glmmTMB(as.formula(paste(i, "~ Position * Day_Sampled + Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[13]] <- summary(DxA)
    P+AxD <- glmmTMB(as.formula(paste(i, "~ Position + Crop_Age_Days * Day_Sampled + (1 | Field)")), data = ModelRich2)
    allsum[[14]] <- summary(DxA)
    PxAxD <- glmmTMB(as.formula(paste(i, "~ Position * Crop_Age_Days * Day_Sampled + (1 | Field)")), data = ModelRich2)
    allsum[[15]] <- summary(DxA)
    
  #collect models
  tempmodlist <- list("null" = null, "P" = P, "A" = A, "D" = D,
                      "P+A" = P+A, "P+D" = P+D, "D+A" = D+A, 
                      "PxA" = PxA, "PxD" = PxD, "DxA" = DxA, 
                      "P+A+D"=P+A+D, "PxA+D"=PxA+D, "PxD+A"=PxD+A, 
                      "P+AxD" = P+AxD, "PxAxD" = PxAxD)
  
  #check which models didn't coverge/had issue to exclude from AIC
  has_na <- sapply(allsum, function(model) is.na(model$AICtab[1]))
  na_models <- allsum[has_na]
  cleaned_models <- tempmodlist[!has_na] #remove models with NA as AIC
  
 Richlist1 [[i]] <- aictab(cleaned_models,modnames = NULL)
 
}
length(Richlist1)



#Sort of working --> only 10 of the richness groups there

#ORGINAL FOR LOOP----

for (i in rich_names) {
  null <- glmmTMB(i ~ 1 + (1 | Field), family = nbinom2, data = ModelRich2)
  P <- glmmTMB(i ~ Position + (1 | Field), family = nbinom2, data = ModelRich2)
  A <- glmmTMB(i ~ Crop_Age_Days + (1 | Field), family = nbinom2, data = ModelRich2)
  D <- glmmTMB(i ~ Day_Sampled + (1 | Field), family = nbinom2, data = ModelRich2)
  P+A <- glmmTMB(i ~ Position + Crop_Age_Days + (1 | Field), family = nbinom2, data = ModelRich2)
  P+D <- glmmTMB(i ~ Position + Day_Sampled + (1 | Field), family = nbinom2, data = ModelRich2)
  D+A <- glmmTMB(i ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = nbinom2, data = ModelRich2)
  PxA <- glmmTMB(i ~ Position * Crop_Age_Days + (1 | Field), family = nbinom2, data = ModelRich2)
  PxD <- glmmTMB(i ~ Position * Day_Sampled + (1 | Field), family = nbinom2, data = ModelRich2)
  DxA <- glmmTMB(i ~ Day_Sampled * Crop_Age_Days + (1 | Field), family = nbinom2, data = ModelRich2)
  P+A+D <- glmmTMB(i ~ Position + Crop_Age_Days + Day_Sampled + (1 | Field), family = nbinom2, data = ModelRich2)
  PxA+D <- glmmTMB(i ~ Position * Crop_Age_Days + Day_Sampled + (1 | Field), family = nbinom2, data = ModelRich2)
  PxD+A <- glmmTMB(i ~ Position * Day_Sampled + Crop_Age_Days + (1 | Field), family = nbinom2, data = ModelRich2)
  P+AxD <- glmmTMB(i ~ Position + Crop_Age_Days * Day_Sampled + (1 | Field), family = nbinom2, data = ModelRich2)
  PxAxD <- glmmTMB(i ~ Position * Crop_Age_Days * Day_Sampled + (1 | Field), family = nbinom2, data = ModelRich2)
  
  tempmodlist <- list(null, P, A, D, P+A, P+D, D+A, PxA, PxD, DxA, P+A+D, PxA+D, PxD+A, P+AxD, PxAxD)
  
  Richlist1 [[t]] <- aictab(tempmodlist)
  
  t <- t+1
}


##Step 2 - Environmental Variables----


##Step 3 - Check for spatial autocorrelation----


##Step 4 - Visualisation----

###Main Figures----

###Supporting Figures----


#Abundance Modelling----

#Diversity Modelling----


#Binomial Modelling----
















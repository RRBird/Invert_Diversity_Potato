#Libraries----


library("AICcmodavg")
library('lme4')
library("glmmTMB")

#Richness Modelling----

##Step 1 - Design Variables----

#position, age and day - all the various combinations of these 

#try to write a loop for this first part - can put the results into a list so then I can do list[1] and see the modelling results for just the first grouping
head(ModelRich2);dim(ModelRich2)

rich_names <- colnames(ModelRich2)[2:19]
Richlist1 <- list()

t <- 1
for (i in rich_names) {
  null <- glmmTMB(as.formula(paste(i, "~ 1 + (1 | Field)")), family = nbinom2, data = ModelRich2)
  P <- glmmTMB(as.formula(paste(i, "~ Position + (1 | Field)")), family = nbinom2, data = ModelRich2)
  A <- glmmTMB(as.formula(paste(i, "~ Crop_Age_Days + (1 | Field)")), family = nbinom2, data = ModelRich2)
  D <- glmmTMB(as.formula(paste(i, "~ Day_Sampled + (1 | Field)")), family = nbinom2, data = ModelRich2)
  
  tempmodlist <- list("null" = null, "P" = P, "A" = A, "D" = D)
  
 Richlist1 [[t]] <- aictab(tempmodlist,modnames = NULL)
 
 t <- t+1
}
warnings()
t
Richlist1[[1]]

glmmTMB(Hawking ~ Position + (1 | Field), family = nbinom2, data = ModelRich2)

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
















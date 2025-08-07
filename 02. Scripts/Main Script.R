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

ModelRich2$Age_Scaled <- scale(ModelRich$Crop_Age_Days)
ModelRich2$Day_Scaled <- scale(ModelRich$Day_Sampled)


rich_names <- colnames(ModelRich2)[2:19]
length(rich_names)
Richlist1 <- list()

for (i in rich_names) {
  allsum <- list()
  
  cat("=== Processing group:", i, "===\n") #track where loop is up to
  
#run models and collect summaries
    null <- glmmTMB(as.formula(paste(i, "~ 1 + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum [[1]] <- summary(null)
    P <- glmmTMB(as.formula(paste(i, "~ Position + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum [[2]] <- summary(P)
    A <- glmmTMB(as.formula(paste(i, "~ Age_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum [[3]] <- summary(A)
    D <- glmmTMB(as.formula(paste(i, "~ Day_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[4]] <- summary(D)
    PA <- glmmTMB(as.formula(paste(i, "~ Position + Age_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[5]] <- summary(PA)
    PD <- glmmTMB(as.formula(paste(i, "~ Position + Day_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[6]] <- summary(PD)
    DA <- glmmTMB(as.formula(paste(i, "~ Day_Scaled + Age_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[7]] <- summary(DA)
    PxA <- glmmTMB(as.formula(paste(i, "~ Position * Age_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[8]] <- summary(PxA)
    PxD <- glmmTMB(as.formula(paste(i, "~ Position * Day_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[9]] <- summary(PxD)
    DxA <- glmmTMB(as.formula(paste(i, "~ Day_Scaled * Age_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[10]] <- summary(DxA)
    PAD <- glmmTMB(as.formula(paste(i, "~ Position + Age_Scaled + Day_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[11]] <- summary(PAD)
    PxAD <- glmmTMB(as.formula(paste(i, "~ Position * Age_Scaled + Day_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[12]] <- summary(PxAD)
    PxDA <- glmmTMB(as.formula(paste(i, "~ Position * Day_Scaled + Age_Scaled + (1 | Field)")), family = nbinom2, data = ModelRich2)
    allsum[[13]] <- summary(PxDA)
    PAxD <- glmmTMB(as.formula(paste(i, "~ Position + Age_Scaled * Day_Scaled + (1 | Field)")), data = ModelRich2,family = nbinom2,)
    allsum[[14]] <- summary(PAxD)
    
  #collect models
  tempmodlist <- list("null" = null, "P" = P, "A" = A, "D" = D,
                      "PA" = PA, "PD" = PD, "DA" = DA, 
                      "PxA" = PxA, "PxD" = PxD, "DxA" = DxA, 
                      "PAD" = PAD, "PxAD" = PxAD, "PxDA" = PxDA, 
                      "PAxD" = PAxD)
  
  has_issues <- sapply(tempmodlist, function(model) {
    if (is.null(model)) return(TRUE) #model failed
    if (model$fit$convergence != 0) return(TRUE) #convergence issues
    if (is.na(AIC(model))) return(TRUE) #does it have an AIC (not having one will cause all AIC in table to be NA)
    return(FALSE)
  }) 
  
  cleaned_models <- tempmodlist[!has_issues]
  cat("Models with issues:", sum(has_issues), "\n")
  cat("Problem models:", names(tempmodlist)[has_issues], "\n")
  
 Richlist1 [[i]] <- aictab(cleaned_models,modnames = NULL)
 
}

length(Richlist1)

Richlist1

#check which models didn't coverge/had issue to exclude from AIC
has_na <- sapply(allsum, function(model) is.na(model$AICtab[1]))
na_models <- allsum[has_na]
cleaned_models <- tempmodlist[!has_na] #remove models with NA as AIC

##Step 2 - Environmental Variables----

rich_names
Rich_design <- list(
  All = "Age_Scaled",
  Predator = "Day_Scaled + Age_Scaled",
  Herbivore = "Position + Day_Scaled",   
  Omnivore = "Age_Scaled",
  Fungivore = "Age_Scaled",     
  Hematophagous = NULL,
  Web = "Position * Day_Scaled + Age_Scaled",          
  Active_Hunting = "Day_Scaled * Age_Scaled",
  Ambush_Hunting = "Position + Age_Scaled",
  Hawking = "Position + Day_Scaled",      
  Size_1 = "Position + Age_Scaled + Day_Scaled",
  Size_2 = NULL,      
  Size_3 = "Position + Age_Scaled * Day_Scaled",
  Size_4 = "Age_Scaled",
  Always_Winged = NULL,
  Develops_Wings = "Age_Scaled",
  Wingless = "Position * Age_Scaled + Day_Scaled",
  Polymorphic = "Day_Scaled + Age_Scaled *")
head(Rich_design)

Richlist2 <- list()


for (i in rich_names) {
  allsum <- list()
  
  design_formula_part <- Rich_design[[i]]
  
  cat("=== Processing group:", i, "===\n") #track where loop is up to
  cat("Design formula part:", design_formula_part, "\n") #making sure it's taken the right formula
  
  #Handle groups without null model as top in step 1
  if (is.null(design_formula_part) || design_formula_part == "") {
    # No design variables - intercept only model
    base_formula <- paste(i, "~ 1 + (1 | Field)")} else {
    # Has design variables
    base_formula <- paste(i, "~", design_formula_part, "+ (1 | Field)")}
  base_model <- glmmTMB(as.formula(base_formula),family = nbinom2, data = ModelRich2)
  allsum[[1]] <- summary(base_model)

  
  if (is.null(design_formula_part) || design_formula_part == "") {
    height_formula <- paste(i, "~ Height + (1 | Field)")} else {
      height_formula <- paste(i, "~", design_formula_part, "+ Height + (1 | Field)")}
  height_model <- glmmTMB(as.formula(height_formula), family = nbinom2, data = ModelRich2)
  allsum[[2]] <- summary(height_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    gc_formula <- paste(i, "~ GC  + (1 | Field)")} else {
      gc_formula <- paste(i, "~", design_formula_part, "+ GC + (1 | Field)")}
  gc_model <- glmmTMB(as.formula(gc_formula), family = nbinom2, data = ModelRich2)
  allsum[[3]] <- summary(gc_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    Fsize_formula <- paste(i, "~ Field_Area_m2  + (1 | Field)")} else {
      Fsize_formula <- paste(i, "~", design_formula_part, "+ Field_Area_m2 + (1 | Field)")}
  Fsize_model <- glmmTMB(as.formula(Fsize_formula), family = nbinom2, data = ModelRich2)
  allsum[[4]] <- summary(Fsize_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    water_formula <- paste(i, "~ X1km_Prop_Water  + (1 | Field)")} else {water_formula <- paste(i, "~", design_formula_part, "+ X1km_Prop_Water + (1 | Field)")}
  water_model <- glmmTMB(as.formula(water_formula), family = nbinom2, data = ModelRich2)
  allsum[[5]] <- summary(water_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    NDVIf_formula <- paste(i, "~ NDVImean_Field  + (1 | Field)")} else {NDVIf_formula <- paste(i, "~", design_formula_part, "+ NDVImean_Field + (1 | Field)")}
  NDVIf_model <- glmmTMB(as.formula(NDVIf_formula), family = nbinom2, data = ModelRich2)
  allsum[[6]] <- summary(NDVIf_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    NDVI1km_formula <- paste(i, "~ NDVIsum_1km  + (1 | Field)")} else {
      NDVI1km_formula <- paste(i, "~", design_formula_part, "+ NDVIsum_1km + (1 | Field)")}
  NDVI1km_model <- glmmTMB(as.formula(NDVI1km_formula), family = nbinom2, data = ModelRich2)
  allsum[[7]] <- summary(NDVI1km_model)
  
  
  #collect models
  tempmodlist <- list("base" = base_model, "height" = height_model, 
                      "GC" = gc_model, "Feild Size" = Fsize_model,
                      "water" = water_model,"Field NDVI" = NDVIf_model,
                      "NDVI 1km" = NDVI1km_model)
  
  has_issues <- sapply(tempmodlist, function(model) {
    if (is.null(model)) return(TRUE) #model failed
    if (model$fit$convergence != 0) return(TRUE) #convergence issues
    if (is.na(AIC(model))) return(TRUE) #does it have an AIC (not having one will cause all AIC in table to be NA)
    return(FALSE)
  }) 
  
  cleaned_models <- tempmodlist[!has_issues]
  cat("Models with issues:", sum(has_issues), "\n")
  cat("Problem models:", names(tempmodlist)[has_issues], "\n")
  
  Richlist2 [[i]] <- aictab(cleaned_models,modnames = NULL)
  
}

length(Richlist2)

Richlist2

###Final Richness models ----

Rich_All <- glmmTMB(All ~ Age_Scaled + Field_Area_m2 + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_All)

Rich_Pred <- glmmTMB(Predator ~ Day_Scaled + Age_Scaled + X1km_Prop_Water + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Pred)

Rich_Herb <- glmmTMB(Herbivore ~ Position + Day_Scaled + X1km_Prop_Water + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Herb)

Rich_Omni <- glmmTMB(Omnivore ~ Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Omni)

Rich_Fung <- glmmTMB(Fungivore ~ Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Fung)

Rich_Hema <- glmmTMB(Hematophagous ~ NDVIsum_1km + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Hema)

Rich_Web <- glmmTMB(Web ~ Position * Day_Scaled + Age_Scaled + GC + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Web)

Rich_Active <- glmmTMB(Active_Hunting ~ Day_Scaled * Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Active)

Rich_Ambush <- glmmTMB(Ambush_Hunting ~ Position + Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Ambush)

Rich_Hawk <- glmmTMB(Hawking ~ Position + Day_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Hawk)

Rich_Size1 <- glmmTMB(Size_1 ~ Position + Age_Scaled + Day_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Size1)

Rich_Size2 <- glmmTMB(Size_2 ~ Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Size2)

Rich_Size3 <- glmmTMB(Size_3 ~ Position + Age_Scaled * Day_Scaled + X1km_Prop_Water + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Size3)

Rich_Size4 <- glmmTMB(Size_4 ~ Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Size4)

Rich_DevW <- glmmTMB(Develops_Wings ~ Age_Scaled + NDVIsum_1km + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_DevW)

Rich_Wless <- glmmTMB(Wingless ~ Position * Age_Scaled + Day_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Wless)

Rich_Poly <- glmmTMB(Polymorphic ~ Day_Scaled + Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Poly)

##Step 3 - Check for spatial autocorrelation----


##Step 4 - Visualisation----

###Main Figures----

###Supporting Figures??----


#Abundance Modelling----

#Diversity Modelling----


#Binomial Modelling----
















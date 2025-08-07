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
  allsum <- list()
  
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

rich_names
Rich_design <- list(
  All = "Day_Sampled * Crop_Age_Days",
  Predator = "Day_Sampled * Crop_Age_Days",
  Herbivore = "Position + Day_Sampled",   
  Omnivore = "Position + Crop_Age_Days * Day_Sampled",
  Fungivore = "Position + Crop_Age_Days * Day_Sampled",     
  Hematophagous = "Position + Crop_Age_Days * Day_Sampled",
  Web = "Position * Day_Sampled + Crop_Age_Days",          
  Active_Hunting = "Day_Sampled",
  Ambush_Hunting = "Position + Crop_Age_Days * Day_Sampled",
  Hawking = "Position + Crop_Age_Days * Day_Sampled",      
  Size_1 = "Day_Sampled + Crop_Age_Days",
  Size_2 = NULL,      
  Size_3 = "Day_Sampled",
  Size_4 = "Position + Crop_Age_Days * Day_Sampled",
  Always_Winged = "Position + Crop_Age_Days * Day_Sampled",
  Develops_Wings = "Crop_Age_Days",
  Wingless = "Position * Crop_Age_Days + Day_Sampled",
  Polymorphic = "Position * Crop_Age_Days * Day_Sampled")
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
  
  #check which models didn't coverge/had issue to exclude from AIC
  has_na <- as.logical(sapply(allsum, function(model) is.na(model$AICtab[1])))
  na_models <- allsum[has_na]
  cleaned_models <- tempmodlist[!has_na] #remove models with NA as AIC
  
  Richlist2 [[i]] <- aictab(cleaned_models,modnames = NULL)
  
}

length(Richlist2)

Richlist2[[1]]


#Loop is working currently but might need to rescale variables to make sure that models will converge??? 
#need to work out what my system is for this?????

##Step 3 - Check for spatial autocorrelation----


##Step 4 - Visualisation----

###Main Figures----

###Supporting Figures----


#Abundance Modelling----

#Diversity Modelling----


#Binomial Modelling----
















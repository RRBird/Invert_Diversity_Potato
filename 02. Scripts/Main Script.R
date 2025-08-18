
#This script has the main analysis for richness, occurance and diversity modelling 

#Libraries----


library("AICcmodavg")
library('lme4')
library("glmmTMB")
library("dplyr")
library("DHARMa")


#RICHNESS MODELLING----

##Step 1 - Design Variables----

#position, age and day - all the various combinations of these 


head(ModelRich2);dim(ModelRich2)
str(ModelRich2)

#R is struggling with the size column names so will just update them
colnames(ModelRich2)[12] <- "Size_1"
colnames(ModelRich2)[13] <- "Size_2"
colnames(ModelRich2)[14] <- "Size_3"
colnames(ModelRich2)[15] <- "Size_4"

ModelRich2$Age_Scaled <- scale(ModelRich$Crop_Age_Days)
ModelRich2$Day_Scaled <- scale(ModelRich$Day_Sampled)

#remove groups not being modelled for richness

ModelRich2 <- ModelRich2 %>% dplyr::select(-Omnivore, -Hematophagous, -Ambush_Hunting , -Hawking, -Size_4, -Always_Winged, -Polymorphic)

head(ModelRich2);dim(ModelRich2)

colnames(ModelRich2)

rich_names <- colnames(ModelRich2)[2:12]
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
  
  #check which had issues and remove from the mod list
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

#Develops wings null didn't converge properly so investigating that

Richlist1$Develops_Wings

#use optimizer to get a convergence for AIC table
null_develop <- glmmTMB(Develops_Wings ~ 1 + (1 | Field), family = nbinom2, data = ModelRich2,control = glmmTMBControl(optimizer = optim))
summary(null_develop)





##Step 2 - Environmental Variables----

rich_names
Rich_design <- list(
  All = "Age_Scaled",
  Predator = "Day_Scaled + Age_Scaled",
  Herbivore = "Position + Day_Scaled",   
  Fungivore = "Age_Scaled",     
  Web = "Position * Day_Scaled + Age_Scaled",          
  Active_Hunting = "Day_Scaled * Age_Scaled",
  Size_1 = "Position + Age_Scaled + Day_Scaled",
  Size_2 = NULL,      
  Size_3 = "Position + Age_Scaled * Day_Scaled",
  Develops_Wings = "Age_Scaled",
  Wingless = "Position * Age_Scaled + Day_Scaled")
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

###Richness models after first two steps----

Rich_All <- glmmTMB(All ~ Age_Scaled + Field_Area_m2 + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_All)

Rich_Pred <- glmmTMB(Predator ~ Day_Scaled + Age_Scaled + X1km_Prop_Water + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Pred)

Rich_Herb <- glmmTMB(Herbivore ~ Position + Day_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Herb)

Rich_Fung <- glmmTMB(Fungivore ~ Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Fung)

Rich_Web <- glmmTMB(Web ~ Position * Day_Scaled + Age_Scaled + GC + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Web)

Rich_Active <- glmmTMB(Active_Hunting ~ Day_Scaled * Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Active)

Rich_Size1 <- glmmTMB(Size_1 ~ Position + Age_Scaled + Day_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Size1)

Rich_Size2 <- glmmTMB(Size_2 ~ Age_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Size2)

Rich_Size3 <- glmmTMB(Size_3 ~ Position + Age_Scaled * Day_Scaled + X1km_Prop_Water + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Size3)

Rich_DevW <- glmmTMB(Develops_Wings ~ Age_Scaled + NDVIsum_1km + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_DevW)

Rich_Wless <- glmmTMB(Wingless ~ Position * Age_Scaled + Day_Scaled + (1|Field), family = nbinom2, data = ModelRich2)
summary(Rich_Wless)

##Step 3 - Check for spatial autocorrelation----

field_numbers <- unique(ModelRich2$ID)

richmods <- list(Rich_All = Rich_All, Rich_Pred = Rich_Pred, 
                 Rich_Herb = Rich_Herb, Rich_Fung = Rich_Fung,
                 Rich_Web = Rich_Web, Rich_Active = Rich_Active, 
                 Rich_Size1 = Rich_Size1, Rich_Size2 = Rich_Size2,
                 Rich_Size3 = Rich_Size3, Rich_DevW = Rich_DevW, 
                 Rich_Wless = Rich_Wless)
richmods[[3]]

rich_spatial_results <- list()

for (i in names(richmods)) {
  
  cat("=== Processing Model:", i, "===\n") #track which model loop is up to
  
  model_residuals <- simulateResiduals(richmods[[i]])
  spatial_result <- data.frame(
    field = rep(NA, length(field_numbers)),
    statistic = rep(NA, length(field_numbers)),
    p_value = rep(NA, length(field_numbers)),
    method = rep(NA_character_, length(field_numbers)),
    stringsAsFactors = FALSE)
  
  s <- 1
  
  for (f in field_numbers) {
    
    cat("Field", f, "\n") #What field is it doing?
    
    #Extracting specific residuals for individual fields
    field_indices <- which(ModelRich2$ID == f)
    field_residuals <- model_residuals
    field_residuals$scaledResiduals <- 
      model_residuals$scaledResiduals[field_indices]
    field_residuals$fittedPredictedResponse <- 
      model_residuals$fittedPredictedResponse[field_indices]
    
    # Test spatial autocorrelation using your grid coordinates
    spatial_test <- testSpatialAutocorrelation(field_residuals, 
               x = ModelRich2$X_Cor[ModelRich2$ID == f], 
               y = ModelRich2$Y_Cor[ModelRich2$ID == f])
    
    
    spatial_result$field [s] <- f
    spatial_result$statistic [s] <- spatial_test$statistic[1] 
    spatial_result$p_value [s] <- spatial_test$p.value
    spatial_result$method [s] <- spatial_test$method
    
    s <- s + 1
    
  }
  
  rich_spatial_results[[i]] <- spatial_result
  
}

#My own notes/explanation of different bits of results
spatial_test$statistic # The test statistic - note the important one here is observed - want this to be close to zero 
spatial_test$p.value # P-value for the test - indicates if there is significant spatial autocorrelation - statistic above indicates the magnitude and the direction (positive or negative)
spatial_test$method #Method used to get statistic


names(rich_spatial_results)
rich_spatial_results$Rich_Wless

#Moran's I
min(rich_spatial_results$Rich_Wless[2])
max(rich_spatial_results$Rich_Wless[2])

#p-value
min(rich_spatial_results$Rich_Wless[3])
max(rich_spatial_results$Rich_Wless[3])

##Step 4 - Predictions----

#variables included in the top models
Predictions_Day <- seq(min(ModelRich2$Day_Scaled),max(ModelRich2$Day_Scaled),length.out=20) 
Predictions_Age <- seq(min(ModelRich2$Age_Scaled),max(ModelRich2$Age_Scaled),length.out=20)
Predictions_Position <- as.factor(c("Outer","Inner"))
Predictions_Area <- seq(min(ModelRich2$Field_Area_m2),max(ModelRich2$Field_Area_m2),length.out=20)
Predictions_Water <- seq(min(ModelRich2$X1km_Prop_Water),max(ModelRich2$X1km_Prop_Water),length.out=20)
Predictions_GC <- seq(min(ModelRich2$GC),max(ModelRich2$GC),length.out=20)
Predictions_1kmNDVI <- seq(min(ModelRich2$NDVIsum_1km),max(ModelRich2$NDVIsum_1km),length.out=20)


###Creating new data for predictions----
colnames(ModelRich2)

richpred_all <- expand.grid(Age_Scaled = Predictions_Age, Field_Area_m2 = Predictions_Area)
head(richpred_all);dim(richpred_all)

richpred_preds <- expand.grid(Day_Scaled = Predictions_Day, Age_Scaled = Predictions_Age, X1km_Prop_Water = Predictions_Water)
head(richpred_preds);dim(richpred_preds)

richpred_herb <- expand.grid(Position = Predictions_Position, Day_Scaled = Predictions_Day)
head(richpred_herb);dim(richpred_herb)

richpred_fung <- data.frame(Age_Scaled = Predictions_Age)
head(richpred_fung);dim(richpred_fung)

richpred_web <- expand.grid(Position = Predictions_Position, Day_Scaled = Predictions_Day, Age_Scaled = Predictions_Age, GC = Predictions_GC)
head(richpred_web);dim(richpred_web)

richpred_active <- expand.grid(Day_Scaled = Predictions_Day, Age_Scaled = Predictions_Age)
head(richpred_active);dim(richpred_active)

richpred_size1 <- expand.grid(Position = Predictions_Position, Age_Scaled = Predictions_Age, Day_Scaled = Predictions_Day)
head(richpred_size1);dim(richpred_size1)

richpred_size2 <- data.frame(Age_Scaled = Predictions_Age)
head(richpred_size2);dim(richpred_size2)

richpred_size3 <- expand.grid(Position = Predictions_Position, Age_Scaled = Predictions_Age, Day_Scaled = Predictions_Day, X1km_Prop_Water = Predictions_Water)
head(richpred_size3);dim(richpred_size3)

richpred_DevW <- expand.grid(Age_Scaled = Predictions_Age, NDVIsum_1km = Predictions_1kmNDVI)
head(richpred_DevW);dim(richpred_DevW)

richpred_Wless <- expand.grid(Position = Predictions_Position, Age_Scaled = Predictions_Age, Day_Scaled = Predictions_Day)
head(richpred_Wless);dim(richpred_Wless)

richpredlist <- list(richpred_all = richpred_all, 
                     richpred_preds = richpred_preds, 
                     richpred_herb = richpred_herb,
                    richpred_fung = richpred_fung, 
                    richpred_web = richpred_web, 
                    richpred_active = richpred_active,
                    richpred_size1 = richpred_size1, 
                    richpred_size2 = richpred_size2, 
                    richpred_size3 = richpred_size3,
                    richpred_DevW = richpred_DevW, 
                    richpred_Wless = richpred_Wless)

richpredlist[["richpred_size3"]]

names(richmods)[1]

###Loop predictions----

predict(richmods[[1]],newdata=richpredlist[[1]],se.fit = T, type = "link",re.form = NA)

richpredresults <- list()
richpredlist[[1]]

for (i in 1:11) {
  
  d <- names(richpredlist)[i]
  
  cat("=== Processing Predictions:", d, "===\n") #track which prediction the loop is up to
  
  prediction1 <- predict(object = richmods[[i]],newdata= richpredlist[[i]],se.fit = T, type = "link",re.form = NA)
  
  prediction2<-data.frame(richpredlist[[i]],fit.link=prediction1$fit,se.link=prediction1$se.fit)
  
  prediction2$lci.link<-prediction2$fit.link-(1.96*prediction2$se.link)
  prediction2$uci.link<-prediction2$fit.link+(1.96*prediction2$se.link)
  
  prediction2$fit<-exp(prediction2$fit.link)
  prediction2$se<-exp(prediction2$se.link)
  prediction2$lci<-exp(prediction2$lci.link)
  prediction2$uci<-exp(prediction2$uci.link)
  
  richpredresults[[d]] <- prediction2
  
}

head(richpredresults[[1]])


#DIVERSITY MODELLING----

##Step 1 - Design Variables----

head(ModelDiv2);dim(ModelDiv2)
str(ModelDiv2)
colnames(ModelDiv2)

colnames(ModelDiv2)[12] <- "Size_1"
colnames(ModelDiv2)[13] <- "Size_2"
colnames(ModelDiv2)[14] <- "Size_3"
colnames(ModelDiv2)[15] <- "Size_4"

ModelDiv2$Age_Scaled <- scale(ModelDiv2$Crop_Age_Days)
ModelDiv2$Day_Scaled <- scale(ModelDiv2$Day_Sampled)

#remove groups not modelling diversity for

ModelDiv2 <- ModelDiv2 %>% dplyr::select(-Omnivore,-Hematophagous,-Ambush_Hunter,-Hawking,-Size_4,-Always_Winged,-Polymorphic)

head(ModelDiv2);dim(ModelDiv2)

colnames(ModelDiv2)

div_names <- colnames(ModelDiv2)[2:12]
length(div_names)
divlist1 <- list()

#LOOP GOES HERE - WHAT MODEL AM I USING?
#TO DO----

#OCCURRENCE MODELLING----

##Step 1 - Design Variables----

head(ModelOccur2);dim(ModelOccur2)
str(ModelOccur2)
colnames(ModelOccur2)

colnames(ModelOccur2)[10] <- "Size_2"
colnames(ModelOccur2)[11] <- "Size_3"
colnames(ModelOccur2)[12] <- "Size_4"

ModelOccur2$Age_Scaled <- scale(ModelOccur2$Crop_Age_Days)
ModelOccur2$Day_Scaled <- scale(ModelOccur2$Day_Sampled)

#remove develops wings

ModelOccur2 <- ModelOccur2 %>% dplyr::select(-Develops_Wings)

head(ModelOccur2);dim(ModelOccur2)

colnames(ModelOccur2)

occur_names <- colnames(ModelOccur2)[2:15]
length(occur_names)
occurlist1 <- list()

for (i in occur_names) {
  allsum <- list()
  
  cat("=== Processing group:", i, "===\n") #track where loop is up to
  
  #run models and collect summaries
  null <- glmer(as.formula(paste(i, "~ 1 + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum [[1]] <- summary(null)
  P <- glmer(as.formula(paste(i, "~ Position + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum [[2]] <- summary(P)
  A <- glmer(as.formula(paste(i, "~ Age_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum [[3]] <- summary(A)
  D <- glmer(as.formula(paste(i, "~ Day_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[4]] <- summary(D)
  PA <- glmer(as.formula(paste(i, "~ Position + Age_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[5]] <- summary(PA)
  PD <- glmer(as.formula(paste(i, "~ Position + Day_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[6]] <- summary(PD)
  DA <- glmer(as.formula(paste(i, "~ Day_Scaled + Age_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[7]] <- summary(DA)
  PxA <- glmer(as.formula(paste(i, "~ Position * Age_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[8]] <- summary(PxA)
  PxD <- glmer(as.formula(paste(i, "~ Position * Day_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[9]] <- summary(PxD)
  DxA <- glmer(as.formula(paste(i, "~ Day_Scaled * Age_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[10]] <- summary(DxA)
  PAD <- glmer(as.formula(paste(i, "~ Position + Age_Scaled + Day_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[11]] <- summary(PAD)
  PxAD <- glmer(as.formula(paste(i, "~ Position * Age_Scaled + Day_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[12]] <- summary(PxAD)
  PxDA <- glmer(as.formula(paste(i, "~ Position * Day_Scaled + Age_Scaled + (1 | Field)")), family = binomial, data = ModelOccur2)
  allsum[[13]] <- summary(PxDA)
  PAxD <- glmer(as.formula(paste(i, "~ Position + Age_Scaled * Day_Scaled + (1 | Field)")), data = ModelOccur2,family = binomial,)
  allsum[[14]] <- summary(PAxD)
  
  #collect models
  tempmodlist <- list("null" = null, "P" = P, "A" = A, "D" = D,
                      "PA" = PA, "PD" = PD, "DA" = DA, 
                      "PxA" = PxA, "PxD" = PxD, "DxA" = DxA, 
                      "PAD" = PAD, "PxAD" = PxAD, "PxDA" = PxDA, 
                      "PAxD" = PAxD)
  
  #check which had issues and remove from the mod list
  has_issues <- sapply(tempmodlist, function(model) {
    if (is.null(model)) return(TRUE) #model failed
    if (model@optinfo$conv$opt != 0) return(TRUE) #convergence issues
    if (is.na(AIC(model))) return(TRUE) #does it have an AIC (not having one will cause all AIC in table to be NA)
    return(FALSE)
  }) 
  
  cleaned_models <- tempmodlist[!has_issues]
  cat("Models with issues:", sum(has_issues), "\n")
  cat("Problem models:", names(tempmodlist)[has_issues], "\n")
  
  occurlist1 [[i]] <- aictab(cleaned_models,modnames = NULL)
  
}

length(occurlist1)

occurlist1

##Step 2 - Environmental Variables----

occur_names
Occur_design <- list(
  Herbivore = "Day_Scaled",
  Omnivore = NULL,
  Fungivore = "Day_Scaled + Age_Scaled",   
  Hematophagous = 'Day_Scaled * Age_Scaled',
  Web = "Day_Scaled + Age_Scaled",          
  Active_Hunting = "Day_Scaled",
  Ambush_Hunter = NULL,          
  Hawking = "Position + Age_Scaled * Day_Scaled",
  Size_2 = NULL,      
  Size_3 = "Day_Scaled * Age_Scaled",
  Size_4 = "Age_Scaled",
  Always_Winged = NULL,
  Wingless = "Position * Age_Scaled + Day_Scaled",
  Polymorphic = "Position * Day_Scaled")
head(Occur_design)

occurlist2 <- list()


for (i in occur_names) {
  allsum <- list()
  
  design_formula_part <- Occur_design[[i]]
  
  cat("=== Processing group:", i, "===\n") #track where loop is up to
  cat("Design formula part:", design_formula_part, "\n") #making sure it's taken the right formula
  
  #Handle groups without null model as top in step 1
  if (is.null(design_formula_part) || design_formula_part == "") {
    # No design variables - intercept only model
    base_formula <- paste(i, "~ 1 + (1 | Field)")} else {
      # Has design variables
      base_formula <- paste(i, "~", design_formula_part, "+ (1 | Field)")}
  base_model <- glmer(as.formula(base_formula),family = binomial, data = ModelOccur2)
  allsum[[1]] <- summary(base_model)
  
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    height_formula <- paste(i, "~ Height + (1 | Field)")} else {
      height_formula <- paste(i, "~", design_formula_part, "+ Height + (1 | Field)")}
  height_model <- glmer(as.formula(height_formula), family = binomial, data = ModelOccur2)
  allsum[[2]] <- summary(height_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    gc_formula <- paste(i, "~ GC  + (1 | Field)")} else {
      gc_formula <- paste(i, "~", design_formula_part, "+ GC + (1 | Field)")}
  gc_model <- glmer(as.formula(gc_formula), family = binomial, data = ModelOccur2)
  allsum[[3]] <- summary(gc_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    Fsize_formula <- paste(i, "~ Field_Area_m2  + (1 | Field)")} else {
      Fsize_formula <- paste(i, "~", design_formula_part, "+ Field_Area_m2 + (1 | Field)")}
  Fsize_model <- glmer(as.formula(Fsize_formula), family = binomial, data = ModelOccur2)
  allsum[[4]] <- summary(Fsize_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    water_formula <- paste(i, "~ X1km_Prop_Water  + (1 | Field)")} else {water_formula <- paste(i, "~", design_formula_part, "+ X1km_Prop_Water + (1 | Field)")}
  water_model <- glmer(as.formula(water_formula), family = binomial, data = ModelOccur2)
  allsum[[5]] <- summary(water_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    NDVIf_formula <- paste(i, "~ NDVImean_Field  + (1 | Field)")} else {NDVIf_formula <- paste(i, "~", design_formula_part, "+ NDVImean_Field + (1 | Field)")}
  NDVIf_model <- glmer(as.formula(NDVIf_formula), family = binomial, data = ModelOccur2)
  allsum[[6]] <- summary(NDVIf_model)
  
  if (is.null(design_formula_part) || design_formula_part == "") {
    NDVI1km_formula <- paste(i, "~ NDVIsum_1km  + (1 | Field)")} else {
      NDVI1km_formula <- paste(i, "~", design_formula_part, "+ NDVIsum_1km + (1 | Field)")}
  NDVI1km_model <- glmer(as.formula(NDVI1km_formula), family = binomial, data = ModelOccur2)
  allsum[[7]] <- summary(NDVI1km_model)
  
  
  #collect models
  tempmodlist <- list("base" = base_model, "height" = height_model, 
                      "GC" = gc_model, "Feild Size" = Fsize_model,
                      "water" = water_model,"Field NDVI" = NDVIf_model,
                      "NDVI 1km" = NDVI1km_model)
  
  has_issues <- sapply(tempmodlist, function(model) {
    if (is.null(model)) return(TRUE) #model failed
    if (model@optinfo$conv$opt != 0) return(TRUE) #convergence issues
    if (is.na(AIC(model))) return(TRUE) #does it have an AIC (not having one will cause all AIC in table to be NA)
    return(FALSE)
  }) 
  
  cleaned_models <- tempmodlist[!has_issues]
  cat("Models with issues:", sum(has_issues), "\n")
  cat("Problem models:", names(tempmodlist)[has_issues], "\n")
  
  occurlist2 [[i]] <- aictab(cleaned_models,modnames = NULL)
  
}

length(occurlist2)

occurlist2

###Occurrence models after first two steps----
occur_names
Age_Scaled

Occur_herb <- glmer(Herbivore ~ Day_Scaled + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_herb)

Occur_omni <- glmer(Omnivore ~ GC + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_omni)

ModelOccur2$FieldSize_Scaled <- scale(ModelOccur2$Field_Area_m2)
#Needed to scale field size for the model to be more identifiable
Occur_fung <- glmer(Fungivore ~ Day_Scaled + Age_Scaled + FieldSize_Scaled + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_fung)

Occur_hema <- glmer(Hematophagous ~ Day_Scaled * Age_Scaled + NDVImean_Field + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_hema)

Occur_web <- glmer(Web ~ Day_Scaled * Age_Scaled + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_web)

Occur_active <- glmer(Active_Hunting ~ Day_Scaled + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_active)

Occur_hawk <- glmer(Hawking ~ Position + Age_Scaled * Day_Scaled + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_hawk)

Occur_size2 <- glmer(Size_2 ~ Height + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_size2)

Occur_size3 <- glmer(Size_3 ~ Day_Scaled * Age_Scaled + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_size3)

Occur_size4 <- glmer(Size_4 ~ Age_Scaled + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_size4)

Occur_AlWing <- glmer(Always_Winged ~ X1km_Prop_Water + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_AlWing)

Occur_Wless <- glmer(Wingless ~ Position * Age_Scaled + Day_Scaled + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_Wless)

Occur_poly <- glmer(Polymorphic ~ Position * Day_Scaled + (1|Field), family = binomial, data = ModelOccur2)
summary(Occur_poly)

##Step 3 - Check for spatial autocorrelation----

field_numbers3 <- unique(ModelOccur2$ID)

occurmods <- list(Occur_herb = Occur_herb, Occur_omni = Occur_omni, 
                 Occur_fung = Occur_fung, Occur_hema = Occur_hema,
                 Occur_web = Occur_web, Occur_active = Occur_active, 
                 Occur_hawk = Occur_hawk, Occur_size2 = Occur_size2,
                 Occur_size3 = Occur_size3, Occur_size4 = Occur_size4, 
                Occur_AlWing = Occur_AlWing, Occur_Wless = Occur_Wless,
                Occur_poly = Occur_poly)
occurmods[[3]]

occur_spatial_results <- list()

for (i in names(occurmods)) {
  
  cat("=== Processing Model:", i, "===\n") #track which model loop is up to
  
  model_residuals <- simulateResiduals(occurmods[[i]])
  spatial_result <- data.frame(
    field = rep(NA, length(field_numbers3)),
    statistic = rep(NA, length(field_numbers3)),
    p_value = rep(NA, length(field_numbers3)),
    method = rep(NA_character_, length(field_numbers3)),
    stringsAsFactors = FALSE)
  
  s <- 1
  
  for (f in field_numbers3) {
    
    cat("Field", f, "\n") #What field is it doing?
    
    #Extracting specific residuals for individual fields
    field_indices <- which(ModelOccur2$ID == f)
    field_residuals <- model_residuals
    field_residuals$scaledResiduals <- 
      model_residuals$scaledResiduals[field_indices]
    field_residuals$fittedPredictedResponse <- 
      model_residuals$fittedPredictedResponse[field_indices]
    
    # Test spatial autocorrelation using your grid coordinates
    spatial_test <- testSpatialAutocorrelation(field_residuals, 
                                               x = ModelOccur2$X_Cor[ModelOccur2$ID == f], 
                                               y = ModelOccur2$Y_Cor[ModelOccur2$ID == f])
    
    
    spatial_result$field [s] <- f
    spatial_result$statistic [s] <- spatial_test$statistic[1] 
    spatial_result$p_value [s] <- spatial_test$p.value
    spatial_result$method [s] <- spatial_test$method
    
    s <- s + 1
    
  }
  
  occur_spatial_results[[i]] <- spatial_result
  
}


names(occur_spatial_results)
occur_spatial_results$Occur_poly

#Moran's I
min(occur_spatial_results$Occur_poly[2])
max(occur_spatial_results$Occur_poly[2])

#p-value
min(occur_spatial_results$Occur_poly[3])
max(occur_spatial_results$Occur_poly[3])


##Step 4 - Predictions----

#variables included in the top models (not already done earlier in process)

Predictions_Height <- seq(min(ModelOccur2$Height),max(ModelOccur2$Height),length.out = 20)
Predictions_Areascale <- seq(min(ModelOccur2$FieldSize_Scaled),max(ModelOccur2$FieldSize_Scaled),length.out = 20)
Predictions_NDVIfield <- seq(min(ModelOccur2$NDVImean_Field ),max(ModelOccur2$NDVImean_Field ),length.out = 20)


###Creating new data for predictions----
colnames(ModelOccur2)

occurpred_herb <- data.frame(Day_Scaled = Predictions_Day)
head(occurpred_herb);dim(occurpred_herb)

occurpred_omni <- data.frame(GC = Predictions_GC)
head(occurpred_omni);dim(occurpred_omni)

occurpred_fung <- expand.grid(Day_Scaled = Predictions_Day, Age_Scaled = Predictions_Age, FieldSize_Scaled  = Predictions_Areascale)
head(occurpred_fung);dim(occurpred_fung)

occurpred_hema <- expand.grid(Day_Scaled = Predictions_Day, Age_Scaled = Predictions_Age, NDVImean_Field  = Predictions_NDVIfield)
head(occurpred_hema);dim(occurpred_hema)

occurpred_web <- expand.grid(Day_Scaled = Predictions_Day, Age_Scaled = Predictions_Age)
head(occurpred_web);dim(occurpred_web)

occurpred_active <- data.frame(Day_Scaled = Predictions_Day)
head(occurpred_active);dim(occurpred_active)

occurpred_hawk <- expand.grid(Position = Predictions_Position, Age_Scaled = Predictions_Age, Age_Scaled = Predictions_Age,Day_Scaled = Predictions_Day)
head(occurpred_hawk);dim(occurpred_hawk)

occurpred_size2 <- data.frame(Height = Predictions_Height)
head(occurpred_size2);dim(occurpred_size2)

occurpred_size3 <- expand.grid(Day_Scaled = Predictions_Day, Age_Scaled = Predictions_Age)
head(occurpred_size3);dim(occurpred_size3)

occurpred_size4 <- data.frame(Age_Scaled = Predictions_Age)
head(occurpred_size4);dim(occurpred_size4)

occurpred_AlWing <- data.frame(X1km_Prop_Water = Predictions_Water)
head(occurpred_AlWing);dim(occurpred_AlWing)

occurpred_Wless <- expand.grid(Position = Predictions_Position, Age_Scaled = Predictions_Age, Day_Scaled = Predictions_Day)
head(occurpred_Wless);dim(occurpred_Wless)

occurpred_poly <- expand.grid(Position = Predictions_Position, Day_Scaled = Predictions_Day)
head(occurpred_poly);dim(occurpred_poly)

occurpredlist <- list(occurpred_herb = occurpred_herb, 
                      occurpred_omni = occurpred_omni, 
                     occurpred_fung = occurpred_fung,
                     occurpred_hema = occurpred_hema, 
                     occurpred_web = occurpred_web, 
                     occurpred_active = occurpred_active,
                     occurpred_hawk = occurpred_hawk, 
                     occurpred_size2 = occurpred_size2, 
                     occurpred_size3 = occurpred_size3,
                     occurpred_size4 = occurpred_size4,
                     occurpred_AlWing = occurpred_AlWing,
                     occurpred_Wless= occurpred_Wless,
                     occurpred_poly = occurpred_poly)

occurpredlist[["occurpred_AlWing"]]

names(occurmods)[1]


###Loop predictions----

predict(occurmods[[1]],newdata=occurpredlist[[1]],se.fit = T, type = "link",re.form = NA)

occurpredresults <- list()
occurpredlist[[1]]
length(occurpredlist)

for (i in 1:13) {
  
  d <- names(occurpredlist)[i]
  
  cat("=== Processing Predictions:", d, "===\n") #track which prediction the loop is up to
  
  prediction1 <- predict(object = occurmods[[i]],newdata= occurpredlist[[i]],se.fit = T, type = "link",re.form = NA)
  
  prediction2<-data.frame(occurpredlist[[i]],fit.link=prediction1$fit,se.link=prediction1$se.fit)
  
  prediction2$lci.link<-prediction2$fit.link-(1.96*prediction2$se.link)
  prediction2$uci.link<-prediction2$fit.link+(1.96*prediction2$se.link)
  
  prediction2$fit<-exp(prediction2$fit.link)
  prediction2$se<-exp(prediction2$se.link)
  prediction2$lci<-exp(prediction2$lci.link)
  prediction2$uci<-exp(prediction2$uci.link)
  
  occurpredresults[[d]] <- prediction2
  
}

head(occurpredresults[[1]])



#END----


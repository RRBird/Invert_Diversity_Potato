
#Code to run for checking spatial autocorrelation on models----

#includes checking all and checking each field individually - use on top models only to check 

library(DHARMa)

# After fitting your GLM/GLMM
# Replace 'your_model' with your actual model object
model_residuals <- simulateResiduals(your_model)

# Test spatial autocorrelation using your grid coordinates
# Replace 'your_data' with your actual dataframe name
spatial_test <- testSpatialAutocorrelation(model_residuals, 
                                           x = your_data$X, 
                                           y = your_data$Y)

# Print the results
print(spatial_test)

#example results
$statistic
[1] 0.234  # This is the Moran's I value

$p.value  
[1] 0.043  # This is your p-value 

$alternative
[1] "greater"  # Tests if autocorrelation > 0

#How to Interpret Moran's I Values
##Moran's I = 0: No spatial autocorrelation (random pattern)
##Moran's I > 0: Positive autocorrelation (similar values cluster together)
##Moran's I < 0: Negative autocorrelation (dissimilar values cluster together - rare in ecology)

#p > 0.05: No significant spatial autocorrelation - your model is fine
#p < 0.05:️ Significant spatial autocorrelation - your model residuals show spatial patterns

# Optional: Plot the residuals spatially to visualize patterns
plot(your_data$X, your_data$Y, 
     col = ifelse(residuals(model_residuals) > 0, "red", "blue"),
     pch = 16, cex = abs(residuals(model_residuals)) + 0.5,
     main = "Spatial pattern of residuals",
     xlab = "X grid coordinate", ylab = "Y grid coordinate")

# Alternative: More detailed DHARMa diagnostic plot
plot(model_residuals)  # This includes multiple diagnostic plots

# If you want to test field-by-field (optional)
# This loops through each field separately
library(dplyr)
for(field in unique(your_data$Field)) {
  field_data <- filter(your_data, Field == field)
  field_residuals <- residuals(model_residuals)[your_data$Field == field]
  
  if(length(unique(field_data$X)) > 1 & length(unique(field_data$Y)) > 1) {
    cat("Testing field:", field, "\n")
    field_test <- testSpatialAutocorrelation(model_residuals[your_data$Field == field], 
                                             x = field_data$X, 
                                             y = field_data$Y)
    print(field_test)
  }
}

#First attempt at loop richness models step 1 ----

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

#Old - Model Trial ----

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


#Old Diversity calculation (diversity by orders instead of functional groups)----


diversity(x = invert[which(invert$Order=='Araneae'),] )
dim(invert[which(invert$Order=='Araneae'),])

###Araneae

Araneae <- invert[which(invert$Order=='Araneae'),]
head(Araneae);dim(Araneae)
Araneae_matrix <- table(Araneae$Site, Araneae$Morphospecies)
head(Araneae_matrix);dim(Araneae_matrix)

Araneae_diversity<-diversity(Araneae_matrix, index = "invsimpson")
Araneae_diversity<-data.frame(Site = names(Araneae_diversity), Diversity_A = as.numeric(Araneae_diversity))
head(Araneae_diversity);dim(Araneae_diversity)
variables <- merge(variables,Araneae_diversity, by = "Site",all.x = T)
head(variables);dim(variables)
variables$Diversity_A[is.na(variables$Diversity_A)] <- 0.000001



#Code I used for finding the speceis rich of all functional groups
#just change out the column names and the level
length(unique(invert$Morphospecies[!is.na(invert$Wings) & invert$Wings == 'Polymorphic']))
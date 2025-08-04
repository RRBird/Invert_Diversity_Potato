
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
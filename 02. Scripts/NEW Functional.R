#From Tutorial: Methods for assessing functional responses to environmental gradients -- Kleyer et al. -- June 22, 2009

#Libraries----

library("ade4")
library("dplyr")
library("tidyverse")
library("vegan")
library('fpc')
library("tidyr")
library("AICcmodavg")
library('lme4')
library("glmmTMB")
library("DHARMa")

Q_matrix$Size <- factor(Q_matrix$Size, 
                        levels = c("0-2.5", "2.5-5", "5-10", 
                                   ">10","No_Size", "Unknown"),
                        ordered = FALSE) #to make it work properly
str(Q_matrix)

R_matrix$Position <- factor(R_matrix$Position, 
                            levels = c("Inner", "Outer"),
                            ordered = FALSE) #to make it work properly
str(R_matrix)

#Species responses to environmental gradients----

coa1 <- dudi.coa(L_matrix, scannf = F)
cca1 <- pcaiv(coa1, R_matrix, scannf = F)

#percentage of variation in species composition explained by enviro
100 * sum(cca1$eig) / sum(coa1$eig) 

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
s.label(cca1$c1, clabel = 0)
par(mar = c(0.1, 0.1, 0.1, 0.1))
s.arrow(cca1$cor[-1,], add.plot=TRUE)

#RLQ Analysis---- 

pca.traits <- dudi.hillsmith(Q_matrix, row.w = coa1$cw, scannf = FALSE)
pca.env <- dudi.hillsmith(R_matrix, row.w = coa1$lw, scannf = FALSE)

rlq1 <- rlq(pca.env, coa1, pca.traits, scannf = FALSE)
summary(rlq1)

dev.new(height=10,width=15,dpi=80,pointsize=14,noRStudioGD = T)
plot(rlq1)

## Percentage of co-Inertia for each axis
100*rlq1$eig/sum(rlq1$eig)

#To interpret the results, correlations can be computed:
## weighted correlations axes / env.
t(pca.env$tab)%*%(diag(pca.env$lw))%*%as.matrix(rlq1$mR)

## weighted correlations axes / traits.
t(pca.traits$tab)%*%(diag(pca.traits$lw))%*%as.matrix(rlq1$mQ)

## correlations traits / env.
rlq1$tab

#Biplot (traits and enviro) ----
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
s.arrow(rlq1$c1, xlim=c(-1,1), boxes = FALSE)
s.label(rlq1$li, add.plot=T, clab=1.5)

#Another plot
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
s.label(rlq1$lQ, clabel = 0)
par(mar = c(0.1, 0.1, 0.1, 0.1))
pointLabel(rlq1$lQ,row.names(rlq1$lQ), cex=0.7) #pointLabel didn't work so need an alternative to add names to this one

#Classifying scores to obtain functional groups----
hc2 <- hclust(dist(rlq1$lQ), method = "ward.D")
dev.new(height=20,width=40,dpi=80,pointsize=14,noRStudioGD = T)
plot(hc2)

#Calinsky-Harabasz criteria to find best partition 
ntest <- 6
res <- rep(0,ntest - 1)

for (i in 2:ntest){
  fac <- cutree(hc2, k = i)
  res[i-1] <- as.numeric(calinhara(rlq1$lQ, fac))
}
#calinski didn't work so needed to change it to calinhara 

#Trait group numbers----
dev.new(height=5,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mfrow=c(1,2))
plot(2:ntest, res, type='b', pch=20, xlab="Number of groups", ylab = "C-H index")
plot(3:ntest, diff(res), type='b', pch=20, xlab="Number of groups", ylab = "Diff in C-H index")

spe.group2 <- as.factor(cutree(hc2, k = which.max(res) +1))
levels(spe.group2) <- c('E',"C","B","D","A")
spe.group2 <- factor(spe.group2, levels=c("A","B","C","D","E"))

#biplot (trait groups + Traits and also with trait groups + enviro)----
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
ade4::s.class(rlq1$lQ, spe.group2, col= 1:nlevels(spe.group2))
s.arrow(rlq1$c1, add.plot = T,clab=0.8)

#here it is without the traits laid over the top
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
ade4::s.class(rlq1$lQ, spe.group2, col= 1:nlevels(spe.group2))

#here it is with the enviro groups
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
ade4::s.class(rlq1$lQ, spe.group2, col= 1:nlevels(spe.group2))
s.arrow(rlq1$l1, add.plot = T, clab = 0.6,boxes = FALSE)

l1_short <- rlq1$l1
rownames(l1_short) <- c("     H", "GC", "In", "Out","Day", "A           ", "FA", "C", "NF", "N1km","R")
l1_labels <- l1_short * 1.30

dev.new(height=10, width=10, dpi=80, pointsize=14, noRStudioGD = T)
ade4::s.class(rlq1$lQ, spe.group2, col = 1:nlevels(spe.group2))
s.arrow(rlq1$l1, add.plot = T, clab = 0)
s.label(l1_labels, add.plot = T, clab = 0.7, boxes = T)

# Adjust these to match your actual environmental variables in order
# Check what they are first:
rownames(rlq1$l1)



#Trait loadings on Axis 1
print(sort(rlq1$c1[,1], decreasing = TRUE))

#Environmental loadings on Axis 1
print(sort(rlq1$l1[,1], decreasing = TRUE))

#group centroids
group_centers <- aggregate(rlq1$lQ, by = list(spe.group2), FUN = mean)


#We can interpret this partition in terms of traits----
eta2 <- ade4::cor.ratio(Q_matrix[,-1], data.frame(spe.group2), weights = rep(1, length(spe.group2)))

# Calculate eta-squared for each trait
eta2 <- sapply(Q_matrix, function(trait) {
  # Total sum of squares
  grand_mean <- mean(as.numeric(trait))
  TSS <- sum((as.numeric(trait) - grand_mean)^2)
  
  # Between-group sum of squares
  group_means <- tapply(as.numeric(trait), spe.group2, mean)
  group_sizes <- table(spe.group2)
  BSS <- sum(group_sizes * (group_means - grand_mean)^2)
  
  # Eta-squared (proportion of variance explained by groups)
  eta_sq <- BSS / TSS
  return(eta_sq)
})

# View the results
print(eta2)

# Sort to see which traits are most associated with the groups
sort(eta2, decreasing = TRUE)


#Visually see the groups----
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mfrow=n2mfrow(ncol(Q_matrix)))

# Plot first trait
plot(table(spe.group2, Q_matrix[,1]), main = names(Q_matrix)[1])

# Plot remaining traits with eta2 values
for(i in 2:ncol(Q_matrix)){
  label <- paste(names(Q_matrix)[i], "(cor.ratio =", round(eta2[i], 3), ")")  # Changed from eta2[i-1] to eta2[i]
  plot(Q_matrix[,i] ~ spe.group2, main = label, border = 1:nlevels(spe.group2))
}


#get a better idea of traits associated with different groups---

#Create a summary table showing trait composition of each group
for(trait in names(Q_matrix)) {
  cat("\n", trait, ":\n")
  print(table(spe.group2, Q_matrix[, trait]))
  cat("\n")
}

#as proportions (easier to compare)
for(trait in names(Q_matrix)) {
  cat("\n", trait, " (proportions):\n")
  print(round(prop.table(table(spe.group2, Q_matrix[, trait]), margin = 1), 2))
  cat("\n")
}


# Function to find modal (most common) value for each group
group_profiles <- data.frame(Group = levels(spe.group2))

for(trait in names(Q_matrix)) {
  modal_values <- tapply(Q_matrix[, trait], spe.group2, function(x) {
    names(sort(table(x), decreasing = TRUE))[1]  # Most common value
  })
  group_profiles[, trait] <- modal_values
}

print(group_profiles)
#Interesting but I think the proportions table (above) and the heat map below give a better idea of the groups

# Create a heatmap of group-trait associations----

# Calculate proportions for each trait-group combination
heatmap_data <- data.frame()
for(trait in names(Q_matrix)) {
  prop_table <- prop.table(table(spe.group2, Q_matrix[, trait]), margin = 1)
  temp_df <- as.data.frame(prop_table)
  names(temp_df) <- c("Group", "Trait_Value", "Proportion")
  temp_df$Trait <- trait
  heatmap_data <- rbind(heatmap_data, temp_df)
}

str(heatmap_data)

heatmap_data$Trait_Value <- gsub("Active_Hunting", "Active", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("Ambush_Hunter", "Ambush", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("0-2.5", "0-2.5mm", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("2.5-5", "2.5-5mm", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("5-10", "5-10mm", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub(">10", ">10mm", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("No_Size", "No Size", heatmap_data$Trait_Value)

heatmap_data$Trait <- gsub("Hunting.Style", "Hunting", heatmap_data$Trait)

colnames(heatmap_data)[2] <- "Traits"
colnames(heatmap_data)
str(heatmap_data)
colnames(heatmap_data)[4] <- "Functional_Group"

#trying to get a better order for each functional group traits

heatmap_data <- heatmap_data %>%
  group_by(Functional_Group) %>%
  mutate(Traits = if(unique(Functional_Group) == "Size") {
    factor(Traits, levels = c(
      "Unknown","No Size","0-2.5mm","2.5-5mm","5-10mm", ">10mm"))
  } else {
    factor(Traits, levels = sort(unique(Traits)))
  }) %>%
  ungroup()
#worked for size and trophic but not the other two

dev.new(height=10,width=15,dpi=80,pointsize=14,noRStudioGD = T)
ggplot(heatmap_data, aes(x = Group, y = Traits, fill = Proportion)) +
  geom_tile() +
  facet_wrap(~Functional_Group, scales = "free_y") +
  scale_fill_gradient(low = "white", high = "darkblue") + 
  theme_minimal(base_size = 16)


#Extract the trait groups----

Trait_Group <- data.frame(Morphospecies = names(spe.group2), trait_group = spe.group2)
head(Trait_Group);dim(Trait_Group)


invert_trait_group <- merge(Trait_Group,invert_filtered, by = "Morphospecies")
head(invert_trait_group);dim(invert_trait_group)
colnames(invert_trait_group)[colnames(invert_trait_group) == "ID"] <- "Site"

#Creating Modelling data and calulating functuional Diversity 

FDModel <- data.frame(Site = variables$Site)

#functional richness
FUNrichness <- aggregate(trait_group ~ Site, data = invert_trait_group, FUN = function(x) length(unique(x)))
FDModel <- merge(FDModel,FUNrichness,by = "Site",all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[2] <- "Fun_Rich"
head(FDModel);dim(FDModel)

FDModel$Fun_Rich[is.na(FDModel$Fun_Rich)] <- 0

#Functional Diversity
FUNdiversity <- aggregate(trait_group ~ Site, data = invert_trait_group, FUN = function(x) diversity(table(x), index = "invsimpson"))
FDModel <- merge(FDModel,FUNdiversity,by = "Site",all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[3] <- "Fun_Div"
head(FDModel);dim(FDModel)

FDModel$Fun_Div[is.na(FDModel$Fun_Div)] <- 0.00001

#Variables and XY coordinates

FDModel <- merge(FDModel,variables,by = "Site")
head(FDModel);dim(FDModel)

FDModel <- merge(FDModel,coords,by = "Site")
head(FDModel);dim(FDModel)



#Functional Richness----
##Step 1: Design Variables----
#position, age and day - all the various combinations of these 

FDModel$Day_Scaled <- scale(FDModel$Day_Sampled)
FDModel$Age_Scaled <- scale(FDModel$Crop_Age_Days)


head(FDModel);dim(FDModel)
str(FDModel)

FUNRich_null <- glmmTMB(Fun_Rich ~ 1 + (1 | Field), family = poisson, data = FDModel)

FUNRich_P <- glmmTMB(Fun_Rich ~ Position + (1 | Field), family = poisson, data = FDModel)
FUNRich_A <- glmmTMB(Fun_Rich ~ Crop_Age_Days + (1 | Field), family = poisson, data = FDModel)
FUNRich_D <- glmmTMB(Fun_Rich ~ Day_Scaled + (1 | Field), family = poisson, data = FDModel) 

FUNRich_PA <- glmmTMB(Fun_Rich ~ Position + Crop_Age_Days + (1 | Field), family = poisson, data = FDModel)
FUNRich_PD <- glmmTMB(Fun_Rich ~ Position + Day_Scaled + (1 | Field), family = poisson, data = FDModel) 
FUNRich_DA <- glmmTMB(Fun_Rich ~ Day_Scaled + Crop_Age_Days + (1 | Field), family = poisson, data = FDModel) 

FUNRich_PxA <- glmmTMB(Fun_Rich ~ Position * Age_Scaled + (1 | Field), family = poisson, data = FDModel) 
FUNRich_PxD <- glmmTMB(Fun_Rich ~ Position * Day_Scaled + (1 | Field), family = poisson, data = FDModel) 
FUNRich_DxA <- glmmTMB(Fun_Rich ~ Day_Scaled * Age_Scaled + (1 | Field), family = poisson, data = FDModel) 

FUNRich_PAD <- glmmTMB(Fun_Rich ~ Position + Age_Scaled + Day_Scaled + (1 | Field), family = poisson, data = FDModel) 
FUNRich_PxAD <- glmmTMB(Fun_Rich ~ Position * Age_Scaled + Day_Scaled + (1 | Field), family = poisson, data = FDModel)
FUNRich_PxDA <- glmmTMB(Fun_Rich ~ Position * Day_Scaled + Age_Scaled + (1 | Field), family = poisson, data = FDModel)
FUNRich_PAxD <- glmmTMB(Fun_Rich ~ Position + Day_Scaled * Age_Scaled + (1 | Field), family = poisson, data = FDModel) 

#collect models
FUNrichmodlist <- list("null" = FUNRich_null, "P" = FUNRich_P, 
                       "A" = FUNRich_A, "D" = FUNRich_D, 
                       "PA" = FUNRich_PA, "PD" = FUNRich_PD, 
                       "DA" = FUNRich_DA,"PxA" = FUNRich_PxA, 
                       "PxD" = FUNRich_PxD,"DxA" = FUNRich_DxA,
                       "PAD" = FUNRich_PAD, "PxAD" = FUNRich_PxAD, 
                       "PxDA" = FUNRich_PxDA, "PAxD" = FUNRich_PAxD)

aictab(FUNrichmodlist)
#Top model is Age

##Step 2: Environmental Variables----

FDModel$NDVI1km_Scaled <- scale(FDModel$NDVIsum_1km)
FDModel$Field_Area_Scaled <- scale(FDModel$Field_Area_m2)

head(FDModel)

FUNRich_Height <- glmmTMB(Fun_Rich ~ Crop_Age_Days + Height + (1 | Field), family = poisson, data = FDModel) 
FUNRich_GC <- glmmTMB(Fun_Rich ~ Crop_Age_Days + GC + (1 | Field), family = poisson, data = FDModel)

FUNRich_FieldArea <- glmmTMB(Fun_Rich ~ Crop_Age_Days + Field_Area_Scaled + (1 | Field), family = poisson, data = FDModel)
FUNRich_NDVIfield <- glmmTMB(Fun_Rich ~ Crop_Age_Days + NDVImean_Field  + (1 | Field), family = poisson, data = FDModel)  

FUNRich_Rip <- glmmTMB(Fun_Rich ~ Crop_Age_Days + X1km_Rip_Prop + (1 | Field), family = poisson, data = FDModel) 
FUNRich_Crops <- glmmTMB(Fun_Rich ~ Crop_Age_Days + X1km_Prop_Crops + (1 | Field), family = poisson, data = FDModel) 
FUNRich_NDVI1km <- glmmTMB(Fun_Rich ~ Crop_Age_Days + NDVI1km_Scaled + (1 | Field), family = poisson, data = FDModel) 


FUNrichmodlist2 <- list("null" = FUNRich_null, 
                        "Height" = FUNRich_Height, 
                        "GC" = FUNRich_GC,
                        "Field Area" = FUNRich_FieldArea, 
                        "Field NDVI" = FUNRich_NDVIfield,
                        "Rip" = FUNRich_Rip,
                        "Crop" = FUNRich_Crops, 
                        "NDVI 1km" = FUNRich_NDVI1km)
aictab(FUNrichmodlist2)
#Rip is within 2 AICc


##Step 3: Check Spatial Autocorrelation----

FRfield_numbers <- unique(FDModel$ID)


FRmodel_residuals <- simulateResiduals(FUNRich_Rip)
FRspatial_result <- data.frame(
  field = rep(NA, length(FRfield_numbers)),
  statistic = rep(NA, length(FRfield_numbers)),
  p_value = rep(NA, length(FRfield_numbers)),
  method = rep(NA_character_, length(FRfield_numbers)),
  stringsAsFactors = FALSE)

s <- 1

for (f in FRfield_numbers) {
  
  cat("Field", f, "\n") #What field is it doing?
  
  #Extracting specific residuals for individual fields
  FRfield_indices <- which(FDModel$ID == f)
  FRfield_residuals <- FRmodel_residuals
  FRfield_residuals$scaledResiduals <- 
    FRmodel_residuals$scaledResiduals[FRfield_indices]
  FRfield_residuals$fittedPredictedResponse <- 
    FRmodel_residuals$fittedPredictedResponse[FRfield_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  FRspatial_test <- testSpatialAutocorrelation(FRfield_residuals, 
                   x = FDModel$X_Cor[FDModel$ID == f], 
                   y = FDModel$Y_Cor[FDModel$ID == f])
  
  
  FRspatial_result$field [s] <- f
  FRspatial_result$statistic [s] <- FRspatial_test$statistic[1] 
  FRspatial_result$p_value [s] <- FRspatial_test$p.value
  FRspatial_result$method [s] <- FRspatial_test$method
  
  s <- s + 1
  
}

head(FRspatial_result);dim(FRspatial_result)
length(unique(FDModel$ID))

FRspatial_result


Spatial_auto_FR_Field <- FRspatial_result
Spatial_auto_FR_Rip <- FRspatial_result

#Spatial Autocorrelation found in one fields/surveys for both models

#Functional Diversity----
##Step 1: Design Variables----

head(FDModel);dim(FDModel)
str(FDModel)


FUNDiv_null <- glmer(Fun_Div ~ 1 + (1 | Field), family = Gamma(link = "log"), data = FDModel)

FUNDiv_P <- glmer(Fun_Div ~ Position + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_A <- glmer(Fun_Div ~ Crop_Age_Days + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_D <- glmer(Fun_Div ~ Day_Sampled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 

FUNDiv_PA <- glmer(Fun_Div ~ Position + Crop_Age_Days+ (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_PD <- glmer(Fun_Div ~ Position + Day_Sampled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_DA <- glmer(Fun_Div ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = Gamma(link = "log"), data = FDModel) 

FUNDiv_PxA <- glmer(Fun_Div ~ Position * Crop_Age_Days + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_PxD <- glmer(Fun_Div ~ Position * Day_Sampled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_DxA <- glmer(Fun_Div ~ Day_Scaled * Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 

FUNDiv_PAD <- glmer(Fun_Div ~ Position + Crop_Age_Days + Day_Sampled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_PxAD <- glmer(Fun_Div ~ Position * Crop_Age_Days + Day_Sampled + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_PxDA <- glmer(Fun_Div ~ Position * Day_Sampled + Crop_Age_Days + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_PAxD <- glmer(Fun_Div ~ Position + Age_Scaled * Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel)



FUndivlist <- list("null" = FUNDiv_null, "P" = FUNDiv_P, 
                   "A" = FUNDiv_A, "D" = FUNDiv_D,"PA" = FUNDiv_PA,
                   "PD" = FUNDiv_PD, "DA" = FUNDiv_DA,
                   "PxA" = FUNDiv_PxA, "PxD" = FUNDiv_PxD, 
                   "DxA" = FUNDiv_DxA, "PAD" = FUNDiv_PAD, 
                   "PxAD" = FUNDiv_PxAD, "PxDA" = FUNDiv_PxDA, 
                   "PAxD" = FUNDiv_PAxD)

aictab(FUndivlist)
#Top Model is Day


##Step 2: Environmental Variables----

head(FDModel)

FUNDiv_Height <- glmer(Fun_Div ~ Day_Sampled + Height + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_GC <- glmer(Fun_Div ~ Day_Sampled + GC + (1 | Field), family = Gamma(link = "log"), data = FDModel)

FUNDiv_FieldArea <- glmer(Fun_Div ~ Day_Sampled + Field_Area_Scaled  + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_FieldNDVI <- glmer(Fun_Div ~ Day_Sampled + NDVImean_Field + (1 | Field), family = Gamma(link = "log"), data = FDModel)

FUNDiv_Rip <- glmer(Fun_Div ~ Day_Sampled + X1km_Rip_Prop + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_Crops <- glmer(Fun_Div ~ Day_Sampled + X1km_Prop_Crops + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_NDVI1km <- glmer(Fun_Div ~ Day_Sampled + NDVIsum_1km + (1 | Field), family = Gamma(link = "log"), data = FDModel) 

FUNdivlist2 <- list("null" = FUNDiv_null, "height" = FUNDiv_Height, 
                 "GC" = FUNDiv_GC, "Field Area" = FUNDiv_FieldArea, 
                 "Field NDVI" = FUNDiv_FieldNDVI, "Rip" = FUNDiv_Rip,
                 "Crop" = FUNDiv_Crops, "NDVI 1km" = FUNDiv_NDVI1km)
aictab(FUNdivlist2)

#top model is GC
#Crops and field area is within 2 AICc's


##Step 3: Check Spatial Autocorrelation----

FDfield_numbers <- unique(FDModel$ID)


FDmodel_residuals <- simulateResiduals(FUNDiv_FieldArea)
FDspatial_result <- data.frame(
  field = rep(NA, length(FDfield_numbers)),
  statistic = rep(NA, length(FDfield_numbers)),
  p_value = rep(NA, length(FDfield_numbers)),
  method = rep(NA_character_, length(FDfield_numbers)),
  stringsAsFactors = FALSE)

s <- 1

for (f in FDfield_numbers) {
  
  cat("Field", f, "\n") #What field is it doing?
  
  #Extracting specific residuals for individual fields
  FDfield_indices <- which(FDModel$ID == f)
  FDfield_residuals <- FDmodel_residuals
  FDfield_residuals$scaledResiduals <- 
    FDmodel_residuals$scaledResiduals[FDfield_indices]
  FDfield_residuals$fittedPredictedResponse <- 
    FDmodel_residuals$fittedPredictedResponse[FDfield_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  FDspatial_test <- testSpatialAutocorrelation(FDfield_residuals, 
                                               x = FDModel$X_Cor[FDModel$ID == f], 
                                               y = FDModel$Y_Cor[FDModel$ID == f])
  
  
  FDspatial_result$field [s] <- f
  FDspatial_result$statistic [s] <- FDspatial_test$statistic[1] 
  FDspatial_result$p_value [s] <- FDspatial_test$p.value
  FDspatial_result$method [s] <- FDspatial_test$method
  
  s <- s + 1
  
}

head(FDspatial_result);dim(FDspatial_result)
length(unique(FDModel$ID))

FDspatial_result


Spatial_auto_FD_GC <- FDspatial_result
Spatial_auto_FD_Crops <- FDspatial_result
Spatial_auto_FD_Field <- FDspatial_result


#No Spatial Autocorrelation found in any fields/surveys

##Step 4: Predictions----

#GC
summary(FUNDiv_GC)

FUNPredictions_Day <- seq(min(FDModel$Day_Sampled),max(FDModel$Day_Sampled),length.out=20)
FUNPredictions_GC <- seq(min(FDModel$GC),max(FDModel$GC),length.out=20)


FUNdivpred <- expand.grid(Day_Sampled = FUNPredictions_Day, GC = FUNPredictions_GC)
head(FUNdivpred);dim(FUNdivpred)

FUNdivpred1 <- predict(object = FUNDiv_GC,newdata= FUNdivpred,se.fit = T, type = "link",re.form = NA)

FUNdivpred2<-data.frame(FUNdivpred,fit.link=FUNdivpred1$fit,se.link=FUNdivpred1$se.fit)

FUNdivpred2$lci.link<-FUNdivpred2$fit.link-(1.96*FUNdivpred2$se.link)
FUNdivpred2$uci.link<-FUNdivpred2$fit.link+(1.96*FUNdivpred2$se.link)

FUNdivpred2$fit<-exp(FUNdivpred2$fit.link)
FUNdivpred2$se<-exp(FUNdivpred2$se.link)
FUNdivpred2$lci<-exp(FUNdivpred2$lci.link)
FUNdivpred2$uci<-exp(FUNdivpred2$uci.link)

head(FUNdivpred2);dim(FUNdivpred2)

#Crops

summary(FUNDiv_Crops)

FUNPredictions_Crop <- seq(min(FDModel$X1km_Prop_Crops),max(FDModel$X1km_Prop_Crops),length.out=20)


FUNdivpred3 <- expand.grid(Day_Sampled = FUNPredictions_Day, X1km_Prop_Crops = FUNPredictions_Crop)
head(FUNdivpred3);dim(FUNdivpred3)

FUNdivpred4 <- predict(object = FUNDiv_Crops,newdata= FUNdivpred3,se.fit = T, type = "link",re.form = NA)

FUNdivpred5<-data.frame(FUNdivpred3,fit.link=FUNdivpred4$fit,se.link=FUNdivpred4$se.fit)

FUNdivpred5$lci.link<-FUNdivpred5$fit.link-(1.96*FUNdivpred5$se.link)
FUNdivpred5$uci.link<-FUNdivpred5$fit.link+(1.96*FUNdivpred5$se.link)

FUNdivpred5$fit<-exp(FUNdivpred5$fit.link)
FUNdivpred5$se<-exp(FUNdivpred5$se.link)
FUNdivpred5$lci<-exp(FUNdivpred5$lci.link)
FUNdivpred5$uci<-exp(FUNdivpred5$uci.link)

head(FUNdivpred5);dim(FUNdivpred5)


#Field

summary(FUNDiv_FieldArea)

FUNPredictions_FieldScaled <- seq(min(FDModel$Field_Area_Scaled),max(FDModel$Field_Area_Scaled),length.out=20)


FUNdivpred6 <- expand.grid(Day_Sampled = FUNPredictions_Day, Field_Area_Scaled = FUNPredictions_FieldScaled)
head(FUNdivpred6);dim(FUNdivpred6)

FUNdivpred7 <- predict(object = FUNDiv_FieldArea,newdata= FUNdivpred6,se.fit = T, type = "link",re.form = NA)

FUNdivpred8<-data.frame(FUNdivpred6,fit.link=FUNdivpred7$fit,se.link=FUNdivpred7$se.fit)

FUNdivpred8$lci.link<-FUNdivpred8$fit.link-(1.96*FUNdivpred8$se.link)
FUNdivpred8$uci.link<-FUNdivpred8$fit.link+(1.96*FUNdivpred8$se.link)

FUNdivpred8$fit<-exp(FUNdivpred8$fit.link)
FUNdivpred8$se<-exp(FUNdivpred8$se.link)
FUNdivpred8$lci<-exp(FUNdivpred8$lci.link)
FUNdivpred8$uci<-exp(FUNdivpred8$uci.link)

head(FUNdivpred8);dim(FUNdivpred8)


##Step 5: Visualisation----

#GC

summary(FUNDiv_GC)
head(FUNdivpred2);dim(FUNdivpred2)


GG <- FUNdivpred2$GC == FUNPredictions_GC[10]
G_G <- FUNdivpred2$Day_Sampled == FUNPredictions_Day[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Day_Sampled,y = FDModel$Fun_Div,xlab = "Day Sampled",ylab = 'Trait Group Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95)
mtext(side=3,line=0,at = 10,'a)',cex=1.1)

polygon(x = c(FUNdivpred2$Day_Sampled[GG],rev(FUNdivpred2$Day_Sampled[GG])), y = c(FUNdivpred2$lci[GG],rev(FUNdivpred2$uci[GG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=FUNdivpred2$Day_Sampled[GG],y = FUNdivpred2$fit[GG],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$GC,y = FDModel$Fun_Div,xlab = "Ground Cover (%)",ylab = 'Trait Group Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 4,'b)',cex=1.1)

polygon(x = c(FUNdivpred2$GC[G_G],rev(FUNdivpred2$GC[G_G])), y = c(FUNdivpred2$lci[G_G],rev(FUNdivpred2$uci[G_G])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=FUNdivpred2$GC[G_G],y = FUNdivpred2$fit[G_G],lwd = 2,col = 'grey30')

#Crops

summary(FUNDiv_Crops)
head(FUNdivpred5);dim(FUNdivpred5)


CC <- FUNdivpred5$X1km_Prop_Crops  == FUNPredictions_Crop[10]
C_C <- FUNdivpred5$Day_Sampled == FUNPredictions_Day[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Day_Sampled,y = FDModel$Fun_Div,xlab = "Day Sampled",ylab = 'Trait Group Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95)
mtext(side=3,line=0,at = 10,'a)',cex=1.1)

polygon(x = c(FUNdivpred5$Day_Sampled[CC],rev(FUNdivpred5$Day_Sampled[CC])), y = c(FUNdivpred5$lci[CC],rev(FUNdivpred5$uci[CC])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=FUNdivpred5$Day_Sampled[CC],y = FUNdivpred5$fit[CC],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$X1km_Prop_Crops,y = FDModel$Fun_Div,xlab = "Crops within 1km (%)",ylab = 'Trait Group Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 79,'b)',cex=1.1)

polygon(x = c(FUNdivpred5$X1km_Prop_Crops[C_C],rev(FUNdivpred5$X1km_Prop_Crops[C_C])), y = c(FUNdivpred5$lci[C_C],rev(FUNdivpred5$uci[C_C])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=FUNdivpred5$X1km_Prop_Crops[C_C],y = FUNdivpred5$fit[C_C],lwd = 2,col = 'grey30')


#Field Area
summary(FUNDiv_FieldArea)
head(FUNdivpred8);dim(FUNdivpred8)


FF <- FUNdivpred8$Field_Area_Scaled == FUNPredictions_FieldScaled[10]
F_F <- FUNdivpred8$Day_Sampled == FUNPredictions_Day[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Day_Sampled,y = FDModel$Fun_Div,xlab = "Day Sampled",ylab = 'Trait Group Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95)
mtext(side=3,line=0,at = 10,'a)',cex=1.1)

polygon(x = c(FUNdivpred8$Day_Sampled[FF],rev(FUNdivpred8$Day_Sampled[FF])), y = c(FUNdivpred8$lci[FF],rev(FUNdivpred8$uci[FF])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=FUNdivpred8$Day_Sampled[FF],y = FUNdivpred8$fit[FF],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$Field_Area_Scaled,y = FDModel$Fun_Div,xlab = "Field Size (ha)",ylab = 'Trait Group Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(FUNdivpred8$Field_Area_Scaled),to=max(FUNdivpred8$Field_Area_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1),cex.axis=1)
mtext(side=3,line=0,at = -2.3,'b)',cex=1.1)

polygon(x = c(FUNdivpred8$Field_Area_Scaled[F_F],rev(FUNdivpred8$Field_Area_Scaled[F_F])), y = c(FUNdivpred8$lci[F_F],rev(FUNdivpred8$uci[F_F])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=FUNdivpred8$Field_Area_Scaled[F_F],y = FUNdivpred8$fit[F_F],lwd = 2,col = 'grey30')


#END----

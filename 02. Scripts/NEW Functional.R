
#Author: Rhiannon Bird
#Written under version R 4.5.1

#Used Methods for assessing functional responses to environmental gradients -- Kleyer et al. -- June 22, 2009

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
                                   ">10","No_Size"),
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

##Percentage of co-Inertia for each axis
100*rlq1$eig/sum(rlq1$eig)

#To interpret the results, correlations can be computed:
## weighted correlations axes / env.
t(pca.env$tab)%*%(diag(pca.env$lw))%*%as.matrix(rlq1$mR)

##weighted correlations axes / traits.
t(pca.traits$tab)%*%(diag(pca.traits$lw))%*%as.matrix(rlq1$mQ)

##correlations traits / env.
rlq1$tab

#Biplot (traits and enviro) ----

#here's a prelim plot
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
s.arrow(rlq1$c1, xlim=c(-1,1), boxes = FALSE)
s.label(rlq1$li, add.plot=T, clab=1.5)


#Classifying scores to obtain functional groups----
hc2 <- hclust(dist(rlq1$lQ), method = "ward.D")
dev.new(height=20,width=40,dpi=80,pointsize=14,noRStudioGD = T)
plot(hc2) #uninterpretable

#Calinsky-Harabasz criteria to find best partition 
#calinski didn't work so dchanged it to calinhara 

ntest <- 8
res <- rep(0,ntest - 1)

for (i in 2:ntest){
  fac <- cutree(hc2, k = i)
  res[i-1] <- as.numeric(calinhara(rlq1$lQ, fac))
}

#Trait group numbers----
dev.new(height=5,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mfrow=c(1,2))
plot(2:ntest, res, type='b', pch=20, xlab="Number of groups", ylab = "C-H index")
mtext(side=3,line=0,at = 1.5,'a)',cex=1.1)
plot(3:ntest, diff(res), type='b', pch=20, xlab="Number of groups", ylab = "Diff in C-H index")
mtext(side=3,line=0,at = 2.7,'b)',cex=1.1)


spe.group2 <- as.factor(cutree(hc2, k = which.max(res) +1))
levels(spe.group2) <- c("G","F",'E',"C","B","D","A")
spe.group2 <- factor(spe.group2, levels=c("A","B","C","D","E","F","G"))

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

#Making it more readable
l1_short <- rlq1$l1
rownames(l1_short) <- c("     H", "GC", "In", "Out","Day", "A           ", "FA", "C", "NF", "N1km","R")
l1_labels <- l1_short * 1.30

dev.new(height=10, width=10, dpi=80, pointsize=14, noRStudioGD = T)
ade4::s.class(rlq1$lQ, spe.group2, col = 1:nlevels(spe.group2))
s.arrow(rlq1$l1, add.plot = T, clab = 0)
s.label(l1_labels, add.plot = T, clab = 0.7, boxes = T)

#Adjust these to match environmental variables in order
#Check what they are first:
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
  #Total sum of squares
  grand_mean <- mean(as.numeric(trait))
  TSS <- sum((as.numeric(trait) - grand_mean)^2)
  
  #Between-group sum of squares
  group_means <- tapply(as.numeric(trait), spe.group2, mean)
  group_sizes <- table(spe.group2)
  BSS <- sum(group_sizes * (group_means - grand_mean)^2)
  
  #Eta-squared (proportion of variance explained by groups)
  eta_sq <- BSS / TSS
  return(eta_sq)
})

#View the results
print(eta2)

#Sort to see which traits are most associated with the groups
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
#this isn't a great way to interpret what each group contains


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


#Find most common value for each group
group_profiles <- data.frame(Group = levels(spe.group2))

for(trait in names(Q_matrix)) {
  modal_values <- tapply(Q_matrix[, trait], spe.group2, function(x) {
    names(sort(table(x), decreasing = TRUE))[1]  # Most common value
  })
  group_profiles[, trait] <- modal_values
}

group_profiles
#Interesting but I think the proportions table (above) and the heat map below give a better idea of the groups as some are a combination of two traits together

#Create a heatmap of group-trait associations----

#Calculate proportions for each trait-group combination
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

#Trying to get a better order for each functional group traits

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


heatmap_labs <- c("Hunting" = "a)",
                  "Order" = "b)",
                  "Size" = "c)",
                  'Trophic' = "d)")

dev.new(height=20,width=15,dpi=80,pointsize=14,noRStudioGD = T)
ggplot(heatmap_data, aes(x = Group, y = Traits, fill = Proportion)) +
  geom_tile() +
  facet_wrap(~Functional_Group, scales = "free_y", 
             nrow = 4, ncol = 1, 
             labeller = as_labeller(heatmap_labs)) +
  scale_fill_gradient(low = "white", high = "black") + 
  theme_minimal(base_size = 16) +
  theme(strip.text = element_text(hjust = 0))


#Extract the trait groups----

Trait_Group <- data.frame(Morphospecies = names(spe.group2), trait_group = spe.group2)
head(Trait_Group);dim(Trait_Group)


invert_trait_group <- merge(Trait_Group,invert_filtered, by = "Morphospecies")
head(invert_trait_group);dim(invert_trait_group)
colnames(invert_trait_group)[colnames(invert_trait_group) == "ID"] <- "Site"

#Creating Modelling data and calculating functuional Diversity 

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
FDModel$Crops_Scaled <- scale(FDModel$X1km_Prop_Crops)
#scaling to allow models to converge

head(FDModel)

FUNRich_Height <- glmmTMB(Fun_Rich ~ Crop_Age_Days + Height + (1 | Field), family = poisson, data = FDModel) 
FUNRich_GC <- glmmTMB(Fun_Rich ~ Crop_Age_Days + GC + (1 | Field), family = poisson, data = FDModel)

FUNRich_FieldArea <- glmmTMB(Fun_Rich ~ Crop_Age_Days + Field_Area_Scaled + (1 | Field), family = poisson, data = FDModel)
FUNRich_NDVIfield <- glmmTMB(Fun_Rich ~ Crop_Age_Days + NDVImean_Field  + (1 | Field), family = poisson, data = FDModel)  

FUNRich_Rip <- glmmTMB(Fun_Rich ~ Crop_Age_Days + X1km_Rip_Prop + (1 | Field), family = poisson, data = FDModel) 
FUNRich_Crops <- glmmTMB(Fun_Rich ~ Crop_Age_Days + Crops_Scaled + (1 | Field), family = poisson, data = FDModel) 
FUNRich_NDVI1km <- glmmTMB(Fun_Rich ~ Crop_Age_Days + NDVI1km_Scaled + (1 | Field), family = poisson, data = FDModel) 


FUNrichmodlist2 <- list("null" = FUNRich_null, 
                        "Height" = FUNRich_Height, 
                        "GC" = FUNRich_GC,
                        "Field Area" = FUNRich_FieldArea, 
                        "Field NDVI" = FUNRich_NDVIfield,
                        "Rip" = FUNRich_Rip,
                        "Crop" = FUNRich_Crops, 
                        "NDVI 1km" = FUNRich_NDVI1km,
                        "A" = FUNRich_A)
aictab(FUNrichmodlist2)
#field area is top model and none within 2 AICcs

##Step 3: Check Spatial Autocorrelation----

FRfield_numbers <- unique(FDModel$ID)


FRmodel_residuals <- simulateResiduals(FUNRich_FieldArea)
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
#Spatial Autocorrelation found in one fields/surveys but autocorrelation is negligible

write.xlsx(Spatial_auto_FR_Field, 'SpatialResult.xlsx')


##Step 4: Predictions----

summary(FUNRich_FieldArea)

FUNrichpred <- expand.grid(Crop_Age_Days = FUNPredictions_Age, Field_Area_Scaled = FUNPredictions_FieldScaled)
head(FUNrichpred);dim(FUNrichpred)

FUNrichpred1 <- predict(object = FUNRich_FieldArea,newdata= FUNrichpred,se.fit = T, type = "link",re.form = NA)

FUNrichpred2<-data.frame(FUNrichpred,fit.link=FUNrichpred1$fit,se.link=FUNrichpred1$se.fit)

FUNrichpred2$lci.link<-FUNrichpred2$fit.link-(1.96*FUNrichpred2$se.link)
FUNrichpred2$uci.link<-FUNrichpred2$fit.link+(1.96*FUNrichpred2$se.link)

FUNrichpred2$fit<-exp(FUNrichpred2$fit.link)
FUNrichpred2$se<-exp(FUNrichpred2$se.link)
FUNrichpred2$lci<-exp(FUNrichpred2$lci.link)
FUNrichpred2$uci<-exp(FUNrichpred2$uci.link)

head(FUNrichpred2);dim(FUNrichpred2)

##Step 5: Visualisation----

summary(FUNDiv_FieldArea)
head(FUNrichpred2);dim(FUNrichpred2)


AA <- FUNrichpred2$Field_Area_Scaled == FUNPredictions_FieldScaled[10]
A_A <- FUNrichpred2$Crop_Age_Days == FUNPredictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$Fun_Div,xlab = "Crop Age",ylab = 'Trait Group Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95)
mtext(side=3,line=0,at = 47,'a)',cex=1.1)

polygon(x = c(FUNrichpred2$Crop_Age_Days[AA],rev(FUNrichpred2$Crop_Age_Days[AA])), y = c(FUNrichpred2$lci[AA],rev(FUNrichpred2$uci[AA])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=FUNrichpred2$Crop_Age_Days[AA],y = FUNrichpred2$fit[AA],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$Field_Area_Scaled,y = FDModel$Fun_Div,xlab = "Field Size (ha)",ylab = 'Trait Group Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt = 'n')
axis(side=1, at=seq(from=min(FUNrichpred2$Field_Area_Scaled),to=max(FUNrichpred2$Field_Area_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1))
mtext(side=3,line=0,at = -2.2,'b)',cex=1.1)

polygon(x = c(FUNrichpred2$Field_Area_Scaled[A_A],rev(FUNrichpred2$Field_Area_Scaled[A_A])), y = c(FUNrichpred2$lci[A_A],rev(FUNrichpred2$uci[A_A])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=FUNrichpred2$Field_Area_Scaled[A_A],y = FUNrichpred2$fit[A_A],lwd = 2,col = 'grey30')



#Functional Diversity----
##Step 1: Design Variables----

head(FDModel);dim(FDModel)
str(FDModel)


FUNDiv_null <- glmer(Fun_Div ~ 1 + (1 | Field), family = Gamma(link = "log"), data = FDModel)

FUNDiv_P <- glmer(Fun_Div ~ Position + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_A <- glmer(Fun_Div ~ Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_D <- glmer(Fun_Div ~ Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) #model failed to converge

FUNDiv_PA <- glmer(Fun_Div ~ Position + Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_PD <- glmer(Fun_Div ~ Position + Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_DA <- glmer(Fun_Div ~ Day_Scaled + Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 

FUNDiv_PxA <- glmer(Fun_Div ~ Position * Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_PxD <- glmer(Fun_Div ~ Position * Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_DxA <- glmer(Fun_Div ~ Day_Scaled * Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 

FUNDiv_PAD <- glmer(Fun_Div ~ Position + Age_Scaled + Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) #failed to converge
FUNDiv_PxAD <- glmer(Fun_Div ~ Position * Age_Scaled + Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) #Failed to converge
FUNDiv_PxDA <- glmer(Fun_Div ~ Position * Day_Scaled + Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 
FUNDiv_PAxD <- glmer(Fun_Div ~ Position + Age_Scaled * Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel)



FUndivlist <- list("null" = FUNDiv_null, "P" = FUNDiv_P, 
                   "A" = FUNDiv_A, "PA" = FUNDiv_PA,
                   "PD" = FUNDiv_PD, "DA" = FUNDiv_DA,
                   "PxA" = FUNDiv_PxA, "PxD" = FUNDiv_PxD, 
                   "DxA" = FUNDiv_DxA,  "PxDA" = FUNDiv_PxDA, 
                   "PAxD" = FUNDiv_PAxD)

aictab(FUndivlist)
#Top Model is Day x Age but Null is within 2 AICc


##Step 2: Environmental Variables----

FDModel$Height_Scaled <- scale(FDModel$Height)
FDModel$GC_Scaled <- scale(FDModel$GC)


head(FDModel)

FUNDiv_Height <- glmer(Fun_Div ~ Height_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_GC <- glmer(Fun_Div ~ GC_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel)

FUNDiv_FieldArea <- glmer(Fun_Div ~ Field_Area_Scaled  + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_FieldNDVI <- glmer(Fun_Div ~ NDVImean_Field + (1 | Field), family = Gamma(link = "log"), data = FDModel)

FUNDiv_Rip <- glmer(Fun_Div ~ X1km_Rip_Prop + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_Crops <- glmer(Fun_Div ~ Crops_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel)
FUNDiv_NDVI1km <- glmer(Fun_Div ~ NDVI1km_Scaled + (1 | Field), family = Gamma(link = "log"), data = FDModel) 

FUNdivlist2 <- list("null" = FUNDiv_null, "height" = FUNDiv_Height, 
                 "GC" = FUNDiv_GC, "Field Area" = FUNDiv_FieldArea, 
                 "Field NDVI" = FUNDiv_FieldNDVI,"Rip" = FUNDiv_Rip,
                 "Crop" = FUNDiv_Crops, "NDVI 1km" = FUNDiv_NDVI1km)
aictab(FUNdivlist2)

#top model is field area but null is within 2 AICcs
#will still visalise this model but put it into the supporting info only


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

Spatial_auto_FD_Field <- FDspatial_result

write.xlsx(Spatial_auto_FD_Field, 'SpatialResult.xlsx')


##Step 4: Predictions----


#Field

summary(FUNDiv_FieldArea)

FUNPredictions_FieldScaled <- seq(min(FDModel$Field_Area_Scaled),max(FDModel$Field_Area_Scaled),length.out=20)

FUNdivpred <- data.frame(Field_Area_Scaled = FUNPredictions_FieldScaled)

FUNdivpred2 <- predict(object = FUNDiv_FieldArea,
                       newdata= FUNdivpred,
                       se.fit = T, type = "link",re.form = NA)

FUNdivpred3<-data.frame(FUNdivpred,fit.link=FUNdivpred2$fit,se.link=FUNdivpred2$se.fit)

FUNdivpred3$lci.link<-FUNdivpred3$fit.link-(1.96*FUNdivpred3$se.link)
FUNdivpred3$uci.link<-FUNdivpred3$fit.link+(1.96*FUNdivpred3$se.link)

FUNdivpred3$fit<-exp(FUNdivpred3$fit.link)
FUNdivpred3$se<-exp(FUNdivpred3$se.link)
FUNdivpred3$lci<-exp(FUNdivpred3$lci.link)
FUNdivpred3$uci<-exp(FUNdivpred3$uci.link)

head(FUNdivpred3);dim(FUNdivpred3)


##Step 5: Visualisation----

#Field Area
summary(FUNDiv_FieldArea)
head(FUNdivpred3);dim(FUNdivpred3)


FF <- FUNdivpred3$Field_Area_Scaled == FUNPredictions_FieldScaled[10]
F_F <- FUNdivpred3$Day_Sampled == FUNPredictions_Day[10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Field_Area_Scaled,y = FDModel$Fun_Div,xlab = "Field Size (ha)",ylab = 'Trait Group Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(FUNdivpred3$Field_Area_Scaled),to=max(FUNdivpred3$Field_Area_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1),cex.axis=1)

polygon(x = c(FUNdivpred3$Field_Area_Scaled,rev(FUNdivpred3$Field_Area_Scaled)), y = c(FUNdivpred3$lci,rev(FUNdivpred3$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=FUNdivpred3$Field_Area_Scaled,y = FUNdivpred3$fit,lwd = 2,col = 'grey30')



#END----

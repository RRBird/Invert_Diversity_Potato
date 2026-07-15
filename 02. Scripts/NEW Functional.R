
#Author: Rhiannon Bird
#Written under version R 4.5.1

#Used Methods for assessing functional responses to environmental gradients -- Kleyer et al. -- June 22, 2009

options(scipen = 999) #So R doesn't use scientific notation

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

#Calinski-Harabasz criteria to find best partition 

ntest <- 10
res <- rep(0,ntest - 1)

for (i in 2:ntest){
  fac <- cutree(hc2, k = i)
  res[i-1] <- as.numeric(calinhara(rlq1$lQ, fac))
}

#Trait group numbers----
dev.new(height=5,width=7,dpi=80,pointsize=14,noRStudioGD = T)
par(mfrow=c(1,2))
plot(2:ntest, res, type='b', pch=20, xlab="Number of groups", ylab = "C-H index")
mtext(side=3,line=0,at = 1.2,'a)',cex=1.1)
plot(3:ntest, diff(res), type='b', pch=20, xlab="Number of groups", ylab = "Diff in C-H index")
mtext(side=3,line=0,at = 2.2,'b)',cex=1.1)

#nine groups as levels off between nine and 10


spe.group2 <- as.factor(cutree(hc2, k = which.max(res) +1))
levels(spe.group2) <- c("A","B","C","D","E","F","G","H","I")
spe.group2 <- factor(spe.group2, levels=c("A","B","C","D","E","F","G","H","I"))
summary(spe.group2)

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
levels(heatmap_data$Trait_Value)

heatmap_data$Trait_Value <- gsub("Active_Hunting", "Active", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("Ambush_Hunter", "Ambush", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("No_Size", "No Size", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("Web_Building", "Web", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("0-2.5", "0-2.5mm", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("2.5-5", "2.5-5mm", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub("5-10", "5-10mm", heatmap_data$Trait_Value)
heatmap_data$Trait_Value <- gsub(">10", ">10mm", heatmap_data$Trait_Value)


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
  theme(strip.text = element_text(hjust = 0, vjust = 1))


#Extract the trait groups----

Trait_Group <- data.frame(Morphospecies = names(spe.group2), trait_group = spe.group2)
head(Trait_Group);dim(Trait_Group)


invert_trait_group <- merge(Trait_Group,invert_filtered, by = "Morphospecies")
head(invert_trait_group);dim(invert_trait_group)

#Creating Modelling data and calculating functional Diversity 

FDModel <- data.frame(Site = variables$Site)
head(FDModel);dim(FDModel)

##Richness----

#Trait Group A
TG_A <- invert_trait_group[invert_trait_group$trait_group=="A",]
head(TG_A);dim(TG_A)

FUNrichness_A <- aggregate(Morphospecies ~ Site, data = TG_A, FUN = function(x) length(unique(x)))
FDModel <- merge(FDModel,FUNrichness_A,by = "Site",all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[2] <- "A_Rich"
head(FDModel);dim(FDModel)

FDModel$A_Rich[is.na(FDModel$A_Rich)] <- 0

#Trait Group B
TG_B <- invert_trait_group[invert_trait_group$trait_group=="B",]
head(TG_B);dim(TG_B)

FUNrichness_B <- aggregate(Morphospecies ~ Site, data = TG_B, FUN = function(x) length(unique(x)))
FDModel <- merge(FDModel,FUNrichness_B,by = "Site",all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[3] <- "B_Rich"
head(FDModel);dim(FDModel)

FDModel$B_Rich[is.na(FDModel$B_Rich)] <- 0

#Trait Group C
TG_C <- invert_trait_group[invert_trait_group$trait_group=="C",]
head(TG_C);dim(TG_C)

FUNrichness_C <- aggregate(Morphospecies ~ Site, data = TG_C, FUN = function(x) length(unique(x)))
FDModel <- merge(FDModel,FUNrichness_C,by = "Site",all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[4] <- "C_Rich"

FDModel$C_Rich[is.na(FDModel$C_Rich)] <- 0
head(FDModel);dim(FDModel)

#Trait Group D
TG_D <- invert_trait_group[invert_trait_group$trait_group=="D",]
head(TG_D);dim(TG_D)

FUNrichness_D <- aggregate(Morphospecies ~ Site, data = TG_D, FUN = function(x) length(unique(x)))
FDModel <- merge(FDModel,FUNrichness_D,by = "Site",all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[5] <- "D_Rich"

FDModel$D_Rich[is.na(FDModel$D_Rich)] <- 0
head(FDModel);dim(FDModel)

#Trait Group E
TG_E <- invert_trait_group[invert_trait_group$trait_group=="E",]
head(TG_E);dim(TG_E)

FUNrichness_E <- aggregate(Morphospecies ~ Site, data = TG_E, FUN = function(x) length(unique(x)))
FDModel <- merge(FDModel,FUNrichness_E,by = "Site",all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[6] <- "E_Rich"

FDModel$E_Rich[is.na(FDModel$E_Rich)] <- 0
head(FDModel);dim(FDModel)

#Trait Group F
TG_F <- invert_trait_group[invert_trait_group$trait_group=="F",]
head(TG_F);dim(TG_F)

FUNrichness_F <- aggregate(Morphospecies ~ Site, data = TG_F, FUN = function(x) length(unique(x)))
FDModel <- merge(FDModel,FUNrichness_F,by = "Site",all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[7] <- "F_Rich"

FDModel$F_Rich[is.na(FDModel$F_Rich)] <- 0
head(FDModel);dim(FDModel)

#Trait Group G
TG_G <- invert_trait_group[invert_trait_group$trait_group=="G",]
head(TG_G);dim(TG_G)

FUNrichness_G <- aggregate(Morphospecies ~ Site, data = TG_G, FUN = function(x) length(unique(x)))
FDModel <- merge(FDModel,FUNrichness_G,by = "Site",all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[8] <- "G_Rich"

FDModel$G_Rich[is.na(FDModel$G_Rich)] <- 0
head(FDModel);dim(FDModel)

##Diversity----
v <- 9
TG <- list("A","B","C","D","E","F","G")

for (k in TG) {
  TG_Div <- invert_trait_group[invert_trait_group$trait_group==k,]
  
  FUNdiversity <- aggregate(Morphospecies ~ Site, data = TG_Div, FUN = function(x) diversity(table(x), index = "invsimpson"))
  
  FDModel <- merge(FDModel,FUNdiversity,by = "Site",all.x = T)
  
  colnames(FDModel)[v] <- k
  
  v <- v+1
}

head(FDModel);dim(FDModel)


FDModel$A[is.na(FDModel$A)] <- 0.00001
FDModel$B[is.na(FDModel$B)] <- 0.00001
FDModel$C[is.na(FDModel$C)] <- 0.00001
FDModel$D[is.na(FDModel$D)] <- 0.00001
FDModel$E[is.na(FDModel$E)] <- 0.00001
FDModel$F[is.na(FDModel$F)] <- 0.00001
FDModel$G[is.na(FDModel$G)] <- 0.00001

colnames(FDModel)[9] <- "A_Div"
colnames(FDModel)[10] <- "B_Div"
colnames(FDModel)[11] <- "C_Div"
colnames(FDModel)[12] <- "D_Div"
colnames(FDModel)[13] <- "E_Div"
colnames(FDModel)[14] <- "F_Div"
colnames(FDModel)[15] <- "G_Div"

head(FDModel);dim(FDModel)

##lets check prop 0's----

lapply(FDModel[,2:8], function(x){length(which(x==0))/length(x)})

#okay looks like they all need to be binomial models rather then richness and diversity models

FDModel <- FDModel[,1:8]

#Variables and XY coordinates

FDModel <- merge(FDModel,variables,by = "Site")
head(FDModel);dim(FDModel)

FDModel <- merge(FDModel,coords,by = "Site")
head(FDModel);dim(FDModel)

#update to pres/abs
FDModel$A_Rich <- ifelse(FDModel$A_Rich>0, yes = 1,no = 0)
FDModel$B_Rich <- ifelse(FDModel$B_Rich>0, yes = 1,no = 0)
FDModel$C_Rich <- ifelse(FDModel$C_Rich>0, yes = 1,no = 0)
FDModel$D_Rich <- ifelse(FDModel$D_Rich>0, yes = 1,no = 0)
FDModel$E_Rich <- ifelse(FDModel$E_Rich>0, yes = 1,no = 0)
FDModel$F_Rich <- ifelse(FDModel$F_Rich>0, yes = 1,no = 0)
FDModel$G_Rich <- ifelse(FDModel$G_Rich>0, yes = 1,no = 0)

colnames(FDModel)[2] <- "TG_A"
colnames(FDModel)[3] <- "TG_B"
colnames(FDModel)[4] <- "TG_C"
colnames(FDModel)[5] <- "TG_D"
colnames(FDModel)[6] <- "TG_E"
colnames(FDModel)[7] <- "TG_F"
colnames(FDModel)[8] <- "TG_G"


head(FDModel);dim(FDModel)

FDModel$Day_Scaled <- scale(FDModel$Day_Sampled)
FDModel$Age_Scaled <- scale(FDModel$Crop_Age_Days)
FDModel$Field_Size_Scaled <- scale(FDModel$Field_Area_m2)

#Trait Group A----
##Step 1: Modelling Design Variables----



head(FDModel);dim(FDModel)
str(FDModel)

TGA_null <- glmmTMB(TG_A  ~ 1 + (1 | Field), family = binomial, data = FDModel)

TGA_P <- glmmTMB(TG_A ~ Position + (1 | Field), family = binomial, data = FDModel)
TGA_A <- glmmTMB(TG_A ~ Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGA_D <- glmmTMB(TG_A ~ Day_Scaled + (1 | Field), family = binomial, data = FDModel)

TGA_PA <- glmmTMB(TG_A ~ Position + Age_Scaled + (1 | Field), family = binomial, data = FDModel)
TGA_PD <- glmmTMB(TG_A ~ Position + Day_Scaled + (1 | Field), family = binomial, data = FDModel) 
TGA_DA <- glmmTMB(TG_A ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGA_PxA <- glmmTMB(TG_A ~ Position * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGA_PxD <- glmmTMB(TG_A ~ Position * Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGA_DxA <- glmmTMB(TG_A ~ Day_Scaled * Age_Scaled + (1 | Field), family = binomial, data = FDModel) 

TGA_PAD <- glmmTMB(TG_A ~ Position + Age_Scaled + Day_Scaled + (1 | Field), family = binomial, data = FDModel) 
TGA_PxAD <- glmmTMB(TG_A ~ Position * Age_Scaled + Day_Scaled + (1 | Field), family = binomial, data = FDModel) #Convergence problems
TGA_PxDA <- glmmTMB(TG_A ~ Position * Day_Scaled + Age_Scaled + (1 | Field), family = binomial, data = FDModel)
TGA_PAxD <- glmmTMB(TG_A ~ Position + Day_Scaled * Age_Scaled + (1 | Field), family = binomial, data = FDModel) 

#collect models
TGA_modlist <- list("null" = TGA_null, "P" = TGA_P, "A" = TGA_A,
                    'D' = TGA_D, "PA" = TGA_PA, "PD" = TGA_PD, 
                    "DA" = TGA_DA, "PxA" = TGA_PxA, "PxD" =TGA_PxD,
                    "DxA" = TGA_DxA,"PAD" = TGA_PAD,
                    'PxAD'=TGA_PxAD,"PxDA" = TGA_PxDA, 
                    "PAxD" = TGA_PAxD)

aictab(TGA_modlist)
#Top model is null

##Step 2: Environmental Variables----


head(FDModel)

TGA_Height <- glmmTMB(TG_A ~ Height + (1 | Field), family = binomial, data = FDModel)
TGA_GC <- glmmTMB(TG_A ~ GC + (1 | Field), family = binomial, data = FDModel)

TGA_FieldArea <- glmmTMB(TG_A ~ Field_Area_m2 + (1 | Field), family = binomial, data = FDModel)
TGA_NDVIfield <- glmmTMB(TG_A ~ NDVImean_Field  + (1 | Field), family = binomial, data = FDModel)  

TGA_Rip <- glmmTMB(TG_A ~ X1km_Rip_Prop + (1 | Field), family = binomial, data = FDModel) 
TGA_Crops <- glmmTMB(TG_A ~ X1km_Prop_Crops + (1 | Field), family = binomial, data = FDModel) 
TGA_NDVI1km <- glmmTMB(TG_A ~ NDVIsum_1km + (1 | Field), family = binomial, data = FDModel) 


TGA_modlist2 <- list("null" = TGA_null,
                     "Height" = TGA_Height,   
                     "GC" = TGA_GC,
                     "Field Size" = TGA_FieldArea,
                     "Field NDVI" = TGA_NDVIfield,
                     "Rip" = TGA_Rip,
                     "Crop" = TGA_Crops,
                     "NDVI 1km" = TGA_NDVI1km)

aictab(TGA_modlist2)
#Top model is Null
#all except GC is within 2 AICc's of null 

#When null is top model or within 2 AICc's of top model exlcuded from rest of analysis

#Trait Group B----
##Step 1: Modelling Design Variables----

head(FDModel);dim(FDModel)
str(FDModel)

TGB_null <- glmmTMB(TG_B  ~ 1 + (1 | Field), family = binomial, data = FDModel)

TGB_P <- glmmTMB(TG_B ~ Position + (1 | Field), family = binomial, data = FDModel)
TGB_A <- glmmTMB(TG_B ~ Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGB_D <- glmmTMB(TG_B ~ Day_Sampled + (1 | Field), family = binomial, data = FDModel) 

TGB_PA <- glmmTMB(TG_B ~ Position + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGB_PD <- glmmTMB(TG_B ~ Position + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGB_DA <- glmmTMB(TG_B ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGB_PxA <- glmmTMB(TG_B ~ Position * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGB_PxD <- glmmTMB(TG_B ~ Position * Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGB_DxA <- glmmTMB(TG_B ~ Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGB_PAD <- glmmTMB(TG_B ~ Position + Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGB_PxAD <- glmmTMB(TG_B ~ Position * Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel)
TGB_PxDA <- glmmTMB(TG_B ~ Position * Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGB_PAxD <- glmmTMB(TG_B ~ Position + Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

#collect models
TGB_modlist <- list("null" = TGB_null, "P" = TGB_P,"A" = TGB_A,
                    "D" = TGB_D, "PA" = TGB_PA, "PD" = TGB_PD, 
                    "DA" = TGB_DA, "PxA" = TGB_PxA, 
                    "PxD" =TGB_PxD,"DxA" = TGB_DxA,
                    "PAD" = TGB_PAD,"PxAD" = TGB_PxAD, 
                    "PxDA" = TGB_PxDA, "PAxD" = TGB_PAxD)

aictab(TGB_modlist)
#Top model is position by null is within 2 AICc

##Step 2: Environmental Variables----

head(FDModel)

TGB_Height <- glmmTMB(TG_B ~ Height + (1 | Field), family = binomial, data = FDModel)
TGB_GC <- glmmTMB(TG_B ~ GC + (1 | Field), family = binomial, data = FDModel)

TGB_FieldArea <- glmmTMB(TG_B ~ Field_Size_Scaled + (1 | Field), family = binomial, data = FDModel)
TGB_NDVIfield <- glmmTMB(TG_B ~ NDVImean_Field  + (1 | Field), family = binomial, data = FDModel)  

TGB_Rip <- glmmTMB(TG_B ~ X1km_Rip_Prop + (1 | Field), family = binomial, data = FDModel) 
TGB_Crops <- glmmTMB(TG_B ~ X1km_Prop_Crops + (1 | Field), family = binomial, data = FDModel) 
TGB_NDVI1km <- glmmTMB(TG_B ~ NDVIsum_1km + (1 | Field), family = binomial, data = FDModel)


TGB_modlist2 <- list("null" = TGB_null,
                     "Height" = TGB_Height,
                     "GC" = TGB_GC,
                     "Field Size" = TGB_FieldArea,
                     "Field NDVI" = TGB_NDVIfield,
                     "Rip" = TGB_Rip,
                     "Crop" = TGB_Crops,
                     "NDVI" = TGB_NDVI1km)
aictab(TGB_modlist2)
#Top model is Field size (None within 2 AICc's)

##Step 3: Check for Spatial Autocorrelation----

summary(TGB_FieldArea)

TGfield_numbers <- unique(FDModel$ID)

#run before loop
TGmodel_residuals <- simulateResiduals(TGB_FieldArea)
TGspatial_result <- data.frame(
  field = rep(NA, length(TGfield_numbers)),
  statistic = rep(NA, length(TGfield_numbers)),
  p_value = rep(NA, length(TGfield_numbers)),
  method = rep(NA_character_, length(TGfield_numbers)),
  stringsAsFactors = FALSE)

s <- 1

for (f in TGfield_numbers) {
  
  cat("Field", f, "\n") #What field is it doing?
  
  #Extracting specific residuals for individual fields
  TGfield_indices <- which(FDModel$ID == f)
  TGfield_residuals <- TGmodel_residuals
  TGfield_residuals$scaledResiduals <- 
    TGmodel_residuals$scaledResiduals[TGfield_indices]
  TGfield_residuals$fittedPredictedResponse <- 
    TGmodel_residuals$fittedPredictedResponse[TGfield_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  TGspatial_test <- testSpatialAutocorrelation(
    TGfield_residuals, 
    x = FDModel$X_Cor[FDModel$ID == f], 
    y = FDModel$Y_Cor[FDModel$ID == f])
  
  
  TGspatial_result$field [s] <- f
  TGspatial_result$statistic [s] <- TGspatial_test$statistic[1] 
  TGspatial_result$p_value [s] <- TGspatial_test$p.value
  TGspatial_result$method [s] <- TGspatial_test$method
  
  s <- s + 1
  
}

head(TGspatial_result);dim(TGspatial_result)
length(unique(FDModel$ID))

TGspatial_result


TG_B_Spatial <- TGspatial_result
#No spatial autocorrelation found 

write.xlsx(TG_B_Spatial, 'TG_B_Spatial.xlsx')

##Step 4: Predictions----


summary(TGB_FieldArea)

TG_Predictions_FieldSize_Scaled <- seq(min(FDModel$Field_Size_Scaled),max(FDModel$Field_Size_Scaled),length.out=20)


TG_B_pred <- data.frame(Field_Size_Scaled =
                          TG_Predictions_FieldSize_Scaled)
head(TG_B_pred);dim(TG_B_pred)

TG_B_pred1 <- predict(object = TGB_FieldArea,newdata= TG_B_pred,se.fit = T, type = "link",re.form = ~0)
#Predict was struggling but I just went back to modelling and scaled predictor then it was all good

TG_B_pred2<-data.frame(TG_B_pred,fit.link=TG_B_pred1$fit,se.link=TG_B_pred1$se.fit)

TG_B_pred2$lci.link<-TG_B_pred2$fit.link-(1.96*TG_B_pred2$se.link)
TG_B_pred2$uci.link<-TG_B_pred2$fit.link+(1.96*TG_B_pred2$se.link)

TG_B_pred2$fit<-plogis(TG_B_pred2$fit.link)
TG_B_pred2$se<-plogis(TG_B_pred2$se.link)
TG_B_pred2$lci<-plogis(TG_B_pred2$lci.link)
TG_B_pred2$uci<-plogis(TG_B_pred2$uci.link)

head(TG_B_pred2);dim(TG_B_pred2)

##Step 5: Visualisation----

summary(TGB_FieldArea)
head(TG_B_pred2);dim(TG_B_pred2)


dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Field_Size_Scaled,y = FDModel$TG_B,xlab = "Field Size (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,xaxt = 'n',ylim = c(0,1))
axis(side=1, at=seq(from=min(TG_B_pred2$Field_Size_Scaled),to=max(TG_B_pred2$Field_Size_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1),cex.axis=1)

polygon(x = c(TG_B_pred2$Field_Size_Scaled,rev(TG_B_pred2$Field_Size_Scaled)), y = c(TG_B_pred2$lci,rev(TG_B_pred2$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_B_pred2$Field_Size_Scaled,y = TG_B_pred2$fit,lwd = 2,col = 'grey30',lty = 1)


#Trait Group C----
##Step 1: Modelling Design Variables----

head(FDModel);dim(FDModel)
str(FDModel)

TGC_null <- glmmTMB(TG_C ~ 1 + (1 | Field), family = binomial, data = FDModel)

TGC_P <- glmmTMB(TG_C ~ Position + (1 | Field), family = binomial, data = FDModel)
TGC_A <- glmmTMB(TG_C ~ Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGC_D <- glmmTMB(TG_C ~ Day_Sampled + (1 | Field), family = binomial, data = FDModel) 

TGC_PA <- glmmTMB(TG_C ~ Position + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGC_PD <- glmmTMB(TG_C ~ Position + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGC_DA <- glmmTMB(TG_C ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)

TGC_PxA <- glmmTMB(TG_C ~ Position * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGC_PxD <- glmmTMB(TG_C ~ Position * Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGC_DxA <- glmmTMB(TG_C ~ Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)

TGC_PAD <- glmmTMB(TG_C ~ Position + Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGC_PxAD <- glmmTMB(TG_C ~ Position * Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel)
TGC_PxDA <- glmmTMB(TG_C ~ Position * Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGC_PAxD <- glmmTMB(TG_C ~ Position + Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

#collect models
TGC_modlist <- list("null" = TGC_null, "P" = TGC_P,"A" = TGC_A,
                    "D" = TGC_D,"PA" = TGC_PA, "PD" = TGC_PD, 
                    "DA" =TGC_DA, "PxA" =TGC_PxA, "PxD" = TGC_PxD,
                    "PAD" = TGC_PAD, "DxA" = TGC_DxA,
                    "PxAD" = TGC_PxAD,"PxDA" = TGC_PxDA, 
                    "PAxD" = TGC_PAxD)


aictab(TGC_modlist)
#Top model is Position + Age + Day

##Step 2: Environmental Variables----

FDModel$Crop_Scaled <- scale(FDModel$X1km_Prop_Crops)

head(FDModel)

TGC_Height <- glmmTMB(TG_C ~ Position + Crop_Age_Days + Day_Sampled + Height + (1 | Field), family = binomial, data = FDModel)
TGC_GC <- glmmTMB(TG_C ~ Position + Crop_Age_Days + Day_Sampled + GC + (1 | Field), family = binomial, data = FDModel)

TGC_FieldArea <- glmmTMB(TG_C ~ Position + Crop_Age_Days + Day_Sampled + Field_Size_Scaled + (1 | Field), family = binomial, data = FDModel)
TGC_NDVIfield <- glmmTMB(TG_C ~ Position + Crop_Age_Days + Day_Sampled + NDVImean_Field  + (1 | Field), family = binomial, data = FDModel)  

TGC_Rip <- glmmTMB(TG_C ~ Position + Crop_Age_Days + Day_Sampled + X1km_Rip_Prop + (1 | Field), family = binomial, data = FDModel)
TGC_Crops <- glmmTMB(TG_C ~ Position + Crop_Age_Days + Day_Sampled + X1km_Prop_Crops + (1 | Field), family = binomial, data = FDModel) 
TGC_NDVI1km <- glmmTMB(TG_C ~ Position + Crop_Age_Days + Day_Sampled + NDVIsum_1km + (1 | Field), family = binomial, data = FDModel)


TGC_modlist2 <- list("null" = TGC_null,
                     "Height" = TGC_Height,
                     "GC" = TGC_GC,
                     "Field Size" = TGC_FieldArea,
                     "Field NDVI" = TGC_NDVIfield,
                     "Rip" = TGC_Rip,
                     "Crop" = TGC_Crops,
                     "NDVI" = TGC_NDVI1km,
                     "PAD" = TGC_PAD)
aictab(TGC_modlist2)
#Top model is Riparian with field size within 2 AICc's

##Step 3: Check for Spatial Autocorrelation----

#running loon two models this for C so make sure to change over between running the loop
summary(TGC_Rip)
summary(TGC_FieldArea)

#run before loop
TGmodel_residuals <- simulateResiduals(TGC_FieldArea)
TGspatial_result <- data.frame(
  field = rep(NA, length(TGfield_numbers)),
  statistic = rep(NA, length(TGfield_numbers)),
  p_value = rep(NA, length(TGfield_numbers)),
  method = rep(NA_character_, length(TGfield_numbers)),
  stringsAsFactors = FALSE)

s <- 1

for (f in TGfield_numbers) {
  
  cat("Field", f, "\n") #What field is it doing?
  
  #Extracting specific residuals for individual fields
  TGfield_indices <- which(FDModel$ID == f)
  TGfield_residuals <- TGmodel_residuals
  TGfield_residuals$scaledResiduals <- 
    TGmodel_residuals$scaledResiduals[TGfield_indices]
  TGfield_residuals$fittedPredictedResponse <- 
    TGmodel_residuals$fittedPredictedResponse[TGfield_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  TGspatial_test <- testSpatialAutocorrelation(
    TGfield_residuals, 
    x = FDModel$X_Cor[FDModel$ID == f], 
    y = FDModel$Y_Cor[FDModel$ID == f])
  
  
  TGspatial_result$field [s] <- f
  TGspatial_result$statistic [s] <- TGspatial_test$statistic[1] 
  TGspatial_result$p_value [s] <- TGspatial_test$p.value
  TGspatial_result$method [s] <- TGspatial_test$method
  
  s <- s + 1
  
}

head(TGspatial_result);dim(TGspatial_result)
length(unique(FDModel$ID))

TGspatial_result


TG_C_Spatial_Rip <- TGspatial_result
#One field significant but spatial autocorrelation is negligable
TG_C_Spatial_Size <- TGspatial_result
#One field significant but spatial autocorrelation is negligable


write.xlsx(TG_C_Spatial_Rip, 'TG_C_Spatial_Rip.xlsx')
write.xlsx(TG_C_Spatial_Size, 'TG_C_Spatial_Size.xlsx')

##Step 4: Predictions----

#Riparian Predictions
summary(TGC_Rip)

TG_Predictions_Age <- seq(min(FDModel$Crop_Age_Days),max(FDModel$Crop_Age_Days),length.out=20)
TG_Predictions_Day <- seq(min(FDModel$Day_Sampled),max(FDModel$Day_Sampled),length.out=20)
TG_Predictions_Rip <- seq(min(FDModel$X1km_Rip_Prop),max(FDModel$X1km_Rip_Prop),length.out=20)


TG_C_pred <- expand.grid(Position  = unique(FDModel$Position),
                         Crop_Age_Days = TG_Predictions_Age,
                         Day_Sampled = TG_Predictions_Day,
                         X1km_Rip_Prop = TG_Predictions_Rip)

head(TG_C_pred);dim(TG_C_pred)

TG_C_pred1 <- predict(object = TGC_Rip,newdata= TG_C_pred,se.fit = T, type = "link",re.form = ~0)

TG_C_pred2<-data.frame(TG_C_pred,fit.link=TG_C_pred1$fit,se.link=TG_C_pred1$se.fit)

TG_C_pred2$lci.link<-TG_C_pred2$fit.link-(1.96*TG_C_pred2$se.link)
TG_C_pred2$uci.link<-TG_C_pred2$fit.link+(1.96*TG_C_pred2$se.link)

TG_C_pred2$fit<-plogis(TG_C_pred2$fit.link)
TG_C_pred2$se<-plogis(TG_C_pred2$se.link)
TG_C_pred2$lci<-plogis(TG_C_pred2$lci.link)
TG_C_pred2$uci<-plogis(TG_C_pred2$uci.link)

head(TG_C_pred2);dim(TG_C_pred2)

#Field Area predictions
summary(TGC_FieldArea)


TG_C_pred3 <- expand.grid(Position  = unique(FDModel$Position),
                         Crop_Age_Days = TG_Predictions_Age,
                         Day_Sampled = TG_Predictions_Day,
              Field_Size_Scaled = TG_Predictions_FieldSize_Scaled)

head(TG_C_pred3);dim(TG_C_pred3)

TG_C_pred4 <- predict(object = TGC_FieldArea,newdata= TG_C_pred3,se.fit = T, type = "link",re.form = ~0)

TG_C_pred5<-data.frame(TG_C_pred3,fit.link=TG_C_pred4$fit,se.link=TG_C_pred4$se.fit)

TG_C_pred5$lci.link<-TG_C_pred5$fit.link-(1.96*TG_C_pred5$se.link)
TG_C_pred5$uci.link<-TG_C_pred5$fit.link+(1.96*TG_C_pred5$se.link)

TG_C_pred5$fit<-plogis(TG_C_pred5$fit.link)
TG_C_pred5$se<-plogis(TG_C_pred5$se.link)
TG_C_pred5$lci<-plogis(TG_C_pred5$lci.link)
TG_C_pred5$uci<-plogis(TG_C_pred5$uci.link)

head(TG_C_pred5);dim(TG_C_pred5)

##Step 5: Visualisation----

summary(TGC_Rip)
head(TG_C_pred2);dim(TG_C_pred2)

raw_x <- ifelse(FDModel$Position == "Outer", 1, 
                ifelse(FDModel$Position == "Inner", 2, NA))

AA <- TG_C_pred2$Crop_Age_Days == TG_Predictions_Age[10] & TG_C_pred2$Day_Sampled == TG_Predictions_Day [10] & TG_C_pred2$X1km_Rip_Prop ==TG_Predictions_Rip[10]
AAA <- TG_C_pred2$Day_Sampled == TG_Predictions_Day [10] & TG_C_pred2$X1km_Rip_Prop ==TG_Predictions_Rip[10] & TG_C_pred2$Position == "Outer"
A_A <- TG_C_pred2$Crop_Age_Days == TG_Predictions_Age [10] & TG_C_pred2$X1km_Rip_Prop ==TG_Predictions_Rip[10] & TG_C_pred2$Position == "Outer"
A_A_ <- TG_C_pred2$Crop_Age_Days == TG_Predictions_Age [10] & TG_C_pred2$Day_Sampled ==TG_Predictions_Day[10] & TG_C_pred2$Position == "Outer"

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = 1:2,y = TG_C_pred2$fit [AA],xlab = " ",ylab = 'Probability of Occurrence', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,1),xaxt = "n",xlim = c(0,3))
axis(side=1,at=1:2,labels=c('Outer','Inner'))

arrows(x0=1:2, y0=TG_C_pred2$lci [AA],x1=1:2, y1=TG_C_pred2$uci[AA],angle=90,length=0.2, code=3, lwd=2,col = "black")

points(x = jitter(raw_x, factor = 1),y = FDModel$TG_C, pch = 16, cex = 0.4, col = "black")


plot(x = jitter(FDModel$Crop_Age_Days,factor = 1),y = FDModel$TG_C,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_C_pred2$Crop_Age_Days[AAA],rev(TG_C_pred2$Crop_Age_Days[AAA])), y = c(TG_C_pred2$lci[AAA],rev(TG_C_pred2$uci[AAA])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred2$Crop_Age_Days[AAA],y = TG_C_pred2$fit[AAA],lwd = 2,col = 'grey30',lty = 1)


plot(x = jitter(FDModel$Day_Sampled,factor = 1),y = FDModel$TG_C,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_C_pred2$Day_Sampled[A_A],rev(TG_C_pred2$Day_Sampled[A_A])), y = c(TG_C_pred2$lci[A_A],rev(TG_C_pred2$uci[A_A])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred2$Day_Sampled[A_A],y = TG_C_pred2$fit[A_A],lwd = 2,col = 'grey30',lty = 1)


plot(x = jitter(FDModel$X1km_Rip_Prop,factor = 1),y = FDModel$TG_C,xlab = "Proprtion Riparian",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_C_pred2$X1km_Rip_Prop[A_A_],rev(TG_C_pred2$X1km_Rip_Prop[A_A_])), y = c(TG_C_pred2$lci[A_A_],rev(TG_C_pred2$uci[A_A_])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred2$X1km_Rip_Prop[A_A_],y = TG_C_pred2$fit[A_A_],lwd = 2,col = 'grey30',lty = 1)

#Field Model

summary(TGC_FieldArea)
head(TG_C_pred5);dim(TG_C_pred5)

CC <- TG_C_pred5$Crop_Age_Days == TG_Predictions_Age[10] & TG_C_pred5$Day_Sampled == TG_Predictions_Day [10] & TG_C_pred5$Field_Size_Scaled ==TG_Predictions_FieldSize_Scaled[10]
CCC <- TG_C_pred5$Day_Sampled == TG_Predictions_Day [10] & TG_C_pred5$Field_Size_Scaled ==TG_Predictions_FieldSize_Scaled[10] & TG_C_pred5$Position == "Outer"
C_C <- TG_C_pred5$Crop_Age_Days == TG_Predictions_Age [10] & TG_C_pred5$Field_Size_Scaled ==TG_Predictions_FieldSize_Scaled[10] & TG_C_pred5$Position == "Outer"
C_C_ <- TG_C_pred5$Crop_Age_Days == TG_Predictions_Age [10] & TG_C_pred5$Day_Sampled ==TG_Predictions_Day[10] & TG_C_pred5$Position == "Outer"

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = 1:2,y = TG_C_pred5$fit [CC],xlab = " ",ylab = 'Probability of Occurrence', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,1),xaxt = "n",xlim = c(0,3))
axis(side=1,at=1:2,labels=c('Outer','Inner'))

arrows(x0=1:2, y0=TG_C_pred5$lci [CC],x1=1:2, y1=TG_C_pred5$uci[CC],angle=90,length=0.2, code=3, lwd=2,col = "black")

points(x = jitter(raw_x, factor = 1),y = FDModel$TG_C, pch = 16, cex = 0.4, col = "black")


plot(x = jitter(FDModel$Crop_Age_Days,factor = 1),y = FDModel$TG_C,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_C_pred5$Crop_Age_Days[CCC],rev(TG_C_pred5$Crop_Age_Days[CCC])), y = c(TG_C_pred5$lci[CCC],rev(TG_C_pred5$uci[CCC])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred5$Crop_Age_Days[CCC],y = TG_C_pred5$fit[CCC],lwd = 2,col = 'grey30',lty = 1)


plot(x = jitter(FDModel$Day_Sampled,factor = 1),y = FDModel$TG_C,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_C_pred5$Day_Sampled[C_C],rev(TG_C_pred5$Day_Sampled[C_C])), y = c(TG_C_pred5$lci[C_C],rev(TG_C_pred5$uci[C_C])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred5$Day_Sampled[C_C],y = TG_C_pred5$fit[C_C],lwd = 2,col = 'grey30',lty = 1)


plot(x = jitter(FDModel$Field_Size_Scaled,factor = 1),y = FDModel$TG_C,xlab = "Field Size (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_C_pred5$Field_Size_Scaled),to=max(TG_C_pred5$Field_Size_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1),cex.axis=1)

polygon(x = c(TG_C_pred5$Field_Size_Scaled[C_C_],rev(TG_C_pred5$Field_Size_Scaled[C_C_])), y = c(TG_C_pred5$lci[C_C_],rev(TG_C_pred5$uci[C_C_])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred5$Field_Size_Scaled[C_C_],y = TG_C_pred5$fit[C_C_],lwd = 2,col = 'grey30',lty = 1)



#Trait Group D----
##Step 1: Modelling Design Variables----

head(FDModel);dim(FDModel)
str(FDModel)

TGD_null <- glmmTMB(TG_D  ~ 1 + (1 | Field), family = binomial, data = FDModel)


TGD_P <- glmmTMB(TG_D ~ Position + (1 | Field), family = binomial, data = FDModel)
TGD_A <- glmmTMB(TG_D ~ Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGD_D <- glmmTMB(TG_D ~ Day_Sampled + (1 | Field), family = binomial, data = FDModel) 

TGD_PA <- glmmTMB(TG_D ~ Position + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGD_PD <- glmmTMB(TG_D ~ Position + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGD_DA <- glmmTMB(TG_D ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGD_PxA <- glmmTMB(TG_D ~ Position * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGD_PxD <- glmmTMB(TG_D ~ Position * Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGD_DxA <- glmmTMB(TG_D ~ Day_Scaled * Age_Scaled + (1 | Field), family = binomial, data = FDModel) 

TGD_PAD <- glmmTMB(TG_D ~ Position + Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGD_PxAD <- glmmTMB(TG_D ~ Position * Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel)
TGD_PxDA <- glmmTMB(TG_D ~ Position * Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGD_PAxD <- glmmTMB(TG_D ~ Position + Day_Scaled * Age_Scaled + (1 | Field), family = binomial, data = FDModel) 

#collect models
TGD_modlist <- list("null"=TGD_null,"P"=TGD_P,"A"=TGD_A,"D"=TGD_D,
                    "PA" = TGD_PA, "PD" = TGD_PD, "DA" = TGD_DA,
                    "PxA" = TGD_PxA, "PxD" =TGD_PxD,"DxA" = TGD_DxA,
                    "PAD" = TGD_PAD,"PxAD" = TGD_PxAD, 
                    "PxDA" = TGD_PxDA, "PAxD" = TGD_PAxD)

aictab(TGD_modlist)
#Top model is Position * Age + Day

##Step 2: Environmental Variables----

head(FDModel)

TGD_Height <- glmmTMB(TG_D ~ Position * Crop_Age_Days + Day_Sampled +Height + (1 | Field), family = binomial, data = FDModel)
TGD_GC <- glmmTMB(TG_D ~ Position * Crop_Age_Days + Day_Sampled +GC + (1 | Field), family = binomial, data = FDModel)

TGD_FieldArea <- glmmTMB(TG_D ~ Position * Crop_Age_Days + Day_Sampled + Field_Size_Scaled + (1 | Field), family = binomial, data = FDModel)
TGD_NDVIfield <- glmmTMB(TG_D ~ Position * Crop_Age_Days + Day_Sampled +NDVImean_Field  + (1 | Field), family = binomial, data = FDModel)  

TGD_Rip <- glmmTMB(TG_D ~ Position * Crop_Age_Days + Day_Sampled + X1km_Rip_Prop + (1 | Field), family = binomial, data = FDModel) 
TGD_Crops <- glmmTMB(TG_D ~ Position * Crop_Age_Days + Day_Sampled +X1km_Prop_Crops + (1 | Field), family = binomial, data = FDModel) 
TGD_NDVI1km <- glmmTMB(TG_D ~ Position * Crop_Age_Days + Day_Sampled +NDVIsum_1km + (1 | Field), family = binomial, data = FDModel)


TGD_modlist2 <- list("null" = TGD_null,
                     "Height" = TGD_Height,
                     "GC" = TGD_GC,
                     "Field Size" = TGD_FieldArea,
                     "Field NDVI" = TGD_NDVIfield,
                     "Rip" = TGD_Rip,
                     "Crop" = TGD_Crops,
                     "NDVI" = TGD_NDVI1km,
                     "PxAD" = TGD_PxAD)
aictab(TGD_modlist2)
#Top model is Posistion * Age + Day (none within 2 AICc)

##Step 3: Check for Spatial Autocorrelation----

#running loon two models this for C so make sure to change over between running the loop
summary(TGD_PxAD)
summary(TGD_NDVIfield)
summary(TGD_GC)
summary(TGD_FieldArea)
summary(TGD_Height)
summary(TGD_Rip)

#run before loop
TGmodel_residuals <- simulateResiduals(TGD_Rip)
TGspatial_result <- data.frame(
  field = rep(NA, length(TGfield_numbers)),
  statistic = rep(NA, length(TGfield_numbers)),
  p_value = rep(NA, length(TGfield_numbers)),
  method = rep(NA_character_, length(TGfield_numbers)),
  stringsAsFactors = FALSE)

s <- 1

for (f in TGfield_numbers) {
  
  cat("Field", f, "\n") #What field is it doing?
  
  #Extracting specific residuals for individual fields
  TGfield_indices <- which(FDModel$ID == f)
  TGfield_residuals <- TGmodel_residuals
  TGfield_residuals$scaledResiduals <- 
    TGmodel_residuals$scaledResiduals[TGfield_indices]
  TGfield_residuals$fittedPredictedResponse <- 
    TGmodel_residuals$fittedPredictedResponse[TGfield_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  TGspatial_test <- testSpatialAutocorrelation(
    TGfield_residuals, 
    x = FDModel$X_Cor[FDModel$ID == f], 
    y = FDModel$Y_Cor[FDModel$ID == f])
  
  
  TGspatial_result$field [s] <- f
  TGspatial_result$statistic [s] <- TGspatial_test$statistic[1] 
  TGspatial_result$p_value [s] <- TGspatial_test$p.value
  TGspatial_result$method [s] <- TGspatial_test$method
  
  s <- s + 1
  
}

head(TGspatial_result);dim(TGspatial_result)
length(unique(FDModel$ID))

TGspatial_result

TG_D_Spatial_PxAD <- TGspatial_result #One field significant but spatial autocorrelation is negligable
TG_D_Spatial_FieldNDVI <- TGspatial_result #One field significant but spatial autocorrelation is negligable
TG_D_Spatial_GC <- TGspatial_result #One field significant but spatial autocorrelation is negligable
TG_D_Spatial_FieldSize <- TGspatial_result #One field significant but spatial autocorrelation is negligable
TG_D_Spatial_Height <- TGspatial_result #One field significant but spatial autocorrelation is negligable
TG_D_Spatial_Rip <- TGspatial_result #One field significant but spatial autocorrelation is negligable


write.xlsx(TG_D_Spatial_PxAD, 'TG_D_Spatial_PxAD.xlsx')
write.xlsx(TG_D_Spatial_FieldNDVI, 'TG_D_Spatial_FieldNDVI.xlsx')
write.xlsx(TG_D_Spatial_GC, 'TG_D_Spatial_GC.xlsx')
write.xlsx(TG_D_Spatial_FieldSize, 'TG_D_Spatial_FieldSize.xlsx')
write.xlsx(TG_D_Spatial_Height, 'TG_D_Spatial_Height.xlsx')
write.xlsx(TG_D_Spatial_Rip, 'TG_D_Spatial_Rip.xlsx')

##Step 4: Predictions----

#Position * Age + Day
summary(TGD_PxAD)

TG_D_pred <- expand.grid(Position  = unique(FDModel$Position),
                         Crop_Age_Days = TG_Predictions_Age,
                         Day_Sampled = TG_Predictions_Day)

head(TG_D_pred);dim(TG_D_pred)

TG_D_pred1 <- predict(object = TGD_PxAD,newdata= TG_D_pred,se.fit = T, type = "link",re.form = ~0)

TG_D_pred2<-data.frame(TG_D_pred,fit.link=TG_D_pred1$fit,se.link=TG_D_pred1$se.fit)

TG_D_pred2$lci.link<-TG_D_pred2$fit.link-(1.96*TG_D_pred2$se.link)
TG_D_pred2$uci.link<-TG_D_pred2$fit.link+(1.96*TG_D_pred2$se.link)

TG_D_pred2$fit<-plogis(TG_D_pred2$fit.link)
TG_D_pred2$se<-plogis(TG_D_pred2$se.link)
TG_D_pred2$lci<-plogis(TG_D_pred2$lci.link)
TG_D_pred2$uci<-plogis(TG_D_pred2$uci.link)

head(TG_D_pred2);dim(TG_D_pred2)

#NDVIfield
summary(TGD_NDVIfield)

TG_D_pred3 <- expand.grid(Position  = unique(FDModel$Position),
                         Crop_Age_Days = TG_Predictions_Age,
                         Day_Sampled = TG_Predictions_Day,
                         NDVImean_Field = TG_Predictions_NDVIField)

head(TG_D_pred3);dim(TG_D_pred3)

TG_D_pred4 <- predict(object = TGD_NDVIfield,newdata= TG_D_pred3,se.fit = T, type = "link",re.form = ~0)

TG_D_pred5<-data.frame(TG_D_pred3,fit.link=TG_D_pred4$fit,se.link=TG_D_pred4$se.fit)

TG_D_pred5$lci.link<-TG_D_pred5$fit.link-(1.96*TG_D_pred5$se.link)
TG_D_pred5$uci.link<-TG_D_pred5$fit.link+(1.96*TG_D_pred5$se.link)

TG_D_pred5$fit<-plogis(TG_D_pred5$fit.link)
TG_D_pred5$se<-plogis(TG_D_pred5$se.link)
TG_D_pred5$lci<-plogis(TG_D_pred5$lci.link)
TG_D_pred5$uci<-plogis(TG_D_pred5$uci.link)

head(TG_D_pred5);dim(TG_D_pred5)

#Ground cover
summary(TGD_GC)

TG_Predictions_GC <- seq(min(FDModel$GC),max(FDModel$GC),length.out=20)


TG_D_pred6 <- expand.grid(Position  = unique(FDModel$Position),
                         Crop_Age_Days = TG_Predictions_Age,
                         Day_Sampled = TG_Predictions_Day,
                         GC = TG_Predictions_GC)

head(TG_D_pred6);dim(TG_D_pred6)

TG_D_pred7 <- predict(object = TGD_GC,newdata= TG_D_pred6,se.fit = T, type = "link",re.form = ~0)

TG_D_pred8<-data.frame(TG_D_pred6,fit.link=TG_D_pred7$fit,se.link=TG_D_pred7$se.fit)

TG_D_pred8$lci.link<-TG_D_pred8$fit.link-(1.96*TG_D_pred8$se.link)
TG_D_pred8$uci.link<-TG_D_pred8$fit.link+(1.96*TG_D_pred8$se.link)

TG_D_pred8$fit<-plogis(TG_D_pred8$fit.link)
TG_D_pred8$se<-plogis(TG_D_pred8$se.link)
TG_D_pred8$lci<-plogis(TG_D_pred8$lci.link)
TG_D_pred8$uci<-plogis(TG_D_pred8$uci.link)

head(TG_D_pred8);dim(TG_D_pred8)

#Field Area
summary(TGD_FieldArea)
TG_D_pred9 <- expand.grid(Position  = unique(FDModel$Position),
                         Crop_Age_Days = TG_Predictions_Age,
                         Day_Sampled = TG_Predictions_Day,
                         Field_Size_Scaled = 
                           TG_Predictions_FieldSize_Scaled)

head(TG_D_pred9);dim(TG_D_pred9)

TG_D_pred10 <- predict(object = TGD_FieldArea,newdata= TG_D_pred9,se.fit = T, type = "link",re.form = ~0)

TG_D_pred11<-data.frame(TG_D_pred9,fit.link=TG_D_pred10$fit,se.link=TG_D_pred10$se.fit)

TG_D_pred11$lci.link<-TG_D_pred11$fit.link-(1.96*TG_D_pred11$se.link)
TG_D_pred11$uci.link<-TG_D_pred11$fit.link+(1.96*TG_D_pred11$se.link)

TG_D_pred11$fit<-plogis(TG_D_pred11$fit.link)
TG_D_pred11$se<-plogis(TG_D_pred11$se.link)
TG_D_pred11$lci<-plogis(TG_D_pred11$lci.link)
TG_D_pred11$uci<-plogis(TG_D_pred11$uci.link)

head(TG_D_pred11);dim(TG_D_pred11)


#Height
summary(TGD_Height)

TG_Predictions_Height <- seq(min(FDModel$Height),max(FDModel$Height),length.out=20)


TG_D_pred12 <- expand.grid(Position  = unique(FDModel$Position),
                         Crop_Age_Days = TG_Predictions_Age,
                         Day_Sampled = TG_Predictions_Day,
                         Height = TG_Predictions_Height)

head(TG_D_pred12);dim(TG_D_pred12)

TG_D_pred13 <- predict(object = TGD_Height,newdata= TG_D_pred12,se.fit = T, type = "link",re.form = ~0)

TG_D_pred14<-data.frame(TG_D_pred12,fit.link=TG_D_pred13$fit,se.link=TG_D_pred13$se.fit)

TG_D_pred14$lci.link<-TG_D_pred14$fit.link-(1.96*TG_D_pred14$se.link)
TG_D_pred14$uci.link<-TG_D_pred14$fit.link+(1.96*TG_D_pred14$se.link)

TG_D_pred14$fit<-plogis(TG_D_pred14$fit.link)
TG_D_pred14$se<-plogis(TG_D_pred14$se.link)
TG_D_pred14$lci<-plogis(TG_D_pred14$lci.link)
TG_D_pred14$uci<-plogis(TG_D_pred14$uci.link)

head(TG_D_pred14);dim(TG_D_pred14)

#Riparian
summary(TGD_Rip)

TG_D_pred15 <- expand.grid(Position  = unique(FDModel$Position),
                         Crop_Age_Days = TG_Predictions_Age,
                         Day_Sampled = TG_Predictions_Day,
                         X1km_Rip_Prop = TG_Predictions_Rip)

head(TG_D_pred15);dim(TG_D_pred15)

TG_D_pred16 <- predict(object = TGD_Rip,newdata= TG_D_pred15,se.fit = T, type = "link",re.form = ~0)

TG_D_pred17<-data.frame(TG_D_pred15,fit.link=TG_D_pred16$fit,se.link=TG_D_pred16$se.fit)

TG_D_pred17$lci.link<-TG_D_pred17$fit.link-(1.96*TG_D_pred17$se.link)
TG_D_pred17$uci.link<-TG_D_pred17$fit.link+(1.96*TG_D_pred17$se.link)

TG_D_pred17$fit<-plogis(TG_D_pred17$fit.link)
TG_D_pred17$se<-plogis(TG_D_pred17$se.link)
TG_D_pred17$lci<-plogis(TG_D_pred17$lci.link)
TG_D_pred17$uci<-plogis(TG_D_pred17$uci.link)

head(TG_D_pred17);dim(TG_D_pred17)


##Step 5: Visualisation----

#PxAD
summary(TGD_PxAD)
head(TG_D_pred2);dim(TG_D_pred2)

DD <- TG_D_pred2$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred2$Position == "Outer"
DDD <- TG_D_pred2$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred2$Position == "Inner"
D_D <- TG_D_pred2$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred2$Position == "Outer"



dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_G,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred2$Crop_Age_Days[DD],rev(TG_D_pred2$Crop_Age_Days[DD])), y = c(TG_D_pred2$lci[DD],rev(TG_D_pred2$uci[DD])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred2$Crop_Age_Days[DD],y = TG_D_pred2$fit[DD],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_D_pred2$Crop_Age_Days[DDD],rev(TG_D_pred2$Crop_Age_Days[DDD])), y = c(TG_D_pred2$lci[DDD],rev(TG_D_pred2$uci[DDD])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred2$Crop_Age_Days[DDD],y = TG_D_pred2$fit[DDD],lwd = 2,col = 'grey30',lty = 2)
legend('topright',legend = c('Outer', "Inner"), lty = c(1,2), col = 'grey30',pt.cex = 1)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred2$Day_Sampled[D_D],rev(TG_D_pred2$Day_Sampled[D_D])), y = c(TG_D_pred2$lci[D_D],rev(TG_D_pred2$uci[D_D])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred2$Day_Sampled[D_D],y = TG_D_pred2$fit[D_D],lwd = 2,col = 'grey30',lty = 1)

#NDVI field

summary(TGD_NDVIfield)
head(TG_D_pred5);dim(TG_D_pred5)

DD.2 <- TG_D_pred5$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred5$Position == "Outer" & TG_D_pred5$NDVImean_Field == TG_Predictions_NDVIField [10]
DDD.2 <- TG_D_pred5$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred5$Position == "Inner" & TG_D_pred5$NDVImean_Field == TG_Predictions_NDVIField [10]
D_D.2 <- TG_D_pred5$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred5$Position == "Outer" & TG_D_pred5$NDVImean_Field == TG_Predictions_NDVIField [10]
D_DD.2 <- TG_D_pred5$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred5$Position == "Outer" & TG_D_pred5$Day_Sampled == TG_Predictions_Day [10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_G,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred5$Crop_Age_Days[DD.2],rev(TG_D_pred5$Crop_Age_Days[DD.2])), y = c(TG_D_pred5$lci[DD.2],rev(TG_D_pred5$uci[DD.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred5$Crop_Age_Days[DD.2],y = TG_D_pred5$fit[DD.2],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_D_pred5$Crop_Age_Days[DDD.2],rev(TG_D_pred5$Crop_Age_Days[DDD.2])), y = c(TG_D_pred5$lci[DDD.2],rev(TG_D_pred5$uci[DDD.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred5$Crop_Age_Days[DDD.2],y = TG_D_pred5$fit[DDD.2],lwd = 2,col = 'grey30',lty = 2)
legend('bottomright',legend = c('Outer', "Inner"), lty = c(1,2), col = 'grey30',pt.cex = 1)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred5$Day_Sampled[D_D.2],rev(TG_D_pred5$Day_Sampled[D_D.2])), y = c(TG_D_pred5$lci[D_D.2],rev(TG_D_pred5$uci[D_D.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred5$Day_Sampled[D_D.2],y = TG_D_pred5$fit[D_D.2],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$NDVImean_Field,y = FDModel$TG_G,xlab = "Field NDVI",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred5$NDVImean_Field[D_DD.2],rev(TG_D_pred5$NDVImean_Field[D_DD.2])), y = c(TG_D_pred5$lci[D_DD.2],rev(TG_D_pred5$uci[D_DD.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred5$NDVImean_Field[D_DD.2],y = TG_D_pred5$fit[D_DD.2],lwd = 2,col = 'grey30',lty = 1)

#Ground cover
summary(TGD_GC)
head(TG_D_pred8);dim(TG_D_pred8)

DD.3 <- TG_D_pred8$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred8$Position == "Outer" & TG_D_pred8$GC == TG_Predictions_GC [10]
DDD.3 <- TG_D_pred8$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred8$Position == "Inner" & TG_D_pred8$GC == TG_Predictions_GC [10]
D_D.3 <- TG_D_pred8$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred8$Position == "Outer" & TG_D_pred8$GC == TG_Predictions_GC [10]
D_DD.3 <- TG_D_pred8$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred8$Position == "Outer" & TG_D_pred8$Day_Sampled == TG_Predictions_Day [10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_G,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred8$Crop_Age_Days[DD.3],rev(TG_D_pred8$Crop_Age_Days[DD.3])), y = c(TG_D_pred8$lci[DD.3],rev(TG_D_pred8$uci[DD.3])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred8$Crop_Age_Days[DD.3],y = TG_D_pred8$fit[DD.3],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_D_pred8$Crop_Age_Days[DDD.3],rev(TG_D_pred8$Crop_Age_Days[DDD.3])), y = c(TG_D_pred8$lci[DDD.3],rev(TG_D_pred8$uci[DDD.3])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred8$Crop_Age_Days[DDD.3],y = TG_D_pred8$fit[DDD.3],lwd = 2,col = 'grey30',lty = 2)
legend('bottomright',legend = c('Outer', "Inner"), lty = c(1,2), col = 'grey30',pt.cex = 1)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred8$Day_Sampled[D_D.3],rev(TG_D_pred8$Day_Sampled[D_D.3])), y = c(TG_D_pred8$lci[D_D.3],rev(TG_D_pred8$uci[D_D.3])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred8$Day_Sampled[D_D.3],y = TG_D_pred8$fit[D_D.3],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$GC,y = FDModel$TG_G,xlab = "Ground Cover (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred8$GC[D_DD.3],rev(TG_D_pred8$GC[D_DD.3])), y = c(TG_D_pred8$lci[D_DD.3],rev(TG_D_pred8$uci[D_DD.3])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred8$GC[D_DD.3],y = TG_D_pred8$fit[D_DD.3],lwd = 2,col = 'grey30',lty = 1)

#Field Size

summary(TGD_FieldArea)
head(TG_D_pred11);dim(TG_D_pred11)

DD.4 <- TG_D_pred11$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred11$Position == "Outer" & TG_D_pred11$Field_Size_Scaled == TG_Predictions_FieldSize_Scaled [10]
DDD.4 <- TG_D_pred11$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred11$Position == "Inner" & TG_D_pred11$Field_Size_Scaled == TG_Predictions_FieldSize_Scaled [10]
D_D.4 <- TG_D_pred11$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred11$Position == "Outer" & TG_D_pred11$Field_Size_Scaled == TG_Predictions_FieldSize_Scaled [10]
D_DD.4 <- TG_D_pred11$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred11$Position == "Outer" & TG_D_pred11$Day_Sampled == TG_Predictions_Day [10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_G,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred11$Crop_Age_Days[DD.4],rev(TG_D_pred11$Crop_Age_Days[DD.4])), y = c(TG_D_pred11$lci[DD.4],rev(TG_D_pred11$uci[DD.4])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred11$Crop_Age_Days[DD.4],y = TG_D_pred11$fit[DD.4],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_D_pred11$Crop_Age_Days[DDD.4],rev(TG_D_pred11$Crop_Age_Days[DDD.4])), y = c(TG_D_pred11$lci[DDD.4],rev(TG_D_pred11$uci[DDD.4])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred11$Crop_Age_Days[DDD.4],y = TG_D_pred11$fit[DDD.4],lwd = 2,col = 'grey30',lty = 2)
legend('bottomright',legend = c('Outer', "Inner"), lty = c(1,2), col = 'grey30',pt.cex = 1)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred11$Day_Sampled[D_D.4],rev(TG_D_pred11$Day_Sampled[D_D.4])), y = c(TG_D_pred11$lci[D_D.4],rev(TG_D_pred11$uci[D_D.4])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred11$Day_Sampled[D_D.4],y = TG_D_pred11$fit[D_D.4],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$Field_Size_Scaled,y = FDModel$TG_G,xlab = "Field Size (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_D_pred11$Field_Size_Scaled),to=max(TG_D_pred11$Field_Size_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1),cex.axis=1)

polygon(x = c(TG_D_pred11$Field_Size_Scaled[D_DD.4],rev(TG_D_pred11$Field_Size_Scaled[D_DD.4])), y = c(TG_D_pred11$lci[D_DD.4],rev(TG_D_pred11$uci[D_DD.4])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred11$Field_Size_Scaled[D_DD.4],y = TG_D_pred11$fit[D_DD.4],lwd = 2,col = 'grey30',lty = 1)

#Height
summary(TGD_Height)
head(TG_D_pred14);dim(TG_D_pred14)

DD.5 <- TG_D_pred14$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred14$Position == "Outer" & TG_D_pred14$Height == TG_Predictions_Height [10]
DDD.5 <- TG_D_pred14$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred14$Position == "Inner" & TG_D_pred14$Height == TG_Predictions_Height [10]
D_D.5 <- TG_D_pred14$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred14$Position == "Outer" & TG_D_pred14$Height == TG_Predictions_Height [10]
D_DD.5 <- TG_D_pred14$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred14$Position == "Outer" & TG_D_pred14$Day_Sampled == TG_Predictions_Day [10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_G,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred14$Crop_Age_Days[DD.5],rev(TG_D_pred14$Crop_Age_Days[DD.5])), y = c(TG_D_pred14$lci[DD.5],rev(TG_D_pred14$uci[DD.5])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred14$Crop_Age_Days[DD.5],y = TG_D_pred14$fit[DD.5],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_D_pred14$Crop_Age_Days[DDD.5],rev(TG_D_pred14$Crop_Age_Days[DDD.5])), y = c(TG_D_pred14$lci[DDD.5],rev(TG_D_pred14$uci[DDD.5])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred14$Crop_Age_Days[DDD.5],y = TG_D_pred14$fit[DDD.5],lwd = 2,col = 'grey30',lty = 2)
legend('bottomright',legend = c('Outer', "Inner"), lty = c(1,2), col = 'grey30',pt.cex = 1)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred14$Day_Sampled[D_D.5],rev(TG_D_pred14$Day_Sampled[D_D.5])), y = c(TG_D_pred14$lci[D_D.5],rev(TG_D_pred14$uci[D_D.5])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred14$Day_Sampled[D_D.5],y = TG_D_pred14$fit[D_D.5],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$Height,y = FDModel$TG_G,xlab = "Crop Height (cm)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred14$Height[D_DD.5],rev(TG_D_pred14$Height[D_DD.5])), y = c(TG_D_pred14$lci[D_DD.5],rev(TG_D_pred14$uci[D_DD.5])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred14$Height[D_DD.5],y = TG_D_pred14$fit[D_DD.5],lwd = 2,col = 'grey30',lty = 1)

#Riparian
summary(TGD_Rip)
head(TG_D_pred17);dim(TG_D_pred17)

DD.6 <- TG_D_pred17$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred17$Position == "Outer" & TG_D_pred17$X1km_Rip_Prop == TG_Predictions_Rip [10]
DDD.6 <- TG_D_pred17$Day_Sampled == TG_Predictions_Day [10] & TG_D_pred17$Position == "Inner" & TG_D_pred17$X1km_Rip_Prop == TG_Predictions_Rip [10]
D_D.6 <- TG_D_pred17$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred17$Position == "Outer" & TG_D_pred17$X1km_Rip_Prop == TG_Predictions_Rip [10]
D_DD.6 <- TG_D_pred17$Crop_Age_Days == TG_Predictions_Age [10] & TG_D_pred17$Position == "Outer" & TG_D_pred17$Day_Sampled == TG_Predictions_Day [10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_G,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred17$Crop_Age_Days[DD.6],rev(TG_D_pred17$Crop_Age_Days[DD.6])), y = c(TG_D_pred17$lci[DD.6],rev(TG_D_pred17$uci[DD.6])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred17$Crop_Age_Days[DD.6],y = TG_D_pred17$fit[DD.6],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_D_pred17$Crop_Age_Days[DDD.6],rev(TG_D_pred17$Crop_Age_Days[DDD.6])), y = c(TG_D_pred17$lci[DDD.6],rev(TG_D_pred17$uci[DDD.6])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred17$Crop_Age_Days[DDD.6],y = TG_D_pred17$fit[DDD.6],lwd = 2,col = 'grey30',lty = 2)
legend('bottomright',legend = c('Outer', "Inner"), lty = c(1,2), col = 'grey30',pt.cex = 1)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred17$Day_Sampled[D_D.6],rev(TG_D_pred17$Day_Sampled[D_D.6])), y = c(TG_D_pred17$lci[D_D.6],rev(TG_D_pred17$uci[D_D.6])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred17$Day_Sampled[D_D.6],y = TG_D_pred17$fit[D_D.6],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$X1km_Rip_Prop,y = FDModel$TG_G,xlab = "Riparian within 1km (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_D_pred17$X1km_Rip_Prop[D_DD.6],rev(TG_D_pred17$X1km_Rip_Prop[D_DD.6])), y = c(TG_D_pred17$lci[D_DD.6],rev(TG_D_pred17$uci[D_DD.6])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred17$X1km_Rip_Prop[D_DD.6],y = TG_D_pred17$fit[D_DD.6],lwd = 2,col = 'grey30',lty = 1)


#Trait Group E----
##Step 1: Modelling Design Variables----

head(FDModel);dim(FDModel)
str(FDModel)

TGE_null <- glmmTMB(TG_E  ~ 1 + (1 | Field), family = binomial, data = FDModel)

TGE_P <- glmmTMB(TG_E ~ Position + (1 | Field), family = binomial, data = FDModel)
TGE_A <- glmmTMB(TG_E ~ Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGE_D <- glmmTMB(TG_E ~ Day_Sampled + (1 | Field), family = binomial, data = FDModel) 

TGE_PA <- glmmTMB(TG_E ~ Position + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGE_PD <- glmmTMB(TG_E ~ Position + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGE_DA <- glmmTMB(TG_E ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGE_PxA <- glmmTMB(TG_E ~ Position * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGE_PxD <- glmmTMB(TG_E ~ Position * Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGE_DxA <- glmmTMB(TG_E ~ Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGE_PAD <- glmmTMB(TG_E ~ Position + Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGE_PxAD <- glmmTMB(TG_E ~ Position * Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel)
TGE_PxDA <- glmmTMB(TG_E ~ Position * Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGE_PAxD <- glmmTMB(TG_E ~ Position + Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

#collect models
TGE_modlist <- list("null" = TGE_null, "P" = TGE_P,"A" = TGE_A,
                    "D" = TGE_D, "PA" = TGE_PA, "PD" = TGE_PD, 
                    "DA" = TGE_DA, "PxA" = TGE_PxA, 
                    "PxD" = TGE_PxD,"DxA" = TGE_DxA,
                    "PAD" = TGE_PAD,"PxAD" = TGE_PxAD, 
                    "PxDA" = TGE_PxDA, "PAxD" = TGE_PAxD)

aictab(TGE_modlist)
#Top model is Day * Age

##Step 2: Environmental Variables----

head(FDModel)

TGE_Height <- glmmTMB(TG_E ~ Day_Sampled * Crop_Age_Days + Height + (1 | Field), family = binomial, data = FDModel)
TGE_GC <- glmmTMB(TG_E ~ Day_Sampled * Crop_Age_Days + GC + (1 | Field), family = binomial, data = FDModel)

TGE_FieldArea <- glmmTMB(TG_E ~ Day_Sampled * Crop_Age_Days + Field_Area_m2 + (1 | Field), family = binomial, data = FDModel)
TGE_NDVIfield <- glmmTMB(TG_E ~ Day_Sampled * Crop_Age_Days + NDVImean_Field  + (1 | Field), family = binomial, data = FDModel)  

TGE_Rip <- glmmTMB(TG_E ~ Day_Sampled * Crop_Age_Days + X1km_Rip_Prop + (1 | Field), family = binomial, data = FDModel) 
TGE_Crops <- glmmTMB(TG_E ~ Day_Sampled * Crop_Age_Days + X1km_Prop_Crops + (1 | Field), family = binomial, data = FDModel) 
TGE_NDVI1km <- glmmTMB(TG_E ~ Day_Sampled * Crop_Age_Days + NDVIsum_1km + (1 | Field), family = binomial, data = FDModel)


TGE_modlist2 <- list("null" = TGE_null,
                     "Height" = TGE_Height,
                     "GC" = TGE_GC,
                     "Field Size" = TGE_FieldArea,
                     "Field NDVI" = TGE_NDVIfield,
                     "Rip" = TGE_Rip,
                     "Crop" = TGE_Crops,
                     "NDVI 1km" = TGE_NDVI1km,
                     "DxA" = TGE_DxA)
aictab(TGE_modlist2)
#Top model is NDVI 1km (None within 2 AICc's)

##Step 3: Check for Spatial Autocorrelation----

summary(TGE_NDVI1km)


#run before loop
TGmodel_residuals <- simulateResiduals(TGE_NDVI1km)
TGspatial_result <- data.frame(
  field = rep(NA, length(TGfield_numbers)),
  statistic = rep(NA, length(TGfield_numbers)),
  p_value = rep(NA, length(TGfield_numbers)),
  method = rep(NA_character_, length(TGfield_numbers)),
  stringsAsFactors = FALSE)

s <- 1

for (f in TGfield_numbers) {
  
  cat("Field", f, "\n") #What field is it doing?
  
  #Extracting specific residuals for individual fields
  TGfield_indices <- which(FDModel$ID == f)
  TGfield_residuals <- TGmodel_residuals
  TGfield_residuals$scaledResiduals <- 
    TGmodel_residuals$scaledResiduals[TGfield_indices]
  TGfield_residuals$fittedPredictedResponse <- 
    TGmodel_residuals$fittedPredictedResponse[TGfield_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  TGspatial_test <- testSpatialAutocorrelation(
    TGfield_residuals, 
    x = FDModel$X_Cor[FDModel$ID == f], 
    y = FDModel$Y_Cor[FDModel$ID == f])
  
  
  TGspatial_result$field [s] <- f
  TGspatial_result$statistic [s] <- TGspatial_test$statistic[1] 
  TGspatial_result$p_value [s] <- TGspatial_test$p.value
  TGspatial_result$method [s] <- TGspatial_test$method
  
  s <- s + 1
  
}

head(TGspatial_result);dim(TGspatial_result)
length(unique(FDModel$ID))

TGspatial_result


TG_E_Spatial <- TGspatial_result
#No spatial autocorrelation found 

write.xlsx(TG_E_Spatial, 'TG_E_Spatial.xlsx')

##Step 4: Predictions----


summary(TGE_NDVI1km)


TG_Predictions_NDVI1km <- seq(min(FDModel$NDVIsum_1km),max(FDModel$NDVIsum_1km),length.out=20)


TG_E_pred <- expand.grid(Day_Sampled = TG_Predictions_Day,
                         Crop_Age_Days = TG_Predictions_Age,
                         NDVIsum_1km = TG_Predictions_NDVI1km)
head(TG_E_pred);dim(TG_E_pred)

TG_E_pred1 <- predict(object = TGE_NDVI1km,newdata= TG_E_pred,se.fit = T, type = "link",re.form = ~0)

TG_E_pred2<-data.frame(TG_E_pred,fit.link=TG_E_pred1$fit,se.link=TG_E_pred1$se.fit)

TG_E_pred2$lci.link<-TG_E_pred2$fit.link-(1.96*TG_E_pred2$se.link)
TG_E_pred2$uci.link<-TG_E_pred2$fit.link+(1.96*TG_E_pred2$se.link)

TG_E_pred2$fit<-plogis(TG_E_pred2$fit.link)
TG_E_pred2$se<-plogis(TG_E_pred2$se.link)
TG_E_pred2$lci<-plogis(TG_E_pred2$lci.link)
TG_E_pred2$uci<-plogis(TG_E_pred2$uci.link)

head(TG_E_pred2);dim(TG_E_pred2)

##Step 5: Visualisation----
TG_Predictions_Day [4] #Winter
TG_Predictions_Day [15] #Spring

summary(TGE_NDVI1km)
head(TG_E_pred2);dim(TG_E_pred2)

EE <- TG_E_pred2$Day_Sampled == TG_Predictions_Day [4] & TG_E_pred2$NDVIsum_1km == TG_Predictions_NDVI1km[10]
EEE <- TG_E_pred2$Day_Sampled == TG_Predictions_Day [15] & TG_E_pred2$NDVIsum_1km == TG_Predictions_NDVI1km[10]
E_E <- TG_E_pred2$Day_Sampled == TG_Predictions_Day [10] & TG_E_pred2$Crop_Age_Days == TG_Predictions_Age[10]



dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_E,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_E_pred2$Crop_Age_Days[EE],rev(TG_E_pred2$Crop_Age_Days[EE])), y = c(TG_E_pred2$lci[EE],rev(TG_E_pred2$uci[EE])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_E_pred2$Crop_Age_Days[EE],y = TG_E_pred2$fit[EE],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_E_pred2$Crop_Age_Days[EEE],rev(TG_E_pred2$Crop_Age_Days[EEE])), y = c(TG_E_pred2$lci[EEE],rev(TG_E_pred2$uci[EEE])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_E_pred2$Crop_Age_Days[EEE],y = TG_E_pred2$fit[EEE],lwd = 2,col = 'grey30',lty = 2)
legend('topright',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

plot(x = FDModel$NDVIsum_1km,y = FDModel$TG_E,xlab = "NDVI within 1km",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_E_pred2$NDVIsum_1km[E_E],rev(TG_E_pred2$NDVIsum_1km[E_E])), y = c(TG_E_pred2$lci[E_E],rev(TG_E_pred2$uci[E_E])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_E_pred2$NDVIsum_1km[E_E],y = TG_E_pred2$fit[E_E],lwd = 2,col = 'grey30',lty = 1)

#Trait Group F----
##Step 1: Modelling Design Variables----

head(FDModel);dim(FDModel)
str(FDModel)

TGF_null <- glmmTMB(TG_F  ~ 1 + (1 | Field), family = binomial, data = FDModel)

TGF_P <- glmmTMB(TG_F ~ Position + (1 | Field), family = binomial, data = FDModel)
TGF_A <- glmmTMB(TG_F ~ Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGF_D <- glmmTMB(TG_F ~ Day_Sampled + (1 | Field), family = binomial, data = FDModel) 

TGF_PA <- glmmTMB(TG_F ~ Position + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGF_PD <- glmmTMB(TG_F ~ Position + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGF_DA <- glmmTMB(TG_F ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGF_PxA <- glmmTMB(TG_F ~ Position * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGF_PxD <- glmmTMB(TG_F ~ Position * Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGF_DxA <- glmmTMB(TG_F ~ Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGF_PAD <- glmmTMB(TG_F ~ Position + Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGF_PxAD <- glmmTMB(TG_F ~ Position * Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel)
TGF_PxDA <- glmmTMB(TG_F ~ Position * Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGF_PAxD <- glmmTMB(TG_F ~ Position + Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

#collect models
TGF_modlist <- list("null" = TGF_null, "P" = TGF_P,"A" = TGF_A,
                    "D" = TGF_D, "PA" = TGF_PA, "PD" = TGF_PD, 
                    "DA" = TGF_DA, "PxA" = TGF_PxA, 
                    "PxD" =TGF_PxD,"DxA" = TGF_DxA,
                    "PAD" = TGF_PAD,"PxAD" = TGF_PxAD, 
                    "PxDA" = TGF_PxDA, "PAxD" = TGF_PAxD)

aictab(TGF_modlist)
#Top model is null

##Step 2: Environmental Variables----

head(FDModel)

TGF_Height <- glmmTMB(TG_F ~ Height + (1 | Field), family = binomial, data = FDModel)
TGF_GC <- glmmTMB(TG_F ~ GC + (1 | Field), family = binomial, data = FDModel)

TGF_FieldArea <- glmmTMB(TG_F ~ Field_Size_Scaled + (1 | Field), family = binomial, data = FDModel)
TGF_NDVIfield <- glmmTMB(TG_F ~ NDVImean_Field  + (1 | Field), family = binomial, data = FDModel)  

TGF_Rip <- glmmTMB(TG_F ~ X1km_Rip_Prop + (1 | Field), family = binomial, data = FDModel) 
TGF_Crops <- glmmTMB(TG_F ~ X1km_Prop_Crops + (1 | Field), family = binomial, data = FDModel) 
TGF_NDVI1km <- glmmTMB(TG_F ~ NDVIsum_1km + (1 | Field), family = binomial, data = FDModel)


TGF_modlist2 <- list("null" = TGF_null,
                     "Height" = TGF_Height,
                     "GC" = TGF_GC,
                     "Field Size" = TGF_FieldArea,
                     "Field NDVI" = TGF_NDVIfield,
                     "Rip" = TGF_Rip,
                     "Crop" = TGF_Crops,
                     "NDVI 1km" = TGF_NDVI1km)
aictab(TGF_modlist2)
#Top model is GC (None within 2 AICc's)


##Step 3: Check for Spatial Autocorrelation----

summary(TGF_GC)


#run before loop
TGmodel_residuals <- simulateResiduals(TGF_GC)
TGspatial_result <- data.frame(
  field = rep(NA, length(TGfield_numbers)),
  statistic = rep(NA, length(TGfield_numbers)),
  p_value = rep(NA, length(TGfield_numbers)),
  method = rep(NA_character_, length(TGfield_numbers)),
  stringsAsFactors = FALSE)

s <- 1

for (f in TGfield_numbers) {
  
  cat("Field", f, "\n") #What field is it doing?
  
  #Extracting specific residuals for individual fields
  TGfield_indices <- which(FDModel$ID == f)
  TGfield_residuals <- TGmodel_residuals
  TGfield_residuals$scaledResiduals <- 
    TGmodel_residuals$scaledResiduals[TGfield_indices]
  TGfield_residuals$fittedPredictedResponse <- 
    TGmodel_residuals$fittedPredictedResponse[TGfield_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  TGspatial_test <- testSpatialAutocorrelation(
    TGfield_residuals, 
    x = FDModel$X_Cor[FDModel$ID == f], 
    y = FDModel$Y_Cor[FDModel$ID == f])
  
  
  TGspatial_result$field [s] <- f
  TGspatial_result$statistic [s] <- TGspatial_test$statistic[1] 
  TGspatial_result$p_value [s] <- TGspatial_test$p.value
  TGspatial_result$method [s] <- TGspatial_test$method
  
  s <- s + 1
  
}

head(TGspatial_result);dim(TGspatial_result)
length(unique(FDModel$ID))

TGspatial_result


TG_F_Spatial <- TGspatial_result
#One field significant but spatial autocorrelation is negligable


write.xlsx(TG_F_Spatial, 'TG_F_Spatial.xlsx')


##Step 4: Predictions----


summary(TGF_GC)


TG_F_pred <- data.frame(GC = TG_Predictions_GC)
head(TG_F_pred);dim(TG_F_pred)

TG_F_pred1 <- predict(object = TGF_GC,newdata= TG_F_pred,se.fit = T, type = "link",re.form = ~0)

TG_F_pred2<-data.frame(TG_F_pred,fit.link=TG_F_pred1$fit,se.link=TG_F_pred1$se.fit)

TG_F_pred2$lci.link<-TG_F_pred2$fit.link-(1.96*TG_F_pred2$se.link)
TG_F_pred2$uci.link<-TG_F_pred2$fit.link+(1.96*TG_F_pred2$se.link)

TG_F_pred2$fit<-plogis(TG_F_pred2$fit.link)
TG_F_pred2$se<-plogis(TG_F_pred2$se.link)
TG_F_pred2$lci<-plogis(TG_F_pred2$lci.link)
TG_F_pred2$uci<-plogis(TG_F_pred2$uci.link)

head(TG_F_pred2);dim(TG_F_pred2)

##Step 5: Visualisation----

summary(TGF_GC)
head(TG_F_pred2);dim(TG_F_pred2)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$GC,y = FDModel$TG_F,xlab = "Ground Cover (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_F_pred2$GC,rev(TG_F_pred2$GC)), y = c(TG_F_pred2$lci,rev(TG_F_pred2$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_F_pred2$GC,y = TG_F_pred2$fit,lwd = 2,col = 'grey30',lty = 1)



#Trait Group G----
##Step 1: Modelling Design Variables----

head(FDModel);dim(FDModel)
str(FDModel)

TGG_null <- glmmTMB(TG_G  ~ 1 + (1 | Field), family = binomial, data = FDModel)

TGG_P <- glmmTMB(TG_G ~ Position + (1 | Field), family = binomial, data = FDModel)
TGG_A <- glmmTMB(TG_G ~ Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGG_D <- glmmTMB(TG_G ~ Day_Sampled + (1 | Field), family = binomial, data = FDModel) 

TGG_PA <- glmmTMB(TG_G ~ Position + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGG_PD <- glmmTMB(TG_G ~ Position + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGG_DA <- glmmTMB(TG_G ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGG_PxA <- glmmTMB(TG_G ~ Position * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 
TGG_PxD <- glmmTMB(TG_G ~ Position * Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGG_DxA <- glmmTMB(TG_G ~ Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

TGG_PAD <- glmmTMB(TG_G ~ Position + Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel) 
TGG_PxAD <- glmmTMB(TG_G ~ Position * Crop_Age_Days + Day_Sampled + (1 | Field), family = binomial, data = FDModel)
TGG_PxDA <- glmmTMB(TG_G ~ Position * Day_Sampled + Crop_Age_Days + (1 | Field), family = binomial, data = FDModel)
TGG_PAxD <- glmmTMB(TG_G ~ Position + Day_Sampled * Crop_Age_Days + (1 | Field), family = binomial, data = FDModel) 

#collect models
TGG_modlist <- list("null" = TGG_null, "P" = TGG_P,"A" = TGG_A,
                    "D" = TGG_D, "PA" = TGG_PA, "PD" = TGG_PD, 
                    "DA" = TGG_DA, "PxA" = TGG_PxA, 
                    "PxD" = TGG_PxD,"DxA" = TGG_DxA,
                    "PAD" = TGG_PAD,"PxAD" = TGG_PxAD, 
                    "PxDA" = TGG_PxDA, "PAxD" = TGG_PAxD)

aictab(TGG_modlist)
#Top model is Day

##Step 2: Environmental Variables----

head(FDModel)

TGG_Height <- glmmTMB(TG_G ~ Day_Sampled + Height + (1 | Field), family = binomial, data = FDModel)
TGG_GC <- glmmTMB(TG_G ~ Day_Sampled + GC + (1 | Field), family = binomial, data = FDModel)

TGG_FieldArea <- glmmTMB(TG_G ~ Day_Sampled + Field_Size_Scaled + (1 | Field), family = binomial, data = FDModel)
TGG_NDVIfield <- glmmTMB(TG_G ~ Day_Sampled + NDVImean_Field  + (1 | Field), family = binomial, data = FDModel)  

TGG_Rip <- glmmTMB(TG_G ~ Day_Sampled + X1km_Rip_Prop + (1 | Field), family = binomial, data = FDModel) 
TGG_Crops <- glmmTMB(TG_G ~ Day_Sampled + X1km_Prop_Crops + (1 | Field), family = binomial, data = FDModel) 
TGG_NDVI1km <- glmmTMB(TG_G ~ Day_Sampled + NDVIsum_1km + (1 | Field), family = binomial, data = FDModel)


TGG_modlist2 <- list("null" = TGG_null,
                     "Height" = TGG_Height,
                     "GC" = TGG_GC,
                     "Field Size" = TGG_FieldArea,
                     "Field NDVI" = TGG_NDVIfield,
                     "Rip" = TGG_Rip,
                     "Crop" = TGG_Crops,
                     "NDVI" = TGG_NDVI1km,
                     "D" = TGG_D)
aictab(TGG_modlist2)
11.41#Top model is Crop (3 within 2 AICc's)
#Day is within AICc's but crop is not considered uninformative because the models LL is improved by it 


##Step 3: Check for Spatial Autocorrelation----

summary(TGG_Crops)
summary(TGG_D)
summary(TGG_NDVIfield)
summary(TGG_FieldArea)
summary(TGG_GC)



#run before loop
TGmodel_residuals <- simulateResiduals(TGG_D)
TGspatial_result <- data.frame(
  field = rep(NA, length(TGfield_numbers)),
  statistic = rep(NA, length(TGfield_numbers)),
  p_value = rep(NA, length(TGfield_numbers)),
  method = rep(NA_character_, length(TGfield_numbers)),
  stringsAsFactors = FALSE)

s <- 1

for (f in TGfield_numbers) {
  
  cat("Field", f, "\n") #What field is it doing?
  
  #Extracting specific residuals for individual fields
  TGfield_indices <- which(FDModel$ID == f)
  TGfield_residuals <- TGmodel_residuals
  TGfield_residuals$scaledResiduals <- 
    TGmodel_residuals$scaledResiduals[TGfield_indices]
  TGfield_residuals$fittedPredictedResponse <- 
    TGmodel_residuals$fittedPredictedResponse[TGfield_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  TGspatial_test <- testSpatialAutocorrelation(
    TGfield_residuals, 
    x = FDModel$X_Cor[FDModel$ID == f], 
    y = FDModel$Y_Cor[FDModel$ID == f])
  
  
  TGspatial_result$field [s] <- f
  TGspatial_result$statistic [s] <- TGspatial_test$statistic[1] 
  TGspatial_result$p_value [s] <- TGspatial_test$p.value
  TGspatial_result$method [s] <- TGspatial_test$method
  
  s <- s + 1
  
}

head(TGspatial_result);dim(TGspatial_result)
length(unique(FDModel$ID))

TGspatial_result


TG_G_Spatial_Crops <- TGspatial_result
#One field significant but spatial autocorrelation is negligable
TG_G_Spatial_FieldNDVI <- TGspatial_result #One field significant but spatial autocorrelation is negligable
TG_G_Spatial_FieldSize <- TGspatial_result #Two field significant but spatial autocorrelation is negligable
TG_G_Spatial_GC <- TGspatial_result #One field significant but spatial autocorrelation is negligable
TG_G_Spatial_D <- TGspatial_result #One field significant but spatial autocorrelation is negligable



write.xlsx(TG_G_Spatial_Crops, 'TG_G_Spatial_Crops.xlsx')
write.xlsx(TG_G_Spatial_FieldNDVI, 'TG_G_Spatial_FieldNDVI.xlsx')
write.xlsx(TG_G_Spatial_FieldSize, 'TG_G_Spatial_FieldSize.xlsx')
write.xlsx(TG_G_Spatial_GC, 'TG_G_Spatial_GC.xlsx')
write.xlsx(TG_G_Spatial_D, 'TG_G_Spatial_D.xlsx')


##Step 4: Predictions----

#crops
summary(TGG_Crops)


TG_Predictions_Crops <- seq(min(FDModel$X1km_Prop_Crops),max(FDModel$X1km_Prop_Crops),length.out=20)


TG_G_pred <- expand.grid(Day_Sampled = TG_Predictions_Day,
                         X1km_Prop_Crops = TG_Predictions_Crops)
head(TG_G_pred);dim(TG_G_pred)

TG_G_pred1 <- predict(object = TGG_Crops,newdata= TG_G_pred,se.fit = T, type = "link",re.form = ~0)

TG_G_pred2<-data.frame(TG_G_pred,fit.link=TG_G_pred1$fit,se.link=TG_G_pred1$se.fit)

TG_G_pred2$lci.link<-TG_G_pred2$fit.link-(1.96*TG_G_pred2$se.link)
TG_G_pred2$uci.link<-TG_G_pred2$fit.link+(1.96*TG_G_pred2$se.link)

TG_G_pred2$fit<-plogis(TG_G_pred2$fit.link)
TG_G_pred2$se<-plogis(TG_G_pred2$se.link)
TG_G_pred2$lci<-plogis(TG_G_pred2$lci.link)
TG_G_pred2$uci<-plogis(TG_G_pred2$uci.link)

head(TG_G_pred2);dim(TG_G_pred2)

#NDVI field
summary(TGG_NDVIfield)

TG_Predictions_NDVIField <- seq(min(FDModel$NDVImean_Field ),max(FDModel$NDVImean_Field ),length.out=20)


TG_G_pred3 <- expand.grid(Day_Sampled = TG_Predictions_Day,
                         NDVImean_Field = TG_Predictions_NDVIField)
head(TG_G_pred3);dim(TG_G_pred3)

TG_G_pred4 <- predict(object = TGG_NDVIfield,newdata= TG_G_pred3,se.fit = T, type = "link",re.form = ~0)

TG_G_pred5<-data.frame(TG_G_pred3,fit.link=TG_G_pred4$fit,se.link=TG_G_pred4$se.fit)

TG_G_pred5$lci.link<-TG_G_pred5$fit.link-(1.96*TG_G_pred5$se.link)
TG_G_pred5$uci.link<-TG_G_pred5$fit.link+(1.96*TG_G_pred5$se.link)

TG_G_pred5$fit<-plogis(TG_G_pred5$fit.link)
TG_G_pred5$se<-plogis(TG_G_pred5$se.link)
TG_G_pred5$lci<-plogis(TG_G_pred5$lci.link)
TG_G_pred5$uci<-plogis(TG_G_pred5$uci.link)

head(TG_G_pred5);dim(TG_G_pred5)

#Field Area
summary(TGG_FieldArea)

TG_Predictions_NDVIField <- seq(min(FDModel$NDVImean_Field ),max(FDModel$NDVImean_Field ),length.out=20)


TG_G_pred6 <- expand.grid(Day_Sampled = TG_Predictions_Day,
                          Field_Size_Scaled =
                            TG_Predictions_FieldSize_Scaled)
head(TG_G_pred6);dim(TG_G_pred6)

TG_G_pred7 <- predict(object = TGG_FieldArea,newdata= TG_G_pred6,se.fit = T, type = "link",re.form = ~0)

TG_G_pred8<-data.frame(TG_G_pred6,fit.link=TG_G_pred7$fit,se.link=TG_G_pred7$se.fit)

TG_G_pred8$lci.link<-TG_G_pred8$fit.link-(1.96*TG_G_pred8$se.link)
TG_G_pred8$uci.link<-TG_G_pred8$fit.link+(1.96*TG_G_pred8$se.link)

TG_G_pred8$fit<-plogis(TG_G_pred8$fit.link)
TG_G_pred8$se<-plogis(TG_G_pred8$se.link)
TG_G_pred8$lci<-plogis(TG_G_pred8$lci.link)
TG_G_pred8$uci<-plogis(TG_G_pred8$uci.link)

head(TG_G_pred8);dim(TG_G_pred8)

#Ground Cover
summary(TGG_GC)

TG_Predictions_NDVIField <- seq(min(FDModel$NDVImean_Field ),max(FDModel$NDVImean_Field ),length.out=20)


TG_G_pred9 <- expand.grid(Day_Sampled = TG_Predictions_Day,
                          GC = TG_Predictions_GC)
head(TG_G_pred9);dim(TG_G_pred9)

TG_G_pred10 <- predict(object = TGG_GC,newdata= TG_G_pred9,se.fit = T, type = "link",re.form = ~0)

TG_G_pred11<-data.frame(TG_G_pred9,fit.link=TG_G_pred10$fit,se.link=TG_G_pred10$se.fit)

TG_G_pred11$lci.link<-TG_G_pred11$fit.link-(1.96*TG_G_pred11$se.link)
TG_G_pred11$uci.link<-TG_G_pred11$fit.link+(1.96*TG_G_pred11$se.link)

TG_G_pred11$fit<-plogis(TG_G_pred11$fit.link)
TG_G_pred11$se<-plogis(TG_G_pred11$se.link)
TG_G_pred11$lci<-plogis(TG_G_pred11$lci.link)
TG_G_pred11$uci<-plogis(TG_G_pred11$uci.link)

head(TG_G_pred11);dim(TG_G_pred11)

#Day

summary(TGG_D)


TG_G_pred12 <- expand.grid(Day_Sampled = TG_Predictions_Day)
head(TG_G_pred12);dim(TG_G_pred12)

TG_G_pred13 <- predict(object = TGG_D,newdata= TG_G_pred12,se.fit = T, type = "link",re.form = ~0)

TG_G_pred14<-data.frame(TG_G_pred12,fit.link=TG_G_pred13$fit,se.link=TG_G_pred13$se.fit)

TG_G_pred14$lci.link<-TG_G_pred14$fit.link-(1.96*TG_G_pred14$se.link)
TG_G_pred14$uci.link<-TG_G_pred14$fit.link+(1.96*TG_G_pred14$se.link)

TG_G_pred14$fit<-plogis(TG_G_pred14$fit.link)
TG_G_pred14$se<-plogis(TG_G_pred14$se.link)
TG_G_pred14$lci<-plogis(TG_G_pred14$lci.link)
TG_G_pred14$uci<-plogis(TG_G_pred14$uci.link)

head(TG_G_pred14);dim(TG_G_pred14)

##Step 5: Visualisation----

#Crops
summary(TGG_Crops)
head(TG_G_pred2);dim(TG_G_pred2)

GG <- TG_G_pred2$X1km_Prop_Crops == TG_Predictions_Crops [10]
GGG <- TG_G_pred2$Day_Sampled == TG_Predictions_Day [10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_G_pred2$Day_Sampled[GG],rev(TG_G_pred2$Day_Sampled[GG])), y = c(TG_G_pred2$lci[GG],rev(TG_G_pred2$uci[GG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred2$Day_Sampled[GG],y = TG_G_pred2$fit[GG],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$X1km_Prop_Crops,y = FDModel$TG_G,xlab = "Crops within 1km",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_G_pred2$X1km_Prop_Crops[GGG],rev(TG_G_pred2$X1km_Prop_Crops[GGG])), y = c(TG_G_pred2$lci[GGG],rev(TG_G_pred2$uci[GGG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred2$X1km_Prop_Crops[GGG],y = TG_G_pred2$fit[GGG],lwd = 2,col = 'grey30',lty = 1)

#NDVI field
summary(TGG_NDVIfield)
head(TG_G_pred5);dim(TG_G_pred5)

GG.2 <- TG_G_pred5$NDVImean_Field == TG_Predictions_NDVIField[10]
GGG.2 <- TG_G_pred5$Day_Sampled == TG_Predictions_Day [10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_G_pred5$Day_Sampled[GG.2],rev(TG_G_pred5$Day_Sampled[GG.2])), y = c(TG_G_pred5$lci[GG.2],rev(TG_G_pred5$uci[GG.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred5$Day_Sampled[GG.2],y = TG_G_pred5$fit[GG.2],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$NDVImean_Field,y = FDModel$TG_G,xlab = "Field NDVI",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_G_pred5$NDVImean_Field[GGG.2],rev(TG_G_pred5$NDVImean_Field[GGG.2])), y = c(TG_G_pred5$lci[GGG.2],rev(TG_G_pred5$uci[GGG.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred5$NDVImean_Field[GGG.2],y = TG_G_pred5$fit[GGG.2],lwd = 2,col = 'grey30',lty = 1)

#Field Area
summary(TGG_FieldArea)
head(TG_G_pred8);dim(TG_G_pred8)

GG.3 <- TG_G_pred8$Field_Size_Scaled == TG_Predictions_FieldSize_Scaled [10]
GGG.3 <- TG_G_pred8$Day_Sampled == TG_Predictions_Day [10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_G_pred8$Day_Sampled[GG.3],rev(TG_G_pred8$Day_Sampled[GG.3])), y = c(TG_G_pred8$lci[GG.3],rev(TG_G_pred8$uci[GG.3])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred8$Day_Sampled[GG.3],y = TG_G_pred8$fit[GG.3],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$Field_Size_Scaled,y = FDModel$TG_G,xlab = "Field Sixe (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_G_pred8$Field_Size_Scaled),to=max(TG_G_pred8$Field_Size_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1),cex.axis=1)

polygon(x = c(TG_G_pred8$Field_Size_Scaled[GGG.3],rev(TG_G_pred8$Field_Size_Scaled[GGG.3])), y = c(TG_G_pred8$lci[GGG.3],rev(TG_G_pred8$uci[GGG.3])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred8$Field_Size_Scaled[GGG.3],y = TG_G_pred8$fit[GGG.3],lwd = 2,col = 'grey30',lty = 1)

#Ground Cover
summary(TGG_GC)
head(TG_G_pred11);dim(TG_G_pred11)

GG.4 <- TG_G_pred11$GC == TG_Predictions_GC [10]
GGG.4 <- TG_G_pred11$Day_Sampled == TG_Predictions_Day [10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_G_pred11$Day_Sampled[GG.4],rev(TG_G_pred11$Day_Sampled[GG.4])), y = c(TG_G_pred11$lci[GG.4],rev(TG_G_pred11$uci[GG.4])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred11$Day_Sampled[GG.4],y = TG_G_pred11$fit[GG.4],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$GC,y = FDModel$TG_G,xlab = "Ground Cover (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_G_pred11$GC[GGG.4],rev(TG_G_pred11$GC[GGG.4])), y = c(TG_G_pred11$lci[GGG.4],rev(TG_G_pred11$uci[GGG.4])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred11$GC[GGG.4],y = TG_G_pred11$fit[GGG.4],lwd = 2,col = 'grey30',lty = 1)

#Day

summary(TGG_D)
head(TG_G_pred14);dim(TG_G_pred14)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))

polygon(x = c(TG_G_pred14$Day_Sampled,rev(TG_G_pred14$Day_Sampled)), y = c(TG_G_pred14$lci,rev(TG_G_pred14$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred14$Day_Sampled,y = TG_G_pred14$fit,lwd = 2,col = 'grey30',lty = 1)

#Figures----

##Main document----

dev.new(height=15,width=20,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(3,4),mgp=c(1.8,0.5,0),xpd = T)

plot(x = 1:2,y = TG_C_pred2$fit [AA],xlab = " ",ylab = 'Probability of Occurrence', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,1),xaxt = "n",xlim = c(0,3))
axis(side=1,at=c(1,2.3),labels=c('Edge','Interior'),cex.axis = 0.9)
mtext(side=3,line=0,at = -0.3,'a)',cex=0.8)

arrows(x0=1:2, y0=TG_C_pred2$lci [AA],x1=1:2, y1=TG_C_pred2$uci[AA],angle=90,length=0.1, code=3, lwd=2,col = "black")
points(x = jitter(raw_x, factor = 1),y = FDModel$TG_C, pch = 16, cex = 0.4, col = "black")

plot(x = jitter(FDModel$X1km_Rip_Prop,factor = 1),y = FDModel$TG_C,xlab = "Proprtion Riparian (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_C_pred2$X1km_Rip_Prop),to=max(TG_C_pred2$X1km_Rip_Prop),length.out=4),labels=round(seq(from=min(FDModel$X1km_Rip_Prop),to=max(FDModel$X1km_Rip_Prop),length.out=4)),cex.axis=1)
mtext(side=3,line=0,at = -0.7,'b)',cex=0.8)

polygon(x = c(TG_C_pred2$X1km_Rip_Prop[A_A_],rev(TG_C_pred2$X1km_Rip_Prop[A_A_])), y = c(TG_C_pred2$lci[A_A_],rev(TG_C_pred2$uci[A_A_])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred2$X1km_Rip_Prop[A_A_],y = TG_C_pred2$fit[A_A_],lwd = 2,col = 'grey30',lty = 1)


plot(x = jitter(FDModel$Crop_Age_Days,factor = 1),y = FDModel$TG_C,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_C_pred2$Crop_Age_Days),to=max(TG_C_pred2$Crop_Age_Days),length.out=4),labels=round(seq(from=min(FDModel$Crop_Age_Days),to=max(FDModel$Crop_Age_Days),length.out=4)),cex.axis=1)
mtext(side=3,line=0,at = 45,'c)',cex=0.8)
mtext(side=3,line=1,at = 20,'----------------------------------Small fungivore beetles and non-web building spiders----------------------------------',cex=0.7, font = 2)


polygon(x = c(TG_C_pred2$Crop_Age_Days[AAA],rev(TG_C_pred2$Crop_Age_Days[AAA])), y = c(TG_C_pred2$lci[AAA],rev(TG_C_pred2$uci[AAA])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred2$Crop_Age_Days[AAA],y = TG_C_pred2$fit[AAA],lwd = 2,col = 'grey30',lty = 1)


plot(x = jitter(FDModel$Day_Sampled,factor = 1),y = FDModel$TG_C,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1), xaxt = 'n')
axis(side=1, at=seq(from=min(TG_C_pred2$Day_Sampled),to=max(TG_C_pred2$Day_Sampled),length.out=4),labels=round(seq(from=min(FDModel$Day_Sampled),to=max(FDModel$Day_Sampled),length.out=4)),cex.axis=1)
mtext(side=3,line=0,at = 10,'d)',cex=0.8)
mtext(side=1,line=2.5,at = 10,'Autumn/Winter',cex=0.6)
mtext(side=1,line=2.5,at = 157,'Spring',cex=0.6)
arrows(60,-0.48,135,-0.48, length =0.05)

polygon(x = c(TG_C_pred2$Day_Sampled[A_A],rev(TG_C_pred2$Day_Sampled[A_A])), y = c(TG_C_pred2$lci[A_A],rev(TG_C_pred2$uci[A_A])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred2$Day_Sampled[A_A],y = TG_C_pred2$fit[A_A],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$NDVIsum_1km,y = FDModel$TG_E,xlab = "NDVI within 1km",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
mtext(side=3,line=0,at = 360,'e)',cex=0.8)
axis(side=1, at=seq(from=min(TG_E_pred2$NDVIsum_1km),to=max(TG_E_pred2$NDVIsum_1km),length.out=3),labels=round(seq(from=min(FDModel$NDVIsum_1km),to=max(FDModel$NDVIsum_1km),length.out=3)),cex.axis=1)

polygon(x = c(TG_E_pred2$NDVIsum_1km[E_E],rev(TG_E_pred2$NDVIsum_1km[E_E])), y = c(TG_E_pred2$lci[E_E],rev(TG_E_pred2$uci[E_E])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_E_pred2$NDVIsum_1km[E_E],y = TG_E_pred2$fit[E_E],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_E,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
mtext(side=3,line=0,at = 45,'f)',cex=0.8)
mtext(side=3,line=1.2,at = 10,'------Predatory beetles, lacewings and dragonflies------',cex=0.7,font=2)
axis(side=1, at=seq(from=min(TG_E_pred2$Crop_Age_Days),to=max(TG_E_pred2$Crop_Age_Days),length.out=4),labels=round(seq(from=min(FDModel$Crop_Age_Days),to=max(FDModel$Crop_Age_Days),length.out=4)),cex.axis=1)

polygon(x = c(TG_E_pred2$Crop_Age_Days[EE],rev(TG_E_pred2$Crop_Age_Days[EE])), y = c(TG_E_pred2$lci[EE],rev(TG_E_pred2$uci[EE])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_E_pred2$Crop_Age_Days[EE],y = TG_E_pred2$fit[EE],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_E_pred2$Crop_Age_Days[EEE],rev(TG_E_pred2$Crop_Age_Days[EEE])), y = c(TG_E_pred2$lci[EEE],rev(TG_E_pred2$uci[EEE])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_E_pred2$Crop_Age_Days[EEE],y = TG_E_pred2$fit[EEE],lwd = 2,col = 'grey30',lty = 2)
legend('left',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)


plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_G,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_D_pred2$Crop_Age_Days),to=max(TG_D_pred2$Crop_Age_Days),length.out=4),labels=round(seq(from=min(FDModel$Crop_Age_Days),to=max(FDModel$Crop_Age_Days),length.out=4)),cex.axis=1)
mtext(side=3,line=0,at = 45,'g)',cex=0.8)

polygon(x = c(TG_D_pred2$Crop_Age_Days[DD],rev(TG_D_pred2$Crop_Age_Days[DD])), y = c(TG_D_pred2$lci[DD],rev(TG_D_pred2$uci[DD])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred2$Crop_Age_Days[DD],y = TG_D_pred2$fit[DD],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_D_pred2$Crop_Age_Days[DDD],rev(TG_D_pred2$Crop_Age_Days[DDD])), y = c(TG_D_pred2$lci[DDD],rev(TG_D_pred2$uci[DDD])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred2$Crop_Age_Days[DDD],y = TG_D_pred2$fit[DDD],lwd = 2,col = 'grey30',lty = 2)
legend('bottomright',legend = c('Edge', "Interior"), lty = c(1,2), col = 'grey30',pt.cex = 1)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
mtext(side=3,line=0,at = 10,'h)',cex=0.8)
axis(side=1, at=seq(from=min(TG_D_pred2$Day_Sampled),to=max(TG_D_pred2$Day_Sampled),length.out=4),labels=round(seq(from=min(FDModel$Day_Sampled),to=max(FDModel$Day_Sampled),length.out=4)),cex.axis=1)
mtext(side=1,line=2.5,at = 10,'Autumn/Winter',cex=0.6)
mtext(side=1,line=2.5,at = 157,'Spring',cex=0.6)
arrows(60,-0.48,135,-0.48, length =0.05)
mtext(side=3,line=1.2,at = -40,'----------------Web building spiders----------------',cex=0.7, font = 2)


polygon(x = c(TG_D_pred2$Day_Sampled[D_D],rev(TG_D_pred2$Day_Sampled[D_D])), y = c(TG_D_pred2$lci[D_D],rev(TG_D_pred2$uci[D_D])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred2$Day_Sampled[D_D],y = TG_D_pred2$fit[D_D],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$GC,y = FDModel$TG_F,xlab = "Ground Cover (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
mtext(side=3,line=0,at = 2,'i)',cex=0.8)
axis(side=1, at=seq(from=min(TG_F_pred2$GC),to=max(TG_F_pred2$GC),length.out=4),labels=round(seq(from=min(FDModel$GC),to=max(FDModel$GC),length.out=4)),cex.axis=1)
mtext(side=3,line=0.9,at = 50,'Hoverflies and omnivorous \nflies, ants and crickets',cex=0.7,font=2)


polygon(x = c(TG_F_pred2$GC,rev(TG_F_pred2$GC)), y = c(TG_F_pred2$lci,rev(TG_F_pred2$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_F_pred2$GC,y = TG_F_pred2$fit,lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$Field_Size_Scaled,y = FDModel$TG_B,xlab = "Field Size (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,xaxt = 'n',ylim = c(0,1))
axis(side=1, at=seq(from=min(TG_B_pred2$Field_Size_Scaled),to=max(TG_B_pred2$Field_Size_Scaled),length.out=4),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=4)/10000,1),cex.axis=1)
mtext(side=3,line=0,at = -2.3,'j)',cex=0.8)
mtext(side=3,line=0.9,at = 0,'Herbivorous grasshoppers, \nHemipterans and beetles',cex=0.7,font=2)


polygon(x = c(TG_B_pred2$Field_Size_Scaled,rev(TG_B_pred2$Field_Size_Scaled)), y = c(TG_B_pred2$lci,rev(TG_B_pred2$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_B_pred2$Field_Size_Scaled,y = TG_B_pred2$fit,lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$X1km_Prop_Crops,y = FDModel$TG_G,xlab = "Crops within 1km",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1), xaxt = 'n')
mtext(side=3,line=0,at = 79,'k)',cex=0.8)
axis(side=1, at=seq(from=min(TG_G_pred2$X1km_Prop_Crops),to=max(TG_G_pred2$X1km_Prop_Crops),length.out=4),labels=round(seq(from=min(FDModel$X1km_Prop_Crops),to=max(FDModel$X1km_Prop_Crops),length.out=4)),cex.axis=1)

polygon(x = c(TG_G_pred2$X1km_Prop_Crops[GGG],rev(TG_G_pred2$X1km_Prop_Crops[GGG])), y = c(TG_G_pred2$lci[GGG],rev(TG_G_pred2$uci[GGG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred2$X1km_Prop_Crops[GGG],y = TG_G_pred2$fit[GGG],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
mtext(side=3,line=0,at = 9,'l)',cex=0.8)
mtext(side=3,line=1.2,at = -40,'--------Small and large herbivorous Hemiptera--------',cex=0.7,font=2)
axis(side=1, at=seq(from=min(TG_G_pred2$Day_Sampled),to=max(TG_G_pred2$Day_Sampled),length.out=4),labels=round(seq(from=min(FDModel$Day_Sampled),to=max(FDModel$Day_Sampled),length.out=4)),cex.axis=1)
mtext(side=1,line=2.5,at = 10,'Autumn/Winter',cex=0.6)
mtext(side=1,line=2.5,at = 157,'Spring',cex=0.6)
arrows(60,-0.48,135,-0.48, length =0.05)


polygon(x = c(TG_G_pred2$Day_Sampled[GG],rev(TG_G_pred2$Day_Sampled[GG])), y = c(TG_G_pred2$lci[GG],rev(TG_G_pred2$uci[GG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred2$Day_Sampled[GG],y = TG_G_pred2$fit[GG],lwd = 2,col = 'grey30',lty = 1)



##Supporting----

###C Group----

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = 1:2,y = TG_C_pred5$fit [CC],xlab = " ",ylab = 'Probability of Occurrence', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,1),xaxt = "n",xlim = c(0,3))
axis(side=1,at=1:2,labels=c('Edge','Interior'))

arrows(x0=1:2, y0=TG_C_pred5$lci [CC],x1=1:2, y1=TG_C_pred5$uci[CC],angle=90,length=0.2, code=3, lwd=2,col = "black")

points(x = jitter(raw_x, factor = 1),y = FDModel$TG_C, pch = 16, cex = 0.4, col = "black")
mtext(side=3,line=0,at = -0.3,'a)',cex=1)

plot(x = jitter(FDModel$Crop_Age_Days,factor = 1),y = FDModel$TG_C,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_C_pred5$Crop_Age_Days),to=max(TG_C_pred5$Crop_Age_Days),length.out=6),labels=round(seq(from=min(FDModel$Crop_Age_Days),to=max(FDModel$Crop_Age_Days),length.out=6)),cex.axis=0.95)
mtext(side=3,line=0,at = 46,'b)',cex=1)

polygon(x = c(TG_C_pred5$Crop_Age_Days[CCC],rev(TG_C_pred5$Crop_Age_Days[CCC])), y = c(TG_C_pred5$lci[CCC],rev(TG_C_pred5$uci[CCC])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred5$Crop_Age_Days[CCC],y = TG_C_pred5$fit[CCC],lwd = 2,col = 'grey30',lty = 1)


plot(x = jitter(FDModel$Day_Sampled,factor = 1),y = FDModel$TG_C,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_C_pred5$Day_Sampled),to=max(TG_C_pred5$Day_Sampled),length.out=6),labels=round(seq(from=min(FDModel$Day_Sampled),to=max(FDModel$Day_Sampled),length.out=6)),cex.axis=0.9)
mtext(side=3,line=0,at = 12,'c)',cex=1)
mtext(side=1,line=3,at = 10,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 157,'Spring',cex=0.8)
arrows(50,-0.42,140,-0.42, length =0.1)

polygon(x = c(TG_C_pred5$Day_Sampled[C_C],rev(TG_C_pred5$Day_Sampled[C_C])), y = c(TG_C_pred5$lci[C_C],rev(TG_C_pred5$uci[C_C])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred5$Day_Sampled[C_C],y = TG_C_pred5$fit[C_C],lwd = 2,col = 'grey30',lty = 1)


plot(x = jitter(FDModel$Field_Size_Scaled,factor = 1),y = FDModel$TG_C,xlab = "Field Size (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_C_pred5$Field_Size_Scaled),to=max(TG_C_pred5$Field_Size_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1),cex.axis=0.95)
mtext(side=3,line=0,at = -2.3,'d)',cex=1)

polygon(x = c(TG_C_pred5$Field_Size_Scaled[C_C_],rev(TG_C_pred5$Field_Size_Scaled[C_C_])), y = c(TG_C_pred5$lci[C_C_],rev(TG_C_pred5$uci[C_C_])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_C_pred5$Field_Size_Scaled[C_C_],y = TG_C_pred5$fit[C_C_],lwd = 2,col = 'grey30',lty = 1)

###D Group----

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(3,3),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$TG_G,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_D_pred5$Crop_Age_Days),to=max(TG_D_pred5$Crop_Age_Days),length.out=4),labels=round(seq(from=min(FDModel$Crop_Age_Days),to=max(FDModel$Crop_Age_Days),length.out=4)),cex.axis=0.95)
mtext(side=3,line=0,at = 45,'a)',cex=0.95)

polygon(x = c(TG_D_pred5$Crop_Age_Days[DD.2],rev(TG_D_pred5$Crop_Age_Days[DD.2])), y = c(TG_D_pred5$lci[DD.2],rev(TG_D_pred5$uci[DD.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred5$Crop_Age_Days[DD.2],y = TG_D_pred5$fit[DD.2],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(TG_D_pred5$Crop_Age_Days[DDD.2],rev(TG_D_pred5$Crop_Age_Days[DDD.2])), y = c(TG_D_pred5$lci[DDD.2],rev(TG_D_pred5$uci[DDD.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred5$Crop_Age_Days[DDD.2],y = TG_D_pred5$fit[DDD.2],lwd = 2,col = 'grey30',lty = 2)
legend('bottomright',legend = c('Edge', "Interior"), lty = c(1,2), col = 'grey30',pt.cex = 1)


plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_D_pred5$Day_Sampled),to=max(TG_D_pred5$Day_Sampled),length.out=4),labels=round(seq(from=min(FDModel$Day_Sampled),to=max(FDModel$Day_Sampled),length.out=4)),cex.axis=0.95)
mtext(side=3,line=0,at = 12,'b)',cex=0.95)
mtext(side=1,line=3,at = 10,'Autumn/Winter',cex=0.6)
mtext(side=1,line=3,at = 157,'Spring',cex=0.6)
arrows(60,-0.55,135,-0.55, length =0.05)

polygon(x = c(TG_D_pred5$Day_Sampled[D_D.2],rev(TG_D_pred5$Day_Sampled[D_D.2])), y = c(TG_D_pred5$lci[D_D.2],rev(TG_D_pred5$uci[D_D.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred5$Day_Sampled[D_D.2],y = TG_D_pred5$fit[D_D.2],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$NDVImean_Field,y = FDModel$TG_G,xlab = "Field NDVI",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_D_pred5$NDVImean_Field),to=max(TG_D_pred5$NDVImean_Field),length.out=4),labels=round(seq(from=min(FDModel$NDVImean_Field),to=max(FDModel$NDVImean_Field),length.out=4),2),cex.axis=0.9)
mtext(side=3,line=0,at = 0.12,'c)',cex=0.95)

polygon(x = c(TG_D_pred5$NDVImean_Field[D_DD.2],rev(TG_D_pred5$NDVImean_Field[D_DD.2])), y = c(TG_D_pred5$lci[D_DD.2],rev(TG_D_pred5$uci[D_DD.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred5$NDVImean_Field[D_DD.2],y = TG_D_pred5$fit[D_DD.2],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$GC,y = FDModel$TG_G,xlab = "Ground Cover (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.9,ylim = c(0,1))
mtext(side=3,line=0,at = 3,'d)',cex=0.95)

polygon(x = c(TG_D_pred8$GC[D_DD.3],rev(TG_D_pred8$GC[D_DD.3])), y = c(TG_D_pred8$lci[D_DD.3],rev(TG_D_pred8$uci[D_DD.3])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred8$GC[D_DD.3],y = TG_D_pred8$fit[D_DD.3],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$Field_Size_Scaled,y = FDModel$TG_G,xlab = "Field Size (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_D_pred11$Field_Size_Scaled),to=max(TG_D_pred11$Field_Size_Scaled),length.out=4),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=4)/10000,1),cex.axis=1)
mtext(side=3,line=0,at = -2.3,'e)',cex=0.95)

polygon(x = c(TG_D_pred11$Field_Size_Scaled[D_DD.4],rev(TG_D_pred11$Field_Size_Scaled[D_DD.4])), y = c(TG_D_pred11$lci[D_DD.4],rev(TG_D_pred11$uci[D_DD.4])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred11$Field_Size_Scaled[D_DD.4],y = TG_D_pred11$fit[D_DD.4],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$Height,y = FDModel$TG_G,xlab = "Crop Height (cm)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_D_pred14$Height),to=max(TG_D_pred14$Height),length.out=5),labels=round(seq(from=min(FDModel$Height),to=max(FDModel$Height),length.out=5)),cex.axis=0.95)
mtext(side=3,line=0,at = 15,'f)',cex=0.95)

polygon(x = c(TG_D_pred14$Height[D_DD.5],rev(TG_D_pred14$Height[D_DD.5])), y = c(TG_D_pred14$lci[D_DD.5],rev(TG_D_pred14$uci[D_DD.5])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred14$Height[D_DD.5],y = TG_D_pred14$fit[D_DD.5],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$X1km_Rip_Prop,y = FDModel$TG_G,xlab = "Riparian within 1km (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_D_pred17$X1km_Rip_Prop),to=max(TG_D_pred17$X1km_Rip_Prop),length.out=5),labels=round(seq(from=min(FDModel$X1km_Rip_Prop),to=max(FDModel$X1km_Rip_Prop),length.out=5)),cex.axis=0.95)
mtext(side=3,line=0,at = -0.7,'g)',cex=0.95)

polygon(x = c(TG_D_pred17$X1km_Rip_Prop[D_DD.6],rev(TG_D_pred17$X1km_Rip_Prop[D_DD.6])), y = c(TG_D_pred17$lci[D_DD.6],rev(TG_D_pred17$uci[D_DD.6])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_D_pred17$X1km_Rip_Prop[D_DD.6],y = TG_D_pred17$fit[D_DD.6],lwd = 2,col = 'grey30',lty = 1)


###G Group----

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Day_Sampled,y = FDModel$TG_G,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_G_pred2$Day_Sampled),to=max(TG_G_pred2$Day_Sampled),length.out=6),labels=round(seq(from=min(FDModel$Day_Sampled),to=max(FDModel$Day_Sampled),length.out=6)),cex.axis=0.9)
mtext(side=3,line=0,at = 12,'a)',cex=1)
mtext(side=1,line=3,at = 10,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 157,'Spring',cex=0.8)
arrows(50,-0.42,140,-0.42, length =0.1)

polygon(x = c(TG_G_pred2$Day_Sampled[GG],rev(TG_G_pred2$Day_Sampled[GG])), y = c(TG_G_pred2$lci[GG],rev(TG_G_pred2$uci[GG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred2$Day_Sampled[GG],y = TG_G_pred2$fit[GG],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$NDVImean_Field,y = FDModel$TG_G,xlab = "Field NDVI",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
mtext(side=3,line=0,at = 0.12,'b)', cex=1)
axis(side=1, at=seq(from=min(TG_G_pred5$NDVImean_Field),to=max(TG_G_pred5$NDVImean_Field),length.out=5),labels=round(seq(from=min(FDModel$NDVImean_Field),to=max(FDModel$NDVImean_Field),length.out=5),2),cex.axis=0.95)

polygon(x = c(TG_G_pred5$NDVImean_Field[GGG.2],rev(TG_G_pred5$NDVImean_Field[GGG.2])), y = c(TG_G_pred5$lci[GGG.2],rev(TG_G_pred5$uci[GGG.2])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred5$NDVImean_Field[GGG.2],y = TG_G_pred5$fit[GGG.2],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$Field_Size_Scaled,y = FDModel$TG_G,xlab = "Field Size (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1),xaxt = 'n')
axis(side=1, at=seq(from=min(TG_G_pred8$Field_Size_Scaled),to=max(TG_G_pred8$Field_Size_Scaled),length.out=5),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=5)/10000,1),cex.axis=1)
mtext(side=3,line=0,at = -2.3,'c)', cex=1)

polygon(x = c(TG_G_pred8$Field_Size_Scaled[GGG.3],rev(TG_G_pred8$Field_Size_Scaled[GGG.3])), y = c(TG_G_pred8$lci[GGG.3],rev(TG_G_pred8$uci[GGG.3])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred8$Field_Size_Scaled[GGG.3],y = TG_G_pred8$fit[GGG.3],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$GC,y = FDModel$TG_G,xlab = "Ground Cover (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95,ylim = c(0,1))
mtext(side=3,line=0,at = 3,'d)', cex=1)

polygon(x = c(TG_G_pred11$GC[GGG.4],rev(TG_G_pred11$GC[GGG.4])), y = c(TG_G_pred11$lci[GGG.4],rev(TG_G_pred11$uci[GGG.4])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=TG_G_pred11$GC[GGG.4],y = TG_G_pred11$fit[GGG.4],lwd = 2,col = 'grey30',lty = 1)


#END----

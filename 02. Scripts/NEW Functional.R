#From Tutorial: Methods for assessing functional responses to environmental gradients -- Kleyer et al. -- June 22, 2009

#Libraries----

library("ade4")
library("dplyr")
library("tidyverse")
library("vegan")
library('fpc')
library("ggplot2")
library("tidyr")


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
rownames(l1_short) <- c("H", "GC", "In   ", "Out","Day", "A      ", "FA", "W", "NF", "NL")
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

dev.new(height=10,width=15,dpi=80,pointsize=14,noRStudioGD = T)
ggplot(heatmap_data, aes(x = Group, y = Trait_Value, fill = Proportion)) +
  geom_tile() +
  facet_wrap(~Trait, scales = "free_y") +
  scale_fill_gradient(low = "white", high = "darkblue") +
  theme_minimal() +
  labs(title = "Trait composition by group")


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

FDModel$Fun_Rich[is.na(FDModel$Fun_Rich)] <- 0.00001

#Functional Evenness 

FUNevenness <- aggregate(trait_group ~ Site, data = invert_trait_group,                            FUN = function(x) {
                           H <- diversity(table(x), index = "shannon")
                           S <- length(unique(x))
                           if(S <= 1) return(NA)
                           H / log(S)  # Pielou's evenness
                         })
FDModel <- merge(FDModel, FUNevenness, by = "Site", all.x = T)
head(FDModel);dim(FDModel)
colnames(FDModel)[4] <- "Fun_Even"
head(FDModel);dim(FDModel)


FDModel <- merge(FDModel,variables,by = "Site")
head(FDModel);dim(FDModel)

FDModel <- merge(FDModel,coords,by = "Site")
head(FDModel);dim(FDModel)




#Functional Richness----
##Step 1: Design Variables
##Step 2: Environmental Variables
##Step 3: Check Spatial Autocorrelation
##Step 4: Predictions
##Step 5: Visalisation


#Functional Diversity----
##Step 1: Design Variables
##Step 2: Environmental Variables
##Step 3: Check Spatial Autocorrelation
##Step 4: Predictions
##Step 5: Visalisation


#Functional Evenness----
##Step 1: Design Variables
##Step 2: Environmental Variables
##Step 3: Check Spatial Autocorrelation
##Step 4: Predictions
##Step 5: Visalisation


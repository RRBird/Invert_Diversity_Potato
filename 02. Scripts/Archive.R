
#Checking coefficents of models-----


#top model is DxA with crops, field NDVI, NDVI 1km and GC witgin 2 AICc within 2 AICc
#all ranked within 2 improved log liklihood
#now going to look at effect size - see if it overlaps 0 (don't want to overlap 0 - means not a strong influence)

#models I'm investigating
ModList <- list(Div_DxA,Div_Crop, Div_FieldNDVI,Div_NDVI1km, Div_GC)
Coefs_list <- list()
m <- 1

for (M in ModList) {
  Coefs <- data.frame(Estimate = 
                        c(summary(M)$coefficients[,1]),
                      SE = 
                        c(summary(M)$coefficients[,2]),
                      Term = rownames(summary(
                        M)$coefficients))
  
  Coefs$lci <- Coefs$Estimate - (Coefs$SE * 1.96)
  Coefs$uci <- Coefs$Estimate + (Coefs$SE * 1.96)
  rownames(Coefs) <- Coefs$Term
  
  Coefs_list[[m]] <- Coefs
  
  m <- m+1
}

Coefs_list

term_replacements <- c(
  "Age_Scaled"    = "Age", "Day_Scaled"    = "Day",
  "Day_Scaled:Age_Scaled" = "Day:Age", "NDVIsum_1km" = "NDVI 1km",
  "NDVImean_Field"  = "NDVI Field", 
  "X1km_Prop_Crops"   = "Crops 1km")


Coefs_list <- lapply(Coefs_list, function(df) {
  matched <- term_replacements[df$Term]
  df$Term <- ifelse(is.na(matched), df$Term, matched)
  df
})

Coefs_list <- lapply(Coefs_list, function(df) {
  rownames(df) <- df$Term
  df
})

#check coefficents try to get the distance between the two rows better (par for each plot??)
#Age crosses 0 but is in interaction so all good - all Enviro models have the extra parameter crossing 0 so just DxA to continue



dev.new(height=10,width=15,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(5,5,2,3),mfrow=c(2,3),mgp=c(2.5,1,0),xpd = T,oma =c(0,0,1,0))

plot(Coefs_list[[1]]$Estimate, rev(1:nrow(Coefs_list[[1]])),xlim = c(min(Coefs_list[[1]]$lci),max(Coefs_list[[1]]$uci)),las =1, cex = 1.8, ylab = "", xlab = expression(bold("Effect Size")),pch = 20, yaxt = "n", col = "black")
axis(side = 2, at = rev(1:nrow(Coefs_list[[1]])),labels=rownames(Coefs_list[[1]]),las =1)
arrows(Coefs_list[[1]]$uci,rev(1:nrow(Coefs_list[[1]])), Coefs_list[[1]]$lci,rev(1:nrow(Coefs_list[[1]])), lwd =0.8,code = 0)
arrows(0,0.9,0,4.1,code = 0, lwd = 0.8)
mtext("a)", line=0.2, at = -0.7,cex = 0.9)

plot(Coefs_list[[2]]$Estimate, rev(1:nrow(Coefs_list[[2]])),xlim = c(min(Coefs_list[[2]]$lci),max(Coefs_list[[2]]$uci)),las =1, cex = 1.8, ylab = "", xlab = expression(bold("Effect Size")),pch = 20, yaxt = "n", col = "black")
axis(side = 2, at = rev(1:nrow(Coefs_list[[2]])),labels=rownames(Coefs_list[[2]]),las =1)
arrows(Coefs_list[[2]]$uci,rev(1:nrow(Coefs_list[[2]])), Coefs_list[[2]]$lci,rev(1:nrow(Coefs_list[[2]])), lwd =0.8,code = 0)
arrows(0,0.8,0,5.15,code = 0, lwd = 0.8)
mtext("b)", line=0.2, at = -0.9,cex = 0.9)

plot(Coefs_list[[3]]$Estimate, rev(1:nrow(Coefs_list[[3]])),xlim = c(min(Coefs_list[[3]]$lci),max(Coefs_list[[3]]$uci)),las =1, cex = 1.8, ylab = "", xlab = expression(bold("Effect Size")),pch = 20, yaxt = "n", col = "black")
axis(side = 2, at = rev(1:nrow(Coefs_list[[3]])),labels=rownames(Coefs_list[[3]]),las =1)
arrows(Coefs_list[[3]]$uci,rev(1:nrow(Coefs_list[[3]])), Coefs_list[[3]]$lci,rev(1:nrow(Coefs_list[[3]])), lwd =0.8,code = 0)
arrows(0,0.8,0,5.15,code = 0, lwd = 0.8)
mtext("c)", line=0.2, at = -3,cex = 0.9)

plot(Coefs_list[[4]]$Estimate, rev(1:nrow(Coefs_list[[4]])),xlim = c(min(Coefs_list[[4]]$lci),max(Coefs_list[[4]]$uci)),las =1, cex = 1.8, ylab = "", xlab = expression(bold("Effect Size")),pch = 20, yaxt = "n", col = "black")
axis(side = 2, at = rev(1:nrow(Coefs_list[[4]])),labels=rownames(Coefs_list[[4]]),las =1)
arrows(Coefs_list[[4]]$uci,rev(1:nrow(Coefs_list[[4]])), Coefs_list[[4]]$lci,rev(1:nrow(Coefs_list[[4]])), lwd =0.8,code = 0)
arrows(0,0.8,0,5.15,code = 0, lwd = 0.8)
mtext("d)", line=0.2, at = -0.7,cex = 0.9)

plot(Coefs_list[[5]]$Estimate, rev(1:nrow(Coefs_list[[5]])),xlim = c(min(Coefs_list[[5]]$lci),max(Coefs_list[[5]]$uci)),las =1, cex = 1.8, ylab = "", xlab = expression(bold("Effect Size")),pch = 20, yaxt = "n", col = "black")
axis(side = 2, at = rev(1:nrow(Coefs_list[[5]])),labels=rownames(Coefs_list[[5]]),las =1)
arrows(Coefs_list[[5]]$uci,rev(1:nrow(Coefs_list[[5]])), Coefs_list[[5]]$lci,rev(1:nrow(Coefs_list[[5]])), lwd =0.8,code = 0)
arrows(0,0.8,0,5.15,code = 0, lwd = 0.8)
mtext("e)", line=0.2, at = -0.7,cex = 0.9)





#OLD FUNCTIONAL (trait group richness and diversity)----
#extracting trait group

Trait_Group <- data.frame(Morphospecies = names(spe.group2), trait_group = spe.group2)
head(Trait_Group);dim(Trait_Group)


invert_trait_group <- merge(Trait_Group,invert_filtered, by = "Morphospecies")
head(invert_trait_group);dim(invert_trait_group)
colnames(invert_trait_group)[colnames(invert_trait_group) == "ID"] <- "Site"

#Creating Modelling data and calculating functional Diversity 

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

plot(x = FDModel$Crop_Age_Days,y = FDModel$Fun_Div,xlab = "Crop Age (Days)",ylab = 'Trait Group Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.axis = 0.95)
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


dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Field_Area_Scaled,y = FDModel$Fun_Div,xlab = "Field Size (ha)",ylab = 'Trait Group Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(FUNdivpred3$Field_Area_Scaled),to=max(FUNdivpred3$Field_Area_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1),cex.axis=1)

polygon(x = c(FUNdivpred3$Field_Area_Scaled,rev(FUNdivpred3$Field_Area_Scaled)), y = c(FUNdivpred3$lci,rev(FUNdivpred3$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=FUNdivpred3$Field_Area_Scaled,y = FUNdivpred3$fit,lwd = 2,col = 'grey30')




#filter for orders included in functional analysis----

invert_filtered <- invert[invert$Order %in% c("Araneae", "Coleoptera", "Hemiptera", "Diptera"), ]

head(invert_filtered);dim(invert_filtered)

#log liklihood not improved by any below top model, so GC remain the top model



#FUN RICH STEPS----
#in case I need them two option exclude models with spatial auto correlation or do something to account for it so these will just stay here for now
##Step 4: Predictions----

summary(FUNRich_FieldArea)

FUNPredictions_Age <- seq(min(FDModel$Crop_Age_Days),max(FDModel$Crop_Age_Days),length.out=20)
FUNPredictions_FieldArea <- seq(min(FDModel$Field_Area_Scaled),max(FDModel$Field_Area_Scaled),length.out=20)


FUNrichpred <- expand.grid(Crop_Age_Days = FUNPredictions_Age, Field_Area_Scaled = FUNPredictions_FieldArea)
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

#Water
summary(FUNRich_Water)

FUNPredictions_Water <- seq(min(FDModel$X1km_Prop_Water),max(FDModel$X1km_Prop_Water),length.out=20)

FUNrichpred3 <- expand.grid(Crop_Age_Days = FUNPredictions_Age, X1km_Prop_Water = FUNPredictions_Water)
head(FUNrichpred);dim(FUNrichpred)

FUNrichpred4 <- predict(object = FUNRich_Water,newdata= FUNrichpred3,se.fit = T, type = "link",re.form = NA)

FUNrichpred5<-data.frame(FUNrichpred3,fit.link=FUNrichpred4$fit,se.link=FUNrichpred4$se.fit)

FUNrichpred5$lci.link<-FUNrichpred5$fit.link-(1.96*FUNrichpred5$se.link)
FUNrichpred5$uci.link<-FUNrichpred5$fit.link+(1.96*FUNrichpred5$se.link)

FUNrichpred5$fit<-exp(FUNrichpred5$fit.link)
FUNrichpred5$se<-exp(FUNrichpred5$se.link)
FUNrichpred5$lci<-exp(FUNrichpred5$lci.link)
FUNrichpred5$uci<-exp(FUNrichpred5$uci.link)

head(FUNrichpred5);dim(FUNrichpred5)


##Step 5: Visualisation----

#Field area
summary(FUNRich_FieldArea)
head(FUNrichpred2);dim(FUNrichpred2)


FF <- FUNrichpred2$Field_Area_Scaled == FUNPredictions_FieldArea[10]
F_F <- FUNrichpred2$Crop_Age_Days == FUNPredictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$Fun_Rich,xlab = "Crop Age (Days)",ylab = 'Trait Group Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 45,'a)',cex=1.1)

polygon(x = c(FUNrichpred2$Crop_Age_Days[FF],rev(FUNrichpred2$Crop_Age_Days[FF])), y = c(FUNrichpred2$lci[FF],rev(FUNrichpred2$uci[FF])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=FUNrichpred2$Crop_Age_Days[FF],y = FUNrichpred2$fit[FF],lwd = 2,col = 'grey30',lty = 1)

plot(x = FDModel$Field_Area_Scaled,y = FDModel$Fun_Rich,xlab = "Field Area (ha)",ylab = 'Trait Group Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(FUNrichpred2$Field_Area_Scaled),to=max(FUNrichpred2$Field_Area_Scaled),length.out=6),labels=round(seq(from=min(FDModel$Field_Area_m2),to=max(FDModel$Field_Area_m2),length.out=6)/10000,1),cex.axis=1)
mtext(side=3,line=0,at = -2.4,'b)',cex=1.1)

polygon(x = c(FUNrichpred2$Field_Area_Scaled[F_F],rev(FUNrichpred2$Field_Area_Scaled[F_F])), y = c(FUNrichpred2$lci[F_F],rev(FUNrichpred2$uci[F_F])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=FUNrichpred2$Field_Area_Scaled[F_F],y = FUNrichpred2$fit[F_F],lwd = 2,col = 'grey30')

#Water

summary(FUNRich_Water)
head(FUNrichpred5);dim(FUNrichpred5)


WW <- FUNrichpred5$X1km_Prop_Water == FUNPredictions_Water[10]
W_W <- FUNrichpred5$Crop_Age_Days == FUNPredictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = FDModel$Crop_Age_Days,y = FDModel$Fun_Rich,xlab = "Crop Age (Days)",ylab = 'Trait Group Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 45,'a)',cex=1.1)

polygon(x = c(FUNrichpred5$Crop_Age_Days[WW],rev(FUNrichpred5$Crop_Age_Days[WW])), y = c(FUNrichpred5$lci[WW],rev(FUNrichpred5$uci[WW])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=FUNrichpred5$Crop_Age_Days[WW],y = FUNrichpred5$fit[WW],lwd = 2,col = 'grey30',lty = 1)


plot(x = FDModel$X1km_Prop_Water,y = FDModel$Fun_Rich,xlab = "Proportion of Water in 1km",ylab = 'Trait Group Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 1,'b)',cex=1.1)

polygon(x = c(FUNrichpred5$X1km_Prop_Water[W_W],rev(FUNrichpred5$X1km_Prop_Water[W_W])), y = c(FUNrichpred5$lci[W_W],rev(FUNrichpred5$uci[W_W])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=FUNrichpred5$X1km_Prop_Water[W_W],y = FUNrichpred5$fit[W_W],lwd = 2,col = 'grey30',lty = 1)



#with updated data the null did converge


#Develops wings null didn't converge properly so investigating that

Richlist1$Develops_Wings

#use optimizer to get a convergence for AIC table
null_develop <- glmmTMB(Develops_Wings ~ 1 + (1 | Field), family = nbinom2, data = ModelRich2,control = glmmTMBControl(optimizer = optim))
summary(null_develop)


#Old fig

##Step 5 - Richness Model Visalisation----

summary(Rich_Water)
head(richpred2);dim(richpred2)

Predictions_Day[4] #winter
Predictions_Day[15] #Spring

RR <- richpred2$Day_Scaled == Predictions_Day[4] & richpred2$X1km_Prop_Water == Predictions_Water[10]
R_R <- richpred2$Day_Scaled == Predictions_Day[15] & richpred2$X1km_Prop_Water == Predictions_Water[10]
RRR <- richpred2$Day_Scaled == Predictions_Day[10] & richpred2$Age_Scaled == Predictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = TaxModel$Age_Scaled,y = TaxModel$Species_Rich,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(richpred2$Age_Scaled),to=max(richpred2$Age_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Crop_Age_Days),to=max(TaxModel$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.3,'a)',cex=1.1)

polygon(x = c(richpred2$Age_Scaled[RR],rev(richpred2$Age_Scaled[RR])), y = c(richpred2$lci[RR],rev(richpred2$uci[RR])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred2$Age_Scaled[RR],y = richpred2$fit[RR],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(richpred2$Age_Scaled[R_R],rev(richpred2$Age_Scaled[R_R])), y = c(richpred2$lci[R_R],rev(richpred2$uci[R_R])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred2$Age_Scaled[R_R],y = richpred2$fit[R_R],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)


plot(x = TaxModel$X1km_Prop_Water,y = TaxModel$Species_Rich,xlab = expression("Proportion Water within 1km"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 0.9,'b)',cex=1.1)

polygon(x = c(richpred2$X1km_Prop_Water[RRR],rev(richpred2$X1km_Prop_Water[RRR])), y = c(richpred2$lci[RRR],rev(richpred2$uci[RRR])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=richpred2$X1km_Prop_Water[RRR],y = richpred2$fit[RRR],lwd = 2,col = 'grey30')



plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Predator,xlab = expression("Day Sampled"),ylab = 'Predator Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.1,'a)',cex=1)
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-0.7,-3.7,1.15,-3.7, length =0.1)

axis(side=1, at=seq(from=min(richpred_herb.1$Day_Scaled),to=max(richpred_herb.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1))

polygon(x = c(richpred_preds.1$Day_Scaled[ADD],rev(richpred_preds.1$Day_Scaled[ADD])), y = c(richpred_preds.1$lci[ADD],rev(richpred_preds.1$uci[ADD])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_preds.1$Day_Scaled[ADD],y = richpred_preds.1$fit[ADD],lwd = 2,col = 'grey30')


Predictions_Age[5]
Predictions_Age[15]
seq(min(ModelRich2$Crop_Age_Days),max(ModelRich2$Crop_Age_Days),length.out=20)[5]
seq(min(ModelRich2$Crop_Age_Days),max(ModelRich2$Crop_Age_Days),length.out=20)[15]


#on line 2 I think
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
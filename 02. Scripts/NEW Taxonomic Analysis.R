
#Author: Rhiannon Bird
#Written under version R 4.5.1

#This script contains the analysis and modelling of taxpnomic richness and diversity


#Libraries----
library("AICcmodavg")
library('lme4')
library("glmmTMB")
library("dplyr")
library("DHARMa")
library("arm")
library("openxlsx")

#SPECIES RICHNESSS----


TaxModel$Day_Scaled <- scale(TaxModel$Day_Sampled)
TaxModel$Age_Scaled <- scale(TaxModel$Crop_Age_Days)


##Step 1 - Design Variables----

#position, age and day - all the various combinations of these 


head(TaxModel);dim(TaxModel)
str(TaxModel)

Rich_null <- glmmTMB(Species_Rich ~ 1 + (1 | Field), family = nbinom2, data = TaxModel)

Rich_P <- glmmTMB(Species_Rich ~ Position + (1 | Field), family = nbinom2, data = TaxModel)
Rich_A <- glmmTMB(Species_Rich ~ Crop_Age_Days + (1 | Field), family = nbinom2, data = TaxModel)
Rich_D <- glmmTMB(Species_Rich ~ Day_Scaled + (1 | Field), family = nbinom2, data = TaxModel)

Rich_PA <- glmmTMB(Species_Rich ~ Position + Crop_Age_Days + (1 | Field), family = nbinom2, data = TaxModel)
Rich_PD <- glmmTMB(Species_Rich ~ Position + Day_Sampled + (1 | Field), family = nbinom2, data = TaxModel)
Rich_DA <- glmmTMB(Species_Rich ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = nbinom2, data = TaxModel)

Rich_PxA <- glmmTMB(Species_Rich ~ Position * Age_Scaled + (1 | Field), family = nbinom2, data = TaxModel) 
Rich_PxD <- glmmTMB(Species_Rich ~ Position * Day_Sampled + (1 | Field), family = nbinom2, data = TaxModel) 
Rich_DxA <- glmmTMB(Species_Rich ~ Day_Sampled * Crop_Age_Days + (1 | Field), family = nbinom2, data = TaxModel) #model convergence problem

Rich_PAD <- glmmTMB(Species_Rich ~ Position + Age_Scaled + Day_Scaled + (1 | Field), family = nbinom2, data = TaxModel) #model convergence problem
Rich_PxAD <- glmmTMB(Species_Rich ~ Position * Age_Scaled + Day_Scaled + (1 | Field), family = nbinom2, data = TaxModel)
Rich_PxDA <- glmmTMB(Species_Rich ~ Position * Day_Scaled + Age_Scaled + (1 | Field), family = nbinom2, data = TaxModel)
Rich_PAxD <- glmmTMB(Species_Rich ~ Position + Day_Scaled * Age_Scaled + (1 | Field), family = nbinom2, data = TaxModel) #model convergence problems

#collect models
richmodlist <- list("null" = Rich_null, "P" = Rich_P, "A" = Rich_A, 
                    "D" = Rich_D, "PA" = Rich_PA, "PD" = Rich_PD, 
                    "DA" = Rich_DA, "PxD" = Rich_PxD,
                    "PxAD" = Rich_PxAD, "PxDA" = Rich_PxDA)

aictab(richmodlist)
#top model is A

##Step 2 - Environmental Variables----

head(TaxModel)

TaxModel$GC_Scaled <- scale(TaxModel$GC)
TaxModel$FieldNDVI_Scaled <- scale(TaxModel$NDVImean_Field)
TaxModel$FieldArea_Scaled <- scale(TaxModel$Field_Area_m2)

#had to scale some of the variables for the models to converge properly 

Rich_Height <- glmmTMB(Species_Rich ~ Crop_Age_Days + Height + (1 | Field), family = nbinom2, data = TaxModel) 

Rich_GC <- glmmTMB(Species_Rich ~ Crop_Age_Days + GC + (1 | Field), family = nbinom2, data = TaxModel)

Rich_FieldArea <- glmmTMB(Species_Rich ~ Age_Scaled + FieldArea_Scaled + (1 | Field), family = nbinom2, data = TaxModel)
Rich_NDVIfield <- glmmTMB(Species_Rich ~ Crop_Age_Days + FieldNDVI_Scaled + (1 | Field), family = nbinom2, data = TaxModel) 

Rich_Rip <- glmmTMB(Species_Rich ~ Age_Scaled + X1km_Rip_Prop + (1 | Field), family = nbinom2, data = TaxModel) 
Rich_Crops <- glmmTMB(Species_Rich ~ Crop_Age_Days + X1km_Prop_Crops + (1 | Field), family = nbinom2, data = TaxModel) 
Rich_NDVI1km <- glmmTMB(Species_Rich ~ Age_Scaled + NDVIsum_1km + (1 | Field), family = nbinom2, data = TaxModel) 


richmodlist2 <- list("null" = Rich_null, "Height" = Rich_Height,
                     "GC" = Rich_GC, "Field Area" = Rich_FieldArea,
                     "Riparian" = Rich_Rip, 'Crops' = Rich_Crops, 
                     "NDVI 1km" = Rich_NDVI1km, 
                     'Field_NDVI' = Rich_NDVIfield,"Age" = Rich_A)
aictab(richmodlist2)
#top model is Field area
#none within 2 AICc


##Step 3 - Check for spatial autocorrelation----

field_numbers <- unique(TaxModel$ID)

  
model_residuals <- simulateResiduals(Rich_FieldArea)
  
spatial_result <- data.frame(
    field = rep(NA, length(field_numbers)),
    statistic = rep(NA, length(field_numbers)),
    p_value = rep(NA, length(field_numbers)),
    method = rep(NA_character_, length(field_numbers)),
    stringsAsFactors = FALSE)
  
  s <- 1
  
  for (f in field_numbers) {
    
    cat("Field", f, "\n") #What field is the loop doing?
    
    #Extracting specific residuals for individual fields
    field_indices <- which(TaxModel$ID == f)
    field_residuals <- model_residuals
    field_residuals$scaledResiduals <- 
      model_residuals$scaledResiduals[field_indices]
    field_residuals$fittedPredictedResponse <- 
      model_residuals$fittedPredictedResponse[field_indices]
    
    #Testing spatial autocorrelation using XY grid coordinates
    spatial_test <- testSpatialAutocorrelation(field_residuals, 
       x = TaxModel$X_Cor[TaxModel$ID == f], 
       y = TaxModel$Y_Cor[TaxModel$ID == f])
    
    
    spatial_result$field [s] <- f
    spatial_result$statistic [s] <- spatial_test$statistic[1] 
    spatial_result$p_value [s] <- spatial_test$p.value
    spatial_result$method [s] <- spatial_test$method
    
    s <- s + 1
    
  }

head(spatial_result);dim(spatial_result)
length(unique(TaxModel$ID))

#My own notes/explanation of different bits of results
spatial_test$statistic # The test statistic - note the important one here is observed - want this to be close to zero 
spatial_test$p.value # P-value for the test - indicates if there is significant spatial autocorrelation - statistic above indicates the magnitude and the direction (positive or negative)
spatial_test$method #Method used to get statistic


spatial_result
#No Spatial Autocorrelation found in any of the fields/surveys

Spatial_auto_TaxRichTop <- spatial_result

write.xlsx(Spatial_auto_TaxRichTop, 'SpatialResult.xlsx')

##Step 4 - Predictions----


summary(Rich_FieldArea)

Predictions_Age <- seq(min(TaxModel$Age_Scaled),max(TaxModel$Age_Scaled),length.out=20)
Predictions_FieldArea <- seq(min(TaxModel$FieldArea_Scaled),max(TaxModel$FieldArea_Scaled),length.out=20)


richpred <- expand.grid(Age_Scaled = Predictions_Age, FieldArea_Scaled = Predictions_FieldArea)
head(richpred);dim(richpred)

richpred1 <- predict(object = Rich_FieldArea,newdata= richpred,se.fit = T, type = "link",re.form = NA)

richpred2<-data.frame(richpred,fit.link=richpred1$fit,se.link=richpred1$se.fit)

richpred2$lci.link<-richpred2$fit.link-(1.96*richpred2$se.link)
richpred2$uci.link<-richpred2$fit.link+(1.96*richpred2$se.link)

richpred2$fit<-exp(richpred2$fit.link)
richpred2$se<-exp(richpred2$se.link)
richpred2$lci<-exp(richpred2$lci.link)
richpred2$uci<-exp(richpred2$uci.link)

head(richpred2);dim(richpred2)


##Step 5 - Richness Model Visalisation----


summary(Rich_FieldArea)
head(richpred2);dim(richpred2)


RR <- richpred2$FieldArea_Scaled == Predictions_FieldArea[10]
R_R <- richpred2$Age_Scaled == Predictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = TaxModel$Age_Scaled,y = TaxModel$Species_Rich,xlab = "Crop Age (Days)",ylab = 'Taxonomic Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(richpred2$Age_Scaled),to=max(richpred2$Age_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Crop_Age_Days),to=max(TaxModel$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.1,'a)',cex=1.1)

polygon(x = c(richpred2$Age_Scaled[RR],rev(richpred2$Age_Scaled[RR])), y = c(richpred2$lci[RR],rev(richpred2$uci[RR])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred2$Age_Scaled[RR],y = richpred2$fit[RR],lwd = 2,col = 'grey30',lty = 1)

plot(x = TaxModel$FieldArea_Scaled,y = TaxModel$Species_Rich,xlab = expression("Field Size (ha)"),ylab = 'Taxonomic Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.3,'b)',cex=1.1)
axis(side=1, at=seq(from=min(TaxModel$FieldArea_Scaled),to=max(TaxModel$FieldArea_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Field_Area_m2),to=max(TaxModel$Field_Area_m2),length.out=6)/10000,1))

polygon(x = c(richpred2$FieldArea_Scaled[R_R],rev(richpred2$FieldArea_Scaled[R_R])), y = c(richpred2$lci[R_R],rev(richpred2$uci[R_R])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=richpred2$FieldArea_Scaled[R_R],y = richpred2$fit[R_R],lwd = 2,col = 'grey30')


#DIVERSITY MODELLING----

##Step 1 - Design Variables----

#position, age and day - all the various combinations of these 

head(TaxModel);dim(TaxModel)
str(TaxModel)

#update all 0 to be 0.000001 for model fitting
TaxModel$Diversity[TaxModel$Diversity==0] <- 0.000001


Div_null <- glmer(Diversity ~ 1 + (1 | Field), family = Gamma(link = "log"), data = TaxModel)

Div_P <- glmer(Diversity ~ Position + (1 | Field), family = Gamma(link = "log"), data = TaxModel) #Model Covergence issues
Div_A <- glmer(Diversity ~ Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 
Div_D <- glmer(Diversity ~ Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 

Div_PA <- glmer(Diversity ~ Position + Age_Scaled+ (1 | Field), family = Gamma(link = "log"), data = TaxModel) 
Div_PD <- glmer(Diversity ~ Position + Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) #Model Convergence issues
Div_DA <- glmer(Diversity ~ Day_Scaled + Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel)

Div_PxA <- glmer(Diversity ~ Position * Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 
Div_PxD <- glmer(Diversity ~ Position * Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 
Div_DxA <- glmer(Diversity ~ Day_Scaled * Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel)

Div_PAD <- glmer(Diversity ~ Position + Age_Scaled + Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_PxAD <- glmer(Diversity ~ Position * Age_Scaled + Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_PxDA <- glmer(Diversity ~ Position * Day_Scaled + Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 
Div_PAxD <- glmer(Diversity ~ Position + Age_Scaled * Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 


divlist <- list("null" = Div_null, "A" = Div_A,"D" = Div_D, 
                "PA" = Div_PA, "DA" = Div_DA, "PxA" = Div_PxA, 
                "PxD" = Div_PxD, "DxA" = Div_DxA, "PAD" = Div_PAD,
                "PxAD" = Div_PxAD, 'PxDA' = Div_PxDA, 
                'PAxD' = Div_PAxD)

aictab(divlist)
#Day x Age is top model

##Step 2 - Environmental Variables----

head(TaxModel)

TaxModel$Field_Area_Scaled <- scale(TaxModel$Field_Area_m2)
TaxModel$NDVI1km_Scaled <- scale(TaxModel$NDVIsum_1km)

Div_Height <- glmer(Diversity ~ Day_Scaled * Age_Scaled + Height + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_GC <- glmer(Diversity ~ Day_Scaled * Age_Scaled + GC + (1 | Field), family = Gamma(link = "log"), data = TaxModel)

Div_FieldArea <- glmer(Diversity ~ Day_Scaled * Age_Scaled + FieldArea_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_FieldNDVI <- glmer(Diversity ~ Day_Scaled * Age_Scaled + NDVImean_Field + (1 | Field), family = Gamma(link = "log"), data = TaxModel)

Div_Rip <- glmer(Diversity ~ Day_Scaled * Age_Scaled + X1km_Rip_Prop + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_Crop <- glmer(Diversity ~ Day_Scaled * Age_Scaled + X1km_Prop_Crops + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_NDVI1km <- glmer(Diversity ~ Day_Scaled * Age_Scaled + NDVIsum_1km  + (1 | Field), family = Gamma(link = "log"), data = TaxModel)

divlist2 <- list("null" = Div_null, "height" = Div_Height, 
                "GC" = Div_GC, "Field Area" = Div_FieldArea, 
                "Field NDVI" = Div_FieldNDVI, "Rip" = Div_Rip, 
                "Crops" = Div_Crop, "NDVI 1km" = Div_NDVI1km,
                "DxA" = Div_DxA)
aictab(divlist2)

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


##Step 3 - Check for spatial autocorrelation----

model_residuals <- simulateResiduals(Div_DxA)
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
  field_indices <- which(TaxModel$ID == f)
  field_residuals <- model_residuals
  field_residuals$scaledResiduals <- 
    model_residuals$scaledResiduals[field_indices]
  field_residuals$fittedPredictedResponse <- 
    model_residuals$fittedPredictedResponse[field_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  spatial_test <- testSpatialAutocorrelation(field_residuals, 
                                             x = TaxModel$X_Cor[TaxModel$ID == f], 
                                             y = TaxModel$Y_Cor[TaxModel$ID == f])
  
  
  spatial_result$field [s] <- f
  spatial_result$statistic [s] <- spatial_test$statistic[1] 
  spatial_result$p_value [s] <- spatial_test$p.value
  spatial_result$method [s] <- spatial_test$method
  
  s <- s + 1
  
}

head(spatial_result);dim(spatial_result)
length(unique(TaxModel$ID))

spatial_result

#saved results

Spatial_auto_TaxDiv <- spatial_result 


#Spaital autocorrelation found in two fields but statistic is negligable
write.xlsx(Spatial_auto_TaxDiv, 'SpatialResult.xlsx')

##Step 4 - Predictions----

summary(Div_DxA)

Divpred <- expand.grid(Age_Scaled = Predictions_Age, Day_Scaled = Predictions_Day)
head(Divpred);dim(Divpred)

Divpred1 <- predict(object = Div_DxA,newdata= Divpred,se.fit = T, type = "link",re.form = NA)

Divpred2<-data.frame(Divpred,fit.link=Divpred1$fit,se.link=Divpred1$se.fit)

Divpred2$lci.link<-Divpred2$fit.link-(1.96*Divpred2$se.link)
Divpred2$uci.link<-Divpred2$fit.link+(1.96*Divpred2$se.link)

Divpred2$fit<-exp(Divpred2$fit.link)
Divpred2$se<-exp(Divpred2$se.link)
Divpred2$lci<-exp(Divpred2$lci.link)
Divpred2$uci<-exp(Divpred2$uci.link)

head(Divpred2);dim(Divpred2)


##Step 5 - Model Visalisation----


summary(Div_DxA)
head(Divpred2);dim(Divpred2)


YY <- Divpred2$Day_Scaled == Predictions_Day[4] #Winter
Y_Y <- Divpred2$Day_Scaled == Predictions_Day[15] #Spring

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = TaxModel$Age_Scaled,y = TaxModel$Species_Rich,xlab = "Crop Age (Days)",ylab = 'Taxonomic Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(richpred2$Age_Scaled),to=max(richpred2$Age_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Crop_Age_Days),to=max(TaxModel$Crop_Age_Days),length.out=6),-1))

polygon(x = c(Divpred2$Age_Scaled[YY],rev(Divpred2$Age_Scaled[YY])), y = c(Divpred2$lci[YY],rev(Divpred2$uci[YY])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred2$Age_Scaled[YY],y = Divpred2$fit[YY],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(Divpred2$Age_Scaled[Y_Y],rev(Divpred2$Age_Scaled[Y_Y])), y = c(Divpred2$lci[Y_Y],rev(Divpred2$uci[Y_Y])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred2$Age_Scaled[Y_Y],y = Divpred2$fit[Y_Y],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#END----
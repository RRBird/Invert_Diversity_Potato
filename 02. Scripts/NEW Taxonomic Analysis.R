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
#top model is D*A (both scaled)

##Step 2 - Environmental Variables----

head(TaxModel)

TaxModel$GC_Scaled <- scale(TaxModel$GC)
TaxModel$FieldNDVI_Scaled <- scale(TaxModel$NDVImean_Field)
TaxModel$FieldArea_Scaled <- scale(TaxModel$Field_Area_m2)

#had to scale a lot of variables for the models to converge properly 

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
                     'Field_NDVI' = Rich_NDVIfield)
aictab(richmodlist2)
#top model is Field area
#NDVI 1km is within 4 AICc none within 2 AICc


##Step 3 - Check for spatial autocorrelation----

#Field area model

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

#My own notes/explanation of different bits of results
spatial_test$statistic # The test statistic - note the important one here is observed - want this to be close to zero 
spatial_test$p.value # P-value for the test - indicates if there is significant spatial autocorrelation - statistic above indicates the magnitude and the direction (positive or negative)
spatial_test$method #Method used to get statistic


spatial_result
#No Spatial Autocorrelation found in any of the fields/surveys

Spatial_auto_TaxRichTop <- spatial_result

#NDVI 1km


model_residuals <- simulateResiduals(Rich_NDVI1km)
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
#No Spatial Autocorrelation found in any of the fields/surveys

Spatial_auto_TaxRichSI <- spatial_result


##Step 4 - Predictions----

#Field Area

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

#NDVI 1km

summary(Rich_NDVI1km)

Predictions_NDVI1km <- seq(min(TaxModel$NDVIsum_1km),max(TaxModel$NDVIsum_1km),length.out=20)


richpred3 <- expand.grid(Age_Scaled = Predictions_Age, NDVIsum_1km = Predictions_NDVI1km)
head(richpred3);dim(richpred3)

richpred4 <- predict(object = Rich_NDVI1km,newdata= richpred3,se.fit = T, type = "link",re.form = NA)

richpred5<-data.frame(richpred3,fit.link=richpred4$fit,se.link=richpred4$se.fit)

richpred5$lci.link<-richpred5$fit.link-(1.96*richpred5$se.link)
richpred5$uci.link<-richpred5$fit.link+(1.96*richpred5$se.link)

richpred5$fit<-exp(richpred5$fit.link)
richpred5$se<-exp(richpred5$se.link)
richpred5$lci<-exp(richpred5$lci.link)
richpred5$uci<-exp(richpred5$uci.link)

head(richpred5);dim(richpred5)


##Step 5 - Richness Model Visalisation----

#Field Area

summary(Rich_FieldArea)
head(richpred2);dim(richpred2)


RR <- richpred2$FieldArea_Scaled == Predictions_FieldArea[10]
R_R <- richpred2$Age_Scaled == Predictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = TaxModel$Age_Scaled,y = TaxModel$Species_Rich,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(richpred2$Age_Scaled),to=max(richpred2$Age_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Crop_Age_Days),to=max(TaxModel$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.1,'a)',cex=1.1)

polygon(x = c(richpred2$Age_Scaled[RR],rev(richpred2$Age_Scaled[RR])), y = c(richpred2$lci[RR],rev(richpred2$uci[RR])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred2$Age_Scaled[RR],y = richpred2$fit[RR],lwd = 2,col = 'grey30',lty = 1)

plot(x = TaxModel$FieldArea_Scaled,y = TaxModel$Species_Rich,xlab = expression("Field Size (ha)"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.3,'b)',cex=1.1)
axis(side=1, at=seq(from=min(TaxModel$FieldArea_Scaled),to=max(TaxModel$FieldArea_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Field_Area_m2),to=max(TaxModel$Field_Area_m2),length.out=6)/10000,1))

polygon(x = c(richpred2$FieldArea_Scaled[R_R],rev(richpred2$FieldArea_Scaled[R_R])), y = c(richpred2$lci[R_R],rev(richpred2$uci[R_R])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=richpred2$FieldArea_Scaled[R_R],y = richpred2$fit[R_R],lwd = 2,col = 'grey30')

#NDVI 1km

summary(Rich_NDVI1km)
head(richpred5);dim(richpred5)


VV <- richpred5$NDVIsum_1km == Predictions_1kmNDVI[10]
V_V <- richpred5$Age_Scaled == Predictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = TaxModel$Age_Scaled,y = TaxModel$Species_Rich,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(richpred5$Age_Scaled),to=max(richpred5$Age_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Crop_Age_Days),to=max(TaxModel$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.1,'a)',cex=1.1)

polygon(x = c(richpred5$Age_Scaled[VV],rev(richpred5$Age_Scaled[VV])), y = c(richpred5$lci[VV],rev(richpred5$uci[VV])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred5$Age_Scaled[VV],y = richpred5$fit[VV],lwd = 2,col = 'grey30',lty = 1)

plot(x = TaxModel$NDVIsum_1km,y = TaxModel$Species_Rich,xlab = expression("NDVI within 1km"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(richpred5$NDVIsum_1km),to=max(richpred5$NDVIsum_1km),length.out=6),labels=round(seq(from=min(TaxModel$NDVIsum_1km),to=max(TaxModel$NDVIsum_1km),length.out=6),-1))
mtext(side=3,line=0,at = 370,'b)',cex=1.1)

polygon(x = c(richpred5$NDVIsum_1km[V_V],rev(richpred5$NDVIsum_1km[V_V])), y = c(richpred5$lci[V_V],rev(richpred5$uci[V_V])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=richpred5$NDVIsum_1km[V_V],y = richpred5$fit[V_V],lwd = 2,col = 'grey30')


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
                "Crops" = Div_Crop, "NDVI 1km" = Div_NDVI1km)
aictab(divlist2)

#top model is crops and all the rest are within 2 AICc

##Step 3 - Check for spatial autocorrelation----

model_residuals <- simulateResiduals(Div_FieldArea)
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

Spatial_auto_TaxDivCrop <- spatial_result 
Spatial_auto_TaxDivNDVIfield <- spatial_result 
Spatial_auto_TaxDivNDVI1km <- spatial_result 
Spatial_auto_TaxDivGC <- spatial_result 
Spatial_auto_TaxDivFieldArea <- spatial_result 
Spatial_auto_TaxDivRip <- spatial_result 
Spatial_auto_TaxDivHeight <- spatial_result 

#Spaital autocorrelation found in all models
write.xlsx(Spatial_auto_TaxDivHeight, 'SpatialResult.xlsx')

#END----
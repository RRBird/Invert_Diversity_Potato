
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

length(which(TaxModel$Species_Rich==0))/
  length(TaxModel$Species_Rich)

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

plot(x = TaxModel$FieldArea_Scaled,y = TaxModel$Species_Rich,xlab = expression("Field Size (ha)"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.3,'a)',cex=1.1)
axis(side=1, at=seq(from=min(TaxModel$FieldArea_Scaled),to=max(TaxModel$FieldArea_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Field_Area_m2),to=max(TaxModel$Field_Area_m2),length.out=6)/10000,1))

polygon(x = c(richpred2$FieldArea_Scaled[R_R],rev(richpred2$FieldArea_Scaled[R_R])), y = c(richpred2$lci[R_R],rev(richpred2$uci[R_R])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=richpred2$FieldArea_Scaled[R_R],y = richpred2$fit[R_R],lwd = 2,col = 'grey30')

plot(x = TaxModel$Age_Scaled,y = TaxModel$Species_Rich,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(richpred2$Age_Scaled),to=max(richpred2$Age_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Crop_Age_Days),to=max(TaxModel$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.1,'b)',cex=1.1)

polygon(x = c(richpred2$Age_Scaled[RR],rev(richpred2$Age_Scaled[RR])), y = c(richpred2$lci[RR],rev(richpred2$uci[RR])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred2$Age_Scaled[RR],y = richpred2$fit[RR],lwd = 2,col = 'grey30',lty = 1)


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


##Step 3 - Check for spatial autocorrelation----

summary(Div_DxA)
summary(Div_Crop)
summary(Div_FieldNDVI)
summary(Div_NDVI1km)
summary(Div_GC)


model_residuals <- simulateResiduals(Div_GC)
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
Spatial_auto_TaxDivCrop <- spatial_result 
Spatial_auto_TaxDivNDVIfield <- spatial_result
Spatial_auto_TaxDivNDVI1km <- spatial_result
Spatial_auto_TaxDivGC <- spatial_result

#Spatial autocorrelation found in two fields but statistic is negligable
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


summary(Div_Crop)
Div_Predictions_Crops <- seq(min(TaxModel$X1km_Prop_Crops),max(TaxModel$X1km_Prop_Crops),length.out=20)


Divpred3 <- expand.grid(Age_Scaled = Predictions_Age, Day_Scaled = Predictions_Day,X1km_Prop_Crops = Div_Predictions_Crops)
head(Divpred3);dim(Divpred3)

Divpred4 <- predict(object = Div_Crop,newdata= Divpred3,se.fit = T, type = "link",re.form = NA)

Divpred5<-data.frame(Divpred3,fit.link=Divpred4$fit,se.link=Divpred4$se.fit)

Divpred5$lci.link<-Divpred5$fit.link-(1.96*Divpred5$se.link)
Divpred5$uci.link<-Divpred5$fit.link+(1.96*Divpred5$se.link)

Divpred5$fit<-exp(Divpred5$fit.link)
Divpred5$se<-exp(Divpred5$se.link)
Divpred5$lci<-exp(Divpred5$lci.link)
Divpred5$uci<-exp(Divpred5$uci.link)

head(Divpred5);dim(Divpred5)


summary(Div_FieldNDVI)
Div_Predictions_FieldNDVI <- seq(min(TaxModel$NDVImean_Field),max(TaxModel$NDVImean_Field),length.out=20)


Divpred6 <- expand.grid(Age_Scaled = Predictions_Age, Day_Scaled = Predictions_Day,NDVImean_Field = Div_Predictions_FieldNDVI)
head(Divpred6);dim(Divpred6)

Divpred7 <- predict(object = Div_FieldNDVI,newdata= Divpred6,se.fit = T, type = "link",re.form = NA)

Divpred8<-data.frame(Divpred6,fit.link=Divpred7$fit,se.link=Divpred7$se.fit)

Divpred8$lci.link<-Divpred8$fit.link-(1.96*Divpred8$se.link)
Divpred8$uci.link<-Divpred8$fit.link+(1.96*Divpred8$se.link)

Divpred8$fit<-exp(Divpred8$fit.link)
Divpred8$se<-exp(Divpred8$se.link)
Divpred8$lci<-exp(Divpred8$lci.link)
Divpred8$uci<-exp(Divpred8$uci.link)

head(Divpred8);dim(Divpred8)


summary(Div_NDVI1km)
Div_Predictions_NDVI1km <- seq(min(TaxModel$NDVIsum_1km),max(TaxModel$NDVIsum_1km),length.out=20)


Divpred9 <- expand.grid(Age_Scaled = Predictions_Age, Day_Scaled = Predictions_Day,NDVIsum_1km = Div_Predictions_NDVI1km)
head(Divpred9);dim(Divpred9)

Divpred10 <- predict(object = Div_NDVI1km,newdata= Divpred9,se.fit = T, type = "link",re.form = NA)

Divpred11<-data.frame(Divpred9,fit.link=Divpred10$fit,se.link=Divpred10$se.fit)

Divpred11$lci.link<-Divpred11$fit.link-(1.96*Divpred11$se.link)
Divpred11$uci.link<-Divpred11$fit.link+(1.96*Divpred11$se.link)

Divpred11$fit<-exp(Divpred11$fit.link)
Divpred11$se<-exp(Divpred11$se.link)
Divpred11$lci<-exp(Divpred11$lci.link)
Divpred11$uci<-exp(Divpred11$uci.link)

head(Divpred11);dim(Divpred11)



summary(Div_GC)

Div_Predictions_GC <- seq(min(TaxModel$GC),max(TaxModel$GC),length.out=20)


Divpred12 <- expand.grid(Age_Scaled = Predictions_Age, Day_Scaled = Predictions_Day,GC = Div_Predictions_GC)
head(Divpred12);dim(Divpred12)

Divpred13 <- predict(object = Div_GC,newdata= Divpred12,se.fit = T, type = "link",re.form = NA)

Divpred14<-data.frame(Divpred12,fit.link=Divpred13$fit,se.link=Divpred13$se.fit)

Divpred14$lci.link<-Divpred14$fit.link-(1.96*Divpred14$se.link)
Divpred14$uci.link<-Divpred14$fit.link+(1.96*Divpred14$se.link)

Divpred14$fit<-exp(Divpred14$fit.link)
Divpred14$se<-exp(Divpred14$se.link)
Divpred14$lci<-exp(Divpred14$lci.link)
Divpred14$uci<-exp(Divpred14$uci.link)

head(Divpred14);dim(Divpred14)


##Step 5 - Model Visalisation----


summary(Div_DxA)
head(Divpred2);dim(Divpred2)


YY <- Divpred2$Day_Scaled == Predictions_Day[4] #Winter
Y_Y <- Divpred2$Day_Scaled == Predictions_Day[15] #Spring

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = TaxModel$Age_Scaled,y = TaxModel$Diversity,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(Divpred2$Age_Scaled),to=max(Divpred2$Age_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Crop_Age_Days),to=max(TaxModel$Crop_Age_Days),length.out=6),-1))

polygon(x = c(Divpred2$Age_Scaled[YY],rev(Divpred2$Age_Scaled[YY])), y = c(Divpred2$lci[YY],rev(Divpred2$uci[YY])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred2$Age_Scaled[YY],y = Divpred2$fit[YY],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(Divpred2$Age_Scaled[Y_Y],rev(Divpred2$Age_Scaled[Y_Y])), y = c(Divpred2$lci[Y_Y],rev(Divpred2$uci[Y_Y])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred2$Age_Scaled[Y_Y],y = Divpred2$fit[Y_Y],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)


##Supporting Info Models----

ZZ <- Divpred5$Day_Scaled == Predictions_Day[10] & Divpred5$Age_Scaled == Predictions_Age[10]
ZZZ <- Divpred8$Day_Scaled == Predictions_Day[10] & Divpred8$Age_Scaled == Predictions_Age[10]
Z_Z <- Divpred11$Day_Scaled == Predictions_Day[10] & Divpred11$Age_Scaled == Predictions_Age[10]
ZZ_ZZ <- Divpred14$Day_Scaled == Predictions_Day[10] & Divpred14$Age_Scaled == Predictions_Age[10]



dev.new(height=10,width=15,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,3),mgp=c(2.5,1,0),xpd = T)

plot(x = TaxModel$Age_Scaled,y = TaxModel$Diversity,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',cex.lab= 1.3,cex.axis=1.2)
axis(side=1, at=seq(from=min(Divpred2$Age_Scaled),to=max(Divpred2$Age_Scaled),length.out=6),labels=round(seq(from=min(TaxModel$Crop_Age_Days),to=max(TaxModel$Crop_Age_Days),length.out=6),-1),cex.axis=1.2)
mtext(side=3,line=0,at = -2.1,'a)',cex=1)

polygon(x = c(Divpred2$Age_Scaled[YY],rev(Divpred2$Age_Scaled[YY])), y = c(Divpred2$lci[YY],rev(Divpred2$uci[YY])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred2$Age_Scaled[YY],y = Divpred2$fit[YY],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(Divpred2$Age_Scaled[Y_Y],rev(Divpred2$Age_Scaled[Y_Y])), y = c(Divpred2$lci[Y_Y],rev(Divpred2$uci[Y_Y])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred2$Age_Scaled[Y_Y],y = Divpred2$fit[Y_Y],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1.3)


plot(x = TaxModel$X1km_Prop_Crops,y = TaxModel$Diversity,xlab = "Crops within 1km",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.lab= 1.3,cex.axis=1.2)
mtext(side=3,line=0,at = 79.2,'b)',cex=1)

polygon(x = c(Divpred5$X1km_Prop_Crops[ZZ],rev(Divpred5$X1km_Prop_Crops[ZZ])), y = c(Divpred5$lci[ZZ],rev(Divpred5$uci[ZZ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred5$X1km_Prop_Crops[ZZ],y = Divpred5$fit[ZZ],lwd = 2,col = 'grey30',lty = 1)


plot(x = TaxModel$NDVImean_Field,y = TaxModel$Diversity,xlab = "Field NDVI",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.lab= 1.3,cex.axis=1.2,xaxt = 'n')
axis(side=1, at=seq(from=min(Divpred8$NDVImean_Field),to=max(Divpred8$NDVImean_Field),length.out=5),labels=round(seq(from=min(TaxModel$NDVImean_Field),to=max(TaxModel$NDVImean_Field),length.out=5),2),cex.axis=1.2)
mtext(side=3,line=0,at = 0.12,'c)',cex=1)

polygon(x = c(Divpred8$NDVImean_Field[ZZZ],rev(Divpred8$NDVImean_Field[ZZZ])), y = c(Divpred8$lci[ZZZ],rev(Divpred8$uci[ZZZ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred8$NDVImean_Field[ZZZ],y = Divpred8$fit[ZZZ],lwd = 2,col = 'grey30',lty = 1)


plot(x = TaxModel$NDVIsum_1km,y = TaxModel$Diversity,xlab = "NDVI within 1km",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.lab= 1.3,cex.axis=1.2,xaxt = 'n')
axis(side=1, at=seq(from=min(Divpred11$NDVIsum_1km),to=max(Divpred11$NDVIsum_1km),length.out=5),labels=round(seq(from=min(TaxModel$NDVIsum_1km),to=max(TaxModel$NDVIsum_1km),length.out=5),-1),cex.axis=1.2)
mtext(side=3,line=0,at = 370,'d)',cex=1)

polygon(x = c(Divpred11$NDVIsum_1km[Z_Z],rev(Divpred11$NDVIsum_1km[Z_Z])), y = c(Divpred11$lci[Z_Z],rev(Divpred11$uci[Z_Z])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred11$NDVIsum_1km[Z_Z],y = Divpred11$fit[Z_Z],lwd = 2,col = 'grey30',lty = 1)


plot(x = TaxModel$GC,y = TaxModel$Diversity,xlab = "Ground Cover (%)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.lab= 1.3,cex.axis=1.2)
mtext(side=3,line=0,at = 3,'e)',cex=1)

polygon(x = c(Divpred14$GC[ZZ_ZZ],rev(Divpred14$GC[ZZ_ZZ])), y = c(Divpred14$lci[ZZ_ZZ],rev(Divpred14$uci[ZZ_ZZ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Divpred14$GC[ZZ_ZZ],y = Divpred14$fit[ZZ_ZZ],lwd = 2,col = 'grey30',lty = 1)





#END----
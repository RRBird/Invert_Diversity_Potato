
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

Div_Crop_cf <- summary(Div_Crop)$coefficients

Div_Crop_cf <- data.frame(term = row.names(Div_Crop_cf),Div_Crop_cf)
row.names(Div_Crop_cf) <- 1:nrow(Div_Crop_cf)

colnames(Div_Crop_cf) <- c("term","est",'std_err',"t","p")

Div_Crop_cf$lci <- Div_Crop_cf$est-(1.96*Div_Crop_cf$std_err)
Div_Crop_cf$uci <- Div_Crop_cf$est+(1.96*Div_Crop_cf$std_err)

head(Div_Crop_cf)

Div_Crop_cf$term <- c("Intercept", "Day", "Age", "Crops in 1km","Day:Age")


summary(Div_FieldNDVI)
Div_FieldNDVI_cf <- summary(Div_FieldNDVI)$coefficients

Div_FieldNDVI_cf <- data.frame(term = row.names(Div_FieldNDVI_cf),Div_FieldNDVI_cf)
row.names(Div_FieldNDVI_cf) <- 1:nrow(Div_FieldNDVI_cf)

colnames(Div_FieldNDVI_cf) <- c("term","est",'std_err',"z","p")

Div_FieldNDVI_cf$lci <- Div_FieldNDVI_cf$est-(1.96*Div_FieldNDVI_cf$std_err)
Div_FieldNDVI_cf$uci <- Div_FieldNDVI_cf$est+(1.96*Div_FieldNDVI_cf$std_err)

head(Div_FieldNDVI_cf)

Div_FieldNDVI_cf$term <- c("Intercept", "Day", "Age", "Field NDVI","Day:Age")


summary(Div_NDVI1km)
Div_NDVI1km_cf <- summary(Div_NDVI1km)$coefficients

Div_NDVI1km_cf <- data.frame(term = row.names(Div_NDVI1km_cf),Div_NDVI1km_cf)
row.names(Div_NDVI1km_cf) <- 1:nrow(Div_NDVI1km_cf)

colnames(Div_NDVI1km_cf) <- c("term","est",'std_err',"t","p")

Div_NDVI1km_cf$lci <- Div_NDVI1km_cf$est-(1.96*Div_NDVI1km_cf$std_err)
Div_NDVI1km_cf$uci <- Div_NDVI1km_cf$est+(1.96*Div_NDVI1km_cf$std_err)

head(Div_NDVI1km_cf)

Div_NDVI1km_cf$term <- c("Intercept", "Day", "Age", "NDVI in 1km","Day:Age")


summary(Div_GC)

Div_GC_cf <- summary(Div_GC)$coefficients

Div_GC_cf <- data.frame(term = row.names(Div_GC_cf),Div_GC_cf)
row.names(Div_GC_cf) <- 1:nrow(Div_GC_cf)

colnames(Div_GC_cf) <- c("term","est",'std_err',"t","p")

Div_GC_cf$lci <- Div_GC_cf$est-(1.96*Div_GC_cf$std_err)
Div_GC_cf$uci <- Div_GC_cf$est+(1.96*Div_GC_cf$std_err)

head(Div_GC_cf)

Div_GC_cf$term <- c("Intercept", "Day", "Age", "Ground Cover","Day:Age")


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
Div_Crop_cf
Div_FieldNDVI_cf
Div_NDVI1km_cf
Div_GC_cf


dev.new(height=8,width=10,dpi=60,pointsize=14,noRStudioGD = T)
par(mar=c(6,10,1,1),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = F)

plot(x=rev(Div_Crop_cf$est),y=1:length(Div_Crop_cf$term),ylab="",xlab="Model Estimate",pch=16,yaxt="n",xlim = c(min(Div_Crop_cf$lci),max(Div_Crop_cf$uci)))
axis(2, at=1:length(Div_Crop_cf$term),labels = rev(Div_Crop_cf$term),las=1)
arrows(x0 = 0,y0 = 0.5,x1 = 0,y1 = 5.5,length = 0)
arrows(x0 = rev(Div_Crop_cf$lci),y0 = 1:length(Div_Crop_cf$term),x1 = rev(Div_Crop_cf$uci),y1 = 1:length(Div_Crop_cf$term),length = 0)
mtext("a)",side = 3,line = -0.5,at=-8,cex = 0.9)

plot(x=rev(Div_FieldNDVI_cf$est),y=1:length(Div_FieldNDVI_cf$term),ylab="",xlab="Model Estimate",pch=16,yaxt="n",xlim = c(min(Div_FieldNDVI_cf$lci),max(Div_FieldNDVI_cf$uci)))
axis(2, at=1:length(Div_FieldNDVI_cf$term),labels = rev(Div_FieldNDVI_cf$term),las=1)
arrows(x0 = 0,y0 = 0.5,x1 = 0,y1 = 5.5,length = 0)
arrows(x0 = rev(Div_FieldNDVI_cf$lci),y0 = 1:length(Div_FieldNDVI_cf$term),x1 = rev(Div_FieldNDVI_cf$uci),y1 = 1:length(Div_FieldNDVI_cf$term),length = 0)
mtext("b)",side = 3,line = -0.5,at=-7,cex = 0.9)

plot(x=rev(Div_NDVI1km_cf$est),y=1:length(Div_NDVI1km_cf$term),ylab="",xlab="Model Estimate",pch=16,yaxt="n",xlim = c(min(Div_NDVI1km_cf$lci),max(Div_NDVI1km_cf$uci)))
axis(2, at=1:length(Div_NDVI1km_cf$term),labels = rev(Div_NDVI1km_cf$term),las=1)
arrows(x0 = 0,y0 = 0.5,x1 = 0,y1 = 5.5,length = 0)
arrows(x0 = rev(Div_NDVI1km_cf$lci),y0 = 1:length(Div_NDVI1km_cf$term),x1 = rev(Div_NDVI1km_cf$uci),y1 = 1:length(Div_NDVI1km_cf$term),length = 0)
mtext("0.0",side = 1,line = 1,at=0,cex = 0.95)
mtext("1.0",side = 1,line = 1,at=1,cex = 0.95)
mtext("c)",side = 3,line = -0.5,at=-2.5,cex = 0.9)

plot(x=rev(Div_GC_cf$est),y=1:length(Div_GC_cf$term),ylab="",xlab="Model Estimate",pch=16,yaxt="n",xlim = c(min(Div_GC_cf$lci),max(Div_GC_cf$uci)))
axis(2, at=1:length(Div_GC_cf$term),labels = rev(Div_GC_cf$term),las=1)
arrows(x0 = 0,y0 = 0.5,x1 = 0,y1 = 5.5,length = 0)
arrows(x0 = rev(Div_GC_cf$lci),y0 = 1:length(Div_GC_cf$term),x1 = rev(Div_GC_cf$uci),y1 = 1:length(Div_GC_cf$term),length = 0)
mtext("d)",side = 3,line = -0.5,at=-2.3,cex = 0.9)

#END----
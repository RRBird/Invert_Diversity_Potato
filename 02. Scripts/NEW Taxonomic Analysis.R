#This script contains the analysis and modelling of taxpnomic richness and diversity

#Libraries----
library("AICcmodavg")
library('lme4')
library("glmmTMB")
library("dplyr")
library("DHARMa")
library("arm")

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
Rich_D <- glmmTMB(Species_Rich ~ Day_Sampled + (1 | Field), family = nbinom2, data = TaxModel)

Rich_PA <- glmmTMB(Species_Rich ~ Position + Crop_Age_Days + (1 | Field), family = nbinom2, data = TaxModel)
Rich_PD <- glmmTMB(Species_Rich ~ Position + Day_Sampled + (1 | Field), family = nbinom2, data = TaxModel)
Rich_DA <- glmmTMB(Species_Rich ~ Day_Sampled + Crop_Age_Days + (1 | Field), family = nbinom2, data = TaxModel)

Rich_PxA <- glmmTMB(Species_Rich ~ Position * Age_Scaled + (1 | Field), family = nbinom2, data = TaxModel) #Model Covergence Issue
Rich_PxD <- glmmTMB(Species_Rich ~ Position * Day_Scaled + (1 | Field), family = nbinom2, data = TaxModel) #warning message but solved this by using scaled Day
Rich_DxA <- glmmTMB(Species_Rich ~ Day_Sampled * Crop_Age_Days + (1 | Field), family = nbinom2, data = TaxModel) 

Rich_PAD <- glmmTMB(Species_Rich ~ Position + Age_Scaled + Day_Scaled + (1 | Field), family = nbinom2, data = TaxModel) #model covergence problem
Rich_PxAD <- glmmTMB(Species_Rich ~ Position * Crop_Age_Days + Day_Sampled + (1 | Field), family = nbinom2, data = TaxModel)
Rich_PxDA <- glmmTMB(Species_Rich ~ Position * Day_Sampled + Crop_Age_Days + (1 | Field), family = nbinom2, data = TaxModel)
Rich_PAxD <- glmmTMB(Species_Rich ~ Position + Day_Sampled * Crop_Age_Days + (1 | Field), family = nbinom2, data = TaxModel)

#collect models
richmodlist <- list("null" = Rich_null, "P" = Rich_P, "A" = Rich_A, 
                    "D" = Rich_D, "PA" = Rich_PA, "PD" = Rich_PD, 
                    "DA" = Rich_DA, "PxD" = Rich_PxD,"DxA" = Rich_DxA,
                    "PxAD" = Rich_PxAD, "PxDA" = Rich_PxDA, 
                    "PAxD" = Rich_PAxD)

aictab(richmodlist)
#top model is D*A (both scaled)

##Step 2 - Environmental Variables----

head(TaxModel)

TaxModel$Height_Scaled <- scale(TaxModel$Height)
TaxModel$GC_Scaled <- scale(TaxModel$GC)
TaxModel$FieldNDVI_Scaled <- scale(TaxModel$NDVImean_Field)

#had to scale a lot of variables for the models to converge properly 

Rich_Height <- glmmTMB(Species_Rich ~ Day_Scaled * Age_Scaled + Height_Scaled + (1 | Field), family = nbinom2, data = TaxModel) #warning message but scaled height fixed it
Rich_GC <- glmmTMB(Species_Rich ~ Day_Scaled * Age_Scaled + GC_Scaled + (1 | Field), family = nbinom2, data = TaxModel) #warning message

Rich_FieldArea <- glmmTMB(Species_Rich ~ Day_Scaled * Age_Scaled + Field_Area_m2 + (1 | Field), family = nbinom2, data = TaxModel) 
Rich_NDVIfield <- glmmTMB(Species_Rich ~ Day_Sampled * Crop_Age_Days + FieldNDVI_Scaled + (1 | Field), family = nbinom2, data = TaxModel)  #Model Convergence issue

Rich_Water <- glmmTMB(Species_Rich ~ Day_Scaled * Age_Scaled + X1km_Prop_Water + (1 | Field), family = nbinom2, data = TaxModel) 
Rich_NDVI1km <- glmmTMB(Species_Rich ~ Day_Scaled * Age_Scaled + NDVIsum_1km + (1 | Field), family = nbinom2, data = TaxModel) 


richmodlist2 <- list("null" = Rich_null, "Height" = Rich_Height,
                     "GC" = Rich_GC, "Field Area" = Rich_FieldArea,
                     "Water" = Rich_Water, "NDVI 1km" = Rich_NDVI1km)
aictab(richmodlist2)
#top model is water
#none within 2 or 4 AICc


##Step 3 - Check for spatial autocorrelation----

field_numbers <- unique(TaxModel$ID)

  
  model_residuals <- simulateResiduals(Rich_Water)
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


##Step 4 - Predictions----

summary(Rich_Water)

Predictions_Day <- seq(min(TaxModel$Day_Scaled),max(TaxModel$Day_Scaled),length.out=20) 
Predictions_Age <- seq(min(TaxModel$Age_Scaled),max(TaxModel$Age_Scaled),length.out=20)
Predictions_Water <- seq(min(TaxModel$X1km_Prop_Water),max(TaxModel$X1km_Prop_Water),length.out=20)


richpred <- expand.grid(Day_Scaled = Predictions_Day, Age_Scaled = Predictions_Age, X1km_Prop_Water = Predictions_Water)
head(richpred);dim(richpred)

richpred1 <- predict(object = Rich_Water,newdata= richpred,se.fit = T, type = "link",re.form = NA)

richpred2<-data.frame(richpred,fit.link=richpred1$fit,se.link=richpred1$se.fit)

richpred2$lci.link<-richpred2$fit.link-(1.96*richpred2$se.link)
richpred2$uci.link<-richpred2$fit.link+(1.96*richpred2$se.link)

richpred2$fit<-exp(richpred2$fit.link)
richpred2$se<-exp(richpred2$se.link)
richpred2$lci<-exp(richpred2$lci.link)
richpred2$uci<-exp(richpred2$uci.link)

head(richpred2);dim(richpred2)


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


#DIVERSITY MODELLING----

##Step 1 - Design Variables----

#position, age and day - all the various combinations of these 

head(TaxModel);dim(TaxModel)
str(TaxModel)

#update all 0 to be 0.000001 for model fitting
TaxModel$Diversity[TaxModel$Diversity==0] <- 0.000001


Div_null <- glmer(Diversity ~ 1 + (1 | Field), family = Gamma(link = "log"), data = TaxModel)

Div_P <- glmer(Diversity ~ Position + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_A <- glmer(Diversity ~ Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 
Div_D <- glmer(Diversity ~ Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 

Div_PA <- glmer(Diversity ~ Position + Age_Scaled+ (1 | Field), family = Gamma(link = "log"), data = TaxModel) #tried scale, more iterations and changing optimiser and model still failed to converge
Div_PD <- glmer(Diversity ~ Position + Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 
Div_DA <- glmer(Diversity ~ Day_Scaled + Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) #model failed to converge

Div_PxA <- glmer(Diversity ~ Position * Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 
Div_PxD <- glmer(Diversity ~ Position * Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) 
Div_DxA <- glmer(Diversity ~ Day_Scaled * Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) #Model failed to converge

Div_PAD <- glmer(Diversity ~ Position + Age_Scaled + Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) #model failed to converge 
Div_PxAD <- glmer(Diversity ~ Position * Age_Scaled + Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_PxDA <- glmer(Diversity ~ Position * Day_Scaled + Age_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) #model failed to converge 
Div_PAxD <- glmer(Diversity ~ Position + Age_Scaled * Day_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) #model failed to converge 



divlist <- list("null" = Div_null, "P" = Div_P, "A" = Div_A, 
                "D" = Div_D, "PD" = Div_PD, "PxA" = Div_PxA, 
                "PxD" = Div_PxD, "PAD" = Div_PAD, "PxAD" = Div_PxAD)

aictab(divlist)
#Day is top model but Null is within 2 AICc

##Step 2 - Environmental Variables----

head(TaxModel)

TaxModel$Field_Area_Scaled <- scale(TaxModel$Field_Area_m2)
TaxModel$NDVI1km_Scaled <- scale(TaxModel$NDVIsum_1km)

Div_Height <- glmer(Diversity ~ Height + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_GC <- glmer(Diversity ~ GC_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel)

Div_FieldArea <- glmer(Diversity ~ Field_Area_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_FieldNDVI <- glmer(Diversity ~ NDVImean_Field + (1 | Field), family = Gamma(link = "log"), data = TaxModel)

Div_Water <- glmer(Diversity ~ X1km_Prop_Water + (1 | Field), family = Gamma(link = "log"), data = TaxModel)
Div_NDVI1km <- glmer(Diversity ~ NDVI1km_Scaled + (1 | Field), family = Gamma(link = "log"), data = TaxModel) #model failed to converge

divlist2 <- list("null" = Div_null, "height" = Div_Height, 
                "GC" = Div_GC, "Field Area" = Div_FieldArea, 
                "Field NDVI" = Div_FieldNDVI, "Water" = Div_Water)
aictab(divlist2)

#top model is Null all other models are within 2 AICcs

#since none of the design variables made it to the step 2 do I do more variations of the environmental variables with interactions to see if that is better than the null??

#END----
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
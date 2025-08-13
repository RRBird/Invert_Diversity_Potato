
#This script contains all the figures from analysis in Main Script

#RICHNESS----

###All Species ----
summary(Rich_All) #age + field size
richpred_all.1 <- richpredresults[[1]]
head(richpred_all.1);dim(richpred_all.1)


dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

#Crop age
plot(x = ModelRich2$Age_Scaled,y = ModelRich2$All,xlab = "Crop Age (Days)",ylab = 'All Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=3,line=0,at = -2.7,'a)',cex=1)

polygon(x = c(richpred_all.1$Age_Scaled[richpred_all.1$Field_Area_m2 == Predictions_Area[10]],rev(richpred_all.1$Age_Scaled[richpred_all.1$Field_Area_m2 == Predictions_Area[10]])), y = c(richpred_all.1$lci[richpred_all.1$Field_Area_m2 == Predictions_Area[10]],rev(richpred_all.1$uci[richpred_all.1$Field_Area_m2 == Predictions_Area[10]])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_all.1$Age_Scaled[richpred_all.1$Field_Area_m2 == Predictions_Area[10]],y = richpred_all.1$fit[richpred_all.1$Field_Area_m2 == Predictions_Area[10]],lwd = 2,col = 'grey30')

axis(side=1, at=seq(from=min(richpred_all.1$Age_Scaled),to=max(richpred_all.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))


#Field Size

plot(x = ModelRich2$Field_Area_m2,y = ModelRich2$All,xlab = expression("Field Size (ha)"),ylab = 'All Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = 3000,'b)',cex=1)

polygon(x = c(richpred_all.1$Field_Area_m2[richpred_all.1$Age_Scaled == Predictions_Age[10]],rev(richpred_all.1$Field_Area_m2[richpred_all.1$Age_Scaled == Predictions_Age[10]])), y = c(richpred_all.1$lci[richpred_all.1$Age_Scaled == Predictions_Age[10]],rev(richpred_all.1$uci[richpred_all.1$Age_Scaled == Predictions_Age[10]])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_all.1$Field_Area_m2[richpred_all.1$Age_Scaled == Predictions_Age[10]],y = richpred_all.1$fit[richpred_all.1$Age_Scaled == Predictions_Age[10]],lwd = 2,col = 'grey30')

axis(1, at = axTicks(1), labels = axTicks(1)/10000)


###Predators----
summary(Rich_Pred)

richpred_preds.1 <- richpredresults[[2]]
head(richpred_preds.1);dim(richpred_preds.1)

###Herbivores----
summary(Rich_Herb) #position + day

richpred_herb.1 <- richpredresults[[3]]
head(richpred_herb.1);dim(richpred_herb.1)

raw_x <- ifelse(ModelRich2$Position == "Inner", 1, 
                ifelse(ModelRich2$Position == "Outer", 2, NA))

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)
#TO DO - Work out what I'm actually doinf colour wise this looks weird----
#position - not that to make sure the order is right for position had to do 2:1
plot(x = 2:1,y = richpred_herb.1$fit [richpred_herb.1$Day_Scaled == Predictions_Day[10]],xlab = " ",ylab = 'Herbivore Species Richness', type = 'p',pch = 16,cex =2.5,col = 'grey30', las = 1, ylim=c(0,2.5),xaxt = "n",xlim = c(0,3))
mtext(side=3,line=0,at = -0.7,'a)',cex=1)
axis(side=1,at=2:1,labels=c('Outer','Inner'))

arrows(x0=2:1, y0=richpred_herb.1$lci [richpred_herb.1$Day_Scaled == Predictions_Day[10]],x1=2:1, y1=richpred_herb.1$uci[richpred_herb.1$Day_Scaled == Predictions_Day[10]],angle=90,length=0.2, code=3, lwd=2,col = rgb(0.5, 0.5, 0.5, 0.5))

points(x = jitter(raw_x, factor = 1),y = ModelRich2$Herbivore, pch = 16, cex = 0.4, col = "black")

#Day

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Herbivore,xlab = expression("Day Sampled"),ylab = 'Herbivore Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.1,'b)',cex=1)
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1,-1.2,1.25,-1.2, length =0.1)

axis(side=1, at=seq(from=min(richpred_herb.1$Day_Scaled),to=max(richpred_herb.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1))

polygon(x = c(richpred_herb.1$Day_Scaled[richpred_herb.1$Position == Predictions_Position[1]],rev(richpred_herb.1$Day_Scaled[richpred_herb.1$Position == Predictions_Position[1]])), y = c(richpred_herb.1$lci[richpred_herb.1$Position == Predictions_Position[2]],rev(richpred_herb.1$uci[richpred_herb.1$Position == Predictions_Position[1]])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_herb.1$Day_Scaled[richpred_herb.1$Position == Predictions_Position[1]],y = richpred_herb.1$fit[richpred_herb.1$Position == Predictions_Position[1]],lwd = 2,col = 'grey30')


###Fungivores ----
summary(Rich_Fung) #age

richpred_fung.1 <- richpredresults[[4]]
head(richpred_fung.1);dim(richpred_fung.1)


dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Fungivore,xlab = "Crop Age (Days)",ylab = 'Fungivore Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')

polygon(x = c(richpred_fung.1$Age_Scaled,rev(richpred_fung.1$Age_Scaled)), y = c(richpred_fung.1$lci,rev(richpred_fung.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_fung.1$Age_Scaled,y = richpred_fung.1$fit,lwd = 2,col = 'grey30')

axis(side=1, at=seq(from=min(richpred_fung.1$Age_Scaled),to=max(richpred_fung.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

###Web Building----
summary(Rich_Web)

###Active Hunting----
summary(Rich_Active)

###Size 1 (0-2.5)----
summary(Rich_Size1)

###Size 2 (2.5-5)----
summary(Rich_Size2)

richpred_size2.1 <- richpredresults[[8]]
head(richpred_size2.1);dim(richpred_size2.1)


dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Size_2,xlab = "Crop Age (Days)",ylab = '2.5-5cm Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')

polygon(x = c(richpred_size2.1$Age_Scaled,rev(richpred_size2.1$Age_Scaled)), y = c(richpred_size2.1$lci,rev(richpred_size2.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_size2.1$Age_Scaled,y = richpred_size2.1$fit,lwd = 2,col = 'grey30')

axis(side=1, at=seq(from=min(richpred_size2.1$Age_Scaled),to=max(richpred_size2.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

###Size 3 (5-10)----
summary(Rich_Size3)

###Develops Wings----
summary(Rich_DevW) #age + NDVI 1km 

richpred_DevW.1 <- richpredresults[[10]]
head(richpred_DevW.1);dim(richpred_DevW.1)


dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

#Crop age
plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Develops_Wings,xlab = "Crop Age (Days)",ylab = 'Develops Wings Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=3,line=0,at = -2.7,'a)',cex=1)

polygon(x = c(richpred_DevW.1$Age_Scaled[richpred_DevW.1$NDVIsum_1km == Predictions_1kmNDVI[10]],rev(richpred_DevW.1$Age_Scaled[richpred_DevW.1$NDVIsum_1km == Predictions_1kmNDVI[10]])), y = c(richpred_DevW.1$lci[richpred_DevW.1$NDVIsum_1km == Predictions_1kmNDVI[10]],rev(richpred_DevW.1$uci[richpred_DevW.1$NDVIsum_1km == Predictions_1kmNDVI[10]])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_DevW.1$Age_Scaled[richpred_DevW.1$NDVIsum_1km == Predictions_1kmNDVI[10]],y = richpred_DevW.1$fit[richpred_DevW.1$NDVIsum_1km == Predictions_1kmNDVI[10]],lwd = 2,col = 'grey30')

axis(side=1, at=seq(from=min(richpred_DevW.1$Age_Scaled),to=max(richpred_DevW.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

#NDVI

plot(x = ModelRich2$NDVIsum_1km,y = ModelRich2$Develops_Wings,xlab = "NDVI (sum) within 1km",ylab = 'Develops Wings Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 3, lwd = 2)
mtext(side=3,line=0,at = 260,'b)',cex=1)

polygon(x = c(richpred_DevW.1$NDVIsum_1km[richpred_DevW.1$Age_Scaled == Predictions_Age[10]],rev(richpred_DevW.1$NDVIsum_1km[richpred_all.1$Age_Scaled == Predictions_Age[10]])), y = c(richpred_DevW.1$lci[richpred_DevW.1$Age_Scaled == Predictions_Age[10]],rev(richpred_DevW.1$uci[richpred_DevW.1$Age_Scaled == Predictions_Age[10]])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_DevW.1$NDVIsum_1km[richpred_DevW.1$Age_Scaled == Predictions_Age[10]],y = richpred_DevW.1$fit[richpred_DevW.1$Age_Scaled == Predictions_Age[10]],lwd = 2,col = 'grey30')


###Wingless----
summary(Rich_Wless)

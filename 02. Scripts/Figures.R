
#This script contains all the figures from analysis in Main Script

#RICHNESS----

###All Species ----
summary(Rich_All) #age + field size
richpred_all.1 <- richpredresults[[1]]
head(richpred_all.1);dim(richpred_all.1)


dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

#Age

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
summary(Rich_Pred) # Day + Age + Water

richpred_preds.1 <- richpredresults[[2]]
head(richpred_preds.1);dim(richpred_preds.1)

AA <- richpred_preds.1$Age_Scaled == Predictions_Age[10] & richpred_preds.1$X1km_Prop_Water == Predictions_Water[10]

AAA <- richpred_preds.1$Day_Scaled == Predictions_Day[10] & richpred_preds.1$X1km_Prop_Water == Predictions_Water[10]

AAAA <- richpred_preds.1$Age_Scaled == Predictions_Age[10] & richpred_preds.1$Day_Scaled == Predictions_Day[10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

#Day

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Predator,xlab = expression("Day Sampled"),ylab = 'Predator Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.1,'a)',cex=1)
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-0.7,-3.7,1.15,-3.7, length =0.1)

axis(side=1, at=seq(from=min(richpred_preds.1$Day_Scaled),to=max(richpred_preds.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_preds.1$Day_Scaled[AA],rev(richpred_preds.1$Day_Scaled[AA])), y = c(richpred_preds.1$lci[AA],rev(richpred_preds.1$uci[AA])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_preds.1$Day_Scaled[AA],y = richpred_preds.1$fit[AA],lwd = 2,col = 'grey30')

#Age

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Predator,xlab = expression("Crop Age (Days)"),ylab = 'Predator Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.8,'b)',cex=1)

axis(side=1, at=seq(from=min(richpred_preds.1$Age_Scaled),to=max(richpred_preds.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_preds.1$Age_Scaled[AAA],rev(richpred_preds.1$Age_Scaled[AAA])), y = c(richpred_preds.1$lci[AAA],rev(richpred_preds.1$uci[AAA])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_preds.1$Age_Scaled[AAA],y = richpred_preds.1$fit[AAA],lwd = 2,col = 'grey30')

#Water

plot(x = ModelRich2$X1km_Prop_Water,y = ModelRich2$Predator,xlab = expression("Proportion Water within 1km"),ylab = 'Predator Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = -0.5,'c)',cex=1)

polygon(x = c(richpred_preds.1$X1km_Prop_Water[AAAA],rev(richpred_preds.1$X1km_Prop_Water[AAAA])), y = c(richpred_preds.1$lci[AAAA],rev(richpred_preds.1$uci[AAAA])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_preds.1$X1km_Prop_Water[AAAA],y = richpred_preds.1$fit[AAAA],lwd = 2,col = 'grey30')

###Herbivores----
summary(Rich_Herb) #position + day

richpred_herb.1 <- richpredresults[[3]]
head(richpred_herb.1);dim(richpred_herb.1)

raw_x <- ifelse(ModelRich2$Position == "Inner", 1, 
                ifelse(ModelRich2$Position == "Outer", 2, NA))

B <- richpred_herb.1$Day_Scaled == Predictions_Day[10]
BB <- richpred_herb.1$Position == Predictions_Position[1]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

#TO DO - Am I keeping the raw data on there----
#position
#To make sure the order is right for position had to do 2:1

plot(x = 2:1,y = richpred_herb.1$fit [B],xlab = " ",ylab = 'Herbivore Species Richness', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,2.5),xaxt = "n",xlim = c(0,3))
mtext(side=3,line=0,at = -0.7,'a)',cex=1)
axis(side=1,at=2:1,labels=c('Outer','Inner'))

arrows(x0=2:1, y0=richpred_herb.1$lci [B],x1=2:1, y1=richpred_herb.1$uci[B],angle=90,length=0.2, code=3, lwd=2,col = "black")

points(x = jitter(raw_x, factor = 1),y = ModelRich2$Herbivore, pch = 16, cex = 0.4, col = "grey")

#Day

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Herbivore,xlab = expression("Day Sampled"),ylab = 'Herbivore Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.1,'b)',cex=1)
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1,-1.2,1.25,-1.2, length =0.1)

axis(side=1, at=seq(from=min(richpred_herb.1$Day_Scaled),to=max(richpred_herb.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1))

polygon(x = c(richpred_herb.1$Day_Scaled[BB],rev(richpred_herb.1$Day_Scaled[BB])), y = c(richpred_herb.1$lci[BB],rev(richpred_herb.1$uci[BB])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_herb.1$Day_Scaled[BB],y = richpred_herb.1$fit[BB],lwd = 2,col = 'grey30')


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

###Web Building TO DO----
summary(Rich_Web)

D

###Active Hunting----
summary(Rich_Active) # Day * Age

#NEED TO COME BACK TO THIS----
richpred_active.1 <- richpredresults[[6]]
head(richpred_active.1);dim(richpred_active.1)

#Trial and error to get these two - if you use min and max the fit predicted value is fine but the uci is at 200 at the max 
#Can I pick these as the two representative things for age and talk about the high uci in the results or discussion?
#Can I pick things that represent the relationship best or do I need to follow a protocol i.e. max and min, one standard deviation from mean etc?

E <- richpred_active.1$Age_Scaled == Predictions_Age[12]
EE <- richpred_active.1$Age_Scaled == Predictions_Age[20]

E <- richpred_active.1$Age_Scaled == min(Predictions_Age)
EE <- richpred_active.1$Age_Scaled == max(Predictions_Age)

seq(min(ModelRich2$Crop_Age_Days),max(ModelRich2$Crop_Age_Days),length.out=20)[12]
seq(min(ModelRich2$Crop_Age_Days),max(ModelRich2$Crop_Age_Days),length.out=20)[20]


dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Active_Hunting, xlab = "Day Sampled",ylab = 'Active Hunting Species Richness',type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',ylim = c(0,15))
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1,-4.5,1.3,-4.5, length =0.1)

axis(side=1, at=seq(from=min(richpred_active.1$Day_Scaled),to=max(richpred_active.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_active.1$Day_Scaled[E],rev(richpred_active.1$Day_Scaled[E])), y = c(richpred_active.1$lci[E],rev(richpred_active.1$uci[E])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_active.1$Day_Scaled[E],y = richpred_active.1$fit[E],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(richpred_active.1$Day_Scaled[EE],rev(richpred_active.1$Day_Scaled[EE])), y = c(richpred_active.1$lci[EE],rev(richpred_active.1$uci[EE])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_active.1$Day_Scaled[EE],y = richpred_active.1$fit[EE],lwd = 2,lty = 2, col = 'grey30')

legend('topleft',legend = c('Crop Age 84 Days', "Crop Age 109 Days"), lty = c(1,2), col = 'grey30',pt.cex = 1)





###Size 1 (0-2.5)----
summary(Rich_Size1) #Position + Age + Day

richpred_size1.1 <- richpredresults[[7]]
head(richpred_size1.1);dim(richpred_size1.1)

FF <- richpred_size1.1$Age_Scaled == Predictions_Age[10] & richpred_size1.1$Day_Scaled == Predictions_Day[10]
FFF <- richpred_size1.1$Day_Scaled == Predictions_Day[10] & richpred_size1.1$Position == Predictions_Position[1]
FFFF <- richpred_size1.1$Age_Scaled == Predictions_Age[10] & richpred_size1.1$Position == Predictions_Position[1]

raw_x2 <- ifelse(ModelRich2$Position == "Inner", 1, 
                ifelse(ModelRich2$Position == "Outer", 2, NA))

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)
#TO DO - Am I keeping the raw data on there?----

#position

plot(x = 1:2,y = richpred_size1.1$fit [FF],xlab = " ",ylab = 'Size 0-2.5cm Species Richness', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,3),xaxt = "n",xlim = c(0,3))
mtext(side=3,line=0,at = -0.7,'a)',cex=1)
axis(side=1,at=1:2,labels=c('Inner','Outer'))

arrows(x0=1:2, y0=richpred_size1.1$lci [FF],x1=1:2, y1=richpred_size1.1$uci[FF],angle=90,length=0.2, code=3, lwd=2,col = "black")

points(x = jitter(raw_x2, factor = 1),y = ModelRich2$Size_1, pch = 16, cex = 0.4, col = "grey")

#Age

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Size_1,xlab = expression("Crop Age (Days)"),ylab = 'Size 0-2.5cm Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.7,'b)',cex=1)

axis(side=1, at=seq(from=min(richpred_size1.1$Age_Scaled),to=max(richpred_size1.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_size1.1$Age_Scaled[FFF],rev(richpred_size1.1$Age_Scaled[FFF])), y = c(richpred_size1.1$lci[FFF],rev(richpred_size1.1$uci[FFF])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_size1.1$Age_Scaled[FFF],y = richpred_size1.1$fit[FFF],lwd = 2,col = 'grey30')

#Day 

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Size_1,xlab = expression("Day Sampled"),ylab = 'Size 0-2.5cm Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.1,'a)',cex=1)
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-0.7,-2.1,1.15,-2.1, length =0.1)

axis(side=1, at=seq(from=min(richpred_size1.1$Day_Scaled),to=max(richpred_size1.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_size1.1$Day_Scaled[FFFF],rev(richpred_size1.1$Day_Scaled[FFFF])), y = c(richpred_size1.1$lci[FFFF],rev(richpred_size1.1$uci[FFFF])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_size1.1$Day_Scaled[FFFF],y = richpred_size1.1$fit[FFFF],lwd = 2,col = 'grey30')

###Size 2 (2.5-5)----
summary(Rich_Size2) #Age

richpred_size2.1 <- richpredresults[[8]]
head(richpred_size2.1);dim(richpred_size2.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Size_2,xlab = "Crop Age (Days)",ylab = '2.5-5cm Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')

polygon(x = c(richpred_size2.1$Age_Scaled,rev(richpred_size2.1$Age_Scaled)), y = c(richpred_size2.1$lci,rev(richpred_size2.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_size2.1$Age_Scaled,y = richpred_size2.1$fit,lwd = 2,col = 'grey30')

axis(side=1, at=seq(from=min(richpred_size2.1$Age_Scaled),to=max(richpred_size2.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

###Size 3 (5-10) TO DO----
summary(Rich_Size3)

G
GG


###Develops Wings----
summary(Rich_DevW) #age + NDVI 1km 

richpred_DevW.1 <- richpredresults[[10]]
head(richpred_DevW.1);dim(richpred_DevW.1)

H <- richpred_DevW.1$NDVIsum_1km == Predictions_1kmNDVI[10]
HH <- richpred_all.1$Age_Scaled == Predictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

#Crop age 
plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Develops_Wings,xlab = "Crop Age (Days)",ylab = 'Develops Wings Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=3,line=0,at = -2.7,'a)',cex=1)

axis(side=1, at=seq(from=min(richpred_DevW.1$Age_Scaled),to=max(richpred_DevW.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(richpred_DevW.1$Age_Scaled[H],rev(richpred_DevW.1$Age_Scaled[H])), y = c(richpred_DevW.1$lci[H],rev(richpred_DevW.1$uci[H])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_DevW.1$Age_Scaled[H],y = richpred_DevW.1$fit[H],lwd = 2,col = 'grey30')


#NDVI

plot(x = ModelRich2$NDVIsum_1km,y = ModelRich2$Develops_Wings,xlab = "NDVI (sum) within 1km",ylab = 'Develops Wings Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 3, lwd = 2)
mtext(side=3,line=0,at = 260,'b)',cex=1)

polygon(x = c(richpred_DevW.1$NDVIsum_1km[HH],rev(richpred_DevW.1$NDVIsum_1km[HH])), y = c(richpred_DevW.1$lci[HH],rev(richpred_DevW.1$uci[HH])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_DevW.1$NDVIsum_1km[HH],y = richpred_DevW.1$fit[HH],lwd = 2,col = 'grey30')


###Wingless TO DO----
summary(Rich_Wless) #postion * age + day

richpred_Wless.1 <- richpredresults[[11]]
head(richpred_Wless.1);dim(richpred_Wless.1)

I
II

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

#Position * Age 

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Wingless,xlab = "Crop Age (Days)",ylab = 'Wingless Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=3,line=0,at = -2.7,'a)',cex=1)

axis(side=1, at=seq(from=min(richpred_Wless.1$Age_Scaled),to=max(richpred_Wless.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))







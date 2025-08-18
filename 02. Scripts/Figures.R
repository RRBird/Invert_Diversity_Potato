
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
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.1,-3.7,1.15,-3.7, length =0.1)

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
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.2,-1.2,1.25,-1.2, length =0.1)

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

###Web Building----
summary(Rich_Web) #Position * Day + Age + GC

richpred_web.1 <- richpredresults[[5]]
head(richpred_web.1);dim(richpred_web.1)

DD <- richpred_web.1$Position == "Inner" & richpred_web.1$Age_Scaled == Predictions_Age[10] & richpred_web.1$GC == Predictions_GC[10]
DDD <- richpred_web.1$Position == "Outer" & richpred_web.1$Age_Scaled == Predictions_Age[10] & richpred_web.1$GC == Predictions_GC[10]
DDDD <- richpred_web.1$Position == "Outer" & richpred_web.1$Day_Scaled == Predictions_Day[10] & richpred_web.1$GC == Predictions_GC[10]
DDDDD <- richpred_web.1$Position == "Outer" & richpred_web.1$Day_Scaled == Predictions_Day[10] & richpred_web.1$Age_Scaled == Predictions_Age[10]


dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

#Position * Day 

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Web, xlab = "Day Sampled",ylab = 'Web Building Species Richness',type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.1,-2.1,1.1,-2.1, length =0.1)
mtext(side=3,line=0,at = -2,'a)',cex=1)

axis(side=1, at=seq(from=min(richpred_web.1$Day_Scaled),to=max(richpred_web.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_web.1$Day_Scaled[DD],rev(richpred_web.1$Day_Scaled[DD])), y = c(richpred_web.1$lci[DD],rev(richpred_web.1$uci[DD])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_web.1$Day_Scaled[DD],y = richpred_web.1$fit[DD],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(richpred_web.1$Day_Scaled[DDD],rev(richpred_web.1$Day_Scaled[DDD])), y = c(richpred_web.1$lci[DDD],rev(richpred_web.1$uci[DDD])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_web.1$Day_Scaled[DDD],y = richpred_web.1$fit[DDD],lwd = 2,lty = 2, col = 'grey30')

legend('topleft',legend = c('Inner', "Outer"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Age

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Web,xlab = "Crop Age (Days)",ylab = 'Web Buliding Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(richpred_web.1$Age_Scaled),to=max(richpred_web.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.6,'b)',cex=1)

polygon(x = c(richpred_web.1$Age_Scaled[DDDD],rev(richpred_web.1$Age_Scaled[DDDD])), y = c(richpred_web.1$lci[DDDD],rev(richpred_web.1$uci[DDDD])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_web.1$Age_Scaled[DDDD],y = richpred_web.1$fit[DDDD],lwd = 2,col = 'grey30')

#GC

plot(x = ModelRich2$GC,y = ModelRich2$Web,xlab = "Ground Cover (%)",ylab = 'Web Buliding Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = -6,'c)',cex=1)

polygon(x = c(richpred_web.1$GC[DDDDD],rev(richpred_web.1$GC[DDDDD])), y = c(richpred_web.1$lci[DDDDD],rev(richpred_web.1$uci[DDDDD])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_web.1$GC[DDDDD],y = richpred_web.1$fit[DDDDD],lwd = 2,col = 'grey30')


###Active Hunting----
summary(Rich_Active) # Day * Age

richpred_active.1 <- richpredresults[[6]]
head(richpred_active.1);dim(richpred_active.1)

#working out days to represent winter and spring
Predictions_Day
seq(min(ModelRich2$Day_Sampled),max(ModelRich2$Day_Sampled),length.out=20) 
Predictions_Day[4] #winter
Predictions_Day[15] #Spring

E <- richpred_active.1$Day_Scaled == Predictions_Day[4]
EE <- richpred_active.1$Day_Scaled == Predictions_Day[15]


dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Active_Hunting, xlab = "Crop Age (Days)",ylab = 'Active Hunting Species Richness',type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',ylim=c(0,7))
axis(side=1, at=seq(from=min(richpred_active.1$Age_Scaled),to=max(richpred_active.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_active.1$Age_Scaled[E],rev(richpred_active.1$Age_Scaled[E])), y = c(richpred_active.1$lci[E],rev(richpred_active.1$uci[E])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_active.1$Age_Scaled[E],y = richpred_active.1$fit[E],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(richpred_active.1$Age_Scaled[EE],rev(richpred_active.1$Age_Scaled[EE])), y = c(richpred_active.1$lci[EE],rev(richpred_active.1$uci[EE])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_active.1$Age_Scaled[EE],y = richpred_active.1$fit[EE],lwd = 2,lty = 2, col = 'grey30')

legend('topright',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)


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

#position

plot(x = 1:2,y = richpred_size1.1$fit [FF],xlab = " ",ylab = 'Size 0-2.5cm Species Richness', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,3),xaxt = "n",xlim = c(0,3))
mtext(side=3,line=0,at = -0.7,'a)',cex=1)
axis(side=1,at=1:2,labels=c('Inner','Outer'))

arrows(x0=1:2, y0=richpred_size1.1$lci [FF],x1=1:2, y1=richpred_size1.1$uci[FF],angle=90,length=0.1, code=3, lwd=2,col = "black")

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
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.1,-2.1,1.15,-2.1, length =0.1)

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

###Size 3 (5-10)----
summary(Rich_Size3) #Position + Age * Day + Water

richpred_size3.1 <- richpredresults[[9]]
head(richpred_size3.1);dim(richpred_size3.1)

#working out days to represent winter and spring
Predictions_Day
seq(min(ModelRich2$Day_Sampled),max(ModelRich2$Day_Sampled),length.out=20) 
Predictions_Day[4] #winter
Predictions_Day[15] #Spring

G <- richpred_size3.1$Age_Scaled == Predictions_Age[10] & richpred_size3.1$Day_Scaled == Predictions_Day[10] & richpred_size3.1$X1km_Prop_Water == Predictions_Water[10]
GG <- richpred_size3.1$Day_Scaled == Predictions_Day[4] & richpred_size3.1$X1km_Prop_Water == Predictions_Water[10] & richpred_size3.1$Position == "Outer"
GGG <- richpred_size3.1$Day_Scaled == Predictions_Day[15] & richpred_size3.1$X1km_Prop_Water == Predictions_Water[10] & richpred_size3.1$Position == "Outer"
GGGG <- richpred_size3.1$Day_Scaled == Predictions_Day[10] & richpred_size3.1$Position == "Outer" & richpred_size3.1$Age_Scaled == Predictions_Age[10] 

raw_x3 <- ifelse(ModelRich2$Position == "Inner", 1, 
                 ifelse(ModelRich2$Position == "Outer", 2, NA))


dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

#position

plot(x = 2:1,y = richpred_size3.1$fit [G],xlab = " ",ylab = 'Size 5-10cm Species Richness', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1,xaxt = "n",ylim = c(0,1), xlim = c(0,3))
mtext(side=3,line=0,at = -0.7,'a)',cex=1)
axis(side=1,at=2:1,labels=c('Outer','Inner'))

arrows(x0=2:1, y0=richpred_size3.1$lci[G],x1=2:1, y1=richpred_size3.1$uci[G],angle=90,length=0.2, code=3, lwd=2,col = "black")

points(x = jitter(raw_x2, factor = 1),y = ModelRich2$Size_3, pch = 16, cex = 0.4, col = "grey")

#Age * Day

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Size_3,xlab = "Crop Age (Days)",ylab = '5-10cm Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',ylim = c(0,7))
mtext(side=3,line=0,at = -2.6,'b)',cex=1)
axis(side=1, at=seq(from=min(richpred_size3.1$Age_Scaled),to=max(richpred_size3.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_size3.1$Age_Scaled[GG],rev(richpred_size3.1$Age_Scaled[GG])), y = c(richpred_size3.1$lci[GG],rev(richpred_size3.1$uci[GG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_size3.1$Age_Scaled[GG],y = richpred_size3.1$fit[GG],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(richpred_size3.1$Age_Scaled[GGG],rev(richpred_size3.1$Age_Scaled[GGG])), y = c(richpred_size3.1$lci[GGG],rev(richpred_size3.1$uci[GGG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_size3.1$Age_Scaled[GGG],y = richpred_size3.1$fit[GGG],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Water

plot(x = ModelRich2$X1km_Prop_Water,y = ModelRich2$Size_3,xlab = expression("Proportion Water within 1km"),ylab = '5-10cm Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = -0.5,'c)',cex=1)

polygon(x = c(richpred_size3.1$X1km_Prop_Water[GGGG],rev(richpred_size3.1$X1km_Prop_Water[GGGG])), y = c(richpred_size3.1$lci[GGGG],rev(richpred_size3.1$uci[GGGG])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=richpred_size3.1$X1km_Prop_Water[GGGG],y = richpred_size3.1$fit[GGGG],lwd = 2,col = 'grey30')


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


###Wingless----
summary(Rich_Wless) #postion * age + day

richpred_Wless.1 <- richpredresults[[11]]
head(richpred_Wless.1);dim(richpred_Wless.1)

I <- richpred_Wless.1$Position == "Inner" & richpred_Wless.1$Day_Scaled == Predictions_Day[10]
II <- richpred_Wless.1$Position == "Outer" & richpred_Wless.1$Day_Scaled == Predictions_Day[10]

III <- richpred_Wless.1$Position == "Outer" & richpred_Wless.1$Age_Scaled == Predictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

#Position * Age 

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Wingless,xlab = "Crop Age (Days)",ylab = 'Wingless Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=3,line=0,at = -2.7,'a)',cex=1)

axis(side=1, at=seq(from=min(richpred_Wless.1$Age_Scaled),to=max(richpred_Wless.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(richpred_Wless.1$Age_Scaled[I],rev(richpred_Wless.1$Age_Scaled[I])), y = c(richpred_Wless.1$lci[I],rev(richpred_Wless.1$uci[I])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_Wless.1$Age_Scaled[I],y = richpred_Wless.1$fit[I],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(richpred_Wless.1$Age_Scaled[II],rev(richpred_Wless.1$Age_Scaled[II])), y = c(richpred_Wless.1$lci[II],rev(richpred_Wless.1$uci[II])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_Wless.1$Age_Scaled[II],y = richpred_Wless.1$fit[II],lwd = 2,lty = 2, col = 'grey30')

legend('topleft',legend = c('Inner', "Outer"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Day

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Wingless, xlab = "Day Sampled",ylab = 'Wingless Species Richness',type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.2,-1.8,1.2,-1.8, length =0.1)

axis(side=1, at=seq(from=min(richpred_Wless.1$Day_Scaled),to=max(richpred_Wless.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_Wless.1$Day_Scaled[III],rev(richpred_Wless.1$Day_Scaled[III])), y = c(richpred_Wless.1$lci[III],rev(richpred_Wless.1$uci[III])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_Wless.1$Day_Scaled[III],y = richpred_Wless.1$fit[III],lwd = 2,lty = 1, col = 'grey30')

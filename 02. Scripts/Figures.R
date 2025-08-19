
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


#OCCURRENCE----

##Herbivore----

summary(Occur_herb) #Day
occurpred_herb.1 <- occurpredresults[[1]]
head(occurpred_herb.1);dim(occurpred_herb.1)


dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$Day_Scaled,y = ModelOccur2$Herbivore,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Herbivore")
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.2,-0.3,1.25,-0.3, length =0.1)
axis(side=1, at=seq(from=min(occurpred_herb.1$Day_Scaled),to=max(occurpred_herb.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Day_Sampled),to=max(ModelOccur2$Day_Sampled),length.out=6),-1))

polygon(x = c(occurpred_herb.1$Day_Scaled,rev(occurpred_herb.1$Day_Scaled)), y = c(occurpred_herb.1$lci,rev(occurpred_herb.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_herb.1$Day_Scaled,y = occurpred_herb.1$fit,lwd = 2,col = 'grey30')

##Omnivore----

summary(Occur_omni) #Ground Cover
occurpred_omni.1 <- occurpredresults[[2]]
head(occurpred_omni.1);dim(occurpred_omni.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$GC,y = ModelOccur2$Omnivore,xlab = "Ground Cover (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,main = "Omnivore")

polygon(x = c(occurpred_omni.1$GC,rev(occurpred_omni.1$GC)), y = c(occurpred_omni.1$lci,rev(occurpred_omni.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_omni.1$GC,y = occurpred_omni.1$fit,lwd = 2,col = 'grey30')

##Fungivore----

summary(Occur_fung) #Day + Age + Field Size
occurpred_fung.1 <- occurpredresults[[3]]
head(occurpred_fung.1);dim(occurpred_fung.1)

J <- occurpred_fung.1$Age_Scaled == Predictions_Age[10] & occurpred_fung.1$FieldSize_Scaled == Predictions_Areascale[10]
JJ <- occurpred_fung.1$Day_Scaled == Predictions_Day[10] & occurpred_fung.1$FieldSize_Scaled == Predictions_Areascale[10]
JJJ <- occurpred_fung.1$Age_Scaled == Predictions_Age[10] & occurpred_fung.1$Day_Scaled == Predictions_Day[10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow = c(2,2))

#Day

plot(x = ModelOccur2$Day_Scaled,y = ModelOccur2$Fungivore,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.1,-0.42,1.1,-0.42, length =0.1)
axis(side=1, at=seq(from=min(occurpred_fung.1$Day_Scaled),to=max(occurpred_fung.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Day_Sampled),to=max(ModelOccur2$Day_Sampled),length.out=6),-1),cex.axis=0.9)
mtext(side=3,line=0.1,at = 2,expression(bold('Fungivore')),cex=1.2)

polygon(x = c(occurpred_fung.1$Day_Scaled[J],rev(occurpred_fung.1$Day_Scaled[J])), y = c(occurpred_fung.1$lci[J],rev(occurpred_fung.1$uci[J])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_fung.1$Day_Scaled[J],y = occurpred_fung.1$fit[J],lwd = 2,col = 'grey30')

#Age

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Fungivore,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
axis(side=1, at=seq(from=min(occurpred_fung.1$Age_Scaled),to=max(occurpred_fung.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(occurpred_fung.1$Age_Scaled[JJ],rev(occurpred_fung.1$Age_Scaled[JJ])), y = c(occurpred_fung.1$lci[JJ],rev(occurpred_fung.1$uci[JJ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_fung.1$Age_Scaled[JJ],y = occurpred_fung.1$fit[JJ],lwd = 2,col = 'grey30')

#Field Size

plot(x = ModelOccur2$FieldSize_Scaled,y = ModelOccur2$Fungivore,xlab = "Field Size (km2)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
axis(side=1, at=seq(from=min(occurpred_fung.1$FieldSize_Scaled),to=max(occurpred_fung.1$FieldSize_Scaled),length.out=5),labels=round(seq(from=min(ModelOccur2$Field_Area_m2),to=max(ModelOccur2$Field_Area_m2),length.out=5)/1000000,2),cex.axis = 0.9)

polygon(x = c(occurpred_fung.1$FieldSize_Scaled[JJJ],rev(occurpred_fung.1$FieldSize_Scaled[JJJ])), y = c(occurpred_fung.1$lci[JJJ],rev(occurpred_fung.1$uci[JJJ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_fung.1$FieldSize_Scaled[JJJ],y = occurpred_fung.1$fit[JJJ],lwd = 2,col = 'grey30')

##Hematophagous----

summary(Occur_hema) #Day * Age + NDVI Field
occurpred_hema.1 <- occurpredresults[[4]]
head(occurpred_hema.1);dim(occurpred_hema.1)

K <- occurpred_hema.1$Day_Scaled == Predictions_Day[4] & occurpred_hema.1$NDVImean_Field == Predictions_NDVIfield[10]
KK <- occurpred_hema.1$Day_Scaled == Predictions_Day[15] & occurpred_hema.1$NDVImean_Field == Predictions_NDVIfield[10]
K_K <- occurpred_hema.1$Day_Scaled == Predictions_Day[10] & occurpred_hema.1$Age_Scaled == Predictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow = c(1,2))

#Day * Age

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Hematophagous,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
axis(side=1, at=seq(from=min(occurpred_hema.1$Age_Scaled),to=max(occurpred_hema.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0.1,at = 2,expression(bold('Hematophagous')),cex=1.2)

polygon(x = c(occurpred_hema.1$Age_Scaled[K],rev(occurpred_hema.1$Age_Scaled[K])), y = c(occurpred_hema.1$lci[K],rev(occurpred_hema.1$uci[K])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hema.1$Age_Scaled[K],y = occurpred_hema.1$fit[K],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(occurpred_hema.1$Age_Scaled[KK],rev(occurpred_hema.1$Age_Scaled[KK])), y = c(occurpred_hema.1$lci[KK],rev(occurpred_hema.1$uci[KK])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hema.1$Age_Scaled[KK],y = occurpred_hema.1$fit[KK],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#NDVI Field

plot(x = ModelOccur2$NDVImean_Field,y = ModelOccur2$Hematophagous,xlab = "Mean Field NDVI",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)

polygon(x = c(occurpred_hema.1$NDVImean_Field[K_K],rev(occurpred_hema.1$NDVImean_Field[K_K])), y = c(occurpred_hema.1$lci[K_K],rev(occurpred_hema.1$uci[K_K])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hema.1$NDVImean_Field[K_K],y = occurpred_hema.1$fit[K_K],lwd = 2,col = 'grey30',lty = 1)



##Web----

summary(Occur_web) #Day * Age
occurpred_web.1 <- occurpredresults[[5]]
head(occurpred_web.1);dim(occurpred_web.1)

Predictions_Day[4] #winter
Predictions_Day[15] #Spring

L <- occurpred_web.1$Day_Scaled == Predictions_Day[4]
LL <- occurpred_web.1$Day_Scaled == Predictions_Day[15]

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Web,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Web Building")
axis(side=1, at=seq(from=min(occurpred_web.1$Age_Scaled),to=max(occurpred_web.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(occurpred_web.1$Age_Scaled[L],rev(occurpred_web.1$Age_Scaled[L])), y = c(occurpred_web.1$lci[L],rev(occurpred_web.1$uci[L])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_web.1$Age_Scaled[L],y = occurpred_web.1$fit[L],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(occurpred_web.1$Age_Scaled[LL],rev(occurpred_web.1$Age_Scaled[LL])), y = c(occurpred_web.1$lci[LL],rev(occurpred_web.1$uci[LL])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_web.1$Age_Scaled[LL],y = occurpred_web.1$fit[LL],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)


##Active Hunting----

summary(Occur_active) #Day 
occurpred_active.1 <- occurpredresults[[6]]
head(occurpred_active.1);dim(occurpred_active.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$Day_Scaled,y = ModelOccur2$Active_Hunting,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Active Hunting")
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.2,-0.3,1.25,-0.3, length =0.1)
axis(side=1, at=seq(from=min(occurpred_active.1$Day_Scaled),to=max(occurpred_active.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Day_Sampled),to=max(ModelOccur2$Day_Sampled),length.out=6),-1))

polygon(x = c(occurpred_active.1$Day_Scaled,rev(occurpred_active.1$Day_Scaled)), y = c(occurpred_active.1$lci,rev(occurpred_active.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_active.1$Day_Scaled,y = occurpred_active.1$fit,lwd = 2,col = 'grey30')

##Hawking----

summary(Occur_hawk) #Position + Age * Day Scaled 
occurpred_hawk.1 <- occurpredresults[[7]]
head(occurpred_hawk.1);dim(occurpred_hawk.1)

Predictions_Day[4] #Winter
Predictions_Day[15] #Spring

M <- occurpred_hawk.1$Age_Scaled == Predictions_Age[10] & occurpred_hawk.1$Day_Scaled == Predictions_Day [10]
MM <- occurpred_hawk.1$Position == "Outer" & occurpred_hawk.1$Day_Scaled == Predictions_Day[4]
MMM <- occurpred_hawk.1$Position == "Outer" & occurpred_hawk.1$Day_Scaled == Predictions_Day[4]

raw_x4 <- ifelse(ModelOccur2$Position == "Inner", 1, 
                ifelse(ModelOccur2$Position == "Outer", 2, NA))

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow = c(1,2))

#Not fully updated yet
#the y= in plot isn't bringing up only two why?

plot(x = 1:2,y = occurpred_hawk.1$fit[M],xlab = " ",ylab = 'Probability of Occurrence', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,2.5),xaxt = "n",xlim = c(0,3))
mtext(side=3,line=0,at = -0.7,'a)',cex=1)
axis(side=1,at=1:2,labels=c('Outer','Inner'))

arrows(x0=1:2, y0=occurpred_hawk.1$lci [M],x1=1:2, y1=richpred_herb.1$uci[M],angle=90,length=0.2, code=3, lwd=2,col = "black")

points(x = jitter(raw_x4, factor = 1),y = ModelOccur2$Hawking, pch = 16, cex = 0.4, col = "grey")



##Size 2 (2.5-5cm)----

summary(Occur_size2) #Height 
occurpred_size2.1 <- occurpredresults[[8]]
head(occurpred_size2.1);dim(occurpred_size2.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$Height,y = ModelOccur2$Size_2,xlab = "Plant Height (cm)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,main = "Size 2 (2.5-5cm)")

polygon(x = c(occurpred_size2.1$Height,rev(occurpred_size2.1$Height)), y = c(occurpred_size2.1$lci,rev(occurpred_size2.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_size2.1$Height,y = occurpred_size2.1$fit,lwd = 2,col = 'grey30')

##Size 3 (5-10cm)----

summary(Occur_size3) #Day * Age 
occurpred_size3.1 <- occurpredresults[[9]]
head(occurpred_size3.1);dim(occurpred_size3.1)

Predictions_Day[4] #winter
Predictions_Day[15] #Spring

N <- occurpred_size3.1$Day_Scaled == Predictions_Day[4]
NN <- occurpred_size3.1$Day_Scaled == Predictions_Day[15]

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Size_3,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Size 3 (5-10cm)")
axis(side=1, at=seq(from=min(occurpred_size3.1$Age_Scaled),to=max(occurpred_size3.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(occurpred_size3.1$Age_Scaled[N],rev(occurpred_size3.1$Age_Scaled[N])), y = c(occurpred_size3.1$lci[N],rev(occurpred_size3.1$uci[N])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_size3.1$Age_Scaled[N],y = occurpred_size3.1$fit[N],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(occurpred_size3.1$Age_Scaled[NN],rev(occurpred_size3.1$Age_Scaled[NN])), y = c(occurpred_size3.1$lci[NN],rev(occurpred_size3.1$uci[NN])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_size3.1$Age_Scaled[NN],y = occurpred_size3.1$fit[NN],lwd = 2,col = 'grey30',lty = 2)

legend('topright',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

##Size 4 (>10cm)----

summary(Occur_size4) #Age 
occurpred_size4.1 <- occurpredresults[[10]]
head(occurpred_size4.1);dim(occurpred_size4.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Size_4,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Size 4 (>10cm)")
axis(side=1, at=seq(from=min(occurpred_size4.1$Age_Scaled),to=max(occurpred_size4.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(occurpred_size4.1$Age_Scaled,rev(occurpred_size4.1$Age_Scaled)), y = c(occurpred_size4.1$lci,rev(occurpred_size4.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_size4.1$Age_Scaled,y = occurpred_size4.1$fit,lwd = 2,col = 'grey30')


##Always Winged----

summary(Occur_AlWing) #Water 
occurpred_AlWing.1 <- occurpredresults[[11]]
head(occurpred_AlWing.1);dim(occurpred_AlWing.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$X1km_Prop_Water,y = ModelOccur2$Always_Winged,xlab = "Water within 1km (%)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,main = "Always Winged")

polygon(x = c(occurpred_AlWing.1$X1km_Prop_Water,rev(occurpred_AlWing.1$X1km_Prop_Water)), y = c(occurpred_AlWing.1$lci,rev(occurpred_AlWing.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_AlWing.1$X1km_Prop_Water,y = occurpred_AlWing.1$fit,lwd = 2,col = 'grey30')

##Wingless----

summary(Occur_Wless) #Position * Age + Day 
occurpred_Wless.1 <- occurpredresults[[12]]
head(occurpred_Wless.1);dim(occurpred_Wless.1)

O
OO
OOO

##Polymorphic----

summary(Occur_poly) #Position * Day 
occurpred_poly.1 <- occurpredresults[[13]]
head(occurpred_poly.1);dim(occurpred_poly.1)

PP <- occurpred_poly.1$Position == "Inner"
PPP <- occurpred_poly.1$Position == "Outer"

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$Day_Scaled,y = ModelOccur2$Polymorphic,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Polymorphic (Wings)")
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.2,-0.3,1.25,-0.3, length =0.1)
axis(side=1, at=seq(from=min(occurpred_poly.1$Day_Scaled),to=max(occurpred_poly.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Day_Sampled),to=max(ModelOccur2$Day_Sampled),length.out=6),-1))

polygon(x = c(occurpred_poly.1$Day_Scaled[PP],rev(occurpred_poly.1$Day_Scaled[PP])), y = c(occurpred_poly.1$lci[PP],rev(occurpred_poly.1$uci[PP])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_poly.1$Day_Scaled[PP],y = occurpred_poly.1$fit[PP],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(occurpred_poly.1$Day_Scaled[PPP],rev(occurpred_poly.1$Day_Scaled[PPP])), y = c(occurpred_poly.1$lci[PPP],rev(occurpred_poly.1$uci[PPP])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_poly.1$Day_Scaled[PPP],y = occurpred_poly.1$fit[PPP],lwd = 2,col = 'grey30',lty = 2)


legend('topleft',legend = c('Inner', "Outer"), lty = c(1,2), col = 'grey30',pt.cex = 1)


#END----
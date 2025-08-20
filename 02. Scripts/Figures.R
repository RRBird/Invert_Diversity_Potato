
#This script contains all the figures from analysis in Main Script

#RICHNESS----

###All Species ----
summary(Rich_All) #age + field size
richpred_all.1 <- richpredresults[[1]]
head(richpred_all.1);dim(richpred_all.1)


dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

#Age

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$All,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=3,line=0,at = -2.3,'a)',cex=1)
mtext(side=3,line=0.1,at = 2,expression(bold('All')),cex=1.2)

polygon(x = c(richpred_all.1$Age_Scaled[richpred_all.1$Field_Area_m2 == Predictions_Area[10]],rev(richpred_all.1$Age_Scaled[richpred_all.1$Field_Area_m2 == Predictions_Area[10]])), y = c(richpred_all.1$lci[richpred_all.1$Field_Area_m2 == Predictions_Area[10]],rev(richpred_all.1$uci[richpred_all.1$Field_Area_m2 == Predictions_Area[10]])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_all.1$Age_Scaled[richpred_all.1$Field_Area_m2 == Predictions_Area[10]],y = richpred_all.1$fit[richpred_all.1$Field_Area_m2 == Predictions_Area[10]],lwd = 2,col = 'grey30')

axis(side=1, at=seq(from=min(richpred_all.1$Age_Scaled),to=max(richpred_all.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

#Field Size

plot(x = ModelRich2$Field_Area_m2,y = ModelRich2$All,xlab = expression("Field Size (ha)"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = 10000,'b)',cex=1)

polygon(x = c(richpred_all.1$Field_Area_m2[richpred_all.1$Age_Scaled == Predictions_Age[10]],rev(richpred_all.1$Field_Area_m2[richpred_all.1$Age_Scaled == Predictions_Age[10]])), y = c(richpred_all.1$lci[richpred_all.1$Age_Scaled == Predictions_Age[10]],rev(richpred_all.1$uci[richpred_all.1$Age_Scaled == Predictions_Age[10]])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_all.1$Field_Area_m2[richpred_all.1$Age_Scaled == Predictions_Age[10]],y = richpred_all.1$fit[richpred_all.1$Age_Scaled == Predictions_Age[10]],lwd = 2,col = 'grey30')

axis(1, at = axTicks(1), labels = axTicks(1)/10000)


###Predator----
summary(Rich_Pred) # Day + Age + Water

richpred_preds.1 <- richpredresults[[2]]
head(richpred_preds.1);dim(richpred_preds.1)

AA <- richpred_preds.1$Age_Scaled == Predictions_Age[10] & richpred_preds.1$X1km_Prop_Water == Predictions_Water[10]

AAA <- richpred_preds.1$Day_Scaled == Predictions_Day[10] & richpred_preds.1$X1km_Prop_Water == Predictions_Water[10]

AAAA <- richpred_preds.1$Age_Scaled == Predictions_Age[10] & richpred_preds.1$Day_Scaled == Predictions_Day[10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

#Day

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Predator,xlab = expression("Day Sampled"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -1.7,'a)',cex=1)
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.1,-3.7,1.15,-3.7, length =0.1)
mtext(side=3,line=0.1,at = 2,expression(bold('Predator')),cex=1.2)

axis(side=1, at=seq(from=min(richpred_preds.1$Day_Scaled),to=max(richpred_preds.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_preds.1$Day_Scaled[AA],rev(richpred_preds.1$Day_Scaled[AA])), y = c(richpred_preds.1$lci[AA],rev(richpred_preds.1$uci[AA])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_preds.1$Day_Scaled[AA],y = richpred_preds.1$fit[AA],lwd = 2,col = 'grey30')

#Age

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Predator,xlab = expression("Crop Age (Days)"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.3,'b)',cex=1)

axis(side=1, at=seq(from=min(richpred_preds.1$Age_Scaled),to=max(richpred_preds.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_preds.1$Age_Scaled[AAA],rev(richpred_preds.1$Age_Scaled[AAA])), y = c(richpred_preds.1$lci[AAA],rev(richpred_preds.1$uci[AAA])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_preds.1$Age_Scaled[AAA],y = richpred_preds.1$fit[AAA],lwd = 2,col = 'grey30')

#Water

plot(x = ModelRich2$X1km_Prop_Water,y = ModelRich2$Predator,xlab = "Water within 1km (%)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 0.8,'c)',cex=1)

polygon(x = c(richpred_preds.1$X1km_Prop_Water[AAAA],rev(richpred_preds.1$X1km_Prop_Water[AAAA])), y = c(richpred_preds.1$lci[AAAA],rev(richpred_preds.1$uci[AAAA])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_preds.1$X1km_Prop_Water[AAAA],y = richpred_preds.1$fit[AAAA],lwd = 2,col = 'grey30')

###Herbivore----
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

plot(x = 2:1,y = richpred_herb.1$fit [B],xlab = " ",ylab = 'Species Richness', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,2.5),xaxt = "n",xlim = c(0,3))
mtext(side=3,line=0,at = -0.25,'a)',cex=1)
axis(side=1,at=2:1,labels=c('Outer','Inner'))
mtext(side=3,line=0.1,at = 3.5,expression(bold('Herbivore')),cex=1.2)

arrows(x0=2:1, y0=richpred_herb.1$lci [B],x1=2:1, y1=richpred_herb.1$uci[B],angle=90,length=0.2, code=3, lwd=2,col = "black")

points(x = jitter(raw_x, factor = 1),y = ModelRich2$Herbivore, pch = 16, cex = 0.4, col = "grey")

#Day

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Herbivore,xlab = expression("Day Sampled"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -1.6,'b)',cex=1)
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.2,-1.2,1.25,-1.2, length =0.1)

axis(side=1, at=seq(from=min(richpred_herb.1$Day_Scaled),to=max(richpred_herb.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1))

polygon(x = c(richpred_herb.1$Day_Scaled[BB],rev(richpred_herb.1$Day_Scaled[BB])), y = c(richpred_herb.1$lci[BB],rev(richpred_herb.1$uci[BB])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_herb.1$Day_Scaled[BB],y = richpred_herb.1$fit[BB],lwd = 2,col = 'grey30')


###Fungivore ----
summary(Rich_Fung) #age

richpred_fung.1 <- richpredresults[[4]]
head(richpred_fung.1);dim(richpred_fung.1)


dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Fungivore,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',main = "Fungivore")
axis(side=1, at=seq(from=min(richpred_fung.1$Age_Scaled),to=max(richpred_fung.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(richpred_fung.1$Age_Scaled,rev(richpred_fung.1$Age_Scaled)), y = c(richpred_fung.1$lci,rev(richpred_fung.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_fung.1$Age_Scaled,y = richpred_fung.1$fit,lwd = 2,col = 'grey30')


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

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Web, xlab = "Day Sampled",ylab = 'Species Richness',type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.1,-2.1,1.1,-2.1, length =0.1)
mtext(side=3,line=0,at = -1.65,'a)',cex=1)
mtext(side=3,line=0.1,at = 1.8,expression(bold('Web Building')),cex=1.2)

axis(side=1, at=seq(from=min(richpred_web.1$Day_Scaled),to=max(richpred_web.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_web.1$Day_Scaled[DD],rev(richpred_web.1$Day_Scaled[DD])), y = c(richpred_web.1$lci[DD],rev(richpred_web.1$uci[DD])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_web.1$Day_Scaled[DD],y = richpred_web.1$fit[DD],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(richpred_web.1$Day_Scaled[DDD],rev(richpred_web.1$Day_Scaled[DDD])), y = c(richpred_web.1$lci[DDD],rev(richpred_web.1$uci[DDD])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_web.1$Day_Scaled[DDD],y = richpred_web.1$fit[DDD],lwd = 2,lty = 2, col = 'grey30')

legend('topleft',legend = c('Inner', "Outer"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Age

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Web,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
axis(side=1, at=seq(from=min(richpred_web.1$Age_Scaled),to=max(richpred_web.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.2,'b)',cex=1)

polygon(x = c(richpred_web.1$Age_Scaled[DDDD],rev(richpred_web.1$Age_Scaled[DDDD])), y = c(richpred_web.1$lci[DDDD],rev(richpred_web.1$uci[DDDD])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_web.1$Age_Scaled[DDDD],y = richpred_web.1$fit[DDDD],lwd = 2,col = 'grey30')

#GC

plot(x = ModelRich2$GC,y = ModelRich2$Web,xlab = "Ground Cover (%)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 3,'c)',cex=1)

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

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Active_Hunting, xlab = "Crop Age (Days)",ylab = 'Species Richness',type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',ylim=c(0,7),main = "Active Hunting")
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

plot(x = 1:2,y = richpred_size1.1$fit [FF],xlab = " ",ylab = 'Species Richness', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,3),xaxt = "n",xlim = c(0,3))
mtext(side=3,line=0,at = -0.3,'a)',cex=1)
axis(side=1,at=1:2,labels=c('Inner','Outer'))
mtext(side=3,line=0.1,at = 3.5,expression(bold('0-2.5cm')),cex=1.2)

arrows(x0=1:2, y0=richpred_size1.1$lci [FF],x1=1:2, y1=richpred_size1.1$uci[FF],angle=90,length=0.1, code=3, lwd=2,col = "black")
points(x = jitter(raw_x2, factor = 1),y = ModelRich2$Size_1, pch = 16, cex = 0.4, col = "grey")

#Age

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Size_1,xlab = expression("Crop Age (Days)"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.2,'b)',cex=1)

axis(side=1, at=seq(from=min(richpred_size1.1$Age_Scaled),to=max(richpred_size1.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_size1.1$Age_Scaled[FFF],rev(richpred_size1.1$Age_Scaled[FFF])), y = c(richpred_size1.1$lci[FFF],rev(richpred_size1.1$uci[FFF])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)

lines(x=richpred_size1.1$Age_Scaled[FFF],y = richpred_size1.1$fit[FFF],lwd = 2,col = 'grey30')

#Day 

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Size_1,xlab = expression("Day Sampled"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -1.7,'c)',cex=1)
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.1,-2.1,1.15,-2.1, length =0.1)

axis(side=1, at=seq(from=min(richpred_size1.1$Day_Scaled),to=max(richpred_size1.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_size1.1$Day_Scaled[FFFF],rev(richpred_size1.1$Day_Scaled[FFFF])), y = c(richpred_size1.1$lci[FFFF],rev(richpred_size1.1$uci[FFFF])),col = rgb(0.5, 0.5, 0.5, 0.5),border= NA)
lines(x=richpred_size1.1$Day_Scaled[FFFF],y = richpred_size1.1$fit[FFFF],lwd = 2,col = 'grey30')

###Size 2 (2.5-5) ?SI----
#Possible for supporting info - it doesn't show much of a change over crop age
summary(Rich_Size2) #Age

richpred_size2.1 <- richpredresults[[8]]
head(richpred_size2.1);dim(richpred_size2.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Size_2,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',main = "2.5-5cm")
axis(side=1, at=seq(from=min(richpred_size2.1$Age_Scaled),to=max(richpred_size2.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(richpred_size2.1$Age_Scaled,rev(richpred_size2.1$Age_Scaled)), y = c(richpred_size2.1$lci,rev(richpred_size2.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_size2.1$Age_Scaled,y = richpred_size2.1$fit,lwd = 2,col = 'grey30')


###Size 3 (5-10) ?SI----
#Possible for supporting info - it doesn't show much and what it does show is very marginal
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

plot(x = 2:1,y = richpred_size3.1$fit [G],xlab = " ",ylab = 'Species Richness', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1,xaxt = "n",ylim = c(0,1), xlim = c(0,3))
mtext(side=3,line=0,at = -0.3,'a)',cex=1)
axis(side=1,at=2:1,labels=c('Outer','Inner'))
mtext(side=3,line=0.1,at = 3.5,expression(bold('5-10cm')),cex=1.2)

arrows(x0=2:1, y0=richpred_size3.1$lci[G],x1=2:1, y1=richpred_size3.1$uci[G],angle=90,length=0.1, code=3, lwd=2,col = "black")
points(x = jitter(raw_x2, factor = 1),y = ModelRich2$Size_3, pch = 16, cex = 0.4, col = "grey")

#Age * Day

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Size_3,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',ylim = c(0,7))
mtext(side=3,line=0,at = -2.3,'b)',cex=1)
axis(side=1, at=seq(from=min(richpred_size3.1$Age_Scaled),to=max(richpred_size3.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_size3.1$Age_Scaled[GG],rev(richpred_size3.1$Age_Scaled[GG])), y = c(richpred_size3.1$lci[GG],rev(richpred_size3.1$uci[GG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_size3.1$Age_Scaled[GG],y = richpred_size3.1$fit[GG],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(richpred_size3.1$Age_Scaled[GGG],rev(richpred_size3.1$Age_Scaled[GGG])), y = c(richpred_size3.1$lci[GGG],rev(richpred_size3.1$uci[GGG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_size3.1$Age_Scaled[GGG],y = richpred_size3.1$fit[GGG],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Water

plot(x = ModelRich2$X1km_Prop_Water,y = ModelRich2$Size_3,xlab = expression("Proportion Water within 1km"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 0.9,'c)',cex=1)

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
plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Develops_Wings,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=3,line=0,at = -2.2,'a)',cex=1)
mtext(side=3,line=0.1,at= 1.5,expression(bold('Develops Wings')),cex=1.2)
axis(side=1, at=seq(from=min(richpred_DevW.1$Age_Scaled),to=max(richpred_DevW.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(richpred_DevW.1$Age_Scaled[H],rev(richpred_DevW.1$Age_Scaled[H])), y = c(richpred_DevW.1$lci[H],rev(richpred_DevW.1$uci[H])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_DevW.1$Age_Scaled[H],y = richpred_DevW.1$fit[H],lwd = 2,col = 'grey30')

#NDVI

plot(x = ModelRich2$NDVIsum_1km,y = ModelRich2$Develops_Wings,xlab = "Total NDVI within 1km",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 3, lwd = 2)
mtext(side=3,line=0,at = 370,'b)',cex=1)

polygon(x = c(richpred_DevW.1$NDVIsum_1km[HH],rev(richpred_DevW.1$NDVIsum_1km[HH])), y = c(richpred_DevW.1$lci[HH],rev(richpred_DevW.1$uci[HH])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_DevW.1$NDVIsum_1km[HH],y = richpred_DevW.1$fit[HH],lwd = 2,col = 'grey30')

###Wingless ?SI----
#Possible for supporting info - is this showing enough to be included in paper? I'm leaning towards no
summary(Rich_Wless) #posiiton * age + day

richpred_Wless.1 <- richpredresults[[11]]
head(richpred_Wless.1);dim(richpred_Wless.1)

I <- richpred_Wless.1$Position == "Inner" & richpred_Wless.1$Day_Scaled == Predictions_Day[10]
II <- richpred_Wless.1$Position == "Outer" & richpred_Wless.1$Day_Scaled == Predictions_Day[10]

III <- richpred_Wless.1$Position == "Outer" & richpred_Wless.1$Age_Scaled == Predictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

#Position * Age 

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Wingless,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=3,line=0,at = -2.2,'a)',cex=1)
mtext(side=3,line=0.1,at = 1.5,expression(bold('Wingless')),cex=1.2)

axis(side=1, at=seq(from=min(richpred_Wless.1$Age_Scaled),to=max(richpred_Wless.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(richpred_Wless.1$Age_Scaled[I],rev(richpred_Wless.1$Age_Scaled[I])), y = c(richpred_Wless.1$lci[I],rev(richpred_Wless.1$uci[I])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_Wless.1$Age_Scaled[I],y = richpred_Wless.1$fit[I],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(richpred_Wless.1$Age_Scaled[II],rev(richpred_Wless.1$Age_Scaled[II])), y = c(richpred_Wless.1$lci[II],rev(richpred_Wless.1$uci[II])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_Wless.1$Age_Scaled[II],y = richpred_Wless.1$fit[II],lwd = 2,lty = 2, col = 'grey30')

legend('topleft',legend = c('Inner', "Outer"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Day

plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Wingless, xlab = "Day Sampled",ylab = 'Species Richness',type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=3,line=0,at = -1.65,'b)',cex=1)
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.2,-1.8,1.2,-1.8, length =0.1)
axis(side=1, at=seq(from=min(richpred_Wless.1$Day_Scaled),to=max(richpred_Wless.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(richpred_Wless.1$Day_Scaled[III],rev(richpred_Wless.1$Day_Scaled[III])), y = c(richpred_Wless.1$lci[III],rev(richpred_Wless.1$uci[III])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_Wless.1$Day_Scaled[III],y = richpred_Wless.1$fit[III],lwd = 2,lty = 1, col = 'grey30')


#DIVERSITY----

##All----

summary(div_all) #Day
divpred_all.1 <- divpredresults[[1]]
head(divpred_all.1);dim(divpred_all.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelDiv2$Day_Scaled,y = ModelDiv2$All,xlab = "Day Sampled",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "All")
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.2,-2.8,1.25,-2.8, length =0.1)
  axis(side=1, at=seq(from=min(divpred_all.1$Day_Scaled),to=max(divpred_all.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Day_Sampled),to=max(ModelDiv2$Day_Sampled),length.out=6),-1))

polygon(x = c(divpred_all.1$Day_Scaled,rev(divpred_all.1$Day_Scaled)), y = c(divpred_all.1$lci,rev(divpred_all.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_all.1$Day_Scaled,y = divpred_all.1$fit,lwd = 2,col = 'grey30')

##Predator----

summary(div_pred) #Day + Age + Water
divpred_pred.1 <- divpredresults[[2]]
head(divpred_pred.1);dim(divpred_pred.1)

Q <- divpred_pred.1$Age_Scaled == Predictions_Age[10] & divpred_pred.1$X1km_Prop_Water == Predictions_Water[10]
QQ <- divpred_pred.1$Day_Scaled == Predictions_Day[10] & divpred_pred.1$X1km_Prop_Water == Predictions_Water[10]
QQQ <- divpred_pred.1$Day_Scaled == Predictions_Day[10] & divpred_pred.1$Age_Scaled == Predictions_Age[10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow=c(2,2))

#Day 

plot(x = ModelDiv2$Day_Scaled,y = ModelDiv2$Predator,xlab = "Day Sampled",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.1,-3,1.1,-3, length =0.1)
axis(side=1, at=seq(from=min(divpred_pred.1$Day_Scaled),to=max(divpred_pred.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Day_Sampled),to=max(ModelDiv2$Day_Sampled),length.out=6),-1),cex.axis = 0.9)
mtext(side=3,line=0,at = -1.7,'a)',cex=1)
mtext(side=3,line=0.1,at = 2,expression(bold('Predator')),cex=1.2)

polygon(x = c(divpred_pred.1$Day_Scaled[Q],rev(divpred_pred.1$Day_Scaled[Q])), y = c(divpred_pred.1$lci[Q],rev(divpred_pred.1$uci[Q])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_pred.1$Day_Scaled[Q],y = divpred_pred.1$fit[Q],lwd = 2,col = 'grey30')

#Age

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Predator,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
axis(side=1, at=seq(from=min(divpred_pred.1$Age_Scaled),to=max(divpred_pred.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)
mtext(side=3,line=0,at = -2.3,'b)',cex=1)

polygon(x = c(divpred_pred.1$Age_Scaled[QQ],rev(divpred_pred.1$Age_Scaled[QQ])), y = c(divpred_pred.1$lci[QQ],rev(divpred_pred.1$uci[QQ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_pred.1$Age_Scaled[QQ],y = divpred_pred.1$fit[QQ],lwd = 2,col = 'grey30',lty = 1)

#Water

plot(x = ModelDiv2$X1km_Prop_Water,y = ModelDiv2$Predator,xlab = "Water within 1km (%)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 0.75,'c)',cex=1)

polygon(x = c(divpred_pred.1$X1km_Prop_Water[QQQ],rev(divpred_pred.1$X1km_Prop_Water[QQQ])), y = c(divpred_pred.1$lci[QQQ],rev(divpred_pred.1$uci[QQQ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_pred.1$X1km_Prop_Water[QQQ],y = divpred_pred.1$fit[QQQ],lwd = 2,col = 'grey30',lty = 1)

##Herbivore----

summary(div_herb) #Position * Day + Water
divpred_herb.1 <- divpredresults[[3]]
head(divpred_herb.1);dim(divpred_herb.1)

R <- divpred_herb.1$Position == "Inner" & divpred_herb.1$X1km_Prop_Water == Predictions_Water[10]
RR <- divpred_herb.1$Position == "Outer" & divpred_herb.1$X1km_Prop_Water == Predictions_Water[10]
RRR <- divpred_herb.1$Position == "Outer" & divpred_herb.1$Day_Scaled == Predictions_Day[10]


dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow = c(1,2))

plot(x = ModelDiv2$Day_Scaled,y = ModelDiv2$Herbivore,xlab = "Day Sampled",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.25,-1.2,1.25,-1.2, length =0.1)
axis(side=1, at=seq(from=min(divpred_herb.1$Day_Scaled),to=max(divpred_herb.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Day_Sampled),to=max(ModelDiv2$Day_Sampled),length.out=6),-1))
mtext(side=3,line=0,at = -1.6,'a)',cex=1)
mtext(side=3,line=0.1,at = 2,expression(bold('Herbivore')),cex=1.2)

polygon(x = c(divpred_herb.1$Day_Scaled[R],rev(divpred_herb.1$Day_Scaled[R])), y = c(divpred_herb.1$lci[R],rev(divpred_herb.1$uci[R])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_herb.1$Day_Scaled[R],y = divpred_herb.1$fit[R],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_herb.1$Day_Scaled[RR],rev(divpred_herb.1$Day_Scaled[RR])), y = c(divpred_herb.1$lci[RR],rev(divpred_herb.1$uci[RR])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_herb.1$Day_Scaled[RR],y = divpred_herb.1$fit[RR],lwd = 2,col = 'grey30',lty = 2)


legend('topleft',legend = c('Inner', "Outer"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Water

plot(x = ModelDiv2$X1km_Prop_Water,y = ModelDiv2$Herbivore,xlab = "Water within 1km (%)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 1,'c)',cex=1)

polygon(x = c(divpred_herb.1$X1km_Prop_Water[RRR],rev(divpred_herb.1$X1km_Prop_Water[RRR])), y = c(divpred_herb.1$lci[RRR],rev(divpred_herb.1$uci[RRR])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_herb.1$X1km_Prop_Water[RRR],y = divpred_herb.1$fit[RRR],lwd = 2,col = 'grey30',lty = 1)

##Fungivore (SI)----
#For the supporting information - despite being top model - this doesn't show much and the CI are massive on some of them as well

summary(div_fung) #Position + Age * Day + NDVI 1km (scaled)
divpred_fung.1 <- divpredresults[[4]]
head(divpred_fung.1);dim(divpred_fung.1)

Predictions_Day[4] #winter
Predictions_Day[15] #Spring

S <- divpred_fung.1$Age_Scaled == Predictions_Age[10] & divpred_fung.1$Day_Scaled == Predictions_Day[10] & divpred_fung.1$NDVI1km_Scaled == Predictions_NDVI1km_Scaled[10]
SS <- divpred_fung.1$Position == "Outer" & divpred_fung.1$Day_Scaled == Predictions_Day[4] & divpred_fung.1$NDVI1km_Scaled == Predictions_NDVI1km_Scaled[10]
SSS <- divpred_fung.1$Position == "Outer" & divpred_fung.1$Day_Scaled == Predictions_Day[15] & divpred_fung.1$NDVI1km_Scaled == Predictions_NDVI1km_Scaled[10]
S_S <- divpred_fung.1$Position == "Outer" & divpred_fung.1$Age_Scaled == Predictions_Age[10] & divpred_fung.1$Day_Scaled == Predictions_Day[10] 

raw_x5 <- ifelse(ModelDiv2$Position == "Inner", 1, ifelse(ModelDiv2$Position == "Outer", 2, NA))

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow = c(2,2))

plot(x = 2:1,y = divpred_fung.1$fit[S],xlab = " ",ylab = 'Diversity', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,2),xaxt = "n",xlim = c(0,3))
mtext(side=3,line=0,at = -0.3,'a)',cex=1)
axis(side=1,at=2:1,labels=c('Outer','Inner'))
mtext(side=3,line=0.1,at = 3.5,expression(bold('Fungivore')),cex=1.2)

arrows(x0=2:1, y0=divpred_fung.1$lci[S],x1=2:1, y1=divpred_fung.1$uci[S],angle=90,length=0.1, code=3, lwd=2,col = "black")

points(x = jitter(raw_x5, factor = 1),y = ModelDiv2$Fungivore, pch = 16, cex = 0.4, col = "grey")

#Age * Day

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Fungivore,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",ylim = c(0,10))
axis(side=1, at=seq(from=min(divpred_fung.1$Age_Scaled),to=max(divpred_fung.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)
mtext(side=3,line=0,at = -2.3,'b)',cex=1)

polygon(x = c(divpred_fung.1$Age_Scaled[SS],rev(divpred_fung.1$Age_Scaled[SS])), y = c(divpred_fung.1$lci[SS],rev(divpred_fung.1$uci[SS])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_fung.1$Age_Scaled[SS],y = divpred_fung.1$fit[SS],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_fung.1$Age_Scaled[SSS],rev(divpred_fung.1$Age_Scaled[SSS])), y = c(divpred_fung.1$lci[SSS],rev(divpred_fung.1$uci[SSS])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_fung.1$Age_Scaled[SSS],y = divpred_fung.1$fit[SSS],lwd = 2,col = 'grey30',lty = 2)

legend('top',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#NDVI 1km

plot(x = ModelDiv2$NDVI1km_Scaled,y = ModelDiv2$Fungivore,xlab = "Total NDVI within 1km",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Develops Wings",ylim=c(0,4))
axis(side=1, at=seq(from=min(divpred_fung.1$NDVI1km_Scaled),to=max(divpred_fung.1$NDVI1km_Scaled),length.out=5),labels=round(seq(from=min(ModelDiv2$NDVIsum_1km),to=max(ModelDiv2$NDVIsum_1km),length.out=5),-1),cex.axis = 0.95)
mtext(side=3,line=0,at = -2.8,'c)',cex=1)

polygon(x = c(divpred_fung.1$NDVI1km_Scaled[S_S],rev(divpred_fung.1$NDVI1km_Scaled[S_S])), y = c(divpred_fung.1$lci[S_S],rev(divpred_fung.1$uci[S_S])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_fung.1$NDVI1km_Scaled[S_S],y = divpred_fung.1$fit[S_S],lwd = 2,col = 'grey30')

##Web Building----

summary(div_web) #Day * Age
divpred_web.1 <- divpredresults[[5]]
head(divpred_web.1);dim(divpred_web.1)

Predictions_Day[4] #winter
Predictions_Day[15] #Spring

TT <- divpred_web.1$Day_Scaled == Predictions_Day[4]
T_T <- divpred_web.1$Day_Scaled == Predictions_Day[15]

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Web,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Web Building")
axis(side=1, at=seq(from=min(divpred_web.1$Age_Scaled),to=max(divpred_web.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(divpred_web.1$Age_Scaled[TT],rev(divpred_web.1$Age_Scaled[TT])), y = c(divpred_web.1$lci[TT],rev(divpred_web.1$uci[TT])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_web.1$Age_Scaled[TT],y = divpred_web.1$fit[TT],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_web.1$Age_Scaled[T_T],rev(divpred_web.1$Age_Scaled[T_T])), y = c(divpred_web.1$lci[T_T],rev(divpred_web.1$uci[T_T])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_web.1$Age_Scaled[T_T],y = divpred_web.1$fit[T_T],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)


##Active Hunting ?SI----
#Possible for supporting info - is this showing enough to be included in paper? I'm not sure

summary(div_active) #Day * Age + Height
divpred_active.1 <- divpredresults[[6]]
head(divpred_active.1);dim(divpred_active.1)

Predictions_Day[4] #winter
Predictions_Day[15] #Spring

U <- divpred_active.1$Day_Scaled == Predictions_Day[4] & divpred_active.1$Height == Predictions_Height[10]
UU <- divpred_active.1$Day_Scaled == Predictions_Day[15] & divpred_active.1$Height == Predictions_Height[10]
UUU <- divpred_active.1$Day_Scaled == Predictions_Day[10] & divpred_active.1$Age_Scaled == Predictions_Age[10]


dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow=c(1,2))

#Day * Age

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Active_Hunting,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",ylim = c(0,7))
axis(side=1, at=seq(from=min(divpred_active.1$Age_Scaled),to=max(divpred_active.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.2,'a)',cex=1)
mtext(side=3,line=0.1,at = 1.6,expression(bold('Active Hunting')),cex=1.2)

polygon(x = c(divpred_active.1$Age_Scaled[U],rev(divpred_active.1$Age_Scaled[U])), y = c(divpred_active.1$lci[U],rev(divpred_active.1$uci[U])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_active.1$Age_Scaled[U],y = divpred_active.1$fit[U],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_active.1$Age_Scaled[UU],rev(divpred_active.1$Age_Scaled[UU])), y = c(divpred_active.1$lci[UU],rev(divpred_active.1$uci[UU])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_active.1$Age_Scaled[UU],y = divpred_active.1$fit[UU],lwd = 2,col = 'grey30',lty = 2)

legend('topright',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Height

plot(x = ModelDiv2$Height,y = ModelDiv2$Active_Hunting,xlab = "Crop Height (cm)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 15,'b)',cex=1)

polygon(x = c(divpred_active.1$Height[UUU],rev(divpred_active.1$Height[UUU])), y = c(divpred_active.1$lci[UUU],rev(divpred_active.1$uci[UUU])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_active.1$Height[UUU],y = divpred_active.1$fit[UUU],lwd = 2,col = 'grey30',lty = 1)

##Size 3 (5-10cm) ?SI----
#Possible for supporting info - is this showing enough to be included in paper? I'm leaning towards no

summary(div_size3) #Day * Age + Water 
divpred_size3.1 <- divpredresults[[7]]
head(divpred_size3.1);dim(divpred_size3.1)

Predictions_Day[4] #winter
Predictions_Day[15] #Spring

V <- divpred_size3.1$Day_Scaled == Predictions_Day[4] & divpred_size3.1$X1km_Prop_Water == Predictions_Water[10]
VV <- divpred_size3.1$Day_Scaled == Predictions_Day[15] & divpred_size3.1$X1km_Prop_Water == Predictions_Water[10]
VVV <- divpred_size3.1$Day_Scaled == Predictions_Day[10] & divpred_size3.1$Age_Scaled == Predictions_Age[10]


dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow=c(1,2))

#Day * Age

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Size_3,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",ylim=c(0,15))
axis(side=1, at=seq(from=min(divpred_size3.1$Age_Scaled),to=max(divpred_size3.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.2,'a)',cex=1)
mtext(side=3,line=0.1,at = 2,expression(bold('5-10cm')),cex=1.2)

polygon(x = c(divpred_size3.1$Age_Scaled[V],rev(divpred_size3.1$Age_Scaled[V])), y = c(divpred_size3.1$lci[V],rev(divpred_size3.1$uci[V])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_size3.1$Age_Scaled[V],y = divpred_size3.1$fit[V],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_size3.1$Age_Scaled[VV],rev(divpred_size3.1$Age_Scaled[VV])), y = c(divpred_size3.1$lci[VV],rev(divpred_size3.1$uci[VV])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_size3.1$Age_Scaled[VV],y = divpred_size3.1$fit[VV],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Water

plot(x = ModelDiv2$X1km_Prop_Water,y = ModelDiv2$Size_3,xlab = "Water within 1km (%)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,ylim=c(0,8))
mtext(side=3,line=0,at = 1,'c)',cex=1)

polygon(x = c(divpred_size3.1$X1km_Prop_Water[VVV],rev(divpred_size3.1$X1km_Prop_Water[VVV])), y = c(divpred_size3.1$lci[VVV],rev(divpred_size3.1$uci[VVV])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_size3.1$X1km_Prop_Water[VVV],y = divpred_size3.1$fit[VVV],lwd = 2,col = 'grey30',lty = 1)


##Develops Wings----

summary(div_DevW) #Field Size (scaled)
divpred_DevW.1 <- divpredresults[[8]]
head(divpred_DevW.1);dim(divpred_DevW.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelDiv2$Fieldsize_Scaled,y = ModelDiv2$Develops_Wings,xlab = "Feild Size (ha)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Develops Wings")
axis(side=1, at=seq(from=min(divpred_DevW.1$Fieldsize_Scaled),to=max(divpred_DevW.1$Fieldsize_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Field_Area_m2),to=max(ModelDiv2$Field_Area_m2),length.out=6)/10000,1))

polygon(x = c(divpred_DevW.1$Fieldsize_Scaled,rev(divpred_DevW.1$Fieldsize_Scaled)), y = c(divpred_DevW.1$lci,rev(divpred_DevW.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_DevW.1$Fieldsize_Scaled,y = divpred_DevW.1$fit,lwd = 2,col = 'grey30')


##Wingless----

summary(div_Wless) #Age
divpred_Wless.1 <- divpredresults[[9]]
head(divpred_Wless.1);dim(divpred_Wless.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Wingless,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Wingless")
axis(side=1, at=seq(from=min(divpred_Wless.1$Age_Scaled),to=max(divpred_Wless.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=6),-1))

polygon(x = c(divpred_Wless.1$Age_Scaled,rev(divpred_Wless.1$Age_Scaled)), y = c(divpred_Wless.1$lci,rev(divpred_Wless.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_Wless.1$Age_Scaled,y = divpred_Wless.1$fit,lwd = 2,col = 'grey30')

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
mtext(side=3,line=0,at = -1.7,'a)',cex=1)
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.1,-0.42,1.1,-0.42, length =0.1)
axis(side=1, at=seq(from=min(occurpred_fung.1$Day_Scaled),to=max(occurpred_fung.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Day_Sampled),to=max(ModelOccur2$Day_Sampled),length.out=6),-1),cex.axis=0.9)
mtext(side=3,line=0.1,at = 2,expression(bold('Fungivore')),cex=1.2)

polygon(x = c(occurpred_fung.1$Day_Scaled[J],rev(occurpred_fung.1$Day_Scaled[J])), y = c(occurpred_fung.1$lci[J],rev(occurpred_fung.1$uci[J])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_fung.1$Day_Scaled[J],y = occurpred_fung.1$fit[J],lwd = 2,col = 'grey30')

#Age

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Fungivore,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.2,'b)',cex=1)
axis(side=1, at=seq(from=min(occurpred_fung.1$Age_Scaled),to=max(occurpred_fung.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1),cex.axis = 0.9)

polygon(x = c(occurpred_fung.1$Age_Scaled[JJ],rev(occurpred_fung.1$Age_Scaled[JJ])), y = c(occurpred_fung.1$lci[JJ],rev(occurpred_fung.1$uci[JJ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_fung.1$Age_Scaled[JJ],y = occurpred_fung.1$fit[JJ],lwd = 2,col = 'grey30')

#Field Size

plot(x = ModelOccur2$FieldSize_Scaled,y = ModelOccur2$Fungivore,xlab = "Field Size (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -2.4,'c)',cex=1)
axis(side=1, at=seq(from=min(occurpred_fung.1$FieldSize_Scaled),to=max(occurpred_fung.1$FieldSize_Scaled),length.out=5),labels=round(seq(from=min(ModelOccur2$Field_Area_m2),to=max(ModelOccur2$Field_Area_m2),length.out=5)/10000,1),cex.axis = 0.9)

polygon(x = c(occurpred_fung.1$FieldSize_Scaled[JJJ],rev(occurpred_fung.1$FieldSize_Scaled[JJJ])), y = c(occurpred_fung.1$lci[JJJ],rev(occurpred_fung.1$uci[JJJ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_fung.1$FieldSize_Scaled[JJJ],y = occurpred_fung.1$fit[JJJ],lwd = 2,col = 'grey30')

##Hematophagous ?SI----
#Possible for supporting info - is this showing enough to be included in paper? I'm leaning towards no

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
mtext(side=3,line=0.1,at = 1.6,expression(bold('Hematophagous')),cex=1.2)
mtext(side=3,line=0,at = -2.2,'a)',cex=1)

polygon(x = c(occurpred_hema.1$Age_Scaled[K],rev(occurpred_hema.1$Age_Scaled[K])), y = c(occurpred_hema.1$lci[K],rev(occurpred_hema.1$uci[K])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hema.1$Age_Scaled[K],y = occurpred_hema.1$fit[K],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(occurpred_hema.1$Age_Scaled[KK],rev(occurpred_hema.1$Age_Scaled[KK])), y = c(occurpred_hema.1$lci[KK],rev(occurpred_hema.1$uci[KK])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hema.1$Age_Scaled[KK],y = occurpred_hema.1$fit[KK],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#NDVI Field

plot(x = ModelOccur2$NDVImean_Field,y = ModelOccur2$Hematophagous,xlab = "Mean Field NDVI",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2)
mtext(side=3,line=0,at = 0.12,'b)',cex=1)

polygon(x = c(occurpred_hema.1$NDVImean_Field[K_K],rev(occurpred_hema.1$NDVImean_Field[K_K])), y = c(occurpred_hema.1$lci[K_K],rev(occurpred_hema.1$uci[K_K])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hema.1$NDVImean_Field[K_K],y = occurpred_hema.1$fit[K_K],lwd = 2,col = 'grey30',lty = 1)



##Web Building----

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

##Hawking (SI)----
#For the supporting information - despite being top model - this doesn't fit well at all, we can see that visually for both panels

summary(Occur_hawk) #Position + Age * Day 
occurpred_hawk.1 <- occurpredresults[[7]]
head(occurpred_hawk.1);dim(occurpred_hawk.1)

Predictions_Day[4] #Winter
Predictions_Day[15] #Spring

M <- occurpred_hawk.1$Age_Scaled == Predictions_Age[10] & occurpred_hawk.1$Day_Scaled == Predictions_Day [10]
MM <- occurpred_hawk.1$Position == "Outer" & occurpred_hawk.1$Day_Scaled == Predictions_Day[4]
MMM <- occurpred_hawk.1$Position == "Outer" & occurpred_hawk.1$Day_Scaled == Predictions_Day[15]


raw_x4 <- ifelse(ModelOccur2$Position == "Inner", 1, 
                ifelse(ModelOccur2$Position == "Outer", 2, NA))

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow = c(1,2))

plot(x = 2:1,y = occurpred_hawk.1$fit[M],xlab = " ",ylab = 'Probability of Occurrence', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,1),xaxt = "n",xlim = c(0,3))
mtext(side=3,line=0,at = -0.3,'a)',cex=1)
axis(side=1,at=2:1,labels=c('Outer','Inner'))
mtext(side=3,line=0.1,at = 3.5,expression(bold('Hawking')),cex=1.2)

arrows(x0=2:1, y0=occurpred_hawk.1$lci[M],x1=2:1, y1=occurpred_hawk.1$uci[M],angle=90,length=0.2, code=3, lwd=2,col = "black")

points(x = jitter(raw_x4, factor = 1),y = ModelOccur2$Hawking, pch = 16, cex = 0.4, col = "grey")

#Age * Day

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Hawking,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
axis(side=1, at=seq(from=min(occurpred_hawk.1$Age_Scaled),to=max(occurpred_hawk.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0,at = -2.2,'b)',cex=1)

polygon(x = c(occurpred_hawk.1$Age_Scaled[MM],rev(occurpred_hawk.1$Age_Scaled[MM])), y = c(occurpred_hawk.1$lci[MM],rev(occurpred_hawk.1$uci[MM])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hawk.1$Age_Scaled[MM],y = occurpred_hawk.1$fit[MM],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(occurpred_hawk.1$Age_Scaled[MMM],rev(occurpred_hawk.1$Age_Scaled[MMM])), y = c(occurpred_hawk.1$lci[MMM],rev(occurpred_hawk.1$uci[MMM])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hawk.1$Age_Scaled[MMM],y = occurpred_hawk.1$fit[MMM],lwd = 2,col = 'grey30',lty = 2)

legend('topright',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)



##Size 2 (2.5-5cm)----

summary(Occur_size2) #Height 
occurpred_size2.1 <- occurpredresults[[8]]
head(occurpred_size2.1);dim(occurpred_size2.1)

dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = ModelOccur2$Height,y = ModelOccur2$Size_2,xlab = "Plant Height (cm)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,main = "2.5-5cm")

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

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Size_3,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "5-10cm")
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

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Size_4,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = ">10cm")
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

O <- occurpred_Wless.1$Position == "Inner" & occurpred_Wless.1$Day_Scaled == Predictions_Day[10]
OO <- occurpred_Wless.1$Position == "Outer" & occurpred_Wless.1$Day_Scaled == Predictions_Day[10]
OOO <- occurpred_Wless.1$Position == "Outer" & occurpred_Wless.1$Age_Scaled == Predictions_Age[10]

dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mgp=c(2.5,1,0),xpd = T,mfrow=c(1,2))

#position * Age

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Wingless,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
axis(side=1, at=seq(from=min(occurpred_Wless.1$Age_Scaled),to=max(occurpred_Wless.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=6),-1))
mtext(side=3,line=0.1,at = 1.8,expression(bold('Wingless')),cex=1.2)
mtext(side=3,line=0,at = -2.2,'a)',cex=1)

polygon(x = c(occurpred_Wless.1$Age_Scaled[O],rev(occurpred_Wless.1$Age_Scaled[O])), y = c(occurpred_Wless.1$lci[O],rev(occurpred_Wless.1$uci[O])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_Wless.1$Age_Scaled[O],y = occurpred_Wless.1$fit[O],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(occurpred_Wless.1$Age_Scaled[OO],rev(occurpred_Wless.1$Age_Scaled[OO])), y = c(occurpred_Wless.1$lci[OO],rev(occurpred_Wless.1$uci[OO])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_Wless.1$Age_Scaled[OO],y = occurpred_Wless.1$fit[OO],lwd = 2,col = 'grey30',lty = 2)

legend('bottomright',legend = c('Inner', "Outer"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Day

plot(x = ModelOccur2$Day_Scaled,y = ModelOccur2$Wingless,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n")
mtext(side=3,line=0,at = -1.65,'b)',cex=1)
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1.2,-0.3,1.25,-0.3, length =0.1)
axis(side=1, at=seq(from=min(occurpred_Wless.1$Day_Scaled),to=max(occurpred_Wless.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelOccur2$Day_Sampled),to=max(ModelOccur2$Day_Sampled),length.out=6),-1))

polygon(x = c(occurpred_Wless.1$Day_Scaled[OOO],rev(occurpred_Wless.1$Day_Scaled[OOO])), y = c(occurpred_Wless.1$lci[OOO],rev(occurpred_Wless.1$uci[OOO])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=occurpred_Wless.1$Day_Scaled[OOO],y = occurpred_Wless.1$fit[OOO],lwd = 2,col = 'grey30')


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
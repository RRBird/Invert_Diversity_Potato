
#RICHNESS FIGURE----

#DIVERSITY FIGURE ----

dev.new(height=15,width=15,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(3,3),mgp=c(2.5,1,0),xpd = T)

#ALL

plot(x = ModelDiv2$Day_Scaled,y = ModelDiv2$All,xlab = "Day Sampled",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "All",cex.lab=1.2)
mtext(side=1,line=3,at = -1.8,'Winter',cex=0.7)
mtext(side=1,line=3,at = 1.8,'Spring',cex=0.7)
arrows(-1.2,-5.3,1.25,-5.3, length =0.1)
mtext(side=3,line=0,at = -1.7,'a)',cex=0.9)
axis(side=1, at=seq(from=min(divpred_all.1$Day_Scaled),to=max(divpred_all.1$Day_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Day_Sampled),to=max(ModelDiv2$Day_Sampled),length.out=4),0))

polygon(x = c(divpred_all.1$Day_Scaled,rev(divpred_all.1$Day_Scaled)), y = c(divpred_all.1$lci,rev(divpred_all.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_all.1$Day_Scaled,y = divpred_all.1$fit,lwd = 2,col = 'grey30')

#PREDATOR

#Day 

plot(x = ModelDiv2$Day_Scaled,y = ModelDiv2$Predator,xlab = "Day Sampled",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Predator (1)",cex.lab=1.2)
mtext(side=1,line=3,at = -1.8,'Winter',cex=0.7)
mtext(side=1,line=3,at = 1.8,'Spring',cex=0.7)
arrows(-1.2,-4.1,1.25,-4.1, length =0.1)
axis(side=1, at=seq(from=min(divpred_pred.1$Day_Scaled),to=max(divpred_pred.1$Day_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Day_Sampled),to=max(ModelDiv2$Day_Sampled),length.out=4),0))
mtext(side=3,line=0,at = -1.7,'b)',cex=0.8)

polygon(x = c(divpred_pred.1$Day_Scaled[Q],rev(divpred_pred.1$Day_Scaled[Q])), y = c(divpred_pred.1$lci[Q],rev(divpred_pred.1$uci[Q])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_pred.1$Day_Scaled[Q],y = divpred_pred.1$fit[Q],lwd = 2,col = 'grey30')

#Age

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Predator,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Predator (2)",cex.lab=1.2)
axis(side=1, at=seq(from=min(divpred_pred.1$Age_Scaled),to=max(divpred_pred.1$Age_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=4),0))
mtext(side=3,line=0,at = -2.3,'c)',cex=0.8)

polygon(x = c(divpred_pred.1$Age_Scaled[QQ],rev(divpred_pred.1$Age_Scaled[QQ])), y = c(divpred_pred.1$lci[QQ],rev(divpred_pred.1$uci[QQ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_pred.1$Age_Scaled[QQ],y = divpred_pred.1$fit[QQ],lwd = 2,col = 'grey30',lty = 1)

#Water

plot(x = ModelDiv2$X1km_Prop_Water,y = ModelDiv2$Predator,xlab = "Water within 1km (%)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,,main = "Predator (3)",cex.lab=1.2)
mtext(side=3,line=0,at = 0.75,'d)',cex=0.8)

polygon(x = c(divpred_pred.1$X1km_Prop_Water[QQQ],rev(divpred_pred.1$X1km_Prop_Water[QQQ])), y = c(divpred_pred.1$lci[QQQ],rev(divpred_pred.1$uci[QQQ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_pred.1$X1km_Prop_Water[QQQ],y = divpred_pred.1$fit[QQQ],lwd = 2,col = 'grey30',lty = 1)

#HERBIVORE

#Position * Day

plot(x = ModelDiv2$Day_Scaled,y = ModelDiv2$Herbivore,xlab = "Day Sampled",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Herbivore (1)",cex.lab=1.2)
mtext(side=1,line=3,at = -1.8,'Winter',cex=0.7)
mtext(side=1,line=3,at = 1.8,'Spring',cex=0.7)
arrows(-1.25,-2.2,1.25,-2.2, length =0.1)
axis(side=1, at=seq(from=min(divpred_herb.1$Day_Scaled),to=max(divpred_herb.1$Day_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Day_Sampled),to=max(ModelDiv2$Day_Sampled),length.out=4),0))
mtext(side=3,line=0,at = -1.7,'e)',cex=0.8)

polygon(x = c(divpred_herb.1$Day_Scaled[R],rev(divpred_herb.1$Day_Scaled[R])), y = c(divpred_herb.1$lci[R],rev(divpred_herb.1$uci[R])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_herb.1$Day_Scaled[R],y = divpred_herb.1$fit[R],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_herb.1$Day_Scaled[RR],rev(divpred_herb.1$Day_Scaled[RR])), y = c(divpred_herb.1$lci[RR],rev(divpred_herb.1$uci[RR])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_herb.1$Day_Scaled[RR],y = divpred_herb.1$fit[RR],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c("Outer",'Inner'), lty = c(2,1), col = 'grey30',pt.cex = 1)

#Water

plot(x = ModelDiv2$X1km_Prop_Water,y = ModelDiv2$Herbivore,xlab = "Water within 1km (%)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,main = "Herbivore (2)",cex.lab=1.2)
mtext(side=3,line=0,at = 1,'f)',cex=0.8)

polygon(x = c(divpred_herb.1$X1km_Prop_Water[RRR],rev(divpred_herb.1$X1km_Prop_Water[RRR])), y = c(divpred_herb.1$lci[RRR],rev(divpred_herb.1$uci[RRR])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_herb.1$X1km_Prop_Water[RRR],y = divpred_herb.1$fit[RRR],lwd = 2,col = 'grey30',lty = 1)

#WEB BUILDING

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Web,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Web Building",cex.lab=1.2)
axis(side=1, at=seq(from=min(divpred_web.1$Age_Scaled),to=max(divpred_web.1$Age_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=4),0))
mtext(side=3,line=0,at = -2.3,'g)',cex=0.8)

polygon(x = c(divpred_web.1$Age_Scaled[TT],rev(divpred_web.1$Age_Scaled[TT])), y = c(divpred_web.1$lci[TT],rev(divpred_web.1$uci[TT])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_web.1$Age_Scaled[TT],y = divpred_web.1$fit[TT],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_web.1$Age_Scaled[T_T],rev(divpred_web.1$Age_Scaled[T_T])), y = c(divpred_web.1$lci[T_T],rev(divpred_web.1$uci[T_T])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_web.1$Age_Scaled[T_T],y = divpred_web.1$fit[T_T],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#DEVELOPS WINGS

plot(x = ModelDiv2$Fieldsize_Scaled,y = ModelDiv2$Develops_Wings,xlab = "Field Size (ha)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Develops Wings",cex.lab=1.2)
axis(side=1, at=seq(from=min(divpred_DevW.1$Fieldsize_Scaled),to=max(divpred_DevW.1$Fieldsize_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Field_Area_m2),to=max(ModelDiv2$Field_Area_m2),length.out=4)/10000,1))
mtext(side=3,line=0,at = -2.5,'h)',cex=0.8)

polygon(x = c(divpred_DevW.1$Fieldsize_Scaled,rev(divpred_DevW.1$Fieldsize_Scaled)), y = c(divpred_DevW.1$lci,rev(divpred_DevW.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_DevW.1$Fieldsize_Scaled,y = divpred_DevW.1$fit,lwd = 2,col = 'grey30')

#WINGLESS

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Wingless,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Wingless",cex.lab=1.2)
axis(side=1, at=seq(from=min(divpred_Wless.1$Age_Scaled),to=max(divpred_Wless.1$Age_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=4),0))
mtext(side=3,line=0,at = -2.3,'I)',cex=0.8)

polygon(x = c(divpred_Wless.1$Age_Scaled,rev(divpred_Wless.1$Age_Scaled)), y = c(divpred_Wless.1$lci,rev(divpred_Wless.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_Wless.1$Age_Scaled,y = divpred_Wless.1$fit,lwd = 2,col = 'grey30')

#OCCURRENCE FIGURE ----

dev.new(height=15,width=25,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(3,5),mgp=c(2.5,1,0),xpd = T)

#HERBIVORE

plot(x = ModelOccur2$Day_Scaled,y = ModelOccur2$Herbivore,xlab = "Day Sampled",ylab = 'Probability of Occuring', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Herbivore",cex.lab=1.2)
mtext(side=1,line=3,at = -1.8,'Winter',cex=0.7)
mtext(side=1,line=3,at = 1.8,'Spring',cex=0.7)
arrows(-1.2,-0.55,1.25,-0.55, length =0.1)
axis(side=1, at=seq(from=min(occurpred_herb.1$Day_Scaled),to=max(occurpred_herb.1$Day_Scaled),length.out=4),labels=round(seq(from=min(ModelOccur2$Day_Sampled),to=max(ModelOccur2$Day_Sampled),length.out=4),0))
mtext(side=3,line=0,at = -1.7,'a)',cex=0.8)

polygon(x = c(occurpred_herb.1$Day_Scaled,rev(occurpred_herb.1$Day_Scaled)), y = c(occurpred_herb.1$lci,rev(occurpred_herb.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_herb.1$Day_Scaled,y = occurpred_herb.1$fit,lwd = 2,col = 'grey30')

#OMNIVORE

plot(x = ModelOccur2$GC,y = ModelOccur2$Omnivore,xlab = "Ground Cover (%)",ylab = 'Probability of Occuring', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,main = "Omnivore",cex.lab=1.2,xaxt="n")
axis(1, at = round(seq(min(ModelOccur2$GC), max(ModelOccur2$GC), length.out = 4),0))
mtext(side=3,line=0,at = 1,'b)',cex=0.8)

polygon(x = c(occurpred_omni.1$GC,rev(occurpred_omni.1$GC)), y = c(occurpred_omni.1$lci,rev(occurpred_omni.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_omni.1$GC,y = occurpred_omni.1$fit,lwd = 2,col = 'grey30')

#FUNGIVORE

#Day

plot(x = ModelOccur2$Day_Scaled,y = ModelOccur2$Fungivore,xlab = "Day Sampled",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Fungivore (1)",cex.lab=1.2)
mtext(side=3,line=0,at = -1.7,'c)',cex=0.8)
mtext(side=1,line=3,at = -1.8,'Winter',cex=0.7)
mtext(side=1,line=3,at = 1.8,'Spring',cex=0.7)
arrows(-1.2,-0.54,1.2,-0.54, length =0.1)
axis(side=1, at=seq(from=min(occurpred_fung.1$Day_Scaled),to=max(occurpred_fung.1$Day_Scaled),length.out=4),labels=round(seq(from=min(ModelOccur2$Day_Sampled),to=max(ModelOccur2$Day_Sampled),length.out=4),0))

polygon(x = c(occurpred_fung.1$Day_Scaled[J],rev(occurpred_fung.1$Day_Scaled[J])), y = c(occurpred_fung.1$lci[J],rev(occurpred_fung.1$uci[J])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_fung.1$Day_Scaled[J],y = occurpred_fung.1$fit[J],lwd = 2,col = 'grey30')

#Age

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Fungivore,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Fungivore (2)",cex.lab=1.2)
mtext(side=3,line=0,at = -2.2,'d)',cex=0.8)
axis(side=1, at=seq(from=min(occurpred_fung.1$Age_Scaled),to=max(occurpred_fung.1$Age_Scaled),length.out=4),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=4),0))

polygon(x = c(occurpred_fung.1$Age_Scaled[JJ],rev(occurpred_fung.1$Age_Scaled[JJ])), y = c(occurpred_fung.1$lci[JJ],rev(occurpred_fung.1$uci[JJ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_fung.1$Age_Scaled[JJ],y = occurpred_fung.1$fit[JJ],lwd = 2,col = 'grey30')

#Field Size

plot(x = ModelOccur2$FieldSize_Scaled,y = ModelOccur2$Fungivore,xlab = "Field Size (ha)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Fungivore (3)",cex.lab=1.2)
mtext(side=3,line=0,at = -2.4,'e)',cex=0.8)
axis(side=1, at=seq(from=min(occurpred_fung.1$FieldSize_Scaled),to=max(occurpred_fung.1$FieldSize_Scaled),length.out=4),labels=round(seq(from=min(ModelOccur2$Field_Area_m2),to=max(ModelOccur2$Field_Area_m2),length.out=4)/10000,1))

polygon(x = c(occurpred_fung.1$FieldSize_Scaled[JJJ],rev(occurpred_fung.1$FieldSize_Scaled[JJJ])), y = c(occurpred_fung.1$lci[JJJ],rev(occurpred_fung.1$uci[JJJ])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_fung.1$FieldSize_Scaled[JJJ],y = occurpred_fung.1$fit[JJJ],lwd = 2,col = 'grey30')

#WEB BUILDING

#SUPPORTING----

=======================================================================
##Richness----
=======================================================================

dev.new(height=10,width=15,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,3),mgp=c(2.5,1,0),xpd = T)

#WINGLESS

##Position * Age 
plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Wingless,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',cex.lab=1.3,cex.axis=1.2, main = "Wingless (1)")
mtext(side=3,line=0,at = -2.2,'a)',cex=1)

axis(side=1, at=seq(from=min(richpred_Wless.1$Age_Scaled),to=max(richpred_Wless.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),0),cex.axis =1.2)

polygon(x = c(richpred_Wless.1$Age_Scaled[I],rev(richpred_Wless.1$Age_Scaled[I])), y = c(richpred_Wless.1$lci[I],rev(richpred_Wless.1$uci[I])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_Wless.1$Age_Scaled[I],y = richpred_Wless.1$fit[I],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(richpred_Wless.1$Age_Scaled[II],rev(richpred_Wless.1$Age_Scaled[II])), y = c(richpred_Wless.1$lci[II],rev(richpred_Wless.1$uci[II])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_Wless.1$Age_Scaled[II],y = richpred_Wless.1$fit[II],lwd = 2,lty = 2, col = 'grey30')

legend('topleft',legend = c('Inner', "Outer"), lty = c(1,2), col = 'grey30',pt.cex = 1,cex = 1.2)

##Day
plot(x = ModelRich2$Day_Scaled,y = ModelRich2$Wingless, xlab = "Day Sampled",ylab = 'Species Richness',type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',cex.lab = 1.3,cex.axis = 1.2,main = "Wingless (2)")
mtext(side=3,line=0,at = -1.65,'b)',cex=1)
mtext(side=1,line=3,at = -1.5,'Winter',cex=0.9)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.9)
arrows(-1.1,-1.8,1.1,-1.8, length =0.1)
axis(side=1, at=seq(from=min(richpred_Wless.1$Day_Scaled),to=max(richpred_Wless.1$Day_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Day_Sampled),to=max(ModelRich2$Day_Sampled),length.out=6),0),cex.axis =1.2)

polygon(x = c(richpred_Wless.1$Day_Scaled[III],rev(richpred_Wless.1$Day_Scaled[III])), y = c(richpred_Wless.1$lci[III],rev(richpred_Wless.1$uci[III])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=richpred_Wless.1$Day_Scaled[III],y = richpred_Wless.1$fit[III],lwd = 2,lty = 1, col = 'grey30')

#SIZE 2

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Size_2,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',cex.lab=1.3,cex.axis=1.2,main = '2.5-5cm')
axis(side=1, at=seq(from=min(richpred_size2.1$Age_Scaled),to=max(richpred_size2.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),0), cex.axis=1.2)
mtext(side=3,line=0,at = -2.2,'c)',cex=1)

polygon(x = c(richpred_size2.1$Age_Scaled,rev(richpred_size2.1$Age_Scaled)), y = c(richpred_size2.1$lci,rev(richpred_size2.1$uci)),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_size2.1$Age_Scaled,y = richpred_size2.1$fit,lwd = 2,col = 'grey30')

#SIZE 3

##position

plot(x = 2:1,y = richpred_size3.1$fit [G],xlab = " ",ylab = 'Species Richness', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1,xaxt = "n",ylim = c(0,1), xlim = c(0,3),cex.lab=1.3,cex.axis=1.2,main = "5-10cm (1)")
mtext(side=3,line=0,at = -0.3,'d)',cex=1)
axis(side=1,at=2:1,labels=c('Outer','Inner'),cex.axis=1.2)

arrows(x0=2:1, y0=richpred_size3.1$lci[G],x1=2:1, y1=richpred_size3.1$uci[G],angle=90,length=0.1, code=3, lwd=2,col = "black")
points(x = jitter(raw_x2, factor = 1),y = ModelRich2$Size_3, pch = 16, cex = 0.4, col = "grey")

##Age * Day

plot(x = ModelRich2$Age_Scaled,y = ModelRich2$Size_3,xlab = "Crop Age (Days)",ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n',ylim = c(0,7),cex.lab=1.3,cex.axis=1.2, main = "5-10cm (2)")
mtext(side=3,line=0,at = -2.3,'e)',cex=1)
axis(side=1, at=seq(from=min(richpred_size3.1$Age_Scaled),to=max(richpred_size3.1$Age_Scaled),length.out=6),labels=round(seq(from=min(ModelRich2$Crop_Age_Days),to=max(ModelRich2$Crop_Age_Days),length.out=6),0),cex.axis = 1.2)

polygon(x = c(richpred_size3.1$Age_Scaled[GG],rev(richpred_size3.1$Age_Scaled[GG])), y = c(richpred_size3.1$lci[GG],rev(richpred_size3.1$uci[GG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_size3.1$Age_Scaled[GG],y = richpred_size3.1$fit[GG],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(richpred_size3.1$Age_Scaled[GGG],rev(richpred_size3.1$Age_Scaled[GGG])), y = c(richpred_size3.1$lci[GGG],rev(richpred_size3.1$uci[GGG])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=richpred_size3.1$Age_Scaled[GGG],y = richpred_size3.1$fit[GGG],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1,cex=1.2)

##Water

plot(x = ModelRich2$X1km_Prop_Water,y = ModelRich2$Size_3,xlab = expression("Proportion Water within 1km"),ylab = 'Species Richness', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,cex.lab=1.3,cex.axis=1.2,main = "5-10cm (3)")
mtext(side=3,line=0,at = 0.9,'f)',cex=1)

polygon(x = c(richpred_size3.1$X1km_Prop_Water[GGGG],rev(richpred_size3.1$X1km_Prop_Water[GGGG])), y = c(richpred_size3.1$lci[GGGG],rev(richpred_size3.1$uci[GGGG])),col = rgb(0.5, 0.5, 0.5, 0.5),border=NA)
lines(x=richpred_size3.1$X1km_Prop_Water[GGGG],y = richpred_size3.1$fit[GGGG],lwd = 2,col = 'grey30')

=======================================================================
##Diversity----
=======================================================================

dev.new(height=15,width=15,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(3,3),mgp=c(2.5,1,0),xpd = T)

#FUNGIVORE

#Position

plot(x = 2:1,y = divpred_fung.1$fit[S],xlab = " ",ylab = 'Diversity', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,2),xaxt = "n",xlim = c(0,3),main = "Fungivore (1)",cex.lab = 1.2)
mtext(side=3,line=0,at = -0.3,'a)',cex=0.8)
axis(side=1,at=2:1,labels=c(' ',' '))
mtext(side=1,line=0.5,at = 0.9,'Inner',cex=0.8)
mtext(side=1,line=0.5,at = 2.1,'Outer',cex=0.8)

arrows(x0=2:1, y0=divpred_fung.1$lci[S],x1=2:1, y1=divpred_fung.1$uci[S],angle=90,length=0.1, code=3, lwd=2,col = "black")
points(x = jitter(raw_x5, factor = 1),y = ModelDiv2$Fungivore, pch = 16, cex = 0.4, col = "grey")

#Age * Day

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Fungivore,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",ylim = c(0,10),main = "Fungivore (2)",cex.lab = 1.2)
axis(side=1, at=seq(from=min(divpred_fung.1$Age_Scaled),to=max(divpred_fung.1$Age_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=4),0))
mtext(side=3,line=0,at = -2.3,'b)',cex=0.8)

polygon(x = c(divpred_fung.1$Age_Scaled[SS],rev(divpred_fung.1$Age_Scaled[SS])), y = c(divpred_fung.1$lci[SS],rev(divpred_fung.1$uci[SS])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_fung.1$Age_Scaled[SS],y = divpred_fung.1$fit[SS],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_fung.1$Age_Scaled[SSS],rev(divpred_fung.1$Age_Scaled[SSS])), y = c(divpred_fung.1$lci[SSS],rev(divpred_fung.1$uci[SSS])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_fung.1$Age_Scaled[SSS],y = divpred_fung.1$fit[SSS],lwd = 2,col = 'grey30',lty = 2)

legend('top',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#NDVI 1km

plot(x = ModelDiv2$NDVI1km_Scaled,y = ModelDiv2$Fungivore,xlab = "Total NDVI within 1km",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",ylim=c(0,4), main = "Fungivore (3)",cex.lab = 1.2)
axis(side=1, at=seq(from=min(divpred_fung.1$NDVI1km_Scaled),to=max(divpred_fung.1$NDVI1km_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$NDVIsum_1km),to=max(ModelDiv2$NDVIsum_1km),length.out=4),0),cex.axis=0.9)
mtext(side=3,line=0,at = -2.8,'c)',cex=0.8)

polygon(x = c(divpred_fung.1$NDVI1km_Scaled[S_S],rev(divpred_fung.1$NDVI1km_Scaled[S_S])), y = c(divpred_fung.1$lci[S_S],rev(divpred_fung.1$uci[S_S])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_fung.1$NDVI1km_Scaled[S_S],y = divpred_fung.1$fit[S_S],lwd = 2,col = 'grey30')

#ACTIVE HUNTING

#Day * Age

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Active_Hunting,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",ylim = c(0,7),main = "Active Hunting (1)",cex.lab = 1.2)
axis(side=1, at=seq(from=min(divpred_active.1$Age_Scaled),to=max(divpred_active.1$Age_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=4),0))
mtext(side=3,line=0,at = -2.7,'d)',cex=0.8)

polygon(x = c(divpred_active.1$Age_Scaled[U],rev(divpred_active.1$Age_Scaled[U])), y = c(divpred_active.1$lci[U],rev(divpred_active.1$uci[U])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_active.1$Age_Scaled[U],y = divpred_active.1$fit[U],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_active.1$Age_Scaled[UU],rev(divpred_active.1$Age_Scaled[UU])), y = c(divpred_active.1$lci[UU],rev(divpred_active.1$uci[UU])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_active.1$Age_Scaled[UU],y = divpred_active.1$fit[UU],lwd = 2,col = 'grey30',lty = 2)

legend('topright',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Height

plot(x = ModelDiv2$Height,y = ModelDiv2$Active_Hunting,xlab = "Crop Height (cm)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,main = "Active Hunting (2)",,cex.lab = 1.2,xaxt = "n")
axis(1, at = round(seq(min(ModelDiv2$Height), max(ModelDiv2$Height), length.out = 4),0))
mtext(side=3,line=0,at = 7,'e)',cex=0.8)

polygon(x = c(divpred_active.1$Height[UUU],rev(divpred_active.1$Height[UUU])), y = c(divpred_active.1$lci[UUU],rev(divpred_active.1$uci[UUU])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_active.1$Height[UUU],y = divpred_active.1$fit[UUU],lwd = 2,col = 'grey30',lty = 1)

#SIZE 3

#Day * Age

plot(x = ModelDiv2$Age_Scaled,y = ModelDiv2$Size_3,xlab = "Crop Age (Days)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",ylim=c(0,15),main ="5-10cm (1)",cex.lab=1.2)
axis(side=1, at=seq(from=min(divpred_size3.1$Age_Scaled),to=max(divpred_size3.1$Age_Scaled),length.out=4),labels=round(seq(from=min(ModelDiv2$Crop_Age_Days),to=max(ModelDiv2$Crop_Age_Days),length.out=4),0))
mtext(side=3,line=0,at = -2.2,'f)',cex=0.8)

polygon(x = c(divpred_size3.1$Age_Scaled[V],rev(divpred_size3.1$Age_Scaled[V])), y = c(divpred_size3.1$lci[V],rev(divpred_size3.1$uci[V])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_size3.1$Age_Scaled[V],y = divpred_size3.1$fit[V],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(divpred_size3.1$Age_Scaled[VV],rev(divpred_size3.1$Age_Scaled[VV])), y = c(divpred_size3.1$lci[VV],rev(divpred_size3.1$uci[VV])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_size3.1$Age_Scaled[VV],y = divpred_size3.1$fit[VV],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#Water

plot(x = ModelDiv2$X1km_Prop_Water,y = ModelDiv2$Size_3,xlab = "Water within 1km (%)",ylab = 'Diversity', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,ylim=c(0,8),main = "5-10cm (2)",cex.lab=1.2)
mtext(side=3,line=0,at = 1,'g)',cex=0.8)

polygon(x = c(divpred_size3.1$X1km_Prop_Water[VVV],rev(divpred_size3.1$X1km_Prop_Water[VVV])), y = c(divpred_size3.1$lci[VVV],rev(divpred_size3.1$uci[VVV])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=divpred_size3.1$X1km_Prop_Water[VVV],y = divpred_size3.1$fit[VVV],lwd = 2,col = 'grey30',lty = 1)

=======================================================================
##Occurrence----
=======================================================================

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

#HEMATOPHAGOUS

#Day * Age

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Hematophagous,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main = "Hematophagous (1)")
axis(side=1, at=seq(from=min(occurpred_hema.1$Age_Scaled),to=max(occurpred_hema.1$Age_Scaled),length.out=5),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=5),0))
mtext(side=3,line=0,at = -2.2,'a)',cex=0.9)

polygon(x = c(occurpred_hema.1$Age_Scaled[K],rev(occurpred_hema.1$Age_Scaled[K])), y = c(occurpred_hema.1$lci[K],rev(occurpred_hema.1$uci[K])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hema.1$Age_Scaled[K],y = occurpred_hema.1$fit[K],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(occurpred_hema.1$Age_Scaled[KK],rev(occurpred_hema.1$Age_Scaled[KK])), y = c(occurpred_hema.1$lci[KK],rev(occurpred_hema.1$uci[KK])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hema.1$Age_Scaled[KK],y = occurpred_hema.1$fit[KK],lwd = 2,col = 'grey30',lty = 2)

legend('topleft',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)

#NDVI Field

plot(x = ModelOccur2$NDVImean_Field,y = ModelOccur2$Hematophagous,xlab = "Mean Field NDVI",ylab = 'Probability of Occurrence',type = 'p', pch = 16,cex =0.2,col='black', las = 1, lwd = 2,main = "Hematophagous (2)", xaxt = "n")
axis(1, at = round(seq(min(ModelOccur2$NDVImean_Field), max(ModelOccur2$NDVImean_Field), length.out = 4),2))
mtext(side=3,line=0,at = 0.12,'b)',cex=0.9)

polygon(x = c(occurpred_hema.1$NDVImean_Field[K_K],rev(occurpred_hema.1$NDVImean_Field[K_K])), y = c(occurpred_hema.1$lci[K_K],rev(occurpred_hema.1$uci[K_K])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hema.1$NDVImean_Field[K_K],y = occurpred_hema.1$fit[K_K],lwd = 2,col = 'grey30',lty = 1)

#HAWKING

#Position

plot(x = 2:1,y = occurpred_hawk.1$fit[M],xlab = " ",ylab = 'Probability of Occurrence', type = 'p',pch = 16,cex =2.5,col = 'black', las = 1, ylim=c(0,1),xaxt = "n",xlim = c(0,3),main = "Hawking (1)")
mtext(side=3,line=0,at = -0.3,'c)',cex=0.9)
axis(side=1,at=2:1,labels=c('Outer','Inner'))

arrows(x0=2:1, y0=occurpred_hawk.1$lci[M],x1=2:1, y1=occurpred_hawk.1$uci[M],angle=90,length=0.2, code=3, lwd=2,col = "black")
points(x = jitter(raw_x4, factor = 1),y = ModelOccur2$Hawking, pch = 16, cex = 0.4, col = "grey")

#Age * Day

plot(x = ModelOccur2$Age_Scaled,y = ModelOccur2$Hawking,xlab = "Crop Age (Days)",ylab = 'Probability of Occurrence', type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,xaxt ="n",main= "Hawking(2)")
axis(side=1, at=seq(from=min(occurpred_hawk.1$Age_Scaled),to=max(occurpred_hawk.1$Age_Scaled),length.out=5),labels=round(seq(from=min(ModelOccur2$Crop_Age_Days),to=max(ModelOccur2$Crop_Age_Days),length.out=5),0))
mtext(side=3,line=0,at = -2.2,'d)',cex=0.9)

polygon(x = c(occurpred_hawk.1$Age_Scaled[MM],rev(occurpred_hawk.1$Age_Scaled[MM])), y = c(occurpred_hawk.1$lci[MM],rev(occurpred_hawk.1$uci[MM])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hawk.1$Age_Scaled[MM],y = occurpred_hawk.1$fit[MM],lwd = 2,col = 'grey30',lty = 1)

polygon(x = c(occurpred_hawk.1$Age_Scaled[MMM],rev(occurpred_hawk.1$Age_Scaled[MMM])), y = c(occurpred_hawk.1$lci[MMM],rev(occurpred_hawk.1$uci[MMM])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=occurpred_hawk.1$Age_Scaled[MMM],y = occurpred_hawk.1$fit[MMM],lwd = 2,col = 'grey30',lty = 2)

legend('topright',legend = c('Winter', "Spring"), lty = c(1,2), col = 'grey30',pt.cex = 1)


#END----
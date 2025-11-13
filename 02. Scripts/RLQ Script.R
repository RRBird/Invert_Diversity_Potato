#Libraries----

library("ade4")
library("dplyr")
library("tidyverse")
library("vegan")

Q_matrix$Size <- factor(Q_matrix$Size, 
                        levels = c("0-2.5", "2.5-5", "5-10", 
                                   ">10","No_Size", "Unknown"),
                        ordered = FALSE) #to make it work properly
str(Q_matrix)

R_matrix$Position <- factor(R_matrix$Position, 
                            levels = c("Inner", "Outer"),
                            ordered = FALSE) #to make it work properly
str(R_matrix)

#Species responses to environmental gradients----

coa1 <- dudi.coa(L_matrix, scannf = F)
cca1 <- pcaiv(coa1, R_matrix, scannf = F)

#percentage of variation in species composition explained by enviro
100 * sum(cca1$eig) / sum(coa1$eig) 

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
s.label(cca1$c1, clabel = 0)
par(mar = c(0.1, 0.1, 0.1, 0.1))
s.arrow(cca1$cor[-1,], add.plot=TRUE)

#RLQ Analysis---- 

pca.traits <- dudi.hillsmith(Q_matrix, row.w = coa1$cw, scannf = FALSE)
pca.env <- dudi.hillsmith(R_matrix, row.w = coa1$lw, scannf = FALSE)

rlq1 <- rlq(pca.env, coa1, pca.traits, scannf = FALSE)
summary(rlq1)

dev.new(height=10,width=15,dpi=80,pointsize=14,noRStudioGD = T)
plot(rlq1)

## Percentage of co-Inertia for each axis
100*rlq1$eig/sum(rlq1$eig)

#To interpret the results, correlations can be computed:
## weighted correlations axes / env.
t(pca.env$tab)%*%(diag(pca.env$lw))%*%as.matrix(rlq1$mR)

## weighted correlations axes / traits.
t(pca.traits$tab)%*%(diag(pca.traits$lw))%*%as.matrix(rlq1$mQ)

## correlations traits / env.
rlq1$tab

#bplot
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
s.arrow(rlq1$c1, xlim=c(-1,1), boxes = FALSE)
s.label(rlq1$li, add.plot=T, clab=1.5)

#Another plot
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
s.label(rlq1$lQ, clabel = 0)
par(mar = c(0.1, 0.1, 0.1, 0.1))
pointLabel(rlq1$lQ,row.names(rlq1$lQ), cex=0.7) #pointLabel didn't work so need an alternative to add names to this one

#Classifying scores to obtain functional groups
hc2 <- hclust(dist(rlq1$lQ), method = "ward.D")
dev.new(height=20,width=40,dpi=80,pointsize=14,noRStudioGD = T)
plot(hc2)

#Calinsky-Harabasz criteria to find best partition 
ntest <- 6
res <- rep(0,ntest - 1)
for (i in 2:ntest){
  fac <- cutree(hc2, k = i)
  res[i-1] <- calinski(tab=rlq1$lQ, fac = fac)[1]
}
#calinski ISN"T WORKING SO NEED TO FIGURE THAT OUT----


dev.new(height=5,width=5,dpi=80,pointsize=14,noRStudioGD = T)
par(mfrow=c(1,2))
plot(2:ntest, res, type='b', pch=20, xlab="Number of groups", ylab = "C-H index")
plot(3:ntest, diff(res), type='b', pch=20, xlab="Number of groups", ylab = "Diff in C-H index")






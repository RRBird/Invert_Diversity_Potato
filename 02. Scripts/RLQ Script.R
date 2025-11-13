#Libraries----

library("ade4")
library("dplyr")
library("tidyverse")
install.packages("vegan")

#Species responses to environmental gradients----

coa1 <- dudi.coa(L_matrix, scannf = F)
cca1 <- pcaiv(coa1, R_matrix, scannf = F)

#percentage of variation in species composition explained by enviro
100 * sum(cca1$eig) / sum(coa1$eig) 

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
s.label(cca1$c1, clabel = 0)
par(mar = c(0.1, 0.1, 0.1, 0.1))
pointLabel(cca1$c1,row.names(cca1$c1), cex=0.7)
s.arrow(cca1$cor[-1,], add.plot=TRUE)



#correspondence analysis on matrix L
coa1 <- dudi.coa(species, scannf = F)
head(coa1)



#Step 1: Correspondence Analysis on L table
dudi_L <- dudi.coa(L_matrix, scannf = FALSE, nf = 2)

#Step 2: Hill-Smith analysis on Q table (handles mixed trait types)
dudi_Q <- dudi.hillsmith(Q_matrix, scannf = FALSE, nf = 2, row.w = dudi_L$cw)

#Step 3: PCA on R table
dudi_R <- dudi.pca(R_matrix, scannf = FALSE, nf = 2, row.w = dudi_L$lw)

#Step 4: RLQ analysis
rlq_result <- rlq(dudi_R, dudi_L, dudi_Q, scannf = FALSE, nf = 2)

pca.traits <- dudi.pca(traits,row.w = coa1$cw,scannf=F)







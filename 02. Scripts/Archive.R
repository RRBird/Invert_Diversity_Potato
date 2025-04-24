

##Beta Diversity

#NOT COMPLETED

#Bray-Curtis - abundance data

###Araneae

Araneae_matrix3 <- as.data.frame.matrix(Araneae_matrix2)
head(Araneae_matrix3);dim(Araneae_matrix3)

Araneae_Beta <- vegdist(Araneae_matrix3, method = 'bray')


#TRIAL BELOW TO LOOK AT BETA - IS THIS ENOUGH??

bray_matrix <- as.matrix(Araneae_Beta)

# Create a heatmap
heatmap(bray_matrix, 
        main = "Bray-Curtis Dissimilarity Matrix", 
        col = heat.colors(256),  # Color palette
        scale = "none")  # No scaling to preserve raw distances

mds <- cmdscale(Araneae_Beta, k = 2)  # k = 2 for 2D MDS plot

# Plot the results
plot(mds, 
     xlab = "Dimension 1", 
     ylab = "Dimension 2", 
     main = "MDS of Bray-Curtis Dissimilarity")


# NMDS (Non-metric Multidimensional Scaling)
nmds <- metaMDS(Araneae_matrix3, distance = "bray", k = 2)
#Warning message:
##In metaMDS(Araneae_matrix3, distance = "bray", k = 2): stress is (nearly) zero: you may have insufficient data

# Plot the NMDS
plot(nmds, main = "NMDS of Bray-Curtis Dissimilarity")

# Hierarchical clustering using Bray-Curtis dissimilarity
hclust_result <- hclust(Araneae_Beta)

# Plot the dendrogram
plot(hclust_result, main = "Hierarchical Clustering of Sites", xlab = "Sites", sub = "")

###Hemieptera

###Coleoptera
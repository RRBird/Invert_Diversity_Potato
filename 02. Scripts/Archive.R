
#Old Diversity calculation (diversity by orders instead of functional groups)


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
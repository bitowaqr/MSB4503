HPV = read.csv(file.choose())
HPV_meta = read.csv(file.choose())

HPV_ID = HPV[,1]
HPV = HPV[,-1]
HPV = HPV[,-length(HPV)]
HPV[1:5,1:5]

HPV = t(apply(HPV,1,scale))

# Which factor is the major driving source of clustering when considering the clusters of samples?

# loop trhugh explanatory factors/labels
i = 1

names(HPV_meta)
i = i + 1
sample = t(HPV[,])
groupCodes <- HPV_meta[,i]
rownames(sample) <- make.unique(as.character(groupCodes))
colorCodes <- rainbow(nlevels(groupCodes))
distSamples <- dist(sample)
hc <- hclust(distSamples)
dend <- as.dendrogram(hc)
library(dendextend)
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

# PLOT IT
dend <- set(dend, "labels_cex", 0.5)
par(mar = c(3,3,3,7))
plot(dend, horiz = TRUE,main=names(HPV_meta)[i])

## A sub tree - so we can see better what we got:
par(cex = 1)
plot(dend[[2]], horiz = TRUE)



library(dendextend)


HPV = read.csv(file.choose()) # raw data
HPV_meta = read.csv(file.choose()) # meta data

HPV_ID = HPV[,1]
HPV = HPV[,-1] # remove ID
HPV = HPV[,-length(HPV)] # the last column is the mean, remove
# HPV[1:5,1:5]  # have a look

HPV = t(apply(HPV,1,scale)) # scale per row

# Q1: Which factor is the major driving source of clustering 
# when considering the clusters of samples?

# MANUALLY loop trhugh explanatory factors/labels

plot.my.dend = function(i = 1, mat = HPV_donadj,...){
  names(HPV_meta)
  sample = t(mat[,])
  groupCodes <- HPV_meta[,i]
  rownames(sample) <- make.unique(as.character(groupCodes))
  colorCodes <- rainbow(nlevels(groupCodes))
  distSamples <- dist(sample)
  hc <- hclust(distSamples)
  dend <- as.dendrogram(hc)
  
  labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
  
  # PLOT IT
  dend <- set(dend, "labels_cex", 0.5)
  # par(mar = c(3,3,3,7))
  plot(dend, horiz = TRUE,...)
}
par(mfrow=c(1,1))
man = 6 # look at i = 2,...6
plot.my.dend(i = man, mat = HPV, main = names(HPV_meta)[man])

# Q2: To try to eliminate the donor effect, 
# for each donor subtract their pre-vaccination 
# media sample from the other samples.  
# Does this eliminate the donor-effect and give you 
# the expected three clusters in the remaining data?


donors = unique(HPV_meta$X.title.)
HPV_donadj = HPV
total_len = length(unique(HPV_meta$X.title.))
counter = 0
for(d in donors){
  counter = counter +1
  cat("\n", round(counter/total_len,4)*100,"%",sep="")
  index_donor = HPV_meta$X.title. == d
  mean_temp = apply(HPV_donadj[,index_donor],1,mean)
  HPV_donadj[,index_donor] = HPV_donadj[,index_donor] - mean_temp
}


# plot: whole.class labels
par(mfrow=c(1,2))
plot.my.dend(6,HPV_donadj,,axes=F,main="adjusted")
plot.my.dend(6,HPV,axes=F,main="raw")[[1]]



# Q3: Find sets of genes which have the same pattern of gene expression, 
# regardless of the magnitude.
HPV_cor = cor(t(HPV))
dist_cor = dist(HPV_cor) # Wow... this takes some time!
hclust_cor = hclust(dist_cor)
dend_cor <- as.dendrogram(hclust_cor)
# PLOT IT
dend_cor <- set(dend_cor, "labels_cex", 0.5)
# par(mar = c(3,3,3,7))
plot(dend_cor, horiz = TRUE,main="Correlation clustering")


# Q4: Find sets of genes which have the same pattern of gene expression, regardless of the magnitude.


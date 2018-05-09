
# PCA
# SOMETHING IS GOING SEVERELY WRONG!


# load the liver-data
data = read.csv(file.choose())

# create data sets for dose and hour labels
# for uploading to website
data_dose = data[-2,]
write.csv(data_dose,"data_dose.csv",row.names = F)
data_hours = data[-1,]
write.csv(data_hours,"data_hours.csv",row.names = F)

#### data manipulatins...
str(data_dose)
t.data_dose = data_dose[-1,-1]
colnames(t.data_dose) = t.data_dose[1,-1]
rownames(t.data_dose) = t.data_dose[-1,1]
t.data_dose = as.matrix(t(t.data_dose))
t.data_dose = apply(t.data_dose,1,as.numeric)
dim(t.data_dose)
t.data_dose = apply(t.data_dose,1,as.numeric)

r.data_dose = data_dose
names(r.data_dose) = colnames(data_dose)
rownames(r.data_dose) = rownames(data_dose)
rownames(r.data_dose) = r.data_dose[,1]
r.data_dose = r.data_dose[,-1]
r.data_dose[-1,] = apply(r.data_dose[-1,],2,function(x)as.numeric(as.character(x)))
data_dose.dose = r.data_dose[1,]
r.data_dose = r.data_dose[-1,]
r.data_dose = apply(r.data_dose,2,function(x)as.numeric(as.character(x)))
rownames(r.data_dose) = data_dose[-1,1]

# filter step 1: removing 40% of variables with lowest IQR
feature_IQR = apply(r.data_dose,1,IQR)
hist(feature_IQR)
value.40 = quantile(feature_IQR[order(feature_IQR,decreasing = T )],0.4)
select.rows = which(feature_IQR >= value.40)
data_dose_filter.60 = r.data_dose[select.rows,]

# plot what happend
par(mfrow=c(1,2))
d <- density(r.data_dose) 
plot(d) 
d.filtered <- density(data_dose_filter.60) 
plot(d.filtered) 

# mean centering per feature
data_dose_filter.60.centered = apply(data_dose_filter.60,1,function(x)scale(x,center = T,scale = F))
data_dose_filter.60.centered = t(data_dose_filter.60.centered)
dim(data_dose_filter.60.centered)
t.data_dose_filter.60.centered = t(data_dose_filter.60.centered)

# plot what happened
par(mfrow=c(1,2))
d <- density(data_dose_filter.60) 
plot(d) 
d.filtered <- density(data_dose_filter.60.centered ) 
plot(d.filtered) 

# boxplotting the first 20 unfiltered / filtered featrues 
par(mfrow=c(1,2))
group = as.factor(1:20)
boxplot.df = r.data_dose[1:20,]
boxplot.df = data.frame(group=group,boxplot.df)
boxplot.df = reshape2::melt(boxplot.df)
boxplot(value~group,boxplot.df,horizontal = T)

boxplot.df = data_dose_filter.60.centered[1:20,]
boxplot.df = data.frame(group=group,boxplot.df)
boxplot.df = reshape2::melt(boxplot.df)
boxplot(value~group,boxplot.df,horizontal = T)

# PCA on filtered data set [something is going wrong....]
pca.data_dose = prcomp(t.data_dose_filter.60.centered,center = T,scale. = T)
plot(pca.data_dose)
xxx = (pca.data_dose$sdev/sum(pca.data_dose$sdev))
plot(xxx*100,type="l")

pca.data_dose$rotation[,1]
PCAs = pca.data_dose$rotation
PCAs = data.frame(PCAs)
PCAs = data.frame(dose = t(data_dose.dose) , PCAs)
plot(PCAs[,2:4],col=PCAs$dose)
plot(pca.data_dose$rotation)

# PCA on UN-filtered data set [something is going wrong....]
c.r.data_dose = apply(r.data_dose,1,function(x)scale(x,center = T,scale = F))
c.r.data_dose = t(c.r.data_dose)
apply(c.r.data_dose,2,mean)
pca.data_dose = prcomp(c.r.data_dose)
plot(pca.data_dose)
xxx = pca.data_dose$sdev/sum(pca.data_dose$sdev)
plot(xxx*100,type="l")
plot(pca.data_dose$rotation[,1:3])
pca.data_dose$rotation[,1]
PCAs = pca.data_dose$rotation
PCAs = data.frame(PCAs)
PCAs = data.frame(dose = t(data_dose.dose) , PCAs)
plot(PCAs[,2:4],col=PCAs$dose)
plot(pca.data_dose$rotation)


#### USING ANOTHER PCA FUNCTION FROM BIOCINDUCTOR
source("https://bioconductor.org/biocLite.R")
biocLite("ropls")
library(ropls)
?ropls

# opls 
tt.data_dose_filter.60.centered = t(t.data_dose_filter.60.centered)
r.t.dose.data = t(r.data_dose)

## raw data opls
pca.data = opls(r.t.dose.data)
summary(pca.data)
dose.col = t(data_dose.dose)
dose.col = as.factor(dose.col)
plot(pca.data, parAsColFcVn = dose.col)

# again, raw data using generic r pca function.... whats going on here?
pca2 = prcomp(r.t.dose.data)
xxx = pca2$sdev/sum(pca2$sdev)
plot(xxx*100,type="l")



##### frustrated end
library(csv)
library(psych)
library(FactoMineR)
library(PoiClaClu)
library(FactoMineR)
library(questionr)
library(lars)
library(PoiClaClu)
library(plot3D)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(extrafont)
library(RColorBrewer)
library(ggplot2)
library(ggcorrplot)
library(factoextra)
library(ggpubr)
library(preprocessCore)
library(Rtsne)
library(caret)
library(umap)

# UNSUPERVISED ANALYSIS

# Is needed transform to factor the variables to predict example: Diagnostic and sample ID
data$Study_Name=as.factor(data$Study_Name)

# if is needed apply tranformation to you data, select the columns of biomarkers, example

elisa = data[,13:100] #columns of biomarkers selected

dataX=log2(elisa+1) #biomarkers trasformed to log2

#The following function will perform a pdf with PCA and clustering of your data
# You could define inside the funcion the threeshold of the clustering

coul=rainbow(4)# choose the colors based on your number of variables to predict
pdf("descriptive statistics.pdf")
acp=PCA(dataX)
plot(as.numeric(acp$ind$coord[,1]),as.numeric(acp$ind$coord[,2]),col=coul[data$Study_Name],ylim=c(-4,4),xlim=c(-2,4),main="Representation of sample Elisa assay", xlab=paste("Dim1: ",round(acp$eig[1,2],2),"%",sep=""),ylab=paste("Dim2: ",round(acp$eig[2,2],2),"%",sep=""),pch=19)
abline(v=0,h=0,lty=3)
legend("topright",legend=levels(data$Study_Name),pch=19,col=coul,bty='n')

distance2=dist(dataX,method="euclidean")
hier2=hclust(distance2,method="ward.D2") 
ColorDendrogram(hier2,y=coul[as.factor(data$Relapse)],labels=paste(as.character(data$Relapse)),
                main="Clustering",xlab="sample",sub="",branchlength =70)
ColorDendrogram(hier2,y=coul[as.factor(data$Relapse)],labels=row.names(data),
                main="Clustering",xlab="sample",sub="",branchlength = 70)
ColorDendrogram(hier2,y=coul[as.factor(data$Relapse)],labels=paste(as.character(data$Disease.activity.since.last.return)),
                main="Clustering",xlab="sample",sub="",branchlength = 70)

dev.off()

# The following function will perform a PDF with PCA plus ellipses with 95% of CI of your data
# Is more useful to visualize outliers

pdf("descriptive statistics elipses.pdf")
fviz_pca_ind(acp, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = data$Study_Name, 
             col.ind = "black", 
             palette = "ggplot2", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Study") +
  geom_text(aes(label = data$Firalis_Sample_ID), size = 3) +
  ggtitle("2D PCA-plot from assay") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


## t-SNE dimension reduction
#mat is the matrix
data0=read.csv(choose.files(),header=T)

bm<-data0[,12:1048]
dataX=log2(bm+1)

set.seed(80) # Set a seed if you want reproducible results
tsne_out <- Rtsne(as.matrix(bm), dims = 3, perplexity = 1) # Run TSNE

#image(t(as.matrix(dist(tsne_out$Y))))
#Show the objects in the 2D tsne representation

plot(tsne_out$Y,col=as.factor(data0$Relapse.status),main="t-SNE:", xlab=paste("Dim1: ",round(tsne_out$Y[1,2],2),"%",sep=""),ylab=paste("Dim2: ",round(tsne_out$Y[2,2],2),"%",sep=""), pch=19)

# create the legend for the study groups
legend("topleft",
       legend=unique(data0$Relapse.status),
       fill =palette("default"),
       border="black",box.col="black")

#Function plot umap

plotumap <- function(x, labels,
                     main="A UMAP visualization of the dataset",
                     colors=c("red", "blue", "#17becf", "green"),
                     pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
                     cex.main=1, cex.legend=0.85) {
  
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  } 
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u = unique(labels)
  legend.pos = "topleft"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomleft"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

#coordinates used to visualize the dataset
plotumap(Elisa.umap, Elisa.labels, main=" Text ")

#configuration objects
umap.defaults

custom.config = umap.defaults
custom.config$random_state = 200

#changed settings
bm.umap.config = umap(biop_relap.data, config=custom.config)
plotumap(bm.umap.config, biop_relap.labels,
         main=" Main text(different seed = 200)")

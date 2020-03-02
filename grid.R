library(geomorph)
data(plethodon) 
plethodon
Y.gpa<-gpagen(plethodon$land)    #GPA-alignment
Y.gpa

gp <- as.factor(paste(plethodon$species, plethodon$site))
col.gp <- rainbow(length(levels(gp)))
names(col.gp) <- levels(gp)
col.gp <- col.gp[match(gp, names(col.gp))]
PCA <- plotTangentSpace(Y.gpa$coords, verbose = T)
PCA$pc.summary$importance
xlab <- paste('Principal Component 1', '(', round(PCA$pc.summary$importance[2,1]*100, 1), '%)', sep='')
ylab <- paste('Principal Component 2', '(', round(PCA$pc.summary$importance[2,2]*100, 1), '%)', sep='')
mat <- matrix(c(4,5,0,1,1,2,1,1,3), 3)
layout(mat, widths=c(1,1,1), heights=c(1,1,0.6))
par(mar=c(4, 4, 1, 1)) # sets the margins
plot(PCA$pc.scores[,1], PCA$pc.scores[,2], pch=21, cex=2, bg=col.gp, xlab=xlab, ylab=ylab, asp=T)
legend(-0.09, 0.07, legend= unique(gp), pch=19,  col=unique(col.gp))
ref <- mshape(Y.gpa$coords)
plotRefToTarget(ref,PCA$pc.shapes$PC1min,outline=plethodon$outline)

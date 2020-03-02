install.packages("MASS")
install.packages("Morpho")



library(colorspace)  
library(cluster)     
library(DiscriMiner) 
library(ellipse)     
library(geomorph)    
library(ggplot2)     
library(mclust)      
library(NbClust)     
library(shapes)      
library(vegan)  
library(plotly)

setwd("C:/Users/grati/OneDrive/LEPAV/Cic's/CIC 2019/Morpho2.0")
source("colLab.R")

setwd("C:/Users/grati/OneDrive/LEPAV/Cic's/CIC 2019/Morpho2.0/dados para analise/intraespecifica/Tri x Tri")

#DEFINE OS GRUPOS DE ESTUDO PARA FAZER A COMPARAÇÃO DAS ANÁLISES E LEITURA DOS DADOS

list_tt<-list.files(pattern="nts")
list_att <- list.files(pattern = "att")
list_btt <- list.files(pattern = "btt")
list_ctt <- list.files(pattern = "ctt")

t_t<-readmulti.nts(list_tt)
att <- readmulti.nts(list_att)
btt <- readmulti.nts(list_btt)
ctt <- readmulti.nts(list_ctt)

#CRIAR UM DIRETÓRIO PARA COLOCAR OS REULTADOS EM UM LUGAR DIFERENTE

write.table(t_t, file = "t_t")
write.table(att,file = "att.txt")
write.table(btt,file = "btt.txt")
write.table(ctt,file = "ctt.txt")

#GENERALIZANDO A ANALISE DOS PROCRUSTES ANALYSIS (GPA) E OS SHAPES DE CONSENSO

groups_intra <- list(t_t=t_t, att=att, btt=btt, ctt=ctt)

alg_intra <- lapply (groups_intra, gpagen)
alg_t_t<-alg_intra$t_t$coords#; writeland.tps(t_t,"alg_t_t")
alg_att<-alg_intra$att$coords#; writeland.tps(alg_att,"alg_região1_t")
alg_btt<-alg_intra$btt$coords#; writeland.tps(alg_btt,"alg_região2_t")
alg_ctt<-alg_intra$ctt$coords#; writeland.tps(alg_ctt,"alg_região2_t")

cons_namess<-list(alg_t_t=alg_t_t, alg_att=alg_att, alg_btt=alg_btt, alg_ctt=alg_ctt)

cons_shapes<-lapply(cons_namess,mshape)
cons_tabanus<-cons_shapes$alg_t_t
cons_att<-cons_shapes$alg_att
cons_btt<-cons_shapes$alg_btt
cons_ctt<-cons_shapes$alg_ctt

plotOutliers(alg_t_t)
outliers <- find.outliers(alg_t_t)
outliers$dist.sort


######VISUALIZAÇÃO DA FORMA
#FORMA MÉDIA 
ref =mshape(shape.sym)


list

## 5. Principal Component Analysis (PCA)
factors_t_occ <-sub(".nts", "", list_tt)
for (i in (rep(c(0,1,2,3,4,5,6,7,8,9),5))) {
  factors_t_occ <- sub(i, "", factors_t_occ)
}
factors_t_occ <- factor(factors_t_occ)
summary(factors_t_occ)

pca_t_occ <- plotTangentSpace(alg_intra$coords, label=NULL, verbose =T, groups = factors_t_occ, warpgrids=F) ## Perform PCA of Diachlorus curvipes
pca_t_occ$pc.summary$importance##saber a parte cumulativa 

label_dc <- list(expression('Região 1'), expression('Região 2'), ('Região 3'))   
a <- cbind(pca_t_occ$pc.scores, factors_t_occ)
a <- data.frame(a)
a$factors_t_occ <- as.factor(a$factors_t_occ)
df_ell <- data.frame()
for(g in levels(a$factors_t_occ)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$factors_t_occ==g,], ellipse(cor(PC1, PC2), 
                                                                                     scale=c(sd(PC1),sd(PC2)), 
                                                                                     centre=c(mean(PC1),mean(PC2))))),factors_t_occ=g))
}
p_curvi <- ggplot(data=a, aes(x=PC1, y=PC2,colour=factors_t_occ)) + geom_point(size=2, shape= 18) 
p_curvi <- p_curvi + geom_path(data=df_ell, aes(x=x, y=y,colour=factors_t_occ), size=0.5, linetype=2) 
p_curvi <- p_curvi + scale_color_manual(name ="Areas", labels= label_dc,values=c("#1b9e77", "darkorange4","#7570b3")) 
p_curvi <- p_curvi + scale_shape_manual(name ="Areas", labels= label_dc, values=c(15, 16, 17,18)) 


p_curvi <- p_curvi + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) 
p_curvi + theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size=16))
ggplotly(p_curvi)


######################################################################################################################################################################################
## 6. Hierarchical clustering
datavalores_t_occ <- data.frame(pca_t_occ$pc.scores)            
distvalores_t_occ <- dist(datavalores_t_occ)                                            
clust_valores_t_occ <- hclust(distvalores_t_occ)                  
labelColors <- c("#1b9e77", "royalblue4", "#7570b3" )     
clusMember <- rep(1,length(rownames(datavalores_t_occ)))
clusMember[grep("att",rownames(datavalores_t_occ))]<-1 
clusMember[grep("btt",rownames(datavalores_t_occ))]<-2 
clusMember[grep("ctt",rownames(datavalores_t_occ))]<-3
names(clusMember) <- rownames(datavalores_t_occ)       
clusDendro_t_occ <- as.dendrogram(as.hclust(clust_valores_t_occ)) 
clusDendro_t_occ <- dendrapply(clusDendro_t_occ, colLab)  

par(mar=c(0.5, 4.5, 1, 1))
par(oma= c(0,0,1,0))
op = par(bg = "gray90")
par(lwd=3)
plot(clusDendro_t_occ,horiz=F,axes=T, ylab= "Height", cex.axis=1.3,cex.lab=1.7)
par(lwd=1)
legend("topright", pch= 21, pt.bg=c("#1b9e77", "royalblue4"), 
       legend=expression("Região 1", "Região 2", "Região 3"), pt.cex = 1)
title(main= (expression(paste( " ", italic("Tabanus triangulum")))), outer = F)

## 7. Define the number of clusters using model-based clustering, "Gap" statistic and Calinski-Harabasz criterion
### Model-based clustering

new_pcs_t_occ <- subset (pca_t_occ$pc.scores, select = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))         
best_model_t_occ <- Mclust(new_pcs_t_occ, G=1:3)                                                                     
summary(best_model_t_occ)                                                                                                 

plot(best_model_t_occ, what = c("BIC"), dimens = c(1,2))
plot(best_model_t_occ, what = c("classification"), dimens = c(1,2))

##"Gap" statistic

gap_t_occ <- clusGap(new_pcs_t_occ, kmeans, 2, B = 1000, verbose = interactive())                                    ## Perform the gap statistic for D. curvipes
gap_t_occ                                                                                                                 ## See the results
plot(gap_t_occ)

### Calinski-Harabasz criterion

fit_t_occ <- cascadeKM(scale(new_pcs_t_occ, center = TRUE,  scale = TRUE), 1, 2, iter = 1000)                          ## Perform analysis using Calinski-Harabasz criterion for D. curvipes
fit_t_occ$results                                                                                                           ## See the results
            
plot(fit_t_occ)

## 8. Discriminant analysis (linear and quadratic)
## discriminant analysis of D. curvipes

dis_lin_t_occ <- linDA(new_pcs_t_occ, factors_t_occ, validation="crossval")
percentage_confusion_dislineal_t_occ <- NULL
for (i in 1:dim(dis_lin_t_occ$confusion)[1]) {
  row_new<- (dis_lin_t_occ$confusion[i, ]*100)/sum(dis_lin_t_occ$confusion[i, ])
  percentage_confusion_dislineal_t_occ <- rbind(percentage_confusion_dislineal_t_occ, row_new)
}
row.names(percentage_confusion_dislineal_t_occ) <- dimnames(percentage_confusion_dislineal_t_occ)[[2]]
write.csv(percentage_confusion_dislineal_t_occ, file="dislineal_t_occ_confusion.csv")
dis_lin_t_occ$error

dis_qua_t_occ <- quaDA(new_pcs_t_occ, factors_t_occ, validation="crossval")
percentage_confusion_disqua_t_occ <- NULL
for (i in 1:dim(dis_qua_t_occ$confusion)[1]) {
  row_new<- (dis_qua_t_occ$confusion[i, ]*100)/sum(dis_qua_t_occ$confusion[i, ])
  percentage_confusion_disqua_t_occ <- rbind(percentage_confusion_disqua_t_occ, row_new)
}
row.names(percentage_confusion_disqua_t_occ) <- dimnames(percentage_confusion_disqua_t_occ)[[2]]
write.csv(percentage_confusion_disqua_t_occ, file="discua_t_occ_confusion.csv")
dis_qua_t_occ$error

######################################################################################################################################################################################
## 9. Riemmanian distances among species and Goodall's F test 

cons_shapes

names_t_occ <- combn(names(cons_shapes), m=2)
comb_t_occ <- combn(cons_shapes, m=2)
r_t_occ <- NULL
nam <- NULL
for (i in 1:dim(comb_t_occ)[2]) {
  n <- paste(names_t_occ[,i][1], "to", names_t_occ[,i][2])  
  r <- riemdist(as.matrix(comb_t_occ[,i][1][[1]]), as.matrix(comb_t_occ[,i][2][[1]]))
  r_t_occ <- rbind(r_t_occ, r)
  nam <- c(nam, n)
}
rownames(r_t_occ) <- nam
colnames(r_t_occ) <- c("Riemman")
write.csv(r_t_occ, "Riemmanian_distances_occ.csv")

alg_names_t_occ <- combn(names(alg_intra), m=2)
comb_alg_t_occ <- combn(alg_intra, m=2)
name <- NULL
t_goodall_t_occ <- NULL
for (i in 1:dim(comb_alg_t_occ)[2]) {
  n <- paste(alg_names_t_occ[,i][1], "and", alg_names_t_occ[,i][2]) 
  t <- testmeanshapes(comb_alg_t_occ[,i][1][1][[1]][[1]], comb_alg_t_occ[,i][2][1][[1]][[1]], resamples = 1000, replace = FALSE, scale= TRUE)
  t_value <- c(t$G, t$G.pvalue)
  t_goodall_t_occ <- rbind(t_goodall_t_occ, t_value)
  name <- c(name, n)
}
rownames(t_goodall_t_occ) <- name
colnames(t_goodall_t_occ) <- c("G", "G_pvalue")
write.csv(t_goodall_t_occ, "Goodall_t_occ.csv")
t_goodall_t_occ


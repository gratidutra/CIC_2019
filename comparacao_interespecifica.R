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

setwd("C:/Users/grati/Dropbox/IC_Grati/Morpho2.0")
source("colLab.R")



library(plotly)
setwd("C:/Users/grati/Dropbox/IC_Grati/Morpho2.0/dados para analise/tabanus")
list_names<-list.files(pattern = "nts")
list_tt <- list.files(pattern = "tt")      
list_tto <- list.files(pattern = "tto")
list_ato <- list.files(pattern = "ato")

tabanus<-readmulti.nts(list_names)
tt <- readmulti.nts(list_tt)
tto<-readmulti.nts(list_tto)
ato<-readmulti.nts(list_ato)

write.table(tabanus,file = "tabanus.txt")
write.table(tt,file = "tt.txt")
write.table(tto,file = "att.txt")
write.table(ato,file = "ato.txt")

groups_inter <- list(tabanus=tabanus, tto=tto, ato=ato, tt=tt)
alg_inter<-lapply(groups_inter,gpagen)
alg_tabanus<-alg_inter$tabanus$coords; writeland.tps(tabanus,"tabanus")
alg_tto<-alg_inter$tto$coords; writeland.tps(alg_tto,"alg_TO_t_occidentalis")
alg_tt<-alg_inter$tt$coords; writeland.tps(alg_tt,"alg_triangulum")
alg_ato<-alg_inter$ato$coords; writeland.tps(alg_ato,"alg_RS_t_occidentalis")



cons_names<-list(alg_tabanus=alg_tabanus, alg_tto=alg_tto, alg_ato=alg_ato, alg_tt=alg_tt)
cons_shapes<-lapply(cons_names,mshape)
cons_tabanus<-mshape(cons_shapes$alg_tabanus);writeland.tps(alg_tto,"cons_tabanus")
cons_ato<-mshape(cons_shapes$alg_ato)
cons_tto<-mshape(cons_shapes$alg_tto)
cons_tt<-mshape(cons_shapes$alg_tt)


factors_tabanus <-sub(".nts", "", list_names)
for (i in (rep(c(0,1,2,3,4,5,6,7,8,9),5))) {
  factors_tabanus <-sub(i, "", factors_tabanus)
}
factors_tabanus
factors_tabanus <- sub("ato","to", factors_tabanus) ; factors_tabanus <- sub("att", "tt", factors_tabanus)
factors_tabanus <- sub("btt", "tt", factors_tabanus) ; factors_tabanus <- sub("ctt", "tt", factors_tabanus)
factors_tabanus <- sub("tto", "to", factors_tabanus);
factors_tabanus <- factor(factors_tabanus)
factors_tabanus
pca_tabanus <- plotTangentSpace(alg_tabanus, label=NULL, verbose =T, groups = factors_tabanus, warpgrids=T) 
pca=pca_tabanus$pc.summary$importance##saber a parte cumulativa 
write.csv(pca, "pca_cumulative.csv")
summary(factors_tabanus)

#NO BOXPLOT

label <- list(expression(italic('T. occidentalis RS')), expression(italic('T. triangulum')),
              expression(italic("T. occidentalis TO"))) ## Label to the plot 
a <- cbind(pca_tabanus$pc.scores, factors_tabanus)                                                                      ## Attach the factor for the names to the data 
a <- data.frame(a) ## Convert class of 'a' to data frame
a$factors_tabanus <- as.factor(a$factors_tabanus)

### Define the ellipses to the groups (95%)

df_ell <- data.frame()
for(g in levels(a$factors_tabanus)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(a[a$factors_tabanus==g,], 
            ellipse(cor(PC1, PC2),scale=c(sd(PC1),sd(PC2)), 
            centre=c(mean(PC1),mean(PC2))))),factors_tabanus=g))
}

p <- ggplot(data=a, aes(x=PC1, y=PC2,colour=factors_tabanus, shape=factors_tabanus)) + geom_point(size=3) 

p <- p + geom_path(data=df_ell, aes(x=x, y=y,colour=factors_tabanus), size=0.5, linetype=2) 

p <- p + scale_colour_brewer(name= "Species", labels= label, palette = "Set1") +
    scale_shape_discrete(name ="Species", labels= label) 

p <- p + theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
      theme(legend.text = element_text(size = 14)) + theme(legend.title = element_text(size=16))

p
ggplotly(p)
#####################################################################################################################################################################################################

datavalores <- data.frame(pca_tabanus$pc.scores)                                
distvalores <- dist(datavalores)                                                   
clust_valores <- hclust(distvalores)                                              
labelColors <- c("#377eb8","#e43a1c","#984ea3") 
clusMember <- rep(1,length(rownames(datavalores)))
clusMember[grep("ato",rownames(datavalores))]<-1                                    
clusMember[grep("tt",rownames(datavalores))]<-2
clusMember[grep("tto",rownames(datavalores))]<-3

names(clusMember) <- rownames(datavalores)                                         
clusDendro <- as.dendrogram(as.hclust(clust_valores))                              
clusDendro <- dendrapply(clusDendro, colLab)                                      
plot(clusDendro)
par(mar=c(0.5, 4.5, 1, 1))
par(oma= c(0,0,0,0))
op = par(bg = "gray90")
par(lwd=2)
plot(clusDendro,horiz=F,axes=T, ylab= "Height", cex.axis=1.3,cex.lab=1.7)
par(lwd=1)
 legend("topright", pch= 21, pt.bg=c("#377eb8","#e43a1c","#984ea3"), 
   legend = expression(italic("T. occidentalis RS"),italic("T. triangulum"),
            italic("T. occidentalis TO")), 
   pt.cex = 1.5, cex=0.55)

 
#############################################################################################################################
## 7. Define the number of clusters using model-based clustering, "Gap" statistic and Calinski-Harabasz criterion

new_pcs_tabanus <- subset (pca_tabanus$pc.scores, select = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8))      ## Make a sub-dataset using the first eight components
new_pcs_tabanus
best_model_tabanus <- Mclust(new_pcs_tabanus, G=1:4)                                                       ## Perform model clustering
summary(best_model_tabanus)                                                                                   ## See the results of the model-based clustering


plot(best_model_tabanus, what = c("BIC"), dimens = c(1,2))
plot(best_model_tabanus, what = c("classification"), dimens = c(1,2))

##################################################################################################################

gap_tabanus <- clusGap(new_pcs_tabanus, kmeans, 4, B = 1000, verbose = interactive())                     ## Perform the gap statistic
gap_tabanus                                                                                                 ## See the results

plot(gap_tabanus)

#################################################################################################################

fit <- cascadeKM(scale(new_pcs_tabanus, center = TRUE,  scale = TRUE), 1, 4, iter = 1000)                    ## Perform analysis using Calinski-Harabasz criterion
fit$results                                                                                                     ## See the results

plot(fit, sortg = TRUE, grpmts.plot = TRUE)

######################################################################################################################################################################################
## 8. Discriminant analysis (linear and quadratic)

dis_lin_tabanus <- linDA(new_pcs_tabanus, factors_tabanus, validation="crossval")
percentage_confusion_dislineal <- NULL
for (i in 1:dim(dis_lin_tabanus$confusion)[1]) {
  row_new<- (dis_lin_tabanus$confusion[i, ]*100)/sum(dis_lin_tabanus$confusion[i, ])
  percentage_confusion_dislineal <- rbind(percentage_confusion_dislineal, row_new)
}
row.names(percentage_confusion_dislineal) <- dimnames(percentage_confusion_dislineal)[[2]]
write.csv(percentage_confusion_dislineal, file="dislineal_tabanus_confusion.csv")
dis_lin_tabanus$error

dis_qua_tabanus <- quaDA(new_pcs_tabanus, factors_tabanus, validation="crossval")
percentage_confusion_discia <- NULL
for (i in 1:dim(dis_qua_tabanus$confusion)[1]) {
  row_new<- (dis_qua_tabanus$confusion[i, ]*100)/sum(dis_qua_tabanus$confusion[i, ])
  percentage_confusion_discia <- rbind(percentage_confusion_discia, row_new)
}
row.names(percentage_confusion_discia) <- dimnames(percentage_confusion_discia)[[2]]
write.csv(percentage_confusion_discia, file="discua_tabanus_confusion.csv")
dis_qua_tabanus$error

##Plot discriminants of Diachlorus in ggplot2
#goal_tabanus_lin <- NULL
#for (i in 1:dim(percentage_confusion_dislineal)[1]) {
#new_element<- percentage_confusion_dislineal[i, i]
#goal_tabanus_lin <- c(goal_tabanus_lin, new_element)
#}
#goal_tabanus_qua <- NULL
#for (i in 1:dim(percentage_confusion_disqua)[1]) {
#new_element<- percentage_confusion_disqua[i, i]
#  goal_tabanus_qua <- c(goal_tabanus_qua, new_element)
#}
#SS <- c("T.occidentalis", "T. triangulum")
#df1 <- data.frame(Species = factor(SS), accuracy_lineal = goal_diachlorus_lin, accuracy_cuadratica = goal_diachlorus_qua) 
#pdf(file="discri_diachlorus_crossval.pdf", width = 11, height = 5.7)
#d <- ggplot(df1,aes(x=Species,y=accuracy_lineal,fill=Species)) +  stat_summary(fun.y=mean,position=position_dodge(), geom="bar") + scale_fill_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#ffff33", "#984ea3", "#ff7f00")) + scale_x_discrete("Species") +  coord_flip() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)) + theme(legend.text = element_text(size = 14, face = "italic")) + theme(legend.title = element_text(size=16))
#d + geom_errorbar(aes(y=accuracy_cuadratica, ymax=accuracy_cuadratica, ymin=accuracy_cuadratica), linetype="dashed", size= 0.5, position=position_dodge()) + scale_y_continuous("Accuracy (%)", breaks = round(seq(0, 100, by = 5),1)) + theme(axis.ticks = element_blank(), axis.text.y = element_blank())
#dev.off()

###########################################################################################################
#####################################################################################################################################################################################
## 9. Riemmanian distances among species and Goodall's F test 

cons_shapes
names_tabanus <- combn(names(cons_shapes[-1]), m=2)
comb_tabanus <- combn(cons_shapes[-1], m=2)
r_distances <- NULL
nam <- NULL
for (i in 1:dim(comb_tabanus)[2]) {
  n <- paste(names_tabanus[,i][1], "to", names_tabanus[,i][2])  
  r <- riemdist(as.matrix(comb_tabanus[,i][1][[1]]), as.matrix(comb_tabanus[,i][2][[1]]))
  r_distances <- rbind(r_distances, r)
  nam <- c(nam, n)
}
rownames(r_distances) <- nam
colnames(r_distances) <- c("Riemman")
write.csv(r_distances, "Riemmanian_interspecific_distances.csv")

alg_names_comb <- combn(names(alg_inter[-1]), m=2)
comb_alg_tabanus <- combn(alg_inter[-1], m=2)
name <- NULL
t_goodall <- NULL
for (i in 1:dim(comb_alg_tabanus)[2]) {
  n <- paste(alg_names_comb[,i][1], "and", alg_names_comb[,i][2]) 
  t <- testmeanshapes(comb_alg_tabanus[,i][1][1][[1]][[1]], comb_alg_tabanus[,i][2][1][[1]][[1]], resamples = 1000, replace = FALSE, scale= TRUE)
  t_value <- c(t$G, t$G.pvalue)
  t_goodall <- rbind(t_goodall, t_value)
  name <- c(name, n)
}
rownames(t_goodall) <- name
colnames(t_goodall) <- c("G", "G_pvalue")
write.csv(t_goodall, "Goodall_interspecific_significance.csv")



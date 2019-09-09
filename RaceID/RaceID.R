library(RaceID)



ctrl.names <- paste("CTRL_",c(1:378), sep = "")
colnames(ctrl.data) <- ctrl.names

trt.names <- paste("TRT_",c(1:495), sep = "")
colnames(trt.data) <- trt.names

x <- merge(ctrl.data,trt.data, by="row.names", all = TRUE)
x[is.na(x)] <- 0

rownames(x) <- x$Row.names
x$Row.names <- NULL

barplot(colSums(x))
abline(h=1000, col = "red")

# name2id <- function(x,id) {
#   ##  id[sub("\\_\\_chr\\w+","",id) %in% x]
#   n <- c()
#   for ( j in x ){ n <- append(n,id[grep(paste(j,"(\\_\\_chr\\w+|$|\\|)",sep=""),id)])
#   }
#   n
# }

# CGenes <- c("Pcna","Mki67","Mir703","Gm44044","Gm22757","Gm4775","Gm17541","Gm8225","Gm8730","Ptma","Actb","Hsp90aa1","Hsp90ab1","Ppia")
# FGenes<- c("Malat1","Xist")

sc <- SCseq(x)
sc <- filterdata(sc,mintotal=2000, CGenes = NULL, FGenes = NULL) 
dim(sc@ndata)

sc <- compdist(sc,metric="pearson")
sc <- clustexp(sc, cln = 0, sat = T, clustnr = 30)
plotsaturation(sc,disp=FALSE)
plotsaturation(sc,disp=TRUE)
plotjaccard(sc)


#sc <- clustexp(sc,cln=14,sat=FALSE, clustnr = 50)
sc <- findoutliers(sc, probthr = 0.0001)
plotbackground(sc)
plotsensitivity(sc)
plotoutlierprobs(sc)


sc <- comptsne(sc)
sc <- compfr(sc,knn=10)
plotmap(sc)
plotmap(sc,fr=TRUE)

types <- sub("(\\_|\\.).+","", colnames(sc@ndata))
plotsymbolsmap(sc,types,fr=F)

plotexpmap(sc,"Pecam1",logsc=TRUE,fr=F)

# genes <- list()
# dg.all <- list()
# for (i in c(1:6)) {
#   dg <- clustdiffgenes(sc,i,pvalue=.05)
#   dg <- dg[grep("Rpl|Rps|Rik|RP|Malat|Jun|Fos|Hsp|Actb|Eef1a1|Ptma", rownames(dg), invert = T),]
#   dg <- dg[order(dg$fc, decreasing = T),]
#   dg.all[[i]] <- dg
#   genes[[i]] <- head(rownames(dg),5)
#   
# }
# 
# genes <- unique(unlist(genes))
# plotmarkergenes(sc,genes,cl=c(1:6),samples=types,order.cells=F,aggr = T, cap = 3, cluster_cols = T)

dg <- clustdiffgenes(sc,6,pvalue=.05)
dg <- dg[order(dg$fc, decreasing = T),]
head(dg,25)
types <- sub("(\\_|\\.).+","", colnames(sc@ndata))
genes <- head(rownames(dg)[dg$fc>1],10)
plotmarkergenes(sc,genes,samples=types)




ltr <- Ltree(sc)
ltr <- compentropy(ltr)
ltr <- projcells(ltr,cthr=15,nmode=T,fr=F)
#ltr <- projback(ltr,pdishuf=500)
ltr <- lineagegraph(ltr)
ltr <- comppvalue(ltr,pthr=0.01)
plotgraph(ltr,scthr=0.6,showCells=FALSE,showTsne=TRUE)
x <- compscore(ltr,scthr=0.6)



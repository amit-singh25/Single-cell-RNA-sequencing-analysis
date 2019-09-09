library(edgeR)
require(scran)
require(scater)
require(SingleCellExperiment)
require(Matrix)
require(Seurat)
library("dplyr")
library("DropletUtils")
library(ggplot2)
library(cowplot)
####10x genomics does not give genes.tsv rather it give features.tsv, there for one has to rename them when laod to seurat objects 
ctrl.data<- Read10X(data.dir = "~/Desktop/single_cell_analysis/con_agg/outs/filtered_feature_bc_matrix/")
trt.data<- Read10X(data.dir = "~/Desktop/single_cell_analysis/treat_agg/outs/filtered_feature_bc_matrix/")
#####
colnames(x = ctrl.data) <- paste('ctrl', colnames(x = ctrl.data), sep = '_')
colnames(x = trt.data) <- paste('trt', colnames(x = trt.data), sep = '_')
###creat seurat object
#for each, create object, add metadata, filter, normalize, and scale

s2=CreateSeuratObject(raw.data=ctrl.data,project="thymus",min.cells=5)
s2@meta.data$group="ctrl"
mito.genes <- grep(pattern = "^MT-", x = rownames(x = s2@data), value = TRUE)


###check the data
slotNames(s2)
head(s2@meta.data)


####filtred Filtering
lt =c(200,-Inf)
ht=c(2500,0.05)
s2=FilterCells(s2,subset.names="nGene",low.thresholds = lt,high.thresholds = ht)

#s2=FilterCells(s2,subset.names="nGene",low.thresholds = 500,high.thresholds = Inf)

low.thresholds =c(200,-Inf)


s2=NormalizeData(s2)
s2=ScaleData(s2,display.progress = T)


################
s4=CreateSeuratObject(raw.data=trt.data,project="thymus",min.cells=5)
s4@meta.data$group="treat"
s4=FilterCells(s4,subset.names="nGene",low.thresholds = 500,high.thresholds = Inf)
s4=NormalizeData(s4)
s4=ScaleData(s4,display.progress = T)


#select variable genes common to both samples
s2=FindVariableGenes(s2,do.plot=T)
s4=FindVariableGenes(s4,do.plot=T)
g.2=head(rownames(s2@hvg.info),1000)
g.4=head(rownames(s4@hvg.info),1000)
genes.use=unique(c(g.2,g.4))
genes.use=intersect(genes.use,rownames(s2@scale.data))
genes.use=intersect(genes.use,rownames(s4@scale.data))
########
agg=RunCCA(s2,s4,genes.use=genes.use,num.cc=30)
######plot 
####CCA space plot and violin plot of abundance per sample.
p1=DimPlot(object=agg,reduction.use="cca",group.by="group",pt.size=0.5,do.return=T)
p2=VlnPlot(object=agg,features.plot="CC1",group.by="group",do.return=T)
plot_grid(p1,p2)

###The PrintDim function outputs the top distinguishing genes in each CCA dimension.
PrintDim(object=agg,reduction.type="cca",dims.print=1:2,genes.print=10)
#####
MetageneBicorPlot(agg,grouping.var="group",dims.eval=1:30,display.progress=T)
####A heatmap is plotted to associate the most variable genes with each cluster.  
##For this example, I only plotted the first 9 CCs.
#####
DimHeatmap(object=agg,reduction.type="cca",cells.use=500,dim.use=1:9,do.balanced=T)
####
#With this you can now align the data to the CCA subspace–choosethe number of CC dimensions that make sense for your sample.  
#Note that each of these dimension reduction steps produces a new set of data under the @dr slot, 
#so you can refer to this for clustering.  After this point, you will have both “cca” and “cca.aligned” under this slot.

agg=AlignSubspace(agg,reduction.type="cca",grouping.var="group",dims.align=1:16)

#####Let’s plot distributions, as violin plots, for each of the first two CC dimensions.

p1=VlnPlot(object=agg,features.plot="ACC1",group.by="group",do.return=T)
p2=VlnPlot(object=agg,features.plot="ACC2",group.by="group",do.return=T)
plot_grid(p1,p2)
#####tSNE plot
agg=RunTSNE(agg,reduction.use="cca.aligned",dims.use=1:16,do.fast=T)
agg=FindClusters(agg,reduction.type="cca.aligned",resolution=0.6,dims.use=1:16)
p1=TSNEPlot(agg,do.return=T,pt.size=0.5,group.by="group")
p2=TSNEPlot(agg,do.label=T,do.return=T,pt.size=0.5)
plot_grid(p1,p2)
#######Principal Components Analysis
agg=RunPCA(agg,pc.genes=agg@var.genes,do.print=T,pcs.print = 1:5,genes.print = 5)
####
VizPCA(agg,pcs.use = 1:2) #plots component of each of top genes per PC
PCAPlot(agg,dim.1=1,dim.2=2) #all cells plotted on first two PCs
PCAPlot(agg,dim.1=1,dim.2=2,group.by="group") #show source samples
PCHeatmap(agg,pc.use=1,cells.use=500,do.balanced = T,label.columns = F) #first PC only
####
#try running more PCs to visualize how many explain variance
PCHeatmap(agg,pc.use=1:6,cells.use = 500,do.balanced = T,label.columns = F,use.full=F)

#try running more PCs to visualize how many explain variance
PCHeatmap(agg,pc.use=1:6,cells.use = 500,do.balanced = T,label.columns = F,use.full=F)
##Looks good.  Next, project all the data onto PC space for differential expression analysis
agg=ProjectPCA(agg,do.print = F)#Differential Expression by Sample

agg=SetAllIdent(agg,id="group")
cell.markers=FindMarkers(agg,ident.1="ctrl",ident.2="treat",test.use="wilcox")

treat.mrkrs=rownames(head(cell.markers[order(-cell.markers$pct.2 + cell.markers$pct.1),],5))
ctrl.mrkrs=rownames(head(cell.markers[order(cell.markers$pct.2 - cell.markers$pct.1),],5))

avg.cells=log1p(AverageExpression(agg,show.progress = F))
avg.cells$gene=rownames(avg.cells)

FeaturePlot(agg,features.plot = c(ctrl.mrkrs,treat.mrkrs),cols.use=c("grey","blue"),reduction.use="pca")

###JackStraw PCA
agg=JackStraw(agg,num.replicate=100,display.progress = T)
JackStrawPlot(agg,PCs=1:18) #to find how many are significant
PCElbowPlot(agg) #another, simpler way to visualize 
all.markers=FindAllMarkers(agg,only.pos=T,min.pct=0.25,thresh.use=0.25)
all.markers %>% group_by(cluster) %>% top_n(5,avg_logFC)
######
top2=all.markers %>% group_by(cluster) %>% top_n(2,avg_logFC)
SplitDotPlotGG(agg,genes.plot=as.character(top2$gene),cols.use = c("blue","red"),x.lab.rot = T,plot.legend = T,dot.scale = 8,do.return = T,grouping.var = "group")
#####
FeaturePlot(agg,features.plot = as.character(top2$gene),cols.use=c("grey","blue"),reduction.use="pca")

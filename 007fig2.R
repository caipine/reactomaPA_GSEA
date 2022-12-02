cd /rsrch3/scratch/lym_myl_rsch/qcai1/RNAseq/PI3K/scRNA
bsub -Is -q medium -W 10:00 -M 200 -R rusage[mem=200] -n 1  -XF  /bin/bash
conda activate R403


library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(monocle3)
library(SeuratWrappers)
library(ggpubr)
library(ggplot2)

library(scales)
library(RColorBrewer)
library(data.table)
library("Nebulosa")
library(viridis)
library(Seurat)
library(scCustomize)

gexwst <- readRDS("20221116.gexws.rds")
files <- dir("Data")
d0 <- readRDS(paste0("Data/",files[1]))


d1 <- readRDS(paste0("Data/",files[1]))
for ( i in 2:length(files)) {
#for ( i in 2:10) {
print(i)
d1 <-  rbind(d1, readRDS(paste0("Data/",files[i])))
}

dim(d1)
saveRDS(d1, file= "Data/2022_1130.d1.rds")

d1 <- readRDS("Data/2022_1130.d1.rds")

gexwst <- gexws
gexwst@meta.data <- cbind(gexwst@meta.data,t(d1))    ####add pathway to meta.data
#rownames(d1)

gexwsg1 <- subset(gexwst, subset = Phase %in% "G1")
library(Signac)
gexwsg1 <- SCTransform(gexwsg1 , method = "glmGamPoi",vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE, variable.features.n = nrow(gexwsg1) )
gexwsg1 <- RunTFIDF(gexwsg1) #normalization
gexwsg1 <- RunSVD(gexwsg1) #lsi
DefaultAssay(gexwsg1) <- "SCT"
gexwsg1 <- RunUMAP(gexwsg1, reduction = "lsi", dims = 1:30)
gexwsg1 <- FindNeighbors(gexwsg1, reduction = "umap", dims = 1:2, force.recalc = T)
gexwsg1 <- FindClusters(gexwsg1  , resolution = 0.05) #7 cluster

gexwsg1@meta.data$cytoTRACE = readRDS("20221111.CytoTRACE.results_gexwsg1.rds")$CytoTRACE
pal <- viridis(n = 10, option = "D")
ggData = data.table(UMAP1 = as.matrix(gexwsg1@reductions$umap@cell.embeddings)[,1], UMAP2 = as.matrix(gexwsg1@reductions$umap@cell.embeddings)[,2],
                    cytoTRACE = readRDS("20221111.CytoTRACE.results_gexwsg1.rds")$CytoTRACE, 
					orig.ident = gexwsg1@meta.data$orig.ident,  
					Cluster = Idents(gexwsg1),
					CD19 = gexwsg1@assays$SCT@data["CD19",]
					#,
					#Phase = gexwsg1@meta.data$Phase,
					#pseudotime = pseudotime(cds),
					#pseudotime2 = pseudotime(cds2)
					)
ggData[ggData$orig.ident %in% "M",]$orig.ident <- "MinoR" 
ggData[ggData$orig.ident %in% "V",]$orig.ident <- "5µM PBN"
ggData[ggData$orig.ident %in% "X",]$orig.ident <- "10µM PBN"

ggData$orig.ident <- factor(ggData$orig.ident, levels = c("MinoR","5µM PBN","10µM PBN")) 
ggData22 <- cbind(ggData,  gexwsg1@meta.data[,c("JAATINEN_HEMATOPOIETIC_STEM_CELL_UP",
                     "JAATINEN_HEMATOPOIETIC_STEM_CELL_DN",
					 "KEGG_PURINE_METABOLISM")]	)
head(ggData22)				 
ggData2 <- 	melt(ggData22, id = c("UMAP1", "UMAP2", "orig.ident", "Cluster")) 
ggData2_CD19 <-	ggData2[ggData2$variable %in% "CD19",]
ggData2_CD19 <- ggData2_CD19[order(ggData2_CD19$value),]
ggData2_CD19[ggData2_CD19$value ==0,]$value <- NA	 
p222(ggData2_CD19,"CD19", "CD19", pal) 


cytoColor = rev(brewer.pal(11, "Spectral")) # Reproduce cytoTRACE colours
cytoColor[6] = "gold"	
s = 12
pt = 0.1
col2 <- rev(brewer.pal(n = 11, name = "RdBu"))



#############				
f2A <- DimPlot_scCustom(seurat_object = gexwsg1, pt.size =pt) + 
theme_classic(base_size = s) + ggtitle(" Cluster") +
theme(  plot.title = element_text(hjust = 0.5))



p222 <- function(ggData2,f1, f2, colx ) {
return(
ggplot(ggData2[ ggData2$variable %in% f1,], aes(UMAP1, UMAP2, color = value)) +
  geom_point(size = pt) + theme_classic(base_size = s) +
  scale_color_gradientn(colors = colorRampPalette(colx)(50)) +
  facet_wrap(~  orig.ident) + 
theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()
) + 
    theme(legend.key.size = unit(0.5, 'cm')) +
	#theme(legend.position="none") +
  theme(strip.background = element_blank())+
ggtitle(f2) +
theme(  plot.title = element_text(hjust = 0.5),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),  
		axis.text.x=element_blank() #,         axis.ticks.x=element_blank() 
        )
)
}

  

p222N <- function(ggData2, f1, f2, colx ) {
return(
ggplot(ggData2[ ggData2$variable %in% f1,], aes(UMAP1, UMAP2, color = value)) +
  geom_point(size = pt) + theme_classic(base_size = s) +
  scale_color_gradientn(colors = colorRampPalette(colx)(50)) +
  facet_wrap(~  orig.ident) + 
theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()
) + 
    theme(legend.key.size = unit(0, 'cm')) +
	theme(legend.position="none") +
  theme(strip.background = element_blank())+
ggtitle(f2) +
theme(  plot.title = element_text(hjust = 0.5),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),  
		axis.text.x=element_blank() #,         axis.ticks.x=element_blank() 
        )
)
}



f2C_t <-  ggarrange(
p222N(ggData2, "cytoTRACE", "cytoTRACE", cytoColor),
p222N(ggData2_CD19,"CD19", "CD19", pal),
p222N(ggData2,"JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", "JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", col2),
p222N(ggData2,"JAATINEN_HEMATOPOIETIC_STEM_CELL_DN", "JAATINEN_HEMATOPOIETIC_STEM_CELL_DN", col2),
p222N(ggData2,"KEGG_PURINE_METABOLISM", "KEGG_PURINE_METABOLISM", col2),
nrow = 5, ncol =1)



f2C_legend <-  ggarrange(
as_ggplot(get_legend(p222(ggData2, "cytoTRACE", "cytoTRACE", cytoColor))),
as_ggplot(get_legend(p222(ggData2_CD19,"CD19", "CD19", pal))),
as_ggplot(get_legend(p222(ggData2,"JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", "JAATINEN_HEMATOPOIETIC_STEM_CELL_UP", col2))),
as_ggplot(get_legend(p222(ggData2,"JAATINEN_HEMATOPOIETIC_STEM_CELL_DN", "JAATINEN_HEMATOPOIETIC_STEM_CELL_DN", col2))),
as_ggplot(get_legend(p222(ggData2,"KEGG_PURINE_METABOLISM", "KEGG_PURINE_METABOLISM", col2))),
nrow = 5, ncol =1)

	
pdf(file = "pdf/20221202.gexwsg1_fig2_1T3.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 21.2) 
ggarrange(
fig2ABCD_T,
nrow = 2, ncol =1, heights = c(12,10)) +
  theme(plot.margin = margin(1,1,2,1, "cm")) 
dev.off()	
	

load("bulkRNA_results/2022_1202.p1-6.rdata")

s1=10
s2 =10
s3 =10
s4 = 6
f2D <- ggarrange(
p2+ theme(text = element_text(size = s1), axis.text.x = element_text(size = s2), axis.title.y = element_text(size = s2)) +
theme(  plot.title = element_text(hjust = 0.5, size=s4)),
p1+ theme(text = element_text(size = s1), axis.text.x = element_text(size = s2), axis.title.y = element_text(size = s2)) +
 theme(  plot.title = element_text(hjust = 0.5, size=s4),
        axis.title.y=element_blank()),
p3+ theme(text = element_text(size = s1), axis.text.x = element_text(size = s2), axis.title.y = element_text(size = s2))+
theme(  plot.title = element_text(hjust = 0.5, size=s4)),
p5+ theme(text = element_text(size = s1), axis.text.x = element_text(size = s2), axis.title.y = element_text(size = s2))+
 theme(  plot.title = element_text(hjust = 0.5, size=s4),
        axis.title.y=element_blank()),
p4+ theme(text = element_text(size = s1), axis.text.x = element_text(size = s2), axis.title.y = element_text(size = s2))+
theme(  plot.title = element_text(hjust = 0.5, size=s4)),
p6+ theme(text = element_text(size = s1), axis.text.x = element_text(size = s2), axis.title.y = element_text(size = s2))+
 theme(  plot.title = element_text(hjust = 0.5, size=s4),
        axis.title.y=element_blank()),
nrow = 3, ncol =2)
f2D


fig2ABCD_T <- ggarrange(
ggarrange(
f2C_t,
f2C_legend,
nrow = 1, ncol =2, widths = c(15,2.5)),

ggarrange(
f2A ,
ggarrange(
 plot_spacer(),
f2D ,
nrow = 1, ncol =2, widths = c(1,2)),

nrow = 2, ncol =1, heights = c(2,3)),

nrow = 1, ncol =2,  widths = c(2,2))	
	
	
pdf(file = "pdf/20221202.gexwsg1_fig2_1T4.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 21.2) 
ggarrange(
fig2ABCD_T,
plot_spacer(),
nrow = 2, ncol =1, heights = c(12,10)) +
  theme(plot.margin = margin(1,1,2,1, "cm")) 
dev.off()	
	

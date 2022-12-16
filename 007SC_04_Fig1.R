rm(list =ls())
library(Seurat)
library(scCustomize)
library(ggplot2)
library(viridis)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(patchwork)

pt = 0.1
s = 12
col2 <- rev(brewer.pal(n = 11, name = "RdBu"))
pal <- viridis(n = 10, option = "D")
cytoColor = rev(brewer.pal(11, "Spectral")) # Reproduce cytoTRACE colours
cytoColor[6] = "gold"	



gexwsg1a <- readRDS( "Data/20221205.gexwsg1.rds")
DimPlot(gexwsg1a)
gexwsg1a@meta.data$CytoTRACE = readRDS("20221111.CytoTRACE.results_gexwsg1.rds")$CytoTRACE
Idents(gexwsg1a)  = dplyr::recode(Idents(gexwsg1a) ,
                                                "1"="C3",
                                                "2"="C4",
                                                "3"="C5",
                                                "4"="C0",
                                                "5"="C2",
                                                "0"="C1"
                                                )
												
Idents(gexwsg1a)  = dplyr::recode(Idents(gexwsg1a) ,    
                                                "C0"="0",
                                                "C1"="1",
                                                "C2"="2",
                                                "C3"="3",
                                                "C4"="4",
                                                "C5"="5"
												)
Idents(gexwsg1a) <- factor(Idents(gexwsg1a), levels = 0:5)

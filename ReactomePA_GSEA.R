
```{r}
rm(list =ls())
```

```{r}
ddsObj <- readRDS( "20220922.ddsObj.delMCL89_delMCL21.rds")
normalizedCounts <- counts(ddsObj, normalized=TRUE) 
shrinkLvV <- ddsShrink <- readRDS("20220922.ddsShrink.delMCL89_delMCL21.rds" )

shrinkLvV$Significant <- "Not_Sig"
shrinkLvV[is.na(shrinkLvV$padj),]$padj <- 1
shrinkLvV[shrinkLvV$padj < 0.05,]$Significant <- "0.01 < FDR < 0.05"
shrinkLvV[shrinkLvV$padj < 0.01,]$Significant <- "FDR < 0.01"

rld  <- readRDS("20220922.rld.delMCL89_delMCL21.rds" )
sinfo <- data.frame(colData(ddsObj))
pathways.all <- readRDS("C:/Users/qcai1/OneDrive - Inside MD Anderson/007analysis/RNAseq/007B1/results/v2022.pathways.all.rds")
shrinkLvV$Symbol <- rownames(shrinkLvV)
res2 <- shrinkLvV %>% 
  dplyr::select(Symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Symbol) %>% 
  summarize(stat=mean(stat))


library(fgsea)
ranks <- deframe(res2)

#pathways.hallmark %>%    head() %>%   lapply(head)
  
fgseaRes <- fgsea(pathways = pathways.all, stats=ranks, nperm=1000)

###all pathways
fgseaResTidy <- fgseaRes %>%   as_tibble() %>%    arrange(desc(NES))

fgseaResTidy$R <- NA
for ( i in 1: nrow(fgseaResTidy))
fgseaResTidy$R[i] <- length(fgseaResTidy$leadingEdge[[i]])/length(intersect( shrinkLvV$Symbol,
                                                                       pathways.all[[fgseaResTidy$pathway[i]]]) )

max(fgseaResTidy$R)
min(fgseaResTidy$R)
fgseaResTidy[order(fgseaResTidy$R, decreasing = T),]
t1 <- fgseaResTidy[fgseaResTidy$padj < 0.05 & fgseaResTidy$NES < 0,]
t1[t1$R < 0.3,]
```


```{r, fig.align="center", fig.width=8, fig.height=8, echo=FALSE, echo=FALSE, warning=FALSE, message=FALSE}
print_voc1 <- function(shrinkLvV, pathid="BOQUEST_STEM_CELL_UP") {
  print(
ggplot(shrinkLvV, aes(x = log2FoldChange, y= -log10(padj) )) + 
  geom_point(aes(colour= Significant), shape=20, size=1)+
  scale_color_manual(values = c("green","red","grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  #geom_hline(yintercept=2, colour="blue", linetype="dashed") +
  #  geom_hline(yintercept=-log10(0.05), colour="blue", linetype="dashed") +
  geom_vline(xintercept = -1, colour="grey", linetype="dashed") +
    geom_vline(xintercept = 1, colour="grey", linetype="dashed") +
   #coord_cartesian(xlim = c(-5,5)) +
    geom_text_repel(data = subset(shrinkLvV, ( Symbol %in% fgseaResTidy[fgseaResTidy$pathway %in% pathid,]$leadingEdge[[1]])), #& pvalue < 0.05
    aes(label = Symbol), size = 3, box.padding = unit(0.1, "lines"),point.padding = unit(0.1, "lines")) +
    ylab("-log10(False Discovery Rate)") +
   coord_cartesian(ylim = c(0,8), xlim = c(-1.5,1.5))
  )
}
```

```{r, fig.align="center", fig.width=8, fig.height=8, echo=FALSE, echo=FALSE, warning=FALSE, message=FALSE}
print_voc2 <- function(shrinkLvV, pathid="BOQUEST_STEM_CELL_UP") {
  print(
    ggplot(shrinkLvV, aes(x = log2FoldChange, y= -log10(padj) )) + 
  geom_point(aes(colour= Significant), shape=20, size=1)+
  scale_color_manual(values = c("green","red","grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  #geom_hline(yintercept=2, colour="blue", linetype="dashed") +
  #  geom_hline(yintercept=-log10(0.05), colour="blue", linetype="dashed") +
  geom_vline(xintercept = -1, colour="grey", linetype="dashed") +
    geom_vline(xintercept = 1, colour="grey", linetype="dashed") +
   #coord_cartesian(xlim = c(-5,5)) +
    geom_text_repel(data = subset(shrinkLvV, ( Symbol %in% pathways.all[[pathid]])), #& pvalue < 0.05
    aes(label = Symbol), size = 3, box.padding = unit(0.1, "lines"),point.padding = unit(0.1, "lines")) +
    ylab("-log10(False Discovery Rate)") +
   coord_cartesian(ylim = c(0,8), xlim = c(-1.5,1.5))
  )
}
```


```{r, fig.align="center", fig.width=8, fig.height=8, echo=FALSE, echo=FALSE, warning=FALSE, message=FALSE}
id <- "PEREZ_TP53_TARGETS"
print_voc1(shrinkLvV, pathid=id)
print_voc2(shrinkLvV, pathid=id)
```



```{r}
pst <- function(dx) {
  #dx <- fgseaResTidy.MinoR_10
  a <- c()
  for (i in 1:nrow(dx)) {
  a <-  c(a,str_c(dx$leadingEdge[[i]], collapse =  "/"))
#  print(str_c(dx$leadingEdge[[i]], collapse =  "/"))
  }
  return(a)
}

library( stringi )


pst_GeneRatio <- function(dx) {
  #dx <- fgseaResTidy.MinoR_10
  a <- c()
  for (i in 1:nrow(dx)) {
  a <-  c(a, paste0( length(dx$leadingEdge[[i]]), "/", length(dx$leadingEdge)))
  }
  return(a)
}

pst_GeneRatio_inpathway <- function(dx) {
a <- c()
for (i in 1:nrow(dx)) {
  #print(length(dx$leadingEdge[[i]])/length(pathways.all[[dx$pathway[i]]]))
  a <-  c(a,
          length(dx$leadingEdge[[i]])/length(pathways.all[[dx$pathway[i]]])
          )
}
return(a)
}
```

shrink = TRUE, 

```{r}
 keyword1 = NULL
is.null
keyword1 == NULL
```

```{r}
build_res_from_feas_pst <- function(shrinkLvV, fgseaResTidy, padj_cut = 0.25,  keyword1 = NULL, keyword2=NULL , UPonly = FALSE, DNonly = FALSE )
    {
  UP_genes <- rownames(shrinkLvV[shrinkLvV$log2FoldChange > 0 &shrinkLvV$padj < 0.05,])
  DN_genes <- rownames(shrinkLvV[shrinkLvV$log2FoldChange < 0 &shrinkLvV$padj < 0.05,])
geneClusters <- list(UP_genes,DN_genes)

names(geneClusters) <- c("UP", "DOWN" )

  
UP <- fgseaResTidy[fgseaResTidy$padj < padj_cut & fgseaResTidy$NES > 0,]
DOWN <- fgseaResTidy[fgseaResTidy$padj < padj_cut & fgseaResTidy$NES < 0,]
#UP$geneID <- pst_more_up(UP, UP_genes = UP_genes )
#DOWN$geneID <- pst_more_dn(DOWN,  DN_genes = DN_genes )
UP$geneID <- pst(UP)
DOWN$geneID <- pst(DOWN)
UP$Cluster <- "UP"
DOWN$Cluster <- "DOWN"

dx <- rbind(UP, DOWN)

if (UPonly) 
  dx <- UP

if (DNonly) 
  dx <- DOWN

if (!is.null(keyword1)) 
dx1 <- dx[grepl(keyword1, dx$pathway),]

if (!is.null(keyword2)) 
dx2 <- dx[grepl(keyword2, dx$pathway),]

if (!is.null(keyword1) & !is.null(keyword2) )
dx <- rbind(dx1, dx2)


if (!is.null(keyword1) & is.null(keyword2) )
dx <- dx1

if (is.null(keyword1) & !is.null(keyword2) )
dx <- dx2



dx$ID <- dx$pathway
dx$Description <- dx$pathway

dx$pvalue <- dx$pval
dx$p.adjust <- dx$padj
dx$Count <- dx$size
dx$GeneRatio <- pst_GeneRatio(dx)
dx$GeneRatio_inpath  <- pst_GeneRatio_inpathway(dx)

#sort(dx$GeneRatio_inpath )
dx <- dx[dx$Count > 5 & dx$GeneRatio_inpath > 0.05, ]
dx <- dx[,c("Cluster" ,    "ID"    ,      "Description" , "GeneRatio",   "pvalue"  ,    "p.adjust"    ,   "geneID"  ,    "Count" )]
#sort(dx$Count)
res <- new("compareClusterResult",
           compareClusterResult = dx, 
         #  geneClusters = geneClusters, 
           .call = match.call(expand.dots = TRUE))
    res@keytype <- "UNKNOWN"
    res@readable <- FALSE
    res@fun <-"enrichX" #paste("enrich",keyword)
    res@termsim <- as.matrix(0,0)
    res@geneClusters <- geneClusters
    return(res)
}
```

```{r}
source("https://raw.githubusercontent.com/YuLab-SMU/enrichplot/c50579c780158330dd8f828c31580bc187d3f744/R/emapplot_utilities.R")
source("https://raw.githubusercontent.com/YuLab-SMU/enrichplot/01dfdd27acc02ca57e7103fe5689a892c06eebd5/R/utilities.R")
library(igraph)
library(ggraph)
library(scatterpie)

emapplot_cqs <- function(x,
showCategory = 30,
layout = NULL,
coords = NULL,
split = NULL,
pie = "equal",
legend_n = 5,
cex_category = 1,
cex_line = 1,
min_edge=0.05,
cex_label_category  = 1 ,
shadowtext = TRUE ,
with_edge = TRUE,
group_category = FALSE,
label_format = 30,
group_legend = FALSE,
node_label  = "category",
label_style = "shadowtext" ,
repel = FALSE ,
cex_label_group = 1,
nWords = 4 ,
clusterFunction = stats::kmeans,
nCluster = NULL ,
cex_pie2axis = 1,
color="p.adjust") {
  has_pairsim(x)
    label_size_category <- 3
    label_group <- 3
    # y <- get_selected_category(showCategory, x, split)
    y <- fortify(x, showCategory = showCategory,
                 includeAll = TRUE, split = split)
    y$Cluster <- sub("\n.*", "", y$Cluster)

    if ("core_enrichment" %in% colnames(y)) { ## for GSEA result
        y$geneID <- y$core_enrichment
    }
    ## Data structure transformation, combining the same ID (Description) genes
    mergedEnrichDf <- merge_compareClusterResult(y)
     enrichDf = y

    
    segment.size <- get_ggrepel_segsize()
    geneSets <- setNames(strsplit(as.character(mergedEnrichDf$geneID), "/",
                              fixed = TRUE), mergedEnrichDf$ID) 
    
  pair_sim = x@termsim
  method = x@method
## get ggraph object and add edge
  
  
  
enrichDf=mergedEnrichDf
#pair_sim = x@termsim
#method = "JC"
    if (!is.numeric(min_edge) | min_edge < 0 | min_edge > 1) {
    	stop('"min_edge" should be a number between 0 and 1.')
    }

    if (is.null(dim(enrichDf)) | nrow(enrichDf) == 1) {  # when just one node
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- as.character(enrichDf$Description)
        V(g)$color <- "red"
        return(g)
    } else {
        w <- pair_sim[as.character(enrichDf$Description), 
            as.character(enrichDf$Description)]
    }

    wd <- melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    # remove NA
    wd <- wd[!is.na(wd[,3]),]
    if (method != "JC") {
        # map id to names
        wd[, 1] <- enrichDf[wd[, 1], "Description"]
        wd[, 2] <- enrichDf[wd[, 2], "Description"]
    }

    g <- graph.data.frame(wd[, -3], directed=FALSE)
    E(g)$width <- sqrt(wd[, 3] * 5) * cex_line
    # Use similarity as the weight(length) of an edge
    E(g)$weight <- wd[, 3]
    g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == enrichDf$Description)))
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt#/sapply(pathways.all[names(cnt)], length) ######################################
    colVar <- enrichDf[idx, color]
    V(g)$color <- colVar
    
    
     if(is.null(dim(enrichDf)) | nrow(enrichDf) == 1) {
        title <- enrichDf$Cluster
        p <- ggraph(g, "tree") + geom_node_point(color="red", size=5 * cex_category) +
            geom_node_text(aes_(label=~name)) + theme_void() +
            labs(title=title)
    #    return(p)
    }

    if(is.null(dim(mergedEnrichDf)) | nrow(mergedEnrichDf) == 1) {
        p <- ggraph(g, "tree")
        ID_Cluster_mat <- prepare_pie_category(enrichDf = enrichDf, pie=pie)

        ID_Cluster_mat <- cbind(ID_Cluster_mat,1,1,0.1*cex_category)
        colnames(ID_Cluster_mat) <- c(colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
            "x", "y", "radius")


        p <- p + geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
                cols=names(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],
                color=NA)+
            xlim(-3,3) + ylim(-3,3) + coord_equal()+
            geom_node_text(aes_(label=~name), repel=TRUE, segment.size = segment.size) +
            theme_void()+labs(fill = "Cluster")
    #    return(p)

    }
    p <- adj_layout(g = g, layout = layout, coords = coords)    
    ## add edges
    if (with_edge & length(E(g)$width) > 0) {
        p <- p + geom_edge_link(alpha=.8, aes_(width=~I(width)),
                                colour='darkgrey')
    }
    
       cnt <- sapply(geneSets[idx], length)
    p$data$name == names(cnt)
    p$data$size <- cnt/sapply(pathways.all[names(cnt)], length)


if (is.null(dim(y)) | nrow(y) == 1 | is.null(dim(mergedEnrichDf)) | nrow(mergedEnrichDf) == 1)
        return(p)

    # ggData <- p$data
    # if show group cricle or group label, Process p$data and assign color to the group label
    if (group_category || node_label == "all" || node_label == "group") {    
        ggData <- groupNode(p = p, enrichDf = y, nWords = nWords, 
            clusterFunction =  clusterFunction, nCluster = nCluster)
        p$data <- ggData
    }      
    ## add circle
    if (group_category) {
        p <- add_ellipse(p = p, group_legend = group_legend, 
            label_style = label_style, ...)
    }
       
    ## then add the pie plot
    ## Get the matrix data for the pie plot
    ID_Cluster_mat <- get_pie_data(enrichDf = y, pie = pie, mergedEnrichDf = mergedEnrichDf, cex_pie2axis = cex_pie2axis, 
                                   p = p, cex_category = cex_category)

    
    get_pie_data2 <-  function(enrichDf, pie, mergedEnrichDf, cex_pie2axis, p, cex_category) {
    ggData <- p$data
    ID_Cluster_mat <- prepare_pie_category(enrichDf = enrichDf, pie=pie) 
    desc <- mergedEnrichDf$Description[match(rownames(ID_Cluster_mat),
                                      mergedEnrichDf$Description)]
    i <- match(desc, ggData$name)
    ID_Cluster_mat$x <- ggData$x[i]
    ID_Cluster_mat$y <- ggData$y[i]
    ID_Cluster_mat$radius <- ggData$size[i]/2#sqrt(ggData$size[i] / sum(ggData$size) * cex_category * cex_pie2axis)
    return(ID_Cluster_mat)
}
    ## then add the pie plot
    ## Get the matrix data for the pie plot
    ID_Cluster_mat <- get_pie_data2(enrichDf = y, pie = pie, mergedEnrichDf = mergedEnrichDf, cex_pie2axis = cex_pie2axis, 
                                   p = p, cex_category = cex_category)

    
    add_pie_node2 <- function(p, ID_Cluster_mat, node_label, 
                         cex_category, cex_pie2axis,
                         cex_label_category,
                         shadowtext, legend_n,
                         label_size_category) {
    color <- NULL
    if(ncol(ID_Cluster_mat) > 4) {
        p <- p + ggnewscale::new_scale_fill() + 
            geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
            cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +
            coord_equal() 

        if (node_label == "all" || node_label == "category") {
            p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
                cex_label_node = cex_label_category, shadowtext = shadowtext)
        }
        
        p <- p + theme_void() +
            geom_scatterpie_legend(ID_Cluster_mat$radius, 
                x=min(ID_Cluster_mat$x), y=min(ID_Cluster_mat$y),
                n = legend_n,
                labeller=function(x)x*2 #round(sum(p$data$size) * x^2 / cex_category/ cex_pie2axis)
                ) +
            labs(fill = "Cluster")
    } else {
        title <- colnames(ID_Cluster_mat)[1]
        p <- p + theme_void() + geom_node_point(aes_(color=~color, size=~size))
        if (node_label == "all" || node_label == "category") {
            p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
                cex_label_node = cex_label_category, shadowtext = shadowtext)
        }
        p <- p + scale_color_continuous(low="red", high="blue", name = color,
                                        guide=guide_colorbar(reverse=TRUE)) +
            scale_size(range=c(3, 8) * cex_category)  +labs(title= title)  
    }  
    return(p)
}
    
    
    ## add dot and node label
    p <- add_pie_node2(p = p, ID_Cluster_mat = ID_Cluster_mat, 
                  node_label = node_label, cex_category = cex_category,
                  cex_pie2axis = cex_pie2axis, 
                  cex_label_category = cex_label_category,
                  shadowtext = shadowtext, legend_n = legend_n,
                  label_size_category = label_size_category)
   

    ## add group label
    if (node_label == "all" || node_label == "group") {   
        label_location <- get_label_location(ggData = ggData, label_format = label_format)
        p <- add_group_label(label_style = label_style, repel = repel, shadowtext = shadowtext, p = p,
            label_location = label_location, label_group = label_group,
            cex_label_group = cex_label_group)
    }    
    
    return(p)
}
```

```{r, fig.width=10, fig.height=10}
library("enrichplot")
res <- build_res_from_feas_pst(shrinkLvV, fgseaResTidy, padj_cut = 0.05)
x_gsea <- pairwise_termsim(res, method = "JC") # method = "JC"
emapplot(x_gsea,  min_edge = 0.05, showCategory = 30,  layout = "nicely")
emapplot_cqs(x_gsea,  min_edge = 0.05, showCategory = 30,  layout = "nicely")
```


```{r, fig.width=8, fig.height=8}
library("enrichplot")
res <- build_res_from_feas_pst(shrinkLvV, fgseaResTidy, padj_cut = 0.05,keyword1 = "DIFF", keyword2 = NULL)
x_gsea <- pairwise_termsim(res, method = "JC") # method = "JC"
emapplot(x_gsea,  min_edge = 0.05, showCategory = 40,  layout = "nicely")
emapplot_cqs(x_gsea,  min_edge = 0.02, showCategory = 40,  layout = "nicely")
```


```{r, fig.width=8, fig.height=8}
library("enrichplot")
res <- build_res_from_feas_pst(shrinkLvV, fgseaResTidy, padj_cut = 0.05,keyword1 = "STEM_CELL", keyword2 = "HSC")
x_gsea <- pairwise_termsim(res, method = "JC") # method = "JC"
emapplot(x_gsea,  min_edge = 0.05, showCategory = 40,  layout = "nicely")
emapplot_cqs(x_gsea,  min_edge = 0.02, showCategory = 40,  layout = "nicely")
```

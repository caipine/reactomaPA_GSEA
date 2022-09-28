library( stringi )
library(igraph)
library(ggraph)
library(scatterpie)
library(ReactomePA)
source("https://raw.githubusercontent.com/YuLab-SMU/enrichplot/c50579c780158330dd8f828c31580bc187d3f744/R/emapplot_utilities.R")
source("https://raw.githubusercontent.com/YuLab-SMU/enrichplot/01dfdd27acc02ca57e7103fe5689a892c06eebd5/R/utilities.R")


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

build_res_from_feas_pst <- function(shrinkLvV, fgseaResTidy, padj_cut = 0.25,  keywords = NULL, ex_keywords= NULL, UPonly = FALSE, DNonly = FALSE )
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

if (!is.null(keywords)) {
  dx1 <- dx[grepl(keywords[1], dx$pathway),]
  if (length(keywords) > 1) {
    for (i in 2:length(keywords)) {
      dx1 <- rbind(dx1,  dx[grepl(keywords[i], dx$pathway),] )
    }
  } 
  dx <- dx1
} 


if (!is.null(ex_keywords)) {
    for (i in 1:length(ex_keywords)) {
      dx <- dx[!grepl(ex_keywords[i], dx$pathway),] 
    }
} 


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


library( stringi )

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

build_res_from_feas_pst <- function(shrinkLvV, fgseaResTidy, padj_cut = 0.25,  keywords = NULL , UPonly = FALSE, DNonly = FALSE )
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



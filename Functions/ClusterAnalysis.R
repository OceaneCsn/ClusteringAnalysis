library(ggplot2)
library(gridExtra)
library(reshape2)

suppressMessages(library(coseq, warn.conflicts = F, quietly = T))


load("./Data/OntologyAllGenes.RData")
load("./Data/filteredData.RData")
load("./Data/NitrateGenes.RData")

# DEGs <- DEGs[["cnF CnF"]]
# genes <- unique(unlist(DEGs))
# cluster <- clustering(genes, data, nb_clusters = 7:14, norm = "TMM")

translate <- function(text){
  res = ""
  if(grepl("c", text)){res = paste0(res, "AmbientCO2, ")}
  else{res = paste0(res, "ElevatedCO2, ")}
  if(grepl("N", text)){res = paste0(res, "Hight Nitrate, ")}
  else{res = paste0(res, "Low Nitrate, ")}
  if(grepl("f", text)){res = paste0(res, "Iron Starvation")}
  else{res = paste0(res, "Iron Supply")}
  return(res)
}

plotProfile <- function(cluster, k="none"){
  # plot all the profiles or the profile of cluster k
  results <- cluster[[2]]
  clusters <- cluster[[1]]
  profiles <- data.frame(results@y_profiles)
  profiles$gene <- rownames(profiles)
  d <- melt(profiles)
  d$group <- str_split_fixed(d$variable, '_', 2)[,1]
  d$cluster <- clusters[match(d$gene, names(clusters))]
  if(k=="none"){
     g <- ggplot(data = d, aes(x=group, y=value)) + geom_boxplot(alpha=0.7, aes( color = group), fill = "grey", outlier.color = "black",outlier.alpha =0.1)  + geom_jitter(width = 0.1, alpha=0.0) +facet_wrap(~cluster, nrow=4) 
  }
  else{
    g <- ggplot(data = d[d$cluster==k,], aes(x=group, y=value)) + geom_boxplot(outlier.alpha =0.2, outlier.color = "black", alpha=0.7, aes(color = group), fill = "grey")  + geom_jitter(width = 0.1, alpha=0.0) 
  }
  g +theme(plot.title = element_text(size=22, face="bold"),strip.text.x = element_text(size = 20),
           legend.title = element_text(size = 25, face="bold"), legend.text = element_text(size=20),
           axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 0, hjust = 0, colour = "grey50"),
           axis.title=element_text(size=24)) + xlab("") + ylab("Normalized expression") + scale_colour_discrete("Condition", labels=sapply(levels(as.factor(d$group)),translate)) +
    stat_summary(fun.y=median, geom="line", aes(group=1), alpha=0.1, size = 1.5) +
    ylim(0, 0.25) 
}
findNitrateGenes <- function(cluster, k="none"){
  if(k=="none"){
    genesK <- names(cluster[[1]])
  }
  else{
    genesK <- names(cluster[[1]][cluster[[1]]==k])
  }
  res <- data.frame(Gene = genesK)
  res$Wang_2004 = ifelse(res$Gene %in% nGenes$Wang_2004, 1, 0)
  res$Marchive_2013 = ifelse(res$Gene %in% nGenes$Marchive_2013_20min, 1, 0)
  res$Widiez_2011 = ifelse(res$Gene %in% nGenes$HN_induced_Widiez_2011, 1, 0)
  res[,c("Name", "Description")] <- ontologies[match(res$Gene, ontologies$ensembl_gene_id),c("external_gene_name", "description")]
  res$NitrateScore <- rowSums(res[,grepl("_", colnames(res))])
  res <- res[order(-res$NitrateScore),]
  res$cluster <- cluster[[1]][match(res$Gene, names(cluster[[1]]))]
  return(res)
}

rankClusters <- function(cluster){
  clusterStats <- c()
  clusterStatsAbs <- c()
  for( k in seq(1:max(cluster[[1]]))){
    tab <- findNitrateGenes(cluster, k)
    
    clusterStats <- c(clusterStats, sum(tab$NitrateScore)/dim(tab)[1]/3)
    
    clusterStatsAbs <- c(clusterStatsAbs, sum(tab$NitrateScore)/3)
  }
  names(clusterStats) <- seq(1:max(cluster[[1]]))
  names(clusterStatsAbs) <- seq(1:max(cluster[[1]]))
  
  
  d <- data.frame(Cluster = names(clusterStats), Enrichment.Rate = clusterStats, Absolute.Enrichment = clusterStatsAbs)
  d <- melt(d)
  ggplot(data=d, aes(x=Cluster, y=value, fill = Cluster)) + geom_bar(stat="identity",alpha=0.3, color="black")+ facet_wrap(~variable, nrow=1, scales="free" )+
    theme(plot.title = element_text(size=22, face="bold"),strip.text.x = element_text(size = 26), legend.position = "none",
          legend.title = element_text(size = 25, face="bold"), legend.text = element_text(size=20),
          axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 20, hjust = 0, colour = "grey50"),
          axis.title=element_text(size=17)) + coord_flip()
}

#findNitrateGenes(5,cluster)

####### generation des ontologies de genes nitrates
# nitrateGenes <- unique(c(as.vector(nGenes$Wang_2004), as.vector(nGenes$Marchive_2013_20min), as.vector(nGenes$HN_induced_Widiez_2011)))
# nitrateGenes <- nitrateGenes[grepl("AT", nitrateGenes)]
# ontologies <- OntologyProfile(ids = AGIToEntrez$ensembl_gene_id)
# save(ontologies, file = "./Data/OntologyAllGenes.RData")

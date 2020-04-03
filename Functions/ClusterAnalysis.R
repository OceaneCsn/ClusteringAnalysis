library(ggplot2)
library(gridExtra)
library(reshape2)

suppressMessages(library(coseq, warn.conflicts = F, quietly = T))

library(plotly)
library(stringr)

load("./Data/OntologyAllGenes.RData")
load("./Data/filteredData.RData")
load("./Data/NitrateGenes.RData")

# load("D:/These/ClusteringAnalysis/Clusterings/CO2NoIronStarv.RData")


translate <- function(text){
  res = ""
  if(grepl("c", text)){res = paste0(res, "AmbientCO2,")}
  else{res = paste0(res, "ElevatedCO2,")}
  if(grepl("N", text)){res = paste0(res, "HightNitrate,")}
  else{res = paste0(res, "LowNitrate,")}
  if(grepl("f", text)){res = paste0(res, "IronStarvation")}
  else{res = paste0(res, "IronSupply")}
  return(res)
}

plotProfile <- function(cluster, k="none", boxplot=T, expression = "profiles"){
  # plot all the profiles or the profile of cluster k
  results <- cluster[[2]]
  clusters <- cluster[[1]]
  if(expression=="profiles"){profiles <- data.frame(results@y_profiles)
  ylab <- "Normalized expression/Mean(Normalized expression)"}
  if(expression=="counts") {
    profiles <- data.frame(log(as.matrix(results@tcounts)+1))
    ylab <- "log(Normalized expresion)"
  }
  profiles$gene <- rownames(profiles)
  d <- melt(profiles)
  d$group <- str_split_fixed(d$variable, '_', 2)[,1]
  d$cluster <- clusters[match(d$gene, names(clusters))]
  d$geneRep <- paste0(d$gene, substr(d$variable,4,5))
  if(k=="none"){
     g <- ggplot(data = d, aes(x=group, y=value))  +facet_wrap(~cluster, nrow=3) 
  }
  else{
    g <- ggplot(data = d[d$cluster==k,], aes(x=group, y=value))
  }
  if(boxplot) g <- g + geom_boxplot(alpha=0.7, lwd=1.2, aes( color = group), fill = "grey", outlier.color = "black",outlier.alpha =0.1)  + geom_jitter(width = 0.1, alpha=0.0015)
  else{
    g <- g+ geom_line(alpha=0.09,lwd=1.2, color="#333366", aes(group=geneRep))
  }
  
  g <- g +theme(plot.title = element_text(size=22, face="bold"),strip.text.x = element_text(size = 20),legend.position="bottom",
           legend.title = element_text(size = 2, face="bold"), legend.text = element_text(size=15, angle=0),
           axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 0, hjust = 0, colour = "grey50"),legend.text.align=1,
           axis.title=element_text(size=24)) + xlab("") + ylab(ylab) + scale_colour_discrete("", labels=sapply(levels(as.factor(d$group)),translate)) +
    stat_summary(fun.y=median, geom="line", aes(group=1), alpha=0.1, size = 1.5) 
  if(expression=="profiles") g <- g + ylim(0, 0.25) 
  g
}


findNitrateGenes <- function(cluster, k="none"){
  if(k=="none"){
    genesK <- names(cluster[[1]])
  }
  else{
    genesK <- names(cluster[[1]][cluster[[1]]==k])
  }
  res <- data.frame(Gene = genesK)
  for(paper in colnames(nGenes)){
    res[,paper] <- ifelse(res$Gene %in% toupper(nGenes[,paper]), 1, 0)
  }
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
  ggplotly(ggplot(data=d, aes(x=Cluster, y=value, fill = Cluster)) + geom_bar(stat="identity",alpha=0.3, color="black")+ facet_wrap(~variable, nrow=1, scales="free" )+
    theme(plot.title = element_text(size=22, face="bold"),strip.text.x = element_text(size = 26), legend.position = "none",
          legend.title = element_text(size = 25, face="bold"), legend.text = element_text(size=20),
          axis.text.y = element_text(size = 18, angle = 30), axis.text.x = element_text(size = 20, hjust = 0, colour = "grey50"),
          axis.title=element_text(size=17)) + coord_flip())
}


glmCluster <- function(DEgenes, normalized.count){
  
  glmData <- melt(round(normalized.count[DEgenes,], 0))
  print(head(glmData))
  colnames(glmData) <- c("Condition", "Counts")
  glmData <- glmData[sample(rownames(glmData)),]
  
  groups <- str_split_fixed(glmData$Condition, "_", 2)[,1]
  glmData$Co2 <- str_split_fixed(groups, "", 3)[,1]
  glmData$nitrate <- str_split_fixed(groups, "", 3)[,2]
  glmData$fer <- str_split_fixed(groups, "", 3)[,3]
  
  
  glmData$Co2 <- as.factor(ifelse(glmData$Co2 == "c", 0, 1))
  glmData$nitrate <- as.factor(ifelse(glmData$nitrate == "N", 0, 1))
  glmData$fer <- as.factor(ifelse(glmData$fer == "F", 0, 1))
  #glmData <- glmData[c("Counts", "Co2", "nitrate", "fer", "gene")]
  
  formula = "Counts ~ "

  for (factor in c("Co2", "nitrate", "fer")){
    if(length(levels(glmData[,factor]))>1){formula <- paste(formula, factor, '*')}
  }
  
  formula <- substr(formula, 1, nchar(formula)-1)
  
  glm <- glm(formula , data = glmData, family = poisson(link="log"))
  return(glm)
  
}

plotGlmCluster <- function(glm){
  coefs <- glm$coefficients[2:length(glm$coefficients)]
  d <- data.frame(Coefficient =str_remove_all(as.character(names(coefs)), "1"), Values = coefs)
  ggplotly(ggplot(d, aes(x=Coefficient, y=Values, fill=Coefficient)) + geom_bar(color = 'black', stat = "identity", alpha = 0.4) +
             theme(axis.text.x = element_text(size=25, angle = 320), legend.position = "none", 
                   axis.title.x = element_blank(), axis.title.y = element_blank(), 
                   plot.title = element_text(size=22, face="bold"), axis.text.y = element_text(size=20))+
             ggtitle("Generalized linear model's coefficients values"), tooltip = c("Values"))
}


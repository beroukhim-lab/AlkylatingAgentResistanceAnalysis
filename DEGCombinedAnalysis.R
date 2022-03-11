#DEGCombinedAnalysis

########## Load Packages ########## 
library(data.table)
library(ggplot2)
library(ggVennDiagram)
library(tools)
library(prob)

directory <- "/Volumes/xchip_beroukhimlab" #personal computer
#directory <-"/xchip/beroukhimlab" #server

######### ggplot Theme ########
theme <- theme(panel.background = element_blank(),
               panel.border=element_rect(color = "black", size = 1, fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.x=element_text(colour="black", size = 15),
               axis.text.y=element_text(colour="black", size = 15),
               axis.title.y = element_text(colour = "black", size = 15),
               axis.title.x = element_text(colour = "black", size = 15),
               axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"),
               legend.title = element_blank(),
               legend.key = element_blank(),
               plot.title = element_text(face = "bold", hjust = 0.5),
               legend.background = element_rect(fill=alpha('blue', 0)),
               legend.text = element_text(size=14))


plotSignificantGenes <- function(DatabaseName)
{
  #DatabaseName <- "PRISM"
  AlkylatingAgents <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Combined Datasets/AlkylatingAgents.csv', directory)) 
  AlkylatingAgents <- AlkylatingAgents[,..DatabaseName]
  AlkylatingAgents <- na.omit(AlkylatingAgents)
  
  for(i in 1:nrow(AlkylatingAgents))
  {
    #i <- 2
    limmaResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Limma Voom Analysis/%s/Table for %s DEG Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    pearsonResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Pearson Correlation Analysis/%s/Table for %s Pearson Correlation Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    spearmanExactResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Spearman (Exact) Correlation Analysis/%s/Table for %s Spearman (Exact) Correlation Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    #spearmanApproxResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Spearman (Approximate) Correlation Analysis/%s/Table for %s Spearman (Approximate) Correlation Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    DESeqResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/DESeq Analysis/%s/Table for %s DEG Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    EdgeRResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/EdgeR Analysis/%s/Table for %s DEG Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    
    limmaSG <- limmaResults[limmaResults$diffexpressed != "No difference"]$Gene
    pearsonSG <- pearsonResults[pearsonResults$diffexpressed != "No difference"]$Gene
    spearmanExactSG <- spearmanExactResults[spearmanExactResults$diffexpressed != "No difference"]$Gene
    #spearmanApproxSG <- spearmanApproxResults[spearmanApproxResults$diffexpressed != "No difference"]$Gene
    DESeqSG <- DESeqResults[DESeqResults$diffexpressed != "No difference"]$Gene
    EdgeRSG <- EdgeRResults[EdgeRResults$diffexpressed != "No difference"]$Gene

    SpearmanVennDiagram <- ggVennDiagram(list(limmaSG, DESeqSG, pearsonSG, spearmanExactSG),
                                  label_alpha = 0, set_color = "#084594", lwd = 0.8, lty = 1,
                                  label = "count",
                                  category.names = c("Limma",
                                                     "DESeq",
                                                     "Pearson",
                                                     "Spearman Exact")) +
       scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
       scale_x_continuous(expand = expansion(mult = .2)) +
       scale_color_brewer(palette = "Blues") +
       theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5))
    
     title <- paste0(sprintf(" %s Significant Genes for %s ", DatabaseName,toTitleCase(tolower(as.character(AlkylatingAgents[i])))))
     print(SpearmanVennDiagram +
             ggtitle(title))
     ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/Significant Gene Diagrams - Pearson/%s.png', directory, title)))

    EdgeRVennDiagram <- ggVennDiagram(list(limmaSG, DESeqSG, EdgeRSG, spearmanExactSG),
                                         label_alpha = 0, set_color = "#084594", lwd = 0.8, lty = 1, 
                                         label = "count", 
                                         category.names = c("Limma",
                                                            "DESeq",
                                                            "EdgeR",
                                                            "Spearman Exact")) + 
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      scale_x_continuous(expand = expansion(mult = .2)) +
      scale_color_brewer(palette = "Blues") +
      theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5))
    
    title <- paste0(sprintf(" %s Significant Genes for %s ", DatabaseName,toTitleCase(tolower(as.character(AlkylatingAgents[i])))))
    print(EdgeRVennDiagram + 
            ggtitle(title))
    ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/Significant Gene Diagrams - EdgeR/%s.png', directory, title)))
    
 
     MaxLength = max(length(limmaSG),length(limmaSG), length(spearmanExactSG), length(DESeqSG),length(EdgeRSG))
     significantGenes <- data.frame("Limma-Voom" = c(limmaSG, rep(NA, MaxLength - length(limmaSG))),
                                    "Pearson" = c(pearsonSG, rep(NA, MaxLength - length(pearsonSG))),
                                    "Spearman" = c(spearmanExactSG, rep(NA, MaxLength - length(spearmanExactSG))),
                                    "DESeq" = c(DESeqSG, rep(NA, MaxLength - length(DESeqSG))),
                                    "EdgeR" = c(EdgeRSG, rep(NA, MaxLength - length(EdgeRSG))))
     
     title <- paste0(sprintf("Significant %s Genes for %s",toTitleCase(tolower(as.character(AlkylatingAgents[i]))), DatabaseName)," DEG Analysis")
     write.csv(significantGenes, sprintf('%s/Isobel/CancerResistance/datasets/Combined Datasets/%s/%s.csv', directory, DatabaseName, title), row.names = FALSE)
  }
  
}

plotSignificantGenes("PRISM")
plotSignificantGenes("sangerGDSC")
plotSignificantGenes("ctd2")

combineResults <- function(DatabaseName)
{
  #DatabaseName <- "PRISM"
  AlkylatingAgents <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Combined Datasets/AlkylatingAgents.csv', directory)) 
  AlkylatingAgents <- AlkylatingAgents[,..DatabaseName]
  AlkylatingAgents <- na.omit(AlkylatingAgents)
  
  combinedDataset <- data.frame("AlkylatingAgent" = NA, "Gene" = NA, "LimmalogFC" = NA,
                                "LimmaFDR" = NA, "DESeqlogFC" = NA, "DESeqFDR" = NA,
                                "PearsonR" = NA, "PearsonFDR" = NA,
                                "SpearmanExactRho" = NA, "SpearmanExactFDR"= NA, 
                                "SpearmanApproxRho" = NA, "SpearmanApproxFDR"= NA)
  
  for(i in 1:nrow(AlkylatingAgents))
  {
    #i <- 2
    limmaResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Limma Voom Analysis/%s/Table for %s DEG Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    pearsonResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Pearson Correlation Analysis/%s/Table for %s Pearson Correlation Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    spearmanExactResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Spearman (Exact) Correlation Analysis/%s/Table for %s Spearman (Exact) Correlation Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    spearmanApproxResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Spearman (Approximate) Correlation Analysis/%s/Table for %s Spearman (Approximate) Correlation Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    DESeqResults <- fread(sprintf('%s/Isobel/CancerResistance/datasets/DESeq Analysis/%s/Table for %s DEG Analysis.csv', directory, DatabaseName, AlkylatingAgents[i])) 
    
    limmaResults <- limmaResults[order(limmaResults$Gene),]
    pearsonResults <- pearsonResults[order(pearsonResults$Gene),]
    spearmanExactResults <- spearmanExactResults[order(spearmanExactResults$Gene),]
    spearmanApproxResults <- spearmanApproxResults[order(spearmanApproxResults$Gene),]
    DESeqResults <- DESeqResults[order(DESeqResults$Gene),]
    
    #Cut genes not analyzed in all three methods - ignores duplicates
    includedGenes <- intersect(limmaResults$Gene, pearsonResults$Gene, spearmanExactResults$Gene, spearmanApproxResults$Gene, DESeqResults$Gene)
    
    limmaResults <- limmaResults[limmaResults$Gene %in% includedGenes]
    pearsonResults <- pearsonResults[pearsonResults$Gene %in% includedGenes]
    spearmanExactResults <- spearmanExactResults[spearmanExactResults$Gene %in% includedGenes]
    spearmanApproxResults <- spearmanApproxResults[spearmanApproxResults$Gene %in% includedGenes]
    DESeqResults <- DESeqResults[DESeqResults$Gene %in% includedGenes]
    
    
    limmaResults <- subset(limmaResults, !duplicated(Gene))
    pearsonResults <- subset(pearsonResults, !duplicated(Gene))
    spearmanExactResults <- subset(spearmanExactResults, !duplicated(Gene))
    spearmanApproxResults <- subset(spearmanApproxResults, !duplicated(Gene))
    DESeqResults <- subset(DESeqResults, !duplicated(Gene))
    
    
    newDataset  <- data.frame("AlkylatingAgent" =  rep(as.character(AlkylatingAgents[i]),nrow(pearsonResults)), 
                              "Gene" = pearsonResults$Gene, "LimmalogFC" = limmaResults$logFC, "LimmaFDR" = limmaResults$adj.P.Val, 
                              "DESeqlogFC" = DESeqResults$logFC, "DESeqFDR" = DESeqResults$FDR,
                              "PearsonR" = pearsonResults$Rvalue, "PearsonFDR" = pearsonResults$FDR,
                              "SpearmanExactRho" = spearmanExactResults$Rvalue, "SpearmanExactFDR" = spearmanExactResults$FDR,
                              "SpearmanApproxRho" = spearmanApproxResults$Rvalue, "SpearmanApproxFDR" = spearmanApproxResults$FDR)
    
    combinedDataset <- rbind(combinedDataset, newDataset)
  }
  write.csv(combinedDataset[-1,], sprintf('%s/Isobel/CancerResistance/datasets/Combined Datasets/%s Combined Dataset.csv', directory, DatabaseName), row.names = FALSE)  
}

combineResults("PRISM")
combineResults("sangerGDSC")
combineResults("ctd2")

plotExactParameter <- function(DatabaseName)
{
  combinedDataset <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Combined Datasets/%s Combined Dataset.csv', directory, DatabaseName))
  
  ExactParameterPlot <- 
    ggplot(data = combinedDataset, aes(x = SpearmanExactFDR, y = SpearmanApproxFDR)) + 
    geom_point() + theme
  print(ExactParameterPlot + 
          ggtitle(paste0(sprintf("%s Spearman FDR Values for Exact = T vs Exact = F", DatabaseName))) +
          labs(x= "Exact = T (Exact)", y = "Exact = F (Approximate)"))
  ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/Spearman Exact Parameter/%s Spearman FDR Values for Exact = T vs Exact = F.pdf', directory, DatabaseName)))
  
}

plotExactParameter("PRISM")
plotExactParameter("sangerGDSC")
plotExactParameter("ctd2")

plotFDR <- function(DatabaseName)
{
  #DatabaseName <- "PRISM"
  combinedDataset <- fread(sprintf('%s/Isobel/CancerResistance/datasets/Combined Datasets/%s Combined Dataset.csv', directory, DatabaseName))
  plotDataset <- data.frame("Limma" = combinedDataset$LimmaFDR, 
                            "SpearmanExact" = combinedDataset$SpearmanExactFDR, 
                            "DESeq" = combinedDataset$DESeqFDR)
  plotDataset$Limma <- -log10(plotDataset$Limma*sign(combinedDataset$LimmalogFC))
  plotDataset$SpearmanExact <- -log10(plotDataset$SpearmanExact*sign(combinedDataset$SpearmanExactRho))
  plotDataset$DESeq <- -log10(plotDataset$DESeq*sign(combinedDataset$DESeqlogFC))
  
  # plotDataset$Limma <- plotDataset$Limma*sign(combinedDataset$LimmalogFC)
  # plotDataset$SpearmanExact <- plotDataset$SpearmanExact*sign(combinedDataset$SpearmanExactRho)
  # plotDataset$DESeq <- plotDataset$DESeq*sign(combinedDataset$DESeqlogFC)
  
  SpearmanLimmaFDRplot <- 
    ggplot(data = plotDataset, aes(x = SpearmanExact, y = Limma)) + 
    geom_point() + theme
  
  title <- paste0(sprintf("%s Limma vs. Spearman Correlation FDRs", DatabaseName))
  print(SpearmanLimmaFDRplot + 
          ggtitle(title) +
          labs(x= "Spearman Exact -log(FDR)", y = "Limma -log(FDR)"))
  ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/Spearman vs. Limma FDR/%s.pdf', directory, title)))
  
  SpearmanDESeqFDRplot <- 
    ggplot(data = plotDataset, aes(x = SpearmanExact, y = DESeq)) + 
    geom_point() + theme
  
  title <- paste0(sprintf("%s DESeq vs. Spearman Correlation FDRs", DatabaseName))
  print(SpearmanDESeqFDRplot + 
          ggtitle(title) +
          labs(x= "Spearman Exact -log(FDR)", y = "DESeq -log(FDR)"))
  ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/Spearman vs. DESeq FDR/%s.pdf', directory, title)))
  
  LimmaDESeqFDRplot <- 
    ggplot(data = plotDataset, aes(x = Limma, y = DESeq)) + 
    geom_point() + theme
  
  title <- paste0(sprintf("%s Limma vs. DESeq FDRs", DatabaseName))
  print(LimmaDESeqFDRplot + 
          ggtitle(title) +
          labs(x= "Limma -log(FDR)", y = "DESeq -log(FDR)"))
  ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/Limma vs. DESeq FDR/%s.pdf', directory, title)))
  
}

plotFDR("PRISM")
plotFDR("sangerGDSC")
plotFDR("ctd2")

# install.packages("ggVennDiagram")
#install.packages('Rcpp')
install.packages('RVenn')




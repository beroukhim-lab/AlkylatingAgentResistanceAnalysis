#CCLE Analysis Rewrite w/ MGMT included

########## Load Packages ########## 
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(edgeR)

########## Load Data ##########
setwd( "/Users/igarrett/Desktop/R Code/ThesisRewrite ")

#From the CCLE 
RNAseqcounts <- fread("./data rewrite/CCLE_RNAseq_genes_counts_20180929 (1).gct") 
RNAseqcountsGlioma <- select(RNAseqcounts, Name, Description, contains("CENTRAL_NERVOUS_SYSTEM"))

#From the CTD^2 Data Portal 
Sensitivity <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt") 
CellLineMeta <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt") 
CompoundMeta <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt") 
ExperimentMeta <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt") 

########## Create EdgeRdataframe ##########
EdgeRdataframe <- data.frame(celllines = colnames(RNAseqcounts)[-(1:2)])
rownames(EdgeRdataframe) <- c(1:nrow(EdgeRdataframe))

Sensitivity$CellLineName <- NA

# Converts experiment_id to master_ccl_id (from ExperimentMeta) to cpd_name 
# (from CellLineMeta)
for(ExperimentID in unique(Sensitivity$experiment_id))
{
  if (length(which(ExperimentID == ExperimentMeta$experiment_id)) == 0)
  {
    cclID <- NA
  }
  
  else
  {
    cclID <- ExperimentMeta$master_ccl_id[ExperimentID == ExperimentMeta$experiment_id]
    ExperimentIDrows <- Sensitivity$experiment_id == ExperimentID
    Sensitivity$CellLineName[ExperimentIDrows] <- CellLineMeta$ccl_name[cclID == CellLineMeta$master_ccl_id]
  }
}

########## Filter out Alkylators ##########
DrugNametoID <- filter(CompoundMeta,grepl("DNA alkylator",target_or_activity_of_compound)) 
AlkylatorSensitivity <- filter(Sensitivity, master_cpd_id %in% DrugNametoID$master_cpd_id) 

########## Add Drug Sensitivity to EdgeRdataframe ##########
for (AlkylatorcpdID in unique(AlkylatorSensitivity$master_cpd_id))
{ 
  Colnumber <- ncol(EdgeRdataframe)
  CurrentAlkylatorSensitivityData <- subset(AlkylatorSensitivity, master_cpd_id == AlkylatorcpdID)
  AlkylatorSensitivityAverage <- aggregate(CurrentAlkylatorSensitivityData$area_under_curve,by = list(CellLineName = CurrentAlkylatorSensitivityData$CellLineName),data = CurrentAlkylatorSensitivityData,FUN = mean)
  
  for(RowNumber in (1:nrow(EdgeRdataframe)))
  { 
    CellLineName <- EdgeRdataframe$celllines[RowNumber] 
    CellLineNumber <- unlist(strsplit(as.character(CellLineName),"_"))[1]
    
    if (length(which(CellLineNumber == AlkylatorSensitivityAverage$CellLineName)) == 0)
    {
      EdgeRdataframe[RowNumber,Colnumber+1] <- NA
    }
    else
    {
      EdgeRdataframe[RowNumber,Colnumber+1] <- AlkylatorSensitivityAverage$x[which(CellLineNumber == AlkylatorSensitivityAverage$CellLineName)]
    } 
  }
  colnames(EdgeRdataframe)[Colnumber+1] <- paste(DrugNametoID$cpd_name[AlkylatorcpdID==DrugNametoID$master_cpd_id],'SensitivityAUC')
}

########## Filter table for only glioma cell lines ##########
EdgeRdataframe <- EdgeRdataframe[grep("CENTRAL_NERVOUS_SYSTEM",EdgeRdataframe$celllines,fixed=TRUE),]

########## Label cell lines as high/low sensitivity for alkylating agents ##########
GroupNumber <- 8
SensitivityStartIndex <- 2

for(Col in (SensitivityStartIndex:ncol(EdgeRdataframe)))
{
  MedianDrug = median(EdgeRdataframe[, Col], na.rm = TRUE)
  EdgeRdataframe[, Col] <- ntile(EdgeRdataframe[, Col], GroupNumber)
  
  for(Row in (1:nrow(EdgeRdataframe)))
  {
    #Ignores cells with NAs and those in the fourth and fifth octiles
    if (is.na(EdgeRdataframe[Row, Col])|| EdgeRdataframe[Row, Col] == GroupNumber/2
        || EdgeRdataframe[Row, Col] == GroupNumber/2 + 1)
    {
      EdgeRdataframe[Row, Col] = NA
    }
    else if (EdgeRdataframe[Row, Col] < MedianDrug)
    {
      EdgeRdataframe[Row, Col] = 'sensitive'
    }
    else
    {
      EdgeRdataframe[Row, Col] = 'resistant'
    }
  }
}

#Graph theme
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

#########     Making the Design Matrix    ######### 
ListofDesignMatrices <- c()
ListofCellLines <- c() 

for(Col in (SensitivityStartIndex:ncol(EdgeRdataframe)))
{
  Sensitivity <- c()
  CellLineNames <- c() 
  
  for(Row in (1:nrow(EdgeRdataframe)))
  { 
    if ( !is.na(EdgeRdataframe[Row, Col]))
    {
      new_sensitivity <- EdgeRdataframe[Row, Col]
      new_cellline <- as.character(EdgeRdataframe[Row, 1]) 
      
      Sensitivity <- c(Sensitivity, new_sensitivity)
      CellLineNames <- c(CellLineNames, new_cellline) 
    }
    
  }
  
  new_designMat <- model.matrix(~0+Sensitivity) 
  nextEmptySlot <- length(ListofDesignMatrices) + 1
  
  ListofDesignMatrices[[nextEmptySlot]] <- new_designMat
  ListofCellLines[[nextEmptySlot]] <- CellLineNames #list of cell line names for each drug in order of ListofDesignMatrices
  
} 

#Name the slots of the lists of design matrices by the drug they refer to
names(ListofDesignMatrices) <- as.vector(sapply(names(EdgeRdataframe)[SensitivityStartIndex:ncol(EdgeRdataframe)], function(x){substr(x, 1, regexpr(" ", x)-1)}))

#########     EdgeR Workflow    ######### 
dgListGliomaList <- c()
All_DEG <- c()

for(DesignMatrixIndex in 1:length(ListofDesignMatrices))
{
  dgListGlioma<- DGEList(counts=select(RNAseqcountsGlioma, unlist(ListofCellLines[DesignMatrixIndex])), 
                         genes=RNAseqcountsGlioma[,2]) 
  dgListGliomaList <- c(dgListGliomaList, dgListGlioma) 
  countsPerMillion <- cpm(dgListGlioma)
  
  #Filtering + Normalization 
  countCheck <- countsPerMillion > 1 
  keep <- which(rowSums(countCheck) >= 2) 
  dgListGlioma <- dgListGlioma[keep,]
  dgListGlioma <- calcNormFactors(dgListGlioma, method="TMM")
  
  # Plot 1: MDS
  Colors<-select(as.data.frame(ListofDesignMatrices[[DesignMatrixIndex]]),1)
  Colors[Colors=="1"]<-"blue" #resistant
  Colors[Colors=="0"]<-"red" #sensitive
  png(filename=paste0(sprintf("MDS plot for %s",names(ListofDesignMatrices)[DesignMatrixIndex]),".png"))
  MDSplots <- plotMDS(dgListGlioma,
                      pch = 20,
                      col=Colors[[1]],
                      ylim = c(-4, 4),
                      main=sprintf("MDS plot for %s",names(ListofDesignMatrices)[DesignMatrixIndex]))
  legend("topleft",c("Resistant", "Sensitive"), pch = c(16,16), col=c("blue","red"))
  dev.off()
  
  # Estimating Dispersons
  dgListGlioma <- estimateDisp(dgListGlioma, design=ListofDesignMatrices[[DesignMatrixIndex]])
  
  # Plot 2: BCV
  png(filename=paste0(sprintf("Dispersions for %s",names(ListofDesignMatrices)[DesignMatrixIndex]),".png"))
  plotBCV(dgListGlioma,
          main=sprintf("Dispersions for %s",names(ListofDesignMatrices)[DesignMatrixIndex]))
  dev.off()
  
  # Differential Expression
  fit <- glmFit(dgListGlioma, ListofDesignMatrices[[DesignMatrixIndex]]) 
  lrt <- glmLRT(fit, contrast=c(1,-1)) 
  edgeR_result <- topTags(lrt, n=Inf, p.value=.001)
  All_DEG <- c(All_DEG,list(edgeR_result$table))
  
  # Plot 3: Volcano
  lrt$table$diffexpressed <- "No difference"
  lrt$table$FDR <- NA
  lrt$table$FDR <- p.adjust(lrt$table$PValue,method="BH")
  lrt$table$diffexpressed[lrt$table$logFC > 1 & lrt$table$FDR < 0.01] <- "Up for resistant"
  lrt$table$diffexpressed[lrt$table$logFC < -1 & lrt$table$FDR < 0.01] <- "Down for resistant"
  
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("Down for resistant", "Up for resistant", "No difference")
  
  lrt$table$delabel <- NA
  lrt$table$delabel[lrt$table$diffexpressed != "No difference"] <- lrt$genes$Description[lrt$table$diffexpressed != "No difference"]

  lrt$table <- lrt$table[order(lrt$table$FDR),]
  MostSignificantFDR <- head(lrt$table[order(lrt$table$FDR),],10)
  MostSignificantGeneNames <- data.frame("delabel" = NA)
  
  Index <- 1
  for (Gene in lrt$table$delabel)
  {
    if(is.na(Gene) || !(Gene %in% MostSignificantFDR$delabel))
    {
      MostSignificantGeneNames[Index,1] <- NA
      
    }
    else
    {
      MostSignificantGeneNames[Index,1] <- Gene
    }
    Index <- Index + 1
  }
  

  volplot <- ggplot(data=lrt$table, mapping = aes(x=logFC, y=-log10(PValue), col=diffexpressed)) +
    geom_point() +
    theme + 
    geom_text_repel(aes(label = MostSignificantGeneNames$delabel), size = 3,
                    max.overlaps=Inf,
                    show.legend  = F) +
    scale_colour_manual(values = mycolors) +
    labs(title=sprintf("Differentially expressed genes for %s",names(ListofDesignMatrices)[DesignMatrixIndex]),
         y = expression("-log"[10]~"(P Value)"))
  ggsave(paste0(sprintf("Differentially expressed genes for %s",names(ListofDesignMatrices)[DesignMatrixIndex]),".png"))
  
  #Plot 4: PCA
  Flipped<-as.data.frame(t(dgListGlioma[["counts"]]))
  png(filename=paste0(sprintf("PCA plot for %s",names(ListofDesignMatrices)[DesignMatrixIndex]),".png"))
  
  autoplot(prcomp(Flipped, scale. = TRUE),
           labels= TRUE,
           colour=Colors[[1]],
           main=sprintf("PCA plot for %s",names(ListofDesignMatrices)[DesignMatrixIndex]))
  dev.off()
  
  for(Row in 1:nrow(lrt$table))
  {
    lrt$table$GeneName[Row] <- RNAseqcountsGlioma$Description[which(rownames(lrt$table)[Row]==rownames(RNAseqcountsGlioma))]
  }
  write.csv(lrt$table,paste0(sprintf("Table for %s",names(ListofDesignMatrices)[DesignMatrixIndex]),".csv"))
}

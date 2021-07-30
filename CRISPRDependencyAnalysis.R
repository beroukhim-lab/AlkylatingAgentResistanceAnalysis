#CRISPR CERES Score Analysis

########## Load Packages ########## 
library(data.table)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(rcompanion)
library(edgeR)

########## Load Data ##########
setwd( "/Users/igarrett/Desktop/R Code/ThesisRewrite ")

#From CCLE
RNAseqcounts <- fread("./data rewrite/CCLE_RNAseq_genes_counts_20180929 (1).gct") 
RNAseqcountsGlioma <- select(RNAseqcounts, Description, contains("CENTRAL_NERVOUS_SYSTEM"))

#From DepMap
CRISPRScores <- fread("./data rewrite/CRISPR_(DepMap_21Q2_Public+Score,_CERES).csv") 
CRISPRScores <- CRISPRScores[which(CRISPRScores$lineage_1 =="Central Nervous System"),]
CRISPRScores <- subset(CRISPRScores, select = -c(1,3,4,5,6))

#From the CTD^2 Data Portal 
Sensitivity <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt") 
CellLineMeta <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt") 
CompoundMeta <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt") 
ExperimentMeta <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt") 

Sensitivity$CellLineName <- NA

########## Create AnalysisDataframe ##########
AnalysisDataframe <- data.frame(celllines = CRISPRScores$cell_line_display_name)

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

########## Add Drug Sensitivity to AnalysisDataframe ##########
for (AlkylatorcpdID in unique(AlkylatorSensitivity$master_cpd_id))
{ 
  Colnumber <- ncol(AnalysisDataframe)
  CurrentAlkylatorSensitivityData <- subset(AlkylatorSensitivity, master_cpd_id == AlkylatorcpdID)
  AlkylatorSensitivityAverage <- aggregate(CurrentAlkylatorSensitivityData$area_under_curve,by = list(CellLineName = CurrentAlkylatorSensitivityData$CellLineName),data = CurrentAlkylatorSensitivityData,FUN = mean)
  
  for(RowNumber in (1:nrow(AnalysisDataframe)))
  { 
    CellLineName <- AnalysisDataframe$celllines[RowNumber] 
    CellLineNumber <- unlist(strsplit(as.character(CellLineName),"_"))[1]
    
    if (length(which(CellLineNumber == AlkylatorSensitivityAverage$CellLineName)) == 0)
    {
      AnalysisDataframe[RowNumber,Colnumber+1] <- NA
    }
    else
    {
      AnalysisDataframe[RowNumber,Colnumber+1] <- AlkylatorSensitivityAverage$x[which(CellLineNumber == AlkylatorSensitivityAverage$CellLineName)]
    } 
  }
  colnames(AnalysisDataframe)[Colnumber+1] <- paste(DrugNametoID$cpd_name[AlkylatorcpdID==DrugNametoID$master_cpd_id],'SensitivityAUC')
}

########## Label cell lines as high/low sensitivity for alkylating agents ##########
SensitivityStartIndex <- 2
SensitivityEndIndex <- 9
CRISPRScoreStartIndex <- 10

for(Col in (SensitivityStartIndex:ncol(AnalysisDataframe)))
{
  MedianDrug = median(AnalysisDataframe[, Col], na.rm = TRUE)
  AnalysisDataframe[, Col] <- ntile(AnalysisDataframe[, Col], GroupNumber)
  
  for(Row in (1:nrow(AnalysisDataframe)))
  {
    #Ignores cells with NAs and those in the fourth and fifth octiles
    if (is.na(AnalysisDataframe[Row, Col])|| AnalysisDataframe[Row, Col] == GroupNumber/2
        || AnalysisDataframe[Row, Col] == GroupNumber/2 + 1)
    {
      AnalysisDataframe[Row, Col] = NA
    }
    else if (AnalysisDataframe[Row, Col] < MedianDrug)
    {
      AnalysisDataframe[Row, Col] = 'sensitive'
    }
    else
    {
      AnalysisDataframe[Row, Col] = 'resistant'
    }
  }
}

#Add Gene data to AnalysisDataframe
AnalysisDataframe <- cbind(AnalysisDataframe, subset(CRISPRScores, select=-c(1)))

#Filter out cell lines 
ListofCellLineNumbers <- c()

for(Index in 2:ncol(RNAseqcountsGlioma))
{
  CellLineName <- colnames(RNAseqcountsGlioma)[Index] 
  CellLineNumber <- unlist(strsplit(as.character(CellLineName),"_"))[1]
  
  ListofCellLineNumbers <- c(ListofCellLineNumbers, CellLineNumber)
}
AnalysisDataframe <- AnalysisDataframe[AnalysisDataframe$celllines %in% ListofCellLineNumbers ,]

#Determine whether scores form a normal distribution
ShapiroTestResults <- data.frame("Genes" = colnames(AnalysisDataframe)[CRISPRScoreStartIndex:ncol(AnalysisDataframe)],"PValue" = NA)

counter = 1
for(GeneIndex in CRISPRScoreStartIndex:ncol(AnalysisDataframe))
{
  ShapiroTestResults[counter,2] <-shapiro.test(AnalysisDataframe[,GeneIndex])$p.value
  counter <- counter + 1
}

######Calculate P values and Effect Size #######
#Note:Multiplies effect size by - so that negative sign means lower CERES score

AlkylatingAgents <- c("chlorambucil", "dacarbazine", "ifosfamide", "temozolomide","bendamustine", "Platin" 
                      ,"cyclophosphamide", "oxaliplatin")

MWTestPValues <- data.frame(matrix(nrow = length(AlkylatingAgents), ncol = length(colnames(AnalysisDataframe)[10:ncol(AnalysisDataframe)]))) 
rownames(MWTestPValues) <- AlkylatingAgents
colnames(MWTestPValues) <- colnames(AnalysisDataframe)[10:ncol(AnalysisDataframe)]

EffectSize <- data.frame(matrix(nrow = length(AlkylatingAgents), ncol = length(colnames(AnalysisDataframe)[10:ncol(AnalysisDataframe)]))) 
rownames(EffectSize) <- AlkylatingAgents
colnames(EffectSize) <- colnames(AnalysisDataframe)[10:ncol(AnalysisDataframe)]

ResultsRow <- 1
for(SensitivityIndex in SensitivityStartIndex:SensitivityEndIndex) #For each alkylating agent
{
  ReducedDataframe <- AnalysisDataframe[which(!is.na(AnalysisDataframe[,SensitivityIndex])),]
 
  ResultsCol <- 1
  for(GeneIndex in CRISPRScoreStartIndex:ncol(AnalysisDataframe)) #For each gene
  {
    SensitivityLabels <- factor(ReducedDataframe[,SensitivityIndex])
    EffectSize[ResultsRow, ResultsCol] <- wilcoxonR(x = ReducedDataframe[,GeneIndex], g = SensitivityLabels )*-1
    LeveneTest <- leveneTest(ReducedDataframe[,GeneIndex]~SensitivityLabels, data = ReducedDataframe)$`Pr(>F)`[1]
    
    if(LeveneTest > 0.05)
    {
      MWTestPValues[ResultsRow, ResultsCol] <- wilcox.test(ReducedDataframe[,GeneIndex]~SensitivityLabels)$p.value
    }
    ResultsCol <- ResultsCol+1
  }
  
  ResultsRow <- ResultsRow+1
}

#FDR-adjust P Values
MWTestFDR <- apply(MWTestPValues,2,p.adjust,method = "BH")

#Calculate Ranking Metric: Sign of Effect Size * -log(P value)
Signs <- apply(EffectSize,2,sign)
RankingMetric <- Signs*-log(MWTestPValues,10)

####### Save CRISPR Dependencies for GSEA #########
for(Row in 1:nrow(RankingMetric)) #For each alkylating agent
{
  DependencyTable <- data.frame("Genes" = colnames(RankingMetric), "RankingMetric" = NA)
  GeneIndex <- 1
  
  for(Col in 1:ncol(RankingMetric)) #For each gene
  {
      DependencyTable[GeneIndex,2]<- RankingMetric[Row,Col]
      GeneIndex <- GeneIndex+1
  }
  DependencyTable<- DependencyTable[which(!(is.na(DependencyTable$RankingMetric))),]
  write.csv(DependencyTable,paste0(sprintf("CRISPR Ranked List for %s",rownames(RankingMetric)[Row]),".csv"), row.names = FALSE)
}

####### Create Volcano plots for CRISPR dependencies #######

#Theme for figures
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


mycolors <- c("blue", "red", "black")
names(mycolors) <- c("Down for resistant", "Up for resistant", "No difference")
PvalueIndex <- 1
FDRcutoff <- 0.01

for(AlkylatingAgentIndex in c(1:nrow(MWTestPValues)))
{
  VolcanoPlotTable <- data.frame("PValue" = t(MWTestPValues[AlkylatingAgentIndex,]), "FDR" = MWTestFDR[AlkylatingAgentIndex,],"EffectSize" = t(EffectSize[AlkylatingAgentIndex,]), "GeneName" = colnames(MWTestPValues))
  rownames(VolcanoPlotTable) <- c(1:nrow(VolcanoPlotTable))
  colnames(VolcanoPlotTable) <- c("PValue", "FDR", "EffectSize", "GeneName")
  
  VolcanoPlotTable$diff <- "No difference"
  VolcanoPlotTable$diff[VolcanoPlotTable$EffectSize > 0 & VolcanoPlotTable$FDR < FDRcutoff] <- "Up for resistant"
  VolcanoPlotTable$diff[VolcanoPlotTable$EffectSize < 0 & VolcanoPlotTable$FDR < FDRcutoff] <- "Down for resistant"

  #Search for most significant FDRs, will label the ten with largest effect size
  MostSignificantFDR <- subset(VolcanoPlotTable, FDR < FDRcutoff)
  MostSignificantFDR <- MostSignificantFDR[order(-abs(MostSignificantFDR$EffectSize)),]
  MostSignificantFDR <- head(MostSignificantFDR,10)
  MostSignificantGeneNames <- data.frame("GeneName" = NA)
  
  Index <- 1
  for (Gene in VolcanoPlotTable$GeneName)
  {
    if(is.na(Gene) || !(Gene %in% MostSignificantFDR$GeneName))
    {
      MostSignificantGeneNames[Index,1] <- NA
      
    }
    else
    {
      MostSignificantGeneNames[Index,1] <- Gene
    }
    Index <- Index + 1
  }
  
  volplot <- ggplot(data = VolcanoPlotTable, mapping = aes(x=EffectSize, y=-log10(PValue), col = diff)) +
    geom_point() + theme + 
    geom_text_repel(aes(label = MostSignificantGeneNames$GeneName), size = 3,
                    max.overlaps=Inf,show.legend  = F) +
    scale_colour_manual(values = mycolors) + geom_jitter(width = .2) +
    labs(title=sprintf("Differential CRISPR dependencies for %s",rownames(MWTestFDR)[AlkylatingAgentIndex]),
         x = "Effect Size", y = expression("-log"[10]~"(P Value)"))
  
  ggsave(paste0(sprintf("Differential CRISPR dependencies for %s",rownames(MWTestFDR)[AlkylatingAgentIndex]),".png"))
  
}

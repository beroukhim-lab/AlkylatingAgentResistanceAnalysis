#CCLE Analysis Rewrite

########## Load Packages ########## 
library(data.table)
library(car)
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

#From the CCLE 
RNAseqcounts <- fread("./data rewrite/CCLE_RNAseq_genes_counts_20180929 (1).gct") 
RNAseqcountsGlioma <- select(RNAseqcounts, Name, Description, contains("CENTRAL_NERVOUS_SYSTEM"))
Methylation <- fread("./data rewrite/CCLE_RRBS_TSS1kb_20181022.txt") 
RNAseqRSEM <- fread("./data rewrite/CCLE_RNAseq_rsem_genes_tpm_20180929.txt") 
MGMTgeneID <- RNAseqcounts[which(RNAseqcounts$Description =="MGMT"),]$Name

#From the CTD^2 Data Portal 
Sensitivity <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt") 
CellLineMeta <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt") 
CompoundMeta <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt") 
ExperimentMeta <- fread("./data rewrite/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt") 

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

########## Create EdgeRdataframe ##########
MGMTmethylation <- as.data.frame(t(Methylation[grep("MGMT", Methylation$locus_id),]))
MGMTmethylation <- MGMTmethylation[-c(1,2,3),,drop=FALSE]

MGMTexpression <- as.data.frame(t(RNAseqRSEM[which(RNAseqRSEM$gene_id == MGMTgeneID),]))
MGMTexpression <- MGMTexpression[-c(1,2),,drop=FALSE]

for(Row in (1: abs(nrow(MGMTexpression)-nrow(MGMTmethylation))))
{
  if(nrow(MGMTexpression) > nrow(MGMTmethylation))
  {
    MGMTmethylation <- rbind(MGMTmethylation, MGMTmethylation[NA,])
  }
  else
  {
    MGMTexpression <- rbind(MGMTexpression, MGMTexpression[NA,])
  }
}

EdgeRdataframe <- data.frame(celllines = colnames(RNAseqcounts)[-(1:2)])
EdgeRdataframe <- cbind(EdgeRdataframe, MGMTmethylation, MGMTexpression)
rownames(EdgeRdataframe) <- c(1:nrow(EdgeRdataframe))
colnames(EdgeRdataframe)[2] <- "MGMTmethylationValues"
colnames(EdgeRdataframe)[3] <- "MGMTexpressionValues"
EdgeRdataframe$MGMTexpressionLabels <- NA
EdgeRdataframe$MGMTmethylationLabels <- NA

EdgeRdataframe$MGMTmethylationValues <- as.numeric(EdgeRdataframe$MGMTmethylationValues)
EdgeRdataframe$MGMTexpressionValues <- as.numeric(EdgeRdataframe$MGMTexpressionValues)

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
plot(Sensitivity$conc_pts_fit,Sensitivity$area_under_curve) 

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

########## Label cell lines as high or low MGMT ##########
ggplot(data = EdgeRdataframe, aes(x = MGMTexpressionValues )) + geom_histogram(binwidth = 3) #200 

ExpressionMedian <- median(EdgeRdataframe$MGMTexpressionValues, na.rm = TRUE) # 0.12 for RSEM data
MethylationMedian <- median(EdgeRdataframe$MGMTmethylationValues, na.rm = TRUE) 

for(Row in (1:nrow(EdgeRdataframe)))
{
  if(!(is.na(EdgeRdataframe$MGMTexpressionValues[Row])))
  {
      if (EdgeRdataframe$MGMTexpressionValues[Row] < ExpressionMedian) 
      {
          EdgeRdataframe$MGMTexpressionLabels[Row] = 'low'
      }
      else
      {
          EdgeRdataframe$MGMTexpressionLabels[Row] = 'high'
      }
  }
  
  if(!(is.na(EdgeRdataframe$MGMTmethylationValues[Row])))
  {
      if (EdgeRdataframe$MGMTmethylationValues[Row] < MethylationMedian) 
      {
         EdgeRdataframe$MGMTmethylationLabels[Row] = 'low'
      }
      else
      {
         EdgeRdataframe$MGMTmethylationLabels[Row] = 'high'
      }
  }
}
length(which(EdgeRdataframe$MGMTexpressionLabels =='high')) #33 for high, 32 for low

#not normally distributed, need non-parametric t test
shapiro.test(EdgeRdataframe$MGMTexpressionValue) 
shapiro.test(EdgeRdataframe$MGMTmethylationValue) 

#MGMT methylation vs. expression Model
Model <- lm(EdgeRdataframe$MGMTexpressionValues~EdgeRdataframe$MGMTmethylationValues, data = EdgeRdataframe)
MGMTplotRsquared <- round(summary(Model)$r.squared,3)
MGMTplotPvalue <- round(summary(Model)$coefficients[2,4],3)

####### MGMT methylation and expression vs. drug sensitivity AUC Models #######
AlkylatingAgents <- c("Bendamustine", "Chlorambucil", "Dacarbazine", "Ifosfamide", "Platin", "Temozolomide"
                      ,"Cyclophosphamide", "Oxaliplatin")

ExpressionModelResults <- data.frame("AlkylatingAgents" = AlkylatingAgents, "AdjustedRSquared" = NA,
                                     "Pvalue" = NA, "LeveneTest" = NA, "WilcoxonTest" = NA, "EffectSize" = NA)
MethylationModelResults <- data.frame("AlkylatingAgents" = AlkylatingAgents,"AdjustedRSquared" = NA,
                                      "Pvalue" = NA, "LeveneTest" = NA, "WilcoxonTest" = NA, "EffectSize" = NA)

SensitivityStartIndex <- 6

ReducedExpressionDataframe <- EdgeRdataframe[which(!is.na(EdgeRdataframe$MGMTexpressionLabels)),]
ReducedMethylationDataframe <- EdgeRdataframe[which(!is.na(EdgeRdataframe$MGMTmethylationLabels)),]

for(Index in 1:length(ExpressionModelResults$AlkylatingAgents))
{
  ReducedExpressionDataframe <- ReducedExpressionDataframe[which(!is.na(ReducedExpressionDataframe[,Index+SensitivityStartIndex-1])),]
  ReducedMethylationDataframe <- ReducedMethylationDataframe[which(!is.na(  ReducedMethylationDataframe[,Index+SensitivityStartIndex-1])),]
  
  MGMTexpressionFactor <- factor(ReducedExpressionDataframe[,4])
  MGMTmethylationFactor <- factor(ReducedMethylationDataframe[,5])
  
  ExpressionModel <- lm(ReducedExpressionDataframe[,Index+SensitivityStartIndex-1]~MGMTexpressionFactor, data = ReducedExpressionDataframe)
  ExpressionModelResults$AdjustedRSquared[Index] <- summary(ExpressionModel)$r.squared
  ExpressionModelResults$Pvalue[Index] <- summary(ExpressionModel)$coefficients[2,4]
  
  MethylationModel <- lm(ReducedMethylationDataframe[,Index+SensitivityStartIndex-1]~MGMTmethylationFactor, data = ReducedMethylationDataframe)
  MethylationModelResults$AdjustedRSquared[Index] <- summary(MethylationModel)$r.squared
  MethylationModelResults$Pvalue[Index] <- summary(MethylationModel)$coefficients[2,4]
  
  #Tests if variance of the samples is the same - needed to perform Mann-Whitney Test  
  ExpressionModelResults$LeveneTest[Index] <- leveneTest(ReducedExpressionDataframe[,Index+SensitivityStartIndex-1]~MGMTexpressionFactor, data = ReducedExpressionDataframe)$`Pr(>F)`[1]
  MethylationModelResults$LeveneTest[Index] <- leveneTest(ReducedMethylationDataframe[,Index+SensitivityStartIndex-1]~MGMTmethylationFactor, data = ReducedMethylationDataframe)$`Pr(>F)`[1]
  
  # independent 2-group Mann-Whitney Test
  if( ExpressionModelResults$LeveneTest[Index] > 0.05) #Can perform MW Test
  {
    ExpressionModelResults$WilcoxonTest[Index] <- wilcox.test(ReducedExpressionDataframe[,Index+SensitivityStartIndex-1]~MGMTexpressionFactor)$p.value
    ExpressionModelResults$EffectSize[Index] <- wilcoxonR(x = ReducedExpressionDataframe[,Index+SensitivityStartIndex-1], g = MGMTexpressionFactor )
  }
  
  if( MethylationModelResults$LeveneTest[Index] > 0.05) #Can perform MW Test
  {
    MethylationModelResults$WilcoxonTest[Index] <- wilcox.test(ReducedMethylationDataframe[,Index+SensitivityStartIndex-1]~MGMTmethylationFactor)$p.value
    MethylationModelResults$EffectSize[Index] <-wilcoxonR(x = ReducedMethylationDataframe[,Index+SensitivityStartIndex-1], g = MGMTmethylationFactor )
  }
}

#MGMT vs. TMZ AUC
MGMTsensitivityTheme <- theme + 
  theme(panel.grid.major = element_line(colour = "azure3", size = 0.5)) +
  theme(panel.grid.minor = element_line(colour="azure3", size=0.5))

MGMTsensitivityPlot <- ggplot(data = EdgeRdataframe[which(!is.na(EdgeRdataframe$`temozolomide SensitivityAUC`)),], aes(x = MGMTexpressionLabels, y = `temozolomide SensitivityAUC`)) + 
  geom_point() + geom_jitter(width = .4) +
  annotate("text", label =  sprintf("paste(p, \" = %s\")",round(ExpressionModelResults$WilcoxonTest[match("Temozolomide",ExpressionModelResults$AlkylatingAgents)],2)), x = "low", y = 16.9, hjust = -1.3, vjust = 1,  parse = TRUE)+
  annotate("text", label =  sprintf("paste(d, \" = %s\")",round(ExpressionModelResults$EffectSize[match("Temozolomide",ExpressionModelResults$AlkylatingAgents)],2)), x = "low", y = 16.7, hjust = -1.2, vjust = 1,  parse = TRUE)
print(MGMTsensitivityPlot +
        ggtitle("MGMT Expression vs. TMZ Drug Sensitivity")+
        labs(y="TMZ Sensitivity (AUC)", x = "MGMT Gene Expression") + MGMTsensitivityTheme)

########## Label cell lines as high/low sensitivity for alkylating agents ##########
GroupNumber <- 8

for(Col in (SensitivityStartIndex:ncol(EdgeRdataframe)))
{
  MedianDrug = median(EdgeRdataframe[, Col], na.rm = TRUE)
  EdgeRdataframe[, Col] <- ntile(EdgeRdataframe[, Col], GroupNumber)
  
  for(Row in (1:nrow(EdgeRdataframe)))
  {
    #Ignores cells with NAs and those in the fourth and fifth octiles
    if (is.na(EdgeRdataframe[Row, Col]) || EdgeRdataframe[Row, Col] == GroupNumber/2
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

# MGMT methylation vs MGMT expression 
MGMTPlot <- ggplot(data = EdgeRdataframe, aes(x = MGMTmethylationValues, y = MGMTexpressionValues)) +
            annotate("text", x = 0.92, y = 19,
                     label = sprintf("paste(p, \" = %s\")",MGMTplotPvalue), parse = TRUE) +
            annotate("text", x = 0.92, y = 21,
                      label = sprintf("paste(R ^ 2, \" = %s\")",MGMTplotRsquared), parse = TRUE) +
            geom_point() + theme
print(MGMTPlot + 
        ggtitle("MGMT Methylation vs. Expression")+
        labs(y="MGMT Gene Expression (TPM)", x = "MGMT Gene Methylation")) + theme

######## Heatmap for MGMT T Tests ########
MGMTdata <- c("Expression", "Methylation")
Order <- c("Oxaliplatin","Platin","Dacarbazine", "Temozolomide", "Cyclophosphamide", 
                              "Ifosfamide", "Chlorambucil", "Bendamustine")
MGMTHeatmap <- expand.grid(MGMTdata =  MGMTdata, AlkylatingAgents = Order)
MGMTHeatmap <- MGMTHeatmap[order(MGMTHeatmap$MGMTdata),]

OrderedExpressionData <- ExpressionModelResults %>% slice(match(Order, AlkylatingAgents))
OrderedMethylationData <- MethylationModelResults %>% slice(match(Order, AlkylatingAgents))
AllEffectSizes <-c(OrderedExpressionData$EffectSize,OrderedMethylationData$EffectSize)
MGMTHeatmap$EffectSize <- AllEffectSizes

AllPvalues <-c(OrderedExpressionData$Pvalue,OrderedMethylationData$Pvalue)
MGMTHeatmap$Pvalue <- AllPvalues
MGMTHeatmap$Stars <- NA

for(Index in 1:nrow(MGMTHeatmap))
{
  if(MGMTHeatmap$Pvalue[Index] < 0.001)
  {
    MGMTHeatmap$Stars[Index] <- "***"
  }
  else if(MGMTHeatmap$Pvalue[Index] < 0.01)
  {
    MGMTHeatmap$Stars[Index] <- "**"
  }
  else if(MGMTHeatmap$Pvalue[Index] < 0.05)
  {
    MGMTHeatmap$Stars[Index] <- "*"
  }
}

HeatmapTheme <- theme + theme(panel.border = element_blank()) + theme(axis.ticks = element_blank()) 
MGMTHeatmapPlot <- ggplot(data = MGMTHeatmap, mapping = aes(x = MGMTdata, 
                y = AlkylatingAgents, fill = EffectSize)) + 
                scale_fill_gradient2(low = "brown3", mid="cornsilk1", high="turquoise4")+
                geom_tile(color = "black") +
                labs(fill="Effect Size") +
                geom_text(aes(label = Stars), color = "black", size=5) +
                HeatmapTheme

print(MGMTHeatmapPlot + 
        ggtitle("Effect Size for Mann-Whitney U Test")+
        labs(y="", x = "MGMT", fill = "Effect Size")) + HeatmapTheme

#########     Making the Design Matrix    ######### 
ListofDesignMatricesMGMThigh <- c()
ListofDesignMatricesMGMTlow <- c()
CellLinesMGMTHigh <- c() 
CellLinesMGMTLow <- c() 

for(Col in (SensitivityStartIndex:ncol(EdgeRdataframe)))
{
  MGMThighSensitivity <- c()
  MGMTlowSensitivity <- c()
  CellLineNamesMGMThigh <- c() 
  CellLineNamesMGMTlow <- c() 
  
  for(Row in (1:nrow(EdgeRdataframe)))
  { 
    if ( !is.na(EdgeRdataframe[Row, Col]))
    {
      new_sensitivity <- EdgeRdataframe[Row, Col]
      new_cellline <- as.character(EdgeRdataframe[Row, 1]) 
      
      if (EdgeRdataframe$MGMTexpressionLabels[Row]=='high')
      {
        MGMThighSensitivity <- c(MGMThighSensitivity, new_sensitivity)
        CellLineNamesMGMThigh <- c(CellLineNamesMGMThigh, new_cellline) 
      }
      else
      {
        MGMTlowSensitivity <- c(MGMTlowSensitivity, new_sensitivity)
        CellLineNamesMGMTlow <- c(CellLineNamesMGMTlow, new_cellline) 
      }
      
    }
  
  }

  new_designMatMGMThigh <- model.matrix(~0+MGMThighSensitivity) 
  new_designMatMGMTlow <- model.matrix(~0+MGMTlowSensitivity)
  
  nextEmptySlot <- length(ListofDesignMatricesMGMThigh) + 1
  
  ListofDesignMatricesMGMThigh[[nextEmptySlot]] <- new_designMatMGMThigh
  CellLinesMGMTHigh[[nextEmptySlot]] <- CellLineNamesMGMThigh #list of cell line names for each drug in order of ListofDesignMatricesMGMThigh
  ListofDesignMatricesMGMTlow[[nextEmptySlot]] <- new_designMatMGMTlow
  CellLinesMGMTLow[[nextEmptySlot]] <- CellLineNamesMGMTlow #list of cell line names for each drug in order of ListofDesignMatricesMGMTlow
} 

#Name the slots of the lists of design matrices by the drug they refer to
names(ListofDesignMatricesMGMThigh) <- as.vector(sapply(names(EdgeRdataframe)[SensitivityStartIndex:ncol(EdgeRdataframe)], function(x){substr(x, 1, regexpr(" ", x)-1)}))
names(ListofDesignMatricesMGMTlow) <- sapply(names(EdgeRdataframe)[SensitivityStartIndex:ncol(EdgeRdataframe)], function(x){substr(x, 1, regexpr(" ", x)-1)})

#########     EdgeR Workflow    ######### 
ListofAllDesignMatrices <- c(ListofDesignMatricesMGMThigh,ListofDesignMatricesMGMTlow)
ListofAllCellLines <- c(CellLinesMGMTHigh,CellLinesMGMTLow)
dgListGliomaList <- c()
All_DEG <- c()

for(DesignMatrixIndex in 1:length(ListofAllDesignMatrices))
{
  #Skip over NAs
  if (DesignMatrixIndex <= length(ListofAllDesignMatrices)/2)
  { 
    MGMTstatus <- "MGMT high"
  }
  else
  {
    MGMTstatus <- "MGMT low"
  }
  
  dgListGlioma<- DGEList(counts=select(RNAseqcountsGlioma, unlist(ListofAllCellLines[DesignMatrixIndex])), 
                         genes=RNAseqcountsGlioma[,2]) 
  dgListGliomaList <- c(dgListGliomaList, dgListGlioma) 
  countsPerMillion <- cpm(dgListGlioma)
  
  #Filtering + Normalization 
  countCheck <- countsPerMillion > 1 
  keep <- which(rowSums(countCheck) >= 2) 
  dgListGlioma <- dgListGlioma[keep,]
  dgListGlioma <- calcNormFactors(dgListGlioma, method="TMM")
  
  # Plot 1: MDS
  Colors<-select(as.data.frame(ListofAllDesignMatrices[[DesignMatrixIndex]]),1)
  Colors[Colors=="1"]<-"blue" #resistant
  Colors[Colors=="0"]<-"red" #sensitive
  png(filename=paste0(sprintf("MDS plot for %s, %s",names(ListofAllDesignMatrices)[DesignMatrixIndex],MGMTstatus),".png"))
  MDSplots <- plotMDS(dgListGlioma,
                      pch = 20,
                      col=Colors[[1]],
                      ylim = c(-4, 4),
                      main=sprintf("MDS plot for %s, %s",names(ListofAllDesignMatrices)[DesignMatrixIndex],MGMTstatus))
  legend("topleft",c("Resistant", "Sensitive"), pch = c(16,16), col=c("blue","red"))
  dev.off()
  
  # Estimating Dispersons
  dgListGlioma <- estimateDisp(dgListGlioma, design=ListofAllDesignMatrices[[DesignMatrixIndex]])
  
  # Plot 2: BCV
  png(filename=paste0(sprintf("Dispersions for %s, %s",names(ListofAllDesignMatrices)[DesignMatrixIndex],MGMTstatus),".png"))
  plotBCV(dgListGlioma,
          main=sprintf("Dispersions for %s, %s",names(ListofAllDesignMatrices)[DesignMatrixIndex],MGMTstatus))
  dev.off()
  
  # Differential Expression
  fit <- glmFit(dgListGlioma, ListofAllDesignMatrices[[DesignMatrixIndex]]) 
  lrt <- glmLRT(fit, contrast=c(1,-1)) 
  edgeR_result <- topTags(lrt, n=Inf, p.value=.001)
  All_DEG <- c(All_DEG,list(edgeR_result$table))
 
  # Plot 3: Volcano
  lrt$table$diffexpressed <- "No difference"
  lrt$table$FDR <- NA
  lrt$table$FDR <- p.adjust(lrt$table$PValue,method="BH")
  lrt$table$diffexpressed[lrt$table$logFC > 1 & lrt$table$FDR < 0.001] <- "Up for resistant"
  lrt$table$diffexpressed[lrt$table$logFC < -1 & lrt$table$FDR < 0.001] <- "Down for resistant"
 
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("Down for resistant", "Up for resistant", "No difference")
  
  lrt$table$delabel <- NA
  lrt$table$delabel[lrt$table$diffexpressed != "No difference"] <- lrt$genes$Description[lrt$table$diffexpressed != "No difference"]
  LastGene <- length(which(lrt$table$diffexpressed=="Down for resistant"))+length(which(lrt$table$diffexpressed=="Up for resistant"))
  Pvalue <- ((lrt$table[order(lrt$table$FDR), ][LastGene,4])+(lrt$table[order(lrt$table$FDR), ][LastGene+1,4]))/2
  
  volplot <- ggplot(data=lrt$table, mapping = aes(x=logFC, y=-log10(PValue), col=diffexpressed)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel(aes(label = lrt$table$delabel), size = 3,
                    max.overlaps=Inf,
                    show.legend  = F) +
    geom_hline(yintercept=-log10(Pvalue), col="red") +
    scale_colour_manual(values = mycolors) +
    labs(title=sprintf("Differentially expressed genes for %s, %s",names(ListofAllDesignMatrices)[DesignMatrixIndex],MGMTstatus))
  ggsave(paste0(sprintf("Differentially expressed genes for %s, %s",names(ListofAllDesignMatrices)[DesignMatrixIndex],MGMTstatus),".png"))
 
  #Plot 4: PCA
  Flipped<-as.data.frame(t(dgListGlioma[["counts"]]))
  png(filename=paste0(sprintf("PCA plot for %s, %s",names(ListofAllDesignMatrices)[DesignMatrixIndex],MGMTstatus),".png"))
  
  autoplot(prcomp(Flipped, scale. = TRUE),
           labels= TRUE,
           colour=Colors[[1]],
           main=sprintf("PCA plot for %s, %s",names(ListofAllDesignMatrices)[DesignMatrixIndex],MGMTstatus))
  dev.off()
  
  for(Row in 1:nrow(lrt$table))
  {
    lrt$table$GeneName[Row] <- RNAseqcountsGlioma$Description[which(rownames(lrt$table)[Row]==rownames(RNAseqcountsGlioma))]
  }
  write.csv(lrt$table,paste0(sprintf("Table for %s, %s",names(ListofAllDesignMatrices)[DesignMatrixIndex],MGMTstatus),".csv"))
}

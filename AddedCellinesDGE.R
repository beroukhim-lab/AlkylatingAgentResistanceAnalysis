#AddedCellinesDGE

########## Load Packages ########## 
library(data.table)
library(dplyr)
library(taigr)
library(ggplot2)
library(ggrepel)
library(edgeR)
library(tools)
library(DESeq2)
library(ggpubr)
library(edgeR)
library(statmod)

directory <- "/Volumes/xchip_beroukhimlab" #personal computer
#directory <-"/xchip/beroukhimlab" #server

######## Load Data #########

#Combined Cell Line Dataset (2012 CCLE, PRISM, Sanger GDSC, CTD2)
celllines <- fread(sprintf('%s/Isobel/CancerResistance/datasets/20211115_combined_drug_response.csv',directory)) 

AlkylatingAgentMoA <- c('DNA ALKYLATOR','INDUCER OF DNA DAMAGE', 
                        'DNA ALKYLATOR; ORGANOPLATINUM REAGENT', 'DNA SYNTHESIS INHIBITOR', 
                        'DNA INHIBITOR','DNA ALKYLATING AGENT, DNA INHIBITOR', 'DNA ALKYLATING AGENT', 
                        'DNA SYNTHESIS INHIBITOR','DNA ALKYLATING AGENT','CYTOCHROME P450 INHIBITOR')

AlkylatingAgentTarget <- c("DNA CROSSLINKER","DNA ALKYLATING AGENT", "ALKYLATING AGENT")

TISSUE_TYPE <- "CENTRAL_NERVOUS_SYSTEM"

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

#Permissive subsetting to make sure all alkylating agents are caught
subsetDatabase <- function(databaseName)
{
  database <-select(celllines, Cell_Line, Tissue_Type, Compound, AUC, Target, MoA, Dataset)
  database <- database[which(database$Tissue_Type == TISSUE_TYPE),]
  database <- database[which(database$Dataset == databaseName),]
  database <- database[order(Compound)]
  
  #Has no MoA data, but Target data sometimes includes MoA
  #Carmustine must be added manually, no target information 
  if(databaseName == "SANGER_GDSC") 
  {
    database  <- rbind(database[which(database$Target %in% AlkylatingAgentTarget),],
                       database[which(database$Compound == "CARMUSTINE")])
    # database  <- database[which(database$Target %in% AlkylatingAgentTarget),]
  }
  else
  {
    database  <- database[which(database$MoA %in% AlkylatingAgentMoA),]
  }
  return(database)
}

ccleNp24 <- subsetDatabase("CCLE_NP24")
prism <- subsetDatabase("PRISM_20Q2")
sangerGDSC <- subsetDatabase("SANGER_GDSC")
ctd2 <- subsetDatabase("CTD2")

#Exclude non-alkylating agents
#CB10-277 is related to dacarbazine
#TH-302 = Evofosfamide
ExtraneousCompounds <- c("BLEOMYCIN A2", "CLOFARABINE", "CYTARABINE HYDROCHLORIDE",
                         "ANISOMYCIN", "ARTEMISININ", "ARTESUNATE", "AZODICARBONAMIDE", 
                         "COBICISTAT","ENOCITABINE", "EPLERENONE", "FLOXURIDINE", "ITRACONAZOLE",
                         "METRONIDAZOLE", "METYRAPONE", "MORIN", "PF-4981517", "PIBENZIMOL",
                         "SANGIVAMYCIN", "TIRAPAZAMINE")
sangerGDSC <- sangerGDSC[!(sangerGDSC$Compound %in% ExtraneousCompounds)]
prism <- prism[!(prism$Compound %in% ExtraneousCompounds)]
ctd2 <- ctd2[!(ctd2$Compound %in% ExtraneousCompounds)]

#Alkylating Agents in each dataset
ccleNp24Drugs <-  unique(ccleNp24$Compound)
prismDrugs <-  unique(prism$Compound)
sangerGDSCDrugs <-  unique(sangerGDSC$Compound)
ctd2Drugs <- unique(ctd2$Compound)

MaxLength = max(length(ccleNp24Drugs),length(prismDrugs), length(sangerGDSCDrugs), length(ctd2Drugs))
AlkylatingAgentsPerDataset <- data.frame("ccleNp24" = c(ccleNp24Drugs, rep(NA, MaxLength - length(ccleNp24Drugs))), 
                                         "prism" = c(prismDrugs, rep(NA, MaxLength - length(prismDrugs))), 
                                         "sangerGDSC" = c(sangerGDSCDrugs, rep(NA, MaxLength - length(sangerGDSCDrugs))), 
                                         "ctd2" = c(ctd2Drugs, rep(NA, MaxLength - length(ctd2Drugs))))

# #No ccleNp24 drugs are alkylating agents
# AlkylatingAgentsPerDataset <- 
#   AlkylatingAgentsPerDataset[, !(colnames(AlkylatingAgentsPerDataset) %in% c("ccleNp24"))]

########## Create AnalysisDataframe ########## 
SensitivityStartIndex <- 2

createDataframe <- function(database)
{
  AnalysisDataframe <- dcast(data.table(database), Cell_Line ~ Compound, value.var = "AUC", fun.aggregate = mean)
  colnames(AnalysisDataframe)[1] <- "celllines"
  
  for(col in SensitivityStartIndex:ncol(AnalysisDataframe))
  {
    colnames(AnalysisDataframe)[col] <- paste(colnames(AnalysisDataframe)[col],'SensitivityAUC')
    AnalysisDataframe[,col] <- apply(AnalysisDataframe[,..col], 1, function(x){ifelse(is.na(x), return(NA), return(x))})
  }
  AnalysisDataframe <- as.data.frame(AnalysisDataframe)
  
  return (AnalysisDataframe)
}

prismDataframe <- createDataframe(prism)
sangerDataframe <- createDataframe(sangerGDSC)
ctd2Dataframe <- createDataframe(ctd2)

######## Data Filtering ########
#Want 10 or more celllines and some variation in AUC values for each Alkylating agent
dataCheck <- function(dataframe, dataframeName)
{
  CelllineNums <- c()
  AlkylatingAgents <- c()
  Range <- c()
  
  for(col in SensitivityStartIndex:ncol(dataframe))
  {
    CelllineNums <- c(CelllineNums,length(na.omit(dataframe[,col])))
    AlkylatingAgents <- c(AlkylatingAgents ,substr(colnames(dataframe)[col], start = 1, stop = 4))
    Range <- c(Range, max(dataframe[,col], na.rm=TRUE) - min(dataframe[,col], na.rm=TRUE))
  }
  
  dataCheckDataframe <- data.frame("Celllines" = CelllineNums, 
                                   "AlkylatingAgents" = AlkylatingAgents, "Range" = Range)
  
  CelllinePlot <- 
    ggplot(data = dataCheckDataframe, aes(x = AlkylatingAgents, y = Celllines)) + 
    geom_col() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme
  print(CelllinePlot + 
          ggtitle(paste0(sprintf("Celllines per Alkylating Agent for %s Dataset", dataframeName)))+
          labs(x= "Alkylating Agents", y = "Number of Celllines"))
  
  AUCrangePlot <- 
    ggplot(data = dataCheckDataframe, aes(x = AlkylatingAgents, y = Range)) + 
    geom_col() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme
  print(AUCrangePlot + 
          ggtitle(paste0(sprintf("AUC Ranges per Alkylating Agent for %s Dataset", dataframeName)))+
          labs(y="Range (AUC)", x = "Alkylating Agents"))
}

dataCheck(prismDataframe, "PRISM")
dataCheck(sangerDataframe, "sangerGDSC")
dataCheck(ctd2Dataframe, "ctd2")

filterLowVariance <- function(dataframe, dataframeName)
{
  dataframe <- prismDataframe
  Mean <- apply(dataframe[,-1], 2, function(x){mean(x, na.rm = TRUE)})
  Variation <- apply(dataframe[,-1], 2, function(x){var(x, na.rm = TRUE)})
  CelllineNum <- apply(dataframe[,-1], 2, function(x){sum(!is.na(x))})
  # DrugLabels <- unlist(strsplit(colnames(dataframe)[-1], " "))
  # DrugLabels <- DrugLabels[!DrugLabels %in% "SensitivityAUC"]
  
  avgVariance <- data.table(AlkylatingAgent = colnames(dataframe)[-1],
                            Mean = Mean,
                            Variation = Variation,
                            Celllines = CelllineNum)
  
  #Remove drugs with variance in the bottom 20% of variances for that mean
  numGroups <- ceiling(ncol(dataframe)/3)
  avgVariance <-  avgVariance[, meanIntoNTiles := ntile(Mean, numGroups)][ 
    , cutoffVarThisMeanNTile := quantile(Variation, 0.2, na.rm = TRUE), by = meanIntoNTiles]
  
  title <- sprintf("%s - %s groups", dataframeName,numGroups)
  
  MeanvsVariance <- 
    ggplot(avgVariance, aes(x = Mean, y = Variation, label = sub(" .*", "", AlkylatingAgent))) + 
    geom_point(alpha = 0.5) + 
    geom_smooth(aes(x = Mean, y = cutoffVarThisMeanNTile), color = "red", se = F, linetype = "longdash") +
    labs(title = title,
         x = "Mean of drug AUCs across celllines",
         y = "Variance of celllines\nwith available AUCs") +
    theme + geom_text_repel(size = 2)
  ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/MeanvsVariance/%s.pdf',directory,title)))
  
  filteredDrugs <- avgVariance[ Celllines >= 10,][
    Variation > cutoffVarThisMeanNTile, ]
  
  filteredAlkylatingAgents <- colnames(dataframe) %in% filteredDrugs$AlkylatingAgent
  dataframeFiltered <- cbind(dataframe[,1], dataframe %>% select(filteredDrugs$AlkylatingAgent))
  colnames(dataframeFiltered)[1] <- "celllines"
  
  #Remove rows of all NAs
  dataframeFiltered <- filter(dataframeFiltered, rowSums(is.na(dataframeFiltered)) != ncol(dataframeFiltered[,-1]))
  
  return (dataframeFiltered)
}

prismDataframe <- filterLowVariance(prismDataframe, "PRISM")
sangerDataframe <- filterLowVariance(sangerDataframe, "sangerGDSC")
ctd2Dataframe <- filterLowVariance(ctd2Dataframe, "ctd2")

remainingDrugs <- function(Dataframe, MaxLength)
{
  DrugLabels <- unlist(strsplit(colnames(Dataframe)[-1], " "))
  DrugLabels <- DrugLabels[!DrugLabels %in% "SensitivityAUC"]
  
  return (data.frame( "Col" = c(DrugLabels, rep(NA, MaxLength - length(DrugLabels)))))
}

MaxLength <- max(ncol(prismDataframe)-1, ncol(sangerDataframe)-1,ncol(ctd2Dataframe)-1, na.rm = T)
AlkylatingAgentsPerDataset <- cbind(remainingDrugs(prismDataframe, MaxLength), remainingDrugs(sangerDataframe, MaxLength),remainingDrugs(ctd2Dataframe, MaxLength))
names(AlkylatingAgentsPerDataset) <- c("PRISM", "sangerGDSC", "ctd2")
write.csv(AlkylatingAgentsPerDataset, sprintf('%s/Isobel/CancerResistance/datasets/Combined Datasets/AlkylatingAgents.csv', directory), row.names = FALSE)  


prismReadCounts <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=21, data.file='public_20Q2_counts')
ccleReadCounts <- load.from.taiga(data.name='ccle-rnaseq-reads-count-per-gene', data.version=1, data.file='data')

secondaryScreenCelllineInfo <- load.from.taiga(data.name='prism-repurposing-20q2-5e60', data.version=5, data.file='secondary-screen-cell-line-info')
secondaryScreenCelllineInfo <- subset(secondaryScreenCelllineInfo, !is.na(secondaryScreenCelllineInfo$depmap_id))

#Passed str profiling?
secondaryScreenCelllineInfo <- secondaryScreenCelllineInfo[secondaryScreenCelllineInfo$passed_str_profiling == T,]

#Still duplicated celllines?
#View(subset(secondaryScreenCelllineInfo, duplicated(depmap_id))) #95 entries
#View(secondaryScreenCelllineInfo[secondaryScreenCelllineInfo$depmap_id == "ACH-000096",])
#View(secondaryScreenCelllineInfo[secondaryScreenCelllineInfo$depmap_id == "ACH-000082",])
#Duplicated entries are exactly the same, can just delete one of the duplicates
secondaryScreenCelllineInfo <- subset(secondaryScreenCelllineInfo, !duplicated(depmap_id))

prismCelllineNames <- c()
for(row in 1:nrow(prismReadCounts))
{
  CelllineName <- secondaryScreenCelllineInfo$ccle_name[which(secondaryScreenCelllineInfo$depmap_id == rownames(prismReadCounts)[row])]
  if(length(CelllineName) != 0)
  {
    newCelllineName <- sub("[_].*", "", CelllineName)
    prismCelllineNames <- c(prismCelllineNames, newCelllineName)
  }
  else
  {
    prismCelllineNames <- c(prismCelllineNames, "NA")
  }
}

rownames(prismReadCounts) <- unlist(prismCelllineNames)
rownames(ccleReadCounts) <- sub("[_].*", "", rownames(ccleReadCounts))

filterCounts <- function(Dataframe, ReadCounts)
{
  colnames(ReadCounts) <- sub(" .*", "", colnames(ReadCounts))
  ReadCounts <- subset(ReadCounts, rownames(ReadCounts) %in%  Dataframe$celllines)
  ReadCounts <- ReadCounts[order(rownames(ReadCounts)),]
  
  return (ReadCounts) 
} 

prismReadCounts <- filterCounts(prismDataframe, prismReadCounts)
sangerReadCounts <- filterCounts(sangerDataframe, prismReadCounts)
ctd2ReadCounts <- filterCounts(ctd2Dataframe, prismReadCounts)

prismDataframe <- subset(prismDataframe, prismDataframe$celllines %in%  rownames(prismReadCounts))
sangerDataframe <- subset(sangerDataframe, sangerDataframe$celllines %in%  rownames(sangerReadCounts))
ctd2Dataframe <- subset(ctd2Dataframe, ctd2Dataframe$celllines %in%  rownames(ctd2ReadCounts))

######   Limma DEG Analysis   ######

# prismDGList <- DGEList(counts=t(prismReadCounts))
# sangerDGList <- DGEList(counts=t(sangerReadCounts))
# ctd2DGList <- DGEList(counts=t(ctd2ReadCounts))

#For design process demo
#write.csv(prismDataframe, sprintf('%s/Isobel/CancerResistance/datasets/prismDataframe.csv', directory), row.names = FALSE)

##### Create Design Matrix #####
makeDesignMatrix <- function(Dataframe)
{
  designMatrix <- Dataframe[,-1, drop=FALSE]
  designMatrix <- as.matrix(designMatrix)
  
  #Save global option 
  current.na.action <- options("na.action") 
  options(na.action = "na.pass") 
  
  designMatrix <- model.matrix(~designMatrix)#designMatrix <- model.matrix(~0+designMatrix)
  rownames(designMatrix) <- Dataframe[,1]
  colnames(designMatrix)[-1] <- colnames(Dataframe[-c(1)]) #  colnames(designMatrix) <- colnames(Dataframe[-c(1)])
  options("na.action" = current.na.action$na.action)
  
  return (designMatrix)
}

prismDesignMatrix <- makeDesignMatrix(prismDataframe)
sangerDesignMatrix <- makeDesignMatrix(sangerDataframe)
ctd2DesignMatrix <- makeDesignMatrix(ctd2Dataframe)

######## Filtering out low counts ######### 
filterLowCountsLimma <- function(DGList, ReadCounts, ThresholdValue)
{
  CPM <- cpm(DGList) 
  
  plot(CPM[,2],DGList$counts[,2],ylim=c(0,30),xlim=c(0,0.8))
  abline(v=ThresholdValue)
  
  #Good threshold is ~0.3 prism, 0.2 sanger, 0.2 ctd2, corresponding to 10-15 counts
  threshold <- CPM > ThresholdValue
  #Need to be "significant" CPM in at least a quarter of the celllines
  keep <- rowSums(threshold) >= nrow(ReadCounts)/4
  DGList <- DGList[keep,,keep.lib.sizes=FALSE]
  
  #Test if major discrepancies in library sizes for each celline
  barplot(DGList$samples$lib.size/1e06, names=colnames(DGList), las=2, ann=FALSE, cex.names=0.75)
  mtext(side = 2, text = "Library size (millions)", line = 3)
  title("Barplot of Library Sizes")
  #A few cell lines have much higher library sizes (100+ million) but most have ~70 million
  
  #TMM Normalization
  #Accounts for library size variation between samples
  return (calcNormFactors(DGList))
}

#Voom: Using Precision Weights - Useful if library sizes are quite varied
#Counts close to 0 or/and "beak shape" mean data not properly filtered
runVoomAnalysis <- function(DesignMatrix, ReadCounts, DatabaseName)
{
  for(row in SensitivityStartIndex:ncol(DesignMatrix))#for(row in 1:ncol(DesignMatrix))
  {
    # row <- 6
    # DesignMatrix <- prismDesignMatrix
    # ReadCounts <- prismReadCounts
    
    ReducedDesignMatrix <- na.omit(DesignMatrix[,c(1,row), drop = FALSE])
    ReducedReadCounts <- ReadCounts[which(rownames(ReadCounts) %in%  rownames(ReducedDesignMatrix)),]
    
    ReduceddgList <- DGEList(counts=t(ReducedReadCounts))
    #ReduceddgList <- filterLowCountsLimma(ReduceddgList, ReducedReadCounts, 0.2)
    
    keep <- filterByExpr(ReduceddgList, ReducedDesignMatrix)
    ReduceddgList <- ReduceddgList[keep,,keep.lib.sizes=FALSE]
    ReduceddgList <- calcNormFactors(ReduceddgList)
    
    v <- voom(ReduceddgList, ReducedDesignMatrix, plot=F)
    
    #Fit the linear model - estimates group means and gene variances
    fit <- lmFit(v,  ReducedDesignMatrix, method = "robust") 
    #fit <- lmFit(v,  ReducedDesignMatrix) 
    
    #eBayes - performs empirical Bayes shrinkage on the variances
    #Estimates moderated t-statistics and the associated p-values
    fit <- eBayes(fit) #fit <- eBayes(fit, robust = T)
    
    #topTable(fit, coef = 2, p.value = 0.05)
    #summary(decideTests(fit))
    DEGresults <- topTable(fit, 2, n = Inf)#DEGresults <- topTable(fit, coef=ncol(ReducedDesignMatrix), n = Inf)#DEGresults <- topTable(fit, coef=row, n = Inf) 
    colnames(DEGresults)[1] <- "Gene"
    testSummary <- decideTests(fit) 
    #summary(testSummary)
    
    # #P-value Distribution
    # ggplot(data = limmaResults, aes(x=P.Value))+ geom_histogram(binwidth=0.01)+
    #   theme
    # ggplot(data = DESeqResults, aes(x=Pvalue))+ geom_histogram(binwidth=0.01)+
    #   theme
    
    AlkylatingAgentName <- unlist(strsplit(colnames(ReducedDesignMatrix), " "))#AlkylatingAgentName <- unlist(strsplit(colnames(DesignMatrix), " "))
    AlkylatingAgentName <-AlkylatingAgentName[!AlkylatingAgentName %in% c("(Intercept)", "SensitivityAUC")]#AlkylatingAgentName <- AlkylatingAgentName[!AlkylatingAgentName %in% "SensitivityAUC"][row]
    
    #Volcano Plot
    DEGresults$diffexpressed <- "No difference"
    DEGresults$diffexpressed[DEGresults$logFC > 1 & DEGresults$adj.P.Val < 0.05] <- "Down for resistant"
    DEGresults$diffexpressed[DEGresults$logFC < -1 & DEGresults$adj.P.Val < 0.05] <- "Up for resistant"
    
    mycolors <- c("blue", "red", "black")
    names(mycolors) <- c("Down for resistant", "Up for resistant", "No difference")
    
    DEGresults$delabel <- NA
    DEGresults$delabel[DEGresults$diffexpressed != "No difference"] <- DEGresults$Gene[DEGresults$diffexpressed != "No difference"]
    
    TopTen <- topTable(fit, 2, n = 10, sort.by = "p")#TopTen <- topTable(fit, coef=row, n = 10, sort.by = "p")#TopTen <- topTable(fit, coef=ncol(ReducedDesignMatrix), n = 10)
    TopTen$diffexpressed <- "No difference"
    TopTen$diffexpressed[TopTen$logFC > 1 & TopTen$adj.P.Val < 0.05] <- "Down for resistant"
    TopTen$diffexpressed[TopTen$logFC < -1 & TopTen$adj.P.Val < 0.05] <- "Up for resistant"
    TopTen <- TopTen[TopTen$diffexpressed != "No difference",]
    
    ggplot(data=DEGresults, mapping = aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) +
      geom_point() +
      theme +
      geom_text_repel(data =TopTen, aes(label = ID), size = 3,
                      max.overlaps=Inf,
                      show.legend  = F) +
      scale_colour_manual(values = mycolors) +
      labs(title=sprintf("Differentially Expressed Genes For %s",AlkylatingAgentName ))
    
    filename <- paste0(sprintf("Differentially Expressed Genes For %s",AlkylatingAgentName))
    ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/Limma Voom Analysis/%s/%s.pdf', directory, DatabaseName, filename)))
    
    title <- paste0(sprintf("Table for %s",AlkylatingAgentName)," DEG Analysis")
    write.csv(DEGresults, sprintf('%s/Isobel/CancerResistance/datasets/Limma Voom Analysis/%s/%s.csv', directory, DatabaseName, title), row.names = FALSE)
  }
}

runVoomAnalysis(prismDesignMatrix, prismReadCounts, "PRISM")
runVoomAnalysis(sangerDesignMatrix, sangerReadCounts, "sangerGDSC")
runVoomAnalysis(ctd2DesignMatrix, ctd2ReadCounts, "ctd2")

####### logCPM Analysis ######
filterLowCountsCPM <- function(ReadCounts, ThresholdValue)
{
  transposedCounts <- t(ReadCounts)
  CPM <- cpm(transposedCounts) 
  
  #Check to see what threshold CPM value corresponds to ~10-15 counts
  plot(CPM[,1],transposedCounts[,1],ylim=c(0,30),xlim=c(0,0.8))
  abline(v=ThresholdValue)
  
  threshold <- CPM > ThresholdValue
  #Need to be "significant" CPM in at least a quarter of the celllines
  keep <- rowSums(threshold) >= nrow(ReadCounts)/4
  transposedCounts <- transposedCounts[keep,]
  
  #Still created a dgList item because cpm(dgList) uses the stored normalization factors internally
  #Slightly different results when cpm(matrix of counts) used
  dgList <- DGEList(counts=transposedCounts) 
  normFactors <- calcNormFactors(dgList)
  logCPM <- cpm(dgList, log = TRUE)
  
  return (t(logCPM)) #More easily aligns with AUC data
  
}

# ###### Normalization Tests #########
# shapiro.test(logCPM [0:5000])$p.value #Not normally distributed
# 
# #Anderson-Darling normality test
# library(nortest)
# ad.test(logCPM)$p.value #Not normally distributed
# 
# #Use non-parametric correlation tests
runCorrelationAnalysis <- function(Dataframe, ReadCounts, Type, Exact, DatabaseName)
{
  for(col in SensitivityStartIndex:ncol(Dataframe))
  {
    # col <- 2
    # Dataframe <- prismDataframe
    # Type <- "spearman"
    # Exact <- T
    
    SensitivityAUC <- subset(Dataframe,select = c(1,col)) #22-Thiotepa
    SensitivityAUC <- na.omit(SensitivityAUC)
    SensitivityAUC <- SensitivityAUC[order(SensitivityAUC$celllines),]
    
    ReducedCPM <- subset(ReadCounts, rownames(ReadCounts) %in%  SensitivityAUC$celllines)
    ReducedCPM <- ReducedCPM[order(rownames(ReducedCPM)),]
    ReducedCPM <- filterLowCountsCPM(ReducedCPM, 0.3)
    
    GeneNames <- c(colnames(ReducedCPM))
    CPMCorrelations <- data.frame("Gene" = GeneNames, "Rvalue" = NA, "Tstatistic" = NA,
                                  "Pvalue" = NA)
    
    for(row in 1:nrow(CPMCorrelations))
    {
      #row <- 1
      Gene <- as.character(CPMCorrelations$Gene[row])
      plotDataset <- cbind(subset(SensitivityAUC,select = c(2)), ReducedCPM[,Gene, drop=FALSE]) 
      Result <- cor.test(plotDataset[,1],plotDataset[,2], method = Type, exact = Exact)
      CPMCorrelations$Rvalue[row] <- Result$estimate
      CPMCorrelations$Pvalue[row] <- Result$p.value
      CPMCorrelations$Tstatistic[row] <- Result$statistic
    }
    OrderedCPMCorrelations <- CPMCorrelations[order(CPMCorrelations$Rvalue),]
    OrderedCPMCorrelations <- na.omit(OrderedCPMCorrelations) #NAs because for columnns of 0s run into 0 STD error
    OrderedCPMCorrelations$FDR <- p.adjust(OrderedCPMCorrelations$Pvalue, method="BH")
    OrderedCPMCorrelations$diffexpressed <- "No difference"
    
    #Resistant for positive correlation, Sensitive for negative correlation
    OrderedCPMCorrelations$diffexpressed[OrderedCPMCorrelations$Rvalue > 0.5 & OrderedCPMCorrelations$FDR < 0.05] <- "Down for resistant"
    OrderedCPMCorrelations$diffexpressed[OrderedCPMCorrelations$Rvalue < -0.5 & OrderedCPMCorrelations$FDR < 0.05] <- "Up for resistant"
    
    # CPMSignificantResults <- OrderedCPMCorrelations[abs(OrderedCPMCorrelations$Rvalue) > 0.5,]
    
    AlkylatingAgentName <- unlist(strsplit(colnames(SensitivityAUC)[2], " "))
    AlkylatingAgentName <- AlkylatingAgentName[!AlkylatingAgentName %in% "SensitivityAUC"]
    
    # mycolors <- c("blue", "red", "black")
    # names(mycolors) <- c("Down for resistant", "Up for resistant", "No difference")
    # 
    # ggplot(data=OrderedCPMCorrelations, mapping = aes(x=Rvalue, y=-log10(Pvalue), col=diffexpressed)) +
    #     geom_point() +
    #     theme +
    #     geom_text_repel(data =CPMSignificantResults, aes(label = Gene), size = 3,
    #                     max.overlaps=Inf,
    #                     show.legend  = F) +
    #     scale_colour_manual(values = mycolors) +
    #     labs(title=sprintf("Differentially Expressed Genes For %s",AlkylatingAgentName ))
    # 
    # filename <- paste0(sprintf("Differentially Expressed Genes For %s",AlkylatingAgentName))
    # ggsave(file = paste0(sprintf('/Volumes/xchip_beroukhimlab/Isobel/CancerResistance/figures/CPM Correlation Analysis/%s.pdf',filename)))
    
    TypeLabel <- Type
    if(Type == "spearman" && Exact == F)
    {
      TypeLabel <- paste(Type, "(Approximate)", sep=" ")
    }
    
    if(Type == "spearman" && Exact == T)
    {
      TypeLabel <- paste(Type, "(Exact)", sep=" ")
    }
    
    title <- paste0(sprintf("Table for %s %s Correlation Analysis",AlkylatingAgentName, toTitleCase(TypeLabel)))
    write.csv(OrderedCPMCorrelations, sprintf('%s/Isobel/CancerResistance/datasets/%s Correlation Analysis/%s/%s.csv', directory, toTitleCase(TypeLabel), DatabaseName, title), row.names = FALSE)  
    
  }
}

runCorrelationAnalysis(prismDataframe, prismReadCounts, "pearson", T, "PRISM")
runCorrelationAnalysis(sangerDataframe, sangerReadCounts, "pearson", T, "sangerGDSC")
runCorrelationAnalysis(ctd2Dataframe, ctd2ReadCounts, "pearson", T, "ctd2")

runCorrelationAnalysis(prismDataframe, prismReadCounts, "spearman", T, "PRISM")
runCorrelationAnalysis(sangerDataframe, sangerReadCounts, "spearman", T, "sangerGDSC")
runCorrelationAnalysis(ctd2Dataframe, ctd2ReadCounts, "spearman", T, "ctd2")

runCorrelationAnalysis(prismDataframe, prismReadCounts, "spearman", F, "PRISM")
runCorrelationAnalysis(sangerDataframe, sangerReadCounts, "spearman", F, "sangerGDSC")
runCorrelationAnalysis(ctd2Dataframe, ctd2ReadCounts, "spearman", F, "ctd2")

###### Example Plots ########
GeneName <- "GREM1"
AlkylatingAgentName <- "Ifosfamide"
Dataframe <- ctd2Dataframe
DatabaseName <- "ctd2"
ReadCounts <- ctd2ReadCounts

AlkylatingAgentNum <- as.numeric(which(colnames(Dataframe) == sprintf("%s SensitivityAUC", toupper(AlkylatingAgentName))))
AUC <- subset(Dataframe,select = c(1,AlkylatingAgentNum)) 
AUC <- na.omit(AUC)
AUC <- AUC[order(AUC$celllines),]

logCPM <- subset(ReadCounts, rownames(ReadCounts) %in%  AUC$celllines)
logCPM <- logCPM [order(rownames(logCPM)),]
logCPM <- filterLowCountsCPM(logCPM, 0.3)

#corPlot <- data.frame("AUC" = AUC[,2], "GALNT4" = logCPM[,GeneName[1]], "TX1" = logCPM[,GeneName[2]])
corPlot <- data.frame("AUC" = AUC[,2], "logCPM" = logCPM[,GeneName])
equation <- corPlot$logCPM ~ corPlot$AUC

#Need adjusted p-value and Limma data
limmaResults <- fread(sprintf('/Volumes/xchip_beroukhimlab/Isobel/CancerResistance/datasets/Limma Voom Analysis/%s/Table for %s DEG Analysis.csv', DatabaseName, toupper(AlkylatingAgentName)))
spearmanResults <- fread(sprintf('/Volumes/xchip_beroukhimlab/Isobel/CancerResistance/datasets/Spearman (Exact) Correlation Analysis/%s/Table for %s Spearman (Exact) Correlation Analysis.csv', DatabaseName, toupper(AlkylatingAgentName))) 
DESeqResults <- fread(sprintf('/Volumes/xchip_beroukhimlab/Isobel/CancerResistance/datasets/DESeq Analysis/%s/Table for %s DEG Analysis.csv', DatabaseName, toupper(AlkylatingAgentName))) 
EdgeRResults <- fread(sprintf('/Volumes/xchip_beroukhimlab/Isobel/CancerResistance/datasets/EdgeR Analysis/%s/Table for %s DEG Analysis.csv', DatabaseName, toupper(AlkylatingAgentName))) 

RvalueLabel <- sprintf("R = %s",round(spearmanResults$Rvalue[spearmanResults$Gene == GeneName],2))
cpmFDRLabel <- sprintf("p = %s",round(spearmanResults$FDR[spearmanResults$Gene == GeneName],3))
limmalogFCLabel <- sprintf("Limma logFC = %s",round(limmaResults$logFC[limmaResults$Gene == GeneName],2))
limmaFDRLabel <- sprintf("p = %s",round(limmaResults$adj.P.Val[limmaResults$Gene == GeneName],6))
DESeqlogFCLabel  <- sprintf("DESeq logFC = %s",round(DESeqResults$logFC[DESeqResults$Gene == GeneName],2))
DESeqFDRLabel <- sprintf("p = %s",round(DESeqResults$FDR[DESeqResults$Gene == GeneName],14))
EdgeRlogFCLabel  <- sprintf("EdgeR logFC = %s",round(EdgeRResults$logFC[EdgeRResults$Gene == GeneName],2))
EdgeRFDRLabel <- sprintf("p = %s",round(EdgeRResults$FDR[EdgeRResults$Gene == GeneName],13))

#May need to adjust (x,y) of coordinates per graph
annotations <- data.frame(
  xpos = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf, -Inf,-Inf),
  ypos =  c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
  annotateText = c(RvalueLabel, cpmFDRLabel, limmalogFCLabel, limmaFDRLabel,DESeqlogFCLabel,DESeqFDRLabel,EdgeRlogFCLabel,EdgeRFDRLabel),
  hjustvar = c(-0.19,-0.36,-0.09,-0.20, -0.09, -0.20,-0.09, -0.20) ,
  vjustvar = c(3,4.8, 8, 9.8, 13, 14.8, 18, 19.8)) 

logCPMplot <- 
  ggplot(data = corPlot, aes(x = AUC, y = logCPM)) + geom_point() + 
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText, parse = TRUE)) +
  #geom_smooth(method = "lm", se=F) +
  #stat_regline_equation(label.y = 1.5, aes(label = ..eq.label..)) + 
  theme 
print(logCPMplot +
        ggtitle(sprintf("%s Sensitivity and %s Gene Expression Correlation", AlkylatingAgentName, GeneName))+
        labs(y=sprintf("%s Expression (logCPM)", GeneName), x = sprintf("%s (AUC)", AlkylatingAgentName)))

cor(corPlot$AUC,corPlot$logCPM, method = "pearson")
cor.test(corPlot$AUC,corPlot$logCPM, method = "spearman", Exact = F)
########### DESeq2 Analysis #############

filterLowCountsDESeq <- function(ReadCounts, ThresholdValue)
{
  transposedCounts <- t(ReadCounts)
  CPM <- cpm(transposedCounts) 
  
  #Check to see what threshold CPM value corresponds to ~10-15 counts
  plot(CPM[,1],transposedCounts[,1],ylim=c(0,30),xlim=c(0,0.8))
  abline(v=ThresholdValue)
  
  threshold <- CPM > ThresholdValue
  #Need to be "significant" CPM in at least a quarter of the celllines
  keep <- rowSums(threshold) >= nrow(ReadCounts)/4
  transposedCounts <- transposedCounts[keep,]
  
  return (transposedCounts) 
  
}

RunDESeqAnalysis <- function(ReadCounts, DesignMatrix, DatabaseName)
{
  
  for(row in SensitivityStartIndex:ncol(DesignMatrix))
  {
    # row <- 6
    # ReadCounts <- prismReadCounts
    # DesignMatrix <- prismDesignMatrix
    
    ReducedDesignMatrix <- na.omit(DesignMatrix[,c(1,row), drop = FALSE])
    
    AlkylatingAgentName <- unlist(strsplit(colnames(ReducedDesignMatrix), " "))
    AlkylatingAgentName <- AlkylatingAgentName[!AlkylatingAgentName %in% c("(Intercept)","SensitivityAUC")]
    
    colnames(ReducedDesignMatrix)[-1] <- "AUC"
    ReducedReadCounts <- ReadCounts[which(rownames(ReadCounts) %in%  rownames(ReducedDesignMatrix)),]
    #ReducedReadCounts <- t(ReducedReadCounts)
    ReducedReadCounts <- filterLowCountsDESeq(ReducedReadCounts, 0.2) #duplicate genes only show up in this? ReducedReadCounts <- t(ReducedReadCounts)

    #Model corrects for library size
    DESeqModel <- DESeqDataSetFromMatrix(countData = round(ReducedReadCounts),
                                         colData = ReducedDesignMatrix,
                                         design = ~AUC) 
    
    # DESeqModel <- estimateSizeFactors(DESeqModel)
    # nc <- counts(DESeqModel, normalized=TRUE)
    # filter <- rowSums(nc >= 10) >= 2
    # DESeqModel <- DESeqModel[filter,]


    DESeqObject <- DESeq(DESeqModel)
    
    #How many adjusted p-values are less than 0.05
    #sum(DESeqResults$padj < 0.05, na.rm=TRUE) #101
    
    #DESeqResults <- results(DESeqObject, alpha=0.05, lfcThreshold = 1)
    DESeqResults <- results(DESeqObject, name = "AUC")
    #DESeqResults <- lfcShrink(DESeqObject, coef=2, type = "apeglm")#results(DESeqObject)
    #DESeqResults <- results(DESeqObject)
    #summary(DESeqResults)
    #sum(DESeqResults$padj < 0.05, na.rm=TRUE) 
    
    DESeqDataframe <- data.frame("Gene" = DESeqResults@rownames, 
                                 "logFC" = DESeqResults$log2FoldChange, 
                                 "Pvalue" = DESeqResults$pvalue, "FDR" = DESeqResults$padj)
    
    #Volcano Plot
    DESeqDataframe$diffexpressed <- "No difference"
    DESeqDataframe$diffexpressed[DESeqDataframe$logFC > 1 & DESeqDataframe$FDR < 0.05] <- "Down for resistant"
    DESeqDataframe$diffexpressed[DESeqDataframe$logFC < -1 & DESeqDataframe$FDR < 0.05] <- "Up for resistant"
    
    MostSignificantGenes <- DESeqDataframe[DESeqDataframe$diffexpressed != "No difference",]
    MostSignificantGenes <- MostSignificantGenes[order(MostSignificantGenes$FDR),]
    MostSignificantGenes <- top_n(MostSignificantGenes, -10, FDR) #Only label most significant
    
    mycolors <- c("blue", "red", "black")
    names(mycolors) <- c("Down for resistant", "Up for resistant", "No difference")
    
    ggplot(data = DESeqDataframe, mapping = aes(x=logFC, y=-log10(Pvalue), col=diffexpressed)) +
      geom_point() +
      theme +
      geom_text_repel(data = MostSignificantGenes, aes(label = Gene), size = 3,
                      max.overlaps=Inf,
                      show.legend  = F) +
      scale_colour_manual(values = mycolors) +
      labs(title=sprintf("Differentially Expressed Genes For %s",AlkylatingAgentName ))
    
    filename <- paste0(sprintf("Differentially Expressed Genes For %s",AlkylatingAgentName))
    ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/DESeq Analysis/%s/%s.pdf', directory, DatabaseName, filename)))
    
    title <- paste0(sprintf("Table for %s",AlkylatingAgentName)," DEG Analysis")
    write.csv(DESeqDataframe, sprintf('%s/Isobel/CancerResistance/datasets/DESeq Analysis/%s/%s.csv', directory, DatabaseName, title), row.names = FALSE)
    
  }
  
}

RunDESeqAnalysis(prismReadCounts, prismDesignMatrix, "PRISM")
RunDESeqAnalysis(sangerReadCounts, sangerDesignMatrix, "sangerGDSC")
RunDESeqAnalysis(ctd2ReadCounts, ctd2DesignMatrix, "ctd2")

###### EdgeR ########

runEdgeRAnalysis <- function(DesignMatrix, ReadCounts, DatabaseName)
{
  All_DEG <- c()
  
  for(row in SensitivityStartIndex:ncol(DesignMatrix))
  {
    # row <- 6
    # DesignMatrix <- prismDesignMatrix
    # ReadCounts <- prismReadCounts
    
    ReducedDesignMatrix <- na.omit(DesignMatrix[,c(1,row), drop = FALSE])
    ReducedReadCounts <- ReadCounts[which(rownames(ReadCounts) %in%  rownames(ReducedDesignMatrix)),]
    
    ReduceddgList <- DGEList(counts=t(ReducedReadCounts))
    #ReduceddgList <- filterLowCountsLimma(ReduceddgList, ReducedReadCounts, 0.2)
    
    #Filtering + Normalization 
    CPM <- cpm(ReduceddgList)
    countCheck <- CPM > 1 
    keep <- which(rowSums(countCheck) >= 2) 
    ReduceddgList <- ReduceddgList[keep,]
    ReduceddgList <- calcNormFactors(ReduceddgList, method="TMM")
    
    ReduceddgList <- estimateDisp( ReduceddgList, design = ReducedDesignMatrix, robust = TRUE)
    #ReduceddgList <- estimateDisp( ReduceddgList, design = ReducedDesignMatrix[,-1])
    
    fit <- glmFit(ReduceddgList , ReducedDesignMatrix) 
    lrt <- glmLRT(fit, 2) 
    DEGresults <- topTags(lrt, n=Inf) # DEGresults <- topTags(lrt, n=Inf, p.value=.001)
    All_DEG <- c(All_DEG,list(DEGresults$table))
    
    AlkylatingAgentName <- unlist(strsplit(colnames(ReducedDesignMatrix), " "))
    AlkylatingAgentName <-AlkylatingAgentName[!AlkylatingAgentName %in% c("(Intercept)", "SensitivityAUC")]
    
    #Volcano Plot
    lrt$table$FDR <- p.adjust(lrt$table$PValue,method="BH")
    lrt$table$diffexpressed <- "No difference"
    lrt$table$diffexpressed[lrt$table$logFC > 1 & lrt$table$FDR < 0.05] <- "Down for resistant"
    lrt$table$diffexpressed[lrt$table$logFC < -1 & lrt$table$FDR < 0.05] <- "Up for resistant"
    #View(lrt$table[lrt$table$diffexpressed != "No difference",])
    
    mycolors <- c("blue", "red", "black")
    names(mycolors) <- c("Down for resistant", "Up for resistant", "No difference")
    
    lrt$table <- lrt$table[order(lrt$table$FDR),]
    lrt$table$Gene <- rownames(lrt$table)
    rownames(lrt$table) <- seq(1:nrow(lrt$table))
    MostSignificantGenes <- lrt$table[lrt$table$diffexpressed != "No difference",]
    MostSignificantGenes <- top_n(MostSignificantGenes, -10, FDR) #Only label most significant

    ggplot(data=lrt$table, mapping = aes(x=logFC, y=-log10(PValue), col=diffexpressed)) +
      geom_point() +
      theme +
      geom_text_repel(data = MostSignificantGenes, aes(label = Gene), size = 3,
                      max.overlaps=Inf,
                      show.legend  = F) +
      scale_colour_manual(values = mycolors) +
      labs(title=sprintf("Differentially Expressed Genes For %s",AlkylatingAgentName ))
    
    filename <- paste0(sprintf("Differentially Expressed Genes For %s",AlkylatingAgentName))
    ggsave(file = paste0(sprintf('%s/Isobel/CancerResistance/figures/EdgeR Analysis/%s/%s.pdf', directory, DatabaseName, filename)))
    
    title <- paste0(sprintf("Table for %s",AlkylatingAgentName)," DEG Analysis")
    write.csv(lrt$table, sprintf('%s/Isobel/CancerResistance/datasets/EdgeR Analysis/%s/%s.csv', directory, DatabaseName, title), row.names = FALSE)
    
   }
}

runEdgeRAnalysis(prismDesignMatrix, prismReadCounts, "PRISM")
runEdgeRAnalysis(sangerDesignMatrix, sangerReadCounts, "sangerGDSC")
runEdgeRAnalysis(ctd2DesignMatrix, ctd2ReadCounts, "ctd2")



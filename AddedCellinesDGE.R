#AddedCellinesDGE

########## Load Packages ########## 
library(data.table)
library(dplyr)
library(taigr)
library(ggplot2)
library(edgeR)

######## Load Data #########
setwd( "/Users/isobelgarrett/Desktop/College/BSRP 2021/R Code/ThesisRewrite ")

#Combined Cell Line Dataset (2012 CCLE, PRISM, Sanger GDSC, CTD2)
celllines <- fread('/Volumes/xchip_beroukhimlab/Isobel/CancerResistance/datasets/20211115_combined_drug_response.csv') 

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

########## Create AnalysisDataframe ########## 
SensitivityStartIndex <- 2

createDataframe <- function(database, AlkylatingAgentsPerDatasetColumn)
{
    AnalysisDataframe <- data.frame("celllines" = unique(database$Cell_Line))

    for(row in (1:nrow(AlkylatingAgentsPerDataset)))
    {
      if(! is.na(AlkylatingAgentsPerDatasetColumn[row]))
      {
          erCol <- row+1
          AnalysisDataframe[,erCol] <- NA
          colnames(AnalysisDataframe)[erCol] <- 
             paste(AlkylatingAgentsPerDatasetColumn[row],'SensitivityAUC') 
      }
    }
    
    for(col in (SensitivityStartIndex:ncol(AnalysisDataframe)))
    {
      
      for(databaseRow in (1:nrow(database)))
      {
        if(paste(database$Compound[databaseRow],'SensitivityAUC') == colnames(AnalysisDataframe)[col])
        {
          erRow = which(AnalysisDataframe$celllines == database$Cell_Line[databaseRow]) 
          AnalysisDataframe[erRow,col] <- database$AUC[databaseRow]
        }
      }
    }
    
    return (AnalysisDataframe)
}

prismDataframe <- createDataframe(prism, AlkylatingAgentsPerDataset$prism)
sangerDataframe <- createDataframe(sangerGDSC, AlkylatingAgentsPerDataset$sangerGDSC)
ctd2Dataframe <- createDataframe(ctd2, AlkylatingAgentsPerDataset$ctd2)


plotTMZ <- function(dataframe, dataframeName)
{
    TemozolomideSensitivityPlot <- 
      ggplot(data = dataframe, aes(x = dataframe$`TEMOZOLOMIDE SensitivityAUC`)) + 
      geom_histogram(binwidth = 0.05) + theme
    print(TemozolomideSensitivityPlot + 
            ggtitle(paste0(sprintf("Temozolomide Sensitivity for %s Dataset", dataframeName)))+
            labs(y="Number of Celllines", x = "Sensitivity (AUC)"))
}

plotTMZ(prismDataframe, "PRISM")
plotTMZ(sangerDataframe, "sangerGDSC")
plotTMZ(ctd2Dataframe, "ctd2")

#########     Making the Design Matrix (for now just for PRISM)    ######### 
DesignMatrix <- prismDataframe[,-c(1), drop=FALSE]
DesignMatrix <- as.matrix(DesignMatrix)
DesignMatrix <- model.matrix.lm(~0+DesignMatrix, na.action = "na.pass")
rownames(DesignMatrix) <- prismDataframe[,1]

#DesignMatrix <- model.matrix(~0+DesignMatrix) #model.matrix by default removes rows with NAs
#DesignMatrix <- model.frame(~0+DesignMatrix,na.action='na.pass')
#DesignMatrix <- model.matrix(DesignMatrix )

######   Limma DEG Analysis   ######

#Limma-Voom: Take read counts, convert to logCPM, model mean-variance relationship using 
#precision weights(voom) (DONT use Bayes, that won't be using voom)
# 1. Need Matrix of Read Counts
#"If you are working with RSEM gene-level expected counts, then you can just pass them to limma as if they were counts"

#PRISMreadCounts from RSEM
#Rows: Broad (Arxspan) cell line IDs
#Columns: genes (HGNC symbol and Ensembl ID)
#"read count data is created using RSEM"
prismReadCounts <- load.from.taiga(data.name='depmap-rnaseq-expression-data-ccd0', data.version=21, data.file='public_20Q2_counts')
dgList <- DGEList(counts=prismReadCounts) #Warning: Library size of zero detected
#Library size: total number of reads that were sequenced in the run or the total number of mapped reads 

#Filter out 0 or low counts
#Can filter using filterByExpr or cpm
#For cpm need prior.count or else log(0) will cause an error
#filterByExpr - "Determine which genes have sufficiently large counts to be retained"
logCPM <- cpm(dgList, log=TRUE, prior.count = 0.3) #ERROR:library sizes should be finite and non-negative
keep <- filterByExpr(dgList, DesignMatrix) #ERROR: NAs in design matrix
keep <- which(rowSums(logCPM) >= 2) 
dgList <- dgeList[keep,,keep.lib.sizes=FALSE]

#Why is there an error when using cpm()?
#any(is.na(dgList$samples$lib.size)) --> FALSE
#any(dgList$samples$lib.size < 0) --> FALSE
#any(dgList$samples$lib.size == 0) --> TRUE
#range(dgList$samples$lib.size) --> 0 1575851541
#Libraries of size 0 will cause an error, so use prior.count --> Currently not working
#"Prior.count = average count to be added to each observation to avoid taking log of zero"

#TMM Normalization
dgList <- calcNormFactors(dgList)

#Voom: Using Precision Weights
v <- voom(dgList, DesignMatrix, plot=TRUE)
Fit <- lmFit(v, DesignMatrix)



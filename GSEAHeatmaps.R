#Creates Heatmaps from GSEA Results

library(ggplot2)
library(ggfortify)
library(ggrepel)

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

AlkylatingAgents = c("Oxaliplatin","Platin","Dacarbazine", "Temozolomide", "Cyclophosphamide", 
                     "Ifosfamide", "Chlorambucil", "Bendamustine") 

#Hardcode NES scores and P values
NESscores <- data.frame(AlkylatingAgents = AlkylatingAgents)
NESscores$DNArepairDGE <- c(-2.56,1.13,-1.46,-1.60,2.10,1.79,-2.29,-1.05)
NESscores$ApoptosisDGE <- c(3.09,1.17,1.15,-1.03,1.76,1.04,2.02,0.86)

NESscores$DNArepairCRISPR <- c(-1.59,-2.50,1.33,-1.40,-0.95,1.59,-1.11,-1.23)
NESscores$ApoptosisCRISPR <- c(1.71,-0.72,1.57,-0.76,0.71,-0.77,1.26,-1.09)
DataType <- c("DGE", "Dependencies")

Pvalues <- data.frame(AlkylatingAgents = AlkylatingAgents)
Pvalues$DNArepairDGE <- c(8.02e-5,0.30,0.081,0.04,0.003,0.014,0.0016,0.38)
Pvalues$ApoptosisDGE <- c(0,0.26,0.27,0.42,0.02,0.39,0.005, 0.61)
Pvalues$DNArepairCRISPR <- c(0.047,3.2e-4,0.14,0.11,0.52,0.047,0.33,0.20)
Pvalues$ApoptosisCRISPR <- c(0.02,0.81,0.066,0.78,0.87,0.74,0.17,0.35)

#################### DNA Repair Heatmap #######################
DNArepairHeatmap <- expand.grid(DataType =  DataType, AlkylatingAgents = AlkylatingAgents)
DNArepairHeatmap <- DNArepairHeatmap[order(DNArepairHeatmap$DataType),]
AllNES <-c(NESscores$DNArepairDGE,NESscores$DNArepairCRISPR)
DNArepairHeatmap$NES <- AllNES

DNArepairHeatmap$Pvalue <- c(Pvalues$DNArepairDGE,Pvalues$DNArepairCRISPR)
DNArepairHeatmap$Stars <- NA

for(Index in 1:nrow(DNArepairHeatmap))
{
  if(DNArepairHeatmap$Pvalue[Index] < 0.001)
  {
    DNArepairHeatmap$Stars[Index] <- "***"
  }
  else if(DNArepairHeatmap$Pvalue[Index] < 0.01)
  {
    DNArepairHeatmap$Stars[Index] <- "**"
  }
  else if(DNArepairHeatmap$Pvalue[Index] < 0.05)
  {
    DNArepairHeatmap$Stars[Index] <- "*"
  }
}

DNArepairHeatmapTheme  <- theme + theme(panel.border = element_blank()) + 
  theme(axis.ticks = element_blank()) + theme(axis.title.x = element_text(margin = margin(t = 10)))#theme(legend.text = element_text(size=12),legend.key.size = unit(0.5, "in")) 
DNArepairHeatmapPlot <- ggplot(data = DNArepairHeatmap, mapping = aes(x = DataType, y = AlkylatingAgents, fill = NES)) + 
  scale_fill_gradient2(low = "brown3", mid="cornsilk1", high="turquoise4",limits = c(-4,4)) +
  geom_tile(color = "black") + geom_text(aes(label = Stars), color = "black", size=5) + 
  DNArepairHeatmapTheme

print(DNArepairHeatmapPlot + 
        ggtitle("DNA Repair GSEA Results")+
        labs(y="", x = "Differential Data")) + DNArepairHeatmapTheme

#################### Apoptosis Heatmap #######################
ApoptosisHeatmap <- expand.grid(DataType =  DataType, AlkylatingAgents = AlkylatingAgents)
ApoptosisHeatmap <- ApoptosisHeatmap[order(ApoptosisHeatmap$DataType),]
AllNES <-c(NESscores$ApoptosisDGE,NESscores$ApoptosisCRISPR)
ApoptosisHeatmap$NES <- AllNES

AllPvalues <-c(Pvalues$ApoptosisDGE,Pvalues$ApoptosisCRISPR)
ApoptosisHeatmap$Pvalue <- AllPvalues
ApoptosisHeatmap$Stars <- NA

for(Index in 1:nrow(ApoptosisHeatmap))
{
  if(ApoptosisHeatmap$Pvalue[Index] < 0.001)
  {
    ApoptosisHeatmap$Stars[Index] <- "***"
  }
  else if(ApoptosisHeatmap$Pvalue[Index] < 0.01)
  {
    ApoptosisHeatmap$Stars[Index] <- "**"
  }
  else if(ApoptosisHeatmap$Pvalue[Index] < 0.05)
  {
    ApoptosisHeatmap$Stars[Index] <- "*"
  }
}

ApoptosisHeatmapTheme  <- theme + theme(panel.border = element_blank()) + 
  theme(axis.ticks = element_blank())+ theme(axis.title.x = element_text(margin = margin(t = 10)))# + theme(legend.text = element_text(size=12),legend.key.size = unit(0.5, "in")) 
ApoptosisHeatmapPlot <- ggplot(data = ApoptosisHeatmap, mapping = aes(x = DataType, y = AlkylatingAgents, fill = NES)) + 
  scale_fill_gradient2(low = "brown3", mid="cornsilk1", high="turquoise4",limits = c(-4,4)) +
  geom_tile(color = "black") + geom_text(aes(label = Stars), color = "black", size=5) + 
  ApoptosisHeatmapTheme

print(ApoptosisHeatmapPlot + 
        ggtitle("Apoptosis GSEA Results")+
        labs(y="", x = "Differential Data")) + ApoptosisHeatmapTheme

library(data.table)

setwd( "/Users/igarrett/Desktop/Thesis Rewrite")

Names <- c("bendamustine", "chlorambucil", "dacarbazine", "ifosfamide", "Platin", "temozolomide"
          ,"cyclophosphamide", "oxaliplatin")

for(Index in 1:length(Names))
{
  LowMGMT <- fread(sprintf("./Table for %s, MGMT low.csv", Names[Index])) 
  HighMGMT <- fread(sprintf("./Table for %s, MGMT high.csv", Names[Index]))
  
  #3 = logCPM, 2 = logFC, 5 = PValue
  RankingMetric = 2
  
  LowMGMT <- subset(LowMGMT, select = c(9,RankingMetric))
  HighMGMT <- subset(HighMGMT, select = c(9,RankingMetric))
  
  write.csv(LowMGMT, sprintf("./RankedListLowMGMTFC%s.csv", Names[Index]), row.names = FALSE)
  write.csv(HighMGMT, sprintf("./RankedListHighMGMTFC%s.csv", Names[Index]), row.names = FALSE)
}

write.csv(ExpressionList, "./ExpressionListAll.csv", row.names = FALSE)
write.csv(Sensitivity, "./SensitivityAll.csv")

# BendamustineLowMGMT <- fread("./Table for bendamustine, MGMT low.csv")
# BendamustineHighMGMT <- fread("./Table for bendamustine, MGMT high.csv")
# 
# BendamustineLowMGMT <- BendamustineLowMGMT[which(BendamustineLowMGMT$diffexpressed !="No difference"),]
# BendamustineHighMGMT <- BendamustineHighMGMT[which(BendamustineHighMGMT$diffexpressed !="No difference"),]

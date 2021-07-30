library(data.table)

setwd( "/Users/igarrett/Desktop/Thesis Rewrite")

Names <- c("bendamustine", "chlorambucil", "dacarbazine", "ifosfamide", "Platin", "temozolomide"
          ,"cyclophosphamide", "oxaliplatin")

#Format Ranked List w/ Ranking Metric: Sign of logFC * -log(P value)
for(Index in 1:length(Names))
{
  GeneList <- fread(sprintf("./Table for %s.csv", Names[Index]))
  GeneList$PValue <- -log(GeneList$PValue,10)
  GeneList$PValue <- GeneList$PValue * sign(GeneList$logFC)
  
  #3 = logCPM, 2 = logFC, 5 = PValue, 9 = GeneName
  GeneList <- subset(GeneList, select = c(9,5))
  
  write.csv(GeneList, sprintf("./RankedListFinal%s.csv", Names[Index]), row.names = FALSE)

}
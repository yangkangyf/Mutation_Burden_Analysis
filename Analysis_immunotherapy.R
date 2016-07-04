#########################
# Adjustable parameters #
#########################
# Choose which data to analyse
Van_Allen = T
Snyder = T
Rizvi = T
TCGA = T

########################
# IMPACT IMMUNOTHERAPY #
########################
if (TCGA == T)
{
  tcga = read.table("TCGA.txt")
  tcga$nonsynonymous <- tcga$V2
  tcga$overall_survival <- tcga$V3
}

nonsynonymous = list()
overall_survival  = list()
# Import and rearrange Van Allen data
if (Van_Allen == T)
{
  total = read.csv("Van_Allen.csv")
  nonsynonymous_Van_Allen = total$nonsynonymous
  overall_survival_Van_Allen = total$overall_survival
}
# Import and rearrange RIZVI data
if (Rizvi == T)
{
  total = read.csv("Rizvi.csv")
  total$nonsynonymous <- total$Nonsyn.
  total$overall_survival <- total$PFS..mos.*365
  nonsynonymous_Rizvi = total$nonsynonymous
  overall_survival_Rizvi = total$overall_survival
}
# Import and rearrange SNYDER data
if (Snyder == T)
{
  S1 = read.table("Snyder1.txt", comment.char = "%", header = T)
  S2 = read.table("Snyder2.txt", comment.char = "%", header = T)
  S3 = read.table("Snyder3.txt", comment.char = "%", header = T)
  bound_S1S2 <- rbind(S1,S2)
  total <- merge(bound_S1S2, S3, by.y = "Study_ID")
  total$nonsynonymous <- total$Mutation
  total$overall_survival <- total$OS.yr. *365
  nonsynonymous_Snyder = total$nonsynonymous
  overall_survival_Snyder = total$overall_survival
}

plot(nonsynonymous_Van_Allen, overall_survival_Van_Allen, col = "red", pch ="*")
points(nonsynonymous_Rizvi, overall_survival_Rizvi, col = "purple", pch ="*")
points(nonsynonymous_Snyder, overall_survival_Snyder, col = "blue", pch ="*")
points(tcga$nonsynonymous, tcga$overall_survival, col = "green", pch ="*")

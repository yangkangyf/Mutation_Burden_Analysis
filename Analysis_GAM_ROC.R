#########################
# Adjustable parameters #
#########################
# Number of simulations for the permutation test
number_of_simulations = 1000
# Discriminant number of mutation to perform the survival analysis
mutation_discriminant = 192
#Replace eight patients from Snyder's data
replace_eight_patients = F

#################################
# Choose which tests to perform #
#################################
# Choose which data to analyse
Van_Allen = T
Snyder = T
Rizvi = T
Hugo = T
TCGA = T
Snyder2 = T

#######################
# Librairies and data #
#######################
# Load libraries
library(survival)
library(mgcv)
library(gdata)
library(pracma)
library(flux)

# Run the analysis for selected data
selected_data = list()
if (Hugo == T)
  selected_data[length(selected_data)+1]="Hugo"
if (Van_Allen == T)
  selected_data[length(selected_data)+1]="Van_Allen"
if (Rizvi == T)
  selected_data[length(selected_data)+1]="Rizvi"
if (Snyder == T)
  selected_data[length(selected_data)+1]="Snyder"
if (Snyder2 == T)
  selected_data[length(selected_data)+1]="Snyder2"

total <- data.frame(
  ticker=character(),
  value=numeric(),
  date = as.Date(character()),
  stringsAsFactors=FALSE
)

for (which_data in selected_data)
{
  if (TCGA == T)
  {
    tcga = read.table("TCGA.txt")
    tcga$nonsynonymous <- tcga$V2
    tcga$overall_survival <- tcga$V3
  }
  
  # Import and rearrange Hugo data
  if (which_data == "Hugo")
  {
    data1= read.table("Hugo1.txt", header = T)
    data2 = read.table("Hugo2.txt", header = T, sep = "\t")
    data <- merge(data1, data2, by.y = "Patient_ID")
    nonsynonymous <- data$TotalNonSyn
    group <- (data$Response == "R")*1
    OS_or_PFS <- "Overall survival"
    age <- data$Age
    overall_survival <- data$Overall_Survival
    dead <- (data$Vital_Status == "Dead")*1
    stage <- data$Disease_Status
    stage1_str = "M0"
    stage2_str = "M1a"
    stage3_str = "M1b"
    stage4_str = "M1c"
    gender <- data$Gender
    male_str = "M"
    female_str = "F"
    rm(data1, data2, data)
  }
  
  # Import and rearrange Van Allen data
  if (which_data == "Van_Allen")
  {
    data = read.csv("Van_Allen.csv")
    nonsynonymous <- data$nonsynonymous
    overall_survival <- data$overall_survival
    group <- (data$group =="response")*1
    dead <- data$dead
    OS_or_PFS <- "Overall survival"
    age <- data$age_start
    stage <- data$M
    stage1_str = "M0"
    stage2_str = "M1a"
    stage3_str = "M1b"
    stage4_str = "M1c"
    gender <- data$gender
    male_str = "male"
    female_str = "female"
    rm(data)
  }
  
  # Import and rearrange RIZVI data
  if (which_data == "Rizvi")
  {
    data = read.csv("Rizvi.csv")
    nonsynonymous <- data$Nonsyn.
    overall_survival <- data$PFS..mos.*365
    group <- (data$Durable.Clinical.Benefit == "DCB")*1
    dead <- data$Event....
    OS_or_PFS <- "Progression free survival"
    age <- data$Age..years.
    gender <- data$Sex
    male_str = "M"
    female_str = "F"
    stage <- data$stage
    rm(data)
  }
  
  # Import and rearrange SNYDER data
  if (which_data == "Snyder")
  {
    S1 = read.table("Snyder1.txt", comment.char = "%", header = T)
    S2 = read.table("Snyder2.txt", comment.char = "%", header = T)
    S3 = read.table("Snyder3.txt", comment.char = "%", header = T)
    bound_S1S2 <- rbind(S1,S2)
    data <- merge(bound_S1S2, S3, by.y = "Study_ID")
    nonsynonymous <- data$Mutation
    overall_survival <- data$OS.yr. *365
    dead <- data$Alive_at_time_of_censure
    eightpatient <- rbind(data[48,],data[22,],data[23,],data[24,],
                          data[25,],data[13,],data[14,],data[15,])
    OS_or_PFS <- "Overall survival"
    age <- data$Age

    group <- (data$Benefit == "1")*1
    gender <- data$Gender
    male_str = "M"
    female_str = "F"
    stage <- data$M_stage
    stage1_str = "IIIc"
    stage2_str = "M1a"
    stage3_str = "M1b"
    stage4_str = "M1c"
    rm(data)
  }
  
  if (which_data == "Snyder2")
  {
    S1 = read.table("Snyder1.txt", comment.char = "%", header = T)
    S2 = read.table("Snyder2.txt", comment.char = "%", header = T)
    S3 = read.table("Snyder3.txt", comment.char = "%", header = T)
    bound_S1S2 <- rbind(S1,S2)
    data <- merge(bound_S1S2, S3, by.y = "Study_ID")
    nonsynonymous <- data$Mutation
    overall_survival <- data$OS.yr. *365
    dead <- data$Alive_at_time_of_censure
    eightpatient <- rbind(data[48,],data[22,],data[23,],data[24,],
                          data[25,],data[13,],data[14,],data[15,])
    OS_or_PFS <- "Overall survival"
    age <- data$Age
    data$Benefit[48] = 1
    data$Benefit[13] = 1
    data$Benefit[14] = 1
    data$Benefit[15] = 1
    data$Benefit[22] = 1
    data$Benefit[23] = 1
    data$Benefit[24] = 1
    data$Benefit[25] = 1
    group <- (data$Benefit == "1")*1
    gender <- data$Gender
    male_str = "M"
    female_str = "F"
    stage <- data$M_stage
    stage1_str = "IIIc"
    stage2_str = "M1a"
    stage3_str = "M1b"
    stage4_str = "M1c"
    rm(data)
  }
  
  total_temp = as.data.frame(cbind(nonsynonymous, group, age,overall_survival,dead, gender, 
                                   as.list(replicate(length(nonsynonymous), which_data))))
  total <- rbind(total_temp, total)
  total <- data.frame(lapply(total, as.character), stringsAsFactors=FALSE)
  
}
  
  # Distinguish groups according to their disease stage
  total_stage1 = total[which(stage == stage1_str),]
  total_stage2 = total[which(stage == stage2_str),]
  total_stage3 = total[which(stage == stage3_str),]
  total_stage4 = total[which(stage == stage4_str),]
  
  # Distinguish groups according to their Gender
  total_male = total[which(gender == male_str),]
  total_female = total[which(gender == female_str),]

##############
# ROC curves #
##############

TPR_list = list()
FPR_list = list()
nonsynonymous_list <- list()
name_list <- list()
legend_list = as.character()
count_legend=1
for (which_data2 in c('Rizvi', 'Snyder', 'Van_Allen', 'Hugo', 'Snyder2'))
{
  k = 1
  TPR_list_temp = list()
  FPR_list_temp = list()
  nonsynonymous_list_temp = list()
  name_list_temp = list()
  for (i in min(na.omit(as.numeric(total$nonsynonymous))):max(na.omit(as.numeric(total$nonsynonymous)))) 
  {
    TP = total[which(na.omit(as.numeric(total$nonsynonymous)) > i &
                       (na.omit(as.numeric(total$group)) == 1) & (na.omit(total$V7) == which_data2)),]
    TPR_list_temp[k] = length(na.omit(as.numeric(TP$nonsynonymous)))/length(which((na.omit(as.numeric(total$group)) == 1) & (na.omit(total$V7) == which_data2)))
    FP = total[which(na.omit(as.numeric(total$nonsynonymous)) > i &
                       (na.omit(as.numeric(total$group)) == 0) & (na.omit(total$V7) == which_data2)),]
    FPR_list_temp[k] = length(na.omit(as.numeric(FP$nonsynonymous)))/length(which((na.omit(as.numeric(total$group)) == 0) & (na.omit(total$V7) == which_data2)))
    nonsynonymous_list_temp[k] = i
    name_list_temp[k] <- which_data2
    k = k+1
  }
  TPR_list <- append(TPR_list, as.numeric(TPR_list_temp)) 
  FPR_list <- append(FPR_list, as.numeric(FPR_list_temp)) 
  nonsynonymous_list <- append(nonsynonymous_list, as.numeric(nonsynonymous_list_temp))
  name_list <- append(name_list, name_list_temp)
  count_legend = count_legend+1
}
data <- as.data.frame(cbind(unlist(TPR_list), unlist(FPR_list), as.numeric(nonsynonymous_list), as.list(name_list),unlist(unlist(FPR_list))^2 + rep(1,length(unlist(TPR_list))) - unlist(unlist(TPR_list))^2))

TPR_list_all = list()
FPR_list_all = list()
nonsynonymous_list_all = list()
k= 1
for (i in min(na.omit(as.numeric(total$nonsynonymous))):max(na.omit(as.numeric(total$nonsynonymous)))) 
{
  TP = total[which(na.omit(as.numeric(total$nonsynonymous)) > i &
                     (na.omit(as.numeric(total$group)) == 1)) ,]
  TPR_list_all[k] = length(na.omit(as.numeric(TP$nonsynonymous)))/length(which((na.omit(as.numeric(total$group)) == 1)))
  FP = total[which(na.omit(as.numeric(total$nonsynonymous)) > i &
                     (na.omit(as.numeric(total$group)) == 0)) ,]
  FPR_list_all[k] = length(na.omit(as.numeric(FP$nonsynonymous)))/length(which((na.omit(as.numeric(total$group)) == 0) ))
  nonsynonymous_list_all[k] = i
  k = k+1
}
data_all <- as.data.frame(cbind(unlist(TPR_list_all), unlist(FPR_list_all), as.numeric(nonsynonymous_list), unlist(unlist(FPR_list_all))^2 + rep(1,length(unlist(FPR_list_all))) - unlist(unlist(TPR_list_all))^2))


plot(FPR_list_all, TPR_list_all, type = "l",
     col = "black", xlab = "(1-Specificity)", ylab = "Sensitivity",
     xlim = c(0,1.1), ylim = c(0,1.1))
title("ROC curve analysis")
points(sort(unlist(data[which(data$V4 == "Rizvi"),]$V2)), sort(unlist(data[which(data$V4 == "Rizvi"),]$V1)), type = "l",
     col = "green")
points(sort(unlist(data[which(data$V4 == "Snyder"),]$V2)), sort(unlist(data[which(data$V4 == "Snyder"),]$V1)), type = "l",
       col = "red")
points(sort(unlist(data[which(data$V4 == "Hugo"),]$V2)), sort(unlist(data[which(data$V4 == "Hugo"),]$V1)), type = "l",
       col = "purple")
points(sort(unlist(data[which(data$V4 == "Van_Allen"),]$V2)), sort(unlist(data[which(data$V4 == "Van_Allen"),]$V1)), type = "l",
       col = "blue")
points(sort(unlist(data[which(data$V4 == "Snyder2"),]$V2)), sort(unlist(data[which(data$V4 == "Snyder2"),]$V1)), type = "l",
       col = "orange")
abline(a = 0, b =1) 
legend_list[1] = paste("All datas (n=246)", ",  AUC = ",as.character(round(auc(unlist(FPR_list), unlist(TPR_list)),3)))
legend_list[2] = paste("Rizvi (n=34)", ",  AUC = ",as.character(round(auc(sort(unlist(data[which(data$V4 == "Rizvi"),]$V2)), sort(unlist(data[which(data$V4 == "Rizvi"),]$V1))),3)))
legend_list[3] = paste("Snyder (n=64)", ",  AUC = ",as.character(round(auc(sort(unlist(data[which(data$V4 == "Snyder"),]$V2)), sort(unlist(data[which(data$V4 == "Snyder"),]$V1))),3)))
legend_list[4] = paste("Hugo (n=38)", ",  AUC = ",as.character(round(auc(sort(unlist(data[which(data$V4 == "Hugo"),]$V2)), sort(unlist(data[which(data$V4 == "Hugo"),]$V1))),3)))
legend_list[5] = paste("Van Allen (n=110)", ",  AUC = ",as.character(round(auc(sort(unlist(data[which(data$V4 == "Van_Allen"),]$V2)), sort(unlist(data[which(data$V4 == "Van_Allen"),]$V1))),3)))
legend_list[6] = paste("Snyder2 (n=64)", ",  AUC = ",as.character(round(auc(sort(unlist(data[which(data$V4 == "Snyder2"),]$V2)), sort(unlist(data[which(data$V4 == "Snyder2"),]$V1))),3)))

legend("bottomright", # places a legend at the appropriate place 
       legend_list, # puts text in the legend
       pch = c(".",".", ".", ".",".","."), # gives the legend appropriate symbols (lines),
       lwd = c(1,1,1,1,1,1),
       col=c('black','green',"blue", "purple", "red", "orange")) # gives the legend lines the correct color and width


cutoff_all <- as.numeric(which.min(unlist(data_all$V4)))
cutoff_Rizvi <- which(data$V4 == "Rizvi")[as.numeric(which.min(unlist(data[which(data$V4 == "Rizvi"),]$V5)))]
cutoff_Snyder <- which(data$V4 == "Snyder")[as.numeric(which.min(unlist(data[which(data$V4 == "Snyder"),]$V5)))]
cutoff_Hugo <- which(data$V4 == "Hugo")[as.numeric(which.min(unlist(data[which(data$V4 == "Hugo"),]$V5)))]
cutoff_VanAllen <- which(data$V4 == "Van_Allen")[as.numeric(which.min(unlist(data[which(data$V4 == "Van_Allen"),]$V5)))]
cutoff_Snyder2 <- which(data$V4 == "Snyder2")[as.numeric(which.min(unlist(data[which(data$V4 == "Snyder2"),]$V5)))]


points(unlist(data_all$V2[cutoff_all]), unlist(data_all$V1[cutoff_all]), pch = "*", cex = 3, col = "black")
text(unlist(data_all$V2[cutoff_all]), unlist(data_all$V1[cutoff_all]) + 0.04, data_all$V3[cutoff_all], col = "black")
points(unlist(data$V2[cutoff_Rizvi]), unlist(data$V1[cutoff_Rizvi]), pch = "*", cex = 3, col = "green")
text(unlist(data$V2[cutoff_Rizvi]), unlist(data$V1[cutoff_Rizvi]) + 0.04, data$V3[cutoff_Rizvi], col = "green")
points(unlist(data$V2[cutoff_Snyder]), unlist(data$V1[cutoff_Snyder]), pch = "*", cex = 3, col = "blue")
text(unlist(data$V2[cutoff_Snyder]), unlist(data$V1[cutoff_Snyder]) + 0.04, data$V3[cutoff_Snyder], col = "blue")
points(unlist(data$V2[cutoff_Hugo]), unlist(data$V1[cutoff_Hugo]), pch = "*", cex = 3, col = "purple")
text(unlist(data$V2[cutoff_Hugo]), unlist(data$V1[cutoff_Hugo]) + 0.04, data$V3[cutoff_Hugo], col = "purple")
points(unlist(data$V2[cutoff_VanAllen]), unlist(data$V1[cutoff_VanAllen]), pch = "*", cex = 3, col = "red")
text(unlist(data$V2[cutoff_VanAllen]), unlist(data$V1[cutoff_VanAllen]) + 0.04, data$V3[cutoff_VanAllen], col = "red")
points(unlist(data$V2[cutoff_Snyder2]), unlist(data$V1[cutoff_Snyder2]), pch = "*", cex = 3, col = "orange")
text(unlist(data$V2[cutoff_Snyder2]), unlist(data$V1[cutoff_Snyder2]) + 0.04, data$V3[cutoff_Snyder2], col = "orange")



##############################
# Generalized additive model #
##############################  

#GAM_obj= gam(log(na.omit(as.numeric(total$overall_survival))) ~ na.omit(as.numeric(total$nonsynonymous)))
#plot(log((as.numeric(total$overall_survival))) ~ (as.numeric(total$nonsynonymous)), 
#     pch = '.', cex = 5)
#points(log(tcga$overall_survival) ~ tcga$nonsynonymous, 
#       pch = '.', col = 'red', cex = 5)

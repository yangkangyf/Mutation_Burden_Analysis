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

#######################
# Librairies and data #
#######################
# Load libraries
library(survival)
library(mgcv)
library(gdata)

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
    
    if (replace_eight_patients == T)
    {
      data$Benefit[48] = 1
      data$Benefit[13] = 1
      data$Benefit[14] = 1
      data$Benefit[15] = 1
      data$Benefit[22] = 1
      data$Benefit[23] = 1
      data$Benefit[24] = 1
      data$Benefit[25] = 1
    }
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
  
  if (which_data == "Snyder")
  {
    eightpatient_benefit = eightpatient[which(eightpatient$group == clinical_benefit_str),]
    eightpatient_nobenefit = eightpatient[which(eightpatient$group == clinical_nobenefit_str),]
  }
}
  
  # Distinguish groups according to their disease stage
  total_stage1 = total[which(stage == stage1_str),]
  total_stage2 = total[which(stage == stage2_str),]
  total_stage3 = total[which(stage == stage3_str),]
  total_stage4 = total[which(stage == stage4_str),]
  
  # Distinguish groups according to their Gender
  total_male = total[which(gender == male_str),]
  total_female = total[which(gender == female_str),]

  
##############################
# Generalized additive model #
##############################  
  
#GAM_obj= gam(log(total$overall_survival) ~ total$nonsynonymous)
#plot(log(total$overall_survival) ~ total$nonsynonymous, 
#     pch = '.', cex = 3)
#points(log(tcga$overall_survival) ~ tcga$nonsynonymous, 
#      pch = '.', col = 'red', cex = 3)

##############
# ROC curves #
##############

TPR_list = list()
FPR_list = list()

j = 0
for (i in min(na.omit(as.numeric(total$nonsynonymous))):max(na.omit(as.numeric(total$nonsynonymous)))) 
{
TP = total[which(na.omit(as.numeric(total$nonsynonymous)) > i &
                   (na.omit(as.numeric(total$group)) == 1)),]
TPR_list[j] = length(na.omit(as.numeric(TP$nonsynonymous)) )/sum((na.omit(as.numeric(total$group)) == 1)*1)
FP = total[which(na.omit(as.numeric(total$nonsynonymous)) > i &
                   (na.omit(as.numeric(total$group)) == 0)),]  
FPR_list[j] = length(na.omit(as.numeric(FP$nonsynonymous)))/sum((na.omit(as.numeric(total$group)) == 0)*1)
j = j+1
}

TPR_list <- as.numeric(TPR_list)
FPR_list <- as.numeric(FPR_list)
plot(FPR_list,TPR_list, type = "l") 
abline(a = 0, b =1) 

rm(TPR_list, FPR_list)

for (which_data2 in c('Rizvi', 'Snyder', 'Van_Allen', 'Hugo'))
{
  if (which_data2 == 'Snyder')
  {col_string = 'blue'}
  if (which_data2 == 'Van_Allen')
  {col_string = 'red'}
  if (which_data2 == 'Rizvi')
  {col_string = 'green'}
  if (which_data2 == 'Hugo')
  {col_string = 'purple'}
  j = 0
  TPR_list = list()
  FPR_list = list()
  for (i in min(na.omit(as.numeric(total$nonsynonymous))):max(na.omit(as.numeric(total$nonsynonymous)))) 
  {
    TP = total[which(na.omit(as.numeric(total$nonsynonymous)) > i &
                       (na.omit(as.numeric(total$group)) == 1) & (na.omit(total$V7) == which_data2)),]
    TPR_list[j] = length(na.omit(as.numeric(TP$nonsynonymous)))/length(which((na.omit(as.numeric(total$group)) == 1) & (na.omit(total$V7) == which_data2)))
    FP = total[which(na.omit(as.numeric(total$nonsynonymous)) > i &
                       (na.omit(as.numeric(total$group)) == 0) & (na.omit(total$V7) == which_data2)),]
    FPR_list[j] = length(na.omit(as.numeric(FP$nonsynonymous)))/length(which((na.omit(as.numeric(total$group)) == 0) & (na.omit(total$V7) == which_data2)))
    j = j+1
    TPR_list <- as.numeric(TPR_list)
    FPR_list <- as.numeric(FPR_list)
  }
  points(FPR_list,TPR_list, type = "l", col = col_string) 
}



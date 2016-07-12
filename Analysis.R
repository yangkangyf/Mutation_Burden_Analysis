params_tests = list()
params_data = list()
params_values = list()

# Path to photos
path = "~/Figures/"

#########################
# Adjustable parameters #
#########################
# Number of simulations for the permutation test
params_values[1] = 10000
# Discriminant number of mutation to perform the survival analysis
params_values[2] = 102

#################################
# Choose which tests to perform #
#################################
# Choose which data to analyse
Van_Allen = T
Snyder = T
Snyder2 = F
Rizvi = T
Hugo = T


# Plot mutational load VS survival
params_tests[1] = F
# Boxplot of mutational loads for each group
params_tests[2]  = F
# Plot permutation test density function for Mann-Whitney test (by shuffling mutational load)
params_tests[4]  = F
# Plot p-value of permutation test density function for log-rank test (by shuffling overall survival) for all discriminant number of mutations
params_tests[6]  = T

# Trimmed analysis
params_tests[7]  = F
# Trimmed survival analysis
params_tests[8]  = F
# Plot K-M curves for the chosen mutation_discriminant
params_tests[3]  = F
# Plot permutation test density function for log-rank test (by shuffling overall survival)
params_tests[5]  = F


if (FALSE)
{
# Plot mutational load VS survival with points size according to age
params[1]  = T
# Plot sorted mutational loads for each group (responders and non-responders)
params[1]  = T

# Plot spline and find optimal number of mutations
params[1]  = F
# Plot mutational load VS survival by distinguishing by stage severity
params[1]  = F
# Plot mutational load VS survival by distinguishing by gender
params[1]  = F
# Test potential pitfalls (age, gender, stage)
params[1]  = F
# Compute the trimmed data analysis
}
#######################
# Librairies and data #
#######################
# Load libraries
library(survival)
library(ggplot2)
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
if (Snyder2 == T)
  selected_data[length(selected_data)+1]="Snyder2"
rm(Snyder2, Snyder, Van_Allen, Rizvi, Hugo)



for (which_data in selected_data)
{
# Import and rearrange Van Allen data
if (which_data == "Van_Allen")
{
  data = read.csv("Van_Allen.csv")
  nonsynonymous <- data$nonsynonymous
  overall_survival <- data$overall_survival
  group <- data$group
  params_data[1] <- "response"
  params_data[2] <- "nonresponse"
  dead <- data$dead
  params_data[3] <- "Overall survival"
  age <- data$age_start
  stage <- data$M
  params_data[4] = "M0"
  params_data[5] = "M1a"
  params_data[6] = "M1b"
  params_data[7] = "M1c"
  gender <- data$gender
  params_data[8] = "male"
  params_data[9] = "female"
  special_patients <- which(data$overall_survival > 2*365 & data$progression_free < 6*30 & data$RECIST == "PD")
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
  params_data[1] <- "1"
  params_data[2] <- "0"
  dead <- data$Alive_at_time_of_censure
  special_patients <- c(48,25,22,23,24,25,13,14,15) 
  params_data[3] <- "Overall survival"
  age <- data$Age
  group <- data$Benefit
  gender <- data$Gender
  params_data[8] = "M"
  params_data[9] = "F"
  stage <- data$M_stage
  params_data[4] = "IIIc"
  params_data[5] = "M1a"
  params_data[6] = "M1b"
  params_data[7] = "M1c"
  rm(data)
}

  
# Import and rearrange RIZVI data
if (which_data == "Rizvi")
{
  data = read.csv("Rizvi.csv")
  nonsynonymous <- data$Nonsyn.
  overall_survival <- data$PFS..mos.*365
  group <- data$Durable.Clinical.Benefit
  params_data[1] <- "DCB"
  params_data[2] <- "NDB"
  dead <- data$Event....
  params_data[3] <- "Progression free survival"
  age <- data$Age..years.
  gender <- data$Sex
  params_data[8]  = "M"
  params_data[9]  = "F"
  stage <- data$stage
  rm(data)
}
  
  
# Import and rearrange Hugo data
if (which_data == "Hugo")
  {
    data1= read.table("Hugo1.txt", header = T)
    data2 = read.table("Hugo2.txt", header = T, sep = "\t")
    data <- merge(data1, data2, by.y = "Patient_ID")
    nonsynonymous <- data$TotalNonSyn
    group <- data$Response
    params_data[1] <- "R"
    params_data[2] <- "NR"
    params_data[3] <- "Overall survival"
    age <- data$Age
    overall_survival <- data$Overall_Survival
    dead <- (data$Vital_Status == "Dead")*1
    stage <- data$Disease_Status
    params_data[4] = "M0"
    params_data[5] = "M1a"
    params_data[6] = "M1b"
    params_data[7] = "M1c"
    gender <- data$Gender
    params_data[8] = "M"
    params_data[9] = "F"
    special_patients <- c(9)
    rm(data1, data2, data)
}

# Import and rearrange SNYDER2 data
if (which_data == "Snyder2")
{
  S1 = read.table("Snyder1.txt", comment.char = "%", header = T)
  S2 = read.table("Snyder2.txt", comment.char = "%", header = T)
  S3 = read.table("Snyder3.txt", comment.char = "%", header = T)
  bound_S1S2 <- rbind(S1,S2)
  data <- merge(bound_S1S2, S3, by.y = "Study_ID")
  nonsynonymous <- data$Mutation
  overall_survival <- data$OS.yr. *365
  clinical_benefit_str <- "1"
  clinical_nobenefit_str <- "0"
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
  group <- data$Benefit
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
  
total = as.data.frame(cbind(as.numeric((nonsynonymous)), as.character(group), 
                            as.numeric((age)),
                            as.numeric((overall_survival)),
                            as.numeric((dead)), 
                            as.character(gender)), as.is = F,
                            stringsAsFactors = FALSE)
total <- drop.levels(total)
total$V1 <- as.numeric(total$V1)
total$V3 <- as.numeric(total$V3)
total$V4 <- as.numeric(total$V4)
total$V5 <- as.numeric(total$V5)
names(total)[1] <- "nonsynonymous"
names(total)[2] <- "group"
names(total)[3] <- "age"
names(total)[4] <- "overall_survival"
names(total)[5] <- "dead"
names(total)[6] <- "gender"
#names(total)[7] <- "stage"

rm(nonsynonymous, group, age,overall_survival,dead, stage, gender)
# Distinguish groups according to clinical benefit  
#total_benefit = total[which(group == clinical_benefit_str),]
#total_nobenefit = total[which(group == clinical_nobenefit_str),]

# Distinguish groups according to their disease stage
#total_stage1 = total[which(stage == stage1_str),]
#total_stage2 = total[which(stage == stage2_str),]
#total_stage3 = total[which(stage == stage3_str),]
#total_stage4 = total[which(stage == stage4_str),]

# Distinguish groups according to their Gender
#total_male = total[which(gender == male_str),]
#total_female = total[which(gender == female_str),]

params_data = unlist(params_data)
params_values <- as.numeric(params_values)

###############
# Basic plots #
###############


if (params_tests[1] ==T)
{
# Plot Survival against number of nonsynonymous mutations
par(mfrow=c(1,1))
namefile = paste(path, "survVSmut_", which_data, ".tiff")
tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
plot(total[which(total$group == params_data[1]),]$nonsynonymous, total[which(total$group == params_data[1]),]$overall_survival,
     xlab= "Number of Nonsynonymous mutations",
     ylab= params_data[3], pch = 16,col = "cyan3", cex = 1, main = which_data,      
     xlim = c(min(total$nonsynonymous),max(total$nonsynonymous)), 
     ylim= c(min(na.omit(total$overall_survival)),max(na.omit(total$overall_survival))))
points(total[which(total$group == params_data[1]),][which(total[which(total$group == params_data[1]),]$dead == 1),]$nonsynonymous, 
       total[which(total$group == params_data[1]),][which(total[which(total$group == params_data[1]),]$dead == 1),]$overall_survival,
       pch = 1, cex = 1)
points(total[which(total$group == params_data[2]),]$nonsynonymous, 
       total[which(total$group == params_data[2]),]$overall_survival, col ="red",
       pch = 16, cex = 1)
points(total[which(total$group == params_data[2]),][which(total[which(total$group == params_data[2]),]$dead == 1),]$nonsynonymous, 
       total[which(total$group == params_data[2]),][which(total[which(total$group == params_data[2]),]$dead == 1),]$overall_survival, xlab= "Number of Nonsynonymous mutations",
       ylab= params_data[3], pch = 1, cex = 1)
legend("topright", # places a legend at the appropriate place 
       c("Responders","Non-responders", "Dead"), # puts text in the legend
       pch = c(16,16,1),
       col=c("cyan3", "red", "black")) # gives the legend lines the correct color and width
if (which_data == "Snyder" || which_data == "Snyder2")
{
  points(total[special_patients,]$nonsynonymous, total[special_patients,]$overall_survival, 
         col ="grey", cex= 1, pch = 16)
  points(total[special_patients,][which(total[special_patients,]$dead == 1),]$nonsynonymous, 
         total[special_patients,][which(total[special_patients,]$dead == 1),]$overall_survival
         , pch = 1, cex = 1)
  legend("topright", # places a legend at the appropriate place 
         c("Responders","Non-responders","Dead", "Reassigned patients"), # puts text in the legend
         pch = c(16,16,1,16), # gives the legend appropriate symbols (lines),
         col=c("cyan", "red","black", "grey")) # gives the legend lines the correct color and width
}
if (which_data == "Van_Allen")
{
  points(total[special_patients,]$nonsynonymous, total[special_patients,]$overall_survival, 
         col ="grey", cex= 1, pch = 16)
  points(total[special_patients,][which(total[special_patients,]$dead == 1),]$nonsynonymous, 
         total[special_patients,][which(total[special_patients,]$dead == 1),]$overall_survival
         , pch = 1, cex = 1)
  legend("topright", # places a legend at the appropriate place 
         c("Responders","Non-responders","Dead", "Subset"), # puts text in the legend
         pch = c(16,16,1,16), # gives the legend appropriate symbols (lines),
         col=c("cyan", "red","black", "grey")) # gives the legend lines the correct color and width
}

dev.off()

# Same figure, but in log scale
par(mfrow=c(1,1))
namefile = paste(path,"survVSlogmut_", which_data, ".tiff")
tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
plot(log(total[which(total$group == params_data[1]),]$nonsynonymous), total[which(total$group == params_data[1]),]$overall_survival,
     xlab= "Log number of Nonsynonymous mutations",
     ylab= params_data[3], pch = 16,col = "cyan3", cex = 1, main = which_data,      
     xlim = c(min(log(total$nonsynonymous)),max(log(total$nonsynonymous))), 
     ylim= c(min(na.omit(total$overall_survival)),max(na.omit(total$overall_survival))))
points(log(total[which(total$group == params_data[1]),][which(total[which(total$group == params_data[1]),]$dead == 1),]$nonsynonymous), 
       total[which(total$group == params_data[1]),][which(total[which(total$group == params_data[1]),]$dead == 1),]$overall_survival,
       pch = 1, cex = 1)
points(log(total[which(total$group == params_data[2]),]$nonsynonymous), 
       total[which(total$group == params_data[2]),]$overall_survival, col ="red",
       pch = 16, cex = 1)
points(log(total[which(total$group == params_data[2]),][which(total[which(total$group == params_data[2]),]$dead == 1),]$nonsynonymous), 
       total[which(total$group == params_data[2]),][which(total[which(total$group == params_data[2]),]$dead == 1),]$overall_survival, xlab= "Number of Nonsynonymous mutations",
       ylab= params_data[3], pch = 1, cex = 1)
legend("topright", # places a legend at the appropriate place 
       c("Responders","Non-responders", "Dead"), # puts text in the legend
       pch = c(16,16,1),
       col=c("cyan3", "red", "black")) # gives the legend lines the correct color and width
if (which_data == "Snyder" || which_data == "Snyder2")
{
  points(log(total[special_patients,]$nonsynonymous), total[special_patients,]$overall_survival, 
         col ="grey", cex= 1, pch = 16)
  points(log(total[special_patients,][which(total[special_patients,]$dead == 1),]$nonsynonymous), 
         total[special_patients,][which(total[special_patients,]$dead == 1),]$overall_survival
         , pch = 1, cex = 1)
  legend("topright", # places a legend at the appropriate place 
         c("Responders","Non-responders","Dead", "Reassigned"), # puts text in the legend
         pch = c(16,16,1,16), # gives the legend appropriate symbols (lines),
         col=c("cyan", "red","black", "grey")) # gives the legend lines the correct color and width
}
if (which_data == "Van_Allen")
{
  points(log(total[special_patients,]$nonsynonymous), total[special_patients,]$overall_survival, 
         col ="grey", cex= 1, pch = 16)
  points(log(total[special_patients,][which(total[special_patients,]$dead == 1),]$nonsynonymous), 
         total[special_patients,][which(total[special_patients,]$dead == 1),]$overall_survival
         , pch = 1, cex = 1)
  legend("topright", # places a legend at the appropriate place 
         c("Responders","Non-responders","Dead", "Subset"), # puts text in the legend
         pch = c(16,16,1,16), # gives the legend appropriate symbols (lines),
         col=c("cyan", "red","black", "grey")) # gives the legend lines the correct color and width
}
dev.off()
}


##############################################################
# Association between Mutational Burden and Clinical Benefit #
##############################################################
if (params_tests[2] == T)
{
### Plot number of nonsynonymous mutations against benefit
par(oma = c(4, 1, 1, 1))
namefile = paste(path,"mutBox", which_data, ".tiff")
tiff(namefile, width = 8, height = 6, units = 'in', res = 500)


if (nlevels(as.factor(total$group)) == 2)
{
  col_ = c("red", "cyan3")
  names_ = c("Non-responders", "Responders")
}

if (nlevels(as.factor(total$group)) == 3)
{
  col_ =  c("grey", "red", "cyan3")
  names_ = c("Subset", "Non-responders", "Responders")
}

boxplot(log(total$nonsynonymous) ~ as.factor(total$group), data = total, lwd = 2, 
         pwcol = 1 + as.numeric(total$dead),offset = .5,
         ylab = 'Log number of nonymomymous mutations',
         xlab = 'Response categories',
         col = col_,
         names = names_)

stripchart(log(total[which(total$dead == 0),]$nonsynonymous) ~ as.factor(total[which(total$dead == 0),]$group), vertical = TRUE, data = total[which(total$dead == 0),], 
           method = "jitter", add = TRUE, pch = 16, col = 'black')
stripchart(log(total[which(total$dead == 1),]$nonsynonymous) ~ as.factor(total[which(total$dead == 1),]$group), vertical = TRUE, data = total[which(total$dead == 1),], 
           method = "jitter", add = TRUE, pch = 1, col = 'black')
wilcox.test(total[which(total$group == params_data[1]),]$nonsynonymous, 
           total[which(total$group == params_data[2]),]$nonsynonymous)
dev.off()  
}


############################
# Permutation Mann-Whitney #
############################
if (params_tests[4] == T)
{
  # Mann_Whitney (benefit vs no benefit)
  real_test <- wilcox.test(na.omit(total[which(total$group == params_data[1]),]$nonsynonymous), 
                           na.omit(total[which(total$group == params_data[2]),]$nonsynonymous))
  real_W <- as.numeric(as.character(real_test$statistic))
  list_i = list()
  list_W= vector()
  j = 1
  try(
    for (i in 1:params_values[1])
    {
      all_nonsynonymous <- append(total[which(total$group == params_data[1]),]$nonsynonymous, 
                                  total[which(total$group == params_data[2]),]$nonsynonymous)
      all_nonsynonymous = sample(all_nonsynonymous)
      total_benefit_shuffled = all_nonsynonymous[1:length(total[which(total$group == params_data[1]),]$nonsynonymous)]
      total_nobenefit_shuffled = all_nonsynonymous[length(total[which(total$group == params_data[1]),]$nonsynonymous)+1:length(all_nonsynonymous)]
      res <- wilcox.test(total_benefit_shuffled, total_nobenefit_shuffled)
      W = as.numeric(as.character(res$statistic))
      list_i[j] = i
      list_W[j] = W
      j = j+1
    })
  namefile = paste(path, "Perm_MW", which_data, ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  h <- hist(list_W, breaks = params_values[1]/50, main = which_data, xlab = "W statisitics", 
            xlim = c(min(list_W),max(real_W + real_W/20, max(list_W))))
  abline(v= real_W, col = 'red')
  pvalue_comp = length(which(list_W > real_W))/params_values[1]
  pvalue = paste("p-value = ",pvalue_comp)
  legend("topright", # places a legend at the appropriate place 
         c(pvalue), # puts text in the legend
         pch = c(".")) # gives the legend lines the correct color and width
  dev.off()
  rm(W, real_test,real_W,list_i,list_W,j,all_nonsynonymous, total_benefit_shuffled,total_nobenefit_shuffled,res,
     pvalue_comp,pvalue )
}

if (params_tests[6] == T)
{
  pvalue_comp = vector()
  mut_chosen_list = vector()
  k = 0
  try(
    for (mut_chosen in (min(total$nonsynonymous)):(max(total$nonsynonymous)))
    {
      real_test <- survdiff(Surv(total$overall_survival, total$dead) ~ total$nonsynonymous > mut_chosen)
      real_Z <- as.numeric(as.character(((real_test$obs[1]-real_test$exp[1])^2)/real_test$exp[1] + ((real_test$obs[2]-real_test$exp[2])^2)/real_test$exp[2]))
      list_i = list()
      list_Z= vector()
      j = 1
      for (i in 1:params_values[1])
      {
        overall_survival_shuffled = sample(total$overall_survival)
        res <- survdiff(Surv(overall_survival_shuffled, total$dead) ~ total$nonsynonymous > mut_chosen)
        Z = as.numeric(as.character(((res$obs[1]-res$exp[1]))/res$exp[1] + ((res$obs[2]-res$exp[2]))/res$exp[2]))
        list_i[j] = i
        list_Z[j] = Z
        j = j+1
      }
      pvalue_comp[k] = length(which(list_Z > real_Z))/params_values[2]
      mut_chosen_list[k] = mut_chosen
      k = k+1
    })
  namefile6 = paste("Permutation_Log-Rank_Multiple",which_data, ".jpg")
  jpeg(namefile6)
  plot(mut_chosen_list, pvalue_comp, xlim = c(min(total$nonsynonymous),max(total$nonsynonymous)), 
       ylim = c(min(pvalue_comp), max(pvalue_comp)), pch = ".",
       xlab = "Discriminant number of mutations", ylab = "p_value")
  lines(mut_chosen_list, pvalue_comp)
  dev.off()
  rm(list_i ,list_Z, real_test, real_Z, pvalue_comp, mut_chosen_list, mut_chosen, overall_survival_shuffled, res, Z)
}
}

##########################
# Log odd ratio analysis #
##########################
#RE = Responders-Exposed
#RN = Responders non exposed
#NE = Non-responders exposed
#NN = Non-responders non-exposed


OR_list = list()
OR_IC_u = list()
OR_IC_l = list()
j = 0
i_list = 0
try(
for (i in sort(unique(total$nonsynonymous)))
{
  RE <- total[which((total$group == as.character(params_data[1])) & (total$nonsynonymous > i)),]$nonsynonymous
  RE = length(RE)
  NE <- total[which((total$group == as.character(params_data[2])) & (total$nonsynonymous > i)),]$nonsynonymous
  NE = length(NE)
  RN <- total[which((total$group == as.character(params_data[1])) & (total$nonsynonymous < i)),]$nonsynonymous
  RN = length(RN)
  NN <- total[which((total$group == as.character(params_data[2])) & (total$nonsynonymous < i)),]$nonsynonymous
  NN = length(NN)
  OR_list[j] = log((RE/RN)/(NE/NN))
  OR_IC_u[j] = log((RE/RN)/(NE/NN)) + 1.96*sqrt(1/RE+1/NE+1/RN+1/NN)
  OR_IC_l[j] = log((RE/RN)/(NE/NN)) - 1.96*sqrt(1/RE+1/NE+1/RN+1/NN)
  i_list[j]=i
  j = j+1
})
plot((as.numeric(OR_list)), pch = ".", xlab = "Discriminant number of mutation",
     ylab= "Odd ratio", cex = 3)

d = NaRV.omit(data.frame(cbind(as.numeric(i_list),
                               as.numeric(OR_list),
                               as.numeric(OR_IC_l),
                               as.numeric(OR_IC_u))))
names(d) = c("x","y
             ","ylo","yhi")
namefile = paste(path,"Forestplot_", which_data, ".tiff")
tiff(namefile, width = 8, height = 12, units = 'in', res = 800)
forestplot(d$x,d$y, d$ylo, d$yhi,lwd.zero=1, reo = "red")
dev.off()

#####################
# Trimming analysis #
#####################
if (params_tests[7] == T)
{
  list_wilcox_trimmed = list()
  total_trimmed <- total[order(total$nonsynonymous),]
  for (i in 1:floor(length(total[which(total$group == params_data[1]),]$nonsynonymous)))
  {
  ### Plot number of nonsynonymous mutations against benefit
  total_trimmed <- head(total_trimmed,-1) #trim tail
  total_trimmed <- tail(total_trimmed,-1) #trim head
  list_wilcox_trimmed[i] <- wilcox.test(total_trimmed[which(total_trimmed$group == params_data[1]),]$nonsynonymous, 
              total_trimmed[which(total_trimmed$group == params_data[2]),]$nonsynonymous)$p.value
  }
  par(mfrow=c(1,1))
  namefile = paste(path,"Trimmed_data_", which_data, ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  plot(as.numeric(list_wilcox_trimmed), type = 'l', 
       xlab = "Number of values trimmed on each side",
       ylab = "P-value (Mann-Whitney test)")
  points(which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE),
         as.numeric(list_wilcox_trimmed)[which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE)],
         pch = "*", col = "red", cex = 2)
  abline(a=0.05, b=0, col = 'red') 
  dev.off()
  
namefile = paste(path,"Trimmed_boxplot_", which_data, ".tiff")
tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
total_trimmed <- total[order(total$nonsynonymous),]
total_trimmed <- head(total_trimmed,-which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE)) #trim tail
total_trimmed <- tail(total_trimmed,-which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE)) #trim head
boxplot(log(total_trimmed$nonsynonymous) ~ as.factor(total_trimmed$group), 
        data = total_trimmed, lwd = 2, 
        pwcol = 1 + as.numeric(total_trimmed$dead),offset = .5,
        ylab = 'Log number of nonymomymous mutations',
        xlab = 'Response categories',
        col = c("grey", "red", "cyan3"),
        names = c("Subset", "Non-responders", "Responders"))
dev.off()

list_wilcox_trimmed = list()
total_trimmed <- total[order(total$nonsynonymous),]
for (i in 1:floor(length(total[which(total$group == params_data[1]),]$nonsynonymous)))
{
  ### Plot number of nonsynonymous mutations against benefit
  total_trimmed <- head(total_trimmed,-1) #trim tail
  list_wilcox_trimmed[i] <- wilcox.test(total_trimmed[which(total_trimmed$group == params_data[1]),]$nonsynonymous, 
                                        total_trimmed[which(total_trimmed$group == params_data[2]),]$nonsynonymous)$p.value
}
par(mfrow=c(1,1))
namefile = paste(path,"Trimmed_right_", which_data, ".tiff")
tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
plot(as.numeric(list_wilcox_trimmed), type = 'l', 
     xlab = "Number of values trimmed on the right side",
     ylab = "P-value (Mann-Whitney test)")
points(which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE),
       as.numeric(list_wilcox_trimmed)[which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE)],
       pch = "*", col = "red", cex = 2)
abline(a=0.05, b=0, col = 'red') 
dev.off()


list_wilcox_trimmed = list()
total_trimmed <- total[order(total$nonsynonymous),]
for (i in 1:floor(length(total[which(total$group == params_data[1]),]$nonsynonymous)))
{
  ### Plot number of nonsynonymous mutations against benefit
  total_trimmed <- tail(total_trimmed,-1) #trim tail
  list_wilcox_trimmed[i] <- wilcox.test(total_trimmed[which(total_trimmed$group == params_data[1]),]$nonsynonymous, 
                                        total_trimmed[which(total_trimmed$group == params_data[2]),]$nonsynonymous)$p.value
}
par(mfrow=c(1,1))
namefile = paste(path,"Trimmed_right_", which_data, ".tiff")
tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
plot(as.numeric(list_wilcox_trimmed), type = 'l', 
     xlab = "Number of values trimmed on the left side",
     ylab = "P-value (Mann-Whitney test)")
points(which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE),
       as.numeric(list_wilcox_trimmed)[which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE)],
       pch = "*", col = "red", cex = 2)
abline(a=0.05, b=0, col = 'red') 
dev.off()


namefile = paste(path,"Trimmed_boxplot_right_", which_data, ".tiff")
tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
total_trimmed <- total[order(total$nonsynonymous),]
total_trimmed <- head(total_trimmed,-which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE)) #trim tail
total_trimmed <- tail(total_trimmed,-which(as.numeric(list_wilcox_trimmed) == min(as.numeric(list_wilcox_trimmed)), arr.ind = TRUE)) #trim head
boxplot(log(total_trimmed$nonsynonymous) ~ as.factor(total_trimmed$group), 
        data = total_trimmed, lwd = 2, 
        pwcol = 1 + as.numeric(total_trimmed$dead),offset = .5,
        ylab = 'Log number of nonymomymous mutations',
        xlab = 'Response categories',
        col = c("grey", "red", "cyan3"),
        names = c("Subset", "Non-responders", "Responders"))
dev.off()
}

#############################
# Trimmed Survival Analysis #
#############################
if (params_tests[8] == T)
{
  difference_survival_list = list()
  total_trimmed <- total[order(total$nonsynonymous),]
  for (i in 1:floor(length(total[which(total$group == params_data[1]),]$nonsynonymous)))
  {
    ### Plot number of nonsynonymous mutations against benefit
    total_trimmed_u <- head(total_trimmed,-1) #trim tail
    total_trimmed_l <- tail(total_trimmed,-1) #trim head
    med_R = mean(total_trimmed_u$overall_survival)
    med_NR = mean(total_trimmed_l$overall_survival)
    difference_survival_list[i] <- med_R - med_NR
  }
  namefile = paste(path,"Difference_mean_survival_trimmed_", which_data, ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  plot(as.numeric(difference_survival_list), type = "l", 
       xlab = "Number of values (NsMs) trimmed on both sides",
       ylab = "Difference in mean survival")
  dev.off()
  
  difference_survival_list = list()
  total_trimmed <- total[order(total$nonsynonymous),]
  for (i in 1:50)
  {
    ### Plot number of nonsynonymous mutations against benefit
    lower_quantile <- head(total_trimmed,i) #trim tail
    lower_quantile <- head(lower_quantile,-5)
    upper_quantile <- tail(total_trimmed,i) #trim head
    upper_quantile <- head(upper_quantile,-5)
    rbind(upper_quantile,lower_quantile)
    med_u = median(upper_quantile$overall_survival)
    med_l = median(lower_quantile$overall_survival)
    if (i<10)
    {difference_survival_list[i] <- NA}
    else
    {difference_survival_list[i] <- med_u - med_l}
  }
  namefile = paste(path,"Difference_mean_survival_quantiles_", which_data, ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  plot(as.numeric(difference_survival_list), type = "l", 
       xlab = "Number of patients in each extreme quantile",
       ylab = "Difference in mean survival")
  dev.off()
  
  difference_survival_list = list()
  total_trimmed <- total[order(total$nonsynonymous),]
  for (i in 1:50)
  {
    ### Plot number of nonsynonymous mutations against benefit
    lower_quantile <- head(total_trimmed,i) #trim tail
    lower_quantile <- head(total_trimmed,i) #trim tail
    upper_quantile <- tail(total_trimmed,i) #trim head
    rbind(upper_quantile,lower_quantile)
    med_u = mean(upper_quantile$overall_survival)
    med_l = mean(lower_quantile$overall_survival)
    if (i<10)
    {difference_survival_list[i] <- NA}
    else
    {difference_survival_list[i] <- med_u - med_l}
  }
  namefile = paste(path,"Difference_mean_survival_quantiles_", which_data, ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  plot(as.numeric(difference_survival_list), type = "l", 
       xlab = "Number of patients in each extreme quantile",
       ylab = "Difference in mean survival")
  dev.off()
  
  namefile = paste(path,"Survival_analysis_quantiles", which_data, ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  
  list_p_values = list()
  for (i in 1:50)
  {
    ### Plot number of nonsynonymous mutations against benefit
    lower_quantile <- head(total_trimmed,i) #trim tail
    upper_quantile <- tail(total_trimmed,i) #trim head
    total_subset <- rbind(upper_quantile,lower_quantile)
    real_test <- survdiff(Surv(total_subset$overall_survival, total_subset$dead) ~ total_subset$nonsynonymous > min(upper_quantile$nonsynonymous))
    list_p_values[i] <- round(1 - pchisq(real_test$chisq, length(real_test$n) - 1),4)
  }

  title_plot = paste("Log-rank p-value = ", p.val)
  total_Mut_over = upper_quantile
  mini.surv2 <- survfit(Surv(total_Mut_over$overall_survival , total_Mut_over$dead)~ 1, conf.type="none")
  plot(mini.surv2, mark.time = T, col = "red", main = which_data, xlab = "Time in days",
       ylab ="Proportion of survivors")
  total_Mut_under = lower_quantile
  mini.surv2 <- survfit(Surv(total_Mut_under$overall_survival , total_Mut_under$dead)~ 1, conf.type="none")
  lines(mini.surv2, mark.time = T, col = "blue")
  legend("topright", # places a legend at the appropriate place 
         c(under_string,over_string, title_plot), # puts text in the legend
         pch = c(".",".", ""), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),col=c("blue", "red", "white")) # gives the legend lines the correct color and width
  dev.off()
  rm(over_string, under_string, total_Mut_over, mini.surv2)
  }


#####################
# Survival Analysis #
#####################
if (params_tests[3] == T)
{
  namefile = paste(path,"Survival_analysis_", which_data, ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  real_test <- survdiff(Surv(total$overall_survival, total$dead) ~ total$nonsynonymous > params_values[2])
  p.val <- round(1 - pchisq(real_test$chisq, length(real_test$n) - 1),4)
  title_plot = paste("Log-rank p-value = ", p.val)
  over_string = paste(">",  params_values[2], " mutations")
  under_string = paste("<",  params_values[2], " mutations")
  total_Mut_over = total[which(total$nonsynonymous >  params_values[2]),] 
  mini.surv2 <- survfit(Surv(total_Mut_over$overall_survival , total_Mut_over$dead)~ 1, conf.type="none")
  plot(mini.surv2, mark.time = T, col = "red", main = which_data, xlab = "Time in days",
       ylab ="Proportion of survivors")
  total_Mut_under = total[which(total$nonsynonymous <=  params_values[2]),] 
  mini.surv2 <- survfit(Surv(total_Mut_under$overall_survival , total_Mut_under$dead)~ 1, conf.type="none")
  lines(mini.surv2, mark.time = T, col = "blue")
  legend("topright", # places a legend at the appropriate place 
         c(under_string,over_string, title_plot), # puts text in the legend
         pch = c(".",".", ""), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),col=c("blue", "red", "white")) # gives the legend lines the correct color and width
  dev.off()
  rm(over_string, under_string, total_Mut_over, mini.surv2)
}


#####################
# Permutation Tests #
#####################

if (params_tests[5] == T)
{
# Log-rank (benefit vs no benefit)
real_test <- survdiff(Surv(total$overall_survival, total$dead) ~ total$nonsynonymous > params_values[2])
real_Z <- as.numeric(as.character(((real_test$obs[1]-real_test$exp[1])^2)/real_test$exp[1] + ((real_test$obs[2]-real_test$exp[2])^2)/real_test$exp[2]))
list_i = list()
list_Z= vector()
j = 1
try(
  for (i in 1:params_values[1])
  {
    overall_survival_shuffled = sample(total$overall_survival)
    res <- survdiff(Surv(overall_survival_shuffled, total$dead) ~ total$nonsynonymous > params_values[2])
    Z = as.numeric(as.character(((res$obs[1]-res$exp[1]))/res$exp[1] + ((res$obs[2]-res$exp[2]))/res$exp[2]))
    list_i[j] = i
    list_Z[j] = Z
    j = j+1
  })
namefile = paste(path, "Permutation_Log_rank", which_data, ".tiff")
tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
h <- hist(list_Z, breaks = params_values[1]/50, main = which_data, xlab = "Z statisitic", 
          xlim = c(min(list_Z),max(real_Z + real_Z/20, max(list_Z))))
#xfit<-seq(min(list_Z),max(real_Z + real_Z/20, max(list_Z)),length=40) 
#yfit<-dnorm(xfit,mean=mean(list_Z),sd=sd(list_Z)) 
#yfit <- yfit*diff(h$mids[1:2])*length(list_Z) 
#lines(xfit, yfit, col="blue", lwd=2)
abline(v= real_Z, col = 'red')
pvalue_comp = length(which(list_Z > real_Z))/params_values[1]
pvalue = paste("p-value = ",pvalue_comp)
legend("topright", # places a legend at the appropriate place 
       c(pvalue), # puts text in the legend
       pch = c(".")) # gives the legend lines the correct color and width
dev.off()
rm(list_i ,list_Z, real_test, real_Z, pvalue_comp, overall_survival_shuffled, res, Z, pvalue)
}


















######################
# Threshold Mutation #
######################
if (threshold_mutation == T)
{
  ## P_value vs threshold of mutations
  par(mfrow=c(1,1))
  list_i = list()
  list_pvalues = list()
  j = 1
  try(
    for (i in min(na.omit(nonsynonymous)):max(na.omit(nonsynonymous))) 
    {
      res <- survdiff(Surv(overall_survival*12, dead) ~ nonsynonymous > i)
      pvalue = pchisq(res$chisq, length(res$n)-1, lower.tail = FALSE)
      list_i[j] = i
      list_pvalues[j] = pvalue
      j = j+1
    })
  namefile3 = paste("Pvalue_VS_MutationThreshold_", which_data, ".jpg")
  jpeg(namefile3)
  plot(list_i, list_pvalues, pch = ".", col = 'purple', cex = 5, xlab = "Discriminant number of mutations",
       ylab = "p-values", main = which_data)
  abline(a=0.05, b=0, col = "purple")
  rm(list_i, list_pvalues, j, res)
  dev.off()
}

#############################
# Optimal cutpoint analysis #
#############################

if (optimal_cutpoint == T)
{
## P_value vs threshold of mutations
par(mfrow=c(1,1))
list_j = list()
list_pvalues = list()
j = 1
try(
  for (i in sort(unique(na.omit(nonsynonymous)))) 
  {
    res <- survdiff(Surv(overall_survival*12, dead) ~ nonsynonymous > i)
    pvalue = pchisq(res$chisq, length(res$n)-1, lower.tail = FALSE)
    list_j[j] = i
    list_pvalues[j] = pvalue
    j = j+1
  })
    namefile7 = paste("Smoothing_spline_", which_data, ".jpg")
    jpeg(namefile7)
    plot(list_j, list_pvalues, pch = ".", col = 'purple', cex = 5, xlab = "Discriminant number of mutations",
         ylab = "p-values", main = which_data)
    spline = smooth.spline(list_j, list_pvalues, df = 5)
    lines(spline, lwd = 5, col = 'red')
    dev.off()
    rm(list_j, list_pvalues, res, pvalue)
}

###########################
# Stage severity analysis #
###########################

if (stage_severity == T && which_data != "Rizvi")
{
  namefile8 = paste("Stage_severity_", which_data, ".jpg")
  jpeg(namefile8)
  plot(total_stage1$nonsynonymous, total_stage1$overall_survival, pch = ".", col = 'yellow', cex = 10,
       xlim = c(min(na.omit(nonsynonymous)), max(na.omit(nonsynonymous))), 
       ylim = c(min(na.omit(overall_survival)), max(na.omit(overall_survival))),
       xlab = "Number of non synonymous mutations", ylab = "Overall survival", main = which_data)
  points(total_stage2$nonsynonymous, total_stage2$overall_survival, pch = ".", col = 'orange', cex = 10)
  points(total_stage3$nonsynonymous, total_stage3$overall_survival, pch = ".", col = 'red', cex = 10)
  points(total_stage4$nonsynonymous, total_stage4$overall_survival, pch = ".", col = 'black', cex = 10)
  legend("topright", # places a legend at the appropriate place 
         c(stage1_str,stage2_str,stage3_str,stage4_str), # puts text in the legend
         pch = c(".",".",".","."), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5,2.5,2.5),col=c("yellow", "orange","red", "black")) # gives the legend lines the correct color and width
  dev.off()
  
  # Same in log scale
  namefile9 = paste("Stage_severity_", which_data, ".jpg")
  jpeg(namefile9)
  plot(log(total_stage1$nonsynonymous), log(total_stage1$overall_survival), pch = ".", col = 'yellow', cex = 10,
       xlim = c(min(na.omit(log(nonsynonymous))), max(na.omit(log(nonsynonymous)))), 
       ylim = c(min(na.omit(log(overall_survival))), max(na.omit(log(overall_survival)))),
       xlab = "Log number of non synonymous mutations", ylab = "Log overall survival", main = which_data)
  points(log(total_stage2$nonsynonymous), log(total_stage2$overall_survival), pch = ".", col = 'orange', cex = 10)
  points(log(total_stage3$nonsynonymous), log(total_stage3$overall_survival), pch = ".", col = 'red', cex = 10)
  points(log(total_stage4$nonsynonymous), log(total_stage4$overall_survival), pch = ".", col = 'black', cex = 10)
  legend("topright", # places a legend at the appropriate place 
         c(stage1_str,stage2_str,stage3_str,stage4_str), # puts text in the legend
         pch = c(".",".",".","."), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5,2.5,2.5),col=c("yellow", "orange","red", "black")) # gives the legend lines the correct color and width
  dev.off()
}

if (gender_analysis == T)
{
  namefile8 = paste("gender_", which_data, ".jpg")
  jpeg(namefile8)
  plot(total_male$nonsynonymous, total_male$overall_survival, pch = ".", col = 'blue', cex = 10,
       xlim = c(min(na.omit(nonsynonymous)), max(na.omit(nonsynonymous))), 
       ylim = c(min(na.omit(overall_survival)), max(na.omit(overall_survival))),
       xlab = "Number of non synonymous mutations", ylab = "Overall survival", main = which_data)
  points(total_female$nonsynonymous, total_female$overall_survival, pch = ".", col = 'pink', cex = 10)
  legend("topright", # places a legend at the appropriate place 
         c("Male", "Female"), # puts text in the legend
         pch = c(".","."), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),col=c("blue", "pink")) # gives the legend lines the correct color and width
  dev.off()
  
  # Same in log scale
  namefile9 = paste("Log_gender_", which_data, ".jpg")
  jpeg(namefile9)
  plot(log(total_male$nonsynonymous), log(total_male$overall_survival), pch = ".", col = 'blue', cex = 10,
       xlim = c(min(na.omit(log(nonsynonymous))), max(na.omit(log(nonsynonymous)))), 
       ylim = c(min(na.omit(log(overall_survival))), max(na.omit(log(overall_survival)))),
       xlab = "Log number of non synonymous mutations", ylab = "Log overall survival", main = which_data)
  points(log(total_female$nonsynonymous), log(total_female$overall_survival), pch = ".", col = 'pink', cex = 10)
  legend("topright", # places a legend at the appropriate place 
         c("Male", "Female"), # puts text in the legend
         pch = c(".","."), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),col=c("blue", "pink")) # gives the legend lines the correct color and width
  dev.off()
}
rm(total)

#################
# Log odd ratio #
#################


###########################
# Stage severity analysis #
###########################

if (pitfalls == T)
{
# Test if male and females have different mutational load
wilcox.test(total_male$nonsynonymous, total_female$nonsynonymous)
# Test if male and females have different overall survival
wilcox.test(total_male$overall_survival, total_female$overall_survival)

# Test if responders and nonresponders have different mutational load
wilcox.test(total_benefit$nonsynonymous, total_nobenefit$nonsynonymous)
# Test if male and females have different overall survival
wilcox.test(total_benefit$overall_survival, total_nobenefit$overall_survival)

# Test if responders and nonresponders have different ages
wilcox.test(total_benefit$age, total_nobenefit$age)

if (which_data != "Rizvi")
  {
  # Test if stage 1 and stage 4 patients have different mutational load
  wilcox.test(total_stage1$nonsynonymous, total_stage4$nonsynonymous)
  # Test if stage 1 and stage 4 patients have different overall survival
  wilcox.test(total_stage1$overall_survival, total_stage4$overall_survival)

  # Test if stage 1 and stage 3 patients have different mutational load
  wilcox.test(total_stage1$nonsynonymous, total_stage3$nonsynonymous)
  # Test if stage 1 and stage 3 patients have different mutational load
  wilcox.test(total_stage1$overall_survival, total_stage3$overall_survival)

  # Test if stage 1 and stage 2 patients have different mutational load
  wilcox.test(total_stage1$nonsynonymous, total_stage2$nonsynonymous)
  # Test if stage 1 and stage 2 patients have different mutational load
  wilcox.test(total_stage1$overall_survival, total_stage2$overall_survival)
  
  # Test if stage 2 and stage 3 patients have different mutational load
  wilcox.test(total_stage2$nonsynonymous, total_stage3$nonsynonymous)
  # Test if stage 2 and stage 3 patients have different mutational load
  wilcox.test(total_stage2$overall_survival, total_stage3$overall_survival)
  
  # Test if stage 2 and stage 4 patients have different mutational load
  wilcox.test(total_stage2$nonsynonymous, total_stage4$nonsynonymous)
  # Test if stage 2 and stage 4 patients have different mutational load
  wilcox.test(total_stage2$overall_survival, total_stage4$overall_survival)
  
  # Test if stage 3 and stage 4 patients have different mutational load
  wilcox.test(total_stage3$nonsynonymous, total_stage4$nonsynonymous)
  # Test if stage 3 and stage 4 patients have different mutational load
  wilcox.test(total_stage3$overall_survival, total_stage4$overall_survival)
}
}
Done_string = paste(which_data, " done!")
print(Done_string)


if(FALSE) {
if (burden_VS_survival_and_age == T)
{
  par(mfrow=c(1,1))
  namefile = paste("Survival_VS_Mutational_load_and_age", which_data, ".tiff")
  tiff(namefile, width = 12, height = 8, units = 'in', res = 500)
  plot(total_benefit$nonsynonymous, total_benefit$overall_survival, xlab= "Number of Nonsynonymous mutations",
       ylab= OS_or_PFS, col ="blue", pch = ".", cex = total_benefit$age*0.1, main = which_data,      
       xlim = c(min(na.omit(nonsynonymous)),max(na.omit(nonsynonymous))), 
       ylim= c(min(na.omit(overall_survival)),max(na.omit(overall_survival))))
  points(total_nobenefit$nonsynonymous, total_nobenefit$overall_survival, 
         cex = total_nobenefit$age*0.1, pch = ".", col = "red")
  legend("topright", # places a legend at the appropriate place 
         c("Responders","Non-responders"), # puts text in the legend
         pch = c(".","."),
         lwd = c(1,1),
         col=c("blue", "red")) # gives the legend lines the correct color and width
  dev.off()
}
par(mfrow=c(1,1))
namefile2 = paste("Benefit_VS_MutationalLoad_", which_data, ".jpg")
jpeg(namefile2)
plot(sort(na.omit(total_benefit$nonsynonymous)), col = "blue", xlab = "Rank", 
     ylab = "Number of Nonsynonymous mutations", pch = ".",
     cex = 5, main = which_data)
legend("topright", # places a legend at the appropriate place 
       c("Responders","Non-responders"), # puts text in the legend
       pch = c(".","."), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),col=c("blue", "red")) # gives the legend lines the correct color and width
abline(a=mean(na.omit(total_benefit$nonsynonymous)), b=0, col = "blue")
points(sort(na.omit(total_nobenefit$nonsynonymous)), col = "red", xlab = "", 
       pch = ".", cex = 5)
abline(a=mean(na.omit(total_nobenefit$nonsynonymous)), b=0, col = "red")
dev.off() }


###############
# Random data #
###############
par(mfrow=c(1, 2))
# create dataset with 1000 normally distributed points
df <- data.frame(x = rnorm(100, mean(log(total$nonsynonymous)), sd(log(total$nonsynonymous))),
                 y = rpois(100,mean(log(total$overall_survival))))
# create a ggplot2 scatterplot
p <- ggplot(df, aes(x, y)) + geom_point() + theme_classic()
# add marginal histograms
ggExtra::ggMarginal(p, type = "histogram")

# create dataset with 1000 normally distributed points
df <- data.frame(x = log(total$nonsynonymous), y = (total$overall_survival))
# create a ggplot2 scatterplot
p <- ggplot(df, aes(x, y)) + geom_point() + theme_classic()
# add marginal histograms
ggExtra::ggMarginal(p, type = "histogram")


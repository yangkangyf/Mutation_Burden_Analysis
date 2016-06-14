#########################
# Adjustable parameters #
#########################
# Number of simulations for the permutation test
number_of_simulations = 10000
# Discriminant number of mutation to perform the survival analysis
mutation_discriminant = 100
# Choose which data to analyse
Van_Allen = F
Snyder = F
Rizvi = F
Hugo = T

# Choose which test to perform
replace_eight_patients = F
burden_VS_survival_and_age = T
burden_VS_survival = T
burden_VS_benefit = T
threshold_mutation = T
survival_analysis = T
permutation_tests_LR = T
permutation_tests_MW = T
permutation_tests_LR_multiple_mut = T

#######################
# Librairies and data #
#######################
# Load libraries
library(survival)

# Run the analysis for selected data
selected_data = list()
if (Van_Allen == T)
  selected_data[length(selected_data)+1]="Van_Allen"
if (Rizvi == T)
  selected_data[length(selected_data)+1]="Rizvi"
if (Snyder == T)
  selected_data[length(selected_data)+1]="Snyder"
if (Hugo == T)
  selected_data[length(selected_data)+1]="Hugo"

for (which_data in selected_data)
{
# Import and rearrange Hugo data
if (which_data == "Hugo")
  {
    total1= read.table("Hugo1.txt", header = T)
    total2 = read.table("Hugo2.txt", header = T, sep = "\t")
    total <- merge(total1, total2, by.y = "Patient_ID")
    total$nonsynonymous <- total$TotalNonSyn
    total$group <- total$Response
    clinical_benefit_str <- "R"
    clinical_nobenefit_str <- "NR"
    OS_or_PFS <- "Overall survival"
    total$age <- total$Age
    total$overall_survival <- total$Overall_Survival
    total$dead <- (total$Vital_Status == "Dead")*1
  }
# Import and rearrange Van Allen data
if (which_data == "Van_Allen")
{
  total = read.csv("Van_Allen.csv")
  total$nonsynonymous <- total$nonsynonymous
  total$overall_survival <- total$overall_survival
  total$group <- total$group
  clinical_benefit_str <- "response"
  clinical_nobenefit_str <- "nonresponse"
  total$dead <- total$dead
  OS_or_PFS <- "Overall survival"
  total$age <- total$age_start
}
# Import and rearrange RIZVI data
if (which_data == "Rizvi")
{
total = read.csv("Rizvi.csv")
total$nonsynonymous <- total$Nonsyn.
total$overall_survival <- total$PFS..mos.*365
total$group <- total$Durable.Clinical.Benefit
clinical_benefit_str <- "DCB"
clinical_nobenefit_str <- "NDB"
total$dead <- total$Event....
OS_or_PFS <- "Progression free survival"
total$age <- total$Age..years.
}
# Import and rearrange SNYDER data
if (which_data == "Snyder")
{
S1 = read.table("Snyder1.txt", comment.char = "%", header = T)
S2 = read.table("Snyder2.txt", comment.char = "%", header = T)
S3 = read.table("Snyder3.txt", comment.char = "%", header = T)
bound_S1S2 <- rbind(S1,S2)
total <- merge(bound_S1S2, S3, by.y = "Study_ID")
total$nonsynonymous <- total$Mutation
total$overall_survival <- total$OS.yr. *365
clinical_benefit_str <- "1"
clinical_nobenefit_str <- "0"
total$dead <- total$Alive_at_time_of_censure
eightpatient <- rbind(total[48,],total[22,],total[23,],total[24,],
                      total[25,],total[13,],total[14,],total[15,])
OS_or_PFS <- "Overall survival"
total$age <- total$Age
if (replace_eight_patients == T)
  {
  total$Benefit[48] = 1
  total$Benefit[13] = 1
  total$Benefit[14] = 1
  total$Benefit[15] = 1
  total$Benefit[22] = 1
  total$Benefit[23] = 1
  total$Benefit[24] = 1
  total$Benefit[25] = 1
  }
total$group <- total$Benefit
}
  
total_benefit = total[which(total$group == clinical_benefit_str),]
total_nobenefit = total[which(total$group == clinical_nobenefit_str),]
if (which_data == "Snyder")
{
  eightpatient_benefit = eightpatient[which(eightpatient$group == clinical_benefit_str),]
  eightpatient_nobenefit = eightpatient[which(eightpatient$group == clinical_nobenefit_str),]
}
  
###############
# Basic plots #
###############

if (burden_VS_survival_and_age == T)
{
  par(mfrow=c(1,1))
  namefile2 = paste("Survival_VS_Mutational_load_and_age", which_data, ".jpg")
  jpeg(namefile2)
  plot(total_benefit$nonsynonymous, total_benefit$overall_survival, xlab= "Number of Nonsynonymous mutations",
        ylab= OS_or_PFS, col ="blue", pch = ".", cex = total_benefit$age*0.1, main = which_data,      
        xlim = c(min(na.omit(total$nonsynonymous)),max(na.omit(total$nonsynonymous))), 
        ylim= c(min(na.omit(total$overall_survival)),max(na.omit(total$overall_survival))))
  points(total_nobenefit$nonsynonymous, total_nobenefit$overall_survival, 
         cex = total_nobenefit$age*0.1, pch = ".", col = "red")
  legend("topright", # places a legend at the appropriate place 
         c("Responders","Non-responders"), # puts text in the legend
         pch = c(".","."),
         lwd = c(1,1),
         col=c("blue", "red")) # gives the legend lines the correct color and width
  dev.off()
}

if (burden_VS_survival==T)
{
# Plot Survival against number of nonsynonymous mutations
par(mfrow=c(1,1))
namefile1 = paste("Survival_VS_Mutational_load_", which_data, ".jpg")
jpeg(namefile1)

plot( total_benefit$nonsynonymous, total_benefit$overall_survival, xlab= "Number of Nonsynonymous mutations",
     ylab= OS_or_PFS, col ="blue", pch = ".", cex = 5, main = which_data,      
     xlim = c(min(na.omit(total$nonsynonymous)),max(na.omit(total$nonsynonymous))), 
     ylim= c(min(na.omit(total$overall_survival)),max(na.omit(total$overall_survival))))
points(total_nobenefit$nonsynonymous, total_nobenefit$overall_survival, col ="red",
  pch = ".", cex = 5)
legend("topright", # places a legend at the appropriate place 
       c("Responders","Non-responders"), # puts text in the legend
       pch = c(".","."),
       lwd = c(1,1),
       col=c("blue", "red")) # gives the legend lines the correct color and width
if (which_data == "Snyder")
{
  points(eightpatient$nonsynonymous, eightpatient$overall_survival, col ="orange", cex= 3)
  legend("topright", # places a legend at the appropriate place 
         c("Responders","Non-responders", "Eight patients"), # puts text in the legend
         pch = c(".",".", "."), # gives the legend appropriate symbols (lines),
         lwd = c(1,1,1),
         col=c("blue", "red", "orange")) # gives the legend lines the correct color and width
}
dev.off()

# Same figure, but in log scale
namefile1bis = paste("Survival_VS_Mutational_load_log", which_data, ".jpg")
jpeg(namefile1bis)
plot(log(total_benefit$nonsynonymous), log(total_benefit$overall_survival), xlab= "Number of Nonsynonymous mutations",
      ylab= "Overall survival", col ="blue", pch = ".", cex = 5, main = which_data,
      xlim = c(min(log(na.omit(total$nonsynonymous))),max(log(na.omit(total$nonsynonymous)))), 
      ylim= c(min(log(na.omit(total$overall_survival))),max(log(na.omit(total$overall_survival)))))
points(log(na.omit(total_nobenefit$nonsynonymous)), log(na.omit(total_nobenefit$overall_survival)),  xlab= "Number of Nonsynonymous mutations",
       ylab= OS_or_PFS, col ="red", pch = ".", cex = 5, main = which_data)
legend("topright", # places a legend at the appropriate place 
       c("Responders","Non-responders"), # puts text in the legend
       pch = c(".","."), # gives the legend appropriate symbols (lines),
       lwd = c(1,1),
       col=c("blue", "red")) # gives the legend lines the correct color and width
if (which_data == "Snyder")
{
  points(log(na.omit(eightpatient$nonsynonymous)), log(na.omit(eightpatient$overall_survival)), col ="orange", cex= 3)
  legend("topright", # places a legend at the appropriate place 
         c("Responders","Non-responders", "Eight patients"), # puts text in the legend
         pch = c(".",".", "."), # gives the legend appropriate symbols (lines),
         lwd = c(1,1,1),
         col=c("blue", "red", "orange")) # gives the legend lines the correct color and width
}
dev.off()
}

##############################################################
# Association between Mutational Burden and Clinical Benefit #
##############################################################
if (burden_VS_benefit == T)
{
### Plot number of nonsynonymous mutations against benefit
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
dev.off()
}

#####################
# Survival Analysis #
#####################
if (survival_analysis == T)
{
  namefile2bis = paste("Survival_Anaysis_", which_data, ".jpg", ".jpg")
  jpeg(namefile2bis)
  over_string = paste(">", mutation_discriminant, " mutations")
  under_string = paste("<", mutation_discriminant, " mutations")
  total_Mut_over = total[which(na.omit(total$nonsynonymous) > mutation_discriminant),] 
  mini.surv2 <- survfit(Surv(total_Mut_over$overall_survival , total_Mut_over$dead)~ 1, conf.type="none")
  plot(mini.surv2, mark.time = T, col = "red", main = which_data, xlab = "Time in days",
       ylab ="Proportion of survivors")
  # Survival in Discovery set with <= 100 mutations
  total_Mut_under = total[which(na.omit(total$nonsynonymous) <= mutation_discriminant),] 
  mini.surv2 <- survfit(Surv(na.omit(total_Mut_under$overall_survival) , na.omit(total_Mut_under$dead))~ 1, conf.type="none")
  lines(mini.surv2, mark.time = T, col = "blue")
  legend("topright", # places a legend at the appropriate place 
         c(under_string,over_string), # puts text in the legend
         pch = c(".","."), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),col=c("blue", "red")) # gives the legend lines the correct color and width
  dev.off()
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
for (i in min(na.omit(total$nonsynonymous)):max(na.omit(total$nonsynonymous))) 
{
  res <- survdiff(Surv(total$overall_survival*12, total$dead) ~ total$nonsynonymous > i)
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
dev.off()
}

#####################
# Permutation Tests #
#####################
if (permutation_tests_MW == T)
{
# Mann_Whitney (benefit vs no benefit)
real_test <- wilcox.test(na.omit(total_benefit$nonsynonymous), na.omit(total_nobenefit$nonsynonymous))
real_W <- as.numeric(as.character(real_test$statistic))
list_i = list()
list_W= vector()
j = 1
try(
for (i in 1:number_of_simulations)
{
  all_nonsynonymous <- append(na.omit(total_benefit$nonsynonymous), na.omit(total_nobenefit$nonsynonymous))
  all_nonsynonymous = sample(all_nonsynonymous)
  total_benefit_shuffled = all_nonsynonymous[1:length(na.omit(total_benefit$nonsynonymous))]
  total_nobenefit_shuffled = na.omit(all_nonsynonymous[length(na.omit(total_benefit$nonsynonymous))+1:length(all_nonsynonymous)])
  res <- wilcox.test(total_benefit_shuffled, total_nobenefit_shuffled)
  W = as.numeric(as.character(res$statistic))
  list_i[j] = i
  list_W[j] = W
  j = j+1
})
namefile4 = paste("Permutation_Mann-Whitney_", which_data, ".jpg")
jpeg(namefile4)
h <- hist(list_W, breaks = number_of_simulations/50, main = which_data, xlab = "W statisitics", 
     xlim = c(min(list_W),max(real_W + real_W/20, max(list_W))))
#xfit<-seq(min(list_W),max(real_W + real_W/20, max(list_W)),length=40) 
#yfit<-dnorm(xfit,mean=mean(list_W),sd=sd(list_W)) 
#yfit <- yfit*diff(h$mids[1:2])*length(list_W) 
#lines(xfit, yfit, col="blue", lwd=2)
abline(v= real_W, col = 'red')
pvalue_comp = length(which(list_W > real_W))/number_of_simulations
pvalue = paste("p-value = ",pvalue_comp)
legend("topright", # places a legend at the appropriate place 
       c(pvalue), # puts text in the legend
       pch = c(".")) # gives the legend lines the correct color and width
dev.off()
}

if (permutation_tests_LR == T)
{
# Log-rank (benefit vs no benefit)
real_test <- survdiff(Surv(total$overall_survival, total$dead) ~ total$nonsynonymous > mutation_discriminant)
real_Z <- as.numeric(as.character(((real_test$obs[1]-real_test$exp[1])^2)/real_test$exp[1] + ((real_test$obs[2]-real_test$exp[2])^2)/real_test$exp[2]))
list_i = list()
list_Z= vector()
j = 1
try(
  for (i in 1:number_of_simulations)
  {
    overall_survival_shuffled = sample(total$overall_survival)
    res <- survdiff(Surv(overall_survival_shuffled, total$dead) ~ total$nonsynonymous > mutation_discriminant)
    Z = as.numeric(as.character(((res$obs[1]-res$exp[1]))/res$exp[1] + ((res$obs[2]-res$exp[2]))/res$exp[2]))
    list_i[j] = i
    list_Z[j] = Z
    j = j+1
  })
namefile5 = paste("Permutation_Log-Rank_",which_data, ".jpg")
jpeg(namefile5)
h <- hist(list_Z, breaks = number_of_simulations/50, main = which_data, xlab = "Z statisitics", 
          xlim = c(min(list_Z),max(real_Z + real_Z/20, max(list_Z))))
#xfit<-seq(min(list_Z),max(real_Z + real_Z/20, max(list_Z)),length=40) 
#yfit<-dnorm(xfit,mean=mean(list_Z),sd=sd(list_Z)) 
#yfit <- yfit*diff(h$mids[1:2])*length(list_Z) 
#lines(xfit, yfit, col="blue", lwd=2)
abline(v= real_Z, col = 'red')
pvalue_comp = length(which(list_Z > real_Z))/number_of_simulations
pvalue = paste("p-value = ",pvalue_comp)
legend("topright", # places a legend at the appropriate place 
       c(pvalue), # puts text in the legend
       pch = c(".")) # gives the legend lines the correct color and width
dev.off()
}

if (permutation_tests_LR_multiple_mut == T)
{
  pvalue_comp = vector()
  mut_chosen_list = vector()
  k = 0
  try(
  for (mut_chosen in (min(na.omit(total$nonsynonymous))):(max(na.omit(total$nonsynonymous))))
    {
    real_test <- survdiff(Surv(total$overall_survival, total$dead) ~ total$nonsynonymous > mut_chosen)
    real_Z <- as.numeric(as.character(((real_test$obs[1]-real_test$exp[1])^2)/real_test$exp[1] + ((real_test$obs[2]-real_test$exp[2])^2)/real_test$exp[2]))
    list_i = list()
    list_Z= vector()
    j = 1
      for (i in 1:number_of_simulations)
      {
        overall_survival_shuffled = sample(total$overall_survival)
        res <- survdiff(Surv(overall_survival_shuffled, total$dead) ~ total$nonsynonymous > mut_chosen)
        Z = as.numeric(as.character(((res$obs[1]-res$exp[1]))/res$exp[1] + ((res$obs[2]-res$exp[2]))/res$exp[2]))
        list_i[j] = i
        list_Z[j] = Z
        j = j+1
      }
    pvalue_comp[k] = length(which(list_Z > real_Z))/number_of_simulations
    mut_chosen_list[k] = mut_chosen
    k = k+1
  })
  namefile6 = paste("Permutation_Log-Rank_Multiple",which_data, ".jpg")
  jpeg(namefile6)
  plot(mut_chosen_list, pvalue_comp, xlim = c(min(na.omit(total$nonsynonymous)),max(na.omit(total$nonsynonymous))), 
       ylim = c(min(pvalue_comp), max(pvalue_comp)), pch = ".",
       xlab = "Discriminant number of mutations", ylab = "p_value")
  lines(mut_chosen_list, pvalue_comp)
  dev.off()
}
}
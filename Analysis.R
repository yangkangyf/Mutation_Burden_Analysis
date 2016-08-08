#######################
# Librairies and data #
#######################
# Load libraries
library(survival)
library(ggplot2)
library(gdata)

#########################
# Adjustable parameters #
#########################
params_values = list()
# Number of simulations for the permutation test
params_values[1] = 10000
# Discriminant number of mutation to perform the survival analysis
params_values[2] = 102
# Redefine clinical classification
redef = F

###############
# Run program #
###############
# Import all datasets
imported = import()
total = as.data.frame(imported[1])
special_patients = as.numeric(unlist(imported[2]))
# Choose clinical classification (0 for original data, 1 for classfcation according to OS)
if (redef == T)
{
  group_temp = total$group
  total$group = total$group2
  total$group2 = group_temp
  path = "~/Figures/Group_1/"
} else
{
  path = "~/Figures/Group_0/"
}
total_hugo = total[which(total$dataset == "mel1"),]
total_snyder = total[which(total$dataset == "mel2"),]
total_vanallen = total[which(total$dataset == "mel3"),]
total_rizvi = total[which(total$dataset == "nslc"),]
# Plot mutational load VS survival
plot_nonlog(list(total_hugo, total_snyder, total_vanallen, total_rizvi))
plot_log(list(total_hugo, total_snyder, total_vanallen, total_rizvi))
# Boxplots
boxplot_MW()
# Permutation MW
permutation_MW(list(total_hugo, total_snyder, total_vanallen))
# Compare effect to random distribution
rand_comp(list(total_hugo, total_snyder, total_vanallen))
# Contingency table to compare with random normal distribution
cont_table(list(total_hugo, total_snyder, total_vanallen))
# To plot random data
rand_plot(list(total_hugo, total_snyder, total_vanallen))
# Analysis by trimming data
trim_analysis(list(total_hugo, total_snyder, total_vanallen))
# Correlation plots
corr_plots(list(total_hugo, total_snyder, total_vanallen))



########################
# Functions to compile #
########################
# Import data from mel1, mel2, mel3 and nslc
import <- function()
{
total <- data.frame(
  ticker=character(),
  value=numeric(),
  date = as.Date(character()),
  stringsAsFactors=FALSE
)
for (which_data in c("mel1","mel2", "mel3", "nslc"))
{
params_data = list()
# Import and rearrange Van Allen data
if (which_data == "mel3")
{
  data = read.csv("Van_Allen.csv")
  nonsynonymous <- data$nonsynonymous
  overall_survival <- data$overall_survival
  group <- data$group
  pre_th <- data$pre_therapies
  params_data[1] <- "response"
  params_data[2] <- "nonresponse"
  dead <- data$dead
  age <- data$age_start
  stage <- data$M
  params_data[4] = "M0"
  params_data[5] = "M1a"
  params_data[6] = "M1b"
  params_data[7] = "M1c"
  gender <- data$gender
  params_data[8] = "male"
  params_data[9] = "female"
  neoantigens = data$neos500
  special_patients <- which(data$overall_survival > 2*365 & data$progression_free < 6*30 & data$RECIST == "PD")
  rm(data)
}
  
# Import and rearrange SNYDER data
if (which_data == "mel2")
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
if (which_data == "nslc")
{
  data = read.csv("Rizvi.csv")
  nonsynonymous <- data$Nonsyn.
  overall_survival <- data$PFS..mos.*365
  group <- data$Durable.Clinical.Benefit
  params_data[1] <- "DCB"
  params_data[2] <- "NDB"
  dead <- data$Event....
  age <- data$Age..years.
  gender <- data$Sex
  params_data[8]  = "M"
  params_data[9]  = "F"
  stage <- replicate(length(as.numeric(nonsynonymous)), "IV")
  rm(data)
}
  
  
# Import and rearrange Hugo data
if (which_data == "mel1")
  {
    data1= read.table("Hugo1.txt", header = T)
    data2 = read.table("Hugo2.txt", header = T, sep = "\t")
    data <- merge(data1, data2, by.y = "Patient_ID")
    nonsynonymous <- data$TotalNonSyn
    overall_survival <- data$Overall_Survival
    group <- data$Response
    params_data[1] <- "R"
    params_data[2] <- "NR"
    "Overall Survival" <- "Overall survival"
    age <- data$Age
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

params_data = unlist(params_data)

group2 = group
group2 = as.character(group2) 
group = as.character(group) 
group[which(group == as.character(params_data[1]))] = "Responders"
group[which(group == as.character(params_data[2]))] = "Nonresponders"
if (which_data == "nslc")
{
group2[which(as.numeric((overall_survival)) >= 2*365)] = "Responders"
group2[which(as.numeric((overall_survival)) < 2*365)] = "Nonresponders"
} else {
group2[which(as.numeric((overall_survival)) >= 365)] = "Responders"
group2[which(as.numeric((overall_survival)) < 365)] = "Nonresponders"
}


total_temp = as.data.frame(cbind(as.numeric((nonsynonymous)), 
                            as.character(group), 
                            as.numeric(age),
                            as.numeric(overall_survival),
                            as.numeric(dead), 
                            as.character(gender), 
                            as.character(stage),
                            as.character(replicate(length(nonsynonymous), which_data)),
                            as.character(group2)),
                            as.is = F,
                            stringsAsFactors = FALSE)

total_temp <- drop.levels(total_temp)
total_temp$V1 <- as.numeric(total_temp$V1)
total_temp$V3 <- as.numeric(total_temp$V3)
total_temp$V4 <- as.numeric(total_temp$V4)
total_temp$V5 <- as.numeric(total_temp$V5)
names(total_temp)[1] <- "nonsynonymous"
names(total_temp)[2] <- "group"
names(total_temp)[3] <- "age"
names(total_temp)[4] <- "overall_survival"
names(total_temp)[5] <- "dead"
names(total_temp)[6] <- "gender"
names(total_temp)[7] <- "stage"
names(total_temp)[8] <- "dataset"
names(total_temp)[9] <- "group2"
total <- rbind(total_temp, total)
}
return(list(na.omit(total), special_patients))
}


###############
# Basic plots #
###############
# Plot survival against number of nonsynonymous mutation

plot_nonlog <- function(list_arg)
  {
  for (data in list_arg)
  {
    par(mfrow=c(1,1))
    which_data = data$dataset[1]
    namefile = paste(path, "survVSmut_", which_data, ".tiff")
    tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
    plot(data[which(data$group == "Responders"),]$nonsynonymous, 
         data[which(data$group == "Responders"),]$overall_survival,
         xlab= "Number of nonsynonymous mutations",
         ylab= "Overall Survival", pch = 16,col = "cyan3", cex = 1, #main = which_data,      
         xlim = c(min(data$nonsynonymous),max(data$nonsynonymous)), 
         ylim= c(min(na.omit(data$overall_survival)),max(na.omit(data$overall_survival))))

    if (which_data = "nslc")
    {
      abline(h=2*365, lty = 2, col = "purple")
      text(max(data$nonsynonymous)/2, 365*2+55, "OS = 2 year", cex = 0.9, col = "purple")
    }
    else
    {
      abline(h=365, lty = 2, col = "purple")
      text(max(data$nonsynonymous)/2, 420, "OS = 1 year", cex = 0.9, col = "purple")
    }
    points(data[which(data$group == "Responders"),][which(data[which(data$group == "Responders"),]$dead == 1),]$nonsynonymous, 
           data[which(data$group == "Responders"),][which(data[which(data$group == "Responders"),]$dead == 1),]$overall_survival,
           pch = 1, cex = 1)
    points(data[which(data$group == "Nonresponders"),]$nonsynonymous, 
           data[which(data$group == "Nonresponders"),]$overall_survival, col ="red",
           pch = 16, cex = 1)
    points(data[which(data$group == "Nonresponders"),][which(data[which(data$group == "Nonresponders"),]$dead == 1),]$nonsynonymous, 
           data[which(data$group == "Nonresponders"),][which(data[which(data$group == "Nonresponders"),]$dead == 1),]$overall_survival, xlab= "Number of Nonsynonymous mutations",
           ylab= "Overall Survival", pch = 1, cex = 1)
    #legend("topright", # places a legend at the appropriate place 
    #       c("Responders","Non-responders", "Dead"), # puts text in the legend
    #       pch = c(16,16,1),
    #       col=c("cyan3", "red", "black")) # gives the legend lines the correct color and width
    if (which_data == "mel2"|| which_data == "Snyder2")
    {
      points(data[special_patients,]$nonsynonymous, data[special_patients,]$overall_survival, 
             col ="grey", cex= 1, pch = 16)
      points(data[special_patients,][which(data[special_patients,]$dead == 1),]$nonsynonymous, 
             data[special_patients,][which(data[special_patients,]$dead == 1),]$overall_survival
             , pch = 1, cex = 1)
      #legend("topright", # places a legend at the appropriate place 
      #       c("Responders","Non-responders","Dead", "Reassigned patients"), # puts text in the legend
      #       pch = c(16,16,1,16), # gives the legend appropriate symbols (lines),
      #       col=c("cyan", "red","black", "grey")) # gives the legend lines the correct color and width
    }
    if (which_data == "mel3")
    {
      points(data[special_patients,]$nonsynonymous, data[special_patients,]$overall_survival, 
             col ="grey", cex= 1, pch = 16)
      points(data[special_patients,][which(data[special_patients,]$dead == 1),]$nonsynonymous, 
             data[special_patients,][which(data[special_patients,]$dead == 1),]$overall_survival
             , pch = 1, cex = 1)
      #  legend("topright", # places a legend at the appropriate place 
      #         c("Responders","Non-responders","Dead", "Subset"), # puts text in the legend
      #         pch = c(16,16,1,16), # gives the legend appropriate symbols (lines),
      #         col=c("cyan", "red","black", "grey")) # gives the legend lines the correct color and width
    #}
    }
    dev.off()
  }
}

# Same figure, but in log scale
plot_log <- function(list_arg)
{
  for (data in list_arg)
  {
  par(mfrow=c(1,1))
  which_data = data$dataset[1]
  print(which_data)
  namefile = paste(path,"survVSlogmut_", which_data, ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  plot(log(data[which(data$group == "Responders"),]$nonsynonymous), data[which(data$group == "Responders"),]$overall_survival,
       xlab= "Log number of nonsynonymous mutations",
       ylab= "Overall Survival", pch = 16,col = "cyan3", cex = 1.2, main = which_data,      
       xlim = c(min(log(data$nonsynonymous)),max(log(data$nonsynonymous))), 
       ylim= c(min(na.omit(data$overall_survival)),max(na.omit(data$overall_survival))))
  if (which_data == "nslc")
  {
    abline(h=2*365, lty = 2, col = "purple")
    text(max(data$nonsynonymous)/2, 365*2+55, "OS = 2 year", cex = 1.2, col = "purple")
  } else   {
    abline(h=365, lty = 2, col = "purple")
    text(max(data$nonsynonymous)/2, 420, "OS = 1 year", cex = 1.2, col = "purple")
  }
  points(log(data[which(data$group == "Responders"),][which(data[which(data$group == "Responders"),]$dead == 1),]$nonsynonymous), 
       data[which(data$group == "Responders"),][which(data[which(data$group == "Responders"),]$dead == 1),]$overall_survival,
       pch = 1, cex = 1.2)
  points(log(data[which(data$group == "Nonresponders"),]$nonsynonymous), 
       data[which(data$group == "Nonresponders"),]$overall_survival, col ="red",
       pch = 16, cex = 1.2)
  points(log(data[which(data$group == "Nonresponders"),][which(data[which(data$group == "Nonresponders"),]$dead == 1),]$nonsynonymous), 
       data[which(data$group == "Nonresponders"),][which(data[which(data$group == "Nonresponders"),]$dead == 1),]$overall_survival, xlab= "Number of Nonsynonymous mutations",
       ylab= "Overall Survival", pch = 1, cex = 1.2)
#legend("topright", # places a legend at the appropriate place 
#       c("Responders","Non-responders", "Dead"), # puts text in the legend
#       pch = c(16,16,1),
#       col=c("cyan3", "red", "black")) # gives the legend lines the correct color and width
  if (which_data == "mel2"|| which_data == "Snyder2")
    {
    points(log(data[special_patients,]$nonsynonymous), data[special_patients,]$overall_survival, 
         col ="grey", cex= 1.2, pch = 16)
    points(log(data[special_patients,][which(data[special_patients,]$dead == 1),]$nonsynonymous), 
         data[special_patients,][which(data[special_patients,]$dead == 1),]$overall_survival
         , pch = 1, cex = 1.2)
  #  legend("topright", # places a legend at the appropriate place 
  #         c("Responders","Non-responders","Dead", "Reassigned"), # puts text in the legend
  #         pch = c(16,16,1,16), # gives the legend appropriate symbols (lines),
  #         col=c("cyan", "red","black", "grey")) # gives the legend lines the correct color and width
  }
  
  if (which_data == "mel3")
    {
    points(log(data[special_patients,]$nonsynonymous), data[special_patients,]$overall_survival, 
           col ="grey", cex= 1.2, pch = 16)
    points(log(data[special_patients,][which(data[special_patients,]$dead == 1),]$nonsynonymous), 
           data[special_patients,][which(data[special_patients,]$dead == 1),]$overall_survival
           , pch = 1, cex = 1.2)
    #  legend("topright", # places a legend at the appropriate place 
    #         c("Responders","Non-responders","Dead", "Subset"), # puts text in the legend
    #         pch = c(16,16,1,16), # gives the legend appropriate symbols (lines),
    #         col=c("cyan", "red","black", "grey")) # gives the legend lines the correct color and width
  }
  dev.off()
  }
}

##############################################################
# Association between Mutational Burden and Clinical Benefit #
##############################################################
# Boxplots of responders and non responders mutational load and 
# associated p-value by Mann-Whitney test

boxplot_MW <- function()
{
  par(mar=c(6, 4, 4, 6), xpd=TRUE)
  namefile = paste(path,"mutBox", ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)

  if (nlevels(as.factor(total$group))== 3)
    {
    col_ =  c("red", "cyan3", "grey") 
    total$group=factor(total$group , levels=levels(factor(total$group))[c(2,3,1)])
    boxplot(log(total$nonsynonymous) ~ interaction(as.factor(total$group), as.factor(total$dataset)),
          data = total, lwd = 2, 
          pwcol = 1 + as.numeric(total$dead),offset = .5,
          ylab = 'Log number of nonsynonymous mutations',
          xaxt = "n",
          na = c('Hugo', 'Snyder', 'Van Allen'),
          col = c(col_,col_, col_),
          ylim = c(min(log(total$nonsynonymous)), max(log(total$nonsynonymous)) + 2),
          outline = F,
          at = c(1,2,3,4,5,6,7,8,9),
          medlwd = 2)
    mtext(c("mel1", "mel2", "mel3"),1 ,line=1,at=c(1.5 , 4.5 , 8))
    pvalue_hugo = wilcox.test(total_hugo[which(total_hugo$group == "Responders"),]$nonsynonymous, 
                          total_hugo[which(total_hugo$group == "Nonresponders"),]$nonsynonymous)
    pvalue_str = paste("p-value = ", round(pvalue_hugo$p.value,4))
    arrows(1, max(log(total_hugo$nonsynonymous)) +1, 2, max(log(total_hugo$nonsynonymous)) +1, 
         lty = 1, cex = 3, col = "black", lwd = 1, angle  = 90, code = 3, length = 0.07)
    text(1.5, max(log(total_hugo$nonsynonymous)) + 1.5, pvalue_str, col = "black", cex = 1)

    pvalue_snyder = wilcox.test(total_snyder[which(total_snyder$group == "Responders"),]$nonsynonymous, 
                            total_snyder[which(total_snyder$group == "Nonresponders"),]$nonsynonymous)
    pvalue_str = paste("p-value = ", round(pvalue_snyder$p.value,4))
    arrows(4, max(log(total_snyder$nonsynonymous)) +1, 5, max(log(total_snyder$nonsynonymous)) +1, 
         lty = 1, cex = 3, col = "black", lwd = 1, angle  = 90, code = 3, length = 0.07)
    text(4.5, max(log(total_snyder$nonsynonymous)) + 1.5, pvalue_str, col = "black", cex = 1)

    pvalue_vanallen = wilcox.test(total_vanallen[which(total_vanallen$group == "Responders"),]$nonsynonymous, 
                              total_vanallen[which(total_vanallen$group == "Nonresponders"),]$nonsynonymous)
    pvalue_str = paste("p-value = ", round(pvalue_vanallen$p.value,4))
    arrows(7, max(log(total_vanallen$nonsynonymous)) +1, 8, max(log(total_vanallen$nonsynonymous)) +1, 
       lty = 1, cex = 3, col = "black", lwd = 1, angle  = 90, code = 3, length = 0.07)
    text(7.5, max(log(total_vanallen$nonsynonymous)) + 1.5, pvalue_str, col = "black", cex = 1)
    #legend("top", # places a legend at the appropriate place 
    #       inset=c(0, -0.15),
    #       c("Responders","Non-responders", "Subset", 'Dead'), # puts text in the legend
    #       pch = c(15,15,15,1,16), # gives the legend appropriate symbols (lines),
    #       col=c("cyan", "red","grey", "black",  "black"),
    #       horiz = T, bty = "n")

    stripchart(log(total_hugo[which(total_hugo$dead == 0),]$nonsynonymous) ~ 
             as.factor(total_hugo[which(total_hugo$dead == 0),]$group),
             vertical = TRUE, data = total_hugo[which(total_hugo$dead == 0),], 
             method = "jitter", add = TRUE, pch = 16, col = 'black', cex = 0.7)

    stripchart(log(total_hugo[which(total_hugo$dead == 1),]$nonsynonymous) ~ 
             as.factor(total_hugo[which(total_hugo$dead == 1),]$group), 
             vertical = TRUE, data = total_hugo[which(total_hugo$dead == 1),], 
             method = "jitter", add = TRUE, pch = 1, col = 'black', cex = 0.7)
  
    stripchart(log(total_snyder[which(total_snyder$dead == 0),]$nonsynonymous) ~ 
             as.factor(total_snyder[which(total_snyder$dead == 0),]$group),
             vertical = TRUE, data = total_snyder[which(total_snyder$dead == 0),], 
             method = "jitter", add = TRUE, pch = 16, col = 'black', at = c(4,5), cex = 0.7)

    stripchart(log(total_snyder[which(total_snyder$dead == 1),]$nonsynonymous) ~ 
             as.factor(total_snyder[which(total_snyder$dead == 1),]$group), 
             vertical = TRUE, data = total_snyder[which(total_snyder$dead == 1),], 
             method = "jitter", add = TRUE, pch = 1, col = 'black', at = c(4,5), cex = 0.7)
  
    stripchart(log(total_vanallen[which(total_vanallen$dead == 0),]$nonsynonymous) ~ 
             as.factor(total_vanallen[which(total_vanallen$dead == 0),]$group),
             vertical = TRUE, data = total_vanallen[which(total_vanallen$dead == 0),], 
             method = "jitter", add = TRUE, pch = 16, col = 'black', at = c(9,7,8), cex = 0.7)

    stripchart(log(total_vanallen[which(total_vanallen$dead == 1),]$nonsynonymous) ~ 
             as.factor(total_vanallen[which(total_vanallen$dead == 1),]$group), 
             vertical = TRUE, data = total_vanallen[which(total_vanallen$dead == 1),], 
             method = "jitter", add = TRUE, pch = 1, col = 'black', at = c(9,7,8), cex = 0.7)
}

  if (nlevels(as.factor((total)$group))== 2)
    {
    col_ =  c("red", "cyan3") 
    total$group=factor(total$group, levels=levels(factor(total$group)))
    boxplot(log(total$nonsynonymous) ~ interaction(as.factor(total$group), 
          as.factor(total$dataset)),
          data = total, lwd = 2, 
          pwcol = 1 + as.numeric(total$dead),offset = .5,
          ylab = 'Log number of nonsynonymous mutations',
          xaxt = "n",
          na = c('Hugo', 'Snyder', 'Van Allen'),
          col = c(col_,col_,col_),
          ylim = c(min(log(total$nonsynonymous)), max(log(total$nonsynonymous)) + 2),
          outline = F,
          at = c(1,2,4,5,7,8),
          medlwd = 2)
    mtext(c("mel1", "mel2", "mel3"),1 ,line=1,at=c(1.5, 4.5 , 7.5))
    pvalue_hugo = wilcox.test(total_hugo[which(total_hugo$group == "Responders"),]$nonsynonymous, 
                            total_hugo[which(total_hugo$group == "Nonresponders"),]$nonsynonymous)
    pvalue_str = paste("p-value = ", round(pvalue_hugo$p.value,4))
    arrows(1, max(log(total_hugo$nonsynonymous)) +1, 2, max(log(total_hugo$nonsynonymous)) +1, 
         lty = 1, cex = 3, col = "black", lwd = 1, angle  = 90, code = 3, length = 0.07)
    text(1.5, max(log(total_hugo$nonsynonymous)) + 1.5, pvalue_str, col = "black", cex = 1)
  
    pvalue_snyder = wilcox.test(total_snyder[which(total_snyder$group == "Responders"),]$nonsynonymous, 
                              total_snyder[which(total_snyder$group == "Nonresponders"),]$nonsynonymous)
    pvalue_str = paste("p-value = ", round(pvalue_snyder$p.value,4))
    arrows(4, max(log(total_snyder$nonsynonymous)) + 2, 5, max(log(total_snyder$nonsynonymous)) +2, 
         lty = 1, cex = 3, col = "black", lwd = 1, angle  = 90, code = 3, length = 0.07)
    text(4.5, max(log(total_snyder$nonsynonymous)) + 2.5, pvalue_str, col = "black", cex = 1)
  
    pvalue_vanallen = wilcox.test(total_vanallen[which(total_vanallen$group == "Responders"),]$nonsynonymous, 
                                total_vanallen[which(total_vanallen$group == "Nonresponders"),]$nonsynonymous)
    pvalue_str = paste("p-value = ", round(pvalue_vanallen$p.value,4))
    arrows(7, max(log(total_vanallen$nonsynonymous)) +1, 8, max(log(total_vanallen$nonsynonymous)) +1, 
         lty = 1, cex = 3, col = "black", lwd = 1, angle  = 90, code = 3, length = 0.07)
    text(7.5, max(log(total_vanallen$nonsynonymous)) + 1.5, pvalue_str, col = "black", cex = 1)
    #  legend("top", # places a legend at the appropriate place 
    #         inset=c(0, -0.1),
    #         c("Responders","Non-responders"), # puts text in the legend
    #         pch = c(15,15), # gives the legend appropriate symbols (lines),
    #         col=c("cyan", "red"),
    #         horiz = T, bty = "n")
    stripchart(log(total_hugo[which(total_hugo$dead == 0),]$nonsynonymous) ~ 
               as.factor(total_hugo[which(total_hugo$dead == 0),]$group),
               vertical = TRUE, data = total_hugo[which(total_hugo$dead == 0),], 
               method = "jitter", add = TRUE, pch = 16, col = 'black', cex = 0.7)
  
    stripchart(log(total_hugo[which(total_hugo$dead == 1),]$nonsynonymous) ~ 
             as.factor(total_hugo[which(total_hugo$dead == 1),]$group), 
             vertical = TRUE, data = total_hugo[which(total_hugo$dead == 1),], 
             method = "jitter", add = TRUE, pch = 1, col = 'black', cex = 0.7)
  
    stripchart(log(total_snyder[which(total_snyder$dead == 0),]$nonsynonymous) ~ 
              as.factor(total_snyder[which(total_snyder$dead == 0),]$group),
              vertical = TRUE, data = total_snyder[which(total_snyder$dead == 0),], 
              method = "jitter", add = TRUE, pch = 16, col = 'black', at = c(5), cex = 0.7)
  
    stripchart(log(total_snyder[which(total_snyder$dead == 1),]$nonsynonymous) ~ 
               as.factor(total_snyder[which(total_snyder$dead == 1),]$group), 
               vertical = TRUE, data = total_snyder[which(total_snyder$dead == 1),], 
               method = "jitter", add = TRUE, pch = 1, col = 'black', at = c(4,5), cex = 0.7)
  
    stripchart(log(total_vanallen[which(total_vanallen$dead == 0),]$nonsynonymous) ~ 
               as.factor(total_vanallen[which(total_vanallen$dead == 0),]$group),
               vertical = TRUE, data = total_vanallen[which(total_vanallen$dead == 0),], 
               method = "jitter", add = TRUE, pch = 16, col = 'black', at = c(8), cex = 0.7)
  
    stripchart(log(total_vanallen[which(total_vanallen$dead == 1),]$nonsynonymous) ~ 
               as.factor(total_vanallen[which(total_vanallen$dead == 1),]$group), 
              vertical = TRUE, data = total_vanallen[which(total_vanallen$dead == 1),], 
              method = "jitter", add = TRUE, pch = 1, col = 'black', at = c(7,8), cex = 0.7)
  }
dev.off()  
}


############################
# Permutation Mann-Whitney #
############################
# Randomly changes labels of response and compute Mann-Whitney p-value
# for each obtained dataset

permutation_MW <- function(list_arg)
{
  for (data in list_arg)
  {
  which_data = data$dataset[1]
  # Mann_Whitney (benefit vs no benefit)
  namefile = paste(path, "Perm_MW", which_data, ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  real_test <- wilcox.test(na.omit(data[which(data$group == "Responders"),]$nonsynonymous), 
                           na.omit(data[which(data$group == "Nonresponders"),]$nonsynonymous))
  real_W <- as.numeric(as.character(real_test$statistic))
  list_i = list()
  list_W= vector()
  j = 1
  try(
    for (i in 1:as.numeric(params_values[1]))
    {
      all_nonsynonymous <- append(data[which(data$group == "Responders"),]$nonsynonymous, 
                                  data[which(data$group == "Nonresponders"),]$nonsynonymous)
      all_nonsynonymous = sample(all_nonsynonymous)
      data_benefit_shuffled = all_nonsynonymous[1:length(data[which(data$group == "Responders"),]$nonsynonymous)]
      data_nobenefit_shuffled = all_nonsynonymous[length(data[which(data$group == "Responders"),]$nonsynonymous)+1:length(all_nonsynonymous)]
      res <- wilcox.test(data_benefit_shuffled, data_nobenefit_shuffled)
      W = as.numeric(as.character(res$statistic))
      list_i[j] = i
      list_W[j] = W
      j = j+1
    })
  h <- hist(list_W, breaks = as.numeric(params_values[1])/50, main = which_data, xlab = "W statisitics", 
            xlim = c(min(list_W),max(real_W + real_W/20, max(list_W))))
  abline(v= real_W, col = 'red')
  pvalue_comp = length(which(list_W > real_W))/as.numeric(params_values[1])
  pvalue = paste("p-value = ",pvalue_comp)
  legend("topright", # places a legend at the appropriate place 
         c(pvalue), # puts text in the legend
         pch = c(".")) # gives the legend lines the correct color and width
  dev.off()
  }
}

#####################
# Contingency table #
#####################

# To compare with random data when benefit is not redefined
cont_table <- function(list_arg)
{
  for (data in list_arg)
  {  
  list_i = list()
  list_p = list()
  h1 <- hist(log(data$nonsynonymous), plot = F)
  for (i in 1:1000)
  {
    x<-rnorm(length(data$overall_survival), 
             mean(log(data$nonsynonymous)), 
             sd(log(data$nonsynonymous)))
    h2 <- hist(x, plot = F)
    breaks_ = sort(unique(append(h1$breaks, h2$breaks)))
    h1 <- hist(log(data$nonsynonymous), breaks = breaks_, plot = F)
    h2 <- hist(x, breaks = breaks_, plot = F)
    cont_table = t(cbind(h1$counts, h2$counts))
    cont_table = cont_table[,which(!apply(cont_table,2,FUN = function(x){all(x == 0)}))]
    list_i[i]= i
    list_p[i]= as.numeric(chisq.test(cont_table, simulate.p.value = T, B = 20000)$p.value)
    }
  hist(sort(as.numeric(list_p)), breaks = 100)
  }
}


#########################################################
# Comparison with random data when benefit is redefined #
#########################################################

rand_comp <- function(list_arg)
{
  for (data in list_arg)
  {  
    list_p = list()
    for (i in 1:as.numeric(params_values[1]))
    {
    y<-data$overall_survival
    x<-rnorm(length(data$overall_survival), 
             mean(log(data$nonsynonymous)), 
             sd(log(data$nonsynonymous)))
    response = y
    response[which(y > 365)] =  "Responders"
    response[which(y < 365)] =  "Nonresponders"
    rand_data = as.data.frame(cbind(y,x,response),
                              as.is = F,
                              stringsAsFactors = FALSE)
    list_p[i] = as.numeric(wilcox.test(as.numeric(rand_data[which(rand_data$response == "Responders"),]$x), 
                as.numeric(rand_data[which(rand_data$response == "Nonresponders"),]$x))$p.value)
    }
    real_p = as.numeric(unlist(wilcox.test(as.numeric(data[which(data$group == "Responders"),]$nonsynonymous),
                as.numeric(data[which(data$group == "Nonresponders"),]$nonsynonymous))$p.value))
    p_value = sum(list_p < 0.05)/as.numeric(params_values[1])
    print(p_value)
    h <- hist(sort(as.numeric(list_p)), breaks = as.numeric(params_values[1])/50, main = 'Mann-Whitney p-values for ')
    abline(v = real_p, col = 'red')
  }
}

####################
# Plot random data #
####################

# Plot overall survival against nonsynonymous mutation load for 
# random data with same y distribution and normal x distribution
# (with same mean and standard deiation as the data)
rand_plot <- function(list_arg)
{
scatterhist = function(x, y, xlab="", ylab=""){
    zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
    xhist = hist(x, plot=FALSE)
    yhist = hist(y, plot=FALSE)
    top = max(c(xhist$counts, yhist$counts))
    par(mar=c(3,3,1,1))
    plot(x,y,xlab = "Log number of nonsynonymous mutations", ylab = "Overall survival (in days)")
    par(mar=c(0,3,1,1))
    barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
    par(mar=c(3,0,1,1))
    barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
    par(oma=c(3,3,0,0))
    mtext(xlab, side=1, line=1, outer=TRUE, adj=0,     at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
    mtext(ylab, side=2, line=1, outer=TRUE, adj=0,     at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
  }
  
  namefile = paste(path, "Random_data", ".tiff")
  tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  y<-data$overall_survival
  x<-rnorm(length(data$overall_survival), 
           mean(log(data$nonsynonymous)), 
           sd(log(data$nonsynonymous)))
  scatterhist(x,y)
  dev.off()
}
  
#####################
# Correlation plots #
#####################
# Spearman and Pearsons correlation coefficient for random data as defined above
# to compare with real data

corr_plots <- function(list_arg)
{ 
  for (data in list_arg)
  {
    which_data = data$dataset[1]
    namefile = paste(path, paste("Correlations_", which_data), ".tiff")
    tiff(namefile, width = 4, height = 8, units = 'in', res = 500)
    par(mfrow=c(1,1))
    spearman_list = list()
    pearson_list = list()
    y<-data$overall_survival
    for (i in 1:as.numeric(params_values[1]))
    {
      x<-rnorm(length(data$overall_survival), 
             mean(log(data$nonsynonymous)), 
             sd(log(data$nonsynonymous)))
      pearson_list[i] = as.numeric(cor.test(x,y, method = 'pearson')$estimate)
      spearman_list[i] = as.numeric(cor.test(x,y, method = 'spearman')$estimate)
    }
    Pearson_data = cbind(as.numeric(pearson_list), 
          replicate(length(as.numeric(pearson_list)), 'Pearson'))
    Spearman_data = cbind(as.numeric(spearman_list), 
                          replicate(length(as.numeric(spearman_list)), 'Spearman'))
    cor_data = as.data.frame(rbind(Spearman_data, Pearson_data))
    boxplot(as.numeric(as.character(cor_data$V1))~(cor_data$V2), col = c('orange', 'pink'),
            names = c('Pearson', 'Spearman'), main = which_data, ylim = c(-1,1),
            ylab = "Correlation coefficient")
    real_points = as.data.frame(cbind(c(as.numeric(cor.test(data$nonsynonymous,y, method = 'pearson')$estimate),
                                        as.numeric(cor.test(data$nonsynonymous,y, method = 'spearman')$estimate)),
                                      c('Pearson', 'Spearman')))
    stripchart(as.numeric(as.character(real_points$V1)) ~ as.factor(real_points$V2),
               vertical = TRUE, method = "stack", add = TRUE, pch = '*', col = 'red', cex = 3)
    dev.off()
  }
}

#####################
# Trimming analysis #
#####################
# Plot Mann-Whitney pvalue when data is trimmed. Also
# plot evolution of the proportion of non-responders
# as the data is trimmed

trim_analysis_head <- function(list_arg)
{ 
  for (data in list_arg)
  {   
  par(mfrow=c(1,1))
  #namefile = paste(path,"Trimmed_data_", which_data, ".tiff")
  #tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
  list_wilcox_trimmed = list()
  data_trimmed <- data[order(data$nonsynonymous),]
  #25
  stack <- data.frame(
    ticker=character(),
    value=numeric(),
    stringsAsFactors=FALSE)  
  for (i in 1:floor(length(data[which(data$group == "Responders"),]$nonsynonymous)))
  {
    #data_trimmed <- head(data_trimmed,-1) #trim tail
    data_trimmed <- tail(data_trimmed,-1) #trim head
    list_wilcox_trimmed[i] <- wilcox.test(data_trimmed[which(data_trimmed$group == "Responders"),]$nonsynonymous, 
                                          data_trimmed[which(data_trimmed$group == "Nonresponders"),]$nonsynonymous)$p.value
    stack = rbind(stack, list(length(data_trimmed[which(data_trimmed$group == "Responders"),]$nonsynonymous),
      length(data_trimmed[which(data_trimmed$group == "Nonresponders"),]$nonsynonymous)))
  }
  prop_tabl = prop.table(t(stack), margin = 2)
  
  par(mar=c(5,4,4,5)+.1)
  plot(as.numeric(list_wilcox_trimmed),
       xlab = "Number of values trimmed on each side",
       ylab = "P-value (Mann-Whitney test)", type = 'l', col = 'black') 
  abline(h=0.05)
  points(as.numeric(list_wilcox_trimmed), pch = 16, cex = 0.7)
  par(new=TRUE)
  plot(prop_tabl[2,],xaxt="n",yaxt="n",xlab="",ylab="", type = 'l', col = 'red')  
  points(prop_tabl[2,], pch = 16, cex = 0.7, col = 'red')  
  mtext("Proportion of non-responders",side=4,line=3, col = 'red')
  axis(4)
  }
}

trim_analysis_tails <- function(list_arg)
{ 
  for (data in list_arg)
  {   
    par(mfrow=c(1,1))
    #namefile = paste(path,"Trimmed_data_", which_data, ".tiff")
    #tiff(namefile, width = 8, height = 6, units = 'in', res = 500)
    list_wilcox_trimmed = list()
    data_trimmed <- data[order(data$nonsynonymous),]
    #25
    stack <- data.frame(
      ticker=character(),
      value=numeric(),
      stringsAsFactors=FALSE)  
    for (i in 1:floor(length(data[which(data$group == "Responders"),]$nonsynonymous)))
    {
      data_trimmed <- head(data_trimmed,-1) #trim tail
      #data_trimmed <- tail(data_trimmed,-1) #trim head
      list_wilcox_trimmed[i] <- wilcox.test(data_trimmed[which(data_trimmed$group == "Responders"),]$nonsynonymous, 
                                            data_trimmed[which(data_trimmed$group == "Nonresponders"),]$nonsynonymous)$p.value
      stack = rbind(stack, list(length(data_trimmed[which(data_trimmed$group == "Responders"),]$nonsynonymous),
                                length(data_trimmed[which(data_trimmed$group == "Nonresponders"),]$nonsynonymous)))
    }
    prop_tabl = prop.table(t(stack), margin = 2)
    par(new=FALSE)
    
    par(mar=c(5,4,4,5)+.1)
    plot(as.numeric(list_wilcox_trimmed),
         xlab = "Number of values trimmed on each side",
         ylab = "P-value (Mann-Whitney test)", type = 'l', col = 'black')
    abline(h = 0.05)
    points(as.numeric(list_wilcox_trimmed), pch = 16, cex = 0.7)
    par(new=TRUE)
    plot(prop_tabl[1,],xaxt="n",yaxt="n",xlab="",ylab="", type = 'l', col = 'cyan3',col.axis = 'cyan3')
    points(prop_tabl[1,], pch = 16, cex = 0.7, col = 'cyan3')
    mtext("Proportion of responders",side=4,line=3, col = 'cyan3')
    axis(4, col = 'cyan3', col.axis = 'cyan3', col.ticks = 'cyan3')
  }
}


#############################
# Threshold Mutation and LR #
#############################
# For each mutation cutoff, compute associated p-value of log-rang test

threshol_mut <- function(list_arg)
{ 
  for (data in list_arg)
  {   
  ## P_value vs threshold of mutations
  par(mfrow=c(1,1))
  list_i = list()
  list_pvalues = list()
  j = 1
  try(
    for (i in sort(unique(data$nonsynonymous)))
    {
      print(i)
      res <- survdiff(Surv(data$overall_survival, data$dead) ~ data$nonsynonymous > i)
      pvalue = pchisq(res$chisq, length(res$n)-1, lower.tail = FALSE)
      list_i[j] = i
      list_pvalues[j] = pvalue
      j = j+1
    })
  #namefile3 = paste("Pvalue_VS_MutationThreshold_", which_data, ".jpg")
  #jpeg(namefile3)
  plot(list_i, list_pvalues, pch = ".", col = 'purple', cex = 5, xlab = "Discriminant number of mutations",
       ylab = "p-values", main = which_data)
  abline(a=0.05, b=0, col = "purple")
  rm(list_i, list_pvalues, j, res)
  dev.off()
}
}

########################################
# Permutation test for median analysis #
########################################

threshold_median <- function(list_arg)
{ 
  for (data in list_arg)
  {   
    ## P_value vs threshold of mutations
    par(mfrow=c(1,1))
    list_i = list()
    list_deltaS = list()
    j = 1
    try(
      for (i in sort(unique(data$nonsynonymous)))
      {
        print(i)
        S1 <- median(data$overall_survival[data$nonsynonymous > i])
        S2 <- median(data$overall_survival[data$nonsynonymous <= i])
        S1[is.na(S1)] <- 0
        deltaS = S1-S2
        list_i[j] = i
        list_deltaS[j] = deltaS
        j = j+1
      })
    #namefile3 = paste("Pvalue_VS_MutationThreshold_", which_data, ".jpg")
    #jpeg(namefile3)
    plot(list_i, list_deltaS, pch = ".", col = 'purple', cex = 5, xlab = "Discriminant number of mutations",
         ylab = "Delta S", main = which_data)
    abline(a=0.05, b=0, col = "purple")
    rm(list_i, list_deltaS, j, res)
    dev.off()
  }
}


########################
# Permutation log-rank #
########################

perm_LR <- function(list_arg)
{ 
  for (data in list_arg)
  { 
  pvalue_comp = vector()
  mut_chosen_list = vector()
  k = 0
  try(
    for (mut_chosen in (min(data$nonsynonymous)):(max(data$nonsynonymous)))
    {
      real_test <- survdiff(Surv(data$overall_survival, data$dead) ~ data$nonsynonymous > mut_chosen)
      real_Z <- as.numeric(as.character(((real_test$obs[1]-real_test$exp[1])^2)/real_test$exp[1] + ((real_test$obs[2]-real_test$exp[2])^2)/real_test$exp[2]))
      list_i = list()
      list_Z= vector()
      j = 1
      for (i in 1:as.numeric(params_values[1]))
      {
        overall_survival_shuffled = sample(data$overall_survival)
        res <- survdiff(Surv(overall_survival_shuffled, data$dead) ~ data$nonsynonymous > mut_chosen)
        Z = as.numeric(as.character(((res$obs[1]-res$exp[1]))/res$exp[1] + ((res$obs[2]-res$exp[2]))/res$exp[2]))
        list_i[j] = i
        list_Z[j] = Z
        j = j+1
      }
      pvalue_comp[k] = length(which(list_Z > real_Z))/as.numeric(params_values[1])
      mut_chosen_list[k] = mut_chosen
      k = k+1
    })
  #namefile6 = paste("Permutation_Log-Rank_Multiple",which_data, ".jpg")
  #jpeg(namefile6)
  plot(mut_chosen_list, pvalue_comp, xlim = c(min(total$nonsynonymous),max(total$nonsynonymous)), 
       ylim = c(min(pvalue_comp), max(pvalue_comp)), pch = ".",
       xlab = "Discriminant number of mutations", ylab = "p_value")
  lines(mut_chosen_list, pvalue_comp)
  #dev.off()
  rm(h,i, list_i ,list_Z, real_test, real_Z, pvalue_comp, mut_chosen_list, mut_chosen, overall_survival_shuffled, res, Z)
  }
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


####################
# Survival and age #
####################

if (params_tests[9]  == T)
{
  par(mfrow=c(1,1))
  namefile = paste(path, "survVSlogmutage_", which_data, ".tiff")
  tiff(namefile, width = 12, height = 8, units = 'in', res = 500)
  plot(log(total[which(total$group == "Responders"),]$nonsynonymous), 
       total[which(total$group == "Responders"),]$overall_survival, xlab= "Number of Nonsynonymous mutations",
       ylab= "Overall Survival", col ="cyan", pch = 16, 
       cex = total[which(total$group == "Responders"),]$age*0.02, main = which_data,      
       xlim = c(min(log(total$nonsynonymous)),max(log(total$nonsynonymous))), 
       ylim= c(min(na.omit(total$overall_survival)),max(na.omit(total$overall_survival))))
  points(log(total[which(total$group == "Nonresponders"),]$nonsynonymous), 
         total[which(total$group == "Nonresponders"),]$overall_survival, 
         cex = total[which(total$group == "Nonresponders"),]$age*0.02, pch = 16, col = "red")
  legend("topright", # places a legend at the appropriate place 
         c("Responders","Non-responders"), # puts text in the legend
         pch = c(".","."),
         lwd = c(1,1),
         col=c("blue", "red")) # gives the legend lines the correct color and width
  dev.off()
}

###########################
# Stage severity analysis #
###########################

if (params_tests[10]  == T && which_data != "Rizvi")
{
  par(mfrow=c(1,1))
  namefile = paste(path, "survVSlogmutstage_", which_data, ".tiff")
  tiff(namefile, width = 12, height = 8, units = 'in', res = 500)
  plot(log(total[which(total$stage == params_data[4]),]$nonsynonymous), 
       total[which(total$stage == params_data[4]),]$overall_survival, 
       pch = 16, col = 'yellow', cex = 1,
       xlim = c(min(log(total$nonsynonymous)), max(log(total$nonsynonymous))), 
       ylim = c(min(na.omit(total$overall_survival)), max(na.omit(total$overall_survival))),
       xlab = "Log number of non synonymous mutations", ylab = "Overall survival", main = which_data)
  points(log(total[which(total$stage == params_data[5]),]$nonsynonymous), 
         total[which(total$stage == params_data[5]),]$overall_survival,
         pch = 16, col = 'orange', cex = 1)
  points(log(total[which(total$stage == params_data[6]),]$nonsynonymous), 
         total[which(total$stage == params_data[6]),]$overall_survival,
         pch = 16, col = 'red', cex = 1)
  points(log(total[which(total$stage == params_data[7]),]$nonsynonymous), 
         total[which(total$stage == params_data[7]),]$overall_survival,
         pch = 16, col = 'black', cex = 1)  
  legend("topright", # places a legend at the appropriate place 
         c(params_data[4],params_data[5],params_data[6],params_data[7]), # puts text in the legend
         pch = c(".",".",".","."), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5,2.5,2.5),col=c("yellow", "orange","red", "black")) # gives the legend lines the correct color and width
  dev.off()
}


if (params_tests[11]  == T)
{
  # Same in log scale
  namefile = paste(path, "survVSlogmutgender_", which_data, ".tiff")
  tiff(namefile, width = 12, height = 8, units = 'in', res = 500)
  plot(log(total[which(total$gender == params_data[8]),]$nonsynonymous), 
       (total[which(total$gender == params_data[8]),]$overall_survival), pch = 16, col = 'blue', cex = 1,
       xlim = c(min(log(total$nonsynonymous)), max(log(total$nonsynonymous))), 
       ylim = c(min(na.omit((total$overall_survival))), max(na.omit((total$overall_survival)))),
       xlab = "Log number of non synonymous mutations", ylab = "Log overall survival", main = which_data)
  points(log(total[which(total$gender == params_data[9]),]$nonsynonymous), 
         (total[which(total$gender == params_data[9]),]$overall_survival), pch = 16, col = 'pink', cex = 1)
  legend("topright", # places a legend at the appropriate place 
         c("Male", "Female"), # puts text in the legend
         pch = c(".","."), # gives the legend appropriate symbols (lines)
         lwd=c(2.5,2.5),col=c("blue", "pink")) # gives the legend lines the correct color and width
  dev.off()
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






library(GeneOverlap)

## Set Working Directory
setwd("/Users/anauntaroiu/Documents/Data")

## Load Essential Reaction Lists for Each Interval
Essential_Rxns_P_10 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_10.csv"))
Essential_Rxns_P_20 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_20.csv"))
Essential_Rxns_P_30 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_30.csv"))
Essential_Rxns_P_40 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_40.csv"))
Essential_Rxns_P_50 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_50.csv"))
Essential_Rxns_P_60 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_60.csv"))
Essential_Rxns_P_70 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_70.csv"))
Essential_Rxns_P_80 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_80.csv"))
Essential_Rxns_P_90 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_90.csv"))
Essential_Rxns_P_100 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_Models.csv"))
Essential_Rxns_Unconstrained <- as.matrix(read.csv("Essential_Rxns_Unconstrained_Model.csv"))

#################################################################
## Remove reactions essential in unconstrained model
#
# INPUT: Essential Reactions from Conditional Models, Essential Reactions from Unconstrained model
#
# OUTPUT: List of conditionally essential reactions

remove_unconstrained_ess_rxn <- function(Essential_Rxns, Essential_Rxns_Unconstrained){
  # Create list to hold conditionally essential reactions
  Conditionally_ess_rxns <- list()
  
  for (i in 1:dim(Essential_Rxns)[2]){ # Loop through each sample of essential rxns
    
    index <- grep('R_', Essential_Rxns[,i]) # Record postions of essential reactions
    sample <- as.matrix(Essential_Rxns[(index),i]) # Keep only ess rxns (remove NAs)
    Conditionally_ess_rxns[[i]] <- setdiff(sample, Essential_Rxns_Unconstrained) # Remove unconstrained ess rxns
  }
  return(Conditionally_ess_rxns)
}
#################################################
# Call remove_unconstrained_ess_rxn to get conditionally essential reactions

Essential_Rxns_P_10 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_10, Essential_Rxns_Unconstrained)
Essential_Rxns_P_20 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_20, Essential_Rxns_Unconstrained)
Essential_Rxns_P_30 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_30, Essential_Rxns_Unconstrained)
Essential_Rxns_P_40 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_40, Essential_Rxns_Unconstrained)
Essential_Rxns_P_50 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_50, Essential_Rxns_Unconstrained)
Essential_Rxns_P_60 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_60, Essential_Rxns_Unconstrained)
Essential_Rxns_P_70 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_70, Essential_Rxns_Unconstrained)
Essential_Rxns_P_80 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_80, Essential_Rxns_Unconstrained)
Essential_Rxns_P_90 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_90, Essential_Rxns_Unconstrained)
Essential_Rxns_P_100 <- remove_unconstrained_ess_rxn(Essential_Rxns_P_100, Essential_Rxns_Unconstrained)

#################################################################
## Get median number of conditionally essential reactions for each subset
#
# INPUT: Essential Reactions from one subset
#
# OUTPUT: Median number of Essential reactions

median_ess_rxn <- function(Essential_Rxns){
  Sample <- c(median(c(length(Essential_Rxns[[1]]), length(Essential_Rxns[[2]]), length(Essential_Rxns[[3]]), length(Essential_Rxns[[4]]), length(Essential_Rxns[[5]]),
                     length(Essential_Rxns[[6]]), length(Essential_Rxns[[7]]), length(Essential_Rxns[[8]]), length(Essential_Rxns[[9]]), length(Essential_Rxns[[10]]),                                                                                           
                     length(Essential_Rxns[[11]]), length(Essential_Rxns[[12]]), length(Essential_Rxns[[13]]), length(Essential_Rxns[[14]]), length(Essential_Rxns[[15]]))),
              median(c(length(Essential_Rxns[[16]]), length(Essential_Rxns[[17]]), length(Essential_Rxns[[18]]), length(Essential_Rxns[[19]]), length(Essential_Rxns[[20]]),
                     length(Essential_Rxns[[21]]), length(Essential_Rxns[[22]]), length(Essential_Rxns[[23]]), length(Essential_Rxns[[24]]), length(Essential_Rxns[[25]]),                                                   
                     length(Essential_Rxns[[26]]), length(Essential_Rxns[[27]]), length(Essential_Rxns[[28]]), length(Essential_Rxns[[29]]), length(Essential_Rxns[[30]]))),
              median(c(length(Essential_Rxns[[31]]), length(Essential_Rxns[[32]]), length(Essential_Rxns[[33]]), length(Essential_Rxns[[34]]), length(Essential_Rxns[[35]]),
                     length(Essential_Rxns[[36]]), length(Essential_Rxns[[37]]), length(Essential_Rxns[[38]]), length(Essential_Rxns[[39]]), length(Essential_Rxns[[40]]),                                                                                           
                     length(Essential_Rxns[[41]]), length(Essential_Rxns[[42]]), length(Essential_Rxns[[43]]), length(Essential_Rxns[[44]]), length(Essential_Rxns[[45]]))),
              median(c(length(Essential_Rxns[[46]]), length(Essential_Rxns[[47]]), length(Essential_Rxns[[48]]), length(Essential_Rxns[[49]]), length(Essential_Rxns[[50]]),
                     length(Essential_Rxns[[51]]), length(Essential_Rxns[[52]]), length(Essential_Rxns[[53]]), length(Essential_Rxns[[54]]), length(Essential_Rxns[[55]]),
                     length(Essential_Rxns[[56]]), length(Essential_Rxns[[57]]), length(Essential_Rxns[[58]]), length(Essential_Rxns[[59]]), length(Essential_Rxns[[60]]))),
              median(c(length(Essential_Rxns[[61]]), length(Essential_Rxns[[62]]), length(Essential_Rxns[[63]]), length(Essential_Rxns[[64]]), length(Essential_Rxns[[65]]),                                                                                           
                     length(Essential_Rxns[[66]]), length(Essential_Rxns[[67]]), length(Essential_Rxns[[68]]), length(Essential_Rxns[[69]]), length(Essential_Rxns[[70]]),
                     length(Essential_Rxns[[71]]), length(Essential_Rxns[[72]]), length(Essential_Rxns[[73]]), length(Essential_Rxns[[74]]), length(Essential_Rxns[[75]]))))
  return(Sample)
}
#################################################
# Call median_ess-rxn function
Sample_10 <- median_ess_rxn(Essential_Rxns_P_10)
Sample_20 <- median_ess_rxn(Essential_Rxns_P_20)
Sample_30 <- median_ess_rxn(Essential_Rxns_P_30)
Sample_40 <- median_ess_rxn(Essential_Rxns_P_40)
Sample_50 <- median_ess_rxn(Essential_Rxns_P_50)
Sample_60 <- median_ess_rxn(Essential_Rxns_P_60)
Sample_70 <- median_ess_rxn(Essential_Rxns_P_70)
Sample_80 <- median_ess_rxn(Essential_Rxns_P_80)
Sample_90 <- median_ess_rxn(Essential_Rxns_P_90)
Sample_100 <- c(median(length(Essential_Rxns_P_100[[1]]), length(Essential_Rxns_P_100[[2]]), length(Essential_Rxns_P_100[[3]])),
                median(length(Essential_Rxns_P_100[[4]]), length(Essential_Rxns_P_100[[5]]),length(Essential_Rxns_P_100[[6]])),
                median(length(Essential_Rxns_P_100[[7]]), length(Essential_Rxns_P_100[[8]]), length(Essential_Rxns_P_100[[9]])),
                median(length(Essential_Rxns_P_100[[10]]), length(Essential_Rxns_P_100[[11]]), length(Essential_Rxns_P_100[[12]])),
                median(length(Essential_Rxns_P_100[[13]]), length(Essential_Rxns_P_100[[14]]), length(Essential_Rxns_P_100[[15]])))

# Organize data into vectors
Sample_I <- c(Sample_10[1], Sample_20[1], Sample_30[1],
              Sample_40[1], Sample_50[1], Sample_60[1],
              Sample_70[1], Sample_80[1], Sample_90[1],
              Sample_100[1])
Sample_II <- c(Sample_10[2], Sample_20[2], Sample_30[2],
               Sample_40[2], Sample_50[2], Sample_60[2],
               Sample_70[2], Sample_80[2], Sample_90[2],
               Sample_100[2])
Sample_III <- c(Sample_10[3], Sample_20[3], Sample_30[3],
                Sample_40[3], Sample_50[3], Sample_60[3],
                Sample_70[3], Sample_80[3], Sample_90[3],
                Sample_100[3])
Sample_IV <- c(Sample_10[4], Sample_20[4], Sample_30[4],
               Sample_40[4], Sample_50[4], Sample_60[4],
               Sample_70[4], Sample_80[4], Sample_90[4],
               Sample_100[4])
Sample_V <- c(Sample_10[5], Sample_20[5], Sample_30[5],
              Sample_40[5], Sample_50[5], Sample_60[5],
              Sample_70[5], Sample_80[5], Sample_90[5],
              Sample_100[5])

Sample_I <- cbind(Sample_I, c(10,20,30,40,50,60,70,80,90,100), rep('I', times = 10))
Sample_II <- cbind(Sample_II, c(10,20,30,40,50,60,70,80,90,100), rep('II', times = 10))
Sample_III <- cbind(Sample_III, c(10,20,30,40,50,60,70,80,90,100), rep('III', times = 10))
Sample_IV <- cbind(Sample_IV, c(10,20,30,40,50,60,70,80,90,100), rep('IV', times = 10))
Sample_V <- cbind(Sample_V, c(10,20,30,40,50,60,70,80,90,100), rep('V', times = 10))
data <- rbind(Sample_I, Sample_II, Sample_III, Sample_IV, Sample_V)
colnames(data) <- c("Conditional_Essential_Reactions", "Percentage_of_Data", "Sample")

data <- data.frame(Essential_Rxns = as.numeric(data[,1]), Percent = as.numeric(data[,2]), Sample = data[,3])
library(ggplot2)
ggplot(data, aes(x=Percent, y=Essential_Rxns, colour = Sample)) + geom_line(lwd=1.1) + geom_point(size=2.4) + theme_bw() +
  xlab("Percentage of Transcriptomics Integrated") + ylab("Conditionally Essential Reactions") +
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
  theme(axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y=element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 12))

## Transcriptomics Data
## Load Essential Reaction Lists for Each Interval
Essential_Rxns_T_10 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_10.csv"))
Essential_Rxns_T_20 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_20.csv"))
Essential_Rxns_T_30 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_30.csv"))
Essential_Rxns_T_40 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_40.csv"))
Essential_Rxns_T_50 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_50.csv"))
Essential_Rxns_T_60 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_60.csv"))
Essential_Rxns_T_70 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_70.csv"))
Essential_Rxns_T_80 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_80.csv"))
Essential_Rxns_T_90 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_90.csv"))
Essential_Rxns_T_100 <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_Models.csv"))
Essential_Rxns_Unconstrained <- as.matrix(read.csv("Essential_Rxns_Unconstrained_Model.csv"))

#################################################################
## Remove reactions essential in unconstrained model
#
# INPUT: Essential Reactions from Conditional Models, Essential Reactions from Unconstrained model
#
# OUTPUT: List of conditionally essential reactions

remove_unconstrained_ess_rxn <- function(Essential_Rxns, Essential_Rxns_Unconstrained){
  # Create list to hold conditionally essential reactions
  Conditionally_ess_rxns <- list()
  
  for (i in 1:dim(Essential_Rxns)[2]){ # Loop through each sample of essential rxns
    
    index <- grep('R_', Essential_Rxns[,i]) # Record postions of essential reactions
    sample <- as.matrix(Essential_Rxns[(index),i]) # Keep only ess rxns (remove NAs)
    Conditionally_ess_rxns[[i]] <- setdiff(sample, Essential_Rxns_Unconstrained) # Remove unconstrained ess rxns
  }
  return(Conditionally_ess_rxns)
}
#################################################
# Call remove_unconstrained_ess_rxn to get conditionally essential reactions

Essential_Rxns_T_10 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_10, Essential_Rxns_Unconstrained)
Essential_Rxns_T_20 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_20, Essential_Rxns_Unconstrained)
Essential_Rxns_T_30 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_30, Essential_Rxns_Unconstrained)
Essential_Rxns_T_40 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_40, Essential_Rxns_Unconstrained)
Essential_Rxns_T_50 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_50, Essential_Rxns_Unconstrained)
Essential_Rxns_T_60 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_60, Essential_Rxns_Unconstrained)
Essential_Rxns_T_70 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_70, Essential_Rxns_Unconstrained)
Essential_Rxns_T_80 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_80, Essential_Rxns_Unconstrained)
Essential_Rxns_T_90 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_90, Essential_Rxns_Unconstrained)
Essential_Rxns_T_100 <- remove_unconstrained_ess_rxn(Essential_Rxns_T_100, Essential_Rxns_Unconstrained)

#################################################################
## Get median number of conditionally essential reactions for each subset
#
# INPUT: Essential Reactions from one subset
#
# OUTPUT: Median number of Essential reactions

median_ess_rxn <- function(Essential_Rxns){
  Sample <- c(median(c(length(Essential_Rxns[[1]]), length(Essential_Rxns[[2]]), length(Essential_Rxns[[3]]), length(Essential_Rxns[[4]]), length(Essential_Rxns[[5]]),
                       length(Essential_Rxns[[6]]), length(Essential_Rxns[[7]]), length(Essential_Rxns[[8]]), length(Essential_Rxns[[9]]), length(Essential_Rxns[[10]]),                                                                                           
                       length(Essential_Rxns[[11]]), length(Essential_Rxns[[12]]), length(Essential_Rxns[[13]]), length(Essential_Rxns[[14]]), length(Essential_Rxns[[15]]))),
              median(c(length(Essential_Rxns[[16]]), length(Essential_Rxns[[17]]), length(Essential_Rxns[[18]]), length(Essential_Rxns[[19]]), length(Essential_Rxns[[20]]),
                       length(Essential_Rxns[[21]]), length(Essential_Rxns[[22]]), length(Essential_Rxns[[23]]), length(Essential_Rxns[[24]]), length(Essential_Rxns[[25]]),                                                   
                       length(Essential_Rxns[[26]]), length(Essential_Rxns[[27]]), length(Essential_Rxns[[28]]), length(Essential_Rxns[[29]]), length(Essential_Rxns[[30]]))),
              median(c(length(Essential_Rxns[[31]]), length(Essential_Rxns[[32]]), length(Essential_Rxns[[33]]), length(Essential_Rxns[[34]]), length(Essential_Rxns[[35]]),
                       length(Essential_Rxns[[36]]), length(Essential_Rxns[[37]]), length(Essential_Rxns[[38]]), length(Essential_Rxns[[39]]), length(Essential_Rxns[[40]]),                                                                                           
                       length(Essential_Rxns[[41]]), length(Essential_Rxns[[42]]), length(Essential_Rxns[[43]]), length(Essential_Rxns[[44]]), length(Essential_Rxns[[45]]))),
              median(c(length(Essential_Rxns[[46]]), length(Essential_Rxns[[47]]), length(Essential_Rxns[[48]]), length(Essential_Rxns[[49]]), length(Essential_Rxns[[50]]),
                       length(Essential_Rxns[[51]]), length(Essential_Rxns[[52]]), length(Essential_Rxns[[53]]), length(Essential_Rxns[[54]]), length(Essential_Rxns[[55]]),
                       length(Essential_Rxns[[56]]), length(Essential_Rxns[[57]]), length(Essential_Rxns[[58]]), length(Essential_Rxns[[59]]), length(Essential_Rxns[[60]]))),
              median(c(length(Essential_Rxns[[61]]), length(Essential_Rxns[[62]]), length(Essential_Rxns[[63]]), length(Essential_Rxns[[64]]), length(Essential_Rxns[[65]]),                                                                                           
                       length(Essential_Rxns[[66]]), length(Essential_Rxns[[67]]), length(Essential_Rxns[[68]]), length(Essential_Rxns[[69]]), length(Essential_Rxns[[70]]),
                       length(Essential_Rxns[[71]]), length(Essential_Rxns[[72]]), length(Essential_Rxns[[73]]), length(Essential_Rxns[[74]]), length(Essential_Rxns[[75]]))))
  return(Sample)
}
#################################################
# Call median_ess-rxn function
Sample_10 <- median_ess_rxn(Essential_Rxns_T_10)
Sample_20 <- median_ess_rxn(Essential_Rxns_T_20)
Sample_30 <- median_ess_rxn(Essential_Rxns_T_30)
Sample_40 <- median_ess_rxn(Essential_Rxns_T_40)
Sample_50 <- median_ess_rxn(Essential_Rxns_T_50)
Sample_60 <- median_ess_rxn(Essential_Rxns_T_60)
Sample_70 <- median_ess_rxn(Essential_Rxns_T_70)
Sample_80 <- median_ess_rxn(Essential_Rxns_T_80)
Sample_90 <- median_ess_rxn(Essential_Rxns_T_90)
Sample_100 <- c(median(length(Essential_Rxns_T_100[[1]]), length(Essential_Rxns_T_100[[2]]), length(Essential_Rxns_T_100[[3]])),
                median(length(Essential_Rxns_T_100[[4]]), length(Essential_Rxns_T_100[[5]]),length(Essential_Rxns_T_100[[6]])),
                median(length(Essential_Rxns_T_100[[7]]), length(Essential_Rxns_T_100[[8]]), length(Essential_Rxns_T_100[[9]])),
                median(length(Essential_Rxns_T_100[[10]]), length(Essential_Rxns_T_100[[11]]), length(Essential_Rxns_T_100[[12]])),
                median(length(Essential_Rxns_T_100[[13]]), length(Essential_Rxns_T_100[[14]]), length(Essential_Rxns_T_100[[15]])))

# Organize data into vectors
Sample_I <- c(Sample_10[1], Sample_20[1], Sample_30[1],
              Sample_40[1], Sample_50[1], Sample_60[1],
              Sample_70[1], Sample_80[1], Sample_90[1],
              Sample_100[1])
Sample_II <- c(Sample_10[2], Sample_20[2], Sample_30[2],
               Sample_40[2], Sample_50[2], Sample_60[2],
               Sample_70[2], Sample_80[2], Sample_90[2],
               Sample_100[2])
Sample_III <- c(Sample_10[3], Sample_20[3], Sample_30[3],
                Sample_40[3], Sample_50[3], Sample_60[3],
                Sample_70[3], Sample_80[3], Sample_90[3],
                Sample_100[3])
Sample_IV <- c(Sample_10[4], Sample_20[4], Sample_30[4],
               Sample_40[4], Sample_50[4], Sample_60[4],
               Sample_70[4], Sample_80[4], Sample_90[4],
               Sample_100[4])
Sample_V <- c(Sample_10[5], Sample_20[5], Sample_30[5],
              Sample_40[5], Sample_50[5], Sample_60[5],
              Sample_70[5], Sample_80[5], Sample_90[5],
              Sample_100[5])

Sample_I <- cbind(Sample_I, c(10,20,30,40,50,60,70,80,90,100), rep('I', times = 10))
Sample_II <- cbind(Sample_II, c(10,20,30,40,50,60,70,80,90,100), rep('II', times = 10))
Sample_III <- cbind(Sample_III, c(10,20,30,40,50,60,70,80,90,100), rep('III', times = 10))
Sample_IV <- cbind(Sample_IV, c(10,20,30,40,50,60,70,80,90,100), rep('IV', times = 10))
Sample_V <- cbind(Sample_V, c(10,20,30,40,50,60,70,80,90,100), rep('V', times = 10))
data <- rbind(Sample_I, Sample_II, Sample_III, Sample_IV, Sample_V)
colnames(data) <- c("Conditional_Essential_Reactions", "Percentage_of_Data", "Sample")

data <- data.frame(Essential_Rxns = as.numeric(data[,1]), Percent = as.numeric(data[,2]), Sample = data[,3])
library(ggplot2)
ggplot(data, aes(x=Percent, y=Essential_Rxns, colour = Sample)) + geom_line(lwd=1.1) + geom_point(size=2.4) + theme_bw() +
  xlab("Percentage of Transcriptomics Integrated") + ylab("Conditionally Essential Reactions") +
  scale_x_continuous(breaks = c(20, 40, 60, 80, 100)) +
  theme(axis.text.x = element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y=element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 12))







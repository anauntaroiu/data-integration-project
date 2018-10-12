
## Set Working Directory to location of data
setwd("C:/Users/Ana/Documents/Capstone/R Code/Essential Reaction Overlap/Complete Data")

## Load Essential Reaction Lists for Each Sample and Model (Protemoics vs. Transcriptomics)
Essential_Rxns_Proteomics <- as.matrix(read.csv("Essential_Rxns_Proteomics_Models.csv"))
Essential_Rxns_Transcriptomics <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_Models.csv"))
Essential_Rxns_Unconstrained <- as.matrix(read.csv("Essential_Rxns_Unconstrained_Model.csv"))

## Load Random Essential Reaction Lists
Essential_Rxns_Proteomics_Random <- as.matrix(read.csv("Essential_Rxns_Proteomics_Random.csv"))
Essential_Rxns_Transcriptomics_Random <- as.matrix(read.csv("Essential_Rxns_Transcriptomics_Random.csv"))

## Set Working Directory to location of data
setwd("C:/Users/Ana/Documents/Capstone/R Code/Essential Reaction Overlap/Relative Data")

## Load relative data
Essential_Rxns_Relative_Proteomics <- as.matrix(read.csv("Essential_Rxns_Relative_Proteomics.csv"))
Essential_Rxns_Relative_Proteomics_Random <- as.matrix(read.csv("Essential_Rxns_Relative_Proteomics_Random.csv"))

## Set Working Directory to location of data
setwd("C:/Users/Ana/Documents/Capstone/R Code/Essential Reaction Overlap/Absolute Data")

## Load absolute data
Essential_Rxns_Absolute_Proteomics <- as.matrix(read.csv("Essential_Rxns_Absolute_Proteomics.csv"))
Essential_Rxns_Absolute_Proteomics_Random <- as.matrix(read.csv("Essential_Rxns_Absolute_Proteomics_Random.csv"))

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

###################################################################
# Call remove_unconstrained_ess_rxn to get conditionally essential reactions

## Complete Data
Essential_Rxns_Proteomics <- remove_unconstrained_ess_rxn(Essential_Rxns_Proteomics, Essential_Rxns_Unconstrained)
Essential_Rxns_Transcriptomics <- remove_unconstrained_ess_rxn(Essential_Rxns_Transcriptomics, Essential_Rxns_Unconstrained)
Essential_Rxns_Proteomics_Random <- remove_unconstrained_ess_rxn(Essential_Rxns_Proteomics_Random, Essential_Rxns_Unconstrained)
Essential_Rxns_Transcriptomics_Random <- remove_unconstrained_ess_rxn(Essential_Rxns_Transcriptomics_Random, Essential_Rxns_Unconstrained)

## Relative Proteomics Data
Essential_Rxns_Relative_Proteomics <- remove_unconstrained_ess_rxn(Essential_Rxns_Relative_Proteomics, Essential_Rxns_Unconstrained)
Essential_Rxns_Relative_Proteomics_Random <- remove_unconstrained_ess_rxn(Essential_Rxns_Relative_Proteomics_Random, Essential_Rxns_Unconstrained)

## Absolute Proteomics Data
Essential_Rxns_Absolute_Proteomics <- remove_unconstrained_ess_rxn(Essential_Rxns_Absolute_Proteomics, Essential_Rxns_Unconstrained)
Essential_Rxns_Absolute_Proteomics_Random <- remove_unconstrained_ess_rxn(Essential_Rxns_Absolute_Proteomics_Random, Essential_Rxns_Unconstrained)

####################################################################
# Function to compare reactions across data samples (I-V)
#
# INPUT: Two matrices of reactions lists
#
# OUTPUT: Matrix showing number of overlaping reactions for each sample comparision

compare_reactions <- function(Reactions_1, Reactions_2){
  # Create a matrix for pair-wise comparisons of samples
  comparisons <- matrix(data = 0, nrow = length(Reactions_1), ncol = length(Reactions_2))
  
  for (i in 1:length(Reactions_1)){ # Loop through each sample in first list of essential rxns
    
    for (j in 1:length(Reactions_2)){ # Loop through each sample in first list of essential rxns
      overlap <- intersect(Reactions_1[[i]], Reactions_2[[j]]) # Find overlapping rxns
      comparisons[i,j] <- length(overlap) # Save number of rxns
    }
  }
  # Create matrix to hold averages in rxns overlap - divide by three since there are three replicates per sample
  overlap_medians <- matrix(data = 0, nrow = length(Reactions_1)/3, ncol = length(Reactions_2)/3)
  
  index <- c(1,4,7,10,13) # Indexes where each sample's replicates begin
  for (i in 1:5){ # Loop through sample I-V (rows)
    
    for (j in 1:5){ # Loop through sample I-V (columns)
      overlap_medians[i, j] <- round(median(comparisons[index[i]:(index[i]+2), index[j]:(index[j]+2)]),
                                     digits = 1)# Record median across nine comparisons
    }
  }
  
  return(overlap_medians)
}

#################################################################
# Call compare_reaction function for each comparison
T_v_P_Complete <- compare_reactions(Essential_Rxns_Transcriptomics, Essential_Rxns_Proteomics)
T_v_P_Relative <- compare_reactions(Essential_Rxns_Transcriptomics, Essential_Rxns_Relative_Proteomics)
T_v_P_Absolute <- compare_reactions(Essential_Rxns_Transcriptomics, Essential_Rxns_Absolute_Proteomics)

##################################################################
# Function to compare reactions across samples against random
#
# INPUT: Two matrices of reactions lists
#
# OUTPUT: Matrix showing number of overlaping reactions for each sample comparision

compare_reactions_against_random <- function(Data_Reactions, Random_Reactions){
  # Create matrix to save pair-wise comparisons of overlaping rxns
  comparisons <- matrix(data = 0, nrow = length(Data_Reactions), ncol = length(Random_Reactions))
  index <- 0 # Index to track random samples (5 random sample for each replicate)
  for (i in 1:length(Data_Reactions)){ # Loop through each sample of essential rxns
    
    for (j in 1:5){ # Loop through each random sample of essential rxns (corresponding to the sample of interest)
      index <- index + 1
      overlap <- intersect(Data_Reactions[[i]], Random_Reactions[[index]]) # Find rxns overlap
      comparisons[i,index] <- length(overlap) # Record number of overlapping rxns
    }
  }
  # Matrix to record medians in overlap for each sample
  overlap_medians <- matrix(data = 0, nrow = 5, ncol = 1)
  r <- 1
  c <- 1
  for (i in 1:5){ # Loop through each sample (I-V)
   
      sample_1 <- median(comparisons[r, c:(c+4)])
      r <- r + 1; c <- c + 5
      sample_2 <- median(comparisons[r, c:(c+4)])
      r <- r + 1; c <- c + 5
      sample_3 <- median(comparisons[r, c:(c+4)])
      r <- r + 1; c <- c + 5
      overlap_medians[i,1] <- round(median(c(sample_1, sample_2, sample_3)), digits = 1)
  }
  
  return(overlap_medians)
}
##################################################################
# Call compare_reactions_against_random function
T_v_R_Complete <- compare_reactions_against_random(Essential_Rxns_Transcriptomics, Essential_Rxns_Transcriptomics_Random)
P_v_R_Complete <- compare_reactions_against_random(Essential_Rxns_Proteomics, Essential_Rxns_Proteomics_Random)
P_v_R_Relative <- compare_reactions_against_random(Essential_Rxns_Relative_Proteomics, Essential_Rxns_Relative_Proteomics_Random)
P_v_R_Absolute <- compare_reactions_against_random(Essential_Rxns_Absolute_Proteomics, Essential_Rxns_Absolute_Proteomics_Random)

##################################################################
## Create matrix to store data for proportions
#
# INPUT: Overlap of rxns between proteomics and transcriptomics; 
#        Overlap of rxns between transcriptomics and random;
#        Overlap of rxns between proteomics and random;
#        List of essential rxns for transcriptomics; List of essential rxns for proteomics
#
# OUTPUT: Matrix of proportions for figure

generate_proportions_matrix <- function(Data, Random_T, Random_P, 
                                        Essential_Rxns_Transcriptomics, Essential_Rxns_Proteomics,
                                        Essential_Rxns_Transcriptomics_Random, 
                                        Essential_Rxns_Proteomics_Random){
  Num_rxns <- matrix(data = 0, nrow = 1, ncol = 5)
  index <- 1
  for (i in 1:5){
    Num_rxns[1,i] <- median(c(length(Essential_Rxns_Transcriptomics[[index]]),
                              length(Essential_Rxns_Transcriptomics[[index+1]]),
                              length(Essential_Rxns_Transcriptomics[[index+2]])))
    index <- index + 3
  }
  
  Data <- t(Data)
  for (i in 1:5){
      Random_T[i,1] <- round(Random_T[i,1] / as.numeric(Num_rxns[1,i]), digits = 2)
      Random_P[i,1] <- round(Random_P[i,1] / as.numeric(Num_rxns[1,i]), digits = 2)
    for(j in 1:5){
      Data[j,i] <- round(as.numeric(Data[j,i]) / as.numeric(Num_rxns[1,i]), digits = 2)  
    }
  }
  
  overlap_random <- matrix(data = 0, nrow = 1, ncol = 150)
  list_length <- matrix(data = 0, nrow = 1, ncol = 150)
  for (i in 1:length(Essential_Rxns_Transcriptomics_Random)){
    overlap_random[1,i] <- length(intersect(Essential_Rxns_Transcriptomics_Random[[i]],
                                            Essential_Rxns_Proteomics_Random[[i]]))
    list_length[1,i] <- length(Essential_Rxns_Transcriptomics_Random[[i]])
  }
  
  
  Data <- rbind(t(Random_T), Data); Data <- cbind(Data, rev(rbind(Random_P, round(
    median(overlap_random)/median(list_length), digits = 2))))
  return(Data)
}

#####################################################################
## Create heatmap for Complete Proteomics and Transcriptomics

# Initialize matrix
data <- generate_proportions_matrix(T_v_P_Complete, T_v_R_Complete, P_v_R_Complete, 
                                    Essential_Rxns_Transcriptomics, Essential_Rxns_Proteomics,
                                    Essential_Rxns_Transcriptomics_Random,
                                    Essential_Rxns_Relative_Proteomics_Random)
colnames(data) <- c('I_T','II_T','III_T','IV_T','V_T', 'R')
rownames(data) <- c('R','V_P','IV_P','III_P','II_P','I_P')

library(reshape2)
melted_cormat <- melt(data, na.rm = TRUE)

# Heatmap for Complete Data
library(ggplot2)
Complete_Data <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = as.numeric(value))) +
    geom_tile(color = "white") + scale_fill_gradient2(low = "blue",
    high = "red",  mid = "grey", midpoint = 0.36, limit = c(0.15,0.56),
    space = "Lab", name="Overlap of Conditionally\nEssential Reactions") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 11, hjust = 1), axis.text.y = element_text(size = 11, vjust = 1, hjust = 1)) + coord_fixed()


Complete_Data + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face="bold", size=13),
    axis.text.y = element_text(face="bold", size=13),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    #legend.justification = c(1, 0),
    #legend.position = c(0.83, 0.76),
    legend.position = "top",
    legend.direction = "horizontal", legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
    guides(fill = guide_colorbar(barwidth = 10, barheight = 1.3,
    title.position = "top", title.hjust = 0.5))

#####################################################################
## Create heatmap for Relative Proteomics and Transcriptomics

# Initialize matrix
data <- generate_proportions_matrix(T_v_P_Relative, T_v_R_Complete, P_v_R_Relative, 
                                    Essential_Rxns_Transcriptomics, Essential_Rxns_Relative_Proteomics,
                                    Essential_Rxns_Transcriptomics_Random,
                                    Essential_Rxns_Relative_Proteomics_Random)
colnames(data) <- c('I_T','II_T','III_T','IV_T','V_T', 'R')
rownames(data) <- c('R','V_P','IV_P','III_P','II_P','I_P')

melted_cormat <- melt(data, na.rm = TRUE)

# Heatmap for Relative Data
Relative_Data <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = as.numeric(value))) +
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "blue", high = "red",  mid = "grey",
                       midpoint = 0.32, limit = c(0.08,0.56), space = "Lab", 
                       name="Overlap of Conditionally\nEssential Reactions") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1), 
                          axis.text.y = element_text(size = 11, vjust = 1, hjust = 1)) + coord_fixed()


Relative_Data + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face="bold", size=13),
    axis.text.y = element_text(face="bold", size=13),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    #legend.justification = c(1, 0),
    #legend.position = c(0.83, 0.76),
    legend.position = "top",
    legend.direction = "horizontal", legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))

#####################################################################
## Create heatmap for Absolute Proteomics and Transcriptomics

# Initialize matrix
data <- generate_proportions_matrix(T_v_P_Absolute, T_v_R_Complete, P_v_R_Absolute, 
                                    Essential_Rxns_Transcriptomics, Essential_Rxns_Absolute_Proteomics,
                                    Essential_Rxns_Transcriptomics_Random,
                                    Essential_Rxns_Relative_Proteomics_Random)
colnames(data) <- c('I_T','II_T','III_T','IV_T','V_T', 'R')
rownames(data) <- c('R','V_P','IV_P','III_P','II_P','I_P')

melted_cormat <- melt(data, na.rm = TRUE)

# Heatmap for Relative Data
Absolute_Data <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = as.numeric(value))) +
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "blue", high = "red",  mid = "grey",
                       midpoint = 0.34, limit = c(0.12,0.56), space = "Lab", 
                       name="Overlap of Conditionally\nEssential Reactions") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1), 
                          axis.text.y = element_text(size = 11, vjust = 1, hjust = 1)) + coord_fixed()


Absolute_Data + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face="bold", size=13),
    axis.text.y = element_text(face="bold", size=13),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    #legend.justification = c(1, 0),
    #legend.position = c(0.83, 0.76),
    legend.position = "top",
    legend.direction = "horizontal", legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))


##############################
## Heatmap with raw numbers ##
##############################

##################################################################
## Create matrix to store data
#
# INPUT: Overlap of rxns between proteomics and transcriptomics; 
#        Overlap of rxns between transcriptomics and random;
#        Overlap of rxns between proteomics and random;
#       
#
# OUTPUT: Matrix for figure

generate_raw_number <- function(Data, Random_T, Random_P,
                                Essential_Rxns_Transcriptomics_Random,
                                Essential_Rxns_Relative_Proteomics_Random){
  
  overlap_random <- matrix(data = 0, nrow = 1, ncol = 150)
  
  for (i in 1:length(Essential_Rxns_Transcriptomics_Random)){
    overlap_random[1,i] <- length(intersect(Essential_Rxns_Transcriptomics_Random[[i]],
                                            Essential_Rxns_Proteomics_Random[[i]]))
  }
  
  Data <- t(Data)
  Data <- rbind(t(Random_T), Data); Data <- cbind(Data, rev(rbind(Random_P, round(median(overlap_random), digits = 1))))
  return(Data)
}

#####################################################################
## Create heatmap for Complete Proteomics and Transcriptomics (Raw Numbers)

# Initialize matrix
data <- generate_raw_number(T_v_P_Complete, T_v_R_Complete, P_v_R_Complete,
                            Essential_Rxns_Transcriptomics_Random,
                            Essential_Rxns_Relative_Proteomics_Random)
colnames(data) <- c('I_T','II_T','III_T','IV_T','V_T', 'R')
rownames(data) <- c('R','V_P','IV_P','III_P','II_P','I_P')

melted_cormat <- melt(data, na.rm = TRUE)

# Heatmap for Complete Data
Complete_Data_Raw <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = as.numeric(value))) +
  geom_tile(color = "white") + scale_fill_gradient2(low = "blue",
                                                    high = "red",  mid = "grey", midpoint = 22, limit = c(7,37),
                                                    space = "Lab", name="Overlap of Conditionally\nEssential Reactions") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                                     size = 11, hjust = 1), axis.text.y = element_text(size = 11, vjust = 1, hjust = 1)) + coord_fixed()


Complete_Data_Raw + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face="bold", size=13),
    axis.text.y = element_text(face="bold", size=13),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    #legend.justification = c(1, 0),
    #legend.position = c(0.83, 0.76),
    legend.position = "top",
    legend.direction = "horizontal", legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))

#####################################################################
## Create heatmap for Relative Proteomics and Transcriptomics (Raw Numbers)

# Initialize matrix
data <- generate_raw_number(T_v_P_Relative, T_v_R_Complete, P_v_R_Relative,
                            Essential_Rxns_Transcriptomics_Random,
                            Essential_Rxns_Relative_Proteomics_Random)
colnames(data) <- c('I_T','II_T','III_T','IV_T','V_T', 'R')
rownames(data) <- c('R','V_P','IV_P','III_P','II_P','I_P')

melted_cormat <- melt(data, na.rm = TRUE)

# Heatmap for Relative Data
Relative_Data_Raw <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = as.numeric(value))) +
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "blue", high = "red",  mid = "grey",
                       midpoint = 21, limit = c(4,37), space = "Lab", 
                       name="Overlap of Conditionally\nEssential Reactions") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1), 
                          axis.text.y = element_text(size = 11, vjust = 1, hjust = 1)) + coord_fixed()


Relative_Data_Raw + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face="bold", size=13),
    axis.text.y = element_text(face="bold", size=13),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    #legend.justification = c(1, 0),
    #legend.position = c(0.83, 0.76),
    legend.position = "top",
    legend.direction = "horizontal", legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))

#####################################################################
## Create heatmap for Absolute Proteomics and Transcriptomics (Raw Number)

# Initialize matrix
data <- generate_raw_number(T_v_P_Absolute, T_v_R_Complete, P_v_R_Absolute,
                            Essential_Rxns_Transcriptomics_Random,
                            Essential_Rxns_Relative_Proteomics_Random)
colnames(data) <- c('I_T','II_T','III_T','IV_T','V_T', 'R')
rownames(data) <- c('R','V_P','IV_P','III_P','II_P','I_P')

melted_cormat <- melt(data, na.rm = TRUE)

# Heatmap for Relative Data
Absolute_Data_Raw <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = as.numeric(value))) +
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "blue", high = "red",  mid = "grey",
                       midpoint = 22, limit = c(7,37), space = "Lab", 
                       name="Overlap of Conditionally\nEssential Reactions") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1), 
                          axis.text.y = element_text(size = 11, vjust = 1, hjust = 1)) + coord_fixed()


Absolute_Data_Raw + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(face="bold", size=13),
    axis.text.y = element_text(face="bold", size=13),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    #legend.justification = c(1, 0),
    #legend.position = c(0.83, 0.76),
    legend.position = "top",
    legend.direction = "horizontal", legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))


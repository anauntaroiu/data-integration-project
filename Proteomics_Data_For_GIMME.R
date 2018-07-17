
setwd("C:/Users/Ana/Documents/Capstone/R Code/Proteomics") # set working directory to were the file can be found

# Load data
Rel_Pro <- as.matrix(read.csv("Combined_Relative_Proteomics.csv"))
Rel_Pro <- Rel_Pro[,2:17] # Delete column with no useful data

Abs_Pro <- as.matrix(read.csv("Combined_Absolute_Proteomics.csv"))
Abs_Pro <- Abs_Pro[,2:17] # Delete column with no useful data

#############################################################
# Split replicates from each sample

## Relative Data
## Sample I
Sample_I_R_I <- Rel_Pro[,1:2]
Sample_I_R_II <- cbind(Rel_Pro[,1], Rel_Pro[,3])
Sample_I_R_III <- cbind(Rel_Pro[,1], Rel_Pro[,4])

## Sample II
Sample_II_R_I <- cbind(Rel_Pro[,1], Rel_Pro[,5])
Sample_II_R_II <- cbind(Rel_Pro[,1], Rel_Pro[,6])
Sample_II_R_III <- cbind(Rel_Pro[,1], Rel_Pro[,7])

## Sample III
Sample_III_R_I <- cbind(Rel_Pro[,1], Rel_Pro[,8])
Sample_III_R_II <- cbind(Rel_Pro[,1], Rel_Pro[,9])
Sample_III_R_III <- cbind(Rel_Pro[,1], Rel_Pro[,10])

## Sample IV
Sample_IV_R_I <- cbind(Rel_Pro[,1], Rel_Pro[,11])
Sample_IV_R_II <- cbind(Rel_Pro[,1], Rel_Pro[,12])
Sample_IV_R_III <- cbind(Rel_Pro[,1], Rel_Pro[,13])

## Sample V
Sample_V_R_I <- cbind(Rel_Pro[,1], Rel_Pro[,14])
Sample_V_R_II <- cbind(Rel_Pro[,1], Rel_Pro[,15])
Sample_V_R_III <- cbind(Rel_Pro[,1], Rel_Pro[,16])

##################################################
# Create function to remove genes with no expression values

## INPUTS: Three replicates from one sample

## OUTPUTS: A list with the three input replicates, where genes with missing expression values are removed

remove_missing_expression <- function(Replicate_I, Replicate_II, Replicate_III ){
  
  ## Replicate 1
  genes_to_remove <- matrix(data = NA, nrow = length(Replicate_I), ncol = 1) # Create matrix for genes to remove
  index <- 1 # Matrix index for genes_to_remove matrix
  
  for (i in 1:dim(Replicate_I)[1]){ # Loop through each gene in sample
    
    if (is.na(Replicate_I[i,2])){ # If no expression value exists
      
      genes_to_remove[index,1] <- i # add gene index to delete it later
      index <- index + 1 # Increment index by one
    }
  }
  Replicate_I <- as.matrix(Replicate_I[-c(genes_to_remove[1:(index-1),1]),]) # Delete genes with no expression values
  
  # Transform log2 format of data into log10 format
  Replicate_I[,2] <- log10(2^as.numeric(Replicate_I[,2]))
  
  ## Replicate 2
  genes_to_remove <- matrix(data = NA, nrow = length(Replicate_II), ncol = 1) # Create matrix for genes to remove
  index <- 1 # Matrix index for genes_to_remove matrix
  
  for (i in 1:dim(Replicate_II)[1]){ # Loop through each gene in sample
    
    if (is.na(Replicate_II[i,2])){ # If no expression value exists
      
      genes_to_remove[index,1] <- i # add gene index to delete it later
      index <- index + 1 # Increment index by one
    }
  }
  Replicate_II <- as.matrix(Replicate_II[-c(genes_to_remove[1:(index-1),1]),]) # Delete genes with no expression values
  
  # Transform log2 format of data into log10 format
  Replicate_II[,2] <- log10(2^as.numeric(Replicate_II[,2]))
  
  ## Replicate 3
  genes_to_remove <- matrix(data = NA, nrow = length(Replicate_III), ncol = 1) # Create matrix for genes to remove
  index <- 1 # Matrix index for genes_to_remove matrix
  
  for (i in 1:dim(Replicate_III)[1]){ # Loope through each gene in sample
    
    if (is.na(Replicate_III[i,2])){ # If no expression value exists
      
      genes_to_remove[index,1] <- i # Add gene index to delete it later
      index <- index + 1 # Increment index by one 
    }
  }
  Replicate_III <- as.matrix(Replicate_III[-c(genes_to_remove[1:(index-1),1]),]) # Delete genes with no expression values
  
  # Transform log2 format of data into log10 format
  Replicate_III[,2] <- log10(2^as.numeric(Replicate_III[,2])) 
  
  ## Store replicates in a list
  All_Replicates <- list(Replicate_I, Replicate_II, Replicate_III)
  
  return(All_Replicates)
}

##################################################
# Invoke remove_missing_expression function for all samples

Sample_I_R <- remove_missing_expression(Sample_I_R_I, Sample_I_R_II, Sample_I_R_III)
rm(Sample_I_R_I, Sample_I_R_II, Sample_I_R_III) # Remove objects that are no longer needed

Sample_II_R <- remove_missing_expression(Sample_II_R_I, Sample_II_R_II, Sample_II_R_III)
rm(Sample_II_R_I, Sample_II_R_II, Sample_II_R_III) # Remove objects that are no longer needed

Sample_III_R <- remove_missing_expression(Sample_III_R_I, Sample_III_R_II, Sample_III_R_III)
rm(Sample_III_R_I, Sample_III_R_II, Sample_III_R_III) # Remove objects that are no longer needed

Sample_IV_R <- remove_missing_expression(Sample_IV_R_I, Sample_IV_R_II, Sample_IV_R_III)
rm(Sample_IV_R_I, Sample_IV_R_II, Sample_IV_R_III) # Remove objects that are no longer needed

Sample_V_R <- remove_missing_expression(Sample_V_R_I, Sample_V_R_II, Sample_V_R_III)
rm(Sample_V_R_I, Sample_V_R_II, Sample_V_R_III) # Remove objects that are no longer needed

###################################################
# Split replicates from each sample

## Absolute Data
## Sample I
Sample_I_A_I <- Abs_Pro[,1:2]
Sample_I_A_II <- cbind(Abs_Pro[,1], Abs_Pro[,3])
Sample_I_A_III <- cbind(Abs_Pro[,1], Abs_Pro[,4])

## Sample II
Sample_II_A_I <- cbind(Abs_Pro[,1], Abs_Pro[,5])
Sample_II_A_II <- cbind(Abs_Pro[,1], Abs_Pro[,6])
Sample_II_A_III <- cbind(Abs_Pro[,1], Abs_Pro[,7])

## Sample III
Sample_III_A_I <- cbind(Abs_Pro[,1], Abs_Pro[,8])
Sample_III_A_II <- cbind(Abs_Pro[,1], Abs_Pro[,9])
Sample_III_A_III <- cbind(Abs_Pro[,1], Abs_Pro[,10])

## Sample IV
Sample_IV_A_I <- cbind(Abs_Pro[,1], Abs_Pro[,11])
Sample_IV_A_II <- cbind(Abs_Pro[,1], Abs_Pro[,12])
Sample_IV_A_III <- cbind(Abs_Pro[,1], Abs_Pro[,13])

## Sample V
Sample_V_A_I <- cbind(Abs_Pro[,1], Abs_Pro[,14])
Sample_V_A_II <- cbind(Abs_Pro[,1], Abs_Pro[,15])
Sample_V_A_III <- cbind(Abs_Pro[,1], Abs_Pro[,16])

##################################################
# Invoke remove_missing_expression function for all samples

Sample_I_A <- remove_missing_expression(Sample_I_A_I, Sample_I_A_II, Sample_I_A_III)
rm(Sample_I_A_I, Sample_I_A_II, Sample_I_A_III) # Remove objects that are no longer needed

Sample_II_A <- remove_missing_expression(Sample_II_A_I, Sample_II_A_II, Sample_II_A_III)
rm(Sample_II_A_I, Sample_II_A_II, Sample_II_A_III) # Remove objects that are no longer needed

Sample_III_A <- remove_missing_expression(Sample_III_A_I, Sample_III_A_II, Sample_III_A_III)
rm(Sample_III_A_I, Sample_III_A_II, Sample_III_A_III) # Remove objects that are no longer needed

Sample_IV_A <- remove_missing_expression(Sample_IV_A_I, Sample_IV_A_II, Sample_IV_A_III)
rm(Sample_IV_A_I, Sample_IV_A_II, Sample_IV_A_III) # Remove objects that are no longer needed

Sample_V_A <- remove_missing_expression(Sample_V_A_I, Sample_V_A_II, Sample_V_A_III)
rm(Sample_V_A_I, Sample_V_A_II, Sample_V_A_III) # Remove objects that are no longer needed

############################################
# Create a function to normalize the data based on Z-scores

# INPUTS: One sample (Absolute/Relative)

# OUTPUTS: Sample normalized by Z-scores

normalize <- function(Sample){
  normalized_sample <- list() # Create list to hold normalized replicates
  
  for(i in 1:3){ # Loop through each replicate
    
    stan_dev <- sd(Sample[[i]][,2])*(sqrt((length(Sample[[i]][,2])-1)/(length(Sample[[i]][,2])))) # Calcluate Pop. Stand. Dev.
    samp_mean <- mean(as.numeric(Sample[[i]][,2])) # Calculate data mean
    z_scores <- ((as.numeric(Sample[[i]][,2]) - samp_mean)/stan_dev) # Calcluate Z-Scores for each data point
    z_scores <- z_scores - min(z_scores) # Shift data to prevent negative numbers
    normalized_sample[[i]] <- cbind(Sample[[i]][,1], z_scores) # Save Z-scores with gene IDs in list
  }
  
  return(normalized_sample)
}

#############################################
# Invoke normalize function for all samples

Sample_I_A <- normalize(Sample_I_A)
Sample_II_A <- normalize(Sample_II_A)
Sample_III_A <- normalize(Sample_III_A)
Sample_IV_A <- normalize(Sample_IV_A)
Sample_V_A <- normalize(Sample_V_A)

Sample_I_R <- normalize(Sample_I_R)
Sample_II_R <- normalize(Sample_II_R)
Sample_III_R <- normalize(Sample_III_R)
Sample_IV_R <- normalize(Sample_IV_R)
Sample_V_R <- normalize(Sample_V_R)

#############################################
## Combine Absolute and Relative Data
## Sample I
Exprs_matrix_I <- list()
for (i in 1:3){
  Exprs_matrix_I[[i]] <- rbind(Sample_I_A[[i]], Sample_I_R[[i]])
}
rm(Sample_I_A, Sample_I_R) # Remove no longer needed variables

## Sample II
Exprs_matrix_II <- list()
for (i in 1:3){
  Exprs_matrix_II[[i]] <- rbind(Sample_II_A[[i]], Sample_II_R[[i]])
}
rm(Sample_II_A, Sample_II_R) # Remove no longer needed variables

## Sample III
Exprs_matrix_III <- list()
for (i in 1:3){
  Exprs_matrix_III[[i]] <- rbind(Sample_III_A[[i]], Sample_III_R[[i]])
}
rm(Sample_III_A, Sample_III_R) # Remove no longer needed variables

## Sample IV
Exprs_matrix_IV <- list()
for (i in 1:3){
  Exprs_matrix_IV[[i]] <- rbind(Sample_IV_A[[i]], Sample_IV_R[[i]])
}
rm(Sample_IV_A, Sample_IV_R) # Remove no longer needed variables

## Sample V
Exprs_matrix_V <- list()
for (i in 1:3){
  Exprs_matrix_V[[i]] <- rbind(Sample_V_A[[i]], Sample_V_R[[i]])
}
rm(Sample_V_A, Sample_V_R) # Remove no longer needed variables

##############################################

## import all reactions in the model and list of genes with matching reactions
Reactions <- as.matrix(read.csv("Reactions.csv")) ## All reactions found in model
genes_and_reactions <- as.matrix(read.csv("genes_and_rxns.csv")) # list of genes with corresponding rxns

# Load required package for generating random distributions
library(fGarch)
library(sn)

## Function to generate numbers for random models
# INPUT: Sample data list (includes all three replicates)

# OUTPUT: Data frame with ten random samples for each replicate

generate_random_samples <- function(sample){
  ## Generate 10 random samples for each replicate

  # Generate empty data frame to store random samples
  random_samples <- matrix(data = NA, ncol = 30, nrow = length(Reactions)) 
  index <- 1 # Create to store index of columns in random_samples matrix

  for (j in 1:3){ # Loop through each replicate
  
  for (i in 1:10){ # Loop to generate 10 random replicates
    
    # Sample skewed random distribution to generate random samples
        ## rsnorm function
        ## n: sample size
        ## xi (location parameter): set as median of data
        ## omega (scale parameter): set as standard deviation
        ## alpha (slant parameter): set as skewness of data
        ## NOTE: take absolute value of data to ensure all data is positive
    random_samples[,index] <- abs(rsn(n = length(Reactions), 
                                     xi = median(as.numeric(sample[[j]][,2])), omega = sd(sample[[j]][,2]),
                                     alpha = skewness(as.numeric(sample[[j]][,2]))))
    
    index <- index + 1 # Increase column index by one every time loop runs
   }
  }
  return(random_samples)
}

# Generate matrix with 10 random samples for each replicate by calling generate_random_samples function
random_samples <- cbind(generate_random_samples(Exprs_matrix_I),
                        generate_random_samples(Exprs_matrix_II),
                        generate_random_samples(Exprs_matrix_III),
                        generate_random_samples(Exprs_matrix_IV),
                        generate_random_samples(Exprs_matrix_V))

# Label columns based on sample number
colnames(random_samples) <- c(rep('I', 30), rep('II', 30), rep('III', 30), rep('IV', 30), rep('V', 30))
  
###########################################################

## Find the corresponding reactions for genes
## in order to integrate data with GIMME

matching_reactions <- matrix(nrow = length(Reactions), ncol = 1) # Matrix to save gene names to rxns from model
for (i in 1:length(Reactions)){
  hit <- grep(Reactions[i], genes_and_reactions[,2]) # Search rxns/gene list to find rxn in model
  if (length(hit) != 0){
    matching_reactions[i] <- genes_and_reactions[hit[1],1] # Save gene name found above
  }else{
    matching_reactions[i] <- -1 # If no gene ID, use '-1' as a placeholder
    for (j in 1:150){
      random_samples[i,j] <- -1 # Add place holder to random samples
    }
  }
}
genes_and_reactions <- cbind(Reactions, matching_reactions)

############################################################
## Function to match proteomics data to reactions 
## using gene and reaction list 

# INPUT: Sample of expression data; Position where samples begin in random matrix for this data sample (ex. 1-150);
#        Matrix of random data; List of genes with matching reactions

# OUTPUT: Matrix of expression data matched to model reaction list

match_expression <- function(sample, index, random_samples, genes_and_reactions){
  # create matrix to hold sorted data
  matching_expression <- matrix(nrow = length(Reactions), ncol = 3)
  
   for (x in 1:3){ # Run loop for each replicate
   
    for (i in 1:dim(genes_and_reactions)[1]){ # Run loop for each reaction in model
    
      if (genes_and_reactions[i,2] == -1){ # If no gene exists for reaction, put placeholder
        
       matching_expression[i, x] <- -1 # Add placeholder
    }else{ # If a matching gene/rxn exists, assign expression value to the rxn ID
      
      hit <- grep(genes_and_reactions[i,2], sample[[x]][,1]) # Find rxn ID in expression dataset
      
      if (length(hit) == 0){ # If expression was not measure for the reaction
        
        # Add placeholder in expression data and random data
        matching_expression[i, x] <- -1
        random_samples[i, index:(index+9)] <- -1 # use index variable to modify correct samples
        
      }else{ # If expression value was measured for rxn, assign value to rxn ID
        
        matching_expression[i, x] <- sample[[x]][hit[1],2]
      }
    }
    }
     # Increase index for the next set of random samples
     index <- index + 10
  }
  # return matched expression and random data as a list
  combined_data <- list(matching_expression, random_samples)
  return(combined_data)
}

#############################################################
## Call match_expression function for all samples
Data_I <- match_expression(Exprs_matrix_I, 1, random_samples, genes_and_reactions)

# Feed returned random samples into match_expression function, so it can be modified
# for every repilcate
Data_II <- match_expression(Exprs_matrix_II, 31, Data_I[[2]], genes_and_reactions)
Data_III <- match_expression(Exprs_matrix_III, 61, Data_II[[2]], genes_and_reactions)
Data_IV <- match_expression(Exprs_matrix_IV, 91, Data_III[[2]], genes_and_reactions)
Data_V <- match_expression(Exprs_matrix_V, 121, Data_IV[[2]], genes_and_reactions)

# Save random data returned from last sample
random_samples <- Data_V[[2]]

# Save only expression data
Data_I <- Data_I[[1]]
Data_II <- Data_II[[1]]
Data_III <- Data_III[[1]]
Data_IV <- Data_IV[[1]]
Data_V <- Data_V[[1]]

#############################################
# Function to find medians of data
give_medians <- function(Data){
  
  # create matrix to store medians for all repicates
  medians <- matrix(data = NA, nrow = 1, ncol = 3)
  
  for (i in 1:3){ # Loop through each replicate
    
    missing_values <- grep('-1', Data[,i]) # Find positions of missing values
    
    # remove missing values and find median
    medians[1,i] <- median(as.numeric(Data[-c(missing_values),i]))
    rm(missing_values) # Clear variables each time loop runs
  }
  return(medians)
}

#############################################
# Call give_medians function for each variable
median_I <- give_medians(Data_I)
median_II <- give_medians(Data_II)
median_III <- give_medians(Data_III)
median_IV <- give_medians(Data_IV)
median_V <- give_medians(Data_V)

# Combine medians for each sample into one variable
medians_data <- rbind(median_I, median_II, median_III, median_IV, median_V)
colnames(medians_data) <- c("Replicate_1", "Replicate_2", "Replicate_3")

#############################################
## Find medians of random data
#
# Create a variable to hold medians
medians_random <- matrix(data = NA, nrow = 150, ncol = 1)
colnames(medians_random) <- "Medians"
rownames(medians_random) <- colnames(random_samples)

for (i in 1:dim(random_samples)[2]){ # Loop through each random sample
  
  # Find positions in sample with no expression value (ex. '-1')
  missing_values <- grep('-1', random_samples[,i])
  
  # Exclude missing values and find medians of expression values
  medians_random[i,1] <- median(random_samples[-c(missing_values),i])
  rm(missing_values) # Remove positions of missing values every time loop runs
}

##############################################
# Save Data (Protein Expression for Corresponding Reactions in the Model)

setwd("C:/Users/Ana/Documents/MATLAB/Capstone/Combined Proteomics") # set working directory to were the file can be found

# Save expression data as CSV files
colnames(Data_I) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_I, file = "Proteomics_Data_For_GIMME_Sample_I.csv", row.names = FALSE)
colnames(Data_II) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_II, file = "Proteomics_Data_For_GIMME_Sample_II.csv", row.names = FALSE)
colnames(Data_III) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_III, file = "Proteomics_Data_For_GIMME_Sample_III.csv", row.names = FALSE)
colnames(Data_IV) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_IV, file = "Proteomics_Data_For_GIMME_Sample_IV.csv", row.names = FALSE)
colnames(Data_V) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_V, file = "Proteomics_Data_For_GIMME_Sample_V.csv", row.names = FALSE)

# Save medians of expression data
write.csv(medians_data, file = "Medians_for_Complete_Proteomics.csv", row.names = FALSE)

# Save data for Random Models
write.csv(random_samples, file = "Random_Data_For_GIMME_Proteomics.csv", row.names = FALSE)

# Save medians of random samples
write.csv(medians_random, file = "Medians_for_Proteomics_Random.csv", row.names = FALSE)


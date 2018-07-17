
## Set working directory to were the file can be found
setwd("C:/Users/Ana/Documents/Capstone/R Code/Transcriptomics/Samples") 

## Load Samples
Sample_I_1 <- as.matrix(read.csv("Sample I - 1.csv"))
Sample_I_2 <- as.matrix(read.csv("Sample I - 2.csv"))
Sample_I_3 <- as.matrix(read.csv("Sample I - 3.csv"))

Sample_II_1 <- as.matrix(read.csv("Sample II - 1.csv"))
Sample_II_2 <- as.matrix(read.csv("Sample II - 2.csv"))
Sample_II_3 <- as.matrix(read.csv("Sample II - 3.csv"))

Sample_III_1 <- as.matrix(read.csv("Sample III - 1.csv"))
Sample_III_2 <- as.matrix(read.csv("Sample III - 2.csv"))
Sample_III_3 <- as.matrix(read.csv("Sample III - 3.csv"))

Sample_IV_1 <- as.matrix(read.csv("Sample IV - 1.csv"))
Sample_IV_2 <- as.matrix(read.csv("Sample IV - 2.csv"))
Sample_IV_3 <- as.matrix(read.csv("Sample IV - 3.csv"))

Sample_V_1 <- as.matrix(read.csv("Sample V - 1.csv"))
Sample_V_2 <- as.matrix(read.csv("Sample V - 2.csv"))
Sample_V_3 <- as.matrix(read.csv("Sample V - 3.csv"))

## Load Probe Information
Probe_names <- as.matrix(read.csv("Platform Information.csv"))

##############################################################
# Unlog expression values and save replicates in one matrix
Expression_Data <- cbind(Probe_names, 10^as.numeric(Sample_I_1),10^as.numeric(Sample_I_2), 
          10^as.numeric(Sample_I_3), 10^as.numeric(Sample_II_1), 10^as.numeric(Sample_II_2),
          10^as.numeric(Sample_II_3), 10^as.numeric(Sample_III_1), 10^as.numeric(Sample_III_2),
          10^as.numeric(Sample_III_3), 10^as.numeric(Sample_IV_1), 10^as.numeric(Sample_IV_2),
          10^as.numeric(Sample_IV_3), 10^as.numeric(Sample_V_1), 10^as.numeric(Sample_V_2), 
          10^as.numeric(Sample_V_3))

# Label columns by sample and replicate number
colnames(Expression_Data) <- c("Gene.ID", "I_1", "I_2", "I_3",
                               "II_1", "II_2", "II_3",
                               "III_1", "III_2", "III_3",
                               "IV_1", "IV_2", "IV_3",
                               "V_1", "V_2", "V_3")

# Remove variables that are no longer needed
rm(Sample_I_1, Sample_I_2, Sample_I_3, Sample_II_1, Sample_II_2,
   Sample_II_3, Sample_III_1, Sample_III_2, Sample_III_3, Sample_IV_1, 
   Sample_IV_2, Sample_IV_3, Sample_V_1, Sample_V_2, Sample_V_3)

###############################################################
# Delete Probes with no gene ID information

index <- matrix() # Create a matrix to save expression values with missing probe IDs
count <- 1 # Used to track positions in index matrix

for (i in 1:dim(Expression_Data)[1]){ # Loop through each probe ID
  
  if (Expression_Data[i,1] == " "){ # If no probe ID exists
    
    index[count] <- i # Record position of missing probe
    count <- count + 1 # Increase counter for index matrix by one
  }
}
Expression_Data <- Expression_Data[-index,] # Remove expression values with no probe IDs

################################################################
## Replace multiple repeating probes with the probe with the largest variance

hits <- matrix() # Records probe IDs that repeat in expression data
count <- 1 # Used to track positions in hits matrix

for (i in 1:length(Probe_names)){ # Loop through each Probe ID
  
  # reset variables for each probe ID
  hits <- matrix()
  count <- 1
  
  # This loop finds the probes that repeat in the list of probe names 
  for (j in 1:dim(Expression_Data)[1]){ # Loop through gene IDs measured in data
    
    if (Probe_names[i] == Expression_Data[j,1]){ # If probe ID is found in expression data
      
      hits[count] <- j # Record the position of the matching probe ID
      count <- count + 1 # Increase counter for hits matrix by one
    }
  }
  
  if (length(hits) > 1){ # If the probe repeats
    # Create a variable to store variance across replicates for each repeat probe
    variances <- matrix(nrow = length(hits), ncol = 1);
    
    # Find the variance of each probe that is repeated 
    for (k in 1:length(hits)){ # Loop through the repeating probes
      
      variances[k,1] <- var(Expression_Data[k,2:16]) # Record variance across replicates
    }
    
    ## Find position of probe with the largest variance
    max_variance <- max(variances) # Record the largest variance
    index <- 1 # Used to store position of probe with largest variance in loop
    
    for (k in 1:length(variances)){ # Loop through each repeating probe
      
      if (max_variance == variances[k,1]){ # If probe has the largest variance
        index <- k # Record position
      }
    }
    
    # Place the expression values for the probe with the highest variance in the first probe in data
    Expression_Data[hits[1],2:16] <- Expression_Data[hits[k],2:16]
    
    # Delete the rest of the repeating probes
    Expression_Data <- Expression_Data[-c(hits[2:length(hits)]),]
  }
}

###################################################################
## Generate random samples for each replicate in data 

## import all reactions in the model and list of genes with matching reactions
Reactions <- as.matrix(read.csv("Reactions.csv")) ## All reactions found in model
genes_and_reactions <- as.matrix(read.csv("genes_and_rxns.csv")) # list of genes with corresponding rxns

# Load required package for generating random distributions
library(fGarch)

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
                                          xi = median(as.numeric(sample[,j])), omega = sd(sample[,j]),
                                          alpha = skewness(as.numeric(sample[,j]))))
      
      index <- index + 1 # Increase column index by one every time loop runs
    }
  }
  return(random_samples)
}

# Generate matrix with 10 random samples for each replicate by calling generate_random_samples function
random_samples <- cbind(generate_random_samples(Expression_Data[,2:4]), # Sample I
                        generate_random_samples(Expression_Data[,5:7]), # Sample II
                        generate_random_samples(Expression_Data[,8:10]), # Sample III
                        generate_random_samples(Expression_Data[,11:13]), # Sample IV
                        generate_random_samples(Expression_Data[,14:16])) # Sample V

# Label columns based on sample number
colnames(random_samples) <- c(rep('I', 30), rep('II', 30), rep('III', 30), rep('IV', 30), rep('V', 30))

#########################################################

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
  
  # Save gene IDs and expression values
  IDs <- sample[,1] # gene IDs
  expression <- sample[,2:4] # Expression values
  
  for (x in 1:3){ # Run loop for each replicate
    
    for (i in 1:dim(genes_and_reactions)[1]){ # Run loop for each reaction in model
      
      if (genes_and_reactions[i,2] == -1){ # If no gene exists for reaction, put placeholder
        
        matching_expression[i, x] <- -1 # Add placeholder
      }else{ # If a matching gene/rxn exists, assign expression value to the rxn ID
        
        hit <- grep(genes_and_reactions[i,2], IDs) # Find rxn ID in expression dataset
        
        if (length(hit) == 0){ # If expression was not measure for the reaction
          
          # Add placeholder in expression data and random data
          matching_expression[i, x] <- -1
          random_samples[i, index:(index+9)] <- -1 # use index variable to modify correct samples
          
        }else{ # If expression value was measured for rxn, assign value to rxn ID
          
          matching_expression[i, x] <- expression[hit[1],x]
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

## Call match_expression function for all samples
Data_I <- match_expression(Expression_Data[,1:4], 1, random_samples, genes_and_reactions)

# Feed returned random samples into match_expression function, so it can be modified
# for every repilcate
Data_II <- match_expression(Expression_Data[,c(1,5:7)], 31, Data_I[[2]], genes_and_reactions)
Data_III <- match_expression(Expression_Data[,c(1,8:10)], 61, Data_II[[2]], genes_and_reactions)
Data_IV <- match_expression(Expression_Data[,c(1,11:13)], 91, Data_III[[2]], genes_and_reactions)
Data_V <- match_expression(Expression_Data[,c(1,14:16)], 121, Data_IV[[2]], genes_and_reactions)

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


setwd("C:/Users/Ana/Documents/MATLAB/Capstone/Transcriptomics") # set working directory to were the file can be found

# Save expression data as CSV files
colnames(Data_I) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_I, file = "Transcriptomics_Data_For_GIMME_Sample_I.csv", row.names = FALSE)
colnames(Data_II) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_II, file = "Transcriptomics_Data_For_GIMME_Sample_II.csv", row.names = FALSE)
colnames(Data_III) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_III, file = "Transcriptomics_Data_For_GIMME_Sample_III.csv", row.names = FALSE)
colnames(Data_IV) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_IV, file = "Transcriptomics_Data_For_GIMME_Sample_IV.csv", row.names = FALSE)
colnames(Data_V) <- c("Replicate_1", "Replicate_2", "Replicate_3")
write.csv(Data_V, file = "Transcriptomics_Data_For_GIMME_Sample_V.csv", row.names = FALSE)

# Save medians of expression data
write.csv(medians_data, file = "Medians_for_Complete_Transcriptomics.csv", row.names = FALSE)

# Save data for Random Models
write.csv(random_samples, file = "Random_Data_For_GIMME_Transcriptomics.csv", row.names = FALSE)

# Save medians of random samples
write.csv(medians_random, file = "Medians_for_Transcriptomics_Random.csv", row.names = FALSE)


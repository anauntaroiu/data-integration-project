
# set working directory to were the file can be found
setwd("C:/Users/Ana/Documents/Capstone/R Code/Proteomics")

# Load data
Cytosol_Data <- as.matrix(read.csv("Relative Proteomics - Cytosol.csv"))
Membrane_Data <- as.matrix(read.csv("Relative Proteomics - Membrane.csv"))

###########################################################
# Add proteins measured in membrane relative proteomics data if not included
# in the cytosol relative proteomics data

Additional_Proteins <- matrix(data = NA, nrow = 1000, ncol = 16) # Generate matrix to store data
index <- 1 # Used to track positions in matrix for data storage
counter <- 0 # Used to track whether protein is found in cytosol data

for (i in 1:dim(Membrane_Data)[1]){ # Loop through each protein measured in membrane data
  
  for (j in 1:dim(Cytosol_Data)[1]){ # Loop through each protein measured in cytosol data
    
    if (agrepl(Membrane_Data[i,1], Cytosol_Data[j,1], max.distance = 0.0)){ # Check if proteins in membrane data are in cytosol data
      counter  <- 1 # Protein is found in membrane and cytosol data
    }
    
    if (j == dim(Cytosol_Data)[1] && counter == 0){ # if you reach the end of membrane data and no match, add protein
      Additional_Proteins[index,1:16] <- Membrane_Data[i,1:16] # Add membrane proteomics data point to data matrix
      index <- index + 1
    }
    
    if (j == dim(Cytosol_Data)[1] && counter == 1){
      counter <- 0 # Reset counter to default before entire loop runs again
    }
  }
}

# Add additional proteins in membrane data to cytosol data
Data <- rbind(Cytosol_Data, Additional_Proteins[1:(index-1),])

##########################################################
# Convert data from log2 scale to linear scale
for (i in 2:dim(Data)[2]){ # Loop through each column in data - skip gene ID column
  
  for (j in 1:dim(Data)[1]){ # Loop through each row in data
    
    if (is.na(Data[j,i])){ # Check if data point was measured
      Data[j,i] <- -1 # Use placeholder for missing values
      
    }else{ # If data point was measure, convert scale
      Data[j,i] <- 2^as.numeric(Data[j,i]) # Convert from log2 scale to linear scale
    }
  }
}

############################################################

## import all reactions in the model and list of genes with matching reactions
Reactions <- as.matrix(read.csv("Reactions.csv")) ## All reactions found in model
genes_and_reactions <- as.matrix(read.csv("genes_and_rxns.csv")) # list of genes with corresponding rxns

# Generate random data for relative data
############################################################
## Generate 10 random samples for each replicate

# Load required packages
library(sn)
library(fGarch)

# Generate empty data frame to store random samples
random_samples <- matrix(data = NA, ncol = 150, nrow = length(Reactions)) 
index <- 1 # Create to store index of columns in random_samples matrix

# Create variable to save medians of replicates/samples
medians <- matrix(data = NA, ncol = 1, nrow = 15)

for (j in 2:16){ # Loop through each replicate
  
  # Extract only measured values
  measured_values <- as.numeric(Data[grep('-1', Data[,j], invert = T), j])
  
  # Save medians of extracted data for GIMME integration
  medians[(j-1), 1] <- median(as.numeric(Data[grep('-1', Data[,j], invert = T), j]))
  
  for (i in 1:10){ # Loop to generate 10 random replicates
    
    # Sample skewed random distribution to generate random samples
        ## rsnorm function
        ## n: sample size
        ## xi (location parameter): set as median of data
        ## omega (scale parameter): set as standard deviation
        ## alpha (slant parameter): set as skewness of data
        ## NOTE: take absolute value of data to ensure all data is positive
    random_samples[,index] <- abs(rsn(n = length(Reactions), 
                                         xi = median(as.numeric(measured_values)), omega = sd(measured_values),
                                         alpha = skewness(as.numeric(measured_values))))
   
    index <- index + 1 # Increase column index by one every time loop runs
  }
}

####################################################################

## Find the corresponding reactions for genes
## in order to integrate data with GIMME

matching_reactions <- matrix(nrow = length(Reactions), ncol = 1) # Matrix to save gene names to rxns from model

for (i in 1:length(Reactions)){ # Loop through reactions in model
  
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
# Combine to make a list of reactions with matching gene IDs
genes_and_reactions <- cbind(Reactions, matching_reactions) 

#######################################################################
## Use gene and reaction list to match proteomics data to model reactions

# Create a matrix to store expression values matched to reactions in model
matching_expression <- matrix(nrow = length(Reactions), ncol = 15)
index <- 1 # matrix index for random samples

for (x in 2:16){ # Run loop for each replicate and sample of data matrix
  
  for (i in 1:dim(genes_and_reactions)[1]){ # Run loop for each reaction in model
    
    if (genes_and_reactions[i,2] == -1){ # If no gene ID exists for reaction, put placeholder
      matching_expression[i, (x-1)] <- -1 # Save placeholder for expression value
      random_samples[i, index:(index + 9)] <- -1 # Put placeholder in random samples data as well
      
    }else{ # If gene ID exists for reaction
      hit <- grep(genes_and_reactions[i,2], Data[,1]) # Determine whether expression was measured for a reaction
      
      if (length(hit) == 0 || Data[hit[1],x] == -1){ # If expression was not measured for a reaction
        
        matching_expression[i, (x-1)] <- -1 # Put placeholder for expression value
        random_samples[i, index:(index + 9)] <- -1 # Put placeholder in random samples data as well
        
      }else{ # Expression values exists
        matching_expression[i, (x-1)] <- as.numeric(Data[hit[1],x]) # Save expression value
      }
    }
  }
  index <- index + 10 # Increase index to next set of random samples for next replicate/sample
}

# Find medians of random sample data, excluding placeholders
medians_random <- matrix(data = NA, ncol = 1, nrow = 150)
for (i in 1:150){
  # Save medians of extracted data for GIMME integration
  medians_random[i, 1] <- median(as.numeric(random_samples[grep('-1', random_samples[,i], invert = T), i]))
}

##################################################################
# Save Data

# Set working directory to where files should be saved
setwd("C:/Users/Ana/Documents/MATLAB/Capstone/Relative Proteomics")

write.csv(matching_expression, file = "Relative Proteomics For GIMME.csv", row.names = F)
write.csv(medians, file = "Medians For Relative Proteomics.csv", row.names = F)
write.csv(random_samples, file = "Random Samples For Relative Proteomics.csv", row.names = F)
write.csv(medians_random, file = "Random Sample Medians For Relative Proteomics.csv", row.names = F)


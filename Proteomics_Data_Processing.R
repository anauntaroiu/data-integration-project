
setwd("C:/Users/Ana/Documents/Capstone/R Code/Proteomics") # set working directory to were the file can be found

# Load data

Cytosol_Data <- as.matrix(read.csv("Relative Proteomics - Cytosol.csv"))
Membrane_Data <- as.matrix(read.csv("Relative Proteomics - Membrane.csv"))
Absolute_Data <- as.matrix(read.csv("Absolute Proteomics.csv"))

########################################################################
# Add proteins measured in membrane relative proteomics data if not included
# in the cytosol relative proteomics data

# Generate variable to hold additional proteins
Additional_Proteins <- matrix(data = NA, nrow = 1000, ncol = 16)
index <- 1 # Tracks the indices in Additional_Proteins matrix
counter <- 0 # Tracks whether protein is present in membrane and cytosol data

for (i in 1:dim(Membrane_Data)[1]){ # Loop through each protein measured in membrane data
  
  for (j in 1:dim(Cytosol_Data)[1]){ # Loop through each protein measured in cytosol data
    
  if (agrepl(Membrane_Data[i,1], Cytosol_Data[j,1], max.distance = 0.0)){ # Check if proteins in membrane data are in cytosol data
    counter  <- 1 # If present, set counter to one
  }
    
  if (j == dim(Cytosol_Data)[1] && counter == 0){ # If you reach the end and no match, add protein
    Additional_Proteins[index,1:16] <- Membrane_Data[i,1:16] # Add protein to Additional_Proteins
    index <- index + 1 # Increase index by one
  }
    
    if (j == dim(Cytosol_Data)[1] && counter == 1){ # If you reach the end and there is a match
      counter <- 0 # Reset counter to zero for the next loop
    }
  }
}

# Combine cytosol data with additional proteins
Data <- rbind(Cytosol_Data, Additional_Proteins[1:(index-1),])

################################################################################
# Add proteins measured in absolute proteomics data if not included
# in the relative proteomics data
#
# Use the relative data as a default (since it has higher protein coverage)
# and add additional proteins measured in absolute data

# Generate a variable to store additional proteins
Additional_Proteins <- matrix(data = NA, nrow = 1000, ncol = 16)
index <- 1 # Tracks the indices of the Additional_Proteins matrix

for (i in 2:dim(Absolute_Data)[1]){ # Loop through protein measured in absolute data
  
  for (j in 1:dim(Data)[1]){ # Loop through each protein in the relative data (cytosol + membrane data)
    
      if (agrepl(Absolute_Data[i,1], Data[j,1], max.distance = 0.0)){ # If protein is found in both absolute and relative data
      break # Exit loop
    }
  
  if (j == dim(Data)[1]){ # If you reach the end of the relative data (meaning no match is found)
    Additional_Proteins[index,1:16] <- Absolute_Data[i,1:16] # Add the protein to the additional proteins
    index <- index + 1 # Increase index by one
  }
 }
}

# Save the absolute data from proteins not found in relative data
Absolute_Proteins <- Additional_Proteins[1:(index-1),]

############################################################################
# Save Data

write.csv(Data, file = "Combined_Relative_Proteomics.csv")
write.csv(Absolute_Proteins, file = "Combined_Absolute_Proteomics.csv")






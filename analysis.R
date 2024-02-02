# Created by Mackenzie Mekler 
# University of Florida Center for Neurogenetics - Swanson lab
# 
# Analyze and group neural stem cell data for neurosphere assay

# import required libraries 
library(tidyverse)

# set directory with the data
workdir <- "//192.168.68.57/SwansonLab/Neural Stem Cell/Data/Neurosphere Assays/DKO"

### Make the combined_data.csv file #############
# make vector of folder names that contain the data
folders <- list.dirs("../7dpp")[-1]

# dataframe of all neurosphere init before for loop
combined_data <- data.frame()

# nested for loops to go through all .csv files in all of the folders
for(folder in folders){
  for(img in list.files(folder)){
    # isolating the .csv files 
    if(str_sub(img, -4, -1) == ".csv"){
      table <- read.csv(paste(folder, "/", img, sep=''))
      # table of all spheres 
      sphere <- subset(table, Area != 0)
      # tally number of spheres
      num_spheres <- length(table$X.1) - length(sphere$X.1)
      # calculate average area
      average_area <- sum(table$Area) / num_spheres
      # save to the combined_data table
      results <- c(folder, img, num_spheres, average_area, sum(table$Area))
      combined_data <- rbind(combined_data, results)
    }
  }
}

# write to a file named combined_data.csv
colnames(combined_data) <- c("folder", "image", "num_neurospheres", "ave_area", "area")
write.csv(combined_data, "../combined_data.csv")



### Plot by cell line over time ############

# initialize the passagelines dataframe which will record data for each line by passage date
passagelines <- data.frame()

# start of a while loop, i will increment in the loop to access data and end the loop
# if the folder has not 3 .csv files for each, changing the number added to i in each instance should allow this to still work
i <- 1
while(i <= nrow(combined_data)){
  temp <- combined_data[i:(i+2),] # in the folder, three .csv files are stored back to back for one line this accesses all of them
  passagenum <- str_sub(temp$folder, -2, -1) # this gets the passage number from the folder name
  # this if-else ensures that all three rows we took are from the same passage date and breaks the loop if they aren't 
  # this serves as a safety net 
  if(any(passagenum == passagenum[1])){
    passagenum <- passagenum[1]
  }  
  else{
    print('Error: Passage Dates are not the same')
    print(passage[1])
    print(passage[2])
    print(passage[3])
    break
  }
  # making the data that will be stored in the dataframe 
  total_neurosphere <- sum(as.numeric(temp$num_neurospheres))
  total_area <- sum(as.numeric(temp$area))
  ave <- total_area / total_neurosphere
  # vectorizing the data and adding to the table
  temp_vec <- c(passagenum, temp[1,]$image, total_neurosphere, total_area, ave)
  passagelines <- rbind(passagelines, temp_vec)
  # increment i accordingly 
  i <- i + 3
}
colnames(passagelines) <- c("passage", "cell_line", "neurosphere_num", "area", "ave_area")



### By genotype  ##########################
# establish the genotypes 
dko <- c("2-2", "10-7")
wt <- c("5-8", "7-3", "12-10")
ko2 <- c("10-4", "12-5")

# initialize dataframe so that it can be used during the loops 
genotype <- data.frame(matrix(data=c('init', 'init', 'init', 'init', 'init'),nrow=1, ncol=5))
colnames(genotype) <- c("passage", 'type', "neurosphere_num", "area", "ave_area")

# reset i to increment during loop
i <- 1
while(i <= length(dko)){
  line <- dko[i] # current line in loop is taken from the dko vector and position i
  target_rows <- passagelines[(grepl(dko[i], passagelines$cell_line)),] # take any row from passagelines that has current cell line
  for(j in 1:nrow(target_rows)){ # iterate through rows in target_rows with j
    row <- target_rows[j,] # current row 
    if(!row$passage[1] %in% genotype$passage){ # add a new row if this passage date is not already in genotype dataframe
      vector <- c(row$passage, "DKO", row$neurosphere_num, row$area, as.numeric(row$area) / as.numeric(row$neurosphere_num)) 
      genotype <- rbind(genotype, vector)
    }
    else if (!any(genotype[genotype$passage == row$passage[1],]$type =='DKO')){ # if the passage date is in but not for DKO add a new row
      vector <- c(row$passage, "DKO", row$neurosphere_num, row$area, as.numeric(row$area) / as.numeric(row$neurosphere_num)) 
      genotype <- rbind(genotype, vector)
    }
    else{ # otherwise data needs to be added to the current row for DKO at that passage date
      passage_vector <- genotype$passage == row$passage
      type_vector <- genotype$type == "DKO"
      selected <- genotype[(passage_vector & type_vector),]
      
      neurosphere_original <- as.numeric(selected$neurosphere_num)
      area_original <- as.numeric(selected$area)
      
      new_neurosphere <- neurosphere_original + as.numeric(row$neurosphere_num)
      new_area <- area_original + as.numeric(row$area)
      # replace old values with new ones locating cells where passage number and genotype are what we expect
      genotype[(passage_vector & type_vector), "neurosphere_num"] <- new_neurosphere
      genotype[(passage_vector & type_vector), "area"] <- new_area
      genotype[(passage_vector & type_vector), "ave_area"] <- new_area / new_neurosphere
    }
  }
  # increase i by one
  i <- i + 1
}

# for wt and ko2 they are the same as above for dko 
i <- 1
while(i <= length(wt)){
  line <- wt[i]
  target_rows <- passagelines[(grepl(wt[i], passagelines$cell_line)),]
  for(j in 1:nrow(target_rows)){
    row <- target_rows[j,]
    if(!row$passage[1] %in% genotype$passage){
      vector <- c(row$passage, "WT", row$neurosphere_num, row$area, as.numeric(row$area) / as.numeric(row$neurosphere_num)) 
      genotype <- rbind(genotype, vector)
    }
    else if (!(any(genotype[genotype$passage == row$passage[1],]$type =='WT'))){
      vector <- c(row$passage, "WT", row$neurosphere_num, row$area, as.numeric(row$area) / as.numeric(row$neurosphere_num)) 
      genotype <- rbind(genotype, vector)
    }
    else{
      passage_vector <- genotype$passage == row$passage
      type_vector <- genotype$type == "WT"
      selected <- genotype[(passage_vector & type_vector),]
      
      neurosphere_original <- as.numeric(selected$neurosphere_num)
      area_original <- as.numeric(selected$area)
      
      new_neurosphere <- neurosphere_original + as.numeric(row$neurosphere_num)
      new_area <- area_original + as.numeric(row$area)
      
      genotype[(passage_vector & type_vector), "neurosphere_num"] <- new_neurosphere
      genotype[(passage_vector & type_vector), "area"] <- new_area
      genotype[(passage_vector & type_vector), "ave_area"] <- new_area / new_neurosphere
    }
  }
  i <- i + 1
}

i <- 1
while(i <= length(ko2)){
  line <- ko2[i]
  target_rows <- passagelines[(grepl(ko2[i], passagelines$cell_line)),]
  for(j in 1:nrow(target_rows)){
    row <- target_rows[j,]
    if(!row$passage[1] %in% genotype$passage){
      vector <- c(row$passage, "2KO", row$neurosphere_num, row$area, as.numeric(row$area) / as.numeric(row$neurosphere_num)) 
      genotype <- rbind(genotype, vector)
    }
    else if (!any(genotype[genotype$passage == row$passage[1],]$type =='2KO')){
      vector <- c(row$passage, "2KO", row$neurosphere_num, row$area, as.numeric(row$area) / as.numeric(row$neurosphere_num)) 
      genotype <- rbind(genotype, vector)
    }
    else{
      passage_vector <- genotype$passage == row$passage
      type_vector <- genotype$type == "2KO"
      selected <- genotype[(passage_vector & type_vector),]
      
      neurosphere_original <- as.numeric(selected$neurosphere_num)
      area_original <- as.numeric(selected$area)
      
      new_neurosphere <- neurosphere_original + as.numeric(row$neurosphere_num)
      new_area <- area_original + as.numeric(row$area)
      
      genotype[(passage_vector & type_vector), "neurosphere_num"] <- new_neurosphere
      genotype[(passage_vector & type_vector), "area"] <- new_area
      genotype[(passage_vector & type_vector), "ave_area"] <- new_area / new_neurosphere
    }
  }
  i <- i + 1
}


# write genotype dataframe to MBNL12_neurosphere.csv after removing the init row
genotype = genotype[genotype$type != 'init',]
write.csv(genotype, "../MBNL12_neurosphere.csv")

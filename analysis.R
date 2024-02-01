library(tidyverse)

# this is the directory with all of the folders
workdir <- "//192.168.68.57/SwansonLab/Neural Stem Cell/Data/Neurosphere Assays/DKO"

# getwd()

# this is a vector of the folders with all of the data
folders <- list.dirs("../7dpp")[-1]

# dataframe init before for loop
final <- data.frame()

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
      # save to a final table
      results <- c(folder, img, num_spheres, average_area, sum(table$Area))
      final <- rbind(final, results)
    }
  }
}

# write to a final.csv file 
colnames(final) <- c("folder", "image", "num_neurospheres", "ave_area", "area")
write.csv(final, "../final.csv")
final

# Cell lines
final[1]

### Plot by cell line over time ############
line10_4 <- final[grepl("10-4", final$image),]
p4104 <- line10_4[1:3,]
p5104 <- line10_4[4:6,]

passagelines <- data.frame()
i <- 1
while(i <= nrow(final)){
  temp <- final[i:(i+2),]
  passagenum <- str_sub(temp$folder, -2, -1)
  if(passagenum[1] == passagenum[2] && passagenum[2] == passagenum [3]){
    passagenum <- passagenum[1]
  }
  else{
    print('Error: Passage Dates are not the same')
    print(passage[1])
    print(passage[2])
    print(passage[3])
    break
  }
  total_neurosphere <- sum(as.numeric(temp$num_neurospheres))
  total_area <- sum(as.numeric(temp$area))
  ave <- total_area / total_neurosphere
  temp_vec <- c(passagenum, temp[1,]$image, total_neurosphere, total_area, ave)
  passagelines <- rbind(passagelines, temp_vec)
  
  i <- i + 3
}
colnames(passagelines) <- c("passage", "cell_line", "neurosphere_num", "area", "ave_area")



### Separate into genotypes ##########################
dko <- c("2-2", "10-7")
wt <- c("5-8", "7-3", "12-10")
ko2 <- c("10-4", "12-5")

i <- 1
genotype <- data.frame(matrix(data=c('init', 'init', 'init', 'init', 'init'),nrow=1, ncol=5))
colnames(genotype) <- c("passage", 'type', "neurosphere_num", "area", "ave_area")
while(i <= length(dko)){
  line <- dko[i]
  target_rows <- passagelines[(grepl(dko[i], passagelines$cell_line)),]
  # need to determine passage number and type 
  # add a new row for passage number and type if it does not already exist and add to existing one if it does
  for(num in 1:nrow(target_rows)){
    print(num)
    row <- target_rows[num,]
    if(!row$passage[1] %in% genotype$passage){
      vector <- c(row$passage, "DKO", row$neurosphere_num, row$area, as.numeric(row$area) / as.numeric(row$neurosphere_num)) 
      genotype <- rbind(genotype, vector)
    }
    else if (!any(genotype[genotype$passage == row$passage[1],]$type =='DKO')){
      vector <- c(row$passage, "DKO", row$neurosphere_num, row$area, as.numeric(row$area) / as.numeric(row$neurosphere_num)) 
      genotype <- rbind(genotype, vector)
    }
    else{
      passage_vector <- genotype$passage == row$passage
      type_vector <- genotype$type == "DKO"
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
while(i <= length(wt)){
  line <- wt[i]
  target_rows <- passagelines[(grepl(wt[i], passagelines$cell_line)),]
  # need to determine passage number and type 
  # add a new row for passage number and type if it does not already exist and add to existing one if it does
  for(num in 1:nrow(target_rows)){
    print(num)
    row <- target_rows[num,]
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
  # need to determine passage number and type 
  # add a new row for passage number and type if it does not already exist and add to existing one if it does
  for(num in 1:nrow(target_rows)){
    print(num)
    row <- target_rows[num,]
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

write.csv(genotype, "../MBNL12_neurosphere.csv")

#header file for shared functions, 11/12/2025
# Alicia Melotik

library(collections)
library(dplyr)
library(readxl)
library(googledrive) #https://googledrive.tidyverse.org/
library(googlesheets4)
library(tidyverse)

read_into_dataframe <- function(raw_data) {
  all_data <- raw_data
  
  #create column names for data frame
  for (i in 3:7) {
    colnames(all_data)[i] = raw_data[2, i]
  }
  for (i in 8:length(all_data[1,])) {
    colnames(all_data)[i] = raw_data[1, i]
    #make all areas numeric too, not sure why they aren't automatically numbers
  }
  
  #remove first 2 rows of data frame (were just set as column names)
  all_data <- data.frame(lapply(all_data, function(x) tail(x, -2)))
  
  #make the row names the trial names
  row.names(all_data) <- all_data$Name
  all_data[, 1:3] <- NULL
  
  #force blank levels to 0, and make all levels numeric instead of char
  all_data <- all_data %>% mutate(Level = if_else(Level == "", 0, as.numeric(Level)))
  
  
  return (all_data)
}

## ---------------------------------------------------------
## processing_header.R 
##
## Purpose: clean input and create shared, global variables
## Author: Alicia Melotik
## Date Created: 11/12/2025
## Date Modified: 2/4/2026
## ---------------------------------------------------------

library(openxlsx)    #https://www.rdocumentation.org/packages/openxlsx/versions/4.2.8.1
library(googledrive) #https://googledrive.tidyverse.org/
library(tidyverse)   #https://github.com/tidyverse/tidyverse

read_into_dataframe <- function(raw_data) {
  ### WARNING: ASSUMING SAME FORMAT FOR ALL RAW DATA ###
  all_data <- raw_data
  
  #create column names for data frame
  for (i in 2:6) {
    if (!is.na(raw_data[1, i]))
      colnames(all_data)[i] = raw_data[1, i]
  }

  #remove first row of data frame (necessary ones were just set as column names)
  all_data <- data.frame(lapply(all_data, function(x) tail(x, -1)))
  
  #make the row names the trial names and remove data irrelevant to calculations
  rownames(all_data) <- all_data$Name
  #remove columns before "level" column and one after
  while (colnames(all_data)[1] != "Level") {
    all_data[, 1] <- NULL
  }
  all_data[, 2] <- NULL
  
  #force blank levels to 0, and make all levels numeric instead of char
  all_data <- all_data %>% mutate(Level = if_else(is.na(Level), 0, if_else(Level == "", 0, as.numeric(Level)) ))
  
  for (i in 2:length(colnames(all_data))) {
    #make all areas numeric too, not sure why they aren't automatically numbers
    all_data[i] <- lapply(all_data[i], as.numeric)

    #also, change the automatic X in front of analyte names to A so they are sorted correctly
    if (startsWith(colnames(all_data)[i], "X")) {
      name <- colnames(all_data)[i]
      colnames(all_data)[i] <- paste0("a", substr(name, 2, nchar(name)))
    }
  }
  
  all_data <- as.data.frame(all_data)
  
  return (all_data)
}

get_shared_vars <- function(all_data) {
  #function to generate vars shared across most functions, make available for global use
  all_col_names <<- colnames(all_data)
  
  #start at 2 in sequence to skip "Level" column
  #since analytes and istds alternate, select every other column for each list
  analyte_cols <<- all_col_names[seq(2, length(all_col_names) - 1, 2)]
  istd_cols <<- all_col_names[seq(3, length(all_col_names), 2)]
  
  #mapping of which analytes correspond to which istds, 
  # since analytes are in the same order as their corresponding istd
  mapping <<- tibble(
    Analyte = analyte_cols,
    ISTD = istd_cols
  )
  mapping <<- arrange(mapping, Analyte)
}


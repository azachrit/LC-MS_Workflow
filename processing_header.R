## ---------------------------------------------------------
## processing_header.R 
##
## Purpose: clean input and create shared, global variables
## Author: Alicia Melotik
## Date Created: 11/12/2025
## Date Modified: 8/31/2026
## ---------------------------------------------------------

library(openxlsx)    #https://www.rdocumentation.org/packages/openxlsx/versions/4.2.8.1
library(googledrive) #https://googledrive.tidyverse.org/
library(tidyverse)   #https://github.com/tidyverse/tidyverse

# finds most recent file (or by name) from folder given as arg in shared Google Drive
download_file <- function(folder, type) {
  if (!(exists("sd_meta"))) {
    sd_meta <- shared_drive_get("SWEL Lab")
  }
  
  name <- NA
  while (is.na(name)) {
    name <- readline(prompt = paste0("Enter the ", type, " name or press enter to use most recent one: "))
    
    # check if name is valid/file exists
    if (name != "") {
      if (!(endsWith(name, ".xlsx"))) {
        name <- paste0(name, ".xlsx")
      }
      
      file <- drive_get(path = name, shared_drive = sd_meta)
      if (count(file) > 1) {
        print("ERROR: Mutliple files with that name found. Please rename the file or 
              specify which folder the file is in. For example, enter Proccessed Data/new_data.xlsx")
        name <- NA
      } else if (count(file) == 0) {
        paste0("ERROR: No file called ", name, " could be found.")
        name <- NA
      }
      id <- unlist(file[[1, "id"]])
    }
    else {
      # get most recently created file in the proper folder
      if (!(exists("LCMS_files"))) {
        LCMS_files <- drive_find(pattern = "LC-MS/MS", shared_drive = "SWEL Lab", type = "folder") %>% filter(name == "LC-MS/MS")
      }
      
      #files <- drive_ls(path = paste0("SWEL Lab/4. Instrumentation/LC-MS/MS/", folder), 
      #                 type = ".xlsx", order_by = "createdTime desc")
      
      files <- drive_ls(path = as_id(LCMS_files[["id"]]), pattern = folder)
      files <- drive_ls(path = as_id(files), orderBy = "createdTime desc")
      
      id <- unlist(files[[1, "id"]])
      name <- unlist(files[1, "name"])
      
      #check if there are any files in the folder
      if ((nrow(files) == 0) || (identical(id, character(0)) )) {
        paste0("ERROR: Please enter your data into an Excel File in the ", folder, " folder (& make sure the folder exists)")
        return()
      }
    }
  }
  
  # Create a temp file path to download to & download it
  temp_file <- tempfile(fileext = ".xlsx")
  drive_download(as_id(id), path = temp_file, overwrite = TRUE)
  return (list(temp_file, id))
}

### WARNING: ASSUMING SAME FORMAT FOR ALL RAW DATA ###
read_into_dataframe <- function(raw_data) {
  all_data <- raw_data
  colnames(all_data)[1:4] <- c("Name", "Data File", "Type", "Level")

  #remove first row of data frame (necessary ones were just set as column names)
  all_data <- data.frame(lapply(all_data, function(x) tail(x, -1)))
  
  #make the row names the trial names and remove data irrelevant to calculations
  rownames(all_data) <- all_data$Name
  #remove columns before "level" column and one after (acq date/time)
  while (colnames(all_data)[1] != "Level") {
    all_data[, 1] <- NULL
  }
  all_data[, 2] <- NULL
  
  #force empty cells to 0, and make all data numeric instead of char
  all_data <- all_data %>%
    rename_with(~ str_replace(., "^X", ""), starts_with("X")) %>%
    mutate(across(everything(), ~ {
      x <- na_if(.x, "")      # Blanks to NA
      x <- as.numeric(x)      # Force numeric
      coalesce(x, 0)          # All NAs to 0
    }))
  
  all_data <- as.data.frame(all_data)
  
  return (all_data)
}

#function to generate vars shared across most functions, make available for global use
get_shared_vars <- function(all_data, sorted=TRUE) {
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
  
  if (sorted) {
    mapping <<- mapping %>% arrange(tolower(Analyte)) 
  }
}


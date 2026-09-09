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

# finds most recent file (or by name) from folders given as arg in shared Google Drive
download_file <- function(outer_folder, inner_folder, type) {
  if (!(exists("sd_meta"))) {
    sd_meta <- shared_drive_get("SWEL Lab")
  }
  
  name <- NA
  while (is.na(name)) {
    name <- readline(prompt = paste0("Enter the ", type, " file name or press enter to use most recent one: "))
    
    # check if name is valid/file exists
    if (name != "") {
      if (!(endsWith(name, ".xlsx"))) { # makes it only work for excel type files (not google sheets)
        name <- paste0(name, ".xlsx")
      }
      
      file <- drive_get(path = name, shared_drive = sd_meta)
      id <- unlist(file[[1, "id"]])
    }
    else {
      # get most recently created file in the proper folder
      if (!exists("LCMS_folder_id")) {
        LCMS_folder_id <- drive_find(pattern = paste0("^", outer_folder, "$"), shared_drive = "SWEL Lab", type = "folder")$id[1]
      }
      
      # (Uses Drive API 'q' parameter to search directly inside LC-MS/MS)
      subfolder <- drive_find(
        q = sprintf("'%s' in parents and name = '%s' and mimeType = 'application/vnd.google-apps.folder'", LCMS_folder_id, inner_folder),
        shared_drive = "SWEL Lab"
      )
      
      # Check if subfolder exists
      if (nrow(subfolder) == 0) {
        stop(paste0("ERROR: Subfolder '", inner_folder, "' does not exist in LC-MS/MS."))
      }
      
      # Get the most recent file inside that subfolder
      files <- drive_ls(path = as_id(subfolder$id[1]), orderBy = "createdTime desc")
      
      # Check if subfolder contains any files
      if (nrow(files) == 0) {
        stop(paste0("ERROR: Please enter your data into an Excel File in the ", inner_folder, " folder."))
      }
      
      id <- files$id[1]
      name <- files$name[1]
    }
  }
  
  # Create a temp file path to download to & download it
  temp_file <- tempfile(fileext = ".xlsx")
  drive_download(as_id(id), path = temp_file, overwrite = TRUE)
  return (list(temp_file, id, name))
}

### WARNING: ASSUMING SAME FORMAT FOR ALL RAW DATA ###
read_into_dataframe <- function(raw_data) {
  all_data <- raw_data
  
  while (all_data[1, 1] == "") {
    all_data[, 1] <- NULL
  }

  i <- 1
  while (all_data[1, i] != "Area") {
    colnames(all_data)[i] <- all_data[1, i]
    i <- i + 1
  }
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


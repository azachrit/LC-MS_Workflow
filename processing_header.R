## ---------------------------------------------------------
## processing_header.R 
##
## Purpose: clean input and create shared, global variables
## Author: Alicia Melotik
## Date Created: 11/12/2025
## Date Modified: 10/28/2026
## ---------------------------------------------------------

library(openxlsx)    #https://www.rdocumentation.org/packages/openxlsx/versions/4.2.8.1
library(googledrive) #https://googledrive.tidyverse.org/
library(tidyverse)   #https://github.com/tidyverse/tidyverse

#' Standardize Analyte Names to Match Headers Across Data Frames
#'
#' Converts chemical names with hyphens/colons (e.g. "3-3 FTCA" or "3:3 FTCA")
#' and R auto-formatted column names (e.g. "X3.3.FTCA.Results") into a matching,
#' consistent key string.
normalize_analyte_names <- function(cols) {
  cols %>%
    # Remove leading 'X' assigned by R name checking on numeric starts
    str_replace("^X", "") %>%
    # Convert colons, hyphens, or multiple consecutive dots to standard single dots
    str_replace_all("[:\\-\\s]+", ".") %>%
    # Collapse duplicate dots created during conversion (e.g., "3..3" -> "3.3")
    str_replace_all("\\.{2,}", ".") %>%
    # Trim trailing/leading dots
    str_replace("^\\.|\\.$", "")
}

#' Download Most Recent or Specified File from Google Drive
#'
#' @param outer_folder Name of the parent folder in Google Drive
#' @param inner_folder Name of the target subfolder in Google Drive
#' @param type Descriptive label for the prompt string
#'
#' @return A named list containing `temp_file`, `id`, and `name`
download_file <- function(outer_folder, inner_folder, type) {
  if (!(exists("sd_meta"))) {
    sd_meta <- shared_drive_get("SWEL Lab")
  }
  
  input_name <- NA
  while (is.na(input_name)) {
    input_name <- readline(prompt = paste0("Enter the ", type, " file name or press enter to use most recent one: "))
    input_name <- trimws(input_name)
    
    # check if name is valid/file exists
    if (input_name != "") {
      if (!(endsWith(tolower(input_name), ".xlsx"))) { # makes it only work for excel type files (not google sheets)
        input_name <- paste0(input_name, ".xlsx")
      }
      
      file_info <- drive_get(path = input_name, shared_drive = sd_meta)
      if (nrow(file_info) == 0) {
        message("File not found in Drive. Please check the name and try again.")
        next
      }
      id <- file_info$id[1]
      name <- file_info$name[1]
    } else {
      # get most recently created file in the proper folder
      if (!exists("LCMS_folder")) {
        LCMS_folder <- drive_find(pattern = paste0("^", outer_folder, "$"), shared_drive = "SWEL Lab", type = "folder")
      }
      
      if (nrow(LCMS_folder) == 0) {
        stop(sprintf("ERROR: Parent folder '%s' not found.", outer_folder))
      }
      
      # (Uses Drive API 'q' parameter to search directly inside LC-MS/MS)
      subfolder <- drive_find(
        q = sprintf("'%s' in parents and name = '%s' and mimeType = 'application/vnd.google-apps.folder'", LCMS_folder$id[1], inner_folder),
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

#' Read and Format Raw Data into Standardized Data Frame
#'
#' @param raw_data Initial data frame imported from raw excel file
#'
#' @return Formatted data frame with cleaned column headers and sample ordering
### WARNING: ASSUMING SAME FORMAT FOR ALL RAW DATA ###
read_into_dataframe <- function(raw_data) {
  all_data <- raw_data
  
  # Drop leading empty columns vectorially
  non_empty_cols <- which(colSums(!is.na(all_data) & all_data != "") > 0)
  if (length(non_empty_cols) > 0) {
    all_data <- all_data[, non_empty_cols[1]:ncol(all_data), drop = FALSE]
  }
  
  # Promote row 1 header metadata up to "Area" marker
  area_col_idx <- match("Area", unlist(all_data[1, ]))
  if (!is.na(area_col_idx) && area_col_idx > 1) {
    colnames(all_data)[1:(area_col_idx - 1)] <- unlist(all_data[1, 1:(area_col_idx - 1)])
  }
  
  #remove first row of data frame (necessary ones were just set as column names)
  all_data <- all_data[-1, , drop = FALSE]
  
  #make the row names the trial names and remove data irrelevant to calculations
  sample_names <- make.unique(as.character(all_data$Name))
  rownames(all_data) <- sample_names
  
  #remove columns before "level" column and the one after (acq date/time)
  level_col_idx <- match("Level", colnames(all_data))
  if (!is.na(level_col_idx)) {
    # Keep Level column (col 1) and skip col 2 (Acquisition Date/Time)
    cols_to_keep <- c(level_col_idx, (level_col_idx + 2):ncol(all_data))
    cols_to_keep <- cols_to_keep[cols_to_keep <= ncol(all_data)]
    all_data <- all_data[, cols_to_keep, drop = FALSE]
  }
  
  #sort rows by putting replicates first, followed by blanks
  #find row indices of 1 NGML samples
  ngml_idx <- grep("1 *.*NGML|1 *.*ppb", rownames(all_data), ignore.case = TRUE)
  blank_idx <- grep("blank", rownames(all_data), ignore.case = TRUE)
  start_idx <- unique(c(ngml_idx, blank_idx))
  remaining_idx <- setdiff(seq_len(nrow(all_data)), start_idx)
  
  all_data <- all_data[c(start_idx, remaining_idx), , drop = FALSE ]
  
  #force empty cells to 0, and make all data numeric instead of char
  all_data <- all_data %>%
    rename_with(normalize_analyte_names) %>%
    mutate(across(everything(), ~ {
      x <- na_if(.x, "")      # Blanks to NA
      x <- as.numeric(x)      # Force numeric
      coalesce(x, 0)          # All NAs to 0
    }))
  
  all_data <- as.data.frame(all_data)
  
  return (all_data)
}

#' Generate Mapping Data Structure for Analytes and Internal Standards
#'
#' @param all_data Data frame processed by `read_into_dataframe`
#' @param sorted Logical flag indicating whether mapping should be sorted alphabetically by analyte name
#'
#' @return Named list containing `analyte_cols`, `istd_cols`, `mapping`, and `all_col_names`
get_shared_vars <- function(all_data, sorted=TRUE) {
  all_col_names <- normalize_analyte_names(colnames(all_data))
  
  #start at 2 in sequence to skip "Level" column
  #since analytes and istds alternate, select every other column for each list
  analyte_cols <- all_col_names[seq(2, length(all_col_names) - 1, 2)]
  istd_cols <- all_col_names[seq(3, length(all_col_names), 2)]
  
  #mapping of which analytes correspond to which istds, 
  # since analytes are in the same order as their corresponding istd
  mapping <- tibble(
    Analyte = analyte_cols,
    ISTD = istd_cols
  )
  
  if (sorted) {
    mapping <- mapping %>% arrange(tolower(Analyte)) 
  }

  return (list(
    all_col_names = all_col_names,
    analyte_cols  = analyte_cols,
    istd_cols     = istd_cols,
    mapping       = mapping
  ))
}


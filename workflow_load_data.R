# LC data file input, 11/3/25
# Alicia Melotik

library(collections)
library(dplyr)
library(readxl)
library(googledrive) #https://googledrive.tidyverse.org/
library(googlesheets4)

#include functions from header file?
#source("my_header.R")

read_into_dataframe <- function(full_path) {
  raw_data <- read.csv(full_path, header=FALSE)
  all_data <- raw_data
  
  #create column names for data frame
  for (i in 3:7) {
    colnames(all_data)[i] = raw_data[2, i]
  }
  for (i in 8:length(raw_data)) {
    if (grepl("^[0-9]", raw_data[1, i])) {
      #if name starts with a number, add an underscore in front
      raw_data[1, i] <- paste0(".", raw_data[1, i])
    }
    col_name <- raw_data[1, i]
    if (col_name %in% colnames(all_data)) {
      col_name <- paste0(col_name, ".", i)
    }
    colnames(all_data)[i] = col_name
    #also make measurements into numbers from strings
    all_data[, col_name] <- c(0, 0, as.numeric(all_data[3:length(all_data[, col_name]), col_name]))
  }
  
  ####### DOES NOT LIKE THE REPEATED COL NAMES #########
  ### fixed that but istd col names still have x13 before names? ###
  
  #remove first 2 rows of data frame (just set as column names)
  all_data <- data.frame(lapply(all_data, function(x) tail(x, -2)))
  
  #remove first 2 columns of data frame
  all_data$V1 <- NULL
  all_data$V2 <- NULL
  
  #make the row names the trial names
  row.names(all_data) <- all_data$Name
  
  #force blank levels to 0
  for (row in 1:length(all_data[, "Level"])){
    all_data[row, "Level"] = as.numeric(all_data[row, "Level"])
    if (is.na(all_data[row, "Level"])) {
      all_data[row, "Level"] = 0
    } 
  }
  
  return (all_data)
}

read_data <- function(full_path) {
  raw_data <- read.csv(full_path, header=FALSE)
  trial_names <- tail(raw_data[, 3], -2)
  
  #extract target concentrations for method val
  target_conc <- tail(raw_data[, 6], -2) %>% as.numeric()
  
  #store ISTD peak areas in list for quick lookup later
  ISTD_areas <- list()
  for (col_index in seq(9, length(raw_data), by = 2) ) {
    name <- raw_data[[col_index]][[1]]
    ##check if ISTD has already been added to list: if not, add it and its area values
    if (grepl("ISTD", name) && !(substr(name, 1, nchar(name) - 2) %in% names(ISTD_areas))) {
        ISTD_areas[[name]] <- raw_data[-(1:2), col_index] %>% as.numeric()
    }
  }
  
  #extract area values for each analyte
  native_areas <- list()
  compare_to <- list()
  for (col_index in seq(8, length(raw_data) - 1, by = 2) ) {
    analyte <- raw_data[[col_index]][[1]]
    if ((analyte == "Name") || (analyte == "")) { #skip first rows, not an analyte
      next
    }
    #extract peak areas for each analyte
    for (row_index in 3:length(raw_data[[col_index]])) { #start at row 3 to skip "area"
      #analyte list has peak areas in order
      native_areas[[analyte]] <- raw_data[-(1:2), col_index] %>% as.numeric()
      
      #add the ISTD comparing to as final element
      #analyte_areas[[analyte]] <- c(analyte_areas[[analyte]], raw_data[[col_index + 1]][[1]])
      compare_to[[analyte]] <- raw_data[[col_index + 1]][[1]]
    }
  }

  return (list(trial_names, ISTD_areas, native_areas, compare_to))
}

QAQC <- function(trial_names, ISTD_areas, native_areas) {
  #find row indices of 1 NGML samples
  ngml_samples <- which(grepl("1 NGML", trial_names, ignore.case = TRUE))
  
  ## THIS COULD DEFINITELY BE DONE BETTER!
  calcs <- function(areas, expected) {
    #calculate avg, stddev peak areas for each compound
    avg_ngml <- vector()
    stddev_ngml <- vector()
    percent_exp <- vector()
    ngml_count <- length(ngml_samples)
  
    for (analyte in names(areas)) {
      ngml_vector <- vector()
      for (index in 1:ngml_count) {
        ngml_vector <- c(ngml_vector, areas[[analyte]][[index]])
      }
      avg_ngml[[analyte]] <- sum(ngml_vector) / ngml_count
      stddev_ngml[[analyte]] <- sd(ngml_vector)
      if (analyte %in% names(expected)){
        percent_exp[[analyte]] <- avg_ngml[[analyte]] / expected[[analyte]] * 100
      }
    }
    return (list(avg_ngml, stddev_ngml, percent_exp))
  }
  
  result <- calcs(native_areas, expected_native)
  #can replace this w saving necessary values, writing to sheet
  avg_native <- result[[1]]
  stddev_native <- result[[2]]
  percent_exp_native <- result[[3]]
  
  result <- calcs(ISTD_areas, expected_ISTD)
  #can replace this w saving necessary values, writing to sheet
  avg_ISTD <- result[[1]]
  stddev_ISTD <- result[[2]]
  percent_ISTD <- result[[3]]
  
  #calculate native/ISTD ratio for avgs
  
  #calculate relative response for each analyte for each ngml run
  
}

import_method_val <- function() {
  ### update function to grab most recent method val ###
  ### prob shouldn't hardcode sheets-make sure it matches how method vals are generated ###
  
  
  # If it's a true shared drive (not just shared with you), use shared_drive argument
  # my_shared_folder <- drive_find(shared_drive = "Your Shared Drive Name", type = "folder")
  #drive_path <- ""
  
  #find most recent method val from "Processed Method Val Files" folder in Google Drive
  #Find by name (if added to My Drive)
  shared_folder <- drive_find(pattern = "Processed Method Val Files", type = "folder")
  
  MV_files <- drive_ls(shared_folder) %>% arrange(desc(name))
  cur_MV <- MV_files[2, ] #CHANGE LATER, CURRENTLY IGNORING TRACKING SHEET
  ss <- cur_MV[["id"]]
  
  # Create a temp file path to download to
  tmp_file <- tempfile(fileext = ".xlsx")
  # Download the file
  drive_download(as_id(ss), path = tmp_file, overwrite = TRUE)
  
  ### currently, slopes are transposed as final row: change this? ###
  slope_sheet <- read_excel(tmp_file, sheet = 10)
  slopes <- as.numeric(slope_sheet[length(slope_sheet[[2]]), ])
  
  #get expected values for QAQC calcs
  MV_native_data <- read_excel(tmp_file, sheet = 4)
  MV_ISTD_data <- read_excel(tmp_file, sheet = 5)
  #unlink tmp file to delete
  unlink(tmp_file)
  
  #find which row the average (expected) values are in
  avg_row <- -1
  for (index in 1:length(MV_native_data[[1]])) {
    if (grepl("AVERAGE", MV_native_data[[1]][index], ignore.case = TRUE)) {
      avg_row <- index
    }
  }
  #check if expected values could be found
  if (avg_row == -1) {
    paste("ERROR: could not find average values in spreadsheet")
    return (FALSE)
  }
  #analyte is row 1, col i: create list with analyte: expected value pairs
  
  # Extract analyte names from first row (excluding first column)
  analyte_names <- as.character(MV_native_data[1, -1])
  # Extract numeric values from row `avg_row` (excluding first column)
  analyte_values <- as.numeric(MV_native_data[avg_row, -1])
  # Create named vector
  expected_native <- setNames(analyte_values, analyte_names)
  
  # Extract analyte names from first row (excluding first column)
  analyte_names <- as.character(MV_ISTD_data[1, -1])
  # Extract numeric values from row `avg_row` (excluding first column)
  analyte_values <- as.numeric(MV_ISTD_data[avg_row, -1])
  # Create named vector
  expected_ISTD <- setNames(analyte_values, analyte_names)
  
  ### CHECK IF THE ANALYTES ARE ALWAYS IN THE SAME ORDER (for slopes)###
  
  return (list(expected_native, expected_ISTD, slopes))
}

main <- function() {
  filename <- readline("File name: ")
  
  ### update how to get filepath, it's hardcoded in for now ###
  ### using goodle drive library (already imported) ###
  
  full_path <- paste0("C:/Users/77ali/OneDrive/Desktop/SWEL R/", filename)
  
  #check if file path is valid
  if (!(file.exists(full_path))) {
    paste("The file", filename, "could not be found in the current directory")
  } else {
    print("File exists! Proceeding to process data...")
    # Extract date from file name
    date_from_name <- as.Date(str_extract(current_MV$name, "^\\d{8}"), "%Y%m%d")
    
    #read data from csv at full_path
    dataframe <- read_into_dataframe(full_path)
    
    all_data <- read_data(full_path)
    trial_names <- all_data[[1]]
    ISTD_areas <- all_data[[2]]
    native_areas <- all_data[[3]]
    compare_to <- all_data[[4]]
    
    #get expected values from relevant (most recent) method val
    method_val_data <- import_method_val()
    
    #
    QAQC(trial_names, ISTD_areas, native_areas)
  }
  
}

main()


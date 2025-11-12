# LC data file input, 11/3/25
# Alicia Melotik

#include functions and libraries from header file
### UPDATE THIS LINE, COULD MAKE IT SOURCE HEADER FROM GITHUB ONCE REPO IS PUBLIC ###
#source("https://raw.githubusercontent.com/azachrit/LC-MS_Workflow/refs/heads/main/processing_header.R")
source("C:\\Users\\77ali\\OneDrive\\Desktop\\SWEL R\\processing_header.R")


QAQC <- function(all_data) {
  #find row indices of 1 NGML samples
  ngml_samples <- which(grepl("1 *NGML", all_data[, "Name"], ignore.case = TRUE))
  
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
  # If it's a true shared drive (not just shared with you), use shared_drive argument
  #find most recent method val from "Processed Method Val Files" folder in Google Drive
  MV_path <- drive_find(pattern = "Processed Method Val Files", shared_drive = "SWEL Lab", type = "folder")
  cur_MV <- unlist(drive_ls(MV_path, orderBy = "createdTime desc")[1, "id"])
  
  # Create a temp file path to download to
  temp_file <- tempfile(fileext = ".xlsx")
  # Download the file
  drive_download(as_id(cur_MV), path = temp_file, overwrite = TRUE)
  
  ####### currently, slopes are transposed as final row: change this? #######
  slope_sheet <- read_excel(temp_file, sheet = 10)
  slopes <- as.numeric(slope_sheet[length(slope_sheet[[2]]), ])
  
  #get expected values for QAQC calcs
  MV_native_data <- read_excel(temp_file, sheet = 4)
  MV_ISTD_data <- read_excel(temp_file, sheet = 5)
  #unlink temp file to delete
  unlink(temp_file)
  
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
  ### ----Authenticate to Google Drive-------- ###
  googledrive::drive_auth()
  
  #using google drive library (already imported), find raw data file
  file_path <- drive_find(pattern = "LC-MS/MS", shared_drive = "SWEL Lab", type="folder")

  ######## currently just using most recently created file in LC-MS/MS folder ########
  cur_files <- drive_ls(file_path, orderBy = "createdTime desc")
  file_name <- unlist(cur_files[1, "name"])
  file_id <- unlist(cur_files[1, "id"])

  # Extract date from file name (most recent file currently doesn't match this format)
  date_from_name <- as.Date(str_extract(file_name, "^\\d{8}"), "%Y%m%d")
  
  #read data from csv at full_path
  temp_file <- tempfile(fileext = ".xlsx")
  drive_download(as_id(file_id), path = temp_file, overwrite = TRUE)
  raw_data <- read_excel(temp_file, sheet = 1) #raw data currently in sheet 1
  all_data <- read_into_dataframe(raw_data)
  unlink(temp_file)
  
  #get expected values from relevant (most recent) method val
  method_val_data <- import_method_val()
  
  ###QAQC(all_data) 
}

main()


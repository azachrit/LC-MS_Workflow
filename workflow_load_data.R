# LC-MS data file input, 11/3/25
# Alicia Melotik

#include functions and libraries from header file
source("https://raw.githubusercontent.com/azachrit/LC-MS_Workflow/refs/heads/main/processing_header.R")


QAQC <- function(all_data, expected_native, expected_ISTD, slopes_df) {
  #find row indices of 1 NGML samples
  cal_data <- which(grepl("1 *NGML", all_data[, "Name"], ignore.case = TRUE))
  
  native_df <- tibble(
    Analyte = analyte_cols
  ) %>%
    rowwise() %>%
    mutate(
      Average = mean(cal_data[[Analyte]], na.rm = TRUE), #(na.rm means remove NA values if any present)
      Std_Dev = sd(cal_data[[Analyte]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100,
      Percent_of_Expected = (Average / expected_native[Analyte]) * 100
    ) %>%
    ungroup()
  
  #repeat with ISTDs
  ISTD_df <- tibble(
    Analyte = analyte_cols
  ) %>%
    rowwise() %>%
    mutate(
      Average = mean(cal_data[[istd_cols[[Analyte]]]], na.rm = TRUE),
      Std_Dev = sd(cal_data[[istd_cols[[Analyte]]]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100,
      Percent_of_Expected = (Average / expected_native[[istd_cols[[Analyte]]]]) * 100
    ) %>%
    ungroup() %>%
    select(-ISTD)
  
  #format tables to have analytes lengthwise
  native_df <- native_df %>% 
    pivot_longer(cols = c(Average, Std_Dev, RSD, Percent_of_Expected), names_to = "Native Peak Areas", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
  
  ISTD_df <- ISTD_df %>% 
    pivot_longer(cols = c(Average, Std_Dev, RSD, Percent_of_Expected), names_to = "ISTD Peak Areas", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)

  
  
  
  ### THESE AREN'T GOING TO WORK BELOW, TRYING TO PUT FULL COLUMN IN ONE SPACE AND DIVIDE LIST BY NUMBER 
  
  #calculate RR for each run
  RR_df <- tibble(
    Analyte = analyte_cols
  ) %>%
    rowwise() %>%
    mutate(
      Relative_Response = cal_data[[Analyte]] / cal_data[[istd_cols[[Analyte]]]]
    ) %>%
    ungroup()  %>% 
    pivot_longer(cols = rownames(cal_data), names_to = "Relative Response", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
  
  #measured concentration = RR / slope
  conc_df <- tibble(
    Analyte = analyte_cols
  ) %>%
    rowwise() %>%
    mutate(
      Measured_Conc = RR_df[[Analyte]] / slopes_df[[Analyte]]
    ) %>%
    ungroup()  %>% 
    pivot_longer(cols = c(rownames(cal_data), Average, Percent_Difference, Accuracy), names_to = "Measured Conc (ng/ml)", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
  
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
  
  #read slopes into dataframe
  sheets <- getSheetNames(temp_file)
  name <- "Slope and R"
  if ("SLOPE and R" %in% sheets)
    name <- "SLOPE and R"
  slope_df <- read.xlsx(temp_file, sheet = name)
  
  if ("Slope" %in% colnames(slope_df)) {
    #remove last two lines (transposed slopes)
    slope_df <- slope_df[1:num_analytes, ]
    slope_df["Slope"] <- data.frame(lapply(slope_df["Slope"], as.numeric))
    ### if trying to let old format work too, must transpose df here ###
  } else {
    rownames(slope_df) <- slope_df[[1]]
    slope_df[1] <- NULL
    slope_df[1, ] <- data.frame(lapply(slope_df[1, ], as.numeric))
  }
    
  #get expected values for QAQC calcs
  name <- "Variations"
  ### if trying to let old format work too, must split btwn native and istd here ###
  var_df <- read.xlsx(temp_file, sheet = name)
  #unlink temp file to delete
  unlink(temp_file)
  
  #separate variation data into native & ISTD std dataframes
  ### lots of hard coding in the sheet names/rows for data here :/
  expected_native <- var_df[1:3, ]
  expected_ISTD <- var_df[5:7, ]

  return (list(expected_native, expected_ISTD, slopes))
}

main <- function() {
  ### ----Authenticate to Google Drive-------- ###
  googledrive::drive_auth()
  
  #get most recent method val raw data
  file_path <- drive_find(pattern = "Processed Data", shared_drive = "SWEL Lab", type="folder")
  files <- drive_ls(file_path, orderBy = "createdTime desc")
  
  #check if there are any files in the folder
  if (nrow(files) == 0) {
    print("ERROR: Please enter your raw data into a new Excel File in the 'Method Validations' folder")
    return()
  } 
  #CURRENTLY GRABBING MOST RECENTLY CREATED FILE IN THE FOLDER (could prompt for filename)
  cur_file <- files[1, ]
  file_id <- unlist(cur_file[["id"]])
  
  #read data from csv at full_path
  temp_file <- tempfile(fileext = ".xlsx")
  drive_download(as_id(file_id), path = temp_file, overwrite = TRUE)
  raw_data <- readWorkbook(temp_file, sheet = 1) #raw data currently in sheet 1
  
  #get expected values from relevant (most recent) method val
  mv_data <- import_method_val()
  expected_native <- mv_data[1]
  expected_ISTD <- mv_data[2]
  slopes_df <- mv_data[3]
  
  ###QAQC(all_data) 
  
  unlink(temp_file)
}

main()


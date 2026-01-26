# LC-MS data file input, 11/3/25
# Updated 1/26/2026
# Alicia Melotik

#include functions and libraries from header file
source("https://raw.githubusercontent.com/azachrit/LC-MS_Workflow/refs/heads/main/processing_header.R")

import_method_val <- function() {
  # If it's a true shared drive (not just shared with you), use shared_drive argument
  #find most recent method val from "Processed Method Val Files" folder in Google Drive
  MV_path <- drive_find(pattern = "Processed Method Val Files", shared_drive = "SWEL Lab", type = "folder")
  files <- drive_ls(MV_path, orderBy = "createdTime desc")
  cur_MV <- unlist(drive_ls(MV_path, orderBy = "createdTime desc")[1, "id"])
  
  # Create a temp file path to download to
  temp_file <- tempfile(fileext = ".xlsx")
  # Download the file
  drive_download(as_id(cur_MV), path = temp_file, overwrite = TRUE)
  
  #read slopes into dataframe
  sheets <- getSheetNames(temp_file)
  name <- "Slope and R"
  if (name %in% sheets) {
    slope_df <- read.xlsx(temp_file, sheet = name)
  } else {
    slope_df <- read.xlsx(temp_file, sheet = sheets[length(sheets)-1])
    slope_df <- slope_df[length(slope_df[1])-2:length(slope_df[1]), ] 
  }
  
  
  if ("Slope" %in% colnames(slope_df)) {
    #remove last two lines (transposed slopes)
    slope_df <- slope_df[1:num_analytes, ]
    slope_df["Slope"] <- data.frame(lapply(slope_df["Slope"], as.numeric))
    ### if trying to let old format work too, must transpose df here ###
  } else {
    #rownames(slope_df) <- slope_df[[1]]
    #slope_df[1] <- NULL
    slope_df[1, 2:length(slope_df[1, ])] <- data.frame(lapply(slope_df[1, 2:length(slope_df[1, ])], as.numeric))
    
    
  }
  
  #get expected values for QAQC calcs
  name <- "Variations"
  if (name %in% sheets) {
    var_df <- read.xlsx(temp_file, sheet = name)
    ### lots of hard coding in the sheet names/rows for data here :/
    expected_native <- var_df[1:3, ]
    expected_ISTD <- var_df[5:7, ]
  } else {
    ### if trying to let old format work too, must split btwn native and istd here ###
    expected_native <- read.xlsx(temp_file, sheet = "NATIVE STD VARIATIONS")
    expected_ISTD <- read.xlsx(temp_file, sheet = "ISTD VARIATIONS")
  }
  
  #unlink temp file to delete
  unlink(temp_file)
  
  return (list(expected_native, expected_ISTD, slope_df))
}

QAQC <- function(all_data, expected_native, expected_ISTD, slope_df) {
  #find row indices of 1 NGML samples
  indices <- grep("1 *NGML", rownames(all_data), ignore.case = TRUE)
  cal_data <- all_data[indices, ]
  tib <- tibble(
    Analyte = analyte_cols
  ) %>% rowwise()
  
  native_df <- tib %>%
    mutate(
      Average = mean(cal_data[[Analyte]], na.rm = TRUE), #(na.rm means remove NA values if any present)
      Std_Dev = sd(cal_data[[Analyte]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100,
      Percent_of_Expected = (Average / as.numeric(expected_native[[Analyte]][1])) * 100
    ) %>%
    ungroup()
  
  #repeat with ISTDs
  ISTD_df <-  tib %>%
    mutate(
      Average = mean(cal_data[[istd_cols[[Analyte]]]], na.rm = TRUE),
      Std_Dev = sd(cal_data[[istd_cols[[Analyte]]]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100,
      Percent_of_Expected = (Average / as.numeric(expected_ISTD[[Analyte]][1])) * 100
    ) %>%
    ungroup()
  
  #format tables to have analytes lengthwise
  native_df <- native_df %>% 
    pivot_longer(cols = c(Average, Std_Dev, RSD, Percent_of_Expected), names_to = "Native Peak Areas", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
  
  ISTD_df <- ISTD_df %>% 
    pivot_longer(cols = c(Average, Std_Dev, RSD, Percent_of_Expected), names_to = "ISTD Peak Areas", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
  
  #calculate RR for each run
  RR_df <- cal_data[, analyte_cols] / cal_data[, istd_cols]

  #calculate concentration for each run using conc = RR/slope
  conc_df <- as.data.frame(mapply('/', RR_df, slope_df[1, 2:length(slope_df[1, ])]))
  conc_df <- cbind(rownames(RR_df), conc_df)
  colnames(conc_df)[1] <- "Measured Conc (ng/ml)"
  
  #df for other included concentration metrics
  conc_metrics <-  tib %>%
    mutate(
      Average = mean(conc_df[[Analyte]], na.rm = TRUE),
      Percent_Difference = (1 - Average) / ((Average + 1) / 2) * 100,
      Accuracy = Average/1 * 100
    ) %>%
    ungroup()  %>% 
    pivot_longer(cols = c(Average, Percent_Difference, Accuracy), names_to = "Measured Conc (ng/ml)", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)

  #combine the calculated conc df with other metrics df
  conc_df <- rbind(conc_df, conc_metrics)
  
  return (list(native_df, ISTD_df, RR_df, conc_df))
}

peak_areas <- function(all_data, native_df, ISTD_df, slope_df, LOD) {
  corrected_areas <- all_data
  # need corrected RR_df and conc_df
  to_ratio <- function(x) as.numeric(x)/100
  native_ratio <- lapply(native_df[4, -1], to_ratio) #the -1 removes the first column (which says percent of expected)
  ISTD_ratio <- lapply(ISTD_df[4, -1], to_ratio)
  
}

main <- function() {
  ### ----Authenticate to Google Drive-------- ###
  googledrive::drive_auth()
  
  #get most recent method val raw data
  file_path <- drive_find(pattern = "Processed Data", shared_drive = "SWEL Lab", type="folder")
  files <- drive_ls(file_path, orderBy = "createdTime desc")
  
  #check if there are any files in the folder
  if (nrow(files) == 0) {
    print("ERROR: Please enter your raw data into an Excel File in the 'Processed Data' folder")
    return()
  } 
  #CURRENTLY GRABBING MOST RECENTLY CREATED FILE IN THE FOLDER (could prompt for filename)
  cur_file <- files[1, ]
  file_id <- unlist(cur_file[["id"]])
  
  #read data from csv at full_path
  temp_file <- tempfile(fileext = ".xlsx")
  drive_download(as_id(file_id), path = temp_file, overwrite = TRUE)
  raw_data <- readWorkbook(temp_file, sheet = 1) #raw data currently in sheet 1
  
  #format data, get needed global variables
  all_data <- read_into_dataframe(raw_data)
  get_shared_vars(all_data)
  
  #get expected values from relevant (most recent) method val
  mv_data <- import_method_val()
  expected_native <- mv_data[1]
  expected_ISTD <- mv_data[2]
  slope_df <- mv_data[3]
  #ALSO NEED LOD and LOB
  
  QAQC_data <- QAQC(all_data, expected_native, expected_ISTD, slope_df)
  native_df <- QAQC_data[1]
  ISTD_df <- QAQC_data[2]
  RR_df <- QAQC_data[3]
  conc_df <- QAQC_data[4]
  
  
  
  unlink(temp_file)
}

main()


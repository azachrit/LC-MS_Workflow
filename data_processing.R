## ---------------------------------------------------------
## data_processing.R 
##
## Purpose: Automate LCMS/MS data processing
## Author: Alicia Melotik
## Date Created: 11/3/2025
## Date Modified: 2/16/2026
## ---------------------------------------------------------

#include functions and libraries from header file
source("https://raw.githubusercontent.com/azachrit/LC-MS_Workflow/refs/heads/main/processing_header.R")

import_method_val <- function() {
  # If it's a true shared drive (not just shared with you), use shared_drive argument
  #find most recent method val from "Processed Method Val Files" folder in Google Drive
  MV_path <- drive_find(pattern = "Processed Method Val Files", shared_drive = "SWEL Lab", type = "folder")
  files <- drive_ls(MV_path, orderBy = "createdTime desc")
  cur_MV <- unlist(drive_ls(MV_path, orderBy = "createdTime desc")[1, "id"])
  cur_MV_name <- files[1, "name"]
  
  # Create a temp file path to download to
  temp_file <- tempfile(fileext = ".xlsx")
  # Download the file
  drive_download(as_id(cur_MV), path = temp_file, overwrite = TRUE)
  
  #read slopes into dataframe
  sheets <- getSheetNames(temp_file)
  name <- "Slope and R"
  if (name %in% sheets) { # NEW FORMAT
    slope_df <- read.xlsx(temp_file, sheet = name)
    slope_df <- slope_df[1, ]
    rownames(slope_df) <- slope_df[[1]]
    slope_df[1] <- NULL
  } 
  
  #get expected values for QAQC calcs
  name <- "Variations"
  if (name %in% sheets) { 
    var_df <- read.xlsx(temp_file, sheet = name)
    ### lots of hard coding in the sheet names/rows for data here :/
    expected_native <- var_df[1, ]
    expected_native[, 1] <- NULL
    expected_native[] <- lapply(expected_native, as.numeric)
    expected_ISTD <- var_df[5, ]
    expected_ISTD[, 1] <- NULL
    expected_ISTD[] <- lapply(expected_ISTD, as.numeric)
  }  
  
  #get LOD, LOB, LOQ
  # DID NOT WORRY ABT OLD FORMAT FOR THIS PART AT ALL
  name <- "LOB LOD"
  if (name %in% sheets) { 
    limits_df <- read.xlsx(temp_file, sheet = name)
    limits_df <- limits_df[-1:-4, ] # remove additional lob data (avg, std) and just keep relevant values
    #convert to numeric
    rownames(limits_df) <- limits_df[, 1]
    limits_df <- limits_df[, -1]
    limits_df[] <- lapply(limits_df, as.numeric)
  }
  
  #unlink temp file to delete
  unlink(temp_file)
  
  return (list(expected_native, expected_ISTD, slope_df, limits_df, cur_MV_name))
}

QAQC <- function(cal_data, expected_native, expected_ISTD, slope_df) {
  native_df <- mapping %>%
    rowwise() %>%
    mutate(
      Average = mean(cal_data[[Analyte]], na.rm = TRUE), #(na.rm means remove NA values if any present)
      Std_Dev = sd(cal_data[[Analyte]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100,
      Percent_of_Expected = (Average / as.numeric(expected_native[[Analyte]])) * 100
    ) %>%
    ungroup() %>%
    select(-ISTD)
  
  #repeat with ISTDs
  ISTD_df <-  mapping %>%
    rowwise() %>%
    mutate(
      Average = mean(cal_data[[ISTD]], na.rm = TRUE),
      Std_Dev = sd(cal_data[[ISTD]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100,
      Percent_of_Expected = (Average / as.numeric(expected_ISTD[[Analyte]])) * 100
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
  
  #calculate RR for each run
  RR_df <- cal_data[, mapping[["Analyte"]]] / cal_data[, mapping[["ISTD"]]]

  #calculate concentration for each run using conc = RR/slope
  conc_df <- as.data.frame(mapply('/', RR_df, slope_df))
  conc_df <- cbind(rownames(RR_df), conc_df)
  colnames(conc_df)[1] <- "Measured Conc (ng/ml)"
  
  #df for other included concentration metrics
  conc_metrics <-  mapping %>%
    rowwise() %>%
    mutate(
      Average = mean(conc_df[[Analyte]], na.rm = TRUE),
      Percent_Difference = (1 - Average) / ((Average + 1) / 2) * 100,
      Accuracy = Average/1 * 100
    ) %>%
    ungroup() %>%
    select(-ISTD)  %>% 
    pivot_longer(cols = c(Average, Percent_Difference, Accuracy), names_to = "Measured Conc (ng/ml)", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)

  #combine the calculated conc df with other metrics df
  conc_df <- rbind(conc_df, conc_metrics)
  
  return (list(native_df, ISTD_df, RR_df, conc_df))
}

peak_areas <- function(all_data, native_df, ISTD_df, slope_df) {
  # need corrected RR_df and conc_df
  to_ratio <- function(x) as.numeric(x)/100
  native_ratio <- lapply(native_df[4, -1], to_ratio) #the -1 removes the first column (which says percent of expected)
  ISTD_ratio <- lapply(ISTD_df[4, -1], to_ratio)
  
  #divide each row by the native or ISTD ratio
  corrected_native <- as.data.frame(mapply('/', all_data[, mapping[["Analyte"]]], native_ratio))
  corrected_ISTD <- as.data.frame(mapply('/', all_data[, mapping[["ISTD"]]], ISTD_ratio))
  
  #calculate new concentrations using corrected areas ratio divided by method val slope
  corrected_RR_df <- corrected_native / corrected_ISTD
  corrected_conc_df <- mapply('/', corrected_RR_df, slope_df)

  rownames(corrected_conc_df) <- rownames(all_data)
  
  return (as.data.frame(corrected_conc_df))
}

blank_subs <- function(corrected_conc_df, LOB) {
  #calc 95(?)% confidence interval for maximal blank concentrations
  if (ncol(corrected_conc_df) != length(LOB)) {
    print("ERROR: LOB has different amount of analytes from the dataframe")
    return
  }
  #force NAs to 0 to not impede calculations
  corrected_conc_df[is.na(corrected_conc_df)] <- 0
  corrected_conc_df <- as.matrix(corrected_conc_df)
  
  # calculate maximal blank (first row)
  blank_vals <- as.numeric(corrected_conc_df[1, , drop = TRUE])
  maximal_blank <- blank_vals + (LOB * 1.895)
  
  #repeat maximal_blank rows to reshape to size of conc_df and perform subtractions
  maximal_blank <- maximal_blank[rep(seq_len(nrow(maximal_blank)), nrow(corrected_conc_df)), ]
  corrected_conc_df <- corrected_conc_df - maximal_blank
  corrected_conc_df[corrected_conc_df < 0] <- 0 #floor at zero
 
  return (as.data.frame(corrected_conc_df))
}

check_limit <- function(conc_df, limit_vec, default) {
  #using LOQ for each analyte, indicate where data is below the threshold
  conc_df[] <- mapply(
    function(conc, limit) ifelse(is.na(conc), default, ifelse(conc < limit, default, as.numeric(conc))),
    conc_df,
    limit_vec,
    SIMPLIFY = TRUE
  )
  return (as.data.frame(conc_df))
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
  #CURRENTLY GRABBING MOST RECENTLY CREATED FILE IN THE FOLDER (could prompt for filename or date)
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
  expected_native <- as.data.frame(mv_data[1])
  expected_ISTD <- as.data.frame(mv_data[2])
  slope_df <- as.data.frame(mv_data[3])
  limits_df <- as.data.frame(mv_data[4])
  mv_name <- mv_data[[5]][["name"]]
  
  #find row indices of 1 NGML samples
  indices <- grep("1 *NGML", rownames(all_data), ignore.case = TRUE)
  cal_data <- all_data[indices, ]
  #remove 1 NGML samples from all_data
  indices <- grep("1 *NGML", rownames(all_data), ignore.case = TRUE, invert = TRUE)
  all_data <- all_data[indices, ]
  
  QAQC_data <- QAQC(cal_data, expected_native, expected_ISTD, slope_df)
  native_df <- as.data.frame(QAQC_data[1])
  ISTD_df <- as.data.frame(QAQC_data[2])
  RR_df <- as.data.frame(QAQC_data[3])
  conc_df <- as.data.frame(QAQC_data[4])
  
  #pull LOD/B/Q vector from limits df
  LOD <- limits_df["LOD", ]
  LOB <- limits_df["LOB", ]
  LOQ <- limits_df["LOQ", ]
  
  #peak areas
  corrected_conc_df <- peak_areas(all_data, native_df, ISTD_df, slope_df)
  ### for exporting? set NAs to strings
  #corrected_conc_df[is.na] <- "< LOD"
  # or
  #corrected_conc_df[] <- mapply(
  #  function(conc) ifelse(is.na(conc), "< LOD", conc),
  #  corrected_conc_df,
  #  SIMPLIFY = TRUE
  #)
  
  #perform blank subtractions x2, first time to standardize, second time to subtract blank vals
  blank_sub_df <- blank_subs(corrected_conc_df, LOB)
  blank_sub_df_2 <- blank_subs(blank_sub_df, LOB)
  #set rownames to be the sample names
  rownames(blank_sub_df) <- rownames(all_data)
  rownames(blank_sub_df_2) <- rownames(all_data)
  
  LOQ_flags <- check_limit(blank_sub_df_2, LOQ, "<LOQ")
  
  
  #### CALCULATIONS DONE, NOW WRITE DATA TO EXCEL SHEET ####
  
  #set up global style formatting for workbook
  wb <- loadWorkbook(temp_file, na.convert = FALSE)
  modifyBaseFont(wb, fontName = "Times New Roman")
  options("openxlsx.numFmt" = "0.000")
  hs1 <- createStyle(halign = "justify", textDecoration = "Bold", fontName = "Calibri", fontColour = "black", fgFill = "#e3e1e5")
  
  #add relevant data to temp file through openxlsx workbook
  addWorksheet(wb, "Condensed Data")
  writeData(wb, "Condensed Data", all_data, rowNames = TRUE, headerStyle = hs1)
  
  addWorksheet(wb, "QA QC")
  writeData(wb, "QA QC", expected_native, rowNames = FALSE, startCol = 2, headerStyle = hs1)
  writeData(wb, "QA QC", expected_ISTD, rowNames = FALSE, colNames = FALSE, startRow = 3, startCol = 2)
  writeData(wb, "QA QC", native_df, rowNames = FALSE, startRow = 5, headerStyle = hs1)
  writeData(wb, "QA QC", ISTD_df, rowNames = FALSE, startRow = nrow(native_df) + 7, headerStyle = hs1)
  writeData(wb, "QA QC", RR_df, rowNames = TRUE, startRow = nrow(native_df) + nrow(ISTD_df) + 9, headerStyle = hs1)
  writeData(wb, "QA QC", conc_df, rowNames = FALSE, startRow = nrow(native_df) + nrow(ISTD_df) + nrow(RR_df) + 11, headerStyle = hs1)
  
  #add table names for clarity at top left of each table
  writeData(wb, "QA QC", paste0("Method val data from: ", mv_name), startCol = 1, startRow = 1)
  writeData(wb, "QA QC", "Expected Native", startRow = 2, startCol = 1)
  writeData(wb, "QA QC", "Expected ISTD", startRow = 3)
  writeData(wb, "QA QC", "Relative Response", startRow = nrow(native_df) + nrow(ISTD_df) + 9)
  
  
  addWorksheet(wb, "Peak Areas")
  writeData(wb, "Peak Areas", corrected_conc_df, rowNames = TRUE, headerStyle = hs1)
  # check against lod and reprint to file
  corrected_conc_df <- check_limit(corrected_conc_df, LOD, "<LOD")
  writeData(wb, "Peak Areas", LOD, rowNames = TRUE, startRow = nrow(corrected_conc_df) + 3, headerStyle = hs1)
  writeData(wb, "Peak Areas", corrected_conc_df, rowNames = TRUE, startRow = nrow(corrected_conc_df) + 6, headerStyle = hs1)
  writeData(wb, "Peak Areas", "Concentration (ng/mL)", startCol = 1, startRow = 1)
  writeData(wb, "Peak Areas", "Concentration w/ LOD flags", startCol = 1, startRow = nrow(corrected_conc_df) + 6)
  
  addWorksheet(wb, "Blank Subtractions")
  writeData(wb, "Blank Subtractions", blank_sub_df, rowNames = TRUE, headerStyle = hs1)
  writeData(wb, "Blank Subtractions", blank_sub_df_2, startRow = nrow(blank_sub_df) + 3, rowNames = TRUE, headerStyle = hs1)
  writeData(wb, "Blank Subtractions", "Round 1", startCol = 1, startRow = 1)
  writeData(wb, "Blank Subtractions", "Round 2", startCol = , startRow = nrow(blank_sub_df) + 3)
  
  addWorksheet(wb, "Limit of Quantification")
  writeData(wb, "Limit of Quantification", LOQ_flags, rowNames = TRUE, headerStyle = hs1)
  writeData(wb, "Limit of Quantification", "Concentration (ng/mL)", startCol = 1, startRow = 1)
  
  #addWorksheet(wb, "Final Conversions")
  #writeData(wb, "Final Conversions", slope_df, rowNames = TRUE, headerStyle = hs1)
  
  #adjust all column widths so they're readable
  #first sheet is always raw data (for now)
  first <- 1
  for (sheet_name in names(wb)) {
    num_cols <- ncol(readWorkbook(wb, sheet = sheet_name, colNames = FALSE, rows = 1)) + 1
    #don't adjust spacing for first sheet with raw input, doesn't work well
    if (first == 1) {
      names(wb)[[1]] <- "Raw Data"
      first <- 0
    } else {
      setColWidths(wb, sheet=sheet_name, cols = 1:num_cols, widths="auto")
    }
  }
  #special column width & wrapping style for method val file name in QA QC tab
  setColWidths(wb, sheet="QA QC", cols = 1, widths=22)
  wrap_style <- createStyle(wrapText = TRUE,fontName = "Calibri", fontColour = "red", fgFill = "#e3e1e5")
  addStyle(wb, sheet = "QA QC", wrap_style, rows = 1, cols = 1, gridExpand = TRUE)
  
  #using google drive library (already imported), update file that raw data was pulled from
  new_name <- paste0(format(Sys.Date(), format = "%Y-%m-%d"), "_analyzed.xlsx")
  
  #save and upload excel file to shared SWEL drive
  saveWorkbook(wb, temp_file, overwrite = TRUE)
  drive_update(file = cur_file, media = temp_file, name = new_name)
  
  
  unlink(temp_file)
}

main()


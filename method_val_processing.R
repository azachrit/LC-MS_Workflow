# Method Val Draft R Script, 11/3/25
# Alicia Melotik

# NOTE ON LINEARITY:
###   always forcing intercept through origin
###   always linear curve with no weight
###   range usually (always?) 0.05-100 ng/mL or 50-100,000 ng/L

####### percent differences in accuracy tab are off by a little bit???


#include functions and libraries from header file
source("https://raw.githubusercontent.com/azachrit/LC-MS_Workflow/refs/heads/main/processing_header.R")

calc_slopes <- function(all_data) {
  #blanks removed data
  caldata <- all_data %>% filter(Level > 0)
  
  #get the unique levels and set as the row names
  unique_levels <- sort(unique(caldata$Level))
  num_levels <- length(unique_levels)

  native_avg_table <- matrix(nrow = num_levels, ncol = num_analytes, dimnames = list(unique_levels, analyte_cols))
  ISTD_avg_table <- matrix(nrow = num_levels, ncol = num_analytes, dimnames = list(unique_levels, istd_cols))
  N2Iratio_table <- matrix(nrow = num_levels, ncol = num_analytes, dimnames = list(unique_levels, analyte_cols))
  slopes <- matrix(nrow = 2, ncol = num_analytes, dimnames = list(c("Slope", "R2"), analyte_cols))
  
  #loop through each analyte to calculate each slope
  for (analyte in analyte_cols) {
    ISTD <- istd_cols[[analyte]]
    
    #find N2IRatio across each sample and create vector to store values for lm use
    #N2Iratio <- caldata[[analyte]] / caldata[[ISTD]]
    df_summary <- caldata %>%
      group_by(Level) %>%
      summarise(
        native_avg = mean(.data[[analyte]], na.rm = TRUE),
        ISTD_avg = mean(.data[[ISTD]], na.rm = TRUE),
        RR = native_avg / ISTD_avg,
        .groups = "drop"
      )
    
    #store relevant info into matrices and update column name to current analyte/ISTD
    native_avg_table[, analyte] <- df_summary$native_avg
    ISTD_avg_table[, ISTD] <- df_summary$ISTD_avg
    N2Iratio_table[, analyte] <- df_summary$RR
    
    #linear model regression from Alison's draft r script
    cal_reg <- lm(RR ~ 0 + Level, data = df_summary) ## added "0 +" forces line through origin 
    slope <- coef(cal_reg)[["Level"]]
    r2 <- summary(cal_reg)$r.squared
    
    # Save slope and r^2 results in list
    slopes[, analyte] <- c(slope, r2)
  }
  
  return (list(slopes, native_avg_table, ISTD_avg_table, N2Iratio_table, num_levels))
}

##variation calculations same for other workflow script (avg, std dev, RSD) 
#(use header file again)
variation_calcs <- function(all_data) {
  #avg already calculated above in slope func, move to somewhere where it can be accessed here?
  replicate_data <- all_data %>% filter(Level == 1)
  num_samples <- length(replicate_data[, 1])
  
  native_var_table <- tibble(
    Analyte = analyte_cols
  ) %>%
    rowwise() %>%
    mutate(
      Average = mean(replicate_data[[Analyte]], na.rm = TRUE), #(na.rm means remove NA values if any present)
      Std_Dev = sd(replicate_data[[Analyte]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100
    ) %>%
    ungroup()
  
  #repeat with ISTDs
  ISTD_var_table <- tibble(
    Analyte = analyte_cols,
    ISTD = istd_cols
  ) %>%
    rowwise() %>%
    mutate(
      Average = mean(replicate_data[[istd_cols[[Analyte]]]], na.rm = TRUE),
      Std_Dev = sd(replicate_data[[istd_cols[[Analyte]]]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100
    ) %>%
    ungroup() %>%
    select(-ISTD)
  
  #format tables to have analytes lengthwise
  native_var_df <- native_var_table %>% 
    pivot_longer(cols = c(Average, Std_Dev, RSD), names_to = "1 ng/ml Native", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)

  ISTD_var_df <- ISTD_var_table %>% 
    pivot_longer(cols = c(Average, Std_Dev, RSD), names_to = "1 ng/ml ISTD", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
  
  return (list(native_var_df, ISTD_var_df))
} 

conc_calcs <- function(all_data, slopes) {
  #measured conc = RR/slope
  #### currently doing it exactly like 7/22 and 8/22 examples: 
  ####       not using avgs and only using blanks and 50 ng/l samples ####
  replicate_data <- all_data %>% filter(Level == 0 | Level == 0.05)

  #for each replicate sample, calculate RR and use prior calculated slope to get conc
  conc_df <- matrix(nrow = nrow(replicate_data), ncol = num_analytes + 1, 
                       dimnames = list(rownames(replicate_data), c("Level", analyte_cols)))
  
  for (name in rownames(conc_df)) {
    conc_df[[name, "Level"]] <- all_data[[name, "Level"]]
  }
  
  for (idx in seq_len(num_analytes)) {
    analyte <- analyte_cols[idx]
    istd    <- istd_cols[analyte]
    RR      <- replicate_data[[analyte]] / replicate_data[[istd]]
    slope   <- slopes[[1, analyte]]
    conc_df[, analyte] <- RR / slope
  }
  
  return (conc_df)
}

accuracy <- function(slopes, native_avg, ISTD_avg) {
  #measured conc and percent differences
  accuracy_table <- tibble(
    Analyte = analyte_cols
  ) %>%
    rowwise() %>%
    mutate(
      Relative_Response = (native_avg[[1, Analyte]] / ISTD_avg[[1, Analyte]]),
      Measured_Conc = (Relative_Response / slopes[[1, Analyte]]),
      Percent_Difference = (Measured_Conc - 1) / ((Measured_Conc + 1)/2) * 100
    ) %>%
    ungroup()
  
  accuracy_df <- accuracy_table %>% 
    pivot_longer(cols = c(Relative_Response, Measured_Conc, Percent_Difference), names_to = "Accuracy", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
  
  return (accuracy_df)
}

LOB <- function(conc_df) {
  #LoB = avg conc. + (1.645 * std dev of blank replicates)
  replicate_data <- conc_df %>% filter(Level == 0)
  LOB_table <- tibble(
    Analyte = analyte_cols
  ) %>%
    rowwise() %>%
    mutate(
      avg  = mean(replicate_data[[Analyte]], na.rm = TRUE), #(na.rm means remove NA values if any present)
      sd   = sd(replicate_data[[Analyte]],   na.rm = TRUE),
      LOB = avg + (1.645 * sd)
    ) %>%
    ungroup()
  #transpose so analytes as column names
  LOB_df <- LOB_table %>% pivot_longer(cols = c(avg, sd, LOB), names_to = "metric", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)

  return (LOB_df)
}

LOD_LOQ <- function(conc_df, LOB_df) {
  #LoD = LoB + (1.645 * std dev of 50 ng/L replicates)
  replicate_data <- conc_df %>% filter(Level == 0.05)
  LOD_table <- tibble(
    Analyte = analyte_cols
  ) %>%
    rowwise() %>%
    mutate(
      LOB  = LOB_df[[3, Analyte]],
      sd   = sd(replicate_data[[Analyte]],   na.rm = TRUE),
      LOD = LOB + (1.645 * sd),
      LOQ = pmax(LOD * 3, 0.025)
    ) %>%
    ungroup()
  #transpose so analytes as column names
  LOD_df <- LOD_table %>% 
    pivot_longer(cols = c(LOB, sd, LOD, LOQ), names_to = "metric", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)

  return (LOD_df)
}

upload_mv <- function() {
  ### ----Authenticate to Google Drive-------- ###
  googledrive::drive_auth()
  
  #get most recent method val raw data
  file_path <- drive_find(pattern = "Method Validations", shared_drive = "SWEL Lab", type="folder")
  files <- drive_ls(file_path, orderBy = "createdTime desc")

  #check if there are any files in the folder
  if (nrow(files) == 0) {
    print("ERROR: Please enter your raw data into a new Excel File in the 'Method Validations' folder")
    return()
  } 
  #CURRENTLY GRABBING MOST RECENTLY CREATED FILE IN METHOD VAL FOLDER
  cur_file <- files[1, ]
  file_id <- unlist(cur_file[["id"]])
  
  #read data from csv at full_path
  temp_file <- tempfile(fileext = ".xlsx")
  drive_download(as_id(file_id), path = temp_file, overwrite = TRUE)
  raw_data <- readWorkbook(temp_file, sheet = 1) #raw data currently in sheet 1
  
  #perform relevant operations on data using functions
  all_data <- read_into_dataframe(raw_data)
  get_shared_vars(all_data)
  
  result <- calc_slopes(all_data)
  slopes_df <- result[[1]]
  native_avg_table <- result[[2]]
  ISTD_avg_table <- result[[3]]
  N2Iratio_table <- result[[4]]
  num_levels <- result[[5]]

  result <- variation_calcs(all_data)
  native_var_df <- result[[1]]
  ISTD_var_df <- result[[2]]
  
  accuracy_df <- accuracy(slopes_df, native_var_df, ISTD_var_df)
  
  conc_df <- as.data.frame(conc_calcs(all_data, slopes_df))
  
  LOB_df <- LOB(conc_df)
  LOD_LOQ_df <- LOD_LOQ(conc_df, LOB_df)
  
  #add relevant data to temp file through openxlsx workbook
  hs1 <- createStyle(halign = "justify", textDecoration = "Bold", border = "TopBottomLeftRight", fontColour = "black")

  wb <- loadWorkbook(temp_file, na.convert = FALSE)
  addWorksheet(wb, "Condensed Data")
  writeData(wb, "Condensed Data", all_data, rowNames = TRUE, headerStyle = hs1)

  addWorksheet(wb, "Peak Areas")
  writeData(wb, "Peak Areas", native_avg_table, rowNames = TRUE, headerStyle = hs1)
  writeData(wb, "Peak Areas", ISTD_avg_table, rowNames = TRUE, startCol = num_analytes + 3, headerStyle = hs1)
  writeData(wb, "Peak Areas", N2Iratio_table, rowNames = TRUE, startRow = num_levels + 3, headerStyle = hs1)
  #add table names for clarity
  writeData(wb, "Peak Areas", "Avg Native Peak Areas", startCol = 1, startRow = 1)
  writeData(wb, "Peak Areas", "Avg ISTD Peak Areas", startCol = num_analytes + 3, startRow = 1)
  writeData(wb, "Peak Areas", "Relative Reponse", startCol = 1, startRow = num_levels + 3)
  
  addWorksheet(wb, "Variations")
  writeData(wb, "Variations", native_var_df, headerStyle = hs1)
  writeData(wb, "Variations", ISTD_var_df, startRow = length(rownames(native_var_df)) + 4, headerStyle = hs1)
  writeData(wb, "Variations", accuracy_df, startRow = length(rownames(native_var_df)) + length(rownames(ISTD_var_df)) + 7, headerStyle = hs1)
  
  addWorksheet(wb, "Conc Calcs")
  writeData(wb, "Conc Calcs", conc_df, rowNames = TRUE, headerStyle = hs1)
  
  addWorksheet(wb, "LOB LOD")
  writeData(wb, "LOB LOD", LOB_df, headerStyle = hs1)
  writeData(wb, "LOB LOD", LOD_LOQ_df, startRow = length(rownames(LOB_df)) + 4, headerStyle = hs1)
  
  addWorksheet(wb, "Slope and R")
  writeData(wb, "Slope and R", slopes_df, rowNames = TRUE, headerStyle = hs1)
  
  #adjust all column widths so they're readable
  for (sheet_name in names(wb)) {
    num_cols <- ncol(readWorkbook(wb, sheet = sheet_name, colNames = FALSE, rows = 1)) + 1
    #don't adjust spacing for first sheet with raw input, doesn't work well
    if ((sheet_name == "RAW") | (sheet_name == "Sheet1")) {
      next
    }
    setColWidths(wb, sheet=sheet_name, cols = 1:num_cols, widths="auto")
  }

  #using google drive library (already imported), update file that raw data was pulled from
  new_name <- paste0(format(Sys.Date(), format = "%Y-%m-%d"), "-Method_Val.xlsx")
  
  #save and upload excel file to shared SWEL drive
  saveWorkbook(wb, temp_file, overwrite = TRUE)
  drive_update(file = cur_file, media = temp_file, name = new_name)

  #delete temp file now that it has been uploaded
  unlink(temp_file)
}

upload_mv()

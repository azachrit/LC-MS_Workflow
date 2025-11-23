# Method Val Draft R Script, 11/3/25
# Alicia Melotik

#include functions and libraries from header file
source("https://raw.githubusercontent.com/azachrit/LC-MS_Workflow/refs/heads/main/processing_header.R")
library(openxlsx)

calc_slopes <- function(all_data) {
  all_col_names <- colnames(all_data) # indices 6 to length of col names / 2 = number of analytes
  slopes <- list()
  row_num <- 0
  
  #get the unique levels and set as the row names
  unique_levels <- sort(unique(caldata$Level))
  num_levels <- length(unique_levels)
  num_analytes <- floor((length(all_col_names) - 5) / 2)

  native_avg_table <- matrix(nrow = num_levels, ncol = num_analytes, dimnames = list(unique_levels, rep("", num_analytes)))
  rownames(native_avg_table) <- unique_levels
  ISTD_avg_table <- native_avg_table
  N2Iratio_table <- native_avg_table
  
  #blanks removed data
  caldata <- all_data %>% filter(Level > 0)
  
  #loop through each analyte to calculate each slope
  for (i in seq(5, length(all_col_names) - 1, 2)) {
    analyte <- all_col_names[i]
    ISTD <- all_col_names[i+1]
    row_num <- row_num + 1
    print(i)
    print(analyte)
    print(ISTD)
    
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
    
    #put relevant info into matrices and update column name to current analyte/ISTD
    native_avg_table[, row_num] <- df_summary$native_avg
    ISTD_avg_table[, row_num] <- df_summary$ISTD_avg
    N2Iratio_table[, row_num] <- df_summary$RR
    colnames(native_avg_table)[row_num] <- analyte
    colnames(ISTD_avg_table)[row_num] <- ISTD
    colnames(N2Iratio_table)[row_num] <- analyte
    
    #linear model regression from Alison's draft r script
    cal_reg <- lm(RR ~ 0 + Level, data = df_summary) ## added "0 +" forces line through origin 
    slope <- coef(cal_reg)[["Level"]]
    r2 <- summary(cal_reg)$r.squared
    
    # Save slope and r^2 results in list
    slopes[[analyte]] <- list(
      slope = slope,
      r2 = r2
    )
  }
  slopes_df <- bind_rows(slopes, .id = "Analyte") %>% 
    pivot_longer(cols = c(slope, r2), names_to = "metric", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
  slopes_sorted <- slopes[order(names(slopes))]
  
  return (c(slopes_df, native_avg_table, ISTD_avg_table, N2Iratio_table))
}

##variation calculations same for other workflow script (avg, std dev, RSD) 
#(use header file again)
variation_calcs <- function(slopes) {
  #avg already calculated above in slope func, move to somewhere where it can be accessed here?
  replicate_data <- all_data %>% filter(Level == 1)
  
  # Identify analyte / ISTD pairs
  analyte_cols <- all_col_names[seq(6, length(all_col_names) - 1, 2)]
  istd_cols <- all_col_names[seq(7, length(all_col_names), 2)]
  
  variation_table <- tibble(
    Analyte = analyte_cols,
    ISTD    = istd_cols
  ) %>%
    rowwise() %>%
    mutate(
      n_avg  = mean(replicate_data[[Analyte]], na.rm = TRUE), #(na.rm means remove NA values if any present)
      n_sd   = sd(replicate_data[[Analyte]],   na.rm = TRUE),
      n_rsd  = (n_sd / n_avg) * 100,
      i_avg  = mean(replicate_data[[ISTD]],    na.rm = TRUE),
      i_sd   = sd(replicate_data[[ISTD]],      na.rm = TRUE),
      i_rsd  = (i_sd / i_avg) * 100
    ) %>%
    ungroup()
} 

conc_calcs <- function(slopes) {
  #measured conc = RR/slope
  #### currently doing it exactly like 7/22 and 8/22 examples: 
  ####       not using avgs and only using blanks and 50 ng/l samples ####
  replicate_data <- all_data %>% filter(Level == 0 | Level == 0.05)
  
  # Identify analyte / ISTD pairs
  analyte_cols <- all_col_names[seq(6, length(all_col_names) - 1, 2)]
  istd_cols <- all_col_names[seq(7, length(all_col_names), 2)]

  #for each replicate sample, calculate RR and use prior calculated slope to get conc
  rows = rownames(replicate_data)
  row_num = 1
  conc_table <- matrix(nrow = length(rows), ncol = num_analytes, dimnames = list(rows, analyte_cols))
  
  for (idx in seq_along(analyte_cols)) {
    analyte <- analyte_cols[idx]
    istd    <- istd_cols[idx]
    for (j in seq_len(nrow(replicate_data))) {
      conc_table[j, idx] <- (replicate_data[[analyte]][j]/replicate_data[[ISTD]][j]) / slopes[[analyte]]$slope
    }
  }
  
}

linearity <- function(slopes, r2) {
  #do we need this?
}

accuracy <- function() {
  #measured conc and percent differences
}

LOB <- function(blanks) {
  #LoB = avg conc. + (1.645 * std dev of blank replicates)
  replicate_data <- all_data %>% filter(Level == 0)
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
  
}

LOD <- function(LOBs) {
  #LoD = LoB + (1.645 * std dev of 50 ng/L replicates)
  replicate_data <- all_data %>% filter(Level == 0.05)
  LOD_table <- tibble(
    Analyte = analyte_cols
  ) %>%
    rowwise() %>%
    mutate(
      LOB  = LOB_df[[Analyte]][3],
      sd   = sd(replicate_data[[Analyte]],   na.rm = TRUE),
      LOD = LOB + (1.645 * sd)
    ) %>%
    ungroup()
  #transpose so analytes as column names
  LOD_df <- LOD_table %>% 
    pivot_longer(cols = c(LOB, sd, LOD), names_to = "metric", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
}

LOQ <- function(LODs) {
  #LoQ = LoD * 3 (or 0.025 if calculation is below that)
  
}

upload_mv <- function() {
  ### ----Authenticate to Google Drive-------- ###
  googledrive::drive_auth()
  
  #get most recent method val data
  ######## DIFF WAY TO INPUT RAW DATA? CURRENTLY PULLING FROM PROCESSED FILES ########
  file_path <- drive_find(pattern = "Processed Method Val Files", shared_drive = "SWEL Lab", type="folder")
  cur_files <- drive_ls(file_path, orderBy = "createdTime desc")
  file_name <- unlist(cur_files[1, "name"])
  file_id <- unlist(cur_files[1, "id"])
  
  #read data from csv at full_path
  temp_file <- tempfile(fileext = ".xlsx")
  drive_download(as_id(file_id), path = temp_file, overwrite = TRUE)
  raw_data <- read_excel(temp_file, sheet = 1) #raw data currently in sheet 1
  unlink(temp_file)
  
  #perform relevant operations on data
  all_data <- read_into_dataframe(raw_data)
  result <- calc_slopes(all_data)
  slopes_df <- result[1]
  native_avg_table <- result[2]
  ISTD_avg_table <- result[3]
  N2Iratio_table <- result[4]
  
  #using google drive library (already imported), find raw data file
  file_path <- drive_find(pattern = "LC-MS/MS", shared_drive = "SWEL Lab", type="folder")
  file_name <- paste0(format(Sys.Date(), format = "%Y-%m-%d"), "-Method_val.xlsx")
  temp_file <- tempfile(fileext = ".xlsx")
  
  #add relevant data to temp file
  #writexl::write_xlsx(list(Sheet1 = raw_data, Sheet2 = all_data, Sheet3 = slopes_df), path = temp_file)
  wb <- createWorkbook()
  addWorksheet(wb, "Raw Data")
  writeData(wb, "Raw Data", raw_data, colNames = FALSE)
  addWorksheet(wb, "Condensed Data")
  writeDataTable(wb, "Condensed Data", all_data, rowNames = TRUE,  tableStyle = "TableStyleLight9")
  addWorksheet(wb, "Peak Areas")
  hs1 <- createStyle(halign = "justify", textDecoration = "Bold", border = "TopBottomLeftRight", fontColour = "black")
  writeData(wb, "Peak Areas", analyte_avg_table, rowNames = TRUE, headerStyle = hs1)
  writeData(wb, "Peak Areas", ISTD_avg_table, rowNames = TRUE, startCol = num_analytes + 3, headerStyle = hs1)
  writeData(wb, "Peak Areas", N2Iratio_table, rowNames = TRUE, startRow = length(levels) + 3, headerStyle = hs1)
  #adjust column widths so they're readable
  for (sheet_name in names(wb)) {
    num_cols <- ncol(readWorkbook(wb, sheet = sheet_name, colNames = FALSE, rows = 1)) + 1
    setColWidths(wb, sheet=sheet_name, cols = 1:num_cols, widths="auto")
  }
  #save and upload excel file to shared SWEL drive
  saveWorkbook(wb, temp_file, overwrite = TRUE)
  drive_upload(temp_file, file_path, name=file_name)
  #delete temp file now that it has been uploaded
  unlink(temp_file)
}

upload_mv()

## ---------------------------------------------------------
## method_val_processing.R 
##
## Purpose: Automate method validation processing
## Author: Alicia Melotik
## Date Created: 11/3/2025
## Date Modified: 2/16/2026
## ---------------------------------------------------------

# NOTE ON LINEARITY:
###   always forcing intercept through origin
###   always linear curve with no weight
###   range usually (always?) 0.05-100 ng/mL or 50-100,000 ng/L

####### percent differences in accuracy tab are off by a little bit???

#include functions and libraries from header file
source("https://raw.githubusercontent.com/azachrit/LC-MS_Workflow/refs/heads/main/processing_header.R")

slope_calcs <- function(all_data) {
  #blanks removed data
  caldata <- all_data %>% filter(Level > 0)
  
  #data frames for desired peak areas and relative response data
  native_summary <- caldata %>%
    group_by(Level) %>%
    summarise(
      across(all_of(mapping[["Analyte"]]), ~ mean(as.numeric(.x), na.rm = TRUE))
    ) %>%
    ungroup()
  istd_summary <- caldata %>%
    group_by(Level) %>%
    summarise(
      across(all_of(mapping[["ISTD"]]), \(x) mean(x, na.rm = TRUE))
    ) %>%
    ungroup()

  RR_summary <- native_summary / istd_summary
  RR_summary[, 1] <- native_summary[, 1] #fix the first column to contain the levels, not level / level
  
  #set up empty matrix to hold slopes
  slopes <- matrix(nrow = 2, ncol = ncol(native_summary) - 1, dimnames = list(c("Slope", "R2"), mapping[["Analyte"]]))
  
  #loop through each analyte to calculate each slope and create plot
  for (analyte in mapping[["Analyte"]]) {
    #extract current data into a tibble
    df_summary <- tibble(
      Level = RR_summary[["Level"]],
      RR = RR_summary[[analyte]]
    )
    
    #linear model regression from Alison's draft r script
    cal_reg <- lm(RR ~ 0 + Level, data = df_summary) ## added "0 +" forces line through origin 
    
    #Plot Cal Curve data & Regression Line
    ggplot(df_summary, aes(Level, RR))+
          geom_point()+
          geom_smooth(method = lm, formula = y ~ 0 + x)+
          labs(title = paste(analyte, "Regression"))+
          theme_classic()
    
    slope <- coef(cal_reg)[["Level"]]
    r2 <- summary(cal_reg)$r.squared
    
    # Save slope and r^2 results in list
    slopes[, analyte] <- c(slope, r2)
  }
  
  return (list(slopes, native_summary, istd_summary, RR_summary))
}

##variation calculations same for other workflow script (avg, std dev, RSD) 
#(use header file again?)
variation_calcs <- function(all_data) {
  #avg already calculated above in slope func, move to somewhere where it can be accessed here?
  replicate_data <- all_data %>% filter(Level == 1)
  num_samples <- nrow(replicate_data)
  
  native_var_df <- tibble(
    Analyte = mapping[["Analyte"]]
  ) %>%
    rowwise() %>%
    mutate(
      Average = mean(replicate_data[[Analyte]], na.rm = TRUE), #(na.rm means remove NA values if any present)
      Std_Dev = sd(replicate_data[[Analyte]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100
    ) %>%
    ungroup()
  
  #repeat with ISTDs
  ISTD_var_df <- tibble(
    Analyte = mapping[["Analyte"]],
    ISTD = mapping[["ISTD"]]
  ) %>%
    rowwise() %>%
    mutate(
      Average = mean(replicate_data[[ISTD]], na.rm = TRUE),
      Std_Dev = sd(replicate_data[[ISTD]], na.rm = TRUE),
      RSD = (Std_Dev / Average) * 100
    ) %>%
    select(-ISTD) %>%
    ungroup()
  
  #format tables to have analytes lengthwise
  native_var_df <- native_var_df %>% 
    pivot_longer(cols = c(Average, Std_Dev, RSD), names_to = "1 ng/ml Native", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)

  ISTD_var_df <- ISTD_var_df %>% 
    pivot_longer(cols = c(Average, Std_Dev, RSD), names_to = "1 ng/ml ISTD", values_to = "value") %>%
    pivot_wider(names_from = Analyte, values_from = value)
  
  return (list(native_var_df, ISTD_var_df))
} 

conc_calcs <- function(all_data, slopes_df, LOD_conc) {
  #measured conc = RR/slope
  #### currently doing it exactly like 7/22 and 8/22 examples: 
  ####       not using avgs and only using blanks and samples of chosen level (sample_conc)###
  replicate_data <- all_data %>% filter(Level == 0 | Level == LOD_conc)
  
  #get slope row & replicate RRs, then divide each by corresponding slope
  slope <- slopes_df[1, ]
  conc_df <- replicate_data[mapping[["Analyte"]]] / replicate_data[mapping[["ISTD"]]]
  conc_df <- sweep(conc_df, 2, slope, '/') #margin 2 means it goes by columns
  
  #insert the levels back into the dataframe as first column
  conc_df <- cbind(replicate_data[, 1], conc_df)
  colnames(conc_df)[1] <- colnames(replicate_data)[1]
  
  return (conc_df)
}

accuracy_calcs <- function(slopes, native_avg, ISTD_avg) {
  #measured conc and percent differences
  accuracy_table <- tibble(
    Analyte = mapping[["Analyte"]]
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

LOB_calcs <- function(conc_df) {
  #LoB = avg conc. + (1.645 * std dev of blank replicates)
  replicate_data <- conc_df %>% filter(Level == 0)
  LOB_table <- tibble(
    Analyte = mapping[["Analyte"]]
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

LOD_LOQ_calcs <- function(conc_df, LOB_df, LOD_conc) {
  #LoD = LoB + (1.645 * std dev of chosen replicates)
  replicate_data <- conc_df %>% filter(Level == LOD_conc)

  LOD_table <- tibble(
    Analyte = mapping[["Analyte"]]
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

main <- function() {
  ### ----Authenticate to Google Drive-------- ###
  googledrive::drive_auth()
  
  #get most recent method val raw data
  file_path <- drive_find(pattern = "Processed Method Val Files", shared_drive = "SWEL Lab", type="folder")
  files <- drive_ls(file_path, orderBy = "createdTime desc")

  #check if there are any files in the folder
  if (nrow(files) == 0) {
    print("ERROR: Please enter your raw data into a new Excel File in the 'Processed Method Val Files' folder")
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
  
  result <- slope_calcs(all_data)
  slopes_df <- result[[1]]
  native_avg_table <- result[[2]]
  ISTD_avg_table <- result[[3]]
  N2Iratio_table <- result[[4]]

  result <- variation_calcs(all_data)
  native_var_df <- result[[1]]
  ISTD_var_df <- result[[2]]
  
  accuracy_df <- accuracy_calcs(slopes_df, native_var_df, ISTD_var_df)
  
  #ASK FOR CONC TO CALCULATE LOD WITH
  LOD_conc <- NA
  while (is.na(LOD_conc)) {
    LOD_conc <- readline(prompt = "Concentration to use for LOD/LOQ calculations: ")
    LOD_conc <- suppressWarnings(as.numeric(LOD_conc))
    if (is.na(LOD_conc))
      next
    
    if (LOD_conc > 1) {
      LOD_conc <- LOD_conc / 1000.0
    } 
    level_options <- unique(unlist(all_data[["Level"]]))
    if (!(LOD_conc %in% level_options)) {
      print("Error: Please enter a concentration/level present in the data")
      print("Choices (in ng/L): ")
      print(unique(unlist(lapply(level_options, function(x) x * 1000))))
      LOD_conc <- NA
    }
  }
  
  conc_df <- as.data.frame(conc_calcs(all_data, slopes_df, LOD_conc))
  
  LOB_df <- LOB_calcs(conc_df)
  LOD_LOQ_df <- LOD_LOQ_calcs(conc_df, LOB_df, LOD_conc)
  
  #for formatting, store matrices shapes dimensions in variables
  num_analytes <- ncol(native_avg_table)
  num_levels <- nrow(native_avg_table)

  #### CALCULATIONS DONE, NOW WRITE DATA TO EXCEL SHEET ####
  
  #set up global style formatting for workbook
  wb <- loadWorkbook(temp_file, na.convert = FALSE)
  modifyBaseFont(wb, fontName = "Times New Roman")
  options("openxlsx.numFmt" = "0.000")
  hs1 <- createStyle(halign = "justify", textDecoration = "Bold", fontName = "Calibri", fontColour = "black", fgFill = "#e3e1e5")
  
  #add relevant data to temp file through openxlsx workbook
  addWorksheet(wb, "Condensed Data")
  writeData(wb, "Condensed Data", all_data, rowNames = TRUE, headerStyle = hs1)

  addWorksheet(wb, "Peak Areas")
  writeData(wb, "Peak Areas", native_avg_table, rowNames = TRUE, headerStyle = hs1)
  writeData(wb, "Peak Areas", ISTD_avg_table, rowNames = TRUE, startCol = num_analytes + 3, headerStyle = hs1)
  writeData(wb, "Peak Areas", N2Iratio_table, rowNames = TRUE, startRow = num_levels + 3, headerStyle = hs1)
  #add table names for clarity at top left of each table
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

  #using google drive library (already imported), update file that raw data was pulled from
  new_name <- paste0(format(Sys.Date(), format = "%Y-%m-%d"), "_Method_Val.xlsx")
  
  #save and upload excel file to shared SWEL drive
  saveWorkbook(wb, temp_file, overwrite = TRUE)
  drive_update(file = cur_file, media = temp_file, name = new_name)

  #delete temp file now that it has been uploaded
  unlink(temp_file)
  
  #store data in tracking sheet
  tracking_sheet <- drive_find(pattern = "Method Validation Metric Tracking", shared_drive = "SWEL Lab")
  
}

main()

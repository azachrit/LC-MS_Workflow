# Method Val Draft R Script, 11/3/25
# Alicia Melotik

### use same code from other script to load data, combine via header file ###
source("processing_header.R")

calc_slopes <- function(all_data) {
  all_col_names <- colnames(all_data) # indices 6 to length of col names / 2 = number of analytes
  slopes <- list()
  row_num <- 0
  
  #get the unique levels and set as the row names
  unique_levels <- sort(unique(caldata$Level))
  num_levels <- length(unique_levels)
  num_analytes <- floor((length(all_col_names) - 5) / 2)

  analyte_avg_table <- matrix(nrow = num_levels, ncol = num_analytes, dimnames = list(unique_levels, rep("", num_analytes)))
  rownames(analyte_avg_table) <- unique_levels
  ISTD_avg_table <- analyte_avg_table
  N2Iratio_table <- analyte_avg_table
  
  #blanks removed data
  caldata <- all_data %>% filter(Level > 0)
  
  #loop through each analyte to calculate each slope
  for (i in seq(6, length(all_col_names) - 1, 2)) {
    analyte <- all_col_names[i]
    ISTD <- all_col_names[i+1]
    row_num <- row_num + 1
    
    #find N2IRatio across each sample and create vector to store values for lm use
    #N2Iratio <- caldata[[analyte]] / caldata[[ISTD]]
    df_summary <- caldata %>%
      group_by(Level) %>%
      summarise(
        analyte_avg = mean(.data[[analyte]], na.rm = TRUE),
        ISTD_avg = mean(.data[[ISTD]], na.rm = TRUE),
        RR = analyte_avg / ISTD_avg,
        .groups = "drop"
      )
    
    #put relevant info into matrices and update column name to current analyte/ISTD
    analyte_avg_table[, row_num] <- df_summary$analyte_avg
    ISTD_avg_table[, row_num] <- df_summary$ISTD_avg
    N2Iratio_table[, row_num] <- df_summary$RR
    colnames(analyte_avg_table)[row_num] <- analyte
    colnames(ISTD_avg_table)[row_num] <- ISTD
    colnames(N2Iratio_table)[row_num] <- analyte
    
    #linear model regression from draft r script
    cal_reg <- lm(RR ~ 0 + Level, data = df_summary) ## added "0 +" forces line through origin 
    slope <- coef(cal_reg)[["Level"]]
    r2 <- summary(cal_reg)$r.squared
    
    # Save slope and r^2 results in list
    slopes[[analyte]] <- list(
      slope = slope,
      r2 = r2
    )
  }
  slopes_sorted <- slopes[order(names(slopes))]
}

##variation calculations same for other workflow script (avg, std dev, RSD) 
#(use header file again)
variation_calcs <- function(slopes) {
  #avg already calculated above in slope func, move to somewhere where it can be accessed here?
  replicate_data <- all_data %>% filter(Level == 1)
  native_info <- matrix(nrow = 3, ncol = num_analytes, dimnames = list(c("Average", "Std Dev", "RSD"), rep("", num_analytes)))
  ISTD_info <- native_info
  row_num <- 0
  for (i in seq(6, length(all_col_names) - 1, 2)) {
    analyte <- all_col_names[i]
    ISTD <- all_col_names[i+1]
    row_num <- row_num + 1
    #perform calcs for this analyte
    n_avg <- mean(replicate_data[[analyte]])
    i_avg <- mean(replicate_data[[ISTD]])
    n_sd <- sd(replicate_data[[analyte]])
    i_sd <- sd(replicate_data[[ISTD]])
    n_rsd <- (n_sd / n_avg) * 100
    i_rsd <- (i_sd / i_avg) * 100
    native_info[row_num, "Average"] <- n_avg
    ISTD_info[row_num, "Average"] <- i_avg
    native_info[row_num, "Std Dev"] <- n_sd
    ISTD_info[row_num, "Std Dev"] <- i_sd
    native_info[row_num, "RSD"] <- n_rsd
    ISTD_info[row_num, "RSD"] <- i_rsd
    colnames(analyte_info)[row_num] <- analyte
    colnames(ISTD_info)[row_num] <- ISTD
    
  }
}

conc_calcs <- function(slopes) {
  #measured conc = RR/slope
  
}

linearity <- function(slopes, r2) {

}

accuracy <- function() {
  #measured conc and percent differences
}

LOB <- function(blanks) {
  #LoB = avg conc. + (1.645 * std dev)
  
}

LOD <- function(LOBs) {
  #LoD = LoB + (1.645 * std dev)

}

LOQ <- function(LODs) {
  #LoQ = LoD * 3 (or 0.025 if calculation is below that)
  
}

upload_mv <- function() {
  
}
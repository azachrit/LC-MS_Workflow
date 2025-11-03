# Method Val Draft R Script, 11/3/25
# Alicia Melotik

### use same code from other script to load data, combine via header files later ###

calc_slopes <- function(all_data) {
  slopes <- list()
  
  #blanks removed data, levels forced as numbers
  caldata <- all_data %>% filter(Level > 0)
  caldata <- caldata %>% mutate(Level = as.numeric(Level))
  
  all_col_names <- colnames(all_data)
  #loop through each analyte to calculate each slope
  for (i in seq(6, length(all_col_names), 2)) {
    analyte <- all_col_names[i]
    ISTD <- all_col_names[i+1]
    
    #find N2IRatio across each sample and create vector to store values for lm use
    #N2Iratio <- caldata[[analyte]] / caldata[[ISTD]]
    
    N2Iratio <- vector()
    levels <- vector()
    #use avg N2I for replicates
    for (j in seq(nrow(caldata))) {
      if (!(caldata[j, "Level"] %in% levels)) {
        #find current level
        cur_level <- caldata[j, "Level"]
        replicates <- caldata %>% filter(Level == cur_level)
        levels <- c(levels, cur_level)
        
        #calculate avg ratio for current replicate level
        N2I <- replicates[[analyte]] / replicates[[ISTD]]
        N2Iratio <- c(N2Iratio, mean(N2I))
      }
    }
    
    #linear model regression from draft r script
    cal_reg <- lm(N2Iratio ~ 0 + levels) ## added "0 +" forces line through origin 
    #slope <- coef(cal_reg)[2]
    #slopes <- c(slopes, slope)
    
    slope <- coef(cal_reg)[["levels"]]
    r2 <- summary(cal_reg)$r.squared
    
    # Save results in list
    slopes[[analyte]] <- list(
      slope = slope,
      r2 = r2
    )
    slopes_sorted <- slopes[order(names(slopes))]
  }
}

calc_slopes(all_data)

peak_areas <- function(trial_names, native_areas, ISTD_areas) {
  #averages for replicates, calculate relative response using averages
  avg_native_peak <- list()
  avg_ISTD_peak <- list()

  #for each conc in ppb, find the avg
  #ASSUMPTIONS:
  ##### all names begin with their concentration number in ng/L
  ##### all replicates of same conc are listed continuously, one after the other
  start_index <- 1
  curr_index <- 2
  while (start_index <= length(trial_names) ){
    curr_conc <- grepl("^[0-9]+", trial_names[start_index])
    paste(curr_conc)
    while (curr_conc == grepl("^[0-9]+", trial_names[curr_index])) {
      curr_index <- curr_index + 1
    }
    #replicates currently being processed are from indices [start_index to curr_index] inclusive
    avg_native_peak[as.numeric(curr_conc)/100] <- sum(native_areas[start_index, curr_index])/(curr_index-start_index + 1)
  }
  
}

##variation calculations same for other workflow script (avg, std dev, RSD) 
#(use header file again)

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
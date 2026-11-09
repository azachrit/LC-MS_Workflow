# LC-MS
R scripts to automate the "pipeline" of LC-MS data processing. <br>
Excel sheet folders where data should be uploaded are located in the [SWEL Lab Google Drive](https://drive.google.com/drive/u/0/folders/1Wo-_9w0MelpFjG8tsKr8pQDiBCjwmqy3). <br><br>
  *If you start getting errors, try to redownload the script and run it again. I may have made changes to the file(s) pulled from this Github.*

 ### Dependencies
  Run the following lines of code to install the necessary libraries before running a script for the first time:
  
    install.packages(openxlsx)
    install.packages(googledrive)
    install.packages(tidyverse)

## Method Val Processing
  ### Running the R script locally
  1. Open R Studio [download [here](https://posit.co/download/rstudio-desktop/)]
  
  2. Make sure necessary libraries are downloaded (listed above in "Dependencies")
      - _copy and paste the "install packages" lines into the terminal at the bottom of the R window_
  3. Download method_val_processing.R from this Github and open the file in R.
  4. Export the raw data into the [Processed Method Val Files](https://drive.google.com/drive/u/0/folders/10QyIkD-MeF_0yvfwBlv3FlPDvKlrNkYF) folder in the Google Drive.<br>
       - *copy/upload an Excel file containing one sheet with the raw data into the above folder*
  5. Source the method_val_processing.R script and follow the prompts to authenticate to Google Drive.
  6. To process the most recently uploaded file, press enter. To choose a specific file to process, type or paste the file name and press enter.
  7. When prompted, enter the replicate level to use for LOD/LOQ calculations (in ng/mL or ng/L).
  8. The uploaded Excel sheet will be updated with the processed data within the same folder it was uploaded to earlier.
<br>

### Required Format of Raw Data for method vals
- Blank samples are determined by where the "Level" column is empty or set to 0. Make sure only samples to be used in LOB calculations have this cell blank.
- The columns should be meta data (name, data file, type, level, acq. date/time) followed by interleaved analytes and ISTDs. Each analyte should be followed by its corresponding ISTD.

## Data Processing
  ### Running the R script locally
  1. Open R Studio [download [here](https://posit.co/download/rstudio-desktop/)]
  
  2. Make sure necessary libraries are downloaded (listed above in "Dependencies")
      - *copy and paste the "install packages" lines into the terminal at the bottom of the R window*
  3. Download data_processing.R from this Github and open the file in R.
  4. Export the raw data into the [Processed Data](https://drive.google.com/drive/u/0/folders/1EAaOfnITWIn1KEl3NG1u-HXidcCxF881) folder in the Google Drive.<br>
       - *or copy/upload an Excel file containing one sheet with the raw data into the above folder*
  5. Source the data_processing.R script and follow the prompts to authenticate to Google Drive.
  6. To use the most recently uploaded data file, press enter. To choose a specific file, type or paste the file name and press enter.
  7. To use the copy vals from the most recent method validation file, press enter. To choose a specific method val file, type or paste the file name and press enter.
  8. The Excel sheet will be updated with the processed data within the same folder it was uploaded to earlier.
<br>

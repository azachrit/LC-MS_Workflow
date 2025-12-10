# LC-MS
R scripts to automate the "pipeline" LC-MS data processing
- Excel sheets to upload the data are located in the [SWEL Google Drive](https://drive.google.com/drive/u/0/folders/0ADEP1hOTEkn_Uk9PVA).

## Method Val Processing
  ### Running the R script locally
  1. Download/open R Studio [download [here](https://posit.co/download/rstudio-desktop/)]
  2. Make sure necessary libraries are downloaded (listed below in "Dependencies")
  3. Download method_val_processing.R from this Github.
  4. Paste the raw data into a new Excel sheet within the [Processed Method Val Files](https://drive.google.com/drive/u/0/folders/10QyIkD-MeF_0yvfwBlv3FlPDvKlrNkYF) folder in the Google Drive.
  5. Source the method_val_processing.R script and follow the prompts to authenticate to Google Drive.
  6. The Excel sheet will be updated (& renamed) with the processed data within the same folder it was uploaded to earlier.
  
  ### Dependencies
    install.packages(openxlsx)
    install.packages(googledrive)
    install.packages(tidyverse)
  

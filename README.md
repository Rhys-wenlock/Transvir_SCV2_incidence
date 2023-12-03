# Transvir SARS-CoV-2 Incidence

### Description 

Analysis for TransVIR Sars-CoV-2 Incidence Paper

## Contents

The directory contains:

  - [:file\_folder: data](/Data): The raw input data files. Stored elsewhere.
  - [:file\_folder: scripts](/Scripts): R scripts to be run in number order.
  - [:file\_folder: output](/Output): All output generated from running the scripts. Included any modified data, figures, statistical analyses.
 

## How to replicate analysis

This analysis requires [R software](https://cloud.r-project.org/) and
 [RStudio Desktop](https://rstudio.com/products/rstudio/download/) to be pre-installed.

1. Clone this repository
2. Open RStudio and create new project associated with this directory [See here for details](https://rpubs.com/Dee_Chiluiza/create_RProject).
3. Open each script in `./scripts` and run in number order

##To replicate incidence based analyses (i.e. Anderson-Gill Modelling)
1. Run 01_Data_cleaning.R
2. Import code to produce PCR negative dates (02_PCR_negative_dates.R)
3. Identify all infection episodes (incidence_data_base.R) - this is the default PCR neg definition
4. Convert into a lexis format 04_lexis_code
5. Run Anderson Gill Modelling 05_AG_Results

To run sensitivity analyses, repeat steps 3-5 with difference incidence_data_XX versions of the code
6. Summaryresults_sensitivity.R provides overall numbers for each scenario. 

##To replicate symptom results
1. Run 01_Data_cleaning.R
2. Import code to produce PCR negative dates (02_PCR_negative_dates.R)
3. Identify all infection episodes (incidence_data_base.R) - this is the default PCR neg definition
4. Convert into a lexis format 04_lexis_code
5. Run symptom analysis - 05_symptomcode_V3.R
6. Produce summary description of symptoms - 06_symptom_summarycode.R

#To replicate cluster analysis
1. Run 01_data_cleaning.R
2. Run 02_cluster_code_final.R to produce cluster data sets for analysis
3. Run either 03_HCIR_analysis.R or 03_SAR_analysis.R for HCIR and SAR analyses respectively. 



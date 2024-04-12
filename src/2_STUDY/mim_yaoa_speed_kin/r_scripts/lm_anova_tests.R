install.packages(c("tidyverse","purrr","R.matlab","readxl"))
library(readxl)
library(R.matlab)
library(purrr)
library(tidyverse)

# excel_dir <- "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\AS_dataset\\_studies\\02012023_subjrec_2bounces_rally_serve_human_epochfix_JS_n1p5-1p5\\_figs\\conn\\dDTF08\\R_data\\";
excel_dir <- "R:\\Ferris-Lab\\jsalminen\\Experiments_Funding\\Experiment_6_MIM_OA\\data_saves\\03232023_MIM_OAN70_antsnormalize_iccREMG0p4_powpow0p3_skull0p01\\cluster\\icrej_5\\12\\spec_data\\spec_of_oh\\psd_calcs\\fooof_kinematics_table.xlsx";

read_excel(excel_dir,sheet="Sheet1")
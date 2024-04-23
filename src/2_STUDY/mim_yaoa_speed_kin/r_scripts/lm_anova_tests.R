# install.packages(c("tidyverse","purrr","R.matlab","readxl","dplyr"))
# install.packages("devtools")
library(readxl)
library(R.matlab)
library(purrr)
library(tidyverse)
library(dplyr)

#%% load table
excel_dir <- "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\MIM_dataset\\_studies\\04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01\\cluster\\icrej_5\\12\\spec_data\\group_spec\\psd_calcs\\fooof_kinematics_table.xlsx";
table_out <- read_excel(excel_dir,sheet="Sheet1")
#%% get unique entries
clusters = unique(table_out$cluster_id)
subjects = unique(table_out$subj_char)
eeg_measures = c('alpha_avg_power','beta_avg_power','theta_avg_power','aperiodic_exp','aperiodic_offset');
#%% get speeds only
table_out <- filter_at(table_out,vars('cond_char'), any_vars(. %in% c('0.25','0.5','0.75','1.0')))
#%% convert speeds & groups to factors
table_out <- mutate(table_out,across(c('cond_char','group_char'),
                                     factor))
#%% LOOP through clusters
for (ci in clusters) {
  for (mi in eeg_measures) {
    tmpt <- filter_at(table_out,vars('cluster_id'), any_vars(. %in% ci))
    tapply(tmpt[[mi]],tmpt$cond_char)
    contr.poly(4)
    contrasts(tmpt$cond_char) = contr.poly(4);
    summary(lm(paste(mi,"~ cond_char"),tmpt))
  }
}
# table_out <- filter_at(table_out,vars('cluster_id'), any_vars(. %in% c('3')))
# table_out <- select(table_out,cond_char,cluster_id,group_char,beta_avg_power)
# #%% convert speeds & groups to factors
# table_out <- mutate(table_out,across(c('cond_char','group_char'),
#                 factor))
# # table_out <- mutate(table_out,across(c('cond_char'),
# #                 as.numeric))
# #%% orthogonal polynomial coding
# # table_out$cond_char<-cut(table_out$cond_char, 4,ordered = TRUE)
# # table(table_out)
# tapply(table_out$beta_avg_power,table_out$cond_char)
# contr.poly(4)
# contrasts(table_out$cond_char) = contr.poly(4);
# summary(lm(beta_avg_power~cond_char,table_out))


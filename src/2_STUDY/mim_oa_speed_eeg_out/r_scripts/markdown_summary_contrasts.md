---
title: "Report"
author: "jacob salminen"
date: "April 19, 2024"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Packages

'''{r} 
library(readxl);
library(purrr);
library(tidyverse);
library(dplyr)
'''

## load table 
'''{r} 
excel_dir \<-"M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\\_data\\MIM_dataset\\\_studies\\04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01\\cluster\\icrej_5\\12\\spec_data\\group_spec\\psd_calcs\\fooof_kinematics_table.xlsx";
table_out \<- read_excel(excel_dir,sheet="Sheet1") 
'''

### get unique
entries 
'''{r} 
clusters = unique(table_out$cluster_id) subjects = unique(table_out$subj_char)
eeg_measures = c('alpha_avg_power','beta_avg_power','theta_avg_power','aperiodic_exp','aperiodic_offset');
'''

### get speeds only 
'''{r} table_out
<- filter_at(table_out,vars('cond_char'), any_vars(. %in% c('0.25','0.5','0.75','1.0'))) 
'''
### convert speeds & groups to factors 
'''{r}
table_out <- mutate(table_out,across(c('cond_char','group_char'), factor)) '''
### LOOP through clusters & get constrast summaries 

'''{r} 
for (ci in clusters) { 
  for (mi in eeg_measures) {
    tmpt <- filter_at(table_out,vars('cluster_id'), any_vars(. %in% ci))
    tapply(tmpt[[mi]],tmpt$cond_char)
    contr.poly(4)
    contrasts(tmpt$cond_char) = contr.poly(4); 
    summary(lm(paste(mi,"\~ cond_char"),tmpt)) 
  } 
}
'''

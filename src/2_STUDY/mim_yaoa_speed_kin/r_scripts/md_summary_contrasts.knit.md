---
title: "md_summary_contrasts"
author: "Jacob Salminen"
date: "2024-04-19"
output: pdf_document
---



## Packages

```r
library(readxl);
# library(purrr)
library(tidyverse);
```

```
## -- Attaching core tidyverse packages ------------------------ tidyverse 2.0.0 --
## v dplyr     1.1.4     v readr     2.1.5
## v forcats   1.0.0     v stringr   1.5.1
## v ggplot2   3.5.0     v tibble    3.2.1
## v lubridate 1.9.3     v tidyr     1.3.1
## v purrr     1.0.2     
## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
## i Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

```r
# library(tibble)
library(knitr);
```

## load table 

```r
excel_dir <-"M:/jsalminen/GitHub/par_EEGProcessing/src/_data/MIM_dataset/_studies/04162024_MIM_YAOAN89_antsnormalize_iccREMG0p4_powpow0p3_skull0p01/cluster/icrej_5/12/spec_data/group_spec/psd_calcs/fooof_kinematics_table.xlsx";
table_out <- read_excel(excel_dir,sheet="Sheet1") 
```

### get unique entries 

```r
clusters = unique(table_out$cluster_id);
subjects = unique(table_out$subj_char);
eeg_measures = c('alpha_avg_power','beta_avg_power','theta_avg_power','aperiodic_exp','aperiodic_offset');
```

### get speeds only 

```r
table_out <- filter_at(table_out,vars('cond_char'), any_vars(. %in% c('0.25','0.5','0.75','1.0'))) 
```

### convert speeds & groups to factors 

```r
table_out <- mutate(table_out,across(c('cond_char','group_char'), factor))
head(table_out)
```

```
## # A tibble: 6 x 121
##   speed_ms subj_id subj_cl_ind subj_char comp_id design_id cond_id cond_char
##      <dbl> <chr>         <dbl> <chr>       <dbl> <chr>     <chr>   <fct>    
## 1     0.86 1                 1 H1002           8 2         1       0.25     
## 2     0.87 2                 2 H1004          11 2         1       0.25     
## 3     0.91 3                 3 H1007           8 2         1       0.25     
## 4     0.67 4                 4 H1009           4 2         1       0.25     
## 5     0.78 5                 5 H1010           1 2         1       0.25     
## 6     0.7  7                 6 H1012           5 2         1       0.25     
## # i 113 more variables: group_id <chr>, cluster_id <chr>, aperiodic_exp <dbl>,
## #   aperiodic_offset <dbl>, central_freq_1 <dbl>, central_freq_2 <dbl>,
## #   central_freq_3 <dbl>, power_1 <dbl>, power_2 <dbl>, power_3 <dbl>,
## #   r_squared <dbl>, theta_avg_power <dbl>, alpha_avg_power <dbl>,
## #   beta_avg_power <dbl>, theta_1 <dbl>, theta_2 <dbl>, theta_3 <dbl>,
## #   theta_4 <dbl>, alpha_1 <dbl>, alpha_2 <dbl>, alpha_3 <dbl>, alpha_4 <dbl>,
## #   alpha_5 <lgl>, alpha_6 <lgl>, beta_1 <dbl>, beta_2 <dbl>, beta_3 <dbl>, ...
```

### LOOP through clusters & get constrast summaries 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -4.6089 -1.5618 -0.1448  1.1882  7.0560  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  3.50425    0.12828  27.317   <2e-16 *** ` 
`. cond_char.L -0.20238    0.25656  -0.789    0.431     ` 
`. cond_char.Q -0.03709    0.25656  -0.145    0.885     ` 
`. cond_char.C -0.19416    0.25656  -0.757    0.450     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 2.053 on 252 degrees of freedom ` 
`. Multiple R-squared:  0.004802,	Adjusted R-squared:  -0.007046  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -2.1886 -0.9394 -0.0452  0.8925  3.3452  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.94496    0.07517  25.873   <2e-16 *** ` 
`. cond_char.L -0.04310    0.15035  -0.287    0.775     ` 
`. cond_char.Q -0.09008    0.15035  -0.599    0.550     ` 
`. cond_char.C -0.03007    0.15035  -0.200    0.842     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.203 on 252 degrees of freedom ` 
`. Multiple R-squared:  0.001906,	Adjusted R-squared:  -0.009976  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.5334 -0.8236 -0.3805  0.4921  4.2116  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.990646   0.071164  13.921   <2e-16 *** ` 
`. cond_char.L -0.068890   0.142327  -0.484    0.629     ` 
`. cond_char.Q  0.002479   0.142327   0.017    0.986     ` 
`. cond_char.C -0.046687   0.142327  -0.328    0.743     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.139 on 252 degrees of freedom ` 
`. Multiple R-squared:  0.001356,	Adjusted R-squared:  -0.01053  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.13365 -0.20313 -0.01167  0.18555  0.96105  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.301262   0.017829  72.987   <2e-16 *** ` 
`. cond_char.L -0.038503   0.035657  -1.080    0.281     ` 
`. cond_char.Q -0.015431   0.035657  -0.433    0.666     ` 
`. cond_char.C  0.009113   0.035657   0.256    0.798     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.2853 on 252 degrees of freedom ` 
`. Multiple R-squared:  0.005598,	Adjusted R-squared:  -0.00624  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.83099 -0.35659  0.01263  0.37449  1.37057  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -0.50797    0.03431 -14.805   <2e-16 *** ` 
`. cond_char.L -0.04983    0.06862  -0.726    0.468     ` 
`. cond_char.Q -0.01002    0.06862  -0.146    0.884     ` 
`. cond_char.C  0.01082    0.06862   0.158    0.875     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.549 on 252 degrees of freedom ` 
`. Multiple R-squared:  0.002271,	Adjusted R-squared:  -0.009607  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -2.6267 -1.1757 -0.4039  0.9326  4.9677  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.415211   0.118209  11.972   <2e-16 *** ` 
`. cond_char.L -0.101481   0.236417  -0.429    0.668     ` 
`. cond_char.Q  0.036880   0.236417   0.156    0.876     ` 
`. cond_char.C  0.002545   0.236417   0.011    0.991     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.672 on 196 degrees of freedom ` 
`. Multiple R-squared:  0.001064,	Adjusted R-squared:  -0.01423  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.9824 -0.9907 -0.3215  0.9033  3.0303  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.59238    0.09009  17.675   <2e-16 *** ` 
`. cond_char.L -0.13367    0.18019  -0.742    0.459     ` 
`. cond_char.Q -0.00447    0.18019  -0.025    0.980     ` 
`. cond_char.C  0.01237    0.18019   0.069    0.945     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.274 on 196 degrees of freedom ` 
`. Multiple R-squared:  0.002827,	Adjusted R-squared:  -0.01244  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.4667 -0.7531 -0.1846  0.4160  2.4820  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.941277   0.069739  13.497   <2e-16 *** ` 
`. cond_char.L  0.189844   0.139478   1.361    0.175     ` 
`. cond_char.Q  0.008159   0.139478   0.058    0.953     ` 
`. cond_char.C -0.094604   0.139478  -0.678    0.498     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.9863 on 196 degrees of freedom ` 
`. Multiple R-squared:  0.01168,	Adjusted R-squared:  -0.003449  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -0.92913 -0.19432 -0.00607  0.19529  0.68847  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.04180    0.02103  49.544   <2e-16 *** ` 
`. cond_char.L  0.02637    0.04206   0.627    0.531     ` 
`. cond_char.Q  0.00444    0.04206   0.106    0.916     ` 
`. cond_char.C  0.02048    0.04206   0.487    0.627     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.2974 on 196 degrees of freedom ` 
`. Multiple R-squared:  0.003263,	Adjusted R-squared:  -0.01199  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -2.12300 -0.43800  0.01736  0.60719  1.14942  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -0.922141   0.053134 -17.355   <2e-16 *** ` 
`. cond_char.L  0.032066   0.106269   0.302    0.763     ` 
`. cond_char.Q  0.008143   0.106269   0.077    0.939     ` 
`. cond_char.C  0.031701   0.106269   0.298    0.766     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.7514 on 196 degrees of freedom ` 
`. Multiple R-squared:  0.0009476,	Adjusted R-squared:  -0.01434  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -3.8751 -2.3419 -0.7304  1.9330  8.3644  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  3.14665    0.17670  17.808   <2e-16 *** ` 
`. cond_char.L -0.22385    0.35340  -0.633    0.527     ` 
`. cond_char.Q  0.18090    0.35340   0.512    0.609     ` 
`. cond_char.C -0.05009    0.35340  -0.142    0.887     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 2.893 on 264 degrees of freedom ` 
`. Multiple R-squared:  0.002582,	Adjusted R-squared:  -0.008753  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -3.1851 -1.1615 -0.0239  1.0980  4.1018  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  2.98700    0.09911  30.137   <2e-16 *** ` 
`. cond_char.L -0.22199    0.19823  -1.120    0.264     ` 
`. cond_char.Q  0.04779    0.19823   0.241    0.810     ` 
`. cond_char.C -0.05633    0.19823  -0.284    0.776     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.623 on 264 degrees of freedom ` 
`. Multiple R-squared:  0.005249,	Adjusted R-squared:  -0.006055  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.2380 -0.4501 -0.1846  0.2051  4.7792  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.29335    0.05061   5.796 1.93e-08 *** ` 
`. cond_char.L  0.07583    0.10122   0.749    0.454     ` 
`. cond_char.Q  0.01817    0.10122   0.180    0.858     ` 
`. cond_char.C -0.02415    0.10122  -0.239    0.812     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.8285 on 264 degrees of freedom ` 
`. Multiple R-squared:  0.002458,	Adjusted R-squared:  -0.008878  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -0.75333 -0.11262  0.00342  0.16416  0.68952  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.003912   0.014400  69.715   <2e-16 *** ` 
`. cond_char.L  0.033491   0.028800   1.163    0.246     ` 
`. cond_char.Q  0.002777   0.028800   0.096    0.923     ` 
`. cond_char.C -0.001150   0.028800  -0.040    0.968     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.2357 on 264 degrees of freedom ` 
`. Multiple R-squared:  0.005137,	Adjusted R-squared:  -0.006168  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.40177 -0.27862  0.03705  0.32141  1.52000  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -1.184408   0.031994 -37.019   <2e-16 *** ` 
`. cond_char.L  0.051984   0.063989   0.812    0.417     ` 
`. cond_char.Q -0.003840   0.063989  -0.060    0.952     ` 
`. cond_char.C  0.006379   0.063989   0.100    0.921     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.5238 on 264 degrees of freedom ` 
`. Multiple R-squared:  0.002545,	Adjusted R-squared:  -0.00879  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.7937 -0.7565 -0.2502  0.7686  2.9918  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.605821   0.078166   7.750 9.09e-13 *** ` 
`. cond_char.L  0.050187   0.156332   0.321    0.749     ` 
`. cond_char.Q  0.108147   0.156332   0.692    0.490     ` 
`. cond_char.C -0.001309   0.156332  -0.008    0.993     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.013 on 164 degrees of freedom ` 
`. Multiple R-squared:  0.003534,	Adjusted R-squared:  -0.01469  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.7938 -0.9369 -0.1956  0.6874  3.3820  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.21074    0.08825  13.719   <2e-16 *** ` 
`. cond_char.L -0.06908    0.17651  -0.391    0.696     ` 
`. cond_char.Q  0.03251    0.17651   0.184    0.854     ` 
`. cond_char.C  0.02515    0.17651   0.142    0.887     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.144 on 164 degrees of freedom ` 
`. Multiple R-squared:  0.001263,	Adjusted R-squared:  -0.01701  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.6526 -0.7140 -0.2188  0.4358  3.0909  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.175927   0.075093  15.660   <2e-16 *** ` 
`. cond_char.L  0.189275   0.150186   1.260    0.209     ` 
`. cond_char.Q  0.007568   0.150186   0.050    0.960     ` 
`. cond_char.C -0.017381   0.150186  -0.116    0.908     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.9733 on 164 degrees of freedom ` 
`. Multiple R-squared:  0.009687,	Adjusted R-squared:  -0.008428  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -0.61439 -0.16884 -0.03504  0.19336  0.58301  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.03816    0.01994  52.056   <2e-16 *** ` 
`. cond_char.L  0.02225    0.03989   0.558    0.578     ` 
`. cond_char.Q -0.00407    0.03989  -0.102    0.919     ` 
`. cond_char.C  0.00261    0.03989   0.065    0.948     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.2585 on 164 degrees of freedom ` 
`. Multiple R-squared:  0.001983,	Adjusted R-squared:  -0.01627  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.59441 -0.53811 -0.03078  0.50116  1.19461  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -1.249288   0.052815 -23.654   <2e-16 *** ` 
`. cond_char.L  0.033665   0.105630   0.319    0.750     ` 
`. cond_char.Q -0.010820   0.105630  -0.102    0.919     ` 
`. cond_char.C  0.008511   0.105630   0.081    0.936     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.6846 on 164 degrees of freedom ` 
`. Multiple R-squared:  0.0007224,	Adjusted R-squared:  -0.01756  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -3.3272 -1.6627 -0.1922  1.2921  6.5179  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  2.50000    0.17455  14.322   <2e-16 *** ` 
`. cond_char.L -0.18129    0.34911  -0.519    0.604     ` 
`. cond_char.Q  0.09616    0.34911   0.275    0.783     ` 
`. cond_char.C -0.07619    0.34911  -0.218    0.828     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 2.152 on 148 degrees of freedom ` 
`. Multiple R-squared:  0.00265,	Adjusted R-squared:  -0.01757  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.4689 -0.6931 -0.1146  0.4690  4.1646  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.87671    0.07827  11.201   <2e-16 *** ` 
`. cond_char.L -0.10078    0.15655  -0.644    0.521     ` 
`. cond_char.Q  0.01869    0.15655   0.119    0.905     ` 
`. cond_char.C  0.01250    0.15655   0.080    0.936     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.965 on 148 degrees of freedom ` 
`. Multiple R-squared:  0.002931,	Adjusted R-squared:  -0.01728  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -0.9914 -0.5079 -0.1903  0.3484  2.8234  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.55069    0.05553   9.917   <2e-16 *** ` 
`. cond_char.L  0.05041    0.11106   0.454    0.651     ` 
`. cond_char.Q -0.03339    0.11106  -0.301    0.764     ` 
`. cond_char.C -0.05876    0.11106  -0.529    0.598     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.6846 on 148 degrees of freedom ` 
`. Multiple R-squared:  0.003879,	Adjusted R-squared:  -0.01631  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -0.97670 -0.18134  0.00827  0.19562  1.11074  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.924768   0.029453  31.398   <2e-16 *** ` 
`. cond_char.L  0.006734   0.058907   0.114    0.909     ` 
`. cond_char.Q -0.016935   0.058907  -0.287    0.774     ` 
`. cond_char.C  0.007659   0.058907   0.130    0.897     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.3631 on 148 degrees of freedom ` 
`. Multiple R-squared:  0.0007604,	Adjusted R-squared:  -0.01949  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -2.2833 -0.4239  0.0235  0.4709  2.2919  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -0.892049   0.066619 -13.390   <2e-16 *** ` 
`. cond_char.L  0.050776   0.133238   0.381    0.704     ` 
`. cond_char.Q -0.006801   0.133238  -0.051    0.959     ` 
`. cond_char.C  0.011761   0.133238   0.088    0.930     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.8213 on 148 degrees of freedom ` 
`. Multiple R-squared:  0.00105,	Adjusted R-squared:  -0.0192  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.    Min     1Q Median     3Q    Max  ` 
`. -4.244 -2.248 -0.048  1.720  5.573  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  3.466438   0.223192  15.531   <2e-16 *** ` 
`. cond_char.L -0.320789   0.446383  -0.719    0.474     ` 
`. cond_char.Q  0.004912   0.446383   0.011    0.991     ` 
`. cond_char.C -0.050378   0.446383  -0.113    0.910     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 2.603 on 132 degrees of freedom ` 
`. Multiple R-squared:  0.003994,	Adjusted R-squared:  -0.01864  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.7918 -0.9559 -0.0059  0.7473  3.2750  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.441483   0.096874  14.880   <2e-16 *** ` 
`. cond_char.L -0.182464   0.193749  -0.942    0.348     ` 
`. cond_char.Q -0.042641   0.193749  -0.220    0.826     ` 
`. cond_char.C -0.007101   0.193749  -0.037    0.971     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.13 on 132 degrees of freedom ` 
`. Multiple R-squared:  0.007046,	Adjusted R-squared:  -0.01552  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.4172 -0.7423 -0.4180  0.4214  4.3444  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.658942   0.094665   6.961 1.43e-10 *** ` 
`. cond_char.L  0.026317   0.189330   0.139    0.890     ` 
`. cond_char.Q  0.006489   0.189330   0.034    0.973     ` 
`. cond_char.C -0.045868   0.189330  -0.242    0.809     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.104 on 132 degrees of freedom ` 
`. Multiple R-squared:  0.0005996,	Adjusted R-squared:  -0.02211  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -0.78302 -0.22794  0.02919  0.23962  0.55783  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.911243   0.026257  34.704   <2e-16 *** ` 
`. cond_char.L  0.046993   0.052514   0.895    0.372     ` 
`. cond_char.Q -0.032970   0.052514  -0.628    0.531     ` 
`. cond_char.C -0.003091   0.052514  -0.059    0.953     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.3062 on 132 degrees of freedom ` 
`. Multiple R-squared:  0.008997,	Adjusted R-squared:  -0.01353  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.27035 -0.30844  0.04821  0.34945  1.08823  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -1.009073   0.046622 -21.644   <2e-16 *** ` 
`. cond_char.L  0.104551   0.093243   1.121    0.264     ` 
`. cond_char.Q -0.023716   0.093243  -0.254    0.800     ` 
`. cond_char.C  0.001758   0.093243   0.019    0.985     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.5437 on 132 degrees of freedom ` 
`. Multiple R-squared:  0.009918,	Adjusted R-squared:  -0.01258  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -3.4620 -1.6035 -0.4735  1.0269  9.2229  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.93431    0.14324  13.504   <2e-16 *** ` 
`. cond_char.L -0.09317    0.28647  -0.325    0.745     ` 
`. cond_char.Q  0.11924    0.28647   0.416    0.678     ` 
`. cond_char.C -0.03353    0.28647  -0.117    0.907     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 2.237 on 240 degrees of freedom ` 
`. Multiple R-squared:  0.001218,	Adjusted R-squared:  -0.01127  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -3.1052 -1.7850 -0.2509  1.7671  4.4124  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  2.64009    0.12928  20.421   <2e-16 *** ` 
`. cond_char.L -0.16628    0.25857  -0.643    0.521     ` 
`. cond_char.Q  0.01472    0.25857   0.057    0.955     ` 
`. cond_char.C -0.04258    0.25857  -0.165    0.869     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 2.019 on 240 degrees of freedom ` 
`. Multiple R-squared:  0.001846,	Adjusted R-squared:  -0.01063  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.2924 -0.5654 -0.2042  0.3955  2.9069  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.51278    0.05127  10.002   <2e-16 *** ` 
`. cond_char.L  0.12564    0.10254   1.225    0.222     ` 
`. cond_char.Q -0.01182    0.10254  -0.115    0.908     ` 
`. cond_char.C -0.06699    0.10254  -0.653    0.514     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.8008 on 240 degrees of freedom ` 
`. Multiple R-squared:  0.008024,	Adjusted R-squared:  -0.004376  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -0.94946 -0.16315  0.00568  0.20084  0.56163  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.96610    0.01662  58.140   <2e-16 *** ` 
`. cond_char.L  0.03076    0.03323   0.925    0.356     ` 
`. cond_char.Q  0.01149    0.03323   0.346    0.730     ` 
`. cond_char.C  0.00578    0.03323   0.174    0.862     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.2596 on 240 degrees of freedom ` 
`. Multiple R-squared:  0.004176,	Adjusted R-squared:  -0.008272  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.72503 -0.30032 -0.00651  0.32857  1.23757  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -1.217442   0.035456 -34.336   <2e-16 *** ` 
`. cond_char.L  0.049329   0.070912   0.696    0.487     ` 
`. cond_char.Q  0.009814   0.070912   0.138    0.890     ` 
`. cond_char.C  0.015076   0.070912   0.213    0.832     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.5538 on 240 degrees of freedom ` 
`. Multiple R-squared:  0.002279,	Adjusted R-squared:  -0.01019  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.    Min     1Q Median     3Q    Max  ` 
`. -3.091 -1.956 -1.087  1.831 10.815  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  2.38867    0.17140  13.936   <2e-16 *** ` 
`. cond_char.L -0.18140    0.34281  -0.529    0.597     ` 
`. cond_char.Q  0.19037    0.34281   0.555    0.579     ` 
`. cond_char.C -0.06183    0.34281  -0.180    0.857     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 2.588 on 224 degrees of freedom ` 
`. Multiple R-squared:  0.002764,	Adjusted R-squared:  -0.01059  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.9303 -1.1311 -0.3960  0.9939  3.7190  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.621747   0.090724  17.876   <2e-16 *** ` 
`. cond_char.L -0.125577   0.181448  -0.692    0.490     ` 
`. cond_char.Q  0.062602   0.181448   0.345    0.730     ` 
`. cond_char.C -0.005501   0.181448  -0.030    0.976     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.37 on 224 degrees of freedom ` 
`. Multiple R-squared:  0.002667,	Adjusted R-squared:  -0.01069  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.0885 -0.5032 -0.1954  0.3018  2.7718  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.62845    0.05075  12.383   <2e-16 *** ` 
`. cond_char.L  0.15506    0.10150   1.528    0.128     ` 
`. cond_char.Q  0.04836    0.10150   0.476    0.634     ` 
`. cond_char.C -0.05440    0.10150  -0.536    0.592     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.7663 on 224 degrees of freedom ` 
`. Multiple R-squared:  0.01256,	Adjusted R-squared:  -0.0006693  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -0.84714 -0.16270  0.01075  0.17573  0.54715  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.082314   0.016864  64.181   <2e-16 *** ` 
`. cond_char.L  0.016502   0.033727   0.489    0.625     ` 
`. cond_char.Q -0.031274   0.033727  -0.927    0.355     ` 
`. cond_char.C  0.005322   0.033727   0.158    0.875     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.2546 on 224 degrees of freedom ` 
`. Multiple R-squared:  0.004993,	Adjusted R-squared:  -0.008333  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.61157 -0.41725  0.06029  0.46271  1.27237  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -0.988252   0.041638 -23.734   <2e-16 *** ` 
`. cond_char.L  0.034734   0.083277   0.417    0.677     ` 
`. cond_char.Q -0.034216   0.083277  -0.411    0.682     ` 
`. cond_char.C  0.006772   0.083277   0.081    0.935     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.6287 on 224 degrees of freedom ` 
`. Multiple R-squared:  0.001557,	Adjusted R-squared:  -0.01181  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -3.3390 -1.6877 -0.6289  1.0797  8.8515  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  2.72672    0.15040  18.129   <2e-16 *** ` 
`. cond_char.L -0.18694    0.30081  -0.621    0.535     ` 
`. cond_char.Q  0.16918    0.30081   0.562    0.574     ` 
`. cond_char.C -0.07329    0.30081  -0.244    0.808     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 2.388 on 248 degrees of freedom ` 
`. Multiple R-squared:  0.003063,	Adjusted R-squared:  -0.008997  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.6842 -0.9436 -0.2799  0.6580  3.8062  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.514698   0.078208  19.368   <2e-16 *** ` 
`. cond_char.L -0.061315   0.156415  -0.392    0.695     ` 
`. cond_char.Q -0.008557   0.156415  -0.055    0.956     ` 
`. cond_char.C -0.002743   0.156415  -0.018    0.986     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.242 on 248 degrees of freedom ` 
`. Multiple R-squared:  0.0006325,	Adjusted R-squared:  -0.01146  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -0.9228 -0.3721 -0.1436  0.2717  2.2511  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.411341   0.034744  11.839   <2e-16 *** ` 
`. cond_char.L  0.115775   0.069489   1.666    0.097 .   ` 
`. cond_char.Q -0.007938   0.069489  -0.114    0.909     ` 
`. cond_char.C -0.055756   0.069489  -0.802    0.423     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.5516 on 248 degrees of freedom ` 
`. Multiple R-squared:  0.01365,	Adjusted R-squared:  0.001721  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.35899 -0.16386  0.07761  0.25302  0.64174  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.970483   0.021867  44.381   <2e-16 *** ` 
`. cond_char.L  0.001538   0.043734   0.035    0.972     ` 
`. cond_char.Q -0.011611   0.043734  -0.265    0.791     ` 
`. cond_char.C  0.023720   0.043734   0.542    0.588     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.3471 on 248 degrees of freedom ` 
`. Multiple R-squared:  0.001473,	Adjusted R-squared:  -0.01061  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.8863 -0.3538  0.1550  0.4576  1.0809  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -0.956446   0.038843 -24.623   <2e-16 *** ` 
`. cond_char.L  0.042333   0.077686   0.545    0.586     ` 
`. cond_char.Q -0.004145   0.077686  -0.053    0.957     ` 
`. cond_char.C  0.030195   0.077686   0.389    0.698     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.6166 on 248 degrees of freedom ` 
`. Multiple R-squared:  0.001815,	Adjusted R-squared:  -0.01026  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -3.8908 -1.8458 -0.7996  1.5477  7.7013  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  3.08150    0.16149  19.082   <2e-16 *** ` 
`. cond_char.L -0.37272    0.32298  -1.154    0.250     ` 
`. cond_char.Q  0.06240    0.32298   0.193    0.847     ` 
`. cond_char.C -0.02296    0.32298  -0.071    0.943     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 2.46 on 228 degrees of freedom ` 
`. Multiple R-squared:  0.005991,	Adjusted R-squared:  -0.007089  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -2.8149 -1.0817 -0.3893  1.0458  3.5329  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  2.667006   0.092961  28.690   <2e-16 *** ` 
`. cond_char.L -0.258997   0.185922  -1.393    0.165     ` 
`. cond_char.Q -0.092108   0.185922  -0.495    0.621     ` 
`. cond_char.C -0.008955   0.185922  -0.048    0.962     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.416 on 228 degrees of freedom ` 
`. Multiple R-squared:  0.009507,	Adjusted R-squared:  -0.003526  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.2930 -0.5046 -0.2282  0.2559  4.4439  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.297848   0.053140   5.605 5.98e-08 *** ` 
`. cond_char.L -0.018315   0.106279  -0.172    0.863     ` 
`. cond_char.Q  0.011600   0.106279   0.109    0.913     ` 
`. cond_char.C  0.002016   0.106279   0.019    0.985     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.8094 on 228 degrees of freedom ` 
`. Multiple R-squared:  0.0001841,	Adjusted R-squared:  -0.01297  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -0.70410 -0.19074 -0.04424  0.21814  0.53794  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.006838   0.017479  57.602   <2e-16 *** ` 
`. cond_char.L  0.017090   0.034958   0.489    0.625     ` 
`. cond_char.Q -0.008597   0.034958  -0.246    0.806     ` 
`. cond_char.C -0.002759   0.034958  -0.079    0.937     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.2662 on 228 degrees of freedom ` 
`. Multiple R-squared:  0.001339,	Adjusted R-squared:  -0.0118  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.33073 -0.33141 -0.02141  0.40346  1.14601  ` 
`.  ` 
`. Coefficients: ` 
`.               Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -1.0443211  0.0362100 -28.841   <2e-16 *** ` 
`. cond_char.L  0.0380368  0.0724200   0.525    0.600     ` 
`. cond_char.Q -0.0012050  0.0724200  -0.017    0.987     ` 
`. cond_char.C -0.0001784  0.0724200  -0.002    0.998     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.5515 on 228 degrees of freedom ` 
`. Multiple R-squared:  0.00121,	Adjusted R-squared:  -0.01193  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.8751 -0.9574 -0.2508  0.7058  3.7208  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.05604    0.09439  11.188   <2e-16 *** ` 
`. cond_char.L -0.03939    0.18877  -0.209    0.835     ` 
`. cond_char.Q  0.03159    0.18877   0.167    0.867     ` 
`. cond_char.C -0.01988    0.18877  -0.105    0.916     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.294 on 184 degrees of freedom ` 
`. Multiple R-squared:  0.0004489,	Adjusted R-squared:  -0.01585  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.2808 -0.6547 -0.3580  0.4323  3.3527  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.917638   0.068262  13.443   <2e-16 *** ` 
`. cond_char.L -0.102779   0.136524  -0.753    0.453     ` 
`. cond_char.Q -0.045185   0.136524  -0.331    0.741     ` 
`. cond_char.C -0.002975   0.136524  -0.022    0.983     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.936 on 184 degrees of freedom ` 
`. Multiple R-squared:  0.003665,	Adjusted R-squared:  -0.01258  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -0.9979 -0.6259 -0.2734  0.2832  3.1889  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.70202    0.06676  10.515   <2e-16 *** ` 
`. cond_char.L  0.16223    0.13353   1.215    0.226     ` 
`. cond_char.Q -0.05828    0.13353  -0.436    0.663     ` 
`. cond_char.C -0.01110    0.13353  -0.083    0.934     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.9154 on 184 degrees of freedom ` 
`. Multiple R-squared:  0.009014,	Adjusted R-squared:  -0.007144  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -0.6635 -0.2026 -0.0186  0.1938  0.6166  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.992916   0.020680  48.013   <2e-16 *** ` 
`. cond_char.L  0.013719   0.041360   0.332    0.740     ` 
`. cond_char.Q -0.004810   0.041360  -0.116    0.908     ` 
`. cond_char.C  0.007082   0.041360   0.171    0.864     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.2836 on 184 degrees of freedom ` 
`. Multiple R-squared:  0.0008301,	Adjusted R-squared:  -0.01546  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.8282 -0.5239  0.1718  0.6693  1.2613  ` 
`.  ` 
`. Coefficients: ` 
`.               Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -0.9576265  0.0555651 -17.234   <2e-16 *** ` 
`. cond_char.L  0.0349061  0.1111302   0.314    0.754     ` 
`. cond_char.Q  0.0003066  0.1111302   0.003    0.998     ` 
`. cond_char.C  0.0130703  0.1111302   0.118    0.907     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.7619 on 184 degrees of freedom ` 
`. Multiple R-squared:  0.000611,	Adjusted R-squared:  -0.01568  ` 
 
##EEG Measure alpha_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = alpha_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -2.2536 -1.1404 -0.6822  0.3421  8.8324  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.26115    0.12670   9.954   <2e-16 *** ` 
`. cond_char.L -0.14967    0.25340  -0.591    0.555     ` 
`. cond_char.Q  0.21033    0.25340   0.830    0.407     ` 
`. cond_char.C  0.00406    0.25340   0.016    0.987     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.896 on 220 degrees of freedom ` 
`. Multiple R-squared:  0.004696,	Adjusted R-squared:  -0.008876  ` 
 
##EEG Measure beta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = beta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -2.3191 -1.3923 -0.6921  1.1020  5.0457  ` 
`.  ` 
`. Coefficients: ` 
`.              Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.988062   0.116510  17.063   <2e-16 *** ` 
`. cond_char.L -0.199644   0.233019  -0.857    0.393     ` 
`. cond_char.Q  0.049118   0.233019   0.211    0.833     ` 
`. cond_char.C  0.005371   0.233019   0.023    0.982     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.744 on 220 degrees of freedom ` 
`. Multiple R-squared:  0.003529,	Adjusted R-squared:  -0.01006  ` 
 
##EEG Measure theta_avg_power 
`.  ` 
`. Call: ` 
`. lm(formula = theta_avg_power ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.     Min      1Q  Median      3Q     Max  ` 
`. -1.3052 -0.6680 -0.3322  0.3904  4.2639  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  0.66505    0.07013   9.483   <2e-16 *** ` 
`. cond_char.L  0.13481    0.14026   0.961    0.338     ` 
`. cond_char.Q  0.05151    0.14026   0.367    0.714     ` 
`. cond_char.C -0.05730    0.14026  -0.409    0.683     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 1.05 on 220 degrees of freedom ` 
`. Multiple R-squared:  0.00554,	Adjusted R-squared:  -0.008021  ` 
 
##EEG Measure aperiodic_exp 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_exp ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -0.67605 -0.12607 -0.00989  0.11325  0.74626  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept)  1.00595    0.01405  71.613   <2e-16 *** ` 
`. cond_char.L  0.03732    0.02809   1.329    0.185     ` 
`. cond_char.Q -0.01037    0.02809  -0.369    0.712     ` 
`. cond_char.C  0.01342    0.02809   0.478    0.633     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.2102 on 220 degrees of freedom ` 
`. Multiple R-squared:  0.009587,	Adjusted R-squared:  -0.003919  ` 
 
##EEG Measure aperiodic_offset 
`.  ` 
`. Call: ` 
`. lm(formula = aperiodic_offset ~ 1 + cond_char, data = tmpt) ` 
`.  ` 
`. Residuals: ` 
`.      Min       1Q   Median       3Q      Max  ` 
`. -1.15679 -0.29471 -0.00893  0.24154  1.54605  ` 
`.  ` 
`. Coefficients: ` 
`.             Estimate Std. Error t value Pr(>|t|)     ` 
`. (Intercept) -1.25566    0.03117 -40.284   <2e-16 *** ` 
`. cond_char.L  0.05713    0.06234   0.916    0.360     ` 
`. cond_char.Q -0.01714    0.06234  -0.275    0.784     ` 
`. cond_char.C  0.01933    0.06234   0.310    0.757     ` 
`. --- ` 
`. Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 ` 
`.  ` 
`. Residual standard error: 0.4665 on 220 degrees of freedom ` 
`. Multiple R-squared:  0.004577,	Adjusted R-squared:  -0.008997  ` 
 
# ```{r}
# for (ci in clusters) { 
#   for (mi in eeg_measures) {
#     tmpt <- filter_at(table_out,vars('cluster_id'), any_vars(. %in% ci))
#     tapply(tmpt[[mi]],tmpt$cond_char)
#     contr.poly(4)
#     contrasts(tmpt$cond_char) = contr.poly(4)
#     fit <- lm(paste(mi,"~ cond_char"),tmpt)
#     tidy(fit)
#     modelsummary(fit, output="kableExtra")
#   } 
# }
# ```
# invisible(sapply(seq(eeg_measures), function(i) {
#   fo <- reformulate("cond_char",response=eeg_measures[i])
#   s <- modelsummary(do.call("lm", list(fo,quote(table_out))))
#   cat("\n##EEG Measure",eeg_measures[i],"\n")
#   sapply(1:length(eeg_measures), function(j)
#     cat(paste0("`", ". ", capture.output(s)[j]), "` \n"))
#   cat(" \n")
# }))

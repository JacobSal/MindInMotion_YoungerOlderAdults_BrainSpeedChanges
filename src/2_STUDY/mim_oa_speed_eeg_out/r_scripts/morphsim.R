##########################################################################################################################################
##################################################################### Loading Behavioral Measure Files
DATA_PATH = "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\2_STUDY\\mim_oa_speed_eeg_out\\r_scripts\\BEHAVI";
STAI = read.csv(file.path(DATA_PATH,"STAI_G_X2.csv"))
STAXI = read.csv(file.path(DATA_PATH,"STAXI.csv"))

##########################################################################################################################################
##################################################################### Loading Necessary R Packages

library(parameters)
library(tidyverse)
library(psych)
library(see)
library(lavaan)
library(performance)
library(semPlot)
library(nFactors)
library(vegan)
library(jsonlite)
library(Rtsne)
library(ggplot2)

##########################################################################################################################################
##################################################################### Subject Curation

## metafile
meta = read.csv(file.path(DATA_PATH,"Meta.csv"))

## subject age
youth = subset(meta, Age == "20-25" | Age == "25-30" | Age == "30-35"

## missing structure images or diffusion images
youth = subset(youth, X != "sub-032339" & X != "sub-032341" & X != "sub-032459" & X != "sub-032370"
               & X != "sub-032466" & X != "sub-032438" & X != "sub-032509")

## selecting subjects with no SKID diagnoses
youth$SKID_Diagnoses = as.factor(youth$SKID_Diagnoses)
youth$SKID_Diagnoses = as.numeric(youth$SKID_Diagnoses)
Hyouth = subset(youth, SKID_Diagnoses == 1 | SKID_Diagnoses == 10)

##leaving only relevant info
Hyouth = Hyouth[, -c(4:14, 16:21)]


## old subjects
senior = subset(meta, Age == "60-65" | Age == "65-70" | Age == "70-75")

## selecting subjects with no SKID diagnoses
senior$SKID_Diagnoses = as.factor(senior$SKID_Diagnoses)
senior$SKID_Diagnoses = as.numeric(senior$SKID_Diagnoses)
senior = subset(senior, SKID_Diagnoses == 1 | SKID_Diagnoses == 6)

## missing structure images or diffusion images or STAI measure
senior = subset(senior, X != "sub-032339" & X != "sub-032341" & X != "sub-032459" & X != "sub-032370" 
            & X != "sub-032466" & X != "sub-032438" & X != "sub-032509" & X != "sub-032392" & X != "sub-032443" & X != "sub-032488")

## leaving only relevant info
senior = senior[, -c(4:14, 16:21)]

##########################################################################################################################################
##################################################################### Principal Component Analysis: Negative Affect

## preparing data
NegAff = merge(STAXI, STAI, by = "X")

NegAffHyouth = merge(Hyouth, NegAff, by = "X")
NegAffHyouth = subset(NegAffHyouth, select = c(X, Gender_.1.female_2.male, Age, Hamilton_Scale, STAI_Trait_Anxiety,
                                               STAXI_Trait_Anger))

NegAffsenior = merge(senior, NegAff, by = "X")
NegAffsenior = subset(NegAffsenior, select = c(X, Gender_.1.female_2.male, Age, Hamilton_Scale, STAI_Trait_Anxiety,
                                               STAXI_Trait_Anger))

## excluding subjects without HAM-D score
NegAffHyouth = NegAffHyouth[-119, ]
NegAffsenior = NegAffsenior[-19, ]

## PCA on young adults
negaffHyouth.pca = prcomp(NegAffHyouth[, -c(1:3)], center = T, scale = T)
summary(negaffHyouth.pca)
print(negaffHyouth.pca)
pcHyouth = parameters::n_factors(NegAffHyouth[, -c(1:3)])
as.data.frame(pcHyouth)
pcHyouth

## PCA on old adults
negaffsenior.pca = prcomp(NegAffsenior[, -c(1:3)], center = T, scale = T)
summary(negaffsenior.pca)
print(negaffsenior.pca)
pcsenior = parameters::n_factors(NegAffsenior[, -c(1:3)])
as.data.frame(pcsenior)
pcsenior

## exporting results
Hyouthpc1 = negaffHyouth.pca$x[, 1]
NegAffHyouth$Hyouthpc1 = Hyouthpc1

seniorpc1 = negaffsenior.pca$x[, 1]
NegAffsenior$seniorpc1 = seniorpc1

##########################################################################################################################################
##################################################################### Merging Behavior Data

Hyouth = merge(Hyouth, STAI, by = "X")
Hyouth$NApc1 = Hyouthpc1

senior = merge(senior, STAI, by = "X")
senior$NApc1 = seniorpc1

##########################################################################################################################################
##################################################################### Preparing Behavior Dissimilarity Matrices

HAnnaKSTAI = data.frame()
for (i in 1:119) {
  for (j in 1:119) {
    HAnnaKSTAI[j, i] = (sum(Hyouth[i, 5], Hyouth[j, 5]))/2
  }
}
seniorAnnaKSTAI = data.frame()
for (i in 1:45) {
  for (j in 1:45) {
    seniorAnnaKSTAI[j, i] = (sum(senior[i, 5], senior[j, 5]))/2
  }
}

##########################################################################################################################################
##################################################################### Loading Tract Data

yh_prob_L = read.csv(file.path(DATA_PATH,"yh_probmap_L.csv"), header = T)
yh_prob_R = read.csv(file.path(DATA_PATH,"yh_probmap_R.csv"), header = T)
oh_prob_L = read.csv(file.path(DATA_PATH,"oh_probmap_L.csv"), header = T)
oh_prob_R = read.csv(file.path(DATA_PATH,"oh_probmap_R.csv"), header = T)
yh_prob_5p_L = read.csv(file.path(DATA_PATH,"yh_probmap_5p_L.csv"), header = T)
yh_prob_5p_R = read.csv(file.path(DATA_PATH,"yh_probmap_5p_R.csv"), header = T)
oh_prob_5p_L = read.csv(file.path(DATA_PATH,"oh_probmap_5p_L.csv"), header = T)
oh_prob_5p_R = read.csv(file.path(DATA_PATH,"oh_probmap_5p_R.csv"), header = T)
yh_prob_10p_L = read.csv(file.path(DATA_PATH,"yh_probmap_10p_L.csv"), header = T)
yh_prob_10p_R = read.csv(file.path(DATA_PATH,"yh_probmap_10p_R.csv"), header = T)
oh_prob_10p_L = read.csv(file.path(DATA_PATH,"oh_probmap_10p_L.csv"), header = T)
oh_prob_10p_R = read.csv(file.path(DATA_PATH,"oh_probmap_10p_R.csv"), header = T)
yh_OINE_prob_L = read.csv(file.path(DATA_PATH,"yh_OccInfNE_probmap_L.csv"), header = T)
yh_OINE_prob_R = read.csv(file.path(DATA_PATH,"yh_OccInfNE_probmap_R.csv"), header = T)
yh_OINE_prob_5p_L = read.csv(file.path(DATA_PATH,"yh_OccInfNE_probmap_5p_L.csv"), header = T)
yh_OINE_prob_5p_R = read.csv(file.path(DATA_PATH,"yh_OccInfNE_probmap_5p_R.csv"), header = T)
yh_OINE_prob_10p_L = read.csv(file.path(DATA_PATH,"yh_OccInfNE_probmap_10p_L.csv"), header = T)
yh_OINE_prob_10p_R = read.csv(file.path(DATA_PATH,"yh_OccInfNE_probmap_10p_R.csv"), header = T)
oh_OINE_prob_L = read.csv(file.path(DATA_PATH,"oh_OccInfNE_probmap_L.csv"), header = T)
oh_OINE_prob_R = read.csv(file.path(DATA_PATH,"oh_OccInfNE_probmap_R.csv"), header = T)
oh_OINE_prob_5p_L = read.csv(file.path(DATA_PATH,"oh_OccInfNE_probmap_5p_L.csv"), header = T)
oh_OINE_prob_5p_R = read.csv(file.path(DATA_PATH,"oh_OccInfNE_probmap_5p_R.csv"), header = T)
oh_OINE_prob_10p_L = read.csv(file.path(DATA_PATH,"oh_OccInfNE_probmap_10p_L.csv"), header = T)
oh_OINE_prob_10p_R = read.csv(file.path(DATA_PATH,"oh_OccInfNE_probmap_10p_R.csv"), header = T)
yh_prob_L = yh_prob_L[, -1]
yh_prob_R = yh_prob_R[, -1]
oh_prob_L = oh_prob_L[, -1]
oh_prob_R = oh_prob_R[, -1]
yh_prob_5p_L = yh_prob_5p_L[, -1]
yh_prob_5p_R = yh_prob_5p_R[, -1]
oh_prob_5p_L = oh_prob_5p_L[, -1]
oh_prob_5p_R = oh_prob_5p_R[, -1]
yh_prob_10p_L = yh_prob_10p_L[, -1]
yh_prob_10p_R = yh_prob_10p_R[, -1]
oh_prob_10p_L = oh_prob_10p_L[, -1]
oh_prob_10p_R = oh_prob_10p_R[, -1]
yh_OINE_prob_L = yh_OINE_prob_L[, -1]
yh_OINE_prob_R = yh_OINE_prob_R[, -1]
yh_OINE_prob_5p_L = yh_OINE_prob_5p_L[, -1]
yh_OINE_prob_5p_R = yh_OINE_prob_5p_R[, -1]
yh_OINE_prob_10p_L = yh_OINE_prob_10p_L[, -1]
yh_OINE_prob_10p_R = yh_OINE_prob_10p_R[, -1]
oh_OINE_prob_L = oh_OINE_prob_L[, -1]
oh_OINE_prob_R = oh_OINE_prob_R[, -1]
oh_OINE_prob_5p_L = oh_OINE_prob_5p_L[, -1]
oh_OINE_prob_5p_R = oh_OINE_prob_5p_R[, -1]
oh_OINE_prob_10p_L = oh_OINE_prob_10p_L[, -1]
oh_OINE_prob_10p_R = oh_OINE_prob_10p_R[, -1]

##########################################################################################################################################
##################################################################### Main Analyses

mantel(HAnnaKSTAI, dist(t(yh_prob_L)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_prob_R)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_prob_5p_L)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_prob_5p_R)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_prob_10p_L)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_prob_10p_R)), method = "spearman", permutations = 10000)

##########################################################################################################################################
##################################################################### Control Analyses

## tract integrity approach
FALy = read.csv(file.path(DATA_PATH,"Hyouth_FA_L.csv"))
FARy = read.csv(file.path(DATA_PATH,"Hyouth_FA_R.csv"))
Hyouth = merge(Hyouth, FALy, by = "X")
Hyouth = merge(Hyouth, FARy, by = "X")

cor.test(Hyouth$FALy, Hyouth$STAI_Trait_Anxiety)
cor.test(Hyouth$FARy, Hyouth$STAI_Trait_Anxiety)

TractIntegrityLy = data.frame()
for (i in 1:119) {
  for (j in 1:119) {
    TractIntegrityLy[j, i] = abs(Hyouth[i, 8]-Hyouth[j, 8])
  }
}
TractIntegrityRy = data.frame()
for (i in 1:119) {
  for (j in 1:119) {
    TractIntegrityRy[j, i] = abs(Hyouth[i, 9]-Hyouth[j, 9])
  }
}
mantel(TractIntegrityLy, HAnnaKSTAI, method = "spearman", permutations = 10000)
mantel(TractIntegrityRy, HAnnaKSTAI, method = "spearman", permutations = 10000)

## ROI dissimilarity
yh_ROImap_L = read.csv(file.path(DATA_PATH,"yh_ROImap_L.csv"), header = T)
yh_ROImap_R = read.csv(file.path(DATA_PATH,"yh_ROImap_R.csv"), header = T)
yh_ROImap_L = yh_ROImap_L[, -1]
yh_ROImap_R = yh_ROImap_R[, -1]
mantel(dist(t(yh_ROImap_L)), HAnnaKSTAI, method = "spearman", permutations = 10000)
mantel(dist(t(yh_ROImap_R)), HAnnaKSTAI, method = "spearman", permutations = 10000)

## magnitude matching
logscaled = yh_prob_5p_L
for (i in 1:119) {
  for (j in 1:14770) {
    if (yh_prob_5p_L[j, i] < 2) {
      next
    } else {
      logscaled[j, i] = log10((yh_prob_5p_L[j, i])-1)+1
    }
  }
}

Rlogscaled = yh_prob_5p_R
for (i in 1:119) {
  for (j in 1:15027) {
    if (yh_prob_5p_R[j, i] < 2) {
      next
    } else {
      Rlogscaled[j, i] = log10((yh_prob_5p_R[j, i])-1)+1
    }
  }
}
mantel(HAnnaKSTAI, dist(t(logscaled)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(Rlogscaled)), method = "spearman", permutations = 10000)

## control tract
mantel(HAnnaKSTAI, dist(t(yh_OINE_prob_L)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_OINE_prob_R)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_OINE_prob_5p_L)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_OINE_prob_5p_R)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_OINE_prob_10p_L)), method = "spearman", permutations = 10000)
mantel(HAnnaKSTAI, dist(t(yh_OINE_prob_10p_R)), method = "spearman", permutations = 10000)

## demographics
yAnnaKSex = data.frame()
for (i in 1:119) {
  for (j in 1:119) {
    yAnnaKSex[j, i] = abs(Hyouth[i, 2]-Hyouth[j, 2])
  }
}
mantel(yAnnaKSex, HAnnaKSTAI, method = "spearman", permutations = 10000)
HAnnaKAgeabs = data.frame()
for (i in 1:119) {
  for (j in 1:119) {
    HAnnaKAgeabs[j, i] = abs(Hyouth[i, 3]-Hyouth[j, 3])
  }
}
mantel(HAnnaKAgeabs, HAnnaKSTAI, method = "spearman", permutations = 10000)

## head movement
young_qc = fromJSON(file.path(DATA_PATH,"QC_young/squad/group_db.json"))
young_qc_param = as.matrix(cbind(young_qc$qc_motion, young_qc$qc_parameters))
mantel(dist(young_qc_param[, c(3:8)]), HAnnaKSTAI, method = "spearman", permutations = 10000)

## negative affect
HAnnaKPmean = data.frame()
for (i in 1:118) {
  for (j in 1:118) {
    HAnnaKPmean[j, i] = mean(NegAffHyouth[i, 7], NegAffHyouth[j, 7])
  }
}
mantel(HAnnaKPmean, dist(t(yh_prob_5p_L[, -119])), method = "spearman", permutations = 10000)
mantel(HAnnaKPmean, dist(t(yh_prob_5p_R[, -119])), method = "spearman", permutations = 10000)

##########################################################################################################################################
##################################################################### Main Analyses: Old Sample

mantel(seniorAnnaKSTAI, dist(t(oh_prob_L)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_prob_R)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_prob_5p_L)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_prob_5p_R)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_prob_10p_L)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_prob_10p_R)), method = "spearman", permutations = 10000)

##########################################################################################################################################
##################################################################### Control Analyses: Old Sample

## tract integrity approach
FALo = read.csv(file.path(DATA_PATH,"senior_FA_L.csv"))
FARo = read.csv(file.path(DATA_PATH,"senior_FA_R.csv"))
senior = merge(senior, FALo, by = "X")
senior = merge(senior, FARo, by = "X")

cor.test(senior$FALo, senior$STAI_Trait_Anxiety)
cor.test(senior$FARo, senior$STAI_Trait_Anxiety)

TractIntegrityLo = data.frame()
for (i in 1:45) {
  for (j in 1:45) {
    TractIntegrityLo[j, i] = abs(senior[i, 8]-senior[j, 8])
  }
}
TractIntegrityRo = data.frame()
for (i in 1:45) {
  for (j in 1:45) {
    TractIntegrityRo[j, i] = abs(senior[i, 9]-senior[j, 9])
  }
}
mantel(TractIntegrityLo, seniorAnnaKSTAI, method = "spearman", permutations = 10000)
mantel(TractIntegrityRo, seniorAnnaKSTAI, method = "spearman", permutations = 10000)

## ROI dissimilarity
oh_ROImap_L = read.csv(file.path(DATA_PATH,"oh_ROImap_L.csv"), header = T)
oh_ROImap_R = read.csv(file.path(DATA_PATH,"oh_ROImap_R.csv"), header = T)
oh_ROImap_L = oh_ROImap_L[, -1]
oh_ROImap_R = oh_ROImap_R[, -1]
mantel(dist(t(oh_ROImap_L)), seniorAnnaKSTAI, method = "spearman", permutations = 10000)
mantel(dist(t(oh_ROImap_R)), seniorAnnaKSTAI, method = "spearman", permutations = 10000)

## magnitude matching
ohlogscaled = as.data.frame(oh_prob_5p_L)
for (i in 1:45) {
  for (j in 1:11662) {
    if (oh_prob_5p_L[i, j] < 2) {
      next
    } else {
      ohlogscaled[i, j] = log10((oh_prob_5p_L[i, j])-1)+1
    }
  }
}

Rohlogscaled = as.data.frame(oh_prob_5p_R)
for (i in 1:45) {
  for (j in 1:12506) {
    if (oh_prob_5p_R[i, j] < 2) {
      next
    } else {
      Rohlogscaled[i, j] = log10((oh_prob_5p_R[i, j])-1)+1
    }
  }
}
mantel(seniorAnnaKSTAI, dist(ohlogscaled), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(Rohlogscaled), method = "spearman", permutations = 10000)

## control tract
mantel(seniorAnnaKSTAI, dist(t(oh_OINE_prob_L)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_OINE_prob_R)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_OINE_prob_5p_L)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_OINE_prob_5p_R)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_OINE_prob_10p_L)), method = "spearman", permutations = 10000)
mantel(seniorAnnaKSTAI, dist(t(oh_OINE_prob_10p_R)), method = "spearman", permutations = 10000)

mantel.partial(seniorAnnaKSTAI, dist(t(oh_prob_5p_L)), dist(t(oh_OINE_prob_5p_L)),
               method = "spearman", permutations = 10000)
mantel(dist(t(oh_prob_5p_L)), dist(t(oh_OINE_prob_5p_L)), method = "spearman", permutations = 10000)

## demographics
seniorAnnaKAgeabs = data.frame()
for (i in 1:45) {
  for (j in 1:45) {
    seniorAnnaKAgeabs[j, i] = abs(senior[i, 3]-senior[j, 3])
  }
}
mantel(seniorAnnaKAgeabs, seniorAnnaKSTAI, method = "spearman", permutations = 10000)
mantel.partial(seniorAnnaKSTAI, dist(t(oh_prob_5p_L)), seniorAnnaKAgeabs,
               method = "spearman", permutations = 10000)
seniorAnnaKSex = data.frame()
for (i in 1:45) {
  for (j in 1:45) {
    seniorAnnaKSex[j, i] = abs(senior[i, 2]-senior[j, 2])
  }
}
mantel(seniorAnnaKSex, seniorAnnaKSTAI, method = "spearman", permutations = 10000)

## head movement
old_qc = fromJSON(file.path(DATA_PATH,"QC_old/squad/group_db.json")
old_qc_param = as.matrix(cbind(old_qc$qc_motion, old_qc$qc_parameters))
mantel(dist(old_qc_param[, c(3:8)]), seniorAnnaKSTAI, method = "spearman", permutations = 10000)

## negative affect
seniorAnnaKPmean = data.frame()
for (i in 1:44) {
  for (j in 1:44) {
    seniorAnnaKPmean[j, i] = mean(NegAffsenior[i, 7], NegAffsenior[j, 7])
  }
}
mantel(seniorAnnaKPmean, dist(t(oh_prob_5p_L)[-19, ]), method = "spearman", permutations = 10000)

##########################################################################################################################################
##################################################################### Visualization: Correlogram

theme_set(theme_classic())

## young
correlog = mantel.correlog(dist(t(yh_prob_5p_L)), HAnnaKSTAI, r.type ="spearman", mult = "fdr")
plot(correlog)
correlog$mantel.res

youngcorrelog = correlog$mantel.res
youngcorrelog = subset(youngcorrelog, youngcorrelog[, 1] < 41)
youngcorrelog = cbind(youngcorrelog, c(2, 2, 2, 2, 1, 1, 1, 1))
youngcorrelog = as.data.frame(youngcorrelog)
p = ggplot(data = youngcorrelog,
           mapping = aes(x = class.index,
                         y = Mantel.cor)) + geom_line(size = 2)  + geom_hline(yintercept = 0) + geom_point(pch = 19,
                                                                                                           color = c("black", "red")[youngcorrelog$V6],
                                                                                                           size = 8) + geom_vline(xintercept = c(25, 30, 35, 40), linetype = 3)
print(p)
write.csv(youngcorrelog, file = file.path(DATA_PATH,"youngcorrelogram.csv"))

## old
correlog2 = mantel.correlog(dist(t(oh_prob_5p_L)), seniorAnnaKSTAI, r.type ="spearman", mult = "fdr")
plot(correlog2)
correlog2$mantel.res

oldcorrelog = correlog2$mantel.res
oldcorrelog = subset(oldcorrelog, oldcorrelog[, 1] < 37)
oldcorrelog = cbind(oldcorrelog, c(2, 2, 1, 1, 1, 1))
oldcorrelog = as.data.frame(oldcorrelog)
p2 = ggplot(data = oldcorrelog,
           mapping = aes(x = class.index,
                         y = Mantel.cor)) + geom_line(size = 2)  + geom_hline(yintercept = 0) + geom_point(pch = 19,
                                                                                                           color = c("black", "red")[oldcorrelog$V6],
                                                                                                           size = 8) + geom_vline(xintercept = c(24, 28, 32, 36), linetype = 3)
print(p2)
write.csv(oldcorrelog, file = file.path(DATA_PATH,"oldcorrelogram.csv"))

##########################################################################################################################################
##################################################################### Visualization: TSNE-plot

theme_set(theme_void(18))

## young
set.seed(119)
tSNE_fit = as.data.frame(as.matrix(dist(t(yh_prob_5p_L)))) %>%
  select(where(is.numeric)) %>%
  Rtsne(is_distance = T, perplexity=33, max_iter=10000)

tSNE_df = tSNE_fit$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(X=Hyouth$X)

tSNE_df = tSNE_df %>%
  inner_join(Hyouth, by = "X")

tSNE_df %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             color = STAI_Trait_Anxiety))+
  geom_point(size = 9)+
  scale_color_stepsn(colors = c("red", "white", "steelblue"), breaks = c(32, 36))+
  theme(legend.position="bottom")
ggsave("tSNE_plot_yh.png")

## old
set.seed(45)
tSNE_fit = as.data.frame(as.matrix(dist(oh_prob_5p_L))) %>%
  select(where(is.numeric)) %>%
  Rtsne(is_distance = T, perplexity=10, max_iter=10000)

tSNE_df = tSNE_fit$Y %>%
  as.data.frame() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(X=senior$X)

tSNE_df = tSNE_df %>%
  inner_join(senior, by = "X")

tSNE_df %>%
  ggplot(aes(x = tSNE1,
             y = tSNE2,
             color = STAI_Trait_Anxiety))+
  geom_point(size=11)+
  scale_color_stepsn(colors = c("red", "white", "steelblue"),
                     breaks = c(28, 31))+
  theme(legend.position="bottom")
ggsave("tSNE_plot_oh.png")


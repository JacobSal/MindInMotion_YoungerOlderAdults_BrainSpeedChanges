# install.packages(c("tidyverse","purrr", "R.matlab","circlize"))
library(circlize)
library(R.matlab)
library(purrr)
library(tidyverse)

# mat_dir <- "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\AS_dataset\\_studies\\12282023_bounces_1h_2bm_JS_n1-0p5\\_figs\\conn\\dDTF08\\R_data\\";
# mat_dir <- "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\AS_dataset\\_studies\\12282023_bounces_1h_2bm_JS_n1-0p5\\_figs\\conn\\S\\R_data\\";
# mat_dir <- "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\AS_dataset\\_studies\\01182023_subjrec_2bounces_1h_2bm_JS_n5-1p5\\_figs\\conn\\S\\R_data\\";
# mat_dir <- "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\AS_dataset\\_studies\\01182023_subjrec_2bounces_1h_2bm_JS_n5-1p5\\_figs\\conn\\dDTF08\\R_data\\";
# mat_dir <- "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\AS_dataset\\_studies\\01252023_subjrec_2bounces_rally_serve_human_JS_n5-1p5\\_figs\\conn\\dDTF08\\R_data\\";
# mat_dir <- "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\AS_dataset\\_studies\\01292023_subjrec_2bounces_rally_serve_human_JS_n0p75-0p75\\_figs\\conn\\dDTF08\\R_data\\";
mat_dir <- "M:\\jsalminen\\GitHub\\par_EEGProcessing\\src\\_data\\AS_dataset\\_studies\\01312023_subjrec_2bounces_rally_serve_human_JS_n0p75-0p75\\_figs\\conn\\dDTF08\\R_data\\";


# CLUSTER ASSIGNMENTS
#- test 1
# cluster_ints =   c(1,2,3,4,5,6,7,8,9,10);
# cluster_names = c('RPPa-Oc','Cuneus','Precuneus','RFrontal','LPPa-Oc','LSM','RSM','LTemp','SuppMotor','LFrontal');
# cluster_names = c("RPPa-Oc","LPPa-Oc","Precuneus","Cuneus")
# sub_ints = c(1,5,3,2)
# sub_ints = c(1,2,3,4,5,6,7,9,10);
#- test 2
cluster_ints =   c(1,2,3,4,5,6,7,8,9);
# cluster_names = c('RPPa-Oc','Cuneus','Precuneus','RFrontal','LPPa-Oc','LSM','RSM','SuppMotor','LFrontal');
# sub_ints = c(1,2,3,4,5,6,7,8,9);
# cluster_names = c('RPPa-Oc','RFrontal','LPPa-Oc','LSM','RSM','LFrontal');
# sub_ints = c(1,4,5,6,7,9);
cluster_names = c('RSM','RPPa-Oc','LPPa-Oc','LSM','LFrontal','RFrontal')
sub_ints = c(7,1,5,6,9,4)
# cluster_names = c("RPPa-Oc","LPPa-Oc","Precuneus","Cuneus")
# sub_ints = c(1,5,3,2)
# (01/24/2024) JS, removing ltemporal from analysis

# CONDITIONS = c("Human_Subject_hit","BM_Subject_hit");
# CONDITIONS = c("2Bounce_BM_Subject_receive","2Bounce_Human_Subject_receive");
# CONDITIONS = c("1Bounce_Human_Subject_hit","Serve_Human_Subject_hit");
CONDITIONS = c("2Bounce_Human_Subject_hit","1Bounce_Human_Subject_hit","Serve_Human_Subject_hit");

#%% PLOT PARAMS
# CIRC_XLIM = c(0,20)
# MULTI_FACTOR = 10^2
CIRC_XLIM = c(0,5)
MULTI_FACTOR = 10^4
CIRC_MAJOR_TICKS = c(0,1,2,3,4,5); #c(0,2.5,5,7.5,10,15,20);
CIRC_MINOR_TICK = 2
ALPHA = 0.05;

#%% SET LOOP PARAMETERS
cli = cluster_ints
clj = cluster_ints
for(cond in CONDITIONS){
  for(plotBand_cort in c("theta","alpha","beta","all")){
    netData_sbjs = readMat(paste(mat_dir,"netVals_",cond,"_aveTime_sbjs.mat",sep=""))
    
    # head(netData_sbjs$netVals[1])
    # df = as.data.frame(netData_sbjs$netVals)
    
    #Determine significant connections above 0
    pvals_theta = matrix(0L, nrow = 1, ncol = length(cli)*length(cli)-length(cli))
    pvals_alpha = matrix(0L, nrow = 1, ncol = length(cli)*length(cli)-length(cli))
    pvals_all = matrix(0L, nrow = 1, ncol = length(cli)*length(cli)-length(cli))
    pvals_beta = matrix(0L, nrow = 1, ncol = length(cli)*length(cli)-length(cli))
    counter = 0
    counter_all = 0
    for(i in cli){
      for(j in clj){
        counter_all = counter_all+1
        if(i!=j){
          counter=counter+1
          pvals_theta[counter] = t.test(netData_sbjs$SAVEVAR[[1]][[counter_all]][[1]])$p.value
          pvals_alpha[counter] = t.test(netData_sbjs$SAVEVAR[[2]][[counter_all]][[1]])$p.value
          pvals_all[counter] = t.test(netData_sbjs$SAVEVAR[[4]][[counter_all]][[1]])$p.value
          pvals_beta[counter] = t.test(netData_sbjs$SAVEVAR[[3]][[counter_all]][[1]])$p.value
        }
      }
    }
    
    #Reshape into 16 x 16
    pvals_theta_box = matrix(1L, nrow = length(cli), ncol = length(cli))
    pvals_alpha_box = matrix(1L, nrow = length(cli), ncol = length(cli))
    pvals_all_box = matrix(1L, nrow = length(cli), ncol = length(cli))
    pvals_beta_box = matrix(1L, nrow = length(cli), ncol = length(cli))
    counter = 0
    for(i in cli){
      for(j in clj){
        if(i!=j){
          counter=counter+1
          pvals_theta_box[j,i] = pvals_theta[counter]
          pvals_alpha_box[j,i] = pvals_alpha[counter]
          pvals_all_box[j,i] = pvals_all[counter]
          pvals_beta_box[j,i] = pvals_beta[counter]
        }
      }
    }
    # ANOVA
    
    
    #FDR correction
    for(i in cli){
      pvals_theta_box[,i]=p.adjust(pvals_theta_box[,i],method='fdr')
      pvals_alpha_box[,i]=p.adjust(pvals_alpha_box[,i],method='fdr')
      pvals_beta_box[,i]=p.adjust(pvals_beta_box[,i],method='fdr')
      pvals_all_box[,i]=p.adjust(pvals_all_box[,i],method='fdr')
    }
    
    #Check for significance
    isSig_theta_box = matrix(0L, nrow = length(cli), ncol = length(cli))
    isSig_alpha_box = matrix(0L, nrow = length(cli), ncol = length(cli))
    isSig_all_box = matrix(0L, nrow = length(cli), ncol = length(cli))
    isSig_beta_box = matrix(0L, nrow = length(cli), ncol = length(cli))
    counter = 0
    for(i in cli){
      for(j in clj){
        if(i!=j){
          counter=counter+1
          if (!is.nan(pvals_theta_box[i,j]) & (pvals_theta_box[i,j]<ALPHA)){
            isSig_theta_box[i,j]=1
          }
          if(!is.nan(pvals_alpha_box[i,j]) & pvals_alpha_box[i,j]<ALPHA){
            isSig_alpha_box[i,j]=1
          }
          if(!is.nan(pvals_all_box[i,j]) & pvals_all_box[i,j]<ALPHA){
            isSig_all_box[i,j]=1
          }
          if(!is.nan(pvals_beta_box[i,j]) & pvals_beta_box[i,j]<ALPHA){
            isSig_beta_box[i,j]=1
          }
        }
      }
    }
    
    #Load average data and mask it by significance
    netData_ave = readMat(paste(mat_dir,"netVals_",cond,"_aveTime.mat",sep=""))
    thetaAve = netData_ave$SAVEVAR[[1]]
    alphaAve = netData_ave$SAVEVAR[[2]]
    allAve = netData_ave$SAVEVAR[[4]]
    betaAve = netData_ave$SAVEVAR[[3]]
    thetaAve[isSig_theta_box==0]=0
    alphaAve[isSig_alpha_box==0]=0
    allAve[isSig_all_box==0]=0
    betaAve[isSig_all_box==0]=0
    #Chord plot (for cortico-cortical data)
    library(circlize)
    
    
    #Rearrange order in matrix
    # 
    cortTheta = thetaAve[1:length(sub_ints),1:length(sub_ints)]
    cortAlpha = alphaAve[1:length(sub_ints),1:length(sub_ints)]
    cortAll = thetaAve[1:length(sub_ints),1:length(sub_ints)]
    cortBeta = thetaAve[1:length(sub_ints),1:length(sub_ints)]
    for(i in 1:length(sub_ints)){
      for(j in 1:length(sub_ints)){
        cortTheta[i,j] = thetaAve[sub_ints[i],sub_ints[j]]
        cortAlpha[i,j] = alphaAve[sub_ints[i],sub_ints[j]]
        cortAll[i,j] = allAve[sub_ints[i],sub_ints[j]]
        cortBeta[i,j] = betaAve[sub_ints[i],sub_ints[j]]
      }
    }
    
    if(plotBand_cort=='theta'){
      dat_in = cortTheta
    }else if(plotBand_cort=='alpha'){
      dat_in = cortAlpha
    }else if(plotBand_cort=='all'){
      dat_in = cortAll
    }else if(plotBand_cort=='beta'){
      dat_in = cortBeta
    }
    #%%
    circos.clear()
    png(paste(mat_dir,cond,"_",plotBand_cort,".jpg",sep=""),bg="transparent",width = 4096, height =3160,res=2)
    color_vals = c('purple','orange','blue','magenta','gold','red','cyan','green','yellow','tan','black','pink')#'#FF00FF','#FF6600','#0000FF','#FF0000','#D4AA00','#800080','#00FF00','#00FFFF')
    color_vals_links = adjustcolor(color_vals, alpha.f = 0.5)
    # circos.par(cell.padding = c(0,0,0,0),gap.degree=1)
    # circos.par(cell.padding = c(0,0,0,0),start.degree = 30)
    circos.par(start.degree = 30,canvas.xlim=c(-1.2,1.2),canvas.ylim=c(-1.2,1.2))
    sign_vals = sign(dat_in)
    A = abs(dat_in)*MULTI_FACTOR
    circos.initialize(cluster_names, xlim = CIRC_XLIM) #cbind(mag_upper, mag_lower)) 
    
    par(family="Times New Roman")
    circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05,
                           #panel.fun for each sector
                           panel.fun = function(x, y) {
                             offset = 24
                             #select details of current sector
                             name = get.cell.meta.data("sector.index")
                             i = get.cell.meta.data("sector.numeric.index")
                             xlim = get.cell.meta.data("xlim")
                             ylim = get.cell.meta.data("ylim")
                             
                             #text direction (dd) and adjusmtents (aa)
                             ang1 = 180; #90;
                             ang2 = 360; #270;
                             
                             theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                             # dd <- ifelse(theta < ang1 || theta > ang2, "clockwise", "reverse.clockwise")
                             dd <- ifelse(theta < ang1 || theta > ang2, "inside", "outside")
                             aa = c(0.5, 0.5)
                             if(theta < ang1 || theta > ang2)  aa = c(0.5, 0.25)
                             
                             #plot cortical labels
                             circos.text(x=mean(xlim),family="serif", y=5, labels=name, facing = dd, cex=600,  adj = aa, col=color_vals[i])
                             
                             #plot main sector
                             circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                         col = color_vals[i], border=color_vals[i])
                             #plot axis
                             circos.axis(labels.cex=300,direction = "outside", major.tick=TRUE, major.tick.length=0.6, major.at=CIRC_MAJOR_TICKS, labels.facing=dd, #seq(from=0,to=floor(12),by=4), 
                                         minor.ticks=CIRC_MINOR_TICK, lwd=9) # , labels.away.percentage = 0.15
                           }) 
    
    # Add nonsignificant links first
    curr_starts_gray=matrix(0L, nrow = 1, ncol = length(sub_ints))
    for(i in 1:length(sub_ints)) {
      for(j in 1:length(sub_ints)){
        if((A[i,j]!=0 | A[j,i]!=0) & (A[j,i]>A[i,j])){
          if(sign_vals[i,j]==0){
            if(plotBand_cort=='theta'){
              data_temp = abs(netData_ave$SAVEVAR[[1]])*MULTI_FACTOR
              circos.link(cluster_names[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), cluster_names[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.8)) 
            }else if(plotBand_cort=='alpha'){
              data_temp = abs(netData_ave$SAVEVAR[[2]])*MULTI_FACTOR
              circos.link(cluster_names[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), cluster_names[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.8)) 
            }else if(plotBand_cort=='all'){
              data_temp = abs(netData_ave$SAVEVAR[[4]])*MULTI_FACTOR
              circos.link(cluster_names[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), cluster_names[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.8)) 
            }else if(plotBand_cort=='beta'){
              data_temp = abs(netData_ave$SAVEVAR[[3]])*MULTI_FACTOR
              circos.link(cluster_names[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), cluster_names[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.8)) 
            }
            curr_starts_gray[i] = curr_starts_gray[i]+data_temp[i,j]
            curr_starts_gray[j] = curr_starts_gray[j]+data_temp[j,i]
          }
        }
      }
    }
    
    curr_starts=matrix(0L, nrow = 1, ncol = length(sub_ints))
    for(i in 1:length(sub_ints)) {
      for(j in 1:length(sub_ints)){
        if((A[i,j]!=0 | A[j,i]!=0) & (A[j,i]>A[i,j])){
          if(sign_vals[i,j]>0 & sign_vals[j,i]>0){
            circos.link(cluster_names[i], c(curr_starts[i],curr_starts[i]+A[i,j]), cluster_names[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('red', alpha.f = 0.7))#color_vals_links[8]) 
            curr_starts[i] = curr_starts[i]+A[i,j]
            curr_starts[j] = curr_starts[j]+A[j,i]
          }else if(sign_vals[i,j]<0 & sign_vals[j,i]<0){
            circos.link(cluster_names[i], c(curr_starts[i],curr_starts[i]+A[i,j]), cluster_names[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('blue', alpha.f = 0.7)) 
            curr_starts[i] = curr_starts[i]+A[i,j]
            curr_starts[j] = curr_starts[j]+A[j,i]
          }else if(sign_vals[i,j]!=0 & sign_vals[j,i]!=0){
            circos.link(cluster_names[i], c(curr_starts[i],curr_starts[i]+A[i,j]), cluster_names[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('magenta', alpha.f = 0.7)) 
            curr_starts[i] = curr_starts[i]+A[i,j]
            curr_starts[j] = curr_starts[j]+A[j,i]
          }
        }
      }
    }
    # circos.xaxis();
    dev.off()
    #%%
  }
}

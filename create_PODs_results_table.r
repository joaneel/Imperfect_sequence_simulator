################################################################################################
# This R script gathers statistics from sequences for all created datasets and calculate
# the first and second moments (mean and variance) of all statistics for each dataset.
# 
# Joane Elleouet joane.elleouet@alumni.ubc.ca
#################################################################################################


#!/usr/bin/env Rscript

args=commandArgs(trailingOnly=TRUE)

dir=args[1]
npods=as.numeric(args[2])
mod=args[3]
nseq=args[4]

params=read.table("/data/user/trueparam.txt",header=T)
params2=params[1:npods,]

means_names = c("m_segs_1","m_segs_2","m_segs","m_segs_1_pr","m_segs_2_pr","m_pi_1","m_pi_2","m_theta_pi","m_theta_w_1","m_theta_w_2","m_theta_w","m_tajimasD_1","m_tajimasD_2","m_tajimasD","m_ZnS_1","m_ZnS_2","m_ZnS","m_perc_shared_1_2","m_perc_private_1_2","m_perc_fixed_dif_1_2","m_pairwise_fst_1_2","m_FayWuH_1","m_FayWuH_2","m_FayWuH","m_dvk","m_dvh","m_thomson_est_1","m_thomson_est_2","m_thomson_est","m_thomson_var_1","m_thomson_var_2","m_thomson_var")
vars_names = c("v_segs_1","v_segs_2","v_segs","v_segs_1_pr","v_segs_2_pr","v_pi_1","v_pi_2","v_theta_pi","v_theta_w_1","v_theta_w_2","v_theta_w","v_tajimasD_1","v_tajimasD_2","v_tajimasD","v_ZnS_1","v_ZnS_2","v_ZnS","v_perc_shared_1_2","v_perc_private_1_2","v_perc_fixed_dif_1_2","v_pairwise_fst_1_2","v_FayWuH_1","v_FayWuH_2","v_FayWuH","v_dvk","v_dvh","v_thomson_est_1","v_thomson_est_2","v_thomson_est","v_thomson_var_1","v_thomson_var_2","v_thomson_var")

tabstats=as.data.frame(matrix(vector(),nrow=npods, ncol=64, dimnames=list(c(),c(means_names,vars_names))))

for (i in 1:npods){
        stats=read.table(paste(dir,"POD",i,"_allmarkers_sumstats.txt",sep=""),header=T)
        tabstats[i,]=stats[1,]
}

final_tab=cbind(params2,tabstats)

write.table(final_tab,file=paste(dir,"PODs_result_table.txt",sep=""),quote=FALSE, row.names=FALSE, sep="\t")

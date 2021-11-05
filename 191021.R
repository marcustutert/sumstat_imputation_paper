
#Run the imputation
fizi_impute = function(){
  system("fizi impute hapgen2_sim_data/genotyped_sumstats_munged hapgen2_sim_data/ref_panel --out hapgen2_sim_data/imputed_sumstats --min-prop 0.0001 --verbose --window-size 1499205  --ridge-term 0.9")
}

#Keep only SNPs in the sumsstats file!
convert_HAPGEN_to_PLINK = function(){
  system("/apps/well/plink/2.00a-20170625/plink2 --extract hapgen2_sim_data/snps_to_extract --haps hapgen2_sim_data/AFR_1000G_subset.haps --legend hapgen2_sim_data/AFR_1000G_subset.legend 22  --allow-extra-chr --make-bed --out hapgen2_sim_data/ref_panel --maf 0.00001")
  #Read in the bim file
  library(data.table)
  bim = fread("hapgen2_sim_data/ref_panel.bim")
  bim$V1 = rep(22,nrow(bim))
  write.table(bim,"hapgen2_sim_data/ref_panel.bim",quote = F,row.names = F, col.names = F)
}

#Run the imputation
fizi_impute = function(){
  system("fizi impute hapgen2_sim_data/genotyped_sumstats_munged hapgen2_sim_data/ref_panel --out hapgen2_sim_data/imputed_sumstats --min-prop 0.0001 --verbose --window-size 1499205  --ridge-term 0.9")
}
evaluate_InferLD_sumstat_accuracy()
evaluate_IMPG_accuracy()


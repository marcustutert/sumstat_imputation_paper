# Write a function that takes in a table of the number of individuals and the ancestries and outputs the HAPGEN VCF files
# Note that for now we will only simulate from chr22 for speed reasons\
# Note that our regions will be EXCLUSIVE (ie: if we have a region we will include the SNPs inside that region and not beyond)

# Take 1000G Plink(?) files and extract the 1000G populations we care about (note this will extract the full panel)
# Also will extract the region that we care about as well (defaults to the first 1Mb chunk in the genome)


#####Source InferLD code#####
# system("R CMD build InferLD")
# install.packages("InferLD_0.1.0.tar.gz", repos = NULL)
# library(InferLD)
#system("R CMD check InferLD")

###### Helper Functions ######
filter_and_match_SNPs = function(reference_haplotypes,
                                 reference_legend,
                                 munged_genotyped_summary_statistics,
                                 inference = T){
  
  #Load in the three files we need to process this data
  reference_haplotypes        = t(fread(reference_haplotypes, header = F)) #Haplotype format (no headers though!)
  nhaps                       = nrow(reference_haplotypes)
  reference_legend            = fread(reference_legend, header = T)        #Legend file (see IMPUTE2)
  genotyped_sumstats          = fread(munged_genotyped_summary_statistics, header = T)    #Munged sumstats
  #Make sure that that genotyped sumstats are ordered correctly
  genotyped_sumstats <- genotyped_sumstats[order(BP),]
  #Next, we need to throw out anything in the reference panel that is at 0% MAF
  non_segregating_snps = which(colMeans(reference_haplotypes) < 0.0000001 | colMeans(reference_haplotypes) > 0.99999 ) #Find all non-segregating SNPs in the reference panel
  
  #Remove non-segregating SNPs (if they are any!)
  if (length(non_segregating_snps) > 0) {
    reference_legend               = reference_legend[-non_segregating_snps,]
    #Remove from reference panel
    reference_haplotypes           = reference_haplotypes[,-non_segregating_snps]
  }
  
  #Remove SNPs in the genotyped sumstats that are not in the reference panel
  snps_in_genotyped_sumstats_not_in_reference_panel = which(!genotyped_sumstats$SNP %in% reference_legend$id)

    if (length(snps_in_genotyped_sumstats_not_in_reference_panel) > 0 ) {
    genotyped_sumstats = genotyped_sumstats[-which(!genotyped_sumstats$SNP %in% reference_legend$id),]
    }
  
  #### In cases when we need the SNPs in reference and sumstats to match exactly (eg. for inference)
  if (inference == T) {
    #Remove SNPs in reference panel and legend that are not in the sumstats
    snps_in_reference_panel_not_in_genotyped_sumstats = which(!reference_legend$id %in%  genotyped_sumstats$SNP)
    if (length(snps_in_reference_panel_not_in_genotyped_sumstats) > 0 ) {
      reference_haplotypes = reference_haplotypes[,-snps_in_reference_panel_not_in_genotyped_sumstats]
      reference_legend     = reference_legend[-snps_in_reference_panel_not_in_genotyped_sumstats,]
    }
  }
  
  #Write out the files with filtered in the title
  write.table(x = genotyped_sumstats,file =  "hapgen2_sim_data/genotyped_sumstats_filtered", quote = F, col.names = T, row.names = F)
  write.table(x = reference_legend, file = "hapgen2_sim_data/reference_legend_filtered", quote = F, col.names = T, row.names = F)
  write.table(x = reference_haplotypes, file = "hapgen2_sim_data/reference_haplotypes_filtered", quote = F, col.names = F, row.names = F)
  
}


subset_ancestry = function(ancestry,
                           region)  # Vector of ancestries (strings): eg: c("EUR","AFR")

{
  library(data.table)
  #Read in the sample file for the 1000G project
  KG_sample = fread(file = "/well/mcvean/mtutert/sumstat_imputation_paper/1000G/Impute2/1000GP_Phase3.sample", header = T,)
  for (i in 1:length(ancestry)) { #Loop through ancestries and subset the 1000G files
    #Extract from the sample file the IDs that we need
    ancestry_id = KG_sample$ID[(KG_sample$GROUP == ancestry[i])]
    #Save the IDs 
    ids = cbind(ancestry_id,ancestry_id)
    names = c("FID","IID")
    data = rbind(names,ids)
    write.table(file = "hapgen2_sim_data/ancestry_id", x = data, quote = F, row.names = F, col.names = F)
    system(command = sprintf("/apps/well/plink/2.00a-20170724/plink2 --pfile /well/mcvean/mtutert/sumstat_imputation_paper/1000G/plink/1000G_chr22_plink --from-bp 1 --to-bp 17100000 --chr 22 --keep hapgen2_sim_data/ancestry_id --export hapslegend --out hapgen2_sim_data/%s_1000G_subset",ancestry))
  }
}

#subset_ancestry("EUR")
#Generate the recombination rate maps given the data
generate_HAPGEN2_recomb_map = function(ancestry,
                                       region  = c(1,17100000)){

  #Read in the .map format from the 1000G data
  library(data.table)
  KG_map_file = fread("/well/mcvean/mtutert/sumstat_imputation_paper/1000G/Impute2/genetic_map_chr22_combined_b37.txt", header = T)
  #Extract the positions that we use in the given region 
  start_position = region[1]
  end_position   = region[2]
  #Now find the closest SNP position to this
  start_map_position = (which(abs(KG_map_file$position - start_position) == min(abs(KG_map_file$position - start_position))))
  end_map_position   = (which(abs(KG_map_file$position - end_position) == min(abs(KG_map_file$position - end_position))))
  #Write out this filtered .map file
  write.table(KG_map_file[start_map_position:end_map_position,], file = sprintf("/well/mcvean/mtutert/sumstat_imputation_paper/sumstat_imputation_paper/hapgen2_sim_data/%s_filtered.map", ancestry),quote = F,row.names = F, col.names = T)
}

#Write a function that generates the HAPGEN2 GWAS Data
run_HAPGEN2 = function(ancestry = "EUR"){
  #Need to choose some SNPs to be causal, read in legend file
  #legend = fread(sprintf("hapgen2_sim_data/full_%s_test.legend", ancestry),quote = F,header = T)
  #random_snp = sample(1:nrow(legend),size = 1)
  snp    = 16054740
  #Write out the causal SNPs so we can check it with FINEMAP accuracy
  system(sprintf("/apps/eb/skylake/software/HAPGEN2/2.2.0/hapgen2 -m ./hapgen2_sim_data/%s_filtered.map -l ./hapgen2_sim_data/%s_1000G_subset.legend -h hapgen2_sim_data/%s_1000G_subset.haps -o hapgen2_sim_data/HAPGEN2_%s_simulated -dl %s 1 1 1 -n 1000 1000", ancestry,ancestry,ancestry,ancestry, snp))
}


#Running SNPTEST (for comparison)
run_SNPTEST = function(){
  system(sprintf("/apps/well/snptest/2.5/snptest -data ./hapgen2_sim_data/HAPGEN2_EUR_simulated.controls.gen ./hapgen2_sim_data/HAPGEN2_EUR_simulated.controls.sample ./hapgen2_sim_data/HAPGEN2_EUR_simulated.cases.gen ./hapgen2_sim_data/HAPGEN2_EUR_simulated.cases.sample -o ./hapgen2_sim_data/EUR_GWAS_SNPTEST_results -frequentist 1 -method score -pheno pheno"))
}

#Mask region
munge_sumstats = function(){
  #Read in full summary statistics from SNPTEST
  full_sumstats = fread("./hapgen2_sim_data/EUR_GWAS_SNPTEST_results", skip = 10, header = T)
  #Make things MUNGEABLE
  #Generate Z-score column
  full_sumstats$Z = full_sumstats$frequentist_add_beta_1/full_sumstats$frequentist_add_se_1
  #Generate CHR column
  full_sumstats$CHR = rep(22,nrow(full_sumstats))
  #Convert Position to BP
  colnames(full_sumstats)[colnames(full_sumstats) == 'position'] <- 'BP'
  colnames(full_sumstats)[colnames(full_sumstats) == 'alleleA'] <- 'A2'
  colnames(full_sumstats)[colnames(full_sumstats) == 'alleleB'] <- 'A1'
  colnames(full_sumstats)[colnames(full_sumstats) == 'frequentist_add_beta_1'] <- 'BETA'
  colnames(full_sumstats)[colnames(full_sumstats) == 'frequentist_add_se_1'] <- 'SE'
  colnames(full_sumstats)[colnames(full_sumstats) == 'frequentist_add_pvalue'] <- 'P'
  #Write out the table
  write.table(full_sumstats,"hapgen2_sim_data/full_sumstats_mungeable",quote = F,row.names = F, col.names = T)
  system("fizi munge hapgen2_sim_data/full_sumstats_mungeable --out hapgen2_sim_data/full_sumstats_munged --N 2000 --ignore Z")
  #Unzip 
  system("gzip -d hapgen2_sim_data/full_sumstats_munged.sumstats.gz")
  #Read in the munged sumstats
  munged_full_sumstats = fread("hapgen2_sim_data/full_sumstats_munged.sumstats", header = T)
  #Generate file of SNPs BP/Position that we will writeout for the reference panel later
  snps_extract = munged_full_sumstats$SNP
  write.table(snps_extract,"hapgen2_sim_data/snps_to_extract", quote = F, row.names = F, col.names = F)
}

#Munge the sumstats together
mask_snps = function(){
  #Read in full summary statistics from SNPTEST
  full_sumstats = fread("hapgen2_sim_data/full_sumstats_munged.sumstats", header = T)
  #Mask some of the data
  set.seed(231)
  genotyped_snps = sample(1:nrow(full_sumstats),size = 0.7*nrow(full_sumstats))
  genotyped_sumstats = full_sumstats[genotyped_snps,]
  reference_legend = fread("hapgen2_sim_data/ref_panel.bim", header = F)
  #Remove SNPs in the genotyped sumstats file that are not in the reference
  genotyped_sumstats = genotyped_sumstats[-which(!genotyped_sumstats$SNP %in% reference_legend$V2),] 
  #Write out table
  write.table(genotyped_sumstats, "hapgen2_sim_data/genotyped_sumstats_munged", quote = F, col.names = T, row.names = F)
}

#Convert VCF to PLINK
#Keep only SNPs in the sumsstats file!
convert_HAPGEN_to_PLINK = function(){
  system("/apps/well/plink/2.00a-20170625/plink2 --extract hapgen2_sim_data/snps_to_extract --data hapgen2_sim_data/HAPGEN2_EUR_simulated.controls --allow-extra-chr --make-bed --out hapgen2_sim_data/ref_panel --maf 0.00001")
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

#Run simulated benchmarking for a single region and single population group 

#Evaluate how we did
evaluate_IMPG_accuracy = function(){
  #Easiest to plot the correlation (r2) for ALL, common (>5%) and rare (<5%) SNPs
  #Load in the unmasked full dataset
  full_sumstats = fread("hapgen2_sim_data/full_sumstats_munged.sumstats", header = T)
  #Load in the imputed stuff from fizie
  fizi_imputed = fread("hapgen2_sim_data/imputed_sumstats.sumstat", header = T)
  #Merge together the data
  comparison_set = inner_join(fizi_imputed,full_sumstats, by = "BP") %>% filter(TYPE == "imputed")
  
  #PLot the r2 and graph for the full dataset
  r2_correlation_full = signif(summary(lm(comparison_set$Z.x~comparison_set$Z.y))$r.squared,4)
  #Plot and mark the r2 correlation)
  all_correlation = plot_ly(x = comparison_set$Z.x, y = comparison_set$Z.y, type = "scatter", mode = "markers", themes="Catherine") %>%
    layout(showlegend = F, annotations = list(text = sprintf("r2 = %s",r2_correlation_full), showarrow = F, x= 2,
                                              y= -2,
                                              xref = "x",
                                              yref = "y",
                                              showarrow = T,
                                              ax = 20,
                                              ay = -40), 
                                              xaxis = list(title = 'Imputed Zscores'), 
                                              yaxis = list(title = 'True (Genotyped) Zscores'),
                                              title = "Correlation True to Imputed Zscores (ALL)")
  
  #PLot the r2 and graph for the common variants
  snptest_results = fread("hapgen2_sim_data/full_sumstats_mungeable", header = T)
  comparison_set_common = inner_join(comparison_set,snptest_results, by = "BP") %>% filter(all_maf > 0.05)
  r2_correlation_common = signif(summary(lm(comparison_set_common$Z.x~comparison_set_common$Z.y))$r.squared,4)
  common_correlation = plot_ly(x = comparison_set_common$Z.x, y = comparison_set_common$Z.y, type = "scatter", mode = "markers", themes="Catherine") %>%
    layout(showlegend = F, annotations = list(text = sprintf("r2 = %s",r2_correlation_common), showarrow = F, x= 2,
                                              y= -2,
                                              xref = "x",
                                              yref = "y",
                                              showarrow = T,
                                              ax = 20,
                                              ay = -40),  
           xaxis = list(title = 'Imputed Zscores'), 
           yaxis = list(title = 'True (Genotyped) Zscores'),
           title = "Correlation True to Imputed Zscores (Common)")
  comparison_set_rare = inner_join(comparison_set,snptest_results, by = "BP") %>% filter(all_maf < 0.05)
  r2_correlation_rare = signif(summary(lm(comparison_set_rare$Z.x~comparison_set_rare$Z.y))$r.squared,4)
  rare_correlation = plot_ly(x = comparison_set_rare$Z.x, y = comparison_set_rare$Z.y, type = "scatter", mode = "markers", themes="Catherine") %>%
    layout(showlegend = F, annotations = list(text = sprintf("r2 = %s",r2_correlation_rare), showarrow = F, x= 2,
                                              y= -2,
                                              xref = "x",
                                              yref = "y",
                                              showarrow = T,
                                              ax = 20,
                                              ay = -40), 
           xaxis = list(title = 'Imputed Zscores'), 
           yaxis = list(title = 'True (Genotyped) Zscores'),
           title = "Correlation True to Imputed Zscores (Rare)")
  correlation = subplot(all_correlation,rare_correlation,common_correlation, shareX = T, shareY = T) %>% 
    layout(title = list(text = "IMP-G Imputation Correlation"), xaxis = list(title = "Imputed Zscores"),yaxis = list(title = "True Zscores"))
  return(correlation)
}

benchmark_sumstat_imputation = function(ancestry,region){
  subset_ancestry(ancestry = ancestry, region = region)
  generate_HAPGEN2_recomb_map(ancestry = ancestry, region = region)
  run_HAPGEN2(ancestry = ancestry)
  run_SNPTEST()
  munge_sumstats()
  mask_snps()
  convert_HAPGEN_to_PLINK()
  fizi_impute()
  evaluate_IMPG_accuracy()
}

sumstat_impute = function(typed_snps,
                          untyped_snps_index,
                          typed_z_scores,
                          LD)
{
  LD_typed_untyped   = LD[typed_snps,untyped_snps_index] #LD of typed - untyped pairs
  inv_LD_typed       = solve(LD[typed_snps,typed_snps] + 0.9*diag(length(typed_snps)))#inverse of LD of typed SNPs
  W                  = LD[untyped_snps_index, typed_snps] %*% inv_LD_typed #these are the weights that turn typed to imputed
  infos              = as.numeric(rowSums(W * LD[untyped_snps_index,typed_snps])) #info measures per each untyped
  z.imp              = (W %*% (typed_z_scores))/sqrt(infos) #use scaling 1/sqrt(infos) to get var=1 for z-scores
  return((z.imp))
}


#This function calculates LD between all SNP's
LD_Matrix = function(haplotypes){
  
  #Takes in haplotype matrix (dim of haps x snps) to calculate various LD metrics (r here)
  haplotypes = t(haplotypes)
  pAB = haplotypes %*% t( haplotypes ) / ncol( haplotypes)
  pA  = rowMeans( haplotypes )
  pB  = rowMeans( haplotypes )
  D  = pAB - pA %*% t( pB )
  denA = sqrt( 1 / pA / ( 1 - pA ) )
  denB = sqrt( 1 / pB / ( 1 - pB ) )
  r  = t( D * denA ) * denB
  return(r)
}

#Run InferLD on the munged sumstats and some associated reference panel
run_InferLD = function(munged_genotyped_summary_statistics,  #These are the genotyped summary statistics in munged format
                       reference_haplotypes,                 #These are the reference haplotypes (no column header)
                       reference_legend,
                       case_control_constant)                     #These is the legend file (IMPUTE2 format)
  
{

  #Read in the filtered files (haplotypes/sumstats/reference)
  munged_genotyped_summary_statistics_filtered = fread("hapgen2_sim_data/genotyped_sumstats_filtered", header = T)
  reference_haplotypes_filtered                = (fread("hapgen2_sim_data/reference_haplotypes_filtered", header = F))
  reference_legend_filtered                    = fread("hapgen2_sim_data/reference_legend_filtered", header = T)
  
  #Get GWAS variance from the munged_genotyped_summary_statistics file
  GWAS_variance = munged_genotyped_summary_statistics_filtered$SE^2
  
  ##Check if we already have saved the RData file of the inference results 
  if (!file.exists("inference_results/inference.RData")) {
    inference_results = LD_inference(ref_panel_haplotypes = reference_haplotypes_filtered,
                                     fst                   = 0.1,
                                     alpha                 = 1e3,
                                     nSamples              = 3,
                                     weights_resolution    = 10,
                                     likelihood_toggle     = T, 
                                     gwas_variance         = GWAS_variance,
                                     case_control_constant = 500,
                                     BurnIn                = 0.9,
                                     debug                 = T)
  saveRDS(inference_results, "inference_results/inference.RData")
  }
  
  #Read in inference results
  x = readRDS("inference_results/inference.RData")
  
  #Get the diagnostics from this inference run
  #Extract the post-burnin weights
  Gibbs_Results_Posterior = x$Gibbs_Array[,floor((0.9*ncol(x$Gibbs_Array))):ncol(x$Gibbs_Array)]
  #Get the average of the weights
  wts = rowMeans(Gibbs_Results_Posterior, na.rm = T)
  # 
  # #Calculate our inferred LD (across all SNPs) with these weights
  # filter_and_match_SNPs(munged_genotyped_summary_statistics = munged_genotyped_summary_statistics,
  #                       reference_haplotypes                = reference_haplotypes,
  #                       reference_legend                    = reference_legend, 
  #                       inference = F)
  # 
  #Read in the filtered files (haplotypes/sumstats/reference)
  munged_genotyped_summary_statistics = fread("hapgen2_sim_data/genotyped_sumstats_filtered", header = T)
  reference_haplotypes                = (fread("hapgen2_sim_data/reference_haplotypes_filtered", header = F))
  reference_legend                    = fread("hapgen2_sim_data/reference_panel.legend", header = T)
  reference_haplotypes_all            = t((fread("hapgen2_sim_data/reference_panel.haps", header = F)))
  
  #Calculate LD given the wts
  inferred_LD = (cov.wt(x = reference_haplotypes_all, wt = wts, cor = T))$cor
  # 
  #Get the index of all typed SNPs (SNPs in the sumstats and in the reference panel)
  typed_snps_index   = which(reference_legend$id %in% munged_genotyped_summary_statistics$SNP)
  untyped_snps_index = setdiff(1:ncol(reference_haplotypes_all),typed_snps_index)
  #Perform sumstat imputation (by hand)
  imputed_zs = sumstat_impute(typed_snps         = typed_snps_index,
                              untyped_snps_index = untyped_snps_index,
                              typed_z_scores     = munged_genotyped_summary_statistics$Z,
                              LD                 = inferred_LD)
  #Extract imputed_zs
  imputed_results = cbind(reference_legend[untyped_snps_index,]$id,imputed_zs)
  colnames(imputed_results) = c("SNP","Z")
  print(dim(imputed_results))
  write.table(imputed_results, file = "imputation_results/imputation_test", quote = F, col.names = T, row.names = F)
  #inference_diagnostics()
  # return()
}

#Function to see how well inference does at recapturing the r2 correlation compared to the reference panel
inference_diagnostics = function(reference_legend) #Read in the full reference legend file
{

  #Read in the inference results
  inference_results = readRDS("inference_results/inference.RData")
  #Calculate AF correlation (post-burnIn)
  Gibbs_Array = inference_results$Gibbs_Array
  Gibbs_Results_Posterior = Gibbs_Array[,(0.9*ncol(Gibbs_Array)):ncol(Gibbs_Array)]
  #Get the AF correlation and plot this
  inferred_af = rowMeans(inference_results$inferred_af_given_weights[,(0.9*ncol(inference_results$inferred_af_given_weights)):ncol(inference_results$inferred_af_given_weights)])
  #Read in the SNPTEST file in order
  SNPTEST     = fread("hapgen2_sim_data/EUR_GWAS_SNPTEST_results", skip = 10)
  #Read in the reference legend (not-filtered)
  reference_legend_filtered = fread("hapgen2_sim_data/reference_legend_filtered", header = T)
  #Extract the SNPs from reference legend file within the SNPTEST file
  true_afs = SNPTEST$all_maf[which(SNPTEST$rsid %in% reference_legend_filtered$id)]
  #Convert AFs to the MAF 
  inferred_af[inferred_af > 0.5] <- 1 - inferred_af[inferred_af > 0.5]
  #Get correlation
  r2_correlation_full = signif(summary(lm(inferred_af~true_afs))$r.squared,3)
  #Plot both the correlation with the true allele frequencies and the imputed allele frequencies (included and withheld SNPs)
  reference_haplotypes_filtered = fread("hapgen2_sim_data/reference_haplotypes_filtered", header = F)
  reference_af         = colMeans(reference_haplotypes_filtered)
  #Convert AFs to the MAF 
  reference_af[reference_af > 0.5] <- 1 - reference_af[reference_af > 0.5]

  #Also find the AF of the imputed SNPs (those not included in InferLD inference process)
  #Find out what SNP IDs those are
  #Read in the non-filtered SNP IDs
  reference_legend_all     = fread("hapgen2_sim_data/reference_panel.legend", header = T)   #All SNPs 
  reference_haplotypes_all = t(fread("hapgen2_sim_data/reference_panel.haps", header = F))  #All SNPs 
  
  #Find SNPs in reference legend that are not in the filtered file
  snps_ids_imputed        = reference_legend_all[which(!reference_legend_all$id %in% reference_legend_filtered$id)]$id
  snp_ids_imputed_index   = which(!reference_legend_all$id %in% reference_legend_filtered$id)
  #Get the inferred allele frequencies from these imputed SNPs
  #Read in inference results
  x = readRDS("inference_results/inference.RData")
  #Get the diagnostics from this inference run
  #Extract the post-burnin weights
  Gibbs_Results_Posterior = x$Gibbs_Array[,floor((0.9*ncol(x$Gibbs_Array))):ncol(x$Gibbs_Array)]
  #Get the average of the weights
  wts = rowMeans(Gibbs_Results_Posterior, na.rm = T)/sum(rowMeans(Gibbs_Results_Posterior, na.rm = T))
  inferred_af_imputed      = colSums(as.matrix(reference_haplotypes_all)[,c(snp_ids_imputed_index)] * wts)
  reference_imputed_af     = colMeans(as.matrix(reference_haplotypes_all)[,c(snp_ids_imputed_index)])
  true_imputed_af          = SNPTEST$all_maf[which(SNPTEST$rsid %in% snps_ids_imputed)]
  #Flip the MAF to less than 0.5
  #Convert AFs to the MAF 
  inferred_af_imputed[inferred_af_imputed > 0.5] <- 1 - inferred_af_imputed[inferred_af_imputed > 0.5]
  reference_imputed_af[reference_imputed_af > 0.5] <- 1 - reference_imputed_af[reference_imputed_af > 0.5]
  true_imputed_af[true_imputed_af > 0.5] <- 1 - true_imputed_af[true_imputed_af > 0.5]
  
  #Get all the necessary r2s #TO-DO
  
  #Plot the results
  all_afs = plot_ly(y = inferred_af, x = true_afs, type = "scatter", mode = "markers", themes="Catherine", name = "Inferred (Genotyped): r2 ") %>%
    layout(showlegend = T, 
           xaxis = list(title = 'True Allele Frequencies'), 
           yaxis = list(title = 'Inferred Allele Frequencies'),
           title = "Allele Frequency Inference")
  all_afs <- all_afs %>% add_trace(x =true_afs,  y = reference_af, name = 'Reference (Genotyped): r2', mode = 'markers')
  all_afs <- all_afs %>% add_trace(x =true_imputed_af,  y = inferred_af_imputed, mode = 'markers', name = 'Inferred (Imputed): r2')
  all_afs <- all_afs %>% add_trace(x =true_imputed_af,  y = reference_imputed_af, mode = 'markers', name = 'Reference (Imputed): r2')
  
  #Get the LD diagnostics, nowe we split up the LD into four groups (eg: typed-untyped/typed-typed/untyped-untyped/typed-untyped)
  
  #Get the reference LD
  reference_LD = LD_Matrix(reference_haplotypes_all)
  #Get the true LD, read in the GWAS haplotypes
  case_GWAS_haplotypes    = t(fread("hapgen2_sim_data/HAPGEN2_EUR_simulated.cases.haps", header = F))
  control_GWAS_haplotypes = t(fread("hapgen2_sim_data/HAPGEN2_EUR_simulated.controls.haps", header = F))
  #Cbind them together
  GWAS_haplotypes         = rbind(case_GWAS_haplotypes,control_GWAS_haplotypes)
  #Read in the GWAS legend
  GWAS_legend       = fread("hapgen2_sim_data/HAPGEN2_EUR_simulated.legend", header = T)
  GWAS_SNPs_index   = which(GWAS_legend$rs %in% reference_legend_all$id)
  GWAS_haplotypes   = GWAS_haplotypes[,c(GWAS_SNPs_index)]
  #Get the GWAS LD
  GWAS_LD = LD_Matrix(GWAS_haplotypes)
  #Get the Inferred LD (across all SNPs)
  inferred_LD = (cov.wt(x = reference_haplotypes_all, wt = wts, cor = T))$cor
  
  #Typed-Typed
  true_LD_typed_typed      = GWAS_LD[-c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  inferred_LD_typed_typed  = inferred_LD[-c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  reference_LD_typed_typed = reference_LD[-c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  #Get indexes of sampled typed-typed
  typed_indexes_subsample = sample(1:nrow(true_LD_typed_typed), size = 500)
  plot(true_LD_typed_typed[typed_indexes_subsample,typed_indexes_subsample],inferred_LD_typed_typed[typed_indexes_subsample,typed_indexes_subsample])
  
  #Typed-Untyped
  true_LD_typed_untyped      = GWAS_LD[c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  inferred_LD_typed_untyped  = inferred_LD[c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  reference_LD_typed_untyped = reference_LD[c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  #Untyped-Untyped
  true_LD_untyped_untyped      = GWAS_LD[-c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  inferred_LD_untyped_untyped  = inferred_LD[-c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  reference_LD_untyped_untyped = reference_LD[-c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  #Untyped-Typed
  true_LD_untyped_typed      = GWAS_LD[-c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  inferred_LD_untyped_typed  = inferred_LD[-c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  reference_LD_untyped_typed = reference_LD[-c(snp_ids_imputed_index),-c(snp_ids_imputed_index)]
  
  #Plot the AF
  return(plotly_build(all_afs))
}

evaluate_InferLD_sumstat_accuracy = function(){
  imputed_snps = fread("imputation_results/imputation_test")
  #Compare this to the true data
  #Load in the unmasked full dataset
  full_sumstats = fread("hapgen2_sim_data/full_sumstats_munged.sumstats", header = T)
  #Load in the imputed stuff from fizie
  #Merge together the data
  comparison_set = inner_join(imputed_snps,full_sumstats, by = "SNP")
  
  #PLot the r2 and graph for the full dataset
  r2_correlation_full = signif(summary(lm(comparison_set$Z.x~comparison_set$Z.y))$r.squared,4)
  #Plot and mark the r2 correlation)
  all_correlation = plot_ly(x = comparison_set$Z.x, y = comparison_set$Z.y, type = "scatter", mode = "markers", themes="Catherine") %>%
    layout(showlegend = F, annotations = list(text = sprintf("r2 = %s",r2_correlation_full), showarrow = F, x= 2,
                                              y= -2,
                                              xref = "x",
                                              yref = "y",
                                              showarrow = T,
                                              ax = 20,
                                              ay = -40), 
           xaxis = list(title = 'Imputed Zscores'), 
           yaxis = list(title = 'True (Genotyped) Zscores'),
           title = "Correlation True to Imputed Zscores (ALL)")
  
  #PLot the r2 and graph for the common variants
  snptest_results = fread("hapgen2_sim_data/full_sumstats_mungeable", header = T)
  comparison_set_common = inner_join(comparison_set,snptest_results, by = "BP") %>% filter(all_maf > 0.05)
  r2_correlation_common = signif(summary(lm(comparison_set_common$Z.x~comparison_set_common$Z.y))$r.squared,4)
  common_correlation = plot_ly(x = comparison_set_common$Z.x, y = comparison_set_common$Z.y, type = "scatter", mode = "markers", themes="Catherine") %>%
    layout(showlegend = F, annotations = list(text = sprintf("r2 = %s",r2_correlation_common), showarrow = F, x= 2,
                                              y= -2,
                                              xref = "x",
                                              yref = "y",
                                              showarrow = T,
                                              ax = 20,
                                              ay = -40),  
           xaxis = list(title = 'Imputed Zscores'), 
           yaxis = list(title = 'True (Genotyped) Zscores'),
           title = "Correlation True to Imputed Zscores (Common)")
  comparison_set_rare = inner_join(comparison_set,snptest_results, by = "BP") %>% filter(all_maf < 0.05)
  r2_correlation_rare = signif(summary(lm(comparison_set_rare$Z.x~comparison_set_rare$Z.y))$r.squared,4)
  rare_correlation = plot_ly(x = comparison_set_rare$Z.x, y = comparison_set_rare$Z.y, type = "scatter", mode = "markers", themes="Catherine") %>%
    layout(showlegend = F, annotations = list(text = sprintf("r2 = %s",r2_correlation_rare), showarrow = F, x= 2,
                                              y= -2,
                                              xref = "x",
                                              yref = "y",
                                              showarrow = T,
                                              ax = 20,
                                              ay = -40), 
           xaxis = list(title = 'Imputed Zscores'), 
           yaxis = list(title = 'True (Genotyped) Zscores'),
           title = "Correlation True to Imputed Zscores (Rare)")
  correlation = subplot(all_correlation,rare_correlation,common_correlation, shareX = T, shareY = T) %>% 
    layout(title = list(text = "Infer LD Sumstat Correlation"), xaxis = list(title = "Imputed Zscores"),yaxis = list(title = "True Zscores"))
  return(correlation)
}


  

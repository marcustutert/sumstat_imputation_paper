#Test running InferLD
#source("simulate_GWAS/simulation_functions.R")
set.seed(234)
x = run_InferLD(munged_genotyped_summary_statistics = "hapgen2_sim_data/genotyped_sumstats_munged",
            reference_haplotypes                    = "hapgen2_sim_data/reference_panel.haps",
            reference_legend                        = "hapgen2_sim_data/reference_panel.legend",
            case_control_constant                   = 500)
x = readRDS("inference_results/inference.RData")
plot(x$log_likelihood)

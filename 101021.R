#MWE for InferLD precision
#Build the
install.packages("InferLD_0.1.0.tar.gz", repos = NULL)
library(InferLD)

run_InferLD(reference_haplotypes = 'reference_haplotypes_filtered')
evaluate_InferLD_sumstat_accuracy()

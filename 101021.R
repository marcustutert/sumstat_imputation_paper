#MWE for InferLD precision
#Build the
install.packages("InferLD_0.1.0.tar.gz", repos = NULL)
library(InferLD)

run_InferLD()
evaluate_IMPG_accuracy()
evaluate_InferLD_sumstat_accuracy()
inference_diagnostics()

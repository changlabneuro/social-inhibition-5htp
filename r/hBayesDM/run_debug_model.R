library(hBayesDM)

data_root <- 'D:\\data\\hwwa\\public'
src_data_file <- 'hBayesDM-saline-nogo_trial.txt'
# src_data_file <- 'hBayesDM-5-htp.txt'
# src_data_file <- 'hBayesDM-saline.txt'
src_data_file <- file.path(data_root, src_data_file)

use_model_regressor <- FALSE

src_data <- read.csv(src_data_file)
gng_model <- gng_m1(data=src_data, niter=2000, nwarmup=1000, 
                    nchain=4, ncore=4, modelRegressor=use_model_regressor)
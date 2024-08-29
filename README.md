# RefinedEL

Source code for manuscript "Empirical Likelihood Based Inference for Functional Mean Models Accounting for Within-Subject Correlation"
by Xiang Wang and Honglang Wang

Xiang Wang and Honglang Wang developed this source code. For questions, comments or remarks about it, please contact Xiang Wang (xiangwangphd@outlook.com)


## This folder contains: 
- Functions.R: contains the main functions for our analysis
- real_data_application: contains the script files: ADNI_AD_data_application.R, ADNI_CN_data_application.R, ADNI_EMCI_data_application.R, ADNI_LMCI_data_application.R, ADNI_SMC_data_application.R, cd4_data_application.R, and pm25_data_application.R that has to be run to reproduce the results of the real data analysis, and the data ADNIMERGE.csv, longiCD4.txt, and pm2.5.csv.
- simulation_program.R: R code to reproduce the result of the simulation study. 



## The code has been executed using
- R version 4.2.1 (2022-06-23)
- Platform: x86_64-pc-linux-gnu (64-bit)
- Running under: Red Hat Enterprise Linux 8.9 (Ootpa)

- Matrix products: default
- BLAS:   /geode2/soft/hps/rhel8/r/gnu/4.2.1_X11/lib64/R/lib/libRblas.so
- LAPACK: /geode2/soft/hps/rhel8/r/gnu/4.2.1_X11/lib64/R/lib/libRlapack.so

- locale:
```r
  [1] LC_CTYPE=en_US.UTF-8        LC_NUMERIC=C
  [3] LC_TIME=en_US.UTF-8         LC_COLLATE=en_US.UTF-8
  [5] LC_MONETARY=en_US.UTF-8     LC_MESSAGES=en_US.UTF-8
  [7] LC_PAPER=en_US.UTF-8        LC_NAME=C
  [9] LC_ADDRESS=C                LC_TELEPHONE=C
  [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
```

- attached base packages:
```r
[1] splines   parallel  stats     graphics  grDevices utils     datasets
[8] methods   base
```
- other attached packages:
```r
 [1] xtable_1.8-4      fdapace_0.5.9     doParallel_1.0.17 iterators_1.0.14
 [5] foreach_1.5.2     emplik_1.3-1      stringr_1.5.0     locpol_0.8.0
 [9] PLRModels_1.4     face_0.1-7
```
- loaded via a namespace (and not attached):
```r
 [1] Rcpp_1.0.10         lattice_0.20-45     deldir_1.0-6
 [4] png_0.1-8           digest_0.6.29       utf8_1.2.2
 [7] R6_2.5.1            backports_1.4.1     MatrixModels_0.5-3
[10] pracma_2.4.2        ggplot2_3.5.0       pillar_1.9.0
[13] rlang_1.1.3         data.table_1.14.6   rstudioapi_0.14
[16] SparseM_1.81        rpart_4.1.16        Matrix_1.6-5
[19] checkmate_2.1.0     foreign_0.8-82      htmlwidgets_1.6.2
[22] munsell_0.5.0       numDeriv_2016.8-1.1 compiler_4.2.1
[25] xfun_0.32           pkgconfig_2.0.3     base64enc_0.1-3
[28] mgcv_1.8-40         htmltools_0.5.6     nnet_7.3-17
[31] tidyselect_1.2.0    tibble_3.2.1        gridExtra_2.3
[34] htmlTable_2.4.1     matrixcalc_1.0-6    Hmisc_4.7-2
[37] codetools_0.2-18    fansi_1.0.3         dplyr_1.1.2
[40] MASS_7.3-57         grid_4.2.1          nlme_3.1-157
[43] gtable_0.3.1        lifecycle_1.0.3     magrittr_2.0.3
[46] scales_1.3.0        cli_3.6.2           stringi_1.7.8
[49] latticeExtra_0.6-30 generics_0.1.3      vctrs_0.6.3
[52] Formula_1.2-4       RColorBrewer_1.1-3  tools_4.2.1
[55] interp_1.1-3        glue_1.6.2          jpeg_0.1-10
[58] fastmap_1.1.0       survival_3.5-7      colorspace_2.0-3
[61] cluster_2.1.3       knitr_1.40          quantreg_5.94
```


## Note that:
- our code for the simulations and the analysis of the real data example was run in parallel on Linux and Windows by packages "doParallel_1.0.17" and "foreach_1.5.2" with the number of cores equal to 20 (Line 25 of simulation_program.R and line 27 in data application R script files).
- Figures and tables of our simulation study and real data applications can be reproduced without running all our analyses. To do so, run the simulation_result_Estimator_CI_Test_TrueCov.R and simulation_result_Estimator_CI_Test_EstimatedCov.R, and real_data_application_results.R in figures folder. 


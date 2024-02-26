## ---- eval=F, include = FALSE-------------------------------------------------
#  knitr::opts_chunk$set(
#    collapse = TRUE,
#    comment = "#>"
#  )

## ----eval=F, setup------------------------------------------------------------
#  library(Unico)
#  library(matrixStats)
#  
#  #For visualization in this vignette
#  install.packages(c("ggplot2","ggpubr","hexbin","egg"))
#  source("https://github.com/cozygene/Unico/raw/main/vignettes/vignetts.utils.r")

## ----eval=F-------------------------------------------------------------------
#  data_path <- "./"
#  if(!file.exists(file.path(data_path, "pbmc.rds"))){
#  	download.file("https://github.com/cozygene/Unico/raw/main/vignettes/pbmc.rds", file.path(data_path,"pbmc.rds"))
#  }
#  sim.data = readRDS(file.path(data_path,"pbmc.rds"))

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  unico.res = list()
#  
#  #parameter learning
#  unico.res$params.hat <- Unico(sim.data$X, sim.data$W, C1 = NULL, C2 = NULL, parallel = TRUE)
#  
#  #tensor
#  unico.res$Z.hat = tensor(sim.data$X, W = sim.data$W, C1 = NULL, C2 = NULL,
#                           unico.res$params.hat, parallel = FALSE)
#  

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  # evaluate tensor performance on features with variation
#  unico.res$Z.corrs = calc_Z_corrs(Z.true = sim.data$Z.scale,
#                                   Z.hat  = unico.res$Z.hat,
#                                   eval.feature.source = sim.data$variable.feature.source)
#  
#  colMedians(unico.res$Z.corrs)

## ----eval=F, echo = T, results = 'hide', fig.show='hide'----------------------
#  low.gene = sample(rownames(sim.data$params$entropies[sim.data$params$entropies < quantile(sim.data$params$entropies, 0.25), ,drop= F]), 1)
#  
#  plot(sim.data$Z[3,low.gene, ], unico.res$Z.hat[3,low.gene, ]) + title (low.gene)

## ----echo = FALSE, out.width = "400px"----------------------------------------
knitr:: include_graphics("tensor.cor.png")

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  data_path <- "./"
#  if(!file.exists(file.path(data_path, "liu.rds"))){
#  	download.file("https://github.com/cozygene/Unico/raw/main/vignettes/liu.rds", file.path(data_path,"liu.rds"))
#  }
#  if(!file.exists(file.path(data_path, "hannum.rds"))){
#  	download.file("https://github.com/cozygene/Unico/raw/main/vignettes/hannum.rds", file.path(data_path,"hannum.rds"))
#  }
#  liu    = readRDS(file.path(data_path,"liu.rds"))
#  hannum = readRDS(file.path(data_path,"hannum.rds"))

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  source.ids   = colnames(liu$W)
#  n = ncol(liu$X)
#  m = nrow(liu$X)
#  k = ncol(liu$W)

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  unico.liu = list()
#  unico.liu$params.hat = Unico(X  = liu$X, W = liu$W,
#                               C1 = liu$cov[, c("age", "sex", "disease","smoking")],
#                               C2 = liu$ctrl_pcs)
#  
#  unico.liu$params.hat = association_parametric(X = liu$X, unico.liu$params.hat)

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  unico.hannum = list()
#  unico.hannum$params.hat = Unico(X  = hannum$X, W = hannum$W,
#                                  C1 = hannum$cov[, c("age", "sex", "ethnicity")],
#                                  C2 = cbind(hannum$ctrl_pcs, hannum$cov[,"plate", drop = F]))
#  
#  unico.hannum$params.hat = association_parametric(X = hannum$X, unico.hannum$params.hat)

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  liu.marg.pvals    = unico.liu$params.hat$parametric$gammas_hat_pvals[, paste0(source.ids, ".age")]
#  hannum.marg.pvals = unico.hannum$params.hat$parametric$gammas_hat_pvals[, paste0(source.ids, ".age")]
#  print(sum(liu.marg.pvals < 0.05/(m*k)))
#  print(hannum.marg.pvals[liu.marg.pvals < 0.05/(m*k)])

## ----eval=F, echo = T, results = 'hide', fig.height=4.5, fig.width=22---------
#  qq_age_g = plot_qq(pvals_mat = liu.marg.pvals,
#                     labels = source.ids,
#                     ggarrange.nrow = 1, ggarrange.ncol = k,
#                     alpha = 0.5, text.size = 20,
#                     title = "Parametric association testing (age) at cell-type resolution")
#  qq_age_g

## ----echo = FALSE, out.width = "675px"----------------------------------------
knitr:: include_graphics("qq_age.png")

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  liu.joint.pvals    = unico.liu$params.hat$parametric$gammas_hat_pvals.joint[, "age"]
#  hannum.joint.pvals = unico.hannum$params.hat$parametric$gammas_hat_pvals.joint[, "age"]
#  print(sum(liu.joint.pvals < 0.05/m))
#  print(sum(hannum.joint.pvals[liu.joint.pvals < 0.05/m] < (0.05/sum(liu.joint.pvals < 0.05/m))))

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  unico.liu$params.hat = association_asymptotic(X = liu$X, unico.liu$params.hat)
#  liu.marg.pvals.asym  = unico.liu$params.hat$asymptotic$gammas_hat_pvals[, paste0(source.ids, ".age")]

## ----eval=F, echo = T, results = 'hide', fig.height=4.5, fig.width=22---------
#  qq_compare_g = plot_qq_compare(pvals_mat1 = liu.marg.pvals,
#                                 pvals_mat2 = liu.marg.pvals.asym,
#                                 labels = source.ids,
#                                 ggarrange.nrow = 1, ggarrange.ncol = k,
#                                 alpha = 0.05, text.size = 20,
#                                 xlab = "Parametric", ylab = "Asymptotic",
#                                 title = "Parametric vs Asymptotic pvals")
#  qq_compare_g

## ----echo = FALSE, out.width = "675px"----------------------------------------
knitr:: include_graphics("qq_compare.png")

## ----eval=F, echo = T, results = 'hide'---------------------------------------
#  C1.shuffle = hannum$cov[, c("age", "sex", "ethnicity")]
#  C1.shuffle[, "age"] = hannum$cov[sample(1:nrow(hannum$cov)), "age"]
#  
#  unico.hannum.shuffle = list()
#  unico.hannum.shuffle$params.hat = Unico(X  = hannum$X, W = hannum$W,
#                                          C1 = C1.shuffle,
#                                          C2 = cbind(hannum$ctrl_pcs, hannum$cov[,"plate", drop = F]))
#  unico.hannum.shuffle$params.hat = association_asymptotic(X = hannum$X, unico.hannum.shuffle$params.hat)
#  unico.marg.pvals.asym.calib = unico.hannum.shuffle$params.hat$asymptotic$gammas_hat_pvals[, paste0(source.ids, ".age")]

## ----eval=F, echo = T, results = 'hide', fig.height=4.5, fig.width=22---------
#  qq_calib_g = plot_qq(pvals_mat = unico.marg.pvals.asym.calib,
#                       labels = source.ids,
#                       ggarrange.nrow = 1, ggarrange.ncol = k,
#                       alpha = 0.5, text.size = 20,
#                       title = "Calibration of asymptotic association testing")
#  qq_calib_g
#  

## ----echo = FALSE, out.width = "675px"----------------------------------------
knitr:: include_graphics("qq_calib.png")


## ----setup, class.source='fold-hide', warning=FALSE---------------------------
#install.packages("knitr")
library(knitr)
opts_chunk$set(echo = TRUE)

# To ensure reproducibility
set.seed(12)

## ----packages_installation, eval=FALSE, class.source='fold-hide'--------------
#  pkgs <- c("rchemo")
#  sapply(pkgs, function(x) {
#    if (!requireNamespace(x, quietly = TRUE)) {
#      install.packages(x)
#    }
#  })

## ----packages_load, warning=FALSE, message=FALSE, class.source='fold-hide'----
library(rchemo)  # to load rchemo

## ----import_Zhang2023---------------------------------------------------------
data(Zhang2023, package = "rchemo")

## ----check_dims, class.source='fold-hide'-------------------------------------
# Check dimension
BlockNames <- names(Zhang2023)
nbrBlocs <- length(BlockNames)
dims <- lapply(X=Zhang2023[BlockNames], FUN=dim)
names(dims) <- BlockNames
dims

# Remove unuseful object for the next steps
rm(nbrBlocs, dims)

## ----check_orders_and_names, class.source='fold-hide'-------------------------
# Check rows names in any order
row_names <- lapply(X=Zhang2023[BlockNames], FUN=rownames)
rns <- do.call(cbind, row_names)
rns.unique <- apply(rns, 1, function(x) length(unique(x)))
if (max(rns.unique) > 1) {
  stop("Rows names are not identical between blocks.")
}

# Check order of samples
check_row_names <- all(sapply(X=row_names, FUN=identical, y = row_names[[1]]))
if (!check_row_names && max(rns.unique) == 1) {
  print("Rows names are not in the same order for all blocks.")
}

# Remove unuseful object for the next steps
rm(row_names, rns, rns.unique, check_row_names)

## ----change_data_format-------------------------------------------------------
GCTOF <- Zhang2023$GCTOF
HILICPOS <- Zhang2023$HILICPOS
sample_metadata <- Zhang2023$metadata
rm(Zhang2023)

## ----PCA_separate_datablocks, class.source='fold-hide'------------------------
# GCTOF with unit variance scaling
pcaGCTOF <- pcanipalsna(X = scale(GCTOF[,1:dim(GCTOF)[2]]), nlv = nrow(GCTOF), 
                         gs = TRUE,
                         tol = .Machine$double.eps^0.5, maxit = 200)

## diagram of explained variance
barplot(summary(pcaGCTOF,X = scale(GCTOF[,1:dim(GCTOF)[2]]))$explvar$pvar * 100, names.arg = 1:nrow(GCTOF), main = "diagram of explained variance - PCA GCTOF")

## score plot
plotxy(X= pcaGCTOF$T, group = sample_metadata$Gender, 
       asp = 0, col = 3:4, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = FALSE,
       legend = TRUE, main = "components - PCA GCTOF", ncol = 1,
       pch=16)
## loading plot
plotxy(X= pcaGCTOF$P, group = NULL, 
       asp = 0, col = NULL, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = TRUE,
       legend = FALSE, main = "loadings - PCA GCTOF", ncol = 1,
       cex=0.8) 

#HILICPOS with unit variance scaling
pcaHILICPOS <- pcanipalsna(X = scale(HILICPOS[,1:dim(HILICPOS)[2]]), nlv = nrow(HILICPOS), 
                         gs = TRUE,
                         tol = .Machine$double.eps^0.5, maxit = 200)

## diagram of explained variance
barplot(summary(pcaHILICPOS,X = scale(HILICPOS[,1:dim(HILICPOS)[2]]))$explvar$pvar * 100, names.arg = 1:nrow(HILICPOS), main = "diagram of explained variance - PCA HILIC POS")

## score plot
plotxy(X= pcaHILICPOS$T, group = sample_metadata$Gender, 
       asp = 0, col = 3:4, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = FALSE,
       legend = TRUE, main = "components - PCA HILIC POS", ncol = 1,
       pch=16)

## loading plot
plotxy(X= pcaHILICPOS$P, group = NULL, 
       asp = 0, col = NULL, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = TRUE,
       legend = FALSE, main = "loadings - PCA HILIC POS", ncol = 1,
       cex=0.8)


# Remove unuseful object for the next steps
rm(pcaGCTOF, pcaHILICPOS)

## ----run_nlvtest--------------------------------------------------------------
nlvtestmbplsda <- mbplsr_mbplsda_allsteps(Xlist = list(GCTOF = GCTOF[,1:dim(GCTOF)[2]], 
                                                       HILICPOS = HILICPOS[,1:dim(HILICPOS)[2]]),
                                          Xnames = c("GCTOF", "HILICPOS"), 
                                          Xscaling = c("none","pareto","sd")[3], 
                                          Y = sample_metadata[,"Gender", drop=FALSE], 
                                          Yscaling = c("none","pareto","sd")[1], 
                                          weights = NULL,
                                          newXlist = NULL, newXnames = NULL,
                                          method = c("mbplsrda","mbplslda","mbplsqda")[2],
                                          prior = c("unif", "prop")[1],
                                          step = c("nlvtest","permutation","model","prediction")[1],
                                          nlv = 4, 
                                          #modeloutput = c("scores","loadings","coef","vip"), 
                                          cvmethod = c("kfolds","loo")[1], 
                                          nbrep = 10, 
                                          seed = 123, 
                                          samplingk = NULL, 
                                          nfolds = 3, 
                                          #npermut = 30, 
                                          
                                          criterion = c("err","rmse")[1], 
                                          selection = c("localmin","globalmin","1std")[1],
                                          
                                          outputfilename = NULL)

nlvtestmbplsda

nlvoptmbplsda <- nlvtestmbplsda[nlvtestmbplsda$optimum==1,"nblv"] # to obtain the optimal number of LV.

# to plot the results of the cross-validation
plot(nlvtestmbplsda$nblv, nlvtestmbplsda$err_mean, xlab = "number of LV", ylab = "CV classification error rate", pch = 16, ylim = c(0,0.6))
segments(nlvtestmbplsda$nblv,nlvtestmbplsda$err_mean-nlvtestmbplsda$err_sd,nlvtestmbplsda$nblv,nlvtestmbplsda$err_mean+nlvtestmbplsda$err_sd)
segments(nlvtestmbplsda$nblv-0.1,nlvtestmbplsda$err_mean-nlvtestmbplsda$err_sd,nlvtestmbplsda$nblv+0.1,nlvtestmbplsda$err_mean-nlvtestmbplsda$err_sd)
segments(nlvtestmbplsda$nblv-0.1,nlvtestmbplsda$err_mean+nlvtestmbplsda$err_sd,nlvtestmbplsda$nblv+0.1,nlvtestmbplsda$err_mean+nlvtestmbplsda$err_sd)

# Remove unuseful object for the next steps
rm(nlvtestmbplsda)

## ----run_permutation----------------------------------------------------------
permutmbplsda <- mbplsr_mbplsda_allsteps(Xlist = list(GCTOF = GCTOF[,1:dim(GCTOF)[2]], 
                                                      HILICPOS = HILICPOS[,1:dim(HILICPOS)[2]]), 
                                         Xnames = c("GCTOF", "HILICPOS"), 
                                         Xscaling = c("none","pareto","sd")[3], 
                                         Y = sample_metadata[,"Gender",drop=FALSE], 
                                         Yscaling = c("none","pareto","sd")[1], weights = NULL,
                                         newXlist = NULL, newXnames = NULL,
                                         
                                         method = c("mbplsrda","mbplslda","mbplsqda")[2],
                                         prior = c("unif", "prop")[1],
                                         
                                         step = c("nlvtest","permutation","model","prediction")[2],
                                         nlv = nlvoptmbplsda, 
                                         modeloutput = c("scores","loadings","coef","vip"), 
                                         
                                         cvmethod = c("kfolds","loo")[1], 
                                         nbrep = 10, 
                                         seed = 123, 
                                         samplingk = NULL, 
                                         nfolds = 3, 
                                         npermut = 10, 
                                         
                                         criterion = c("err","rmse")[1], 
                                         # selection = c("localmin","globalmin","1std")[1],
                                         
                                         import = c("R","ChemFlow","W4M")[1],
                                         outputfilename = NULL)

#plot of the results
plot(permutmbplsda, pch = 16, ylab = "CV classification error rate", xlab = "dyssimilarity Y-Ypermuted")

# Remove unuseful object for the next steps
rm(permutmbplsda)

## ----run_model----------------------------------------------------------------
modelmbplsda <- mbplsr_mbplsda_allsteps(Xlist = list(GCTOF = GCTOF[,1:dim(GCTOF)[2]], 
                                                     HILICPOS = HILICPOS[,1:dim(HILICPOS)[2]]), 
                                        Xnames = c("GCTOF", "HILICPOS"), 
                                        Xscaling = c("none","pareto","sd")[3], 
                                        Y = sample_metadata[,"Gender",drop=FALSE], 
                                        Yscaling = c("none","pareto","sd")[1], 
                                        weights = NULL,
                                        newXlist = NULL, newXnames = NULL,
                                        
                                        method = c("mbplsrda","mbplslda","mbplsqda")[2],
                                        prior = c("unif", "prop")[1],
                                        
                                        step = c("nlvtest","permutation","model","prediction")[3],
                                        nlv = nlvoptmbplsda, 
                                        modeloutput = c("scores","loadings","coef","vip"), 
                                        
                                        cvmethod = c("kfolds","loo")[1], 
                                        # nbrep = 30, 
                                        # seed = 123, 
                                        # samplingk = NULL, 
                                        # nfolds = 5, 
                                        # npermut = 30, 
                                        
                                        # criterion = c("err","rmse")[1], 
                                        # selection = c("localmin","globalmin","1std")[1],
                                        
                                        import = c("R","ChemFlow","W4M")[1],
                                        outputfilename = NULL)

# score plot
plotxy(X= modelmbplsda$scores, group = sample_metadata$Gender, 
       asp = 0, col = 3:4, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = FALSE,
       legend = TRUE, main = "scores - MB PLS DA", ncol = 1,
       pch=16)

# loading plot
plotxy(X=  modelmbplsda$loadings, group = substr(rownames(modelmbplsda$loadings),1,6), 
       asp = 0, col = NULL, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = FALSE,
       legend = TRUE, main = "loadings - MB PLS DA", ncol = 1,
       cex=0.8, pch = 16)

# VIP curve
plot(modelmbplsda$vip[order(modelmbplsda$vip[,nlvoptmbplsda], decreasing = TRUE),nlvoptmbplsda], pch = 16,cex = 0.8,
     col = as.numeric(as.factor(substr(rownames(modelmbplsda$vip[order(modelmbplsda$vip[,nlvoptmbplsda], decreasing = TRUE),nlvoptmbplsda, drop=FALSE]),1,6))), ylab = "VIP value",
     main = "VIP curve MB PLS DA")
legend("topright", legend = c("GCTOF", "HILICPOS"), pch = 16, col = 1:2)

# Remove unuseful object for the next steps
rm(modelmbplsda)

## ----run_prediction-----------------------------------------------------------
predmbplsda <- mbplsr_mbplsda_allsteps(Xlist = list(GCTOF = GCTOF[,1:dim(GCTOF)[2]], 
                                                    HILICPOS = HILICPOS[,1:dim(HILICPOS)[2]]), 
                                        Xnames = c("GCTOF", "HILICPOS"), 
                                        Xscaling = c("none","pareto","sd")[3], 
                                        Y = sample_metadata[,"Gender",drop=FALSE], 
                                       Yscaling = c("none","pareto","sd")[1], 
                                       weights = NULL,
                                        newXlist = NULL, newXnames = NULL,
                                        
                                        method = c("mbplsrda","mbplslda","mbplsqda")[2],
                                        prior = c("unif", "prop")[1],
                                        
                                        step = c("nlvtest","permutation","model","prediction")[4],
                                        nlv = nlvoptmbplsda, 
                                        # modeloutput = c("scores","loadings","coef","vip"), 
                                        # 
                                        # cvmethod = c("kfolds","loo")[1], 
                                        # nbrep = 30, 
                                        # seed = 123, 
                                        # samplingk = NULL, 
                                        # nfolds = 5, 
                                        # npermut = 30, 
                                        # 
                                        # criterion = c("err","rmse")[1], 
                                        # selection = c("localmin","globalmin","1std")[1],
                                        
                                        import = c("R","ChemFlow","W4M")[1],
                                        outputfilename = NULL)

predmbplsda
# Remove unuseful object for the next steps
rm(predmbplsda, nlvoptmbplsda)

## ----reproducibility----------------------------------------------------------
sessionInfo()


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
HILICNEG <- Zhang2023$HILICNEG
sample_metadata <- Zhang2023$metadata
rm(Zhang2023)

## ----PCA_separate_datablocks, class.source='fold-hide'------------------------
# GCTOF with unit variance scaling
pcaGCTOF <- pcanipalsna(X = scale(GCTOF[,1:dim(GCTOF)[2]]), nlv = nrow(GCTOF), 
                         gs = TRUE,
                         tol = .Machine$double.eps^0.5, maxit = 200)

## diagram of explained variance
barplot(summary(pcaGCTOF,X = scale(GCTOF[,1:dim(GCTOF)[2]]))$explvar$pvar * 100, names.arg = 1:nrow(summary(pcaGCTOF,X = scale(GCTOF[,1:dim(GCTOF)[2]]))$explvar), main = "diagram of explained variance - PCA GCTOF")

## score plot
plotxy(X= pcaGCTOF$T, group = sample_metadata$Gender, 
       asp = 0, col = 3:4, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = FALSE,
       legend = TRUE, main = "scores - PCA GCTOF", ncol = 1,
       pch=16)
## loading plot
plotxy(X= pcaGCTOF$P, group = NULL, 
       asp = 0, col = NULL, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = TRUE,
       legend = FALSE, main = "loadings - PCA GCTOF", ncol = 1,
       cex=0.8) 

#HILICNEG with unit variance scaling
pcaHILICNEG <- pcanipalsna(X = scale(HILICNEG[,1:dim(HILICNEG)[2]]), nlv = nrow(HILICNEG), 
                         gs = TRUE,
                         tol = .Machine$double.eps^0.5, maxit = 200)

## diagram of explained variance
barplot(summary(pcaHILICNEG,X = scale(HILICNEG[,1:dim(HILICNEG)[2]]))$explvar$pvar * 100, names.arg = 1:nrow(summary(pcaHILICNEG,X = scale(HILICNEG[,1:dim(HILICNEG)[2]]))$explvar), main = "diagram of explained variance - PCA HILIC NEG")

## score plot
plotxy(X= pcaHILICNEG$T, group = sample_metadata$Gender, 
       asp = 0, col = 3:4, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = FALSE,
       legend = TRUE, main = "scores - PCA HILIC NEG", ncol = 1,
       pch=16)

## loading plot
plotxy(X= pcaHILICNEG$P, group = NULL, 
       asp = 0, col = NULL, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = TRUE,
       legend = FALSE, main = "loadings - PCA HILIC NEG", ncol = 1,
       cex=0.8)


# Remove unuseful object for the next steps
rm(pcaGCTOF, pcaHILICNEG)

## ----run_nlvtest--------------------------------------------------------------
nlvtestsoplsda <- soplsr_soplsda_allsteps(Xlist = list(GCTOF = GCTOF[,1:dim(GCTOF)[2]], 
                                                       HILICNEG = HILICNEG[,1:dim(HILICNEG)[2]]), 
                                          Xnames = c("GCTOF", "HILICNEG"), 
                                          Xscaling = c("none","pareto","sd")[3], 
                                          Y = sample_metadata[,"Gender",drop=FALSE], 
                                          Yscaling = c("none","pareto","sd")[1], weights = NULL,
                                          #newXlist = NULL, newXnames = NULL,
                                          
                                          method = c("soplsrda","soplslda","soplsqda")[2],
                                          prior = c("unif", "prop")[1],
                                          
                                          step = c("nlvtest","permutation","model","prediction")[1],
                                          # nlv = c(),
                                          nlvlist = list(0:3, 0:3), 
                                          # modeloutput = c("scores","loadings","coef","vip"), 
                                          
                                          cvmethod = c("kfolds","loo")[1], 
                                          nbrep = 10, 
                                          seed = 123, 
                                          samplingk = NULL, 
                                          nfolds = 3, 
                                          # npermut = 30, 
                                          
                                          optimisation = c("global","sequential")[1],
                                          criterion = c("err","rmse")[1], 
                                          selection = c("localmin","globalmin","1std")[1],
                                          
                                          #import = c("R","ChemFlow","W4M")[1],
                                          outputfilename = NULL)


# to plot the results of the cross-validation
plot(nlvtestsoplsda[,"nlvsum"],nlvtestsoplsda[,"mean"], type = "n", xlab = "nlv sum", ylab = "classification error rate", 
     main = "error rates obtained from the different numbers \n of latent variable combinaisons of the 2 datablocks")
text(x=nlvtestsoplsda[,"nlvsum"], y=nlvtestsoplsda[,"mean"], 
     labels=paste0(nlvtestsoplsda[,"Xlist1"],",",nlvtestsoplsda[,"Xlist2"]))

nlvoptsoplsda <- nlvtestsoplsda[which(nlvtestsoplsda[,"optimum"]==1),c("Xlist1","Xlist2")] # to obtain the optimal number of LV.

# Remove unuseful object for the next steps
rm(nlvtestsoplsda)

## ----run_permutation----------------------------------------------------------
permutsoplsda <- soplsr_soplsda_allsteps(Xlist = list(GCTOF = GCTOF[,1:dim(GCTOF)[2]], 
                                                      HILICNEG = HILICNEG[,1:dim(HILICNEG)[2]]), 
                                         Xnames = c("GCTOF", "HILICNEG"), 
                                         Xscaling = c("none","pareto","sd")[3], 
                                         Y = sample_metadata[,"Gender",drop=FALSE], 
                                         Yscaling = c("none","pareto","sd")[1], weights = NULL,
                                         #newXlist = NULL, newXnames = NULL,
                                         
                                         method = c("soplsrda","soplslda","soplsqda")[2],
                                         prior = c("unif", "prop")[1],
                                         
                                         step = c("nlvtest","permutation","model","prediction")[2],
                                         nlv = nlvoptsoplsda,
                                         #nlvlist = list(),
                                         #modeloutput = c("scores","loadings","coef","vip"), 
                                         
                                         cvmethod = c("kfolds","loo")[1], 
                                         nbrep = 10, 
                                         seed = 123, 
                                         samplingk = NULL, 
                                         nfolds = 3, 
                                         npermut = 10, 
                                         
                                         # optimisation = c("global","sequential")[1],
                                         criterion = c("err","rmse")[1], 
                                         # selection = c("localmin","globalmin","1std")[1],
                                         
                                         #import = c("R","ChemFlow","W4M")[1],
                                         outputfilename = NULL)

#plot of the results
plot(permutsoplsda, pch = 16, ylab = "CV classification error rate", xlab = "dyssimilarity Y-Ypermuted")

# Remove unuseful object for the next steps
rm(permutsoplsda)

## ----run_model----------------------------------------------------------------
modelsoplsda <- soplsr_soplsda_allsteps(Xlist = list(GCTOF = GCTOF[,1:dim(GCTOF)[2]], 
                                                     HILICNEG = HILICNEG[,1:dim(HILICNEG)[2]]), 
                                         Xnames = c("GCTOF", "HILICNEG"), 
                                         Xscaling = c("none","pareto","sd")[3], 
                                         Y = sample_metadata[,"Gender",drop=FALSE], 
                                         Yscaling = c("none","pareto","sd")[1], weights = NULL,
                                         #newXlist = NULL, newXnames = NULL,
                                         
                                         method = c("soplsrda","soplslda","soplsqda")[2],
                                         prior = c("unif", "prop")[1],
                                         
                                         step = c("nlvtest","permutation","model","prediction")[3],
                                         nlv = nlvoptsoplsda,
                                         #nlvlist = list(), 
                                         modeloutput = c("scores","loadings","coef","vip"), 
                                         
                                         # cvmethod = c("kfolds","loo")[1], 
                                         # nbrep = 30, 
                                         # seed = 123, 
                                         # samplingk = NULL, 
                                         # nfolds = 5, 
                                         # npermut = 30, 
                                         # 
                                         # criterion = c("err","rmse")[1], 
                                         # selection = c("localmin","globalmin","1std")[1],
                                         
                                         #import = c("R","ChemFlow","W4M")[1],
                                         outputfilename = NULL)
# score plot

plotxy(do.call("cbind",modelsoplsda$scores)[,1:2], group = sample_metadata[,"Gender"],
       asp = 0, col = 3:4, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = FALSE,
       legend = TRUE, main = "scores - SO PLS DA \n first model (GCTOF variables)", ncol = 1,
       pch=16,
       xlab="lv1 GCTOF",
       ylab="lv2 GCTOF")

plotxy(do.call("cbind",modelsoplsda$scores)[,3:4], group = sample_metadata[,"Gender"],
       asp = 0, col = 3:4, alpha.f = .8,
       zeroes = TRUE, circle = FALSE, ellipse = FALSE,
       labels = FALSE,
       legend = TRUE, main = "scores - SO PLS DA \n 2nd model (HILIC NEG variables)", ncol = 1,
       pch=16,
       xlab="lv1 HILICNEG",
       ylab="lv2 HILICNEG")

# loading plot
plotxy(modelsoplsda$loadings[[1]], pch = 16, xlab="lv1 GCTOF", ylab="lv2 GCTOF",
       main = "loadings - SO PLS DA \n first model (GCTOF variables)")

plotxy(modelsoplsda$loadings[[2]], pch = 16, xlab="lv1 HILIC NEG", ylab="lv2 HILIC NEG",
       main = "loadings - SO PLS DA \n 2nd model (HILIC NEG variables)")

# vip curve
plot(modelsoplsda$vip[[1]][order(modelsoplsda$vip[[1]][,2], decreasing = TRUE),1,drop=FALSE], pch = 16, ylab = "VIP", xlab="GCTOF variables", main = "VIP curve - SO PLS DA \n first model (GCTOF variables)")

plot(modelsoplsda$vip[[2]][order(modelsoplsda$vip[[2]][,2], decreasing = TRUE),1,drop=FALSE], pch = 16, ylab = "VIP", xlab="HILIC NEG variables", main = "VIP curve - SO PLS DA \n 2nd model (HILIC NEG variables)")

# Remove unuseful object for the next steps
rm(modelsoplsda)

## ----run_prediction-----------------------------------------------------------
predsoplsda <- soplsr_soplsda_allsteps(Xlist = list(GCTOF = GCTOF[,1:dim(GCTOF)[2]], 
                                                    HILICNEG = HILICNEG[,1:dim(HILICNEG)[2]]), 
                                        Xnames = c("GCTOF", "HILICNEG"), 
                                        Xscaling = c("none","pareto","sd")[3], 
                                        Y = sample_metadata[,"Gender",drop=FALSE], 
                                        Yscaling = c("none","pareto","sd")[1], 
                                        weights = NULL,
                                        newXlist = list(GCTOF = GCTOF[,1:dim(GCTOF)[2]], 
                                                        HILICNEG = HILICNEG[,1:dim(HILICNEG)[2]]),
                                       newXnames =  c("GCTOF", "HILICNEG"),
                                        
                                        method = c("soplsrda","soplslda","soplsqda")[2],
                                        prior = c("unif", "prop")[1],
                                        
                                        step = c("nlvtest","permutation","model","prediction")[4],
                                        nlv = nlvoptsoplsda, 
                                        # nlvlist = list(),
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
                                        # 
                                        #import = c("R","ChemFlow","W4M")[1],
                                        outputfilename = NULL)

predsoplsda

# Remove unuseful object for the next steps
rm(predsoplsda, nlvoptsoplsda)

## ----reproducibility----------------------------------------------------------
sessionInfo()


## ffbase Provides support for data.frame like objects that connect to
## files stored out of memory
library("ffbase")

## bigmemory provides the big.matrix object which links to matrices stored
## out of memory
library("bigmemory")

## biglm implements a bounded-memory fitter for glms
library("biglm")

## ggplot2 for figure plotting
library("ggplot2")

## sgd is for stochastic gradient descent methods
library("sgd")

## library required for streaming data
library("stream")

mydir <- "/Volumes/Jie Hao/BGSE/Big_Regression"
HiggsDir <- paste0(mydir,"/", "HIGGS","/", "HIGGSffdf")
load.ffdf(HiggsDir)

objectsDir <- paste0(mydir,"/","R_objects")

## Problems arose when I saved and loaded traindat_ffdf and testdat_ffdf
## so I subset them manually each time I load the data

traindat_ffdf <- subset(Higgs_ffdf, test == 0)
testdat_ffdf <- subset(Higgs_ffdf, test == 1)


## Higgs formulae
HiggsFormula1 <- formula(paste("signal ~",
                               paste0("feature", 1:21, collapse = " + "))) 
HiggsFormula2 <- formula(paste("signal ~",
                               paste0("HLfeature", 1:7, collapse = " + ")))
HiggsFormula3 <- formula(paste("signal ~",
                               paste0("feature", 1:21, collapse = " + "), "+",
                               paste0("HLfeature", 1:7, collapse = " + ")))

##------------------------------------------------------------------------------------------------------

## set directory to the directory with stored R objects
setwd(objectsDir)

## loading results

## comparisons of different methods for explicit sgd updates
results1 <- readRDS("results1") # comparison between bigglm, "adagrad", "one-dim" and "rmsprop" for HiggsFormula1
results2 <- readRDS("results2") # HiggsFormula2
results3 <- readRDS("results3") # HiggsFormula3

mod_bigglm1 <- readRDS("mod_bigglm1") # bigglm results identical to the ones presented in lecture
mod_bigglm2 <- readRDS("mod_bigglm2")
mod_bigglm3 <- readRDS("mod_bigglm3")

# first version of sgd chunking function
sgd_chunk <- function(x, chunk_size=NULL, formula, sgd_model="glm", sgd_modelcontrol=list(family=binomial(logit)),
                      sgd_method="sgd", lrvalue="adagrad", no_passes=3,...) {
  
  ## if no chunk size is chosen, we let the chunk function decide based on RAM availability
  higgs_chunks <- chunk.ffdf(x, from=1, to=nrow(x), by=chunk_size) 
  
  ## we loop through the chunks and apply the sgd function to it
  for (i in higgs_chunks) {
    sgd.theta <-sgd(formula, data = x[i,], model = "glm", model.control=sgd_modelcontrol, 
                    sgd.control = list(method=sgd_method, lr=lrvalue, npasses=no_passes, ...))
  }
  
  return(sgd.theta)
  
}

# load objects
a_1 <- readRDS("a_1")
a_2 <- readRDS("a_2")
a_3 <- readRDS("a_3")
a_4 <- readRDS("a_4")
experiment1 <- readRDS("experiment1") # comparison of coefficients

# experiment 1 for chunk size- I have saved these objects for reference
a_1 <- sgd_chunk(traindat_ffdf, formula=HiggsFormula1, lr="adagrad")
a_2 <- sgd_chunk(traindat_ffdf, chunk_size=1E6, formula=HiggsFormula1, lr="adagrad") #converged
a_3 <- sgd_chunk(traindat_ffdf, chunk_size=3E6, formula=HiggsFormula1, lr="adagrad")
a_4 <- sgd_chunk(traindat_ffdf, chunk_size=5E5, formula=HiggsFormula1, lr="adagrad")

## sample code for euclidean distance from in-memory solution
sqrt(t(as.vector(experiment1$`computer allocated`-results1$adagrad))%*%as.vector(experiment1$`computer allocated`-results1$adagrad))




##------------------------------------------------------------------------------------------------------




# revised sgd chunking function
sgd_chunk_revised <- function(x, chunk_size=NULL, formula, sgd_model="glm", sgd_modelcontrol=list(family=binomial(logit)),
                              sgd_method="sgd", lrvalue="adagrad", no_passes=3,...) {
  
  ## if no chunk size is chosen, we let the chunk function decide based on RAM availability
  higgs_chunks <- chunk.ffdf(x, from=1, to=nrow(x), by=chunk_size) 
  
  ## we loop through the chunks and apply the sgd function to it
  for (j in 1:no_passes) {
    for (i in higgs_chunks) {
      sgd.theta <-sgd(formula, data = x[i,], model = "glm", model.control=sgd_modelcontrol, 
                      sgd.control = list(method=sgd_method, lr=lrvalue, npasses=1, ...))
    }
  }
  
  return(sgd.theta)
  
}

## read objects
b_1 <- readRDS("b_1")
b_2 <- readRDS("b_2")
b_3 <- readRDS("b_3")
b_4 <- readRDS("b_4")
experiment2 <- readRDS("experiment2") # comparison of coefficients

## experiment 2 - I have included these objects for reference
b_1 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula1, lr="adagrad")
b_2 <- sgd_chunk_revised(traindat_ffdf, chunk_size=1E6, formula=HiggsFormula1, lr="adagrad") #converged
b_3 <- sgd_chunk_revised(traindat_ffdf, chunk_size=3E6, formula=HiggsFormula1, lr="adagrad")
b_4 <- sgd_chunk_revised(traindat_ffdf, chunk_size=5E5, formula=HiggsFormula1, lr="adagrad")


## load objects
morepasses <- readRDS("morepasses") # 6 passes with chunk size 3E6, failed to converge
morepasses2 <- readRDS("morepasses2") # 10 passes with chunk size 3E6, failed to converge
morepasses3 <- readRDS("morepasses3") # 6 passes with automatically allocated chunk size, failed to converge
morepasses4 <- readRDS("morepasses4") # 10 passes with automatically allocated chunk size, failed to converge


## with chunk size 3E6
morepasses <- sgd_chunk_revised(traindat_ffdf, chunk_size=3E6, formula=HiggsFormula1,no_passes=6)
morepasses2 <- sgd_chunk_revised(traindat_ffdf, chunk_size=3E6, formula=HiggsFormula1,no_passes=10)

## with computer allocated chunk size
morepasses3 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula1,no_passes=6)
morepasses4 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula1,no_passes=10)


## sample code for euclidean distance from in-memory solution 
sqrt(t(as.vector(experiment2$`computer allocated`-results1$adagrad))%*%as.vector(experiment2$`computer allocated`-results1$adagrad))


# load objects
mod_sgd1 <- readRDS("mod_sgd1") 
mod_sgd2 <- readRDS("mod_sgd2")
mod_sgd3 <-readRDS("mod_sgd3")


# choice of hyperparameters for function
mod_sgd1 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula1, lr="adagrad")
mod_sgd2 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula2, lr="adagrad")
mod_sgd3 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula3, lr="adagrad")


##------------------------------------------------------------------------------------------------------

## using data streaming function
sgd_stream <- function(x, formula, sgd_model = "glm", 
                       sgd_modelcontrol = list(family="binomial"), sgd_method="sgd", ...) {
  
  ## defensive programming for invalid function arguments
  
  if (!is.language(formula)){
    return("Formula is not valid")
  }
  
  if (!is.ffdf(x) & !is.data.frame(x)){
    return("Input data not valid")
  }
  
  # create the datastream to be used
  
  higgs_stream <- DSD_Memory(x, k = NA, loop = TRUE)
  
  ## stochastic gradient descent on input data according to input formula
  
  theta <-sgd(formula, data = get_points(higgs_stream,n=nrow(x)), model = "glm", 
              model.control = sgd_modelcontrol, sgd.control = list(method = sgd_method,...))
  
  return(theta)
  
}

##------------------------------------------------------------------------------------------------------

## Calculate the fitted values for the test data from each model, taken from slides
testdat_ffdf$fitted1 <- as.ff(predict(mod_bigglm1, testdat_ffdf,
                                      type = "response")[,1])
testdat_ffdf$fitted2 <- as.ff(predict(mod_bigglm2, testdat_ffdf,
                                      type = "response")[,1])
testdat_ffdf$fitted3 <- as.ff(predict(mod_bigglm3, testdat_ffdf,
                                      type = "response")[,1])
testdat_ffdf$fitted4 <- as.ff(as.numeric(NA), length = nrow(testdat_ffdf))
testdat_ffdf$fitted5 <- as.ff(as.numeric(NA), length = nrow(testdat_ffdf))
testdat_ffdf$fitted6 <- as.ff(as.numeric(NA), length = nrow(testdat_ffdf))

## chunk the ffdf object for use with the predict function
chunks <- chunk.ffdf(testdat_ffdf)    

for(chunkrangeindex in chunks){
  testdat_ffdf$fitted4[chunkrangeindex,] <- 
    predict(mod_sgd1, cbind(1,data.matrix(as.data.frame(testdat_ffdf[chunkrangeindex, ]))[, all.vars(HiggsFormula1[[3]])]), 
                  type = "response")[,1]
}     

for(chunkrangeindex in chunks){
  testdat_ffdf$fitted5[chunkrangeindex,] <- 
    predict(mod_sgd2, cbind(1,data.matrix(as.data.frame(testdat_ffdf[chunkrangeindex, ]))[, all.vars(HiggsFormula2[[3]])]), 
            type = "response")[,1]
}    

for(chunkrangeindex in chunks){
  testdat_ffdf$fitted6[chunkrangeindex,] <- 
    predict(mod_sgd3, cbind(1,data.matrix(as.data.frame(testdat_ffdf[chunkrangeindex, ]))[, all.vars(HiggsFormula3[[3]])]), 
            type = "response")[,1]
}    

##------------------------------------------------------------------------------------------------------

# ROC function from slides
roc <- function(fit, truth, thresholds,
                plot = FALSE, add = FALSE, col = "black") {
  truepos <- sapply(thresholds, function(pthres) {
    sum((fit > pthres) & (truth == 1))/sum(truth == 1)})
  falsepos <- sapply(thresholds, function(pthres) {
    sum((fit > pthres) & (truth == 0))/sum(truth == 0)})
  correctclass <- sapply(thresholds, function(pthres) {
    (sum((fit > pthres) & (truth == 0)) +
       sum((fit > pthres) & (truth == 1)))/length(fit)})
  if (plot) {
    if (add) points(falsepos, truepos, type = "l", col = col)
    else {
      plot(falsepos, truepos, type = "l", col = col)
      abline(0, 1) }
  }
  invisible(data.frame(thresholds, truepos = truepos,
                       falsepos = falsepos,
                       correctclass = correctclass))
}

## i have saved all roc1 to roc6 for easy loading

## set directory

setwd(objectsDir)

## load objects
roc1 <- readRDS("roc1")
roc2 <- readRDS("roc2")
roc3 <- readRDS("roc3")
roc4 <- readRDS("roc4")
roc5 <- readRDS("roc5")
roc6 <- readRDS("roc6")

## code for roc objects
roc1 <- roc(testdat_ffdf$fitted1, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc1$curve <- "bigglm 1"
roc2 <- roc(testdat_ffdf$fitted2, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc2$curve <- "bigglm 2"
roc3 <- roc(testdat_ffdf$fitted3, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc3$curve <- "bigglm 3"
roc4 <- roc(testdat_ffdf$fitted4, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc4$curve <- "sgd 1"
roc5 <- roc(testdat_ffdf$fitted5, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc5$curve <- "sgd 2"
roc6 <- roc(testdat_ffdf$fitted6, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc6$curve <- "sgd 3"

roc1$seq <- seq(0,1,length=100)
roc2$seq <- seq(0,1,length=100)
roc3$seq <- seq(0,1,length=100)
roc4$seq <- seq(0,1,length=100)
roc5$seq <- seq(0,1,length=100)
roc6$seq <- seq(0,1,length=100)

## ROC curves for the bigglm models on the test data
rocs <- rbind(roc1,roc2,roc3)
rocPlot <- ggplot(data = rocs) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "Bigglm") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot)

## ROC curves for the sgd models on the test data
rocs2 <- rbind(roc4,roc5,roc6)
rocPlot2 <- ggplot(data = rocs2) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "sgd") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot2)

## combined plot
rocs3 <- rbind(roc1,roc2,roc3,roc4,roc5,roc6)
rocPlot3 <- ggplot(data = rocs3) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "sgd and bigglm comparison") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot3)


## Comparing ROC curves for each formula
rocs_f1 <- rbind(roc1,roc4)
rocPlot_f1 <- ggplot(data = rocs_f1) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "HiggsFormula1") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot_f1)

## formula 2
rocs_f2 <- rbind(roc2,roc5)
rocPlot_f2 <- ggplot(data = rocs_f2) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "HiggsFormula2") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot_f2)

## formula 3
rocs_f3 <- rbind(roc3,roc6)
rocPlot_f3 <- ggplot(data = rocs_f3) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "HiggsFormula3") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot_f3)

## sample plot to look at rate of sensitivity, specificity, and classification rate
rocPlot_all <- ggplot(data = roc6) +
  geom_line(aes(x = seq, y = truepos, col = "Sensitivity")) +
  geom_line(aes(x = seq, y = 1-falsepos, col = "Specificity")) +
  geom_line(aes(x = seq, y = correctclass, col = "Classification Rate")) +
  labs(x = "Cutoffs", y = "Values", title = "Values of sensitivity, specificity, and classification rate") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot_all)

##------------------------------------------------------------------------------------------------------

## to get the roc curves for other methods, one-dim and rmsprop

## set directory
setwd(objectsDir)

## load results
sgd_onedim1 <- readRDS("sgd_onedim1")
sgd_onedim2 <- readRDS("sgd_onedim2")
sgd_onedim3 <- readRDS("sgd_onedim3")
sgd_rmsprop1 <- readRDS("sgd_rmsprop1")
sgd_rmsprop2 <- readRDS("sgd_rmsprop2")
sgd_rmsprop3 <- readRDS("sgd_rmsprop3")

## one-dim sgd
sgd_onedim1 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula1, lr="one-dim")
sgd_onedim2 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula2, lr="one-dim")
sgd_onedim3 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula3, lr="one-dim")

## rmsprop sgd
sgd_rmsprop1 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula1, lr="rmsprop")
sgd_rmsprop2 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula2, lr="rmsprop")
sgd_rmsprop3 <- sgd_chunk_revised(traindat_ffdf, formula=HiggsFormula3, lr="rmsprop")

## load roc values for one-dim and rmsprop
roc7 <- readRDS("roc7")
roc8 <- readRDS("roc8")
roc9 <- readRDS("roc9")
roc10 <- readRDS("roc10")
roc11 <- readRDS("roc11")
roc12 <- readRDS("roc12")

## create columns for new predicted values
testdat_ffdf$fitted7 <- as.ff(as.numeric(NA), length = nrow(testdat_ffdf))
testdat_ffdf$fitted8 <- as.ff(as.numeric(NA), length = nrow(testdat_ffdf))
testdat_ffdf$fitted9 <- as.ff(as.numeric(NA), length = nrow(testdat_ffdf))
testdat_ffdf$fitted10 <- as.ff(as.numeric(NA), length = nrow(testdat_ffdf))
testdat_ffdf$fitted11 <- as.ff(as.numeric(NA), length = nrow(testdat_ffdf))
testdat_ffdf$fitted12 <- as.ff(as.numeric(NA), length = nrow(testdat_ffdf))

## generate predicted values
for(chunkrangeindex in chunks){
  testdat_ffdf$fitted7[chunkrangeindex,] <- 
    predict(sgd_onedim1, cbind(1,data.matrix(as.data.frame(testdat_ffdf[chunkrangeindex, ]))[, all.vars(HiggsFormula1[[3]])]), 
            type = "response")[,1]
}     

for(chunkrangeindex in chunks){
  testdat_ffdf$fitted8[chunkrangeindex,] <- 
    predict(sgd_onedim2, cbind(1,data.matrix(as.data.frame(testdat_ffdf[chunkrangeindex, ]))[, all.vars(HiggsFormula2[[3]])]), 
            type = "response")[,1]
}    

for(chunkrangeindex in chunks){
  testdat_ffdf$fitted9[chunkrangeindex,] <- 
    predict(sgd_onedim3, cbind(1,data.matrix(as.data.frame(testdat_ffdf[chunkrangeindex, ]))[, all.vars(HiggsFormula3[[3]])]), 
            type = "response")[,1]
}

for(chunkrangeindex in chunks){
  testdat_ffdf$fitted10[chunkrangeindex,] <- 
    predict(sgd_rmsprop1, cbind(1,data.matrix(as.data.frame(testdat_ffdf[chunkrangeindex, ]))[, all.vars(HiggsFormula1[[3]])]), 
            type = "response")[,1]
}     

for(chunkrangeindex in chunks){
  testdat_ffdf$fitted11[chunkrangeindex,] <- 
    predict(sgd_rmsprop2, cbind(1,data.matrix(as.data.frame(testdat_ffdf[chunkrangeindex, ]))[, all.vars(HiggsFormula2[[3]])]), 
            type = "response")[,1]
}    

for(chunkrangeindex in chunks){
  testdat_ffdf$fitted12[chunkrangeindex,] <- 
    predict(sgd_rmsprop3, cbind(1,data.matrix(as.data.frame(testdat_ffdf[chunkrangeindex, ]))[, all.vars(HiggsFormula3[[3]])]), 
            type = "response")[,1]
}

## generate roc functions for the one-dim and rmsprop learning rates
roc7 <- roc(testdat_ffdf$fitted7, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc7$curve <- "one-dim 1"
roc8 <- roc(testdat_ffdf$fitted8, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc8$curve <- "one-dim 2"
roc9 <- roc(testdat_ffdf$fitted9, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc9$curve <- "one-dim 3"
roc10 <- roc(testdat_ffdf$fitted10, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc10$curve <- "rmsprop 1"
roc11 <- roc(testdat_ffdf$fitted11, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc11$curve <- "rmsprop 2"
roc12 <- roc(testdat_ffdf$fitted12, testdat_ffdf$signal,
            seq(0, 1, length = 100), col = "blue")
roc12$curve <- "rmsprop 3"

# adding a sequence for plots comparing classification rate, truepos and falsepos
roc7$seq <- seq(0,1,length=100)
roc8$seq <- seq(0,1,length=100)
roc9$seq <- seq(0,1,length=100)
roc10$seq <- seq(0,1,length=100)
roc11$seq <- seq(0,1,length=100)
roc12$seq <- seq(0,1,length=100)


# roc plots for one-dim for all three formulae
rocs_onedim <- rbind(roc7,roc8,roc9)
rocPlot_onedim <- ggplot(data = rocs_onedim) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "One-dim") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot_onedim)

# roc plots for rmsprop for all three formulae
rocs_rmsprop <- rbind(roc10,roc11,roc12)
rocPlot_rmsprop <- ggplot(data = rocs_rmsprop) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "Rmsprop") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot_rmsprop)

# comparing adagrad and one-dim
rocs_adagrad_onedim <- rbind(roc4,roc5,roc6,roc7,roc8,roc9)
rocPlot_compare_ada_one <- ggplot(data = rocs_adagrad_onedim) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "Adagrad vs one-dim") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot_compare_ada_one)

# comparing one-dim and bigglm
rocs_onedim_bigglm <- rbind(roc1,roc2,roc3,roc7,roc8,roc9)
rocPlot_compare_one_glm <- ggplot(data = rocs_onedim_bigglm) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "bigglm vs one-dim") +
  theme_bw() + theme(legend.position = "top")
print(rocPlot_compare_one_glm)

##------------------------------------------------------------------------------------------------------

## Taken from presentation code chunks
## Estimate log-odds from many independent Bernoullis
set.seed(123)
N <- 10000
y <- rbinom(N, 1, 0.4)
## MLE
betaMLE <- qlogis(mean(y))
## Explicit SVG
betaX <- 0
betaXVec <- numeric(N)
gammaFun <- function(i) 100/i ## try e.g. 0.1/i, 1/sqrt(i), 50/i, 10/i
CFun <- function(beta) 1
for (i in seq.int(N)) {
  prob <- plogis(betaX)
  betaX <- betaX + gammaFun(i) * CFun(betaX) * (y[i] - prob)
  betaXVec[i] <- betaX
}
## Plotting
plot(c(1, N), c(-2, 2), type = "n", xlab = "observation", ylab = "estimate")
points(seq.int(N), betaXVec, type = "l")
abline(h = betaMLE)

# changing the value and running through the figure plotting code
gammaFun <- function(i) 1/i
gammaFun <- function(i) 1000/i
gammaFun <- function(i) 100/sqrt(i)
gammaFun <- function(i) 100/i^2

##------------------------------------------------------------------------------------------------------

## to save space, I include only the roc objects to be loaded

## set directory
setwd(objectsDir)

## load objects
roc_MLE <- readRDS("roc_MLE")
roc_MLE_2 <- readRDS("roc_MLE_2")

roc6$curve <- "sgd | start from zero"
roc3$curve <- "MLE"

roc_check = rbind(roc_MLE,roc_MLE_2,roc3,roc6)
roc_checks <- ggplot(data = roc_check) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "Sensitivity on Starting Values") +
  theme_bw() + theme(legend.position = "top")
print(roc_checks)

##------------------------------------------------------------------------------------------------------

## to save space, I repeat the previous process except with one-dim and store the roc curves
setwd(objectsDir)

## load objects
roc_MLE_onedim <- readRDS("roc_MLE_onedim") # one-dim starting from MLE
roc_MLE_onedim_2 <- readRDS("roc_MLE_onedim_2") # one-dim starting from half of MLE

roc7$curve <- "one-dim | start from zero"
roc3$curve <- "MLE"

roc_check_onedim = rbind(roc_MLE,roc_MLE_2,roc3,roc7)
roc_checks_onedim <- ggplot(data = roc_check_onedim) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "Sensitivity on Starting Values") +
  theme_bw() + theme(legend.position = "top")
print(roc_checks_onedim)

##------------------------------------------------------------------------------------------------------

setwd(objectsDir)

## load objects
attempt1 <- readRDS("attempt1")
attempt2 <- readRDS("attempt2")
attempt3 <- readRDS("attempt3")
attempt4 <- readRDS("attempt4")

## load roc data to save space
roc_attempt1 <- readRDS("roc_attempt1")
roc_attempt2 <- readRDS("roc_attempt2")
roc_attempt3 <- readRDS("roc_attempt3")
roc_attempt4 <- readRDS("roc_attempt4")

# playing with lr.control values for adagrad
attempt1 <- sgd_chunk_revised(formula=HiggsFormula3, x = traindat_ffdf, lr.control=(c(eta=1.1, epsilon=1e-5)))
attempt2 <- sgd_chunk_revised(formula=HiggsFormula3, x = traindat_ffdf, lr.control=(c(eta=0.9, epsilon=1e-5)))
attempt3 <- sgd_chunk_revised(formula=HiggsFormula3, x = traindat_ffdf, lr.control=(c(eta=1.3, epsilon=1e-5)))
attempt4 <- sgd_chunk_revised(formula=HiggsFormula3, x = traindat_ffdf, lr.control=(c(eta=0.7, epsilon=1e-5)))


# plot and compare with value eta=1.1
roc6$curve <- "eta = 1"
roc_eta <- rbind(roc_attempt1,roc_attempt2,roc_attempt3,roc_attempt4,roc6)
roc_eta_compare <- ggplot(data = roc_eta) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "Comparison of eta values for adagrad") +
  theme_bw() + theme(legend.position = "top")
print(roc_eta_compare)

# compare with MLE
roc_eta_bigglm <- rbind(roc_attempt4,roc3)
roc_eta_mle <- ggplot(data = roc_eta_bigglm) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "Comparison of eta=0.7 adagrad with bigglm") +
  theme_bw() + theme(legend.position = "top")
print(roc_eta_mle)

##------------------------------------------------------------------------------------------------------

setwd(objectsDir)

## load objects
attempt5 <- readRDS("attempt5")
attempt6 <- readRDS("attempt6")
attempt7 <- readRDS("attempt7")
attempt8 <- readRDS("attempt8")

## to save space, I load only the roc objects I have created
roc_attempt5 <- readRDS("roc_attempt5")
roc_attempt6 <- readRDS("roc_attempt6")
roc_attempt7 <- readRDS("roc_attempt7")
roc_attempt8 <- readRDS("roc_attempt8")

## comparison with HiggsFormula 1 and 2 for eta=0.7
attempt5 <- sgd_chunk_revised(formula=HiggsFormula1, x = traindat_ffdf, lr.control=(c(eta=0.7, epsilon=1e-5)))
attempt6 <- sgd_chunk_revised(formula=HiggsFormula2, x = traindat_ffdf, lr.control=(c(eta=0.7, epsilon=1e-5)))
attempt7 <- sgd_chunk_revised(formula=HiggsFormula3, x = traindat_ffdf, lr.control=(c(eta=0.8, epsilon=1e-5)))
attempt8 <- sgd_chunk_revised(formula=HiggsFormula3, x = traindat_ffdf, lr.control=(c(eta=0.71, epsilon=1e-5)))


# sgd_chunk_revised with "adagrad" and eta=0.7 and bigglm
roc3$curve <- "bigglm1"
roc_attempt4$curve <- "eta=0.7 | formula 3"
roc_eta_07 <- rbind(roc_attempt4,roc_attempt5,roc_attempt6, roc1,roc2,roc3)
rocs_eta_077 <- ggplot(data = roc_eta_07) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "eta=0.7 performance") +
  theme_bw() + theme(legend.position = "top")
print(rocs_eta_077)

# sgd_chunk_revised with "adagrad" and eta=0.7 and bigglm
roc_attempt8$curve <- "eta=0.71 | formula 3"
roc_attempt7$curve <- "eta=0.8 | formula3"
roc_078 <- rbind(roc_attempt4, roc_attempt8, roc_attempt7)
rocs_078 <- ggplot(data = roc_078) +
  geom_line(aes(x = falsepos, y = truepos, col = curve)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(x = "False positive", y = "True positive", title = "Sensitivity to eta values") +
  theme_bw() + theme(legend.position = "top")
print(rocs_078)


##------------------------------------------------------------------------------------------------------

setwd(objectsDir)

## load objects
one_pass <- readRDS("one_pass") # sgd_chunk_function with npasses=1
two_pass <- readRDS("two_pass")
four_pass <- readRDS("four_pass")
five_pass <- readRDS("five_pass")
ten_pass <- readRDS("ten_pass")
passes <- readRDS("passes") # data frame containing coefficients from different number of passes

## code sample
one_pass <- sgd_chunk_revised(traindat_ffdf, chunk_size=NULL, formula=HiggsFormula1, no_passes=1)
two_pass <- sgd_chunk_revised(traindat_ffdf, chunk_size=NULL, formula=HiggsFormula1, no_passes=2)
four_pass <- sgd_chunk_revised(traindat_ffdf, chunk_size=NULL, formula=HiggsFormula1, no_passes=4)
five_pass <- sgd_chunk_revised(traindat_ffdf, chunk_size=NULL, formula=HiggsFormula1, no_passes=5)
ten_pass <- sgd_chunk_revised(traindat_ffdf, chunk_size=NULL, formula=HiggsFormula1, no_passes=10)

passes <- data.frame(one_pass$coefficients,two_pass$coefficients,mod_sgd1$coefficients,four_pass$coefficients,
                     five_pass$coefficients,ten_pass$coefficients)
names(passes) <- c("one pass","two pass","three pass","four pass","five pass","ten pass")

# euclidean distance of each solution from each other
distance_mat <- data.frame(matrix(0,ncol(passes),ncol(passes)))
for (i in 1:ncol(passes)){
  for (j in 1:ncol(passes)){
    distance_mat[i,j] <- sqrt(t(as.matrix(passes[i]))%*%as.matrix(passes[j]))
  }
}
names(distance_mat) <- names(passes)
rownames(distance_mat) <- names(passes)

passes <- cbind(passes,seq(1,nrow(passes)))
names(passes)[ncol(passes)] <- "coefficient number"
passes <- melt(passes, id.vars="coefficient number")
coeffplot <- ggplot(passes, aes(`coefficient number`,value, col=variable)) + 
  geom_point(size=1) + labs(title="Variation in coefficient values for different number of passes over data")
coeffplot

##------------------------------------------------------------------------------------------------------

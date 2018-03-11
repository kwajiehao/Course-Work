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
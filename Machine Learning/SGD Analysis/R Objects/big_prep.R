# Update the installed packages first
update.packages(ask = FALSE, repos = 'http://cran.rstudio.com/') # Install the packages for this tutorial
PFpackages <- c('biglm', 'ffbase', 'ggplot2', 'sgd', 'bigmemory',
                'glmnet', 'enrichwith', 'car')
install.packages(PFpackages, repos = 'http://cran.rstudio.com/')

# ffbase Provides support for data.frame like objects that connect to
# files stored out of memory
library("ffbase")
# bigmemory provides the big.matrix object which links to matrices stored # out of memory
library("bigmemory")
# biglm implements a bounded-memory fitter for glms library("biglm")
# ggplot2 is for flexible plotting library("ggplot2")
# sgd is for stochastic gradient descent methods library("sgd")
# car for data sets and tools
library("car")


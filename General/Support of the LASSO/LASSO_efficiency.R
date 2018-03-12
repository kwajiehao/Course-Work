library(glmnet)
library(MASS)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)

generator <- function(a, rho, n){
  count <- 0
  for (i in 1:100){
    cov <- c(1,rho,rho,1)
    cov <- matrix(cov,2,2)
    x <- mvrnorm(n, mu=rep(0,2), Sigma=cov)
    
    beta <- c(a,0)
    
    y <- beta %*% t(x) + rnorm(n)
    tryCatch({
      fit = coef(cv.glmnet(x,y,intercept=FALSE))
      if (a == 0){
        if (fit[2] == 0 & fit[3] == 0){
          count <- count + 1
        }
      } else{
        if (fit[2] != 0 & fit[3] == 0){
          count <- count + 1
        }
      }
    },error=function(e){})
  }
  return(count)
}

generator(0.01, 0.5,100)

# running cv.glmnet with the same input does not produce the same output

#################################################

# for rho, it would be interesting to see how values change for increments 
# in rho

# choice of n was also arbitrary over the range

# my intuition says that larger values of a will make it easier to 
# distinguish between a and zero ( the second parameter). I am curious to 
# see if there is a threshold beyond which glmnet recovers the correct
# support

rho_list <- seq(0.01, 0.99, length=10)
n <- c(20, 30, 40, 50, 80, 100, 200, 300, 500, 1000)
a <- c(-50, -20, -1, -0.5, -0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1, 0.5, 1, 20, 50)


#################################################

probability_table <- function(ns, rhos, as){
  table <- array(rep(0, length(ns)*length(rhos)*length(as)), 
                 c(length(ns), length(rhos), length(as)))
  
  for (i in 1:length(ns)){
    for (j in 1:length(rhos)){
      for (k in 1:length(as)){
        table[i,j,k] <- generator(as[k],rhos[j],ns[i])
        print(c(i,j,k))
      }
    }
  }
  return(table)
}
  
results4 <- probability_table(n,rho_list,a)

# settings for results1
rho_list <- seq(0.01, 0.99, length=11)
n <- c(10, 20, 30, 40, 50, 80, 100, 200, 300, 500, 1000)
a <- c(-100,-20, -1, -0.05, -0.01, 0, 0.01, 0.05, 1, 20, 100)

# settings for results2
rho_list <- seq(0.01, 0.99, length=11)
n <- c(10, 20, 30, 40, 50, 80, 100, 200, 300, 500, 1000)
a <- c(-50, -20, -1, -0.5, -0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1, 0.5, 1, 20, 50)

# settings for results3 (with lambda.min)
rho_list <- seq(0.01, 40.99, length=10)
n <- c(20, 30, 40, 50, 80, 100, 200, 300, 500, 1000)
a <- c(-50, -20, -1, -0.5, -0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1, 0.5, 1, 20, 50)

# settings for results4 (with double condition in generator)
rho_list <- seq(0.01, 0.99, length=10)
n <- c(20, 30, 40, 50, 80, 100, 200, 300, 500, 1000)
a <- c(-50, -20, -1, -0.5, -0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1, 0.5, 1, 20, 50)

getwd()
path <- "/Users/kwajiehao/Documents/Github/sup/1/Deterministic Models/convex optimization/hw/"
setwd(path)
save(results, file="results1.Rda")
save(results, file="results2.Rda")
save(results3, file="results3.Rda")
save(results4, file="results4.Rda")
load("results1.Rda")
load("results2.Rda")
load("results3.Rda")
load("results4.Rda")

#################################################

# sample heatmap plot

a1 <- as.data.frame(results4[,11,])
names(a1) <- a
rownames(a1) <- n
a1 <- cbind(a1, n)
a1$n <- as.factor(a1$n)
a1 <- melt(a1, id.vars = "n")

plot1 <- ggplot(a1, aes(variable, n)) + geom_tile(aes(fill = value),colour = "white") + 
  scale_fill_gradient(low = "white", high = "orangered1") + labs(title=bquote(rho))

plot1

# create a function that generates heatmaps for different values of rho

generate_map <- function(data, index){
  
  # convert output into dataframe
  temp <- as.data.frame(data[,index,])
  names(temp) <- a
  rownames(temp) <- n
  temp <- cbind(temp, n)
  # make n a factor so that we can reshape temp
  temp$n <- as.factor(temp$n)
  temp <- melt(temp, id.vars="n")
  
  titl <- bquote(rho == .(rho_list[index]))
  
  plot_ <- ggplot(temp, aes(variable, n)) + geom_tile(aes(fill = value),colour = "white") + 
    scale_fill_gradient(low = "white", high = "blue") + 
    labs(title = titl, x="a")
  
  return(plot_)
}

generate_map(results4,1)

# loop through the rho_list

for(i in 1:length(rho_list)){
  assign(paste0("df",i), generate_map(results4,i))
}

grid.arrange()

# The smaller the lambda, the larger the constraint area and the greater 
# possibility for non-zero estimates of parameters

#################################################

# testing

cov <- c(1,0.01,0.01,1)
cov <- matrix(cov,2,2)
x <- mvrnorm(200, mu=rep(0,2), Sigma=cov)

beta <- c(0.1,0)

y <- beta %*% t(x) + rnorm(200)
vv <- cv.glmnet(x,y,intercept=FALSE)
coef(vv, s="lambda.min")
coef(vv)

# rmarkdown::render('/Users/kwajiehao/Documents/Github/sup/1/Deterministic Models/convex optimization/hw/HW.Rmd')

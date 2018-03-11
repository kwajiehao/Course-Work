library(FNN)
require(parallel)
library(gridExtra)
library(grid)
library(ggplot2)
library(reshape2)

point_generator <- function(n,d){
  
  if (d==1){
    # generates n points from the d-dimensional uniform distribution 
    sample <- replicate(n,runif(d))
    
    # obtaining the averages of the components of X so we get the conditional y distribution
    cond_y <- sample
  }
  else{
    sample <- t(replicate(n,runif(d)))
    cond_y <- apply(sample,1,mean)
  }
  
  # sample the conditional y given x
  sample <- cbind(sample, runif(n))
  sample[,d+1] <- ifelse(sample[,d+1] > cond_y, 1, 0)
  
  return(sample)
}


# Step 1: Generate training data for given n and d
# Step 2: Redraw samples 100 times and average the risk
# Step 3: Store results in array

risk_calculator_V1 <- function(nn, n, d){ # note! this is extremely inefficient
  
  sample <- point_generator(n,d)
  labels <- rep(0,n)
  
  for (i in 1:n){
    
    # make a copy of data that we update
    copy2 <- data.frame(cbind(rep(0,n),sample[,d+1]))
    names(copy2) <- c("dist","label")
    
    for (j in 1:d){
      
      #summing distance over d-dimensions
      copy2[,1] <- copy2[,1] + (sample[,j] - sample[i,j])^2
    }
    
    copy2 <- data.frame(copy2)
    copy2 <- copy2[with(copy2,order(dist)),]
    
    # we go from index 2 because the nearest neighbor of a point will be itself
    labels[i] <- ifelse(sum(copy2[,2][2:nn+1])> (nn/2) , 1, 0)
  }
  
  indicator <- ifelse(labels != sample[,d+1], 1, 0)
  return(mean(indicator))
}

# function to find results

risk <- function(nvalue, dvalue, nnvalue, iter){
  
  # generate an array to store data
  table <- array(rep(0, length(nvalue)*length(dvalue)*length(nnvalue)), 
                 c(length(nvalue), length(dvalue), length(nnvalue)))
  
  for (a in 1:length(nvalue)){
    for (b in 1:length(dvalue)){
      for (c in 1:length(nnvalue)){
        
        # generate a training set
        train <- point_generator(nvalue[a],dvalue[b])
        
        # separate the class for use with FNN package's knn function
        train_class <- train[,(dvalue[b]+1)]
        train<-train[,1:dvalue[b]]
        
        # store the risk calculated in each iteration
        risks <- rep(0,iter)
        
        for (i in 1:iter){
          
          # generate test data
          test <- point_generator(nvalue[a],dvalue[b])
          test_class <- test[,(dvalue[b]+1)]
          test<-test[,1:dvalue[b]]
          
          # apply knn function
          nearest <- knn(as.matrix(train), as.matrix(test), cl = factor(train_class), k = nnvalue[c])
          nearest_label <- as.numeric(levels(nearest))[nearest]
          
          # compare with test labels, assigning a loss of one if misclassified
          indicator <- ifelse(test_class != nearest_label, 1, 0)
          
          # store results
          risks[iter] <- mean(indicator)
        }
        
        # update table
        table[a,b,c] <- mean(risks)
      }
    }
  }
  
  return(table)
}

# choose values of n and d

n <- c(50,100,1000,5000,10000)
d <- c(1,2,3,10,50)
nn <- c(1,3,5,7,9)
iter <- 1000

system.time(results <- risk(n, d, nn, iter))

# test
a <- apply(test, c(3), appfunc)


# generate training data to be used
train <- point_generator(50,2)
train_class <- train[,(2+1)]
train<-train[,1:2]


# create a function that uses knn, then use apply
appfunc <- function(x,nn){
  return(knn(as.matrix(train),as.matrix(x),cl=factor(train_class),k=nn))
}

#knn(as.matrix(train),as.matrix(test[,1]),cl=factor(train_class),k=abc)

# test data
test <- replicate(100,point_generator(50,2))
test_class <- test[,3,]
test <- test[,1:2,]

abc <- 1

# parallelize
nclusters<-detectCores()-1
cl<-makeCluster(nclusters)
clusterExport(cl=cl,c("appfunc", "train", "train_class","abc"))
clusterEvalQ(cl, library(FNN)) 
system.time(results <- parApply(cl=cl,X = test,MARGIN = c(3),FUN = appfunc, nn=abc))
stopCluster(cl)

# convert results
results2 <- as.numeric(results)
results2 <- matrix(results2, dim(results)[1], dim(results)[2])
mean(results2 != test_class)

####################################################################################


## Version 2 with Apply function ##


# create a function that uses knn, then use apply

appfunc <- function(x, nn){
  return(knn(as.matrix(train),as.matrix(x),cl=factor(train_class),k=nn))
}

risk2 <- function(nvalue, dvalue, nnvalue, iter){
  
  # generate an array to store data
  table <- array(rep(0, length(nvalue)*length(dvalue)*length(nnvalue)), 
                 c(length(nvalue), length(dvalue), length(nnvalue)))
  
  for (a in 1:length(nvalue)){
    for (b in 1:length(dvalue)){
        for (c in 1:length(nnvalue)){
          if (dvalue[b]==1){
              
              # generate a training set
              train <- point_generator(nvalue[a],dvalue[b])
              
              # separate the class for use with FNN package's knn function
              train_class <- train[,(dvalue[b]+1)]
              train<-train[,1:dvalue[b]]
              
              # generate test data
              test <- replicate(iter,point_generator(nvalue[a],dvalue[b]))
              test_class <- test[,(dvalue[b]+1),]
              test<-test[,1:dvalue[b],]
              
              # store nn-value
              abc <- nnvalue[c]
              
              # apply knn function and parallelize
              nclusters<-detectCores()-1
              cl<-makeCluster(nclusters)
              
              # export relevant data
              clusterExport(cl=cl,c("appfunc", "train", "train_class","abc"),envir=environment())
              
              # export relevant libraries
              clusterEvalQ(cl, library(FNN)) 
              
              # apply parApply
              print(list(a,b,c)) 
              results <- parApply(cl=cl,X = test, MARGIN = 2,FUN = appfunc, nn=abc) 
              stopCluster(cl)
              
              # compare with test labels, assigning a loss of one if misclassified
              results2 <- as.numeric(results)
              results2 <- matrix(results2, dim(results)[1], dim(results)[2])
              
              # update table
              table[a,b,c] <- mean(results2 != test_class)
        }
      else{
        # generate a training set
        train <- point_generator(nvalue[a],dvalue[b])
        
        # separate the class for use with FNN package's knn function
        train_class <- train[,(dvalue[b]+1)]
        train<-train[,1:dvalue[b]]
        
        # generate test data
        test <- replicate(iter,point_generator(nvalue[a],dvalue[b]))
        test_class <- test[,(dvalue[b]+1),]
        test<-test[,1:dvalue[b],]
        
        # store nn-value
        abc <- nnvalue[c]
        
        # apply knn function and parallelize
        nclusters<-detectCores()-1
        cl<-makeCluster(nclusters)
        
        # export relevant data
        clusterExport(cl=cl,c("appfunc", "train", "train_class","abc"),envir=environment())
        
        # export relevant libraries
        clusterEvalQ(cl, library(FNN)) 
        
        # apply parApply
        print(list(a,b,c)) 
        results <- parApply(cl=cl,X = test, MARGIN = 3,FUN = appfunc, nn=abc) 
        stopCluster(cl)
        
        # compare with test labels, assigning a loss of one if misclassified
        results2 <- as.numeric(results)
        results2 <- matrix(results2, dim(results)[1], dim(results)[2])
        
        # update table
        table[a,b,c] <- mean(results2 != test_class)
      }
      }

    }
  }
  
  return(table)
}

n <- c(50,100,1000,5000,10000)
d <- c(1,2,3,10)
nn <- c(1,3,5,7,9)
iter <- 1000

# somehow, the global environment's train variable interferes with the function call
system.time(results <- risk2(n, d, nn, iter))

dim(results)

####################################################################################


## Version 3 after Misunderstanding ##

risk3 <- function(nvalue, dvalue, nnvalue, iter){
  
  # generate an array to store data
  table <- array(rep(0, length(nvalue)*length(dvalue)*length(nnvalue)), 
                 c(length(nvalue), length(dvalue), length(nnvalue)))
  
  for (a in 1:length(nvalue)){
    for (b in 1:length(dvalue)){
      for (c in 1:length(nnvalue)){
          # generate a training set
          train <- point_generator(nvalue[a],dvalue[b])
          
          # separate the class for use with FNN package's knn function
          train_class <- train[,(dvalue[b]+1)]
          train<-train[,1:dvalue[b]]
          
          # generate test data
          test <- point_generator(iter,dvalue[b])
          test_class <- test[,(dvalue[b]+1)]
          test<-test[,1:dvalue[b]]
          
          # apply knn
          nearest <- knn(as.matrix(train), as.matrix(test), cl = factor(train_class), k = nnvalue[c])
          nearest_label <- as.numeric(levels(nearest))[nearest]
          
          # compare with test labels, assigning a loss of one if misclassified
          indicator <- ifelse(test_class != nearest_label, 1, 0)
          
          # update table
          table[a,b,c] <- mean(indicator)
      }
      
    }
  }
  
  return(table)
}

n <- c(50,100,1000,5000,10000)
d <- c(1,2,3,10,50)
nn <- c(1,3,5,7,9)
iter <- 10000 

system.time(results2 <- risk3(n, d, nn, iter))

dim(results2)


####################################################################################

# heatmap generator for 3d array
generate_map <- function(data, index, scale_lower, scale_upper){
  
  # convert output into dataframe
    temp <- as.data.frame(data[,,index])
    names(temp) <- d
    rownames(temp) <- n
    temp <- cbind(temp, n)
    
    temp$n <- as.factor(temp$n)
    temp <- melt(temp, id.vars="n")
    
    plot_ <- ggplot(temp, aes(variable, n)) + geom_tile(aes(fill = value),colour = "white") + 
      scale_fill_gradient(low = "white", high = "brown", limits=c(scale_lower,scale_upper)) + 
      labs(title = paste0(nn[index],"-nearest neighbor"), y="train size", x="dimensions")
  
  return(plot_)
}

generate_map(results2,1,0.250, 0.5)
generate_map(results2,2,0.250, 0.5)
generate_map(results2,3,0.250, 0.5)
generate_map(results2,4,0.250, 0.5)
generate_map(results2,5,0.250, 0.5)

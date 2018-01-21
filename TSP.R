## Problem 5 TSP Computational Assignment)

Visit the website: http://www.math.uwaterloo.ca/tsp/world/countries.html
Solve the Traveling Salesman Problem for Uruguay based on the dataset provided. You can use your favorite programming language and solution method for the TSP. Provide a printout of your code with detailed documentation, and compare the optimal solution you obtain to the one available at the website.

### Answer

As we've discussed in class, it is problematic to apply dynamic programming to solve this problem as it would require calculating and comparing all possible paths, which happens to be $\frac{(n-1)!}{2}$ total paths. To give an idea of how big that is, inputting that value into Google returns "infinity" as an answer. This is exacerbated by the fact that every node is accessible from every node, giving a complete graph. Even label-correcting methods don't work well since it takes a long time to find an upper bound and also the sheer number of possibilities. We begin with a nearest-neigbor approach: starting from one point, the decision rule is to move to the nearest point to the current point. For this, we need to compute a distance matrix between cities using the location data.

```{r, echo = TRUE}
# clean environment
rm(list=ls())

# packages
library(ggplot2)

# load file
directory <- "/Users/kwajiehao/Documents/GitHub/sup/2/Stochastic Optimization/2/"
file <-"uy734.tsp.txt"
cities <- read.table(paste0(directory,file), header=TRUE, sep=" ")

# calculate pairwise distances - distance analogue of adjacency matrix. however,
# if there were more cities than 734, would probably need to use adjacency lists
# since matrix would be too big.

no_cities <- nrow(cities)
dists <- matrix(0, nrow=no_cities, ncol=no_cities)

# populate the matrix for pairwise distances

for (i in 1:no_cities){
  for (j in 1:no_cities){
    if (i==j){
      dists[i,j] <- Inf
    }
    else{
      dists[i,j] <- sqrt((cities[i,2]-cities[j,2])^2 + (cities[i,3]-cities[j,3])^2) 
    }
  }
}
```

Next, we implement the greedy algorithm of moving to the nearest neighbor.

```{r, echo = TRUE}

greedy_track <- function(){
  
  # records the best paths for each starting point
  best <- rep(0,nrow(cities))
  path <- data.frame(matrix(0,nrow(cities),nrow(cities)))
  colnames(path) <- seq(nrow(cities))
  
  for (i in 1:nrow(cities)){
    
    # S gives a list of numbers from which we remove to find the 
    # available decision space
    S <- seq(nrow(cities))
    S <- S[-i]
    
    # sol tracks the current path length
    sol <- 0
    
    # pos tracks the current position 
    pos <- i
    
    # record nodes traveled
    traveled <- rep(0,nrow(cities))
    traveled[1] <- i
    
    for (j in 2:nrow(cities)){
      
      G <- S
      
      # we are trying to find the nearest neighbor and track its 
      # position so we cycle through the eligible neighbors 
      current_min <- Inf
      min_position <- 0
      
      for (k in 1:length(G)){
        if (dists[pos,G[k]] < current_min){
          current_min <- dists[pos,G[k]]
          min_position <- G[k]
        }
      }
      
      # update current position to the nearest neighbor's position
      pos <- min_position
      traveled[j] <- pos
      
      # remove current position from decision pool
      S <- S[S!=pos]
      
      # add to path length
      sol <- sol + current_min
    }
    
    # update the best path and distance for given start point
    best[i] <- sol + dists[pos,i]
    path[,i] <- traveled
    
  }
  
  return(list(a=best, b=path))
}

# storing the solution
sol2 <- greedy_track()
which.min(sol2$a) # this gives the starting point which gives 
# the shortest path, which is node 266 with a
# distance of 96514.81

# making a copy of the solution
solution <- sol2$b
best <- sol2$a 
```

Next, we improve upon what we already have by running the two-opt algorithm repeatedly on the optimal paths we've found for each starting point. From position $1$ to $n$, the algorithm calculates if it's possible to reduce the distance of the path by swapping edges at two points $j$ and $k$.

```{r, echo = TRUE, eval=FALSE}

# referencing the number of cities
len <- nrow(cities)

twoopt <- function(pathh,minim){ #takes a path and a distance as input
  
  path <- pathh
  mini <- minim
  
  for (j in 2:(len-1)){
    
    # comparing differences in edge costs for node j and k
    # in other words, if we exchange the position of node j with
    # node k, how does the total path distance change? If the change is
    # negative, we switch the positions of the nodes as we have found a 
    # smaller distance
    
    for (k in (j+1):(len)){
      if (k==(len)){ # there is a special case when the second point is
        # the last node as it has to be connected to the 
        # starting node
        edge1 <- dists[path[j],path[j-1]]
        edge2 <- dists[path[k],path[1]]
        edge3 <- dists[path[j-1],path[k]]
        edge4 <- dists[path[j],path[1]]
        
        if (edge1 + edge2 > edge3 +edge4){ # modifying the path 
          path <- c(path[1:(j-1)], rev(path[j:k]))
          mini <- mini - edge1 - edge2 + edge3 + edge4 # updating distance
        }
      }
      else{ # the general case when the second node is not the ending node
        edge1 <- dists[path[j],path[j-1]]
        edge2 <- dists[path[k],path[k+1]]
        edge3 <- dists[path[j-1],path[k]]
        edge4 <- dists[path[j],path[k+1]]
        
        if ((edge1 + edge2) > (edge3 + edge4)){
          path <- c(path[1:(j-1)], rev(path[j:k]), path[(k+1):len])
          mini <- mini - edge1 - edge2 + edge3 + edge4
        }
      }
    }
  }
  return(list(a=path, b=mini))
} 
```

Next, we run the code and find the best paths.

```{r, echo = TRUE,eval=FALSE}
random_search <-  seq(734) 
# this tells the proceeding loop which paths (which starting points)
# to conduct two-opt on. To conduct it on every starting point, we
# set this as the vector with elements from 1 to 734.

system.time(for (i in 1:length(random_search)){
  
  # store the details of the path being investigated 
  path <- as.vector(solution[,random_search[i]])
  mini <- best[random_search[i]]
  
  # conduct two-opt five times on each given path
  a <- twoopt(path,mini)
  b <- twoopt(a$a, a$b)
  c <- twoopt(b$a, b$b)
  d <- twoopt(c$a, c$b)
  e <- twoopt(d$a, d$b)
  
  # updating results
  solution[,random_search[i]] <- e$a
  best[random_search[i]] <- e$b
}) 
```


For two rounds of two-opt, the code took 297 seconds to run  returning a minimum distance of 85424.85 (starting point: 269).
For three rounds, the code took 405 seconds to run  returning a minimum distance of 84431.12 (starting point: 269).
For four rounds, the code took 554 seconds to run  returning a minimum distance of 84063.39 (starting point: 269).
For four rounds, the code took 675 seconds to run  returning a minimum distance of 84024.81 (starting point: 117).

For reference, the best solution provided has a distance of 79114. 

Instead of trying all 754 paths we obtained from the nearest neighbor search, we can try a random subset of starting points multiple times. We obtain this random subset by changing the "random_search".

```{r, echo = TRUE, eval=FALSE}

# for example, the code below includes a random sample of 50 paths
# from the nearest neighbor output and includes an additional path
# starting from point 266
random_search <-  c(266,sample(seq(734)[-266],50))

```

With a random sample of 50 points, it only took 46 seconds. The optimal output in this case was the same as was found above, point 117 with a distance of 84024.81, which is not too far of the best found answer of 79114. In this case, we were lucky to get the globally optimum path, but we can run a random search again with a different sample of 50 points (and even repeat this more than once) to obtain a very good solution without going through all the paths. Below, we've attached the optimal path as an image.

```{r, echo = FALSE, fig.align="center", out.width = "600px"}
knitr::include_graphics("shortestpath.png")
```
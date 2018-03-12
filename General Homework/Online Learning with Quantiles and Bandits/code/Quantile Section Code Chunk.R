library(ggplot2)

################################# Initial Exploration ##################################

# generate points according to distribution
point_generator <- function(dist,n){
  if (dist=="poisson"){
    return(rpois(n,25000))
  }
  
  if (dist=="normal"){
    return(rnorm(n,25000,1000))
  }
  
}

# split into train and test
df <- data.frame(X = point_generator("normal",2000))
train <- data.frame(X= df[1:500,])
test <- data.frame(X= df[501:nrow(df),])

# plot hist

ggplot(train,aes(X)) + geom_histogram(binwidth= 300, col="red", fill="blue", alpha = .2) + 
  theme(axis.title.x=element_blank())


################################# Run Algorithm ##################################

library(gridExtra)

mean <- mean(train$X)
variance <- var(train$X)

# calculate delta
delta2 <- dnorm(qnorm(0.9,mean,sqrt(variance)),mean,sqrt(variance))

# run algorithm
algorithm <- function(phi, tau, epsilon, delta, horizon, data){
  
  info <- data[,1]
  
  # store the values of n_star
  current_n_star <- rep(0, horizon)
  
  # store the values of epsilon
  epsilon_list <- rep(0, horizon)
  
  # store the values of the quantiles
  quantiles <- rep(0, horizon)
  
  for (i in 1:horizon){
    
    if (i==1){
      # gather information in first step
      current_n_star[i] <- 0
      epsilon_list[i] <- 0
      # print(i)
    }
    
    else{
      
      if (i==2){
        # gather information on quantiles but no stopping condition yet
        quantiles[i] <- quantile(info[1:i], tau, type = 8)
        epsilon_list[i] <- epsilon * quantiles[i]
        current_n_star[i] <- log(1/phi)/(2*(delta*epsilon_list[i])^2)
        # print(i)
      }
      
      else{
        if (i==horizon){
          return(list(offer = info[i], epsilons = epsilon_list, n_s = current_n_star, quant = quantiles, stop=(i-1)))
        }
        
        else{
          
          if (i > current_n_star[i-1]){
            if (info[i] > quantiles[i-1]){
              return(list(offer = info[i], epsilons = epsilon_list, n_s = current_n_star, quant = quantiles, stop=(i-1)))
            }
          }
          quantiles[i] <- quantile(info[1:i], tau, type = 8)
          epsilon_list[i] <- epsilon * quantiles[i]
          current_n_star[i] <- log(1/phi)/(2*(delta*epsilon_list[i])^2)
          # print(i)
        }
      }
    }
    
  }
}

a <- algorithm(phi=0.01, tau=0.9, epsilon=0.02, delta=delta2, horizon=1500, data=test)

n <- data.frame(X = a$n_s[2:length(a$n_s)])
epsilon <- data.frame(X=a$epsilons[2:length(a$n_s)])
quants <- data.frame(X=a$quant[2:length(a$n_s)])

# plot sequence of n_s
plot1 <- ggplot(n,aes(y=X, x=seq(1,nrow(n)))) + geom_line(color="blue") + geom_point(size=0.5, position="jitter") + labs(x="Time Period", y="n*") +
  scale_x_continuous(limits = c(0,a$stop)) + scale_y_continuous(limits=c(200,350))

# plot sequence of quantiles
plot2 <- ggplot(quants,aes(y=X, x=seq(1,nrow(n)))) + geom_line(color="blue") + geom_point(size=0.5, position="jitter") + labs(x="Time Period", y=expression(hat(Q))) + scale_x_continuous(limits = c(0,a$stop)) + scale_y_continuous(limits=c(24000,27000))

grid.arrange(plot1,plot2, ncol=2)

################################# Comparison Test ##################################


# odds algorithm
odds <- function(data){
  info <- data[,1]
  
  # find the ratio n/e
  n <- floor(length(info)/exp(1))
  
  # store the max value in first n/e values
  to_beat <- max(info[1:n])
  
  for (i in n+1:length(info)){
    if (i==length(info)){
      return(info[length(info)])
    }
    else{
      if (info[i] > to_beat){
        return(info[i])
      }
    }
    
  }
}

# doing it 1000 times

results <- replicate(2,rep(0,1000))

# once we have learnt the distribution we don't have to re-estimate
for (i in 1:1000){
  test <- data.frame(X = point_generator("normal",1500))
  b <- algorithm(phi=0.01, tau=0.9, epsilon=0.02, delta=delta2, horizon=1500, data=test)
  results[i,1] <- b$offer
  results[i,2] <- odds(test)
}

df2 <- data.frame(results)
names(df2) <- c("quant","odds")

plot3 <- ggplot(df2,aes(quant)) + geom_histogram(binwidth= 100, col="red", fill="blue", alpha = .2) + theme(axis.title.x=element_blank()) + geom_vline(xintercept = qnorm(0.9,25000,1000)) + labs(title="Quantile Solution")

plot4 <- ggplot(df2,aes(odds)) + geom_histogram(binwidth= 100, col="red", fill="blue", alpha = .2) + theme(axis.title.x=element_blank()) + geom_vline(xintercept = qnorm(0.9,25000,1000)) + labs(title="General Solution")

grid.arrange(plot3,plot4, ncol=1)



################################# Poisson Example ##################################

df2 <- data.frame(X = point_generator("poisson",2000))
train2 <- data.frame(X= df2[1:500,])
test2 <- data.frame(X= df2[501:nrow(df2),])

# use 500
ggplot(train2,aes(X)) + geom_histogram(binwidth= 100, col="red", fill="blue", alpha = .2) + 
  theme(axis.title.x=element_blank())

# we see that we have a symmetric distribution similar to the normal distribution.

# next we estimate the moments

mean2 <- mean(train2$X)
variance2 <- var(train2$X)

# calculate delta
delta3 <- dnorm(qnorm(0.9,mean2,sqrt(variance2)),mean2,sqrt(variance2))

b <- algorithm(phi=0.01, tau=0.9, epsilon=0.02, delta=delta3, horizon=1500, data=test2)

n2 <- data.frame(X = b$n_s[2:length(b$n_s)])
epsilon2 <- data.frame(X=b$epsilons[2:length(b$n_s)])
quants2 <- data.frame(X=b$quant[2:length(b$n_s)])
b$offer # true value of quantile is qpois(0.9,lambda = 25000)

# plot sequence of n_s
ggplot(n2,aes(y=X, x=seq(1,nrow(n)))) + geom_line(color="blue") + geom_point(size=0.5, position="jitter") + labs(x="Time Period") +
  scale_x_continuous(limits = c(0,b$stop)) + labs(y="n*") + scale_y_continuous(limits=c(7,9))

# plot sequence of quantiles
ggplot(quants2,aes(y=X, x=seq(1,nrow(n)))) + geom_line(color="blue") + geom_point(size=0.5, position="jitter") + labs(x="Time Period") +
  scale_x_continuous(limits = c(0,b$stop)) + labs(y=expression(hat(Q))) + 
  scale_y_continuous(limits=c(24000,26000))

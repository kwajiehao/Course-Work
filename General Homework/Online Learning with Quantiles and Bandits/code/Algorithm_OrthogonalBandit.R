## Setup
rm(list=ls())
setwd("/home/jan/Documents/Courses/T2_SMO/Project/MultiArmedBandit")

## Load packages
require(tidyr)
require(stats) # for complete.case to efficiently remove NA rows

#############################################################################################
## Get data
#############################################################################################

## Load data (use christian brownlees data) and transform into matrix with returns of series in columns

files <- list.files(path=paste0("/home/jan/Documents/Courses/T2_SMO/Project/",
                       "MultiArmedBandit/Data"))  # load list of files in directory
df_files <- as.data.frame(files) # transform list into dataframe and clean ".csv" part
df_files <- separate(data = df_files, col = files, into = c("series"), sep = "\\.csv")


# Load individual csv series into dataframe containing only date and log return
# derived from adjusted closing price

for (series in df_files$series){
  df <- read.table(paste0("./Data/",series,".csv"), header=TRUE, sep=',')
  D <- df[,c("timestamp","adjusted_close")] # select only relevant columns
  T <- nrow(D)  # number of observations in data set
  D <- D[T:1, ] # reverse order to be temporally increasing
  D[,1] <- as.Date(as.character(D[,1]), '%Y-%m-%d' ) # transform into date format
  D[,3] <- c(NA,D[2:dim(D)[1],2]/D[1:(dim(D)[1]-1),2]) # calculate gross returns
  D_final <- subset(D, select = c(timestamp, V3)) # only take time and log return
  colnames(D_final) <- c("date",series) # change column names
  assign(series,D_final) # assign dataframe to var with series name
}

## Merge all log returns of series with full obs (4528) into one dataframe

all_ret <- merge(x=eval(parse(text=df_files$series[2])), #series[1] has only 293
                 y=eval(parse(text=df_files$series[3])),
                 by="date",all=TRUE) # start by merging second and third columns
for (series in df_files$series[4:length(df_files$series)]){ # merge other columns to all_ret
  if (dim(eval(parse(text=series)))[1] == 4528) {
    all_ret <- merge(x=all_ret, y=eval(parse(text=series)), by="date",all=TRUE)
  }
}

#############################################################################################
## COMPARE ALGORITHM FOR DIFFERENT VALUES
#############################################################################################

m_vals <- c(10,50,100,200,500,1000)
tau_vals <- c(10,50,100,200,500,1000)
c_vals <- c(0,1000,2000)

max(m_vals)+max(tau_vals)+max(c_vals)<dim(Rk)[1] # Check whether nothing runs out of bound

comp_mat <- matrix(rep(0,6*length(m_vals)*length(tau_vals)*length(c_vals)),
                      length(m_vals)*length(tau_vals)*length(c_vals),6)

colnames(comp_mat) <- c("m","tau","c","OBP","EWP","OBP-EWP")

counter <- 1
for (i in m_vals){
  for (j in tau_vals){
    for (k in c_vals){
      comp_mat[counter,1] <- i
      comp_mat[counter,2] <- j
      comp_mat[counter,3] <- k
      comp_mat[counter,4] <- OBP_Maker(m=i,tau=j, c=k)[3][[1]]
      comp_mat[counter,5] <- EWP_Maker(m=i,tau=j, c=k)[3][[1]]
      comp_mat[counter,6] <- comp_mat[counter,4] - comp_mat[counter,5]
      counter <- counter + 1
      print(counter)
    }
  }
}
print(comp_mat)

write.csv(comp_mat, file="/home/jan/Documents/Courses/T2_SMO/Project/MultiArmedBandit/Results/comp_mat.csv")
print(mean(comp_mat[,6]))
tail(comp_mat[,6])

###################################################################################
## Functions
###################################################################################

JS_shrink_est_mean <- function(series, average_all, eigenvalues){
  # In:
    # series (vec): series of interest 
    # average_all (vec): average of all series at each point in time
    # eigenvalues (vec): all eigenvalues of the matrix of series over time
  # Out:
    # JS_mean (double): mean of series using JS shrinkage estimation
  a <- (1/length(series)) * (sum(eigenvalues) - 2*max(eigenvalues))/
    (t(series-average_all)%*%(series-average_all))
  JS_series <- (1-a[1,1])*series + a[1,1]*average_all
  JS_mean <- mean(JS_series)
  return(JS_mean)
}

OBP_Maker <- function(m=100, l=4, dt=1, 
                      Rk=as.data.frame(all_ret[2:dim(all_ret)[1],2:dim(all_ret)[2]]),
                      tau=((dim(all_ret)[2]-1)-m), c=0){
  # IN
  ##### m (int): number of evaluations (h steps ahead)
  ##### l (int): number of factors, i.e. systemic movement in the market (usually 3-5)
  ##### dt (int): <- length of the time interval (currently not used)
  ##### Rk (matrix): gross return matrix
  ##### tau (int): size of training data
  ##### c (int): starting point in time series 
  # OUT
  ##### w_all (matrix): weightings at each point in time k
  ##### mu_all (vector): average return at each point in time k
  ##### cumret_OBP (double): comulative average after m time periods
  
  ## Initialize temporary matrices that are updated in loop
  k_i <- rep(0,dim(Rk)[2]) # chosen portfolios frequency
  
  ## Prepare result matrices
  w_all <- matrix(rep(0,dim(Rk)[2]*m),dim(Rk)[2],m)  # initialize weights (w)
  mu_all <- rep(0,m) # initially average returns (mu)
  
  for (k in (tau+1+c):(tau+m+c)){
    ## 1. Estimate average return and covariance matrix
    ERk <- apply(X=Rk[(k-tau):(k-1),], MARGIN = 2,mean) # expected returns of series
    # ERk <- apply(X=Rk[(k-tau):(k-1),], MARGIN = 2,mean) # expected returns of series
    Sk <- cov(Rk[(k-tau):(k-1),], use="complete.obs") # covariance matrix
    
    ## 2. Implement principal component decomposition
    eigen_obj <- eigen(x=Sk) # calculate eigenvalues and -vectors
    Lk <- diag(eigen_obj$values) # diagonal matrix with eigenvalues in decr. order
    Hk <- eigen_obj$vectors # orthogonal matrix with eigenvectors in columns
    
    ## 3. Compute renormalized similarity matrices and eigenvalues
    Hk_colsums <- apply(X=Hk, MARGIN = 2, FUN = sum) # calculate col sums
    Hk_t <- t(t(Hk)/Hk_colsums) # normalize each column to calculate Hk tilde
    Sk_t <- (t(Hk_t)%*%Sk%*%Hk_t) # sigma tilde (renormalized similarity matrices)
    Lk_t <- Sk_t # lambda tilde (renormalized eigenvalues)        
    
    ## 4. Compute Sharpe ratio of each arm
    SRk <- (t(Hk)%*%(ERk))/sqrt(eigen_obj$values)  # compute expected retuen over sqrt(eigval) ratio
    
    ## 5. Compute adjusted reward function of each arm
    reward_func_vals <- SRk + sqrt((2*log(k+tau))/(tau+k_i)) # calculalte reward
    
    ## 6. Select optimal arms based on adjusted reward function
    i_opt_l <- which.max(reward_func_vals[1:l]) # find index of maximum value in first l arms
    i_opt_nl <- which.max(reward_func_vals[l+1:length(reward_func_vals)]) # find index of maximum 
    # value in last n-l arms
    k_i[i_opt_l] <- k_i[i_opt_l] + 1  # update selection counter to check how often i was selected
    k_i[i_opt_nl] <- k_i[i_opt_nl] + 1  # update selection counter to check how often i was selected
    
    ## 7. Compute optimal mixture weight (eq 12)
    theta_k <- Sk_t[i_opt_l,i_opt_l]/(Sk_t[i_opt_l,i_opt_l] + Sk_t[i_opt_nl,i_opt_nl])
    
    ## 8. Compute optimal portfolio weight
    w_k <- (1-theta_k) * Hk_t[,i_opt_l] + theta_k * Hk_t[,i_opt_nl]
    
    ## Return: weight vectors and portfolio returns
    w_all[,k-(tau+c)] <- w_k
    mu_all[k-(tau+c)] <- (as.matrix(Rk[k,], nrows = 1)  %*% w_k) - 1
    
  }
  
  cumret_OBP <- 1
  cumret_OBP_series <- rep(0,length(mu_all))
  for (i in 1: length(mu_all)){
    cumret_OBP <- cumret_OBP * (1+mu_all[i])
    cumret_OBP_series[i] <- cumret_OBP
  }
  cumret_OBP
  return(list(w_all, mu_all,cumret_OBP,cumret_OBP_series))
}

EWP_Maker <- function(m=100, l=4, dt=1, 
                      Rk=as.data.frame(all_ret[2:dim(all_ret)[1],2:dim(all_ret)[2]]),
                      tau=((dim(all_ret)[2]-1)-m), c=0){
  # IN
  ##### m (int): number of evaluations (h steps ahead)
  ##### l (int): number of factors, i.e. systemic movement in the market (usually 3-5)
  ##### dt (int): <- length of the time interval (currently not used)
  ##### Rk (matrix): gross return matrix
  ##### tau (int): size of training data
  ##### c (int): starting point in time series 
  # OUT
  ##### w_all (matrix): weightings at each point in time k
  ##### mu_ewp (vector): average return at each point in time k
  ##### cumret_EWP (double): cumulative average after m time periods
  
  mu_ewp <- rep(0,m) # initially average returns (mu)
  w_all <- rep(1/dim(Rk)[2],dim(Rk)[2])
  for (k in (tau+1+c):(tau+m+c)){
    ERk <- apply(X=Rk[(k-tau):(k-1),], MARGIN = 2,mean) # expected returns of series
    mu_ewp[k-(tau+c)] <- mean(ERk)
  }
  cumret_EWP <- 1
  cumret_EWP_series <- rep(0,length(mu_ewp))
  for (i in 1: length(mu_ewp)){
    cumret_EWP <- cumret_EWP * (mu_ewp[i])
    cumret_EWP_series[i] <- cumret_EWP
  }
  cumret_EWP
  return(list(w_all, mu_ewp,cumret_EWP, cumret_EWP_series))
}
  
EWP_Maker()

###################################################################################
## Graphics
###################################################################################

tn <- 1000
cn <- 0

plot(y= EWP_Maker(m=1000, tau= tn, c=cn)[4][[1]], x=D[(cn+1):(cn+1000),1],
     type="l",col = "red", ylab = "Cumulative Return", xlab = "",
     ylim = c(-1,3), main = "Results for different choices of learning rates (tau)")
lines(y = OBP_Maker(m=1000,tau=1000, c=cn, l=1)[4][[1]], 
      x=D[(cn+1):(cn+1000),1], col="#252525")
lines(y = OBP_Maker(m=1000,tau=500, c=cn, l=1)[4][[1]], 
      x=D[(cn+1):(cn+1000),1], col="#525252")
lines(y = OBP_Maker(m=1000,tau=100, c=cn, l=1)[4][[1]], 
      x=D[(cn+1):(cn+1000),1], col="#737373")
legend("topleft",
       legend = c("EWP", "OPB tau=1000", "OPB tau=500", "OPB tau=100"),
       fill = c("red","#252525", "#525252", "#737373"))

plot(y= EWP_Maker(m=1000, tau=500, c=cn)[4][[1]], x=D[(cn+1):(cn+1000),1],
     type="l",col = "red", ylab = "Cumulative Return", xlab = "",
     ylim = c(-1,3), main = "Results for different sizes of passive portfolios (l)")
lines(y = OBP_Maker(m=1000,tau=500, c=cn, l=1)[4][[1]], 
      x=D[(cn+1):(cn+1000),1], col="#252525")
lines(y = OBP_Maker(m=1000,tau=500, c=cn, l=2)[4][[1]], 
      x=D[(cn+1):(cn+1000),1], col="#525252")
lines(y = OBP_Maker(m=1000,tau=500, c=cn, l=3)[4][[1]], 
      x=D[(cn+1):(cn+1000),1], col="#737373")
lines(y = OBP_Maker(m=1000,tau=500, c=cn, l=4)[4][[1]], 
      x=D[(cn+1):(cn+1000),1], col="#d9d9d9")
legend("topleft",
       legend = c("EWP", "OPB l=1", "OPB l=2", "OPB l=3", "OPB l=4"),
       fill = c("red","#252525", "#525252", "#737373", "#d9d9d9"))

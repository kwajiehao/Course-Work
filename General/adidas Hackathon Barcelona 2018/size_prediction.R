directory <- "/Users/kwajiehao/Documents/GitHub/sup/3/Miscellaneous/adidas Hackathon/"
file <- "dataset.csv"
df <- read.csv(paste0(directory,file))
df <- subset(df, select=-c(shoulder,chest,waist,hip, X, Unnamed..0, unique_customer,shirt, XS, S, M, L, XL))
str(df)

df$pref <- factor(df$pref, levels=c(1,2,3), labels=c('tight','regular','loose'))
names(df)

# introduce 2.5% measurement error to body measurements

for (i in 1:nrow(df)){
  for (j in 1:4){
    df[i,22+j] = runif(1,0.975, 1.025)*df[i,22+j]
  }
}

names(df)[23:length(names(df))] <- c("shoulder","chest","waist","hip")

write.csv(df, file="dataset_sizes.csv")

# multinomial logistic regression
library(caret)
library(foreign)
library(nnet)

# train test split
indices.training <- createDataPartition(y = df$label, p = 0.7, list = FALSE)
training <- df[indices.training[,1],]
test  <- df[-indices.training[,1],]

# check class split
a <- summary(df$label)
b <- summary(training$label)
c <- summary(test$label)

a <- a/sum(a)
b <- b/sum(b)
c <- c/sum(c)

# fit labels
fit <- multinom(label~. , data = training)
labels <- predict(fit,test,type="probs")

# compare with test labels
pred <- as.vector(apply(labels, MARGIN = 1, FUN = which.max))
test_labels <- rep(0, length(pred))

for (i in 1:length(pred)){
  if (test$label[i]=="L"){
    test_labels[i] = 1
  }
  if (test$label[i]=="M"){
    test_labels[i] = 2
  }
  if (test$label[i]=="S"){
    test_labels[i] = 3
  }
  if (test$label[i]=="XL"){
    test_labels[i] = 4
  }
  if (test$label[i] =="XS"){
    test_labels[i] = 5
  }
}

mean(pred==test_labels) # without measurement error, 97%, with, we have around 83%

# bar graph
library(reshape2)
df2 <- data.frame(Match = c(95,88), Mismatch = c(5,12), id=c("Without error","With measurement error"))
df_2 <- melt(df2, id.vars = "id")

plot <-ggplot(df_2, aes(x=id, y=value, fill = variable)) + 
  geom_bar(stat="Identity", position="dodge") + 
  labs(y="percentage", title = "Accuracy of prediction", x="") + 
  scale_y_continuous(limits = c(0,100)) +
  scale_fill_discrete(name = "Outcome") 

plot

which(pred!=test_labels)

# knn
library(FNN)
names(df)

knn_train <- training[,c(1:21, 23:26)]
dummies_train <- model.matrix(~ pref-1,data=knn_train)
knn_train <- cbind(dummies_train[,1:2], knn_train[,c(2:25)])

knn_classes <- training[,22]
knn_test <- test[,c(1:21, 23:26)]
dummies_test <- model.matrix(~ pref-1,data=knn_test)
knn_test <- cbind(dummies_test[,1:2], knn_test[,c(2:25)])

a <- knn(knn_train, knn_test, knn_classes, k=2)
mean(a == test$label)

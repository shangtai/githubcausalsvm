library(kernlab)
source("causalsvm.r")

grid_flag = TRUE


causalsvm_example<- function(filename){
    # This code take in a csv file
    # The first column is a binary indicator, 1 if it is assigned to treatment group
    # The second column is y observed.
    # The third column is rnd
    # The rest are the features data.

    raw_Data <- read.csv(filename, header= TRUE)

    # belows are some examples of ways to call the function

    #vanilla <- vanilladot()
    vanilla <- polydot(degree = 1, scale = 0.1, offset = 0)
    rbf <- rbfdot(sigma = 0.1)
    poly2 <- polydot(degree = 2, scale = 0.025, offset = 0)
    poly3 <- polydot(degree = 3, scale = 0.025, offset = 0)
    
    gamma_list <- 10^(-10:1)#c(0.00001, 0.0001, 0.001,0.01,0.1,1)#, seq(0.01,0.1,0.01))
    scale_list <- c(0.1)
    
    #split half the data as training and testing data separately
    smp_size <- floor(0.5* nrow(raw_Data))
    train_ind <- sample(seq_len(nrow(raw_Data)), size= smp_size)
    train_Data <- raw_Data[train_ind, ]
    test_Data <-  raw_Data[-train_ind, ]

    # first I will to sort the data so that treatment data is at the top
    # control data is at the bottom. this is to facilitate optimization coding during training phase
    train_Data <- train_Data[order(-train_Data[,1]),]

    treatment_vector_train <- train_Data[,1]
    y_obs_train <- train_Data[,2]
    rnd_train <- train_Data[,3]
    X_train <- train_Data[,4:ncol(raw_Data)]
    X_train<-as.data.frame(X_train)

    treatment_vector_test <- test_Data[,1]
    y_obs_test <- test_Data[,2]
    rnd_test <- test_Data[,3]
    X_test <- test_Data[,4:ncol(test_Data)]
    X_test<-as.data.frame(X_test)
    
    output <- caulsalsvm_predict(treatment_vector_train, y_obs_train, X_train, X_test, rnd_train, rbf, 10^(-6))
    
    return(output)
    #print(predicted)
 }

predicted <- causalsvm_example("langspiralwithnoise.csv")
print(predicted)

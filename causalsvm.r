library(gurobi)
# this software depends on gurobi optimization package
# this software implements the causal SVM paper described in (link)
# scroll down to the very end to see how to typically call such a function
# to call function within this file, setwd() to a directory containing this file and type source('causalsvm.r') at the top of the code

hinge_loss <- function(x){
   # this function returns the hinge loss function
   return(max(1+x,0))
}

evaluate_objective_directly <- function(y_treatment_train, y_control_train, predicted_h, treatment_vector_train, rnd_train){
    # given predicted_h, the predicted value this function evaluates the surrogate hinge loss function directly
    # y_treatment_train is a n x 1 vector, consist of 1 (effective) or -1 (not effective) if it is in the training group, and otherwise hte value doesn't matter
    # y_control_train is an n x 1 vector, consists of 1 (effective) or -1 (not effective) if it is in the training group, and otherwise hte value doesn't matter
    # treatment_vector_train consists of 1 and -1, it it is 1, the correponding index belongs to the treatment group and -1 for the control group
    # rnd_train is the corresponding radon-nykodym derivative that maps the control group to the treatment group

    # treatment_index is a vector that stores indices of the treatment data point
    treatment_index <- which(treatment_vector_train == 1)

    # control_index is a vector that stores indices of the control data point
    control_index <- which(treatment_vector_train == -1)

    # initialize variable to 0
    treatment_average <- 0
    control_average <- 0

    # compute the number of treatment data points and control data points
    treatment_length <- length(treatment_index)
    control_length <- length(control_index)

    # going through each data points, to compute the weighted hinge loss. lapply can be used to speed up the code when necessary.
    for (i in seq(1:length(treatment_vector_train))){
       if (treatment_vector_train[i] == 1){
           treatment_average <- treatment_average + hinge_loss(-predicted_h[i]*y_treatment_train[i])
       } else{
           control_average <- control_average + rnd_train[i] * hinge_loss(predicted_h[i]* y_control_train[i])   
       }
    }

    #average the number of points
    treatment_average <- treatment_average / treatment_length
    control_average <- control_average / control_length
    return(max(treatment_average, control_average))
}

performanceZeroOneLoss <- function(yTreatment, yControl, predictedDiff){
    # given predictedDiff, this function evealuates the quality of the prediction based on 0-1 loss function
    # all the inputs are lists of length n, the number of sample data
    # yTreatment stores the ground truth value if this unit is being assigned to the treatment group. It takes value 1 or -1.
    # yTreatment stores the ground truth value if this unit is being assigned to the treatment group. It takes value 1 or -1.
    # this function is only useful in simulations when the ground truth is known.
    n <- length(yTreatment)
    loss <- 0
    for (i in 1:n){
       if (yTreatment[i] > yControl[i]){
           if (predictedDiff[i] <= 0){
               loss <- loss + 1
           }
       } else if (yTreatment[i] < yControl[i]){
           if (predictedDiff[i] >= 0){
               loss <- loss+1
           }
       } else{
          if (abs(predictedDiff[i]) >= 1){
              loss <- loss +1
           }
       }
}
loss <- loss/n
return(loss)
}

performanceThetaLoss <- function(yTreatment, yControl, predictedDiff){
    # this function return the theta-loss function when the theta value is chosen to be 0.01, 0.1, and 0.5, here theta will be mapped to the quantile
    # such threshold is set in the variable percentile_vector 
    # this function is designed to handle the problem when 1 is not the natural threshold.
    n <- length(yTreatment)
    percentile_vector <- c(0.01,0.1,0.5)

    # map theta to the quantile of the absolute value of the predicted difference
    theta <- quantile(abs(predictedDiff), percentile_vector)

    # initialize the loss vector
    loss <-rep(0, length(percentile_vector))

    # setting various threshold
    for (j in c(1: length(percentile_vector))){
        #going through each data point to update the loss, lapply can help to speed things up here
        for (i in 1:n){
           if (yTreatment[i] > yControl[i]){
               if (predictedDiff[i] <= - theta[j]){
                   loss[j] <- loss[j] + 1
               }
            } else if (yTreatment[i] < yControl[i]){
               if (predictedDiff[i] >= theta[j]){
                   loss[j] <- loss[j] + 1
               }
           } else{
               if (abs(predictedDiff[i]) >= theta[j]){
                   loss[j] <- loss[j] +1
               }
           }
        }
        loss[j] <- loss[j]/n
    }
    return(loss)
}


compute_projected_diff <- function(predictedDiff){
    # this function performs truncation at 1 and -1
    # if a value is positive, defined as >= 1, then it is mapped to 2
    # if a value is negative, defined as <= -1, it is mapped to -2
    # if a value is neutral, it is mapped to 0
    # it returns a vector of 0, 2, and -2

    n <- length(predictedDiff)
    projected_diff <- rep(0,n)
    for (i in 1:n){
       if (predictedDiff[i]>=1){
           projected_diff[i] <- 2
       }
       else if (predictedDiff[i] <= -1){
           projected_diff[i]<- -2
       }
    }
    return(projected_diff)
}


caulsalsvm_predict <- function(treatment_vector_train, y_obs_train, X_train, X_test, rnd_train, kernelused, gamma){
   # this is the main function that we want to call.
   # treatment_vector_train records if a sample point is a treatment or control data_point
   # while the code is written such that it is possible to run simulation on it, the prediction function does not depend on the data that is not available in practice
   # y_obs_train is the y_obs of the training data
   # X_train is the matrix of feature data
   # X_test is the feature for the data that we will test our model on.
   # rnd_train is the radon-nikodym derivative for the training data where we map the control group to the treatment group
   # kernelused is a kernel object
   # gamma is the regularizer for the causal SVM model to prevent overfitting.
   
   # the following two vectors store the indices for the treatment and control training data respectively.
   treatment_vector_index <- which(treatment_vector_train == 1)
   control_vector_index <- which(treatment_vector_train != 1)

   # we append the index together for the convenience of formulating the optimization problem. We then use this to rearrange our training data such that the treatment group and control group are not mixed.
   sorted_index <- c(treatment_vector_index, control_vector_index)
   treatment_vector_train <- treatment_vector_train[sorted_index]
   y_obs_train <- y_obs_train[sorted_index]
   X_train <- X_train[sorted_index,]
    
   # we compute the number of treatment data
   n_treatment <- length(treatment_vector_train[treatment_vector_train == 1])

   # we compute the total number of training data
   length_train_ind <- length(y_obs_train)

   # we compute the number of control data
   n_control <- length_train_ind - n_treatment

   # we create a variable y_perturb_obs that copies y_obs_train but change the sign if it is from the control group
   y_perturb_obs <- y_obs_train
   y_perturb_obs[(n_treatment+1):(n_treatment+n_control)]<- -y_perturb_obs[(n_treatment+1):(n_treatment+n_control)]
   
   # prepare the signed Gram matrix, linear algebra trick is used to map the right sign to each component
   K <- kernelMatrix(kernel=kernelused, as.matrix(X_train))
   K <- (y_perturb_obs %o% y_perturb_obs) * K
   
   # preparing the dual optimization problem in Gurobi
   model <- list()
   model$A <- matrix(0, length_train_ind+2, length_train_ind+2)
   model$A[1:length_train_ind, 1: length_train_ind] <- diag(length_train_ind)
   model$A[1:n_treatment, (length_train_ind+1)] <- -1/n_treatment
   model$A[(n_treatment+1):length_train_ind, (length_train_ind+2)] <- (-1/n_control)*(rnd_train[(n_treatment+1):length_train_ind])
   model$A[length_train_ind+1, (length_train_ind+1):(length_train_ind+2)] <- c(1,1)
   model$A[length_train_ind+2,1:length_train_ind] <- t(y_perturb_obs)

   model$obj <- rep(0, ncol(model$A))
   model$obj[1:length_train_ind] <- rep(1, length_train_ind)

   model$sense<-rep('<=',(length_train_ind+2))
   model$sense[(length_train_ind+1):(length_train_ind+2)] <- c('=','=')

   model$rhs <- rep(0,nrow(model$A))
   model$rhs[length_train_ind+1] <- 1

   model$modelsense <- 'max'

   model$Q <- matrix(0,ncol(model$A),ncol(model$A))
   model$Q[1:length_train_ind,1:length_train_ind] <- -K/(4*gamma)

   # perturb the matrix to avoid singularity
   model$Q<-model$Q-0.00001*diag(ncol(model$A))

   params <- list(Method=2, ResultFile='model.mps')

   # we solve the optimization problem by calling the solver
   result <- gurobi(model, params)

   # extract the maximization optimal solution
   alpha <- result$x[ncol(model$A)-1]
   beta <- result$x[ncol(model$A)]
   solutionVector <- result$x

   # obtain the dual variable, i.e. a solution to the primal problem
   dualFromGurobi <- result$pi
   # extract the numerical intercept, now we have sufficient information to perform prediction.
   w0 <-dualFromGurobi[length(dualFromGurobi)]
   # perform prediction for the train data
    #predicted_train <- w0+solutionVector[1:length_train_ind]%*%sweep(kernelMatrix(kernel= kernelused, as.matrix(X_train),as.matrix(X_train)), MARGIN=1, y_perturb_obs,"*")/(2*gamma) 
   # perform prediction for the test data
   predicted_test <- w0+solutionVector[1:length_train_ind]%*%sweep(kernelMatrix(kernel=kernelused, as.matrix(X_train),as.matrix(X_test)), MARGIN=1, y_perturb_obs,"*")/(2*gamma)

   return(predicted_test)
}
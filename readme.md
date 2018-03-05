Supplement for "A Minimax Surrogate Loss Approach to Conditional Difference Estimation"
=================================================================================================


We present a new machine learning approach to estimate whether a treatment has an effect on an individual, in the setting of the classical potential outcomes framework with binary outcomes. To overcome the problem that both treatment and control outcomes for the same unit are required for supervised learning, we propose surrogate loss functions that incorporate both treatment and control data. The new surrogates yield tighter bounds than the sum of losses for treatment and control groups. A specific choice of loss function, namely a type of hinge loss, yields a minimax support vector machine formulation. The resulting optimization problem requires the solution to only a single convex optimization problem,
incorporating both treatment and control units, and it enables the kernel trick to be used to handle nonlinear (also non-parametric) estimation. Statistical learning bounds are also presented for the framework, and experimental results.

Assuming the label to take 1 (effective) and -1 (not effective) value only, the possible treatment effects are 0 (neutral) , 2 (effective), and -2 (not effective). Our method uses an SVM method, possibly with RBF kernel. The prediction outcome take continous value. 

___________________________
**To run the algorithms**

We provided the R implementation of the code. The current code depends on the Gurobi solver. If one wishes to write their own code, a quadratic programming solver which returns the dual result is recommended. 

The main file is causalsvm.r and the main function is caulsalsvm_predict() which take in the following argument.

*  the first input, treatment_vector_train: It is a vector consisting of 1 and -1. it takes value 1 if it is in the treatment group and -1 if it in the control group.
*  the second input, y_obs_train: It is a vector consisting of 1 and -1. 
*  the third input, X_train, is the feature matrix for training data set. Each row is a sample and each column represents a feature.
*  the fourth input, X_test, is the feature matrix for the test data set. Each row is a sample and each column represents a feature. 
*  the fifth input, rnd_train, is the radon-nikodym derivative for the training data. One possible algorithm to obtain this is unconstrained Least-Squares Importance Fitting (uLSIF) which is included in densratio package in R software.
*  the sisth input is the kernel used. Examples are included in the example file to illustrate how can one construct the kernel object.
*  the seventh input is gamma. This is a parameter in the causal SVM method which regularize the objective function. 
caulsalsvm_predict <- function(treatment_vector_train, y_obs_train, X_train, X_test, rnd_train, kernelused, gamma, filename = "", print_option = FALSE){

An example has been provided in the example.r file to illustrate how the function causalsvm_predict can be called. To run the file, just type source('example.r').
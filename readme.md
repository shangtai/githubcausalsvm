Supplement for "A Minimax Surrogate Loss Approach to Conditional Difference Estimation"
=================================================================================================


We present a new machine learning approach to estimate whether a treatment has an effect on an individual, in the setting of the classical potential outcomes framework with binary outcomes. To overcome the problem that both treatment and control outcomes for the same unit are required for supervised learning, we propose surrogate loss functions that incorporate both treatment and control data. The new surrogates yield tighter bounds than the sum of losses for treatment and control groups. A specific choice of loss function, namely a type of hinge loss, yields a minimax support vector machine formulation. The resulting optimization problem requires the solution to only a single convex optimization problem,
incorporating both treatment and control units, and it enables the kernel trick to be used to handle nonlinear (also non-parametric) estimation. 

The label takes values 1 (effective) and -1 (not effective). The possible treatment effects are 0 (neutral) , 2 (effective), and -2 (not effective). Our method uses an SVM-based method, possibly with an RBF kernel. The prediction outcomes take continuous value. 

___________________________
**To run the algorithm**

We provide an R implementation of the algorithm. The current code uses the Gurobi solver. If one wishes to write their own code, a quadratic programming solver which returns the dual result is recommended. 

The main file is causalsvm.r and the main function is caulsalsvm_predict() which takes in the following arguments.

*  the first input, treatment_vector_train: It is a vector consisting of 1 and -1. It takes value 1 if it is in the treatment group and -1 if it is in the control group.
*  the second input, y_obs_train: It is a vector consisting of 1 and -1. 
*  the third input, X_train, is the feature matrix for the training data set. Each row is a sample and each column represents a feature.
*  the fourth input, X_test, is the feature matrix for the test data set. Each row is a sample and each column represents a feature. 
*  the fifth input, rnd_train, is the Radon-Nikodym derivative for the training data. One possible algorithm to obtain this  quantity is the unconstrained Least-Squares Importance Fitting (uLSIF) method which has been implemented in densratio package in the R software.
*  the sixth input is the kernel used. Examples are included in the example file to illustrate how can one construct the kernel object.
*  the seventh input is gamma. This is a parameter in the causal SVM method which regularizes the objective function. 

An example has been provided in the example.r file to illustrate how the function causalsvm_predict can be called. To run the file, just type source('example.r').
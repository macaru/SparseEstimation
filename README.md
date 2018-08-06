This SAS program is the macro to obtain one-step estimator that was suggested by Zou and Li (2008) using SAS v9.4. 
The program include 4 macros as follows. 
The macro "ose_scad" returns one-step estimator that was suggested by Zou and Li (2008).
The macro "lasso_shot" returns LASSO estimator using shooting algorithm that was suggested by Fu (1998).
The macros "diff" and "cv"returns the value of K folds cross validation (CV) for estimated coefficients.
The macro "optlam" searchs the optimal smoothing parameter based on K folds CV and return the optimal one-step estimates.
 
Example for usage of macro: %optlam(dat, y, col1 - col12, 0.1, 1, 30);

Here, dat is the name of dataset, y is the name of response variable, col1 - col12 are the names of predict variables.
0.1, 1 and 30 in %optlam macro means the min of smoothing parameter, max of smoothing parameter for grid search and the length of grid of smoothing parameter respectively. Implementing the macro %optlam, you could get 3 datasets optcoef, optlambda and seqcoef. The dataset optcoef is estimated coefficient, the dataset seqcoef and optlambda show the tarnsitions of coefficients and values of K folds CV when lambda is increased, respectively.
 
[References]

[1] Fu (1998). Penalized regression: The bridge versus the lasso.
   Journal of Computational and Graphical Statistics 7: 397 -416.

[2] Zou and Li (2008). One-step sparse estimates in nonconacve penalized likelihood models.
   The Annals of Statistics 36: 267-288.
 
Last changed: 06 AUG 2018 by Masaru Kanba

# SAS-BCAR: Spatial Adaptive Selection using Binary Conditional Autoregressive Model with Application to Brain-Computer Interface



![radius_simulation](https://github.com/zikaiLin/sas-bcar/blob/main/radius_simulation.png)



## Abstarct

In medical imaging studies, statistical inferences on scalar-on-image regression are challenging due to the limited sample size and the high-dimensionality of datasets. Also, the imaging predictors often exhibit spatially heterogeneous activation patterns and the complex nonlinear associations with the response variable. To address these issues, we propose a novel prior model for Bayesian scalar-on-image regression, the Spatial Adaptive Selection using Binary Conditional Autoregressive Model (SAS-BCAR) prior.

The proposed SAS-BCAR prior employs a binary conditional autoregressive model to address spatial dependencies among feature selection indicators, effectively identifying spatially structured sparsity patterns in image datasets. Moreover, SAS-BCAR allows for adaptive feature selection to the varying spatial dependencies across different image regions, leading to a more precise and robust feature selection process for image analysis.
We have developed an efficient posterior computation algorithm for SAS-BCAR and  demonstrated the advantages of SAS-BCAR in image classification tasks with the benchmark computer vision datasets and in analysis of electroencephalography data in Brain computer interface applications.



## Code Example



Here's an 

```R
  ## ----load data------##
  i = 1
  mnist_data_training = readRDS(sprintf("./training_150_%i.rds", i))
  training_y = as.integer(mnist_data_training$mnist_data_training_Y_sub == 7)
  training_X = matrix(NA, nrow = nrow(mnist_data_training$mnist_data_training_X_sub), ncol = 28*28)
  for(s in 1:nrow(mnist_data_training$mnist_data_training_X_sub)){
    training_X[s,] = c(mnist_data_training$mnist_data_training_X_sub[s,,])
  }

```



Normalize the imaging data:

```R
  training_X = training_X/255
```



Fit the model:

```R
sasbcar_model <- sasbcar::fit_model(X = training_X, y = training_y)
```



Prediction:

```R
mnist_data_testing = readRDS("testing.rds")
testing_X = matrix(NA, nrow = nrow(mnist_data_testing$mnist_data_testing_X), ncol = 28*28)
testing_y = as.integer(mnist_data_testing$mnist_data_testing_Y == 9)

for(s in 1:nrow(mnist_data_testing$mnist_data_testing_X)){
  testing_X[s,] = c(mnist_data_testing$mnist_data_testing_X[s,,])
}
testing_X = testing_X/255


prediction <- sasbcar::predict_model(model = sasbcar_model, new_data = testing_X)
```


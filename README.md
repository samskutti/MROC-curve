# MROC-curve
Multivariate Receiver Operating Characteristic Curve

This app demonstrates the Multivariate Receiver Operating Characteristic (MROC) model.
Various types of MROC models can be seen: full model, stepwise model and custom model
Full model includes all the variables in the model
Stepwise model is a variable selection procedure that identifies the significant variables and allows them into the model
Custom model allows the user to pick and choose the variables that they would like to make a part of the model.
The measures of the MROC curve like AUC, sensitivity, Specificity are displayed along with their confidence intervals. The measures AUC and sensitivity are tested at against 50% to see if the model is better than random classification
The second tab showcases the coefficients of the model along with the model precision
The third tab displays the MROC curve for the considered model. The closer the curve to the top left corner, better the classification.

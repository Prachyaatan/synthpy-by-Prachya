# Synthpy Package
Created by Prachya Tantivatana.

# Description: 
This package provide function for computation and analysis of synthetic control method using framework laid out by Abadie, Diamond & Haimueller (2015). Synthetic control implemented in this package is used for measuring a treatment effect on a particular outcome variable of a single treated unit using a weight-averaged of a the same variable of control units to create a "counter-factual" of the treated unit when the treatment did not occur. The weight assigned to control is sum up to one and is between (0,1). Characteristic variable can also be included as predictor for the computation of weight. Additionally, this package include visualization function, computation of confidence interval using leave-one-out and hypothesistest function from placebo treatment.

# Details
Suppose that the outcome variable of the treated unit at time t is $ Y_{1t} $ and for control unit it is $ Y_{0t} $. A synthetic control is a weighted average of the characteristics of control unit such that it best resemble the treated unit during pretreatment period. In other words, it is as follows,

  $$ W = argmin_w \Sigma_{m=1}^{k}v_m(X_{1m}-X_{0m}w)^2 $$
  
Where $ X_{1m} $ and $ X_{0m} $ denotes m "predictor" for treated and control unit respectively, the predictors is a scarlar that represent each characteristic of a particular unit which the outcome variable can be included as well. Normally, predictor should be normal accross all characteritic. $ v_m $ represents variable weight that give weights to each characteristic. Intutitively, the more important the variable, the more weight is given. To select  $ v_m $ , pretreatment data is split into training and validation period where the optimal set of V is the one that minimize the root mean square prediction error (RMSPE) in the validation set. The RMSPE is given as follows:

  $$ V = argmin_v  \left( \frac{1}{T_0}\Sigma_{t=1}^{T_0} \left( Y_{1t}- \Sigma_{j=1}^JW^{\`}Y_{0t}\right)^2 \right)^{1/2} $$

Notes:
This is part of UCLA Master of Quantitative Economics program (MQE), capstone project for winter 2023.
Capstone Title: The effect of Thailand 2014 regime changes: evidence using synthetic control.

Reference:
Abadie, A., Diamond, A., & Hainmueller, J. (2015). “Comparative Politics and the Synthetic Control” Method. American Journal of Political Science, 59(2), 495–510. DOI:10.2139/ssrn.1950298

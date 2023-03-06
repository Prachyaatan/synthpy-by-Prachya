# Synthpy Package
Created by Prachya Tantivatana.
Email: prachyaatan@gmail.com

# Description: 
This package provide function for computation and analysis of synthetic control method using framework laid out by Abadie, Diamond & Haimueller (2015). Synthetic control implemented in this package is used for measuring a treatment effect on a particular outcome variable of a single treated unit using a weight-averaged of a the same variable of control units to create a "counter-factual" of the treated unit when the treatment did not occur. The weight assigned to control is sum up to one and is between (0,1). Characteristic variable can also be included as predictor for the computation of weight. Additionally, this package include visualization function, computation of confidence interval using leave-one-out and hypothesistest function from placebo treatment.

# Details
Suppose that the outcome variable of the treated unit at time t is $Y_{1t}$ and for control unit it is $Y_{0t}$. A synthetic control is a weighted average of the characteristics of control unit such that it best resemble the treated unit during pretreatment period. In other words, it is as follows,

  $$ W = argmin_w \Sigma_{m=1}^{k}v_m(X_{1m}-X_{0m}w)^2 $$
  
Where $X_{1m}$ and $X_{0m}$ denotes m "predictor" for treated and control unit respectively, the predictors is a scarlar that represent each characteristic of a particular unit which the outcome variable can be included as well. The default setting for predictors is the mean, averging the variable over pretreatment period and these predictors are standardized accross all characteritic. $v_m$ represents variable weight that give weights to each characteristic. Intutitively, the more important the variable, the more weight is given. To select  $v_m$ , pretreatment data is split into training and validation period where the optimal set of V is the one that minimize the root mean square prediction error (RMSPE) in the validation set. Let the validation period indexed from 1 to $T_0$ , then RMSPE is given as follows:

  $$ V = argmin_v  \left( \frac{1}{T_0}\Sigma_{t=1}^{T_0} \left( Y_{1t}- \Sigma_{j=1}^JW^{\`}Y_{0t}\right)^2 \right)^{1/2} $$

## Technical details 
1) the computation of synthetic control is a nested optimization problem (with W, V), hence,  W is optimized using [cxvpy](https://www.cvxpy.org/) package and V is optimized using [scipy optimize](https://docs.scipy.org/doc/scipy/reference/optimize.html#module-scipy.optimize) function with Nelder-Mead as default method.

2) Confidencial interval is constructed by aggregating the treatment effect (difference between actual and synthetic) over all leave-on-out (LOO) samples, a baseline analysis but with one control unit dropped. With the LOO samples, its mean and standard error will be computed to create the interval with $\alpha$ confidence level

3) Two-tail interval for hypothesis test is constructed using a permutation of place which is created by assigning the treatment to one control unit in period prior to the actual treatment, (e.g. for 10 control units and 5 place year there will be 50 permutation of placebos). The treatment effect will then be tested against the interval from these placbos, intuititively, treatment effect is significant if it is statistically differents from the placebo effect.


````
````

**Notes**:
This is a part of UCLA Master of Quantitative Economics program (MQE), capstone project for winter 2023.
Capstone Title: The effect of Thailand 2014 regime changes: evidence using synthetic control.


**Reference**:
Abadie, A., Diamond, A., & Hainmueller, J. (2015). “Comparative Politics and the Synthetic Control” Method. American Journal of Political Science, 59(2), 495–510. DOI:10.2139/ssrn.1950298

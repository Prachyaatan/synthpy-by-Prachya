# Synthpy Package
Created by Prachya Tantivatana.

# Description: This package provide function for computation and analysis of synthetic control method using framework laid out by Abadie, Diamond & Haimueller (2015). Synthetic control implemented in this package is used for measuring a treatment effect on a particular outcome variable of a single treated unit using a weight-averaged of a the same variable of control units to create a "counter-factual" of the treated unit when the treatment did not occur. The weight assigned to control is sum up to one and is between (0,1). Characteristic variable can also be included as predictor for the computation of weight. Additionally, this package include visualization function, computation of confidence interval using leave-one-out and hypothesistest function from placebo treatment.

Notes:
This is part of UCLA Master of Quantitative Economics program (MQE), capstone project for winter 2023.
Capstone Title: The effect of Thailand 2014 regime changes: evidence using synthetic control.

Reference:
Abadie, A., Diamond, A., & Hainmueller, J. (2015). “Comparative Politics and the Synthetic Control” Method. American Journal of Political Science, 59(2), 495–510. DOI:10.2139/ssrn.1950298

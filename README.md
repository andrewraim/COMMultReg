The `COMMultReg` package demonstrates the Conway-Maxwell Multinomial (CMM)
distribution and regression model from the paper:

> Darcy Steeg Morris, Andrew M. Raim, and Kimberly F. Sellers (2020). A
> Conway-Maxwell-multinomial distribution for flexible modeling of clustered
> categorical data. Journal of Multivariate Analysis, 179:104651. 
> <https://doi.org/10.1016/j.jmva.2020.104651>.

To install the package from Github, run:

```
library(devtools)
install_github("andrewraim/COMMultReg")
```

Some examples are provided in the `inst` folder, including:

1. A Shiny app to visualize the Conway-Maxwell Binomial distribution.
1. A Shiny app to visualize trinomial CMM, as in Figure 2 of the paper.
1. Fitting CMM to the pollen dataset, as in section 6 of the paper.
1. Fitting CMM to the gator dataset, as in the supplement of the paper.


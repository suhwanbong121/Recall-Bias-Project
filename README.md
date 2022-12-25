Accounting for Recall Bias in Observational Studies
---
These are the R code files to apply the estimation methods presented in Bong, Lee and Dominici's *"Differential recall bias in self-reported risk factors in observational studies."*

---
#### Overview 

1. Simulation\
This part aims to replicate the simulation results in Table 2 in Section 5.

2. Application\
This part aims to replicate the results from the Wisconsin Longitudinal Study (WLS) dataset discussed in Section 6. We are unable to share the original data, but it can be downloaded at [here](https://www.ssc.wisc.edu/wlsresearch/data/downloads/). We include an R code file called *"WLS_to_dataset.R"* that can create the same dataset that we used in the paper. Figures 1, 2, 3, and Table 3 can be replicated with this dataset. 

---
#### Description of Directories

There are two directories *Rfunctions* and *Rcode*.

* *Rfunction* - contains one R script.
  * *"basic_functions.R"* - contains all the R functions that implement the estimation methods discussed in the paper. The maximum likelihood, stratification, blocking, and nearest neighborhood combination methods are implemented.

* *Rcode* - contains three R scripts.
  *	*"simulation.R"* - shows R codes of simulation study discussed in Section 5 that can replicate Table 2.
  * *"WLS_to_dataset.R"* - shows how to create a dataset used in Section 6 from the WLS data.
  *	*"WLS_analysis.R"* - contains all the statistical analysis. The two estimation methods are implemented here with the dataset obtained from "WLS_to_dataset.R". Also, we illustrate how readers can reproduce Figures 1, 2, 3, and Table 3.

---

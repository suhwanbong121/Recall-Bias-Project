Differential Recall Bias in Estimating Treatment Effects in Observational Studies
---
These are the R code files to apply the estimation methods presented in Bong, Lee and Dominici's *"Differential recall bias in estimating treatment effects in observational studies."*

---
#### Overview 

1. Simulation\
This part aims to replicate the simulation results in Table 2 in Section 5 and Table 4 in Supplementary Materials.

2. Application\
This part aims to replicate the results from the Wisconsin Longitudinal Study (WLS) dataset discussed in Section 6. We are unable to share the original data, but it can be downloaded at [here](https://www.ssc.wisc.edu/wlsresearch/data/downloads/). We include an R code file called *"WLS_to_dataset.R"* that can create the same dataset that we used in the paper. Figures 1, 2 and Table 3 can be replicated with this dataset. 

3. Supplementary Materials\
This part aims to replicate the results from the Wisconsin Longitudinal Study (WLS) dataset discussed in Supplementary Materials. We use the same data as in the Application section. Figures 1, 2, 3, and Tables 1, 2, and 3 in Supplementary Materials can be replicated with this dataset.

---
#### Description of Directories

There are two directories *Rfunctions* and *Rcode*.

* *Rfunction* - contains one R script.
  * *"basic_functions.R"* - contains all the R functions that implement the estimation methods discussed in the paper. The maximum likelihood, stratification, blocking, nearest neighborhood combination, and covariate balance adjustment methods are implemented. 

* *Rcode* - contains four R scripts.
  *	*"simulation.R"* - shows the R codes of the simulation study discussed in Section 5, which were used to generate the results presented in Table 2. We also show the R codes of double score stratification discussed in Supplementary Materials, which were used to generate the Table 4 in Supplementary Materials.
  * *"WLS_to_dataset.R"* - shows how to create a dataset used in Section 6 from the WLS data.
  *	*"WLS_analysis.R"* - contains all the statistical analysis provided in the main manuscript. The estimation methods are implemented here with the dataset obtained from "WLS_to_dataset.R". Also, we illustrate how readers can reproduce Figures 1, 2, and Table 3 in the main manuscript.
  * *"WLS_supplementary.R"* - contains all the statistical analysis provided in Supplementary Materials. We use the same dataset as in the Application section. We illustrate how readers can reproduce Figures 1, 2, 3, and Tables 1, 2, and 3 in Supplementary Materials.

---

# SHet_python

This is a python implementation of the SHet method proposed in Zhou et al. 2015. This is a meta-analysis-basd method for multivariate GWAS based on summary statistics from univariate tests.

## Getting Started

This repository contains:

The module SHet.py with two functions:
* Truncated: calculates the SHet statistic from the Wald summary statistics of the univariate tests and the correlation matrix between the statistics.
* EstimateGamma: estimates the empirical gamma distribution of SHet that can be used to calculate p-values.

The directory test_SHet with:
* 10 sample quantitative traits for 197 genetically variant samples (10_tierpsy_features.csv).
* Summary statistics from univariate GWAS run for these 10 traits for a given SNP (summary_stats_SNP=1:717.csv).
* The script test_SHet.py that will run the SHet method for this data.

The sample data come from C. elegans quantitative phenotyping.


### Prerequisites

You will need the following python packages:

```
conda install numpy pandas scipy time
```

## Running the tests

You can run the tests by downloading or cloning the data and running the script test_SHet.py.

You can then use your own data to run the functions in the SHet module.


## Authors

* **Eleni Minga** - *Initial work* - [em812](https://github.com/em812)

## Acknowledgments

This python implementation is based on the R implementation available at https://ctg.cncr.nl/software/multivariate_gwas. This R code is a companion to the review of Vroom at al. 2019 on multivariate GWAS methods.


## References
César-Reyer Vroom, Christiaan de Leeuw, Danielle Posthuma, Conor V. Dolan, Sophie van der Sluis, 2019, The more the merrier? Multivariate approaches to genome-wide association analysis. bioRxiv 610287; doi: https://doi.org/10.1101/610287

Zhu, X., Feng, T., Tayo, B.O., Liang, J., Young, J.H., Franceschini, N., Smith, J.A., Yanek, L.R., Sun, Y.V., Edwards, T.L. and Chen, W., 2015. Meta-analysis of correlated traits via summary statistics from GWASs with an application in hypertension. The American Journal of Human Genetics, 96(1), pp.21-36.
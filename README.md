# Multi-Attribute Subset Selection enables prediction of representative phenotypes across microbial populations

This repository serves as the codebase for our preprint: https://doi.org/10.1101/2022.06.20.496733.

Multi-Attribute Subset Selection (MASS) is an algorithm based on mixed-integer linear programming to optimally split a set of attributes (or variables) into predictor and response attributes. For this, the response attributes are modeled as linear combinations of predictor attributes.

We used this algorithm to identify the most descriptive growth conditions in microbial phenotyping experiments. This has implications for designing high-throughput phenotyping efforts as it can help to profile only the most informative experimental conditions and therefore reduces the experimental effort. We further observed, that the top predictor conditions selected by MASS comprise meaningful metabolic axis.

The general nature of the algorithm should make it applicable in other areas as well.

## Content
 - DATASET 1: marine bacteria (DEMO)
     - MASS: [MASS-marine.ipynb](./Code/MASS-marine.ipynb) (jupyter notebook)
     - RF: [RF_evalutation_w_baseline-marine.ipynb](./Code/RF_evalutation_w_baseline-marine.ipynb) (jupyter notebook)
 - DATASET 2: fermentation (BacDive1)
     - data download and preparation: [download_bacdive_traits.Rmd](./Code/download_bacdive_traits.Rmd) (Rmarkdown notebook)
     - MASS: [MASS-bacdive1.ipynb](./Code/MASS-bacdive1.ipynb) (jupyter notebook)
     - RF: [RF_evalutation_w_baseline-bacdive1.ipynb](./Code/RF_evalutation_w_baseline-bacdive1.ipynb) (jupyter notebook)
 - DATASET 3: yeast species
     - MASS: [MASS-yeast.ipynb](./Code/MASS-yeast.ipynb) (jupyter notebook)
     - RF: [RF_evalutation_w_baseline-yeast.ipynb](./Code/RF_evalutation_w_baseline-yeast.ipynb) (jupyter notebook)
 - manuscript figures: [mass_figures.Rmd](./Code/mass_figures.Rmd)  (Rmarkdown notebook)

## System Requirements
### For jupyter notebooks
```
Python (v3.10.6)

IPython (v8.12.0)
jupyter_core (v5.3.0)
gurobi (v9.1.1)
iterative-stratification (v0.1.6)
numpy (v1.19.2)
matplotlib (v3.3.2)
pandas (v1.1.3)
seaborn (v0.11.0)
scikit-learn (v0.23.2)
```

### For Rmarkdown notebooks
```
R (v4.2.1)

pheatmap (v1.0.12)
BacDive (v0.8.0) # installed using `install.packages("BacDive", repos="http://R-Forge.R-project.org")` as noted on https://r-forge.r-project.org/R/?group_id=1573
forcats (v0.5.2)
stringr (v1.5.0)
dplyr (v1.1.0)
purrr (v1.0.1)
readr (v2.1.2)
tidyr (v1.3.0)
tibble (v3.1.8)
ggplot2 (v3.4.1)
tidyverse (v1.3.2)
ggpubr (v0.6.0)
patchwork (v1.1.2)
viridis (v0.6.2)
viridisLite (v0.4.1)
```

## Installation
Gurobi can be downloaded and installed from: http://www.gurobi.com/downloads/gurobi-optimizer-eula/. The installation of Gurobi requires obtaining a license. Free academic license can be obtained from: https://www.gurobi.com/downloads/end-user-license-agreement-academic/.

Typically the installation should be accomplishable within less of a working day.

## Demo
The smallest dataset (DATASET 1) can serve as demo. Please refer to the respective jupyter notebooks to get idea how to run a MASS analysis and what output to expect. A MASS run on this dataset should not last more than a minute on a normal desktop computer.

MASS application on DATASET 2 and 3 are more resource intensive and can run a couple of days.

## Citation
If you find it useful, please cite our paper as follows:

```
@article{forchielli2022prediction,
  title={Prediction of representative phenotypes using Multi-Attribute Subset Selection (MASS)},
  author={Herbst, Konrad and Wang, Taiyao and Forchielli, Elena and Thommes, Meghan and Paschalidis, Ioannis Ch and Segre, Daniel},
  journal={bioRxiv},
  year={2023},
  doi={https://doi.org/10.1101/2022.06.20.496733},
  publisher={Cold Spring Harbor Laboratory}
}
```

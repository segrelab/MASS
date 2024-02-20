# Multi-Attribute Subset Selection enables prediction of representative phenotypes across microbial populations

This repository serves as the codebase for our preprint: https://doi.org/10.1101/2022.06.20.496733.

Multi-Attribute Subset Selection (MASS) is an algorithm based on mixed-integer linear programming to optimally split a set of attributes (or variables) into predictor and response attributes. For this, the response attributes are modeled as linear combinations of predictor attributes.

We used this algorithm to identify the most descriptive growth conditions in microbial phenotyping experiments. This has implications for designing high-throughput phenotyping efforts as it can help to profile only the most informative experimental conditions and therefore reduces the experimental effort. We further observed, that the top predictor conditions selected by MASS comprise meaningful metabolic axis.

The general nature of the algorithm should make it applicable in other areas as well.

## Content
 - DATASET 1: marine
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

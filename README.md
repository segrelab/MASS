# Prediction of representative phenotypes using Multi-Attribute Subset Selection (MASS)

## Content
 - DATASET 1: marine
     - MASS: [MASS-marine.ipynb](./Code/MASS-marine.ipynb) (jupyter notebook)
     - RF: [RF_evalutation_w_baseline-marine.ipynb](./Code/RF_evalutation_w_baseline-marine.ipynb) (jupyter notebook)
 - DATASET 2: fermentation (BacDive1)
     - data download and preparation: [download_bacdive_traits.Rmd](./Code/download_bacdive_traits.Rmd) (Rmarkdown notebook)
     - MASS: _MISSING_ (jupyter notebook)
     - RF: [RF_evalutation_w_baseline-bacdive1.ipynb](./Code/RF_evalutation_w_baseline-bacdive1.ipynb) (jupyter notebook)
 - DATASET 3: yeast species
     - MASS: [MASS-yeast.ipynb](./Code/MASS-yeast.ipynb) (jupyter notebook)
     - RF: [RF_evalutation_w_baseline-yeast.ipynb](./Code/RF_evalutation_w_baseline-yeast.ipynb) (jupyter notebook)
 - manuscript figures: [mass_figures.Rmd](./Code/mass_figures.Rmd)  (Rmarkdown notebook)
 
## TODO 
 - [x] move RF notebooks
 - [x] clean-up paths in RF notebooks
     - [x] DS1
     - [x] DS2
     - [x] DS3
 - [ ] MASS notebooks
     - [x] DS1
     - [ ] DS2
     - [x] DS3
  - [x] clean-up repo (delete?? GIT_milp_multiclass_11_65.py , MIP_classify_RF_f1_score_agg.csv, compare_classify_RF_f1_score_ag, MIP_classify_f1_scores_greedy2_2,MIP_classify_z_all76_greedy2.csv )
 - [ ] update readme (abstract, usage etc.)

## Citation
If you find it useful, please cite our paper as follows:

```
@article{forchielli2022prediction,
  title={Prediction of representative phenotypes using multi-output subset selection},
  author={Forchielli, Elena and Wang, Taiyao and Thommes, Meghan and Paschalidis, Ioannis Ch and Segre, Daniel},
  journal={bioRxiv},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}
```

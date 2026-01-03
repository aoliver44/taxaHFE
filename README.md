![test workflow](https://github.com/aoliver44/taxaHFE/actions/workflows/test.yml/badge.svg)
# **TaxaHFE, TaxaHFE-ML, and DietML** <a><img src='pictures/logo.png' align="right" height="150" /></a>
<a><img src='pictures/dietML_logo_Nov18.png' align="right" height="120" /></a>
```taxaHFE``` is a program for performing hierarchical feature engineering on data with a taxonomic organization (e.g., microbiome data, dietary data). ```taxaHFE-ML``` is a variation of ```taxaHFE``` that splits the input data into training and test data, and then performs hierarchical feature engineering on the training data. It then tests the predictivness of those features on the left out test data. If your goal is to reduce a set of hierarchically organized features, use ```taxaHFE```. If your goal is to use hierarchically organized features in machine learning models, we recommend using ```taxaHFE-ML```.

Finally, ```dietML``` is a machine learning pipeline which takes in a set of features and creates a model to predict a features of interest. ```dietML``` is what is used in ```taxaHFE-ML```, but it can also be used with non-hierarchical features too. 

### [Learn more about TaxaHFE and TaxaHFE-ML](./taxahfe.md)

### [Learn more about DietML](./dietml.md)

### [Instructions for developers](./developers.md)

## **Acknowledgments**

Special thanks to Stephanie M.G. Wilson for the logos.
</br>

## **Contribute**

Feel free to raise an issue, contribute with a pull request, or reach out!

------------------------------
## **Citation**

If you use TaxaHFE or TaxaHFE-ML, please cite:

Oliver A, Kay M, Lemay DG. TaxaHFE: a machine learning approach to collapse microbiome datasets using taxonomic structure. *Bioinformatics Advances* 2023;3:, https://doi.org/10.1093/bioadv/vbad165


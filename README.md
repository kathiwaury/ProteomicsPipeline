# ProteomicsPipeline

The pipeline is based on the [drake](https://github.com/ropensci/drake) package (William Michael Landau, (2018). The drake R package: a pipeline toolkit for reproducibility and high-performance computing. Journal of Open Source Software, 3(21), 550, https://doi.org/10.21105/joss.00550) and can be used for the analysis of proteomics data that is labeled for binary classification.

Following steps are performed:
- Wilcoxon rank sum test for univariate analysis
- Fold change calculation
- Builing of a classification model based on Random Forest
- Ranking and selection of proteins (based on mean Gini index of Random Forest model)

The selected proteins are analysed by:
- Gene ontology term enrichment analysis
- Pathway enrichment analysis
- Disease ontology enrichment analysis

## Data requirements

The pipeline requires two input files in csv-format.
1. The **proteomics data set** contains the protein measurements. The first column must contain a unique sample identifier. The second column has to contain the label used for classification (e.g. "healthy" and "affected"). Only two different labels are allowed and the class of effect (e.g. with treatment/mutation/symptoms) should be listed first if the column is converted into a factor, i.e. the class of effect should be first alphabetically, to allow for the correct baseline for the fold change calculations.
2. The **protein information data set** should contain two columns. The first column will contain a unique description of all proteins/targets measured. The second column will contain the associated gene name of the target. For every target used in the proteomics data set (column names), there should be an entry in the first column of the protein information data set. Duplicated gene names are allowed.

Are any of the requirments not met, it will lead to the termination of the pipeline before any analysis steps.

## Set up

1. Download the files of the repository.
2. Move the two required csv files into the same directory where you saved the downloaded files.
3. Edit plan.R by replacing "path/to/csvfile" with the file name of the two csv files in data_in and protein_info_in respectively.
4. Set the directory containing the repository files as your working directory.
5. Run the following code
```
install.packages("drake")
library(drake)
new_cache()
```
A folder to save your cache has been created. Your directory should have the following structure now:
- .drake 
- R
  - functions.R
  - packages.R 
  - plan.R
- make.R
- Report.Rmd
- your_proteomics_dataset.csv
- your_protein_information_dataset.csv
6. Execute make.R in R to run the pipeline analysis. The results will be in the form of an HTML report that will be created in the same directory.

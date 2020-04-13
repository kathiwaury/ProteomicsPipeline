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

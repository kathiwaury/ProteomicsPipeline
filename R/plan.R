plan <- drake_plan(
  data = data_check(file_in("correctDataGENFI.csv")),
  proteinInfo = protein_info_check(file_in("correctProteinInfoGENFI.csv"), data),
  wilcoxonResults = wilcoxon_test(data, proteinInfo),
  foldChangeResults = fold_change(data, proteinInfo),
  plotResults = volcano_plot(wilcoxonResults, foldChangeResults),
  randomForestResults = random_forest(data),
  proteinsRanked = protein_ranking(randomForestResults, proteinInfo),
  proteinSetUnique = protein_selection(proteinsRanked, wilcoxonResults),
  proteinID = protein_ID(proteinSetUnique),
  GOAnalysisResults = GO_enrichment_analysis(proteinID),
  DOAnalysisResults = DO_enrichment_analysis(proteinID),
  pathwayAnalysisResults = pathway_analysis(proteinID),
  report = rmarkdown::render(
    knitr_in("Report.Rmd"),
    output_file = file_out("Report.html"),
    quiet = TRUE
  )
)

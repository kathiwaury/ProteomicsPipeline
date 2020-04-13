plan <- drake_plan(
  data_in = data_check(file_in("path/to/csvfile")),
  protein_info_in = protein_info_check(file_in("path/to/csvfile"), data_in),
  wilcox_results = wilcox_test(data_in, protein_info_in),
  fold_change_results = fold_change(data_in, protein_info_in),
  plot_results = volcano_plot(wilcox_results, fold_change_results),
  random_forest_results = random_forest(data_in),
  proteins_ranked = protein_ranking(random_forest_results, protein_info_in),
  unique_protein_list = protein_selection(proteins_ranked, wilcox_results),
  protein_ID_list = clusterProfiler::bitr(unique_protein_list$Gene, fromType = "SYMBOL", 
    toType = "ENTREZID", OrgDb = "org.Hs.eg.db"),
  #uniprot_analysis_results = uniprot_analysis(protein_ID_list),
  GO_term_analysis_results = clusterProfiler::enrichGO(gene = protein_ID_list$ENTREZID, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID"),
  pathway_analysis_results = enrichPathway(protein_ID_list$ENTREZID, organism = "human", readable = TRUE),
  disease_analysis_results = DOSE::enrichDGN(protein_ID_list$ENTREZID, readable = TRUE),
  report = rmarkdown::render(
    knitr_in("Report.Rmd"),
    output_file = file_out("Report.html"),
    quiet = TRUE
  )
)


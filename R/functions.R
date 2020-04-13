data_check <- function(file) {
  #check file
  if (file.exists(file) == FALSE) {
    stop("File wasn't loaded. Couldn't find a file of that name.")
  }
  if (reader::get.ext(file) != "csv") {
    stop("File wasn't loaded. The file needs to be in csv format.")
  }
  print("File check successful!")
  #load file as data frame into R
  data <- read.csv(file = file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  print("File read!")
  #check for empty  or NULL data object
  if (is.null(data) == TRUE || nrow(data) == 0) {
    stop("Loaded object is NULL or empty.")
  }
  #check for NA values in data frame
  if (any(is.na(data)) == TRUE) {
    stop("NA values detected.")
  }
  #check if number of columns is sufficient
  if (length(data) <= 3) {
    stop("Less than 4 columns detected.")
  }
  #check correct type of ID and class columns
  N <- length(data)
  if (typeof(data[, 1]) != "character")  {
    stop("The first column must be of type character.")
  }
  if (typeof(data[, 2]) != "character") {
    stop("The second column must be of type character.")
  }
  #check correct type of proteomics data
  for (i in 3:N) {
    if (typeof(data[, i]) != "integer" && typeof(data[, i]) != "double") {
      stop("The columns are not of the required type.")
    }
  }
  data[, 2] <- factor(data[, 2])
  #check for only two classes in the classification column
  if (length(levels(factor(data[, 2]))) > 2) {
    stop("More than two classes for classification detected.")
  }
  if (length(levels(factor(data[, 2]))) < 2) {
    stop("Less than two classes for classification detected.")
  }
  print("Data check successful!")
  
  data
}

protein_info_check <- function(file, table) {
  #check file
  if (file.exists(file) == FALSE) {
    stop("Protein info file wasn't loaded. Couldn't find a file of that name.")
  }
  if (get.ext(file) != "csv") {
    stop("Protein info file wasn't loaded. The file needs to be in csv format.")
  }
  print("Protein info file check successful!")
  #load file as data frame into R
  data <- read.csv(file = file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  #check for more than two columns
  print("Protein info file read!")
  if (length(data) > 2) {
    print("More than two columns detected. Only the first two will be used.")
    data <- data[, 1:2]
  }
  #check for empty  or NULL data object
  if (is.null(data) == TRUE || nrow(data) == 0) {
    stop("Loaded object is NULL or empty.")
  }
  #check for NA values in data frame
  if (any(is.na(data)) == TRUE) {
    stop("NA values detected.")
  }
  #check correct type of label and class columns
  if (typeof(data[, 1]) != "character") {
    stop("The first column is not of type character.")
  }
  if (typeof(data[, 2]) != "character") {
    stop("The second column is not of type character.")
  }
  
  #check that there is a match for every feature in the proteomics data
  features <- names(table)[3:length(table)]
  if (all(features %in% data[, 1]) == FALSE) {
    stop("Protein info is not provided for all features in the data set.")
  }
  
  names(data)[1] <- "Description"
  names(data)[2] <- "Gene"
  
  print("Protein info data check successful!")
  
  data
}

wilcox_test <- function(data, protein_info) {
  N <- length(data)
  #create subsets of data based on two classes
  class.1 <- levels(data[, 2])[1]
  class.2 <- levels(data[, 2])[2]
  subset.1 <- data[which(data[, 2] == class.1), ]
  subset.2 <- data[which(data[, 2] == class.2), ]
  
  #create data frame to save p-values
  wilcoxResults <- data.frame(Feature = names(data)[3:N], 
                              Pvalue = rep(NA, N-2),
                              Pvalue.adj = rep(NA, N-2),
                              stringsAsFactors = FALSE)
  
  #loop over all features
  for (i in 3:N) {
    #check if wilcox test throws error
    try_value <- try(wilcox.test(subset.1[, i], subset.2[, i]))
    if (class(try_value) == "try-error") {
      print("Failed wilcox test.")
    } else {
      #perform wilcox rank-sum test 
      w <- wilcox.test(subset.1[, i], subset.2[, i])
      wilcoxResults$Pvalue[i-2] <- w$p.value
    }
  }
  
  #calculate adjusted p-value for every variable in dataset
  for (i in 1:(N-2)) {
    wilcoxResults$Pvalue.adj <- p.adjust(wilcoxResults$Pvalue, method="fdr")
  }
  
  wilcoxResults <- merge(protein_info, wilcoxResults, by.x = names(protein_info)[1], 
                         by.y = "Feature")
  
  #filter for significant p-values
  wilcoxResultsSignificant <- wilcoxResults %>% 
    filter(Pvalue.adj < 0.05)
  
  if (all(is.na(wilcoxResults$Pvalue))) {
    print("Wilcox test failed for all features.")
  }
  
  #check for number of significant results
  if (nrow(wilcoxResultsSignificant) == 0) {
    stop("No significant features were found. Analysis is terminated.")
  } else {
    print(paste0(nrow(wilcoxResultsSignificant), " significant features were found."))
  }
  
  wilcoxResults
}

fold_change <- function(data, protein_info) {
  N <- length(data)
  #create subsets of data based on two classes
  class.1 <- levels(data[, 2])[1] #must be affected group
  class.2 <- levels(data[, 2])[2] #must be control group
  subset.1 <- data[which(data[, 2] == class.1), ]
  subset.2 <- data[which(data[, 2] == class.2), ]
  
  #create data frame to save fold changes
  foldChangeResults <- data.frame(Feature = names(data)[3:N], 
                                  foldChange = rep(NA, N-2),
                                  stringsAsFactors = FALSE)
  
  #loop over all features
  for (i in 3:N) {
    #check if wilcox test throws error
    try_value <- try(mean(subset.1[, i]) - mean(subset.2[, i]))
    if (class(try_value) == "try-error") {
      print("Failed fold change calculation.")
    } else {
      #calculate fold change
      FC <- mean(subset.1[, i]) - mean(subset.2[, i])
      foldChangeResults$foldChange[i-2] <- FC
    }
  }
  
  #check results data object
  if (nrow(foldChangeResults) == 0) {
    print("Table of fold change results could not be created.")
  }
  
  foldChangeResults <- merge(protein_info, foldChangeResults, by.x = names(protein_info)[1], 
                             by.y = "Feature")
  
  foldChangeResults
}

volcano_plot <- function(wilcox_results, fold_change_results) {
  
  #create data frame of combined wilcox and fold change results
  results <- merge(wilcox_results, fold_change_results[, -2], by.x = names(wilcox_results)[1],
                   by.y = names(fold_change_results)[1])
  
  results <- na.omit(results)
  
  if (nrow(results) < 20) {
    return("Too many NA values in univariate analysis and fold change calculation.")
  }
  
  #save names and indices of interesting proteins
  indSignificant <- which((results$Pvalue.adj < 0.05 & results$foldChange > 0.25) |
                            (results$Pvalue.adj < 0.05 & results$foldChange < -0.25))
  resultsSignificant <- results[indSignificant, ]
  
  plot <- ggplot() +
    geom_point(data = results, mapping = aes(x = foldChange, y = -log10(Pvalue.adj))) +
    geom_point(data = resultsSignificant, mapping =
                 aes(x = foldChange, y = -log10(Pvalue.adj)), color = "dodgerblue4") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted", size = 1) +
    geom_hline(    yintercept = -log10(0.05), linetype = "dotted", size = 1) +
    labs(x = "Log(2) fold change", y = "-Log10(adjusted p-value)",
         title = "Volcano plot of P-value and fold change")
  
  plot
}

random_forest <- function(data) {
  #save length of dataset
  N <- length(data)
  #create subsets of data based on two classes
  class.1 <- levels(factor(data[, 2]))[1]
  class.2 <- levels(factor(data[, 2]))[2]
  
  #set number of iterations and number of trees per random forest
  M <- 10
  trees <- 6000
  
  #initiate empty list
  randomForestResults <- list()
  
  for (i in 1:M){
    #create training (70%) and test (30%) data split 
    set.seed(i)
    trainInds <- data[, 2] %>%
      caret::createDataPartition(p = 0.7, list = FALSE)
    trainData <- data[trainInds, ]
    
    #bootstrap dataset for balanced classes by sampling with replacement 100 times
    set.seed(1)
    indsClass.1 <- sample(which(trainData[, 2] == class.1), 100, replace = TRUE)
    set.seed(1)
    indsClass.2 <- sample(which(trainData[, 2] == class.2), 100, replace = TRUE)
    #combine bootstraped dataset
    trainData <- rbind(trainData[indsClass.1, ], trainData[indsClass.2, ])
    testData <- data[-trainInds, ]
    trainDataVariables <- as.matrix(trainData[, -c(1:2)])
    trainDataClasses <- as.factor(trainData[, 2])
    testDataVariables <- as.matrix(testData[, -c(1:2)])
    testDataClasses <- as.factor(testData[, 2])
    
    #perform random forest
    fit.RF <- randomForest(x = trainDataVariables, y = trainDataClasses, 
                           xtest = testDataVariables, ytest = testDataClasses, ntree = trees)
    randomForestResults[[length(randomForestResults)+1]] <- fit.RF
  }
  
  if (length(randomForestResults) == 0){
    stop("Random Forest could not be performed.")
  } 
  
  randomForestResults
}

protein_ranking <- function(random_forest_results, protein_info) {
  
  # initiate empty vector to save Gini Index results
  featureImportance <- vector()
  
  for (i in 1:length(random_forest_results)) {
    #subset one random forest object
    fit.RF <- random_forest_results[[i]]
    #extract feature importance ranking and append to results vector
    featureImportance <- cbind(featureImportance, importance(fit.RF))
  }
  #calculate mean MeanDecreaseGini for every feature
  meanFeatureImportance <- data.frame("Importance" = apply(featureImportance, 1, mean)) %>%
    tibble::rownames_to_column("Feature") 
  #add gene name to data frame 
  meanFeatureImportance <- merge(protein_info, meanFeatureImportance, by.x = names(protein_info)[1], 
                                 by.y = "Feature") 
  #sort by mean Gini Index
  meanFeatureImportance <- arrange(meanFeatureImportance, desc(Importance)) #%>%
  #tibble::rowid_to_column("Rank")
  
  meanFeatureImportance
}

protein_selection <- function(proteinsRanked, wilcoxResults) {
  
  wilcoxResultsSignificant <- wilcoxResults %>% 
    #filter significant p-values
    filter(Pvalue.adj < 0.05) %>% arrange(Pvalue.adj) %>%
    #remove duplicate genes
    distinct(Gene, .keep_all = TRUE) %>%
    #keep first two columns
    dplyr::select(names(wilcoxResults)[1], names(wilcoxResults)[2])
  
  #calculate highest importance
  maxImportance <- proteinsRanked$Importance[1]
  #calculate minimum importnace realtive to maximum importance
  minImportance <- maxImportance * 0.33
  proteinsRankedSignificant <- proteinsRanked %>% 
    #filter proteins with high enough importance
    filter(Importance > minImportance) %>%
    #remove duplicate genes
    distinct(Gene, .keep_all = TRUE) %>%
    #keep first two columns
    dplyr::select(names(proteinsRanked)[1], names(proteinsRanked)[2])
  
  N <- nrow(proteinsRankedSignificant)
  #if number of proteins is too low, expand with Wilcoxon results
  if (N < 30) {
    if (nrow(wilcoxResultsSignificant) >= (30 - N)) {
      #filter wilcoxon resuls for genes not yet in ranked protein list
      uniqueInds <- which(!(wilcoxResultsSignificant$Gene %in% proteinsRankedSignificant$Gene))
      wilcoxUnique <- wilcoxResultsSignificant[uniqueInds, ]
      #add unique wilcox proteins to table
      proteinsRankedSignificant <- rbind(proteinsRankedSignificant, wilcoxUnique[1:(30-N), ])
    } else {
      stop("Less than 30 significant proteins found. Analysis is terminated.")
    }
  }
  
  #check for only unique genes
  if (any(duplicated(proteinsRankedSignificant$Gene) == TRUE)) {
    proteinsRankedSignificant <- distinct(proteinsRankedSignificant, Gene, .keep_all = TRUE)
  }
  #add rank to data frame
  proteinsRankedSignificant <- tibble::rowid_to_column(proteinsRankedSignificant, "Rank")
  
  proteinsRankedSignificant
}

uniprot_analysis <- function(protein_ID_list) {
  protein_ID_list
}

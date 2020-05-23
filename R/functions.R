#funtion 1
data_check <- function(file) {
  #check file exists
  if (file.exists(file) == FALSE) {
    stop("File wasn't loaded. Couldn't find a file of that name.")
  }
  
  #check file type
  if (reader::get.ext(file) != "csv") {
    stop("File wasn't loaded. The file needs to be in csv format.")
  }
  
  #load file as data frame into R
  data <- read.csv(file = file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  #check for empty  or NULL data object
  if (is.null(data) == TRUE || nrow(data) == 0) {
    stop("Loaded object is NULL or empty.")
  }
  
  #check for NA values in data frame
  if (any(is.na(data)) == TRUE) {
    stop("NA values detected.")
  }
  
  N <- length(data)
  #check if number of columns is sufficient
  if (N <= 3) {
    stop("Less than 4 columns detected.")
  }
  #rename first two columns
  names(data)[1] <- "SampleID"
  names(data)[2] <- "Class"
  
  #sample ID is unique
  if (any(duplicated(data$SampleID) == TRUE)) {
    stop("Sample description is not unique.")
  }
  
  #check correct type of proteomics data
  for (i in 3:N) {
    if (typeof(data[, i]) != "integer" && typeof(data[, i]) != "double") {
      stop("The proteomics columns are not of the numeric type.")
    }
  }
  
  #convert class column into factor
  data$Class <- factor(data$Class)
  
  #set case class
  inputCounter <- 0
  case <- readline(prompt = "Enter the name of the case group: ")
  #check if class name exists
  while (!(case %in% levels(data$Class))) {
    case <- readline(prompt = "Could not find that label. Try again to enter the name of the case group: ")
    #track number of failed attempts
    inputCounter <- inputCounter + 1
    if (inputCounter > 5) {
      stop("Repeatedly failed attempts to assign class. Pipeline is terminated.")
    }
  }
  
  #set control class
  inputCounter <- 0
  control <- readline(prompt = "Enter the name of the control group: ")
  #check if class name exists and is different from case class
  while (!(control %in% levels(data$Class)) || case == control) {
    if (case == control) {
      control <- readline(prompt = "Cannot use the same label for control and case class. Try again to enter the name of the control group: ")
    } else {
      control <- readline(prompt = "Could not find that label. Try again to enter the name of the control group: ")      
    }
    #track number of failed attempts
    inputCounter <- inputCounter + 1
    if (inputCounter > 5) {
      stop("Repeatedly  failed attempts to assign class. Pipeline is terminated.")
    }
  }
  
  #check for only two classes in the classification column
  if (length(levels(data$Class)) > 2) {
    print("More than two labels in class column detected. Samples that are not part of case and control class will be removed.")
    #remove all samples that are not part of control or case class
    data <- data[which(data$Class %in% c(case, control)), ]
  }
  
  #set order of levels
  data$Class <- factor(data$Class, levels = c(case, control))
  
  print("Data check successful!")
  
  data
}

#function 2
protein_info_check <- function(file, table) {
  #check file
  if (file.exists(file) == FALSE) {
    stop("Protein info file wasn't loaded. Couldn't find a file of that name.")
  }
  if (get.ext(file) != "csv") {
    stop("Protein info file wasn't loaded. The file needs to be in csv format.")
  }
  
  #load file as data frame into R
  data <- read.csv(file = file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  #check for more than two columns
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
  
  #check correct type of description column
  if (typeof(data[, 1]) != "character") {
    stop("The first column is not of type character.")
  }
  
  #check correct type of gene name column
  if (typeof(data[, 2]) != "character") {
    stop("The second column is not of type character.")
  }
  
  #check that there is a match for every feature in the proteomics data
  features <- names(table)[3:length(table)]
  
  if (all(features %in% data[, 1]) == FALSE) {
    stop("Protein info is not provided for all features in the data set.")
  }
  
  #rename first two columns
  names(data)[1] <- "Description"
  names(data)[2] <- "Gene"
  
  print("Protein info data check successful!")
  
  data
}

#function 3
wilcoxon_test <- function(data, proteinInfo) {
  
  N <- length(data)
  #create subsets of data based on two classes
  case <- levels(data[, 2])[1]
  control <- levels(data[, 2])[2]
  subsetCase <- data[which(data[, 2] == case), ]
  subsetControl <- data[which(data[, 2] == control), ]
  
  #create data frame to save p-values
  wilcoxonResults <- data.frame(Feature = names(data)[3:N], 
                                Pvalue = rep(NA, N-2),
                                Pvalue.adj = rep(NA, N-2),
                                stringsAsFactors = FALSE)
  
  failCounter <- 0
  #loop over all features
  for (i in 3:N) {
    #check if wilcox test throws error
    try_value <- try(wilcox.test(subsetCase[, i], subsetControl[, i]), silent=TRUE)
    if (class(try_value) == "try-error") {
      failCounter <- failCounter + 1
    } else {
      #perform wilcox rank-sum test 
      w <- wilcox.test(subsetCase[, i], subsetControl[, i])
      #save p-value in data frame
      wilcoxonResults$Pvalue[i-2] <- w$p.value
    }
  }
  
  if (all(is.na(wilcoxonResults$Pvalue))) {
    warning("Wilcoxon test failed for all features.")
  } else {
    
    if (failCounter > 0) {
      if (failCounter == 1) {
        print("1 Wilcoxon test failed.")
      } else {
        print(paste0(failCounter, " Wilcoxon tests failed."))
      }
      
    }
    
    #calculate adjusted p-value for every variable in dataset
    for (i in 1:(N-2)) {
      wilcoxonResults$Pvalue.adj <- p.adjust(wilcoxonResults$Pvalue, method="fdr")
    }
    
    #filter for significant p-values
    wilcoxonResultsSignificant <- wilcoxonResults %>% 
      filter(Pvalue.adj < 0.05)
    
    #check for number of significant results
    if (nrow(wilcoxonResultsSignificant) == 0) {
      stop("No significant features were found. Analysis is terminated.")
    } else {
      print(paste0(nrow(wilcoxonResultsSignificant), " significant features were found."))
    }
  }
  
  #add gene name column
  wilcoxonResults <- merge(proteinInfo, wilcoxonResults, by.x = "Description", 
                           by.y = "Feature")
  
  
  wilcoxonResults
}

#function 4
fold_change <- function(data, proteinInfo) {
  
  N <- length(data)
  #create subsets of data based on two classes
  case <- levels(data[, 2])[1]
  control <- levels(data[, 2])[2] 
  subsetCase <- data[which(data[, 2] == case), ]
  subsetControl <- data[which(data[, 2] == control), ]
  
  #create data frame to save fold changes
  foldChangeResults <- data.frame(Feature = names(data)[3:N], 
                                  foldChange = rep(NA, N-2),
                                  stringsAsFactors = FALSE)
  
  failCounter <- 0
  #loop over all features
  for (i in 3:N) {
    #check if calculation throws error
    try_value <- try(mean(subsetCase[, i]) - mean(subsetControl[, i]))
    if (class(try_value) == "try-error" || is.na(try_value)) {
      failCounter <- failCounter + 1
    } else {
      #calculate fold change
      FC <- mean(subsetCase[, i]) - mean(subsetControl[, i])
      foldChangeResults$foldChange[i-2] <- FC
    }
  }
  
  if (all(is.na(foldChangeResults$foldChange))) {
    warning("Fold change calculation failed for all features.")
  } else {
    if (failCounter > 0) {
      if (failCounter == 1) {
        print("1 fold change could not be calculated.")
      } else {
        print(paste0(failCounter, " fold changes could not be calculated."))
      }
      
    }
  }
  
  foldChangeResults <- merge(proteinInfo, foldChangeResults, by.x = "Description", 
                             by.y = "Feature")
  
  foldChangeResults
}

#function 5
volcano_plot <- function(wilcoxonResults, foldChangeResults) {
  #create data frame of combined wilcoxon and fold change results
  results <- merge(wilcoxonResults, foldChangeResults[, -2], by.x = "Description",
                   by.y = "Description")
  
  #remove features with NAs
  results <- na.omit(results)
  #check if enough features are available
  if (nrow(results) < 20) {
    return("Too many NA values in wilcoxon results and fold change calculation.")
  }
  
  #save names and indices of interesting proteins
  indSignificant <- which((results$Pvalue.adj < 0.05 & results$foldChange > 0.25) |
                            (results$Pvalue.adj < 0.05 & results$foldChange < -0.25))
  resultsSignificant <- results[indSignificant, ]
  
  #create plot
  plot <- ggplot() +
    #all data points
    geom_point(data = results, mapping = aes(x = foldChange, y = -log10(Pvalue.adj))) +
    #significant data points
    geom_point(data = resultsSignificant, mapping =
                 aes(x = foldChange, y = -log10(Pvalue.adj), name = Gene), color = "dodgerblue4") +
    #add threshold lines
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted", size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", size = 1) +
    labs(x = "Log(2) fold change", y = "-Log10(adjusted p-value)",
         title = "Volcano plot of P-value and fold change")
  
  #make plot interactive
  interactivePlot <- ggplotly(plot)
  
  interactivePlot
}

#function 6
random_forest <- function(data) {
  #save length of dataset
  N <- length(data)
  #create subsets of data based on two classes
  case <- levels(factor(data[, 2]))[1]
  control <- levels(factor(data[, 2]))[2]
  
  #set number of iterations and number of trees per random forest
  M <- 20
  Ntrees <- 6000
  
  #initiate empty list
  randomForestResults <- list()
  
  for (i in 1:M){
    #create train (70%) and test (30%) data split 
    set.seed(i)
    trainInds <- data[, 2] %>%
      caret::createDataPartition(p = 0.7, list = FALSE)
    trainData <- data[trainInds, ]
    testData <- data[-trainInds, ]
    
    #bootstrap dataset for balanced classes by sampling both classes with replacement 100 times
    set.seed(1)
    indsCase <- sample(which(trainData[, 2] == case), 100, replace = TRUE)
    set.seed(1)
    indsControl <- sample(which(trainData[, 2] == control), 100, replace = TRUE)
    #combine bootstraped dataset
    trainData <- rbind(trainData[indsCase, ], trainData[indsControl, ])
    
    #create predictor and outcome sets for train and test data
    trainDataVariables <- as.matrix(trainData[, -c(1:2)])
    trainDataClasses <- as.factor(trainData[, 2])
    testDataVariables <- as.matrix(testData[, -c(1:2)])
    testDataClasses <- as.factor(testData[, 2])
    
    #perform random forest
    fitRF <- randomForest(x = trainDataVariables, y = trainDataClasses, 
                          xtest = testDataVariables, ytest = testDataClasses, ntree = Ntrees)
    
    #save fitted model into results list
    randomForestResults[[length(randomForestResults)+1]] <- fitRF
  }
  
  #check for empty list
  if (length(randomForestResults) == 0){
    stop("Random Forest could not be performed.")
  } 
  
  randomForestResults
}

#function 7
protein_ranking <- function(randomForestResults, proteinInfo) {
  
  # initiate empty vector to save Gini Index results
  featureImportance <- vector()
  
  for (i in 1:length(randomForestResults)) {
    #subset one random forest object
    fitRF <- randomForestResults[[i]]
    #extract feature importance ranking and append to results vector
    featureImportance <- cbind(featureImportance, importance(fitRF))
  }
  
  #calculate mean importance of all models for every feature
  meanFeatureImportance <- data.frame("Importance" = apply(featureImportance, 1, mean)) %>%
    tibble::rownames_to_column("Feature")
  
  #add gene name to data frame
  meanFeatureImportance <- merge(proteinInfo, meanFeatureImportance, by.x = "Description",
                                 by.y = "Feature")
  
  #sort by importance
  meanFeatureImportance <- arrange(meanFeatureImportance, desc(Importance))
  
  meanFeatureImportance
}

#function 8
protein_selection <- function(proteinsRanked, wilcoxonResults) {
  
  wilcoxonResultsSignificant <- wilcoxonResults %>% 
    #remove failed wilcoxon tests
    na.omit(wilcoxonResults) %>%
    #filter significant p-values
    filter(Pvalue.adj < 0.05) %>% 
    #order by p-value
    arrange(Pvalue.adj) %>%
    #remove duplicate genes
    distinct(Gene, .keep_all = TRUE) %>%
    #keep first two columns
    dplyr::select("Description", "Gene")
  
  #calculate highest importance
  maxImportance <- proteinsRanked$Importance[1]
  #calculate importance threshold relative to maximum importance
  minImportance <- maxImportance * 0.33
  
  proteinsRankedSignificant <- proteinsRanked %>% 
    #filter proteins above importance threshold
    filter(Importance > minImportance) %>%
    #remove duplicate genes
    distinct(Gene, .keep_all = TRUE) %>%
    #keep first two columns
    dplyr::select("Description", "Gene")
  
  N <- nrow(proteinsRankedSignificant)
  
  #add Source column and set source as random forest
  proteinsRankedSignificant$Source <- rep("Random Forest", N)
  
  #if number of proteins is too low, expand with wilcoxon results
  if (N < 50) {
    if (nrow(wilcoxonResultsSignificant) >= (50 - N)) {
      #filter wilcoxon resuls for genes not yet in ranked protein list
      uniqueInds <- which(!(wilcoxonResultsSignificant$Gene %in% proteinsRankedSignificant$Gene))
      wilcoxonUnique <- wilcoxonResultsSignificant[uniqueInds, ]
      #set source as wilcoxon test
      wilcoxonUnique$Source <- rep("Wilcoxon Test", nrow(wilcoxonUnique))
      #add unique wilcox proteins to table
      proteinsRankedSignificant <- rbind(proteinsRankedSignificant, wilcoxonUnique[1:(50-N), ])
    } else {
      stop("Less than 50 significant proteins found. Analysis is terminated.")
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

#function 9
protein_ID <- function(proteinSetUnique) {
  #extract Entrez ID for every gene symbol in protein list
  entrezList <- clusterProfiler::bitr(proteinSetUnique$Gene, fromType = "SYMBOL", toType = 
                                        c("ENTREZID"), OrgDb = "org.Hs.eg.db")
  mergedList <- merge(proteinSetUnique, entrezList, by.x = "Gene", by.y = "SYMBOL")
  mergedList <- mergedList[, c("Rank", "Gene", "Description", "Source", "ENTREZID")] %>% 
    arrange(Rank)
  
  mergedList
}

#function 10
GO_enrichment_analysis <- function(proteinID) {
  GOAnalysisResults <- clusterProfiler::enrichGO(gene = proteinID$ENTREZID,
    OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = TRUE)
  
  GOAnalysisResults
}

#function 11
DO_enrichment_analysis <- function(proteinID) {
  DOAnalysisResults <- DOSE::enrichDO(proteinID$ENTREZID, readable = TRUE)
  
  DOAnalysisResults
}

#function 12
pathway_analysis <- function(proteinID) {
  pathwayAnalysisResults <- ReactomePA::enrichPathway(proteinID$ENTREZID, 
    organism = "human", readable = TRUE)
  
  pathwayAnalysisResults
}


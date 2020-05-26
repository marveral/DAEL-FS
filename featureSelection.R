library(limma)
library(httr)
library(rlist)
library(stringr)
library(KnowSeq)


#' @param data The data parameter is an expression matrix or data.frame that contains the genes in the columns and the samples in the rows.
#' @param labels A vector or factor that contains the labels for each samples in data parameter.
#' @param vars_selected The genes selected to use in the feature selection process. It can be the final DEGs extracted with the function \code{\link{limmaDEGsExtraction}} or a custom vector of genes.
#' @param mode The algorithm used to calculate the genes ranking. The possibilities are two: mrmr, rf and da.
#' @param disease The name of a disease in order to calculate the Disease Association ranking by using the DEGs indicated in the vars_selected parameter.
#' @return A vector that contains the ranking of genes.
daFeatureSelection <-function(data,labels,mode="da",disease="",maxGenes=30,returnEvidences=TRUE){

  if(!is.data.frame(data) && !is.matrix(data)){
    
    stop("The data argument must be a dataframe or a matrix.")
    
  }
  if(dim(data)[1] != length(labels)){
    
    stop("The length of the rows of the argument data must be the same than the length of the labels. Please, ensures that the rows are the samples and the columns are the variables.")
    
  }
  
  if(!is.character(labels)  && !is.factor(labels)){stop("The class of the labels parameter must be character vector or factor.")}
  if(is.character(labels)){ labels <- as.factor(labels) }
  
  vars_selected <- colnames(data)
  
  if(disease == ""){
    stop("Please, indicate a disease name to acquire the Disease Association Score and Feature selection.")
  }
  if(mode == "da"){
    cat("Calculating ranking of biological relevant genes by using DA implementation...\n")
  }else if(mode == 'daRed'){
    cat("Calculating ranking of biological relevant genes by using DA-Red implementation...\n")
  }
  
  data <- as.data.frame(apply(data,2,as.double))
    

  overallRanking <- rep(0,length(vars_selected))
  names(overallRanking) <- vars_selected

  disease_ <- str_replace_all(disease,' ','-')
  
  r_Ensembl <- GET(paste("https://api.opentargets.io/v3/platform/public/search?q=",disease_,"&size=1&filter=disease",sep = ""))
  respon <- httr::content(r_Ensembl)
  
  if ( 'size' %in% names(respon) && respon$size == 0) stop('Disease not found')
  
  # Look for related genes
  disease.id <- respon$data[[1]]$id
  url  <- paste("https://api.opentargets.io/v3/platform/public/association/filter?disease=",disease.id,"&size=10000",sep='')
  response <- GET(url)
  response <- httr::content(response)
  found.symbols <- unlist(list.map(response$data,target$gene_info$symbol))
  scores <- unlist(list.map(response$data,association_score$overall))
  
  # Keep related genes and their scores
  found.index <- which(found.symbols %in% vars_selected)
  overallRanking[found.symbols[found.index]] <- scores[found.index]
  names(overallRanking) <- vars_selected
  
  # Sort by scores
  overallRanking <- sort(overallRanking,decreasing = TRUE)
  cat("Disease Association ranking: ")
  cat(names(overallRanking)[1:maxGenes])
  cat("\n")

  if (mode=='da'){
    overallRanking <- overallRanking[1:maxGenes]
    if (returnEvidences) return(list('ranking'=overallRanking,evidences=DEGsEvidences(names(overallRanking), disease, size=1000)))
    else return(overallRanking)
  }
  
  # Look for genes evidences
  evidences <- DEGsEvidences(names(overallRanking), disease, size=1000)

  # If mode == daLOD, calculate LOD scores
  if(grepl('daLOD',mode)){
    cat("Calculating genes LOD scores...\n")
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    limma.table.aux <- limmaDEGsExtraction(t(data),labels)$Table

    if(length(levels(labels)) == 2) limma.table <- data.frame('lods'=range01(limma.table.aux$B))
    else if(length(levels(labels)) > 2) limma.table <- data.frame('lods'=range01(rowMeans(limma.table.aux$lods)))
    rownames(limma.table) <- rownames(limma.table.aux)
  }
  
  cat("Calculating genes scores...\n")
  
  # Output: list of selected genes
  selected.genes <- list()
  
  genes <- names(overallRanking)
  
  # Select the gene with maximun relevance
  if(grepl('daRed',mode)) selected.genes[[names(overallRanking)[1]]] = overallRanking[1]
  else{
    tmp.ranking <- overallRanking[rownames(limma.table)]
    tmp <- as.numeric(tmp.ranking)  * limma.table$lods
    tmp.max  <- which.max(tmp)[1]
    selected.genes[[names(tmp.ranking)[tmp.max]]] = tmp[tmp.max]
  }

  # Create empty redundance matrix
  redundances <- matrix(-1,ncol=maxGenes,nrow=length(overallRanking))
  rownames(redundances) = names(overallRanking)
  colnames(redundances) <- rep('',maxGenes)
  colnames(redundances)[1] <- names(selected.genes[1])

  # Iter over all genes
  for (i in seq(maxGenes-1)){
    # actual max score
    max.score <- -1000
    
    # act.genes contains genes to select (genes that are not already selected)
    act.genes <- genes[ ! genes %in% names(selected.genes) ]
    
    # Iter in genes to select
    for ( gen in act.genes){
      # Gene relevance is it's score
      rel <- as.numeric(overallRanking[[gen]])
      # Redundance begin as 0
      red <- 0

      if (rel < max.score) break

      # Calculate redundance between actual gene and selected genes
      for ( sel in names(selected.genes)){
        # Calculate and save redudance 
        if (redundances[gen,sel] == -1){
          redundances[gen,sel] = 0
          if ( class(evidences[[gen]]) == 'list' && class(evidences[[sel]]) == 'list'){
            gen.nevs <- 0
            for (type in names(evidences[[gen]])){
              if (type %in% names(evidences[[sel]])){
                # type.total contains found coincidences for actual type of evidences
                type.total  <- 0
                ncol <- length(evidences[[gen]][[type]][[1]]$evidence)
                # Iter on gen evidences
                for (row1 in evidences[[gen]][[type]]){
                  gen.nevs <- gen.nevs + 1

                  # Boolean matrix that contains coincidences between gen and sel
                  act.total.matrix <- matrix(0,nrow=length(evidences[[sel]][[type]]),ncol=ncol)
                  # Iter on gen evidences
                  for (row2 in evidences[[sel]][[type]]){
                    # Add row to act.total.matrix with boolean values
                    act.total.matrix <- rbind(act.total.matrix, unlist(row1$evidence) == unlist(row2$evidence))
                    # This row fully coincide with actual evidencie, so we stop searching
                    if (rowSums(tail(act.total.matrix,1) == ncol))  break
                  }
                  # If there are any row that fully coincide add 1 (this evidences is fully contained in sel data)
                  if ( any(rowSums(act.total.matrix) == ncol)) type.total = type.total + 1
                  else{
                    # Check which row is the most similar to the actual evidence
                    # Add the percentage of coincidence ( num. coincidences / ncol)
                    act.total.field <- colSums(act.total.matrix)
                    type.total = type.total + length(which(colSums(act.total.matrix) >= 1))/ncol
                  }
                }
                # Add score of type evidences
                redundances[gen,sel] = redundances[gen,sel] + type.total
              }
            }
            # Normalize score dividing by the number of evidences
            if (gen.nevs > 0 )redundances[gen,sel] = redundances[gen,sel] / gen.nevs
          }
        }
        red <- red + redundances[gen,sel]
      }

      # Score = relevance - relevance * redundance/ num of selected genes
      if (grepl('daRed',mode)) score <- rel - red/length(selected.genes) * rel
      else if(grepl('daLOD',mode)){
        if  (gen %in% rownames(limma.table)) score <- rel * limma.table[gen,'lods'] - red/length(selected.genes) * rel
        else {
          score <-  score <- 0 - red/length(selected.genes) * rel
        }
      }else{
        stop('ERROR IN MODE')
      }
      # If actual score is max keep actual gene
      if (score > max.score){
        max.score <- score
        act.selected <- gen
        act.red <- red
      }
    }
    # Save gene with maximun score
    selected.genes[[act.selected]] <- max.score
    colnames(redundances)[length(selected.genes)] <- act.selected
  }
  if (returnEvidences) return(list('ranking'=selected.genes,'evidences'=evidences[match(names(selected.genes),names(evidences))]))
  else return(selected.genes)
}

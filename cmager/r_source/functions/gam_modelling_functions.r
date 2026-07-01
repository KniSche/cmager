######################################################
#### old gene mapping function
#require(HGNChelper)
#map_genes_to_model_features = function(model_features, input_dataset){
#        gene_map_df = data.frame(
#            "og" = model_features
#            )
#        
#        og_updated = checkGeneSymbols(gene_map_df$og)
#        gene_map_df$og_updated = og_updated$Suggested.Symbol
#        gene_map_df[is.na(gene_map_df$og_updated),"og_updated"] = gene_map_df[is.na(gene_map_df$og_updated),"og"]
#        rownames(gene_map_df) = gene_map_df$og
#        
#        input_genes = colnames(input_dataset)
#        input_genes_updated = checkGeneSymbols(input_genes)

#    
#        # trim the values, removing NAs and updating NAs to original input(s)
#        input_genes_updated[is.na(input_genes_updated[,"Suggested.Symbol"]),"Suggested.Symbol"] = input_genes_updated[is.na(input_genes_updated[,"Suggested.Symbol"]),"x"]
#        input_genes_updated = input_genes_updated[input_genes_updated[,"Suggested.Symbol"] %in% gene_map_df[,"og_updated"],]

#    
#        # match up the previous names with the updated og training dataset labels
#        gene_map_df$input_matched_og = NA
#        gene_map_df$input_matched_updated = NA
#        
#        # add in the input dataset gene names original
#        gene_map_df$input_matched_og[match(input_genes_updated[,"Suggested.Symbol"], gene_map_df[,"og_updated"])] = input_genes_updated[,"x"]
#        # add in the input dataset gene names updated
#        gene_map_df$input_matched_updated[match(input_genes_updated[,"Suggested.Symbol"], gene_map_df[,"og_updated"])] = input_genes_updated[,"Suggested.Symbol"]



#        # subset the input dataset by the overlapping genes
#        input_dataset = input_dataset[,colnames(input_dataset) %in% gene_map_df[,"input_matched_og"]]
#        
#        # add old names back in where there is ambiguity
#        gene_map_df[grep("///", gene_map_df$og_updated),"og_updated"] = gene_map_df[grep("///", gene_map_df$og_updated),"og"]
#       
#        # update the new input data matrix colnames
#        colnames(input_dataset) = gene_map_df[,"og"][match(colnames(input_dataset), gene_map_df[,"input_matched_og"])]

#        return(input_dataset)
#    }



process_data_for_gam2 <- function(
    predict_matrix, 
    model_features,
    modality,
    null_data,
    is_already_cp10k_normalised = FALSE,
    library_sizes,
    log_transform = TRUE,
    k
    ){
          if(missing(modality) | sum(modality %in% c("cells", "nuclei")) != length(modality)){
            return(cat('ERROR: unrecognised modality, use either "cells" or "nuclei"'))
        } else {
            modality_vector = as.numeric(c("nuclei"= 0, "cells" = 1)[modality])
        }
        if(length(modality_vector) != nrow(predict_matrix)){
            if(length(modality_vector) == 1){
                modality_vector = rep(modality_vector, nrow(predict_matrix))
            } else {
                return(cat('ERROR: modality vector is the wrong length'))
            }
        }

        # predict_matrix to matrix
            predict_matrix = as.matrix(predict_matrix)
        # extract obs and vars names
            barcodes = rownames(predict_matrix)
            genes    = colnames(predict_matrix)
            
        # feature organisation ... mapping external dataset features to the model features here
            features_in = model_features[model_features %in% genes]
            features_out = model_features[!model_features %in% genes]
            
            cat("    ", k, ": matched ", length(features_in), " out of ", length(model_features), " gene names. Filling in values with training data mean to ", length(features_out), " genes (", paste(features_out, collapse=", "), ") ...", "\n", sep="")

        # Create the output matrix
         output_matrix = matrix(
            0,
            nrow(predict_matrix),
            length(model_features)
            )
        rownames(output_matrix) = rownames(predict_matrix)
        colnames(output_matrix) = model_features

         # fill in new matrix with the input_data values
            output_matrix[barcodes,features_in] = predict_matrix[barcodes,features_in]
        # library size normalise
            if(is_already_cp10k_normalised == FALSE){
                if(sum(is.na(library_sizes)) > 1 | sum(barcodes %in% names(library_sizes)) != length(barcodes)){
                        return(cat("Please check library sizes argument, ensure that the vector corresponds to the barcode rows of the input_data matrix"))
                    } else {
                        output_matrix = (output_matrix / library_sizes[barcodes])*1e4
                }
            }

        # log transform?
            if(log_transform != FALSE){output_matrix = log(output_matrix+1)}

        # add the train mean data, filling missing features:
             if(length(features_out) > 0){
            # add in mean values for genes (based on training data)
                for(gene in features_out){output_matrix[,gene] = null_data[gene, "train_mean"]}
            }
        # Multiplying out interactions by modality
             # now add interactions columns, and the modality columns
            interactions_matrix = output_matrix*modality_vector
                colnames(interactions_matrix) = paste0(colnames(output_matrix), "_interaction")
                rownames(interactions_matrix) = barcodes
        
        # construct the final output matrix
            output_matrix = cbind(output_matrix, interactions_matrix, modality_vector) ### check here
        
            colnames(output_matrix) = c(model_features, colnames(interactions_matrix), "modality")
        return(output_matrix)
    }


predict_age_fromCT_results =  function(
        input_data,
        experiment_predictions,
        model_list,
        granularity,
        modality = "nuclei",
        good_barcodes= NULL,
        is_already_cp10k_normalised = FALSE,
        library_sizes = NULL,
        log_transform = TRUE,
        add_to_experimentRDS = FALSE,
        verbose=TRUE
    ){
        
# Step 1
    # pull out essential variables
        if(is.null(library_sizes)){
            library_sizes = rowSums(input_data)
        }
        names(library_sizes) = rownames(input_data)

        if(length(modality) > 1){
             names(modality) = rownames(input_data)
        }
        
# Step 2 - not required in with the python updated gene mapping implementation

    # re-map genes back to model hgnc symbols, subsetting new dataset to the model space where possible.
        
        all_model_genes <- unique(unlist(lapply(names(model_list), function(celltype){
          if(length(model_list[[celltype]]) > 1){
              return(model_list[[celltype]][["features"]])
          }
        })))
        
        
#   
#     query_genes = colnames(input_data)
#        features_in = query_genes[query_genes %in% all_model_genes]
#        features_out = query_genes[query_genes %in% all_model_genes]
#                if(verbose == TRUE){
#                    cat("    ", "all_model_genes", ": matched ", length(features_in), " out of ", length(all_model_genes), " gene names." )
#                    cat("    ", "\n", sep="")
#                }
#    
#        input_data = map_genes_to_model_features(model_features = all_model_genes, input_dataset = input_data)



# Step 3
    # filter out bad barcodes (optional)
    input_data_barcodes = rownames(input_data)
        if(!is.null(good_barcodes)){
            input_data = input_data[input_data_barcodes[input_data_barcodes %in% good_barcodes],,drop=F]
            library_sizes = library_sizes[input_data_barcodes[input_data_barcodes %in% good_barcodes]]
            input_data_barcodes = rownames(input_data)
            }

        
# Step 4
    if(add_to_experimentRDS != TRUE){
        cell_type_age_predictions = data.frame(
            "barcode" = character(),
            "predicted_age" = numeric()  
        )
    }
    # Pull up cell type models and check for cell type labels in the dataset
     for(k in names(model_list)){
        k_barcodes = rownames(experiment_predictions[[granularity]][experiment_predictions[[granularity]]$predictions == k,]) 
        k_barcodes = k_barcodes[k_barcodes %in% input_data_barcodes]
        
        if(length(modality) > 1){
            modality_k = modality[k_barcodes]
        } else {
            modality_k = modality
        }
         
        #k = "AtrialCardiomyocytes"
      # load model
        gamsel.model = model_list[[k]]
        
      # CHECK FOR VIABLE MODEL      
        #if(gamsel.model == "no significant features"){gamsel.model = "no model"}
        if(class(gamsel.model) %in% c("NULL", "character")){
         if(verbose == TRUE){
            cat("    ", k, ": No model found... Check model generation steps. Likely due to sparse data in training data or lack of informative genes. Skipping age prediction...", "\n", sep="")  
            cat("    ", "\n", sep="")
            cat("    ", "\n", sep="")
        }
            if(add_to_experimentRDS != TRUE){
                cell_type_age_predictions = rbind(
                            cell_type_age_predictions,
                            data.frame(
                                "barcode" = k_barcodes,
                                "predicted_age" = rep(NA, length(k_barcodes)
                            )
                        )
                    )
                }
        } else {

            
        # pick up model parameters:    
            lambda.index = gamsel.model$cv_gamsel$index.1se
            features = gamsel.model$features
            optimal_features = gamsel.model$optimal_features
            
            if(granularity == "fine_grain"){
                null_data = as.data.frame(rowSums((cell_adjust[all_model_genes,c("cell.mean", "nuclei.mean")])))
                  #underlying distributions are non-log - for the sake of the below back-transformation, we will log this mean as an approximation
                    null_data[,1] = log1p(null_data[,1])
                    
                colnames(null_data) = "train_mean"
            } else {   
                null_data = gamsel.model[["gene_averages"]]
            }
        # identify the barcodes where we have the k "celltype" label predicted:
            #predict_subset = experiment_predictions[[granularity]][experiment_predictions[[granularity]]$predictions == k,]
            predict_matrix = input_data[input_data_barcodes[input_data_barcodes %in% k_barcodes],,drop=F]

        # Catch the zero cells found situation
           if(nrow(predict_matrix) > 0){
                   if(verbose == TRUE){
                     cat("    ", k, ": found ", nrow(predict_matrix), " predicted cells...", "\n", sep="")
                    }
                # identify feature space overlap
                  #   features_in = features[features %in% colnames(predict_matrix)]
                   #  features_out = features[!features %in% colnames(predict_matrix)]
                    predict.matrix = process_data_for_gam2(
                        predict_matrix = predict_matrix, 
                        model_features = features,
                        modality = modality_k,
                        null_data = null_data,
                        is_already_cp10k_normalised = is_already_cp10k_normalised,
                        library_sizes = library_sizes,
                        log_transform = log_transform,
                        k=k
                    )
                    if(granularity == "fine_grain"){
                          predict.matrix = expm1(predict.matrix[,features,drop=FALSE])
                              #the fine_grain model was trained on non-log-transformed data... v 0.5.0
                              #to be addressed in later versions
                              
                    }
             # Predict from model:
                   if(verbose == TRUE){
                    cat("    ", k, ": Predicting cell age.............", "\n", sep="")
                   }
                        prediction_k = as.vector(predict(gamsel.model$cv_gamsel$gamsel.fit, newdata = predict.matrix, index=lambda.index))
                        names(prediction_k) = rownames(predict.matrix)
             # 
                 if(add_to_experimentRDS == TRUE){
                     if(verbose == TRUE){
                         cat("    ", k, ":.................... adding to experiment_predictions dataframe...", "\n", sep="")
                     }
                        if(!"predicted_age" %in% colnames(experiment_predictions[[granularity]])){
                          experiment_predictions[[granularity]]$predicted_age = 0
                        }
                          experiment_predictions[[granularity]][names(prediction_k),"predicted_age"] = prediction_k
                     
                        } else {
                            cell_type_age_predictions = rbind(
                                cell_type_age_predictions,
                                data.frame(
                                    "barcode" = names(prediction_k),
                                    "predicted_age" = prediction_k
                                )
                            )
                        }
               } else {
               if(verbose == TRUE){
                    cat("    ", k, ": No cells found in dataset. Skipping age prediction as there is nothing to predict.", "\n", sep="")  
                   }
               }
                if(verbose == TRUE){
                    cat("    ", k, ": DONE.", "\n", sep="")
                    cat("    ", "\n", sep="")
                }
            }
         flush.console()
        }
        if(add_to_experimentRDS == FALSE){
            return(cell_type_age_predictions)
        } else {
            return(experiment_predictions)
        }
    }



############################


  
      

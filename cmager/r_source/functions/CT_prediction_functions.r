#!/usr/bin/env Rscript

# read the CellTypist predicted files
  read.ct.res = function(filepath){
    x = read.csv(filepath, head=T)
    x[,1] = gsub("\t.*$", "", x[,1])
       rownames(x) = x[,1]
    return(x)
  }
  




# create and format a named prediction vector
  label.ct.predictions = function(so, predictions.df){
      labels = predictions.df[colnames(so),"predictions"]
        names(labels) = colnames(so)
      return(labels)
    }



# process predictions into dataframe
    prediction.df = function(
        prediction_dataframe_list,
        probability_matrix_list,
        thresholds,
        cluster_grouping
        ){
          reference_factors = list(
            "coarse_grain" = levels(meta.data[,"coarse_grain"]),
            "mid_grain" = levels(meta.data[,"mid_grain"]),
            "fine_grain" = levels(meta.data[,"fine_grain"])
          )
        
        
        granularity = "coarse_grain"
        coarse_predictions.df = data.frame(
            "barcode" = prediction_dataframe_list[[granularity]][,1], 
            # no reference as it's test data: "reference" = meta.data[CT.test.predictions_coarse[,1],"coarse_grain"],
            "predictions" = factor(prediction_dataframe_list[[granularity]][,2], levels=reference_factors[[granularity]]),
            "probability" = rep(0, dim(prediction_dataframe_list[[granularity]])[1])
          )
            for(k in reference_factors[[granularity]]){
              coarse_predictions.df[coarse_predictions.df$predictions == k,"probability"] = 
              probability_matrix_list[[granularity]][coarse_predictions.df$predictions == k,k]
            }
          rownames(coarse_predictions.df) = coarse_predictions.df[,1]
    
      # odds ratio between the max prob and next-best prob
      coarse_predictions.df$odds_ratio = sapply(
        1:dim(probability_matrix_list[[granularity]])[1], 
        function(i){
          p1 = probability_matrix_list[[granularity]][i,colnames(probability_matrix_list[[granularity]][,-1])[order(-as.numeric(probability_matrix_list[[granularity]][i,-1]))][1]]
          p2 = probability_matrix_list[[granularity]][i,colnames(probability_matrix_list[[granularity]][,-1])[order(-as.numeric(probability_matrix_list[[granularity]][i,-1]))][2]]
          odds_ratio.i = (p1 * (1-p2)) / (p2 * (1 - p1))
          return(odds_ratio.i)
        }
      )
        
#    
    
    granularity = "mid_grain"  
     mid_predictions.df = data.frame(
            "barcode" = prediction_dataframe_list[[granularity]][,1], 
            # no reference as it's test data: "reference" = meta.data[CT.test.predictions_coarse[,1],"coarse_grain"],
            "predictions" = factor(prediction_dataframe_list[[granularity]][,2], levels=reference_factors[[granularity]]),
            "probability" = rep(0, dim(prediction_dataframe_list[[granularity]])[1])
          )
#            for(k in reference_factors[[granularity]]){
#              mid_predictions.df[mid_predictions.df$predictions == k,"probability"] = 
#              probability_matrix_list[[2]][mid_predictions.df$predictions == k,k]
#            }
          rownames(mid_predictions.df) = mid_predictions.df[,1]
     # heirarchichally label the cells 
     #mid_hPredict = mclapply(
     mid_hPredict = lapply(
          names(cluster_grouping),
          function(k){
          
             temp.df = probability_matrix_list[[granularity]][coarse_predictions.df$predictions == k,colnames(probability_matrix_list[[granularity]]) %in% cluster_grouping[[k]],drop=F]
             
             if(nrow(temp.df) > 0){
               k.predictions = sapply(
                1:nrow(temp.df),
                function(l){
                  return(paste(names(which.max(temp.df[l,,drop=F])),"###",rownames(temp.df)[l],sep=""))
                  #return(names(which.max(temp.df[l,])))
                }
              )
            } else {
            k.predictions = character()
            }
             return(k.predictions)
           }#,
           #mc.cores=7
         )
      
   # unlist the coarse -> mid top probability predictions
   # replace all predictions with them
     mid_predictions.df[gsub("^.*###", "", unlist(mid_hPredict)), "predictions"] = gsub("###.*$", "", unlist(mid_hPredict))
   
   # Retrieve the probability scores from the CT matrix
      for(k in reference_factors[[granularity]]){
        mid_predictions.df[mid_predictions.df$predictions == k,"probability"] = 
        probability_matrix_list[[granularity]][mid_predictions.df$predictions == k,k]
      }
   
  # odds ratio between the max prob and next-best prob
    mid_predictions.df$odds_ratio = sapply(
      1:dim(probability_matrix_list[[granularity]])[1], 
      function(i){
        p1 = probability_matrix_list[[granularity]][i,colnames(probability_matrix_list[[granularity]][,-1])[order(-as.numeric(probability_matrix_list[[granularity]][i,-1]))][1]]
        p2 = probability_matrix_list[[granularity]][i,colnames(probability_matrix_list[[granularity]][,-1])[order(-as.numeric(probability_matrix_list[[granularity]][i,-1]))][2]]
        odds_ratio.i = (p1 * (1-p2)) / (p2 * (1 - p1))
        return(odds_ratio.i)
      }
    )
      
    
    
    granularity = "fine_grain"  
     fine_predictions.df = data.frame(
            "barcode" = prediction_dataframe_list[[granularity]][,1], 
            # no reference as it's test data: "reference" = meta.data[CT.test.predictions_coarse[,1],"coarse_grain"],
            "predictions" = factor(prediction_dataframe_list[[granularity]][,2], levels=reference_factors[[granularity]]),
            "probability" = rep(0, dim(prediction_dataframe_list[[granularity]])[1])
          )
#            for(k in reference_factors[[granularity]]){
#              fine_predictions.df[fine_predictions.df$predictions == k,"probability"] = 
#              probability_matrix_list[[3]][fine_predictions.df$predictions == k,k]
#            }
          rownames(fine_predictions.df) = fine_predictions.df[,1]
          
     # heirarchichally label the cells 
    # fine_hPredict = mclapply(
      fine_hPredict = lapply(
          names(cluster_grouping),
          function(k){
          
             temp.df = probability_matrix_list[[granularity]][coarse_predictions.df$predictions == k,colnames(probability_matrix_list[[granularity]]) %in% cluster_grouping[[k]],drop=F]
             if(nrow(temp.df) > 0){
               k.predictions = sapply(
                1:nrow(temp.df),
                function(l){
                  return(paste(names(which.max(temp.df[l,,drop=F])),"###",rownames(temp.df)[l],sep=""))
                  #return(names(which.max(temp.df[l,])))
                }
              )
            } else {
              k.predictions = character()
            }
             return(k.predictions)
           }#,
           #mc.cores=7
         )
      
   # unlist the coarse -> fine top probability predictions
   # replace all predictions with them
     fine_predictions.df[gsub("^.*###", "", unlist(fine_hPredict)), "predictions"] = gsub("###.*$", "", unlist(fine_hPredict))
   
     weird.cells = grep("###NA\\.", unlist(fine_hPredict))
   
     
   # Retrieve the probability scores from the CT matrix
      for(k in reference_factors[[granularity]]){
        if(k %in% fine_predictions.df$predictions){
          fine_predictions.df[fine_predictions.df$predictions == k,"probability"] = 
          probability_matrix_list[[granularity]][fine_predictions.df$predictions == k,k]
          #order(-probability_matrix_list[[granularity]][fine_predictions.df$predictions == k,-1])[1]
        }
      }
    
  # odds ratio between the max prob and next-best prob
    fine_predictions.df$odds_ratio = sapply(
      1:dim(probability_matrix_list[[granularity]])[1], 
      function(i){
        p1 = probability_matrix_list[[granularity]][i,colnames(probability_matrix_list[[granularity]][,-1])[order(-as.numeric(probability_matrix_list[[granularity]][i,-1]))][1]]
        p2 = probability_matrix_list[[granularity]][i,colnames(probability_matrix_list[[granularity]][,-1])[order(-as.numeric(probability_matrix_list[[granularity]][i,-1]))][2]]
        odds_ratio.i = (p1 * (1-p2)) / (p2 * (1 - p1))
        return(odds_ratio.i)
      }
    )
    
    
    return(
      # "hi"
      list(
        "coarse_grain" = coarse_predictions.df,
        "mid_grain" =  mid_predictions.df,
        "fine_grain" =  fine_predictions.df
      )
    )
  }
 
 
# Probability thresholding function
filter_low_probability = function(
    prediction_dataframe_list,
    thresholds,
    odds_ratio = F,
    apply_hierarchy = T
  ){
       
    reference_factors = list(
        "coarse_grain" = levels(meta.data[,"coarse_grain"]),
        "mid_grain" = levels(meta.data[,"mid_grain"]),
        "fine_grain" = levels(meta.data[,"fine_grain"])
      )
    if(odds_ratio == F){  
      if(class(thresholds) == "numeric"){
      
           for(granularity in names(reference_factors)){
              prediction_dataframe = prediction_dataframe_list[[granularity]]
              prediction_dataframe$predictions2 = factor(prediction_dataframe$predictions, levels=c(levels(prediction_dataframe$predictions), "unknown"))
              for(k in reference_factors[[granularity]]){
                  prediction_dataframe[which(prediction_dataframe$predictions == k)[
                    prediction_dataframe[prediction_dataframe$predictions == k,"probability"] <= thresholds
                    ],"predictions2"] = "unknown"
               }
              prediction_dataframe_list[[granularity]] = prediction_dataframe
            }
            
      } else {

        for(granularity in names(reference_factors)){
          prediction_dataframe = prediction_dataframe_list[[granularity]]
          prediction_dataframe$predictions2 = factor(prediction_dataframe$predictions, levels=c(levels(prediction_dataframe$predictions), "unknown"))
          for(k in reference_factors[[granularity]]){
              prediction_dataframe[which(prediction_dataframe$predictions == k)[
                prediction_dataframe[prediction_dataframe$predictions == k,"probability"] <= thresholds[[granularity]][k,2]
                ],"predictions2"] = "unknown"
           }
          prediction_dataframe_list[[granularity]] = prediction_dataframe
        }
      }
    } else {
         if(class(thresholds) == "numeric"){
      
           for(granularity in names(reference_factors)){
              prediction_dataframe = prediction_dataframe_list[[granularity]]
              prediction_dataframe$predictions2 = factor(prediction_dataframe$predictions, levels=c(levels(prediction_dataframe$predictions), "unknown"))
              for(k in reference_factors[[granularity]]){
                  prediction_dataframe[which(prediction_dataframe$predictions == k)[
                    prediction_dataframe[prediction_dataframe$predictions == k,"odds_ratio"] <= thresholds
                    ],"predictions2"] = "unknown"
               }
              prediction_dataframe_list[[granularity]] = prediction_dataframe
            }
            
      } else {

        for(granularity in names(reference_factors)){
          prediction_dataframe = prediction_dataframe_list[[granularity]]
          prediction_dataframe$predictions2 = factor(prediction_dataframe$predictions, levels=c(levels(prediction_dataframe$predictions), "unknown"))
          for(k in reference_factors[[granularity]]){
              prediction_dataframe[which(prediction_dataframe$predictions == k)[
                prediction_dataframe[prediction_dataframe$predictions == k,"odds_ratio"] <= thresholds[[granularity]][k,2]
                ],"predictions2"] = "unknown"
           }
          prediction_dataframe_list[[granularity]] = prediction_dataframe
        }
      }
    }
    
 
  
    if(apply_hierarchy == T){

      hierarchy_list = list(
        "Cardiomyocytes" = list(
           "AtrialCardiomyocytes" = list(
              "AtrialCardiomyocytesLeft",
              "AtrialCardiomyocytesRight",
              "AtrialCardiomyocytesCycling"
            ),
            "VentricularCardiomyocytes" = list(
              "VentricularCardiomyocytesLeftCompact",
              "VentricularCardiomyocytesRightCompact",
              "VentricularCardiomyocytesLeftTrabeculated",
              "VentricularCardiomyocytesRightTrabeculated",
              "VentricularCardiomyocytesCycling"
            ),
            "CardiacConductionSystem" = list(
              "SinoatrialNodePacemakerCells",
              "AtrioventricularNodePacemakerCells",
              "VentricularConductionSystemProximal",
              "VentricularConductionSystemDistal"
            )
        ),
        "Mesenchymal" = list(
          "Fibroblasts" = list(
            "GreatVesselAdventitialFibroblasts",  
            "CoronaryVesselAdventitialFibroblasts",
            "MyocardialInterstitialFibroblasts",
            "SubEpicardialFibroblasts",   
            "Myofibroblasts",                            
            "LymphNodeFibroblasticReticularCells",       
            "ValveInterstitialCells"
          ),                     
          "MuralCells" = list(
            "GreatVesselSmoothMuscleCells",              
            "CoronarySmoothMuscleCells",                 
            "DuctusArteriosusSmoothMuscleCells",         
            "CoronaryPericytes" 
          ),
          "PericardialCells" = list(
            "PericardialCellsIntermediate",              
            "PericardialCellsFibrous",                   
            "PericardialCellsParietal" 
          )
        ),
        "Endothelium" = list(
          "BloodVesselEndothelialCells" = list(
            "GreatVesselArterialEndothelialCells",       
            "GreatVesselVenousEndothelialCells",         
            "CoronaryArterialEndothelialCells",          
            "CoronaryVenousEndothelialCells",            
            "CoronaryCapillaryEndothelialCells" 
          ),
          "EndocardialCells" = list(
            "EndocardialCells",                  
            "EndocardialCushionCells",                   
            "ValveEndothelialCells"
          ),
          "LymphaticEndothelialCells" = list(
            "LymphaticEndothelialCells"
          )
        ),
        "Epicardium" = list(
          "EpicardialCells" = list(
            "MesothelialEpicardialCells",                
            "EpicardiumDerivedCells"
          )
        ),
        "Neural" = list(
          "Neurons" = list(
            "NeuronPrecursors",                          
            "ChromaffinCells",                           
            "SympatheticNeurons",                        
            "ParasympatheticNeurons"
          ),
          "Glia" = list(
            "SchwannCellPrecursors",                     
            "SchwannCells"
          )
        ),
        "Leukocytes" = list(
          "MyeloidCells" = list(
            "MonocytesMPOpos",                           
            "Monocytes",                                 
            "MonocyteDerivedCells",                      
            "MacrophagesCX3CR1pos",                      
            "MacrophagesTIMD4pos",                       
            "MacrophagesLYVE1pos",                       
            "MacrophagesATF3pos",                        
            "DendriticCellsType1",                       
            "DendriticCellsMature",                      
            "PlasmacytoidDendriticCells",                
            "MastCells",                                 
            "Megakaryocytes"
          ),
          "LymphoidCells" = list(
            "TCellsCD4pos",                              
            "TCellsCD8pos",                              
            "TregsCD4pos",                               
            "ProBCells",                                 
            "BCells",                                    
            "BCellsMS4A1pos",                            
            "NaturalKillerCells",                        
            "InnateLymphoidCells"
          )
        )
      )

     # extract unique hierarchy labels
       cell_type_levels = paste(gsub("\\.|[0-9]$", "_", names(unlist(hierarchy_list))),unlist(hierarchy_list),sep="")
     
     # label cell_type across hierarchy
       h_predictions = paste(
        prediction_dataframe_list[["coarse_grain"]]$predictions2,
        prediction_dataframe_list[["mid_grain"]]$predictions2,
        prediction_dataframe_list[["fine_grain"]]$predictions2,
        sep="_"
       )
       
       
      # re-label unknown according to hierarchy - if higher is unknown... then lower is also unknown.
        h_predictions[grep("^unknown_|^unknown_unknown_", h_predictions)] = "unknown_unknown_unknown"
      
    
        prediction_dataframe_list[["h_predictions"]] = h_predictions
        
    }
    return(prediction_dataframe_list)
  } 
 
 
    hierarchy_list = list(
        "Cardiomyocytes" = list(
           "AtrialCardiomyocytes" = list(
              "AtrialCardiomyocytesLeft",
              "AtrialCardiomyocytesRight",
              "AtrialCardiomyocytesCycling"
            ),
            "VentricularCardiomyocytes" = list(
              "VentricularCardiomyocytesLeftCompact",
              "VentricularCardiomyocytesRightCompact",
              "VentricularCardiomyocytesLeftTrabeculated",
              "VentricularCardiomyocytesRightTrabeculated",
              "VentricularCardiomyocytesCycling"
            ),
            "CardiacConductionSystem" = list(
              "SinoatrialNodePacemakerCells",
              "AtrioventricularNodePacemakerCells",
              "VentricularConductionSystemProximal",
              "VentricularConductionSystemDistal"
            )
        ),
        "Mesenchymal" = list(
          "Fibroblasts" = list(
            "GreatVesselAdventitialFibroblasts",  
            "CoronaryVesselAdventitialFibroblasts",
            "MyocardialInterstitialFibroblasts",
            "SubEpicardialFibroblasts",   
            "Myofibroblasts",                            
            "LymphNodeFibroblasticReticularCells",       
            "ValveInterstitialCells"
          ),                     
          "MuralCells" = list(
            "GreatVesselSmoothMuscleCells",              
            "CoronarySmoothMuscleCells",                 
            "DuctusArteriosusSmoothMuscleCells",         
            "CoronaryPericytes" 
          ),
          "PericardialCells" = list(
            "PericardialCellsIntermediate",              
            "PericardialCellsFibrous",                   
            "PericardialCellsParietal" 
          )
        ),
        "Endothelium" = list(
          "BloodVesselEndothelialCells" = list(
            "GreatVesselArterialEndothelialCells",       
            "GreatVesselVenousEndothelialCells",         
            "CoronaryArterialEndothelialCells",          
            "CoronaryVenousEndothelialCells",            
            "CoronaryCapillaryEndothelialCells" 
          ),
          "EndocardialCells" = list(
            "EndocardialCells",                  
            "EndocardialCushionCells",                   
            "ValveEndothelialCells"
          ),
          "LymphaticEndothelialCells" = list(
            "LymphaticEndothelialCells"
          )
        ),
        "Epicardium" = list(
          "EpicardialCells" = list(
            "MesothelialEpicardialCells",                
            "EpicardiumDerivedCells"
          )
        ),
        "Neural" = list(
          "Neurons" = list(
            "NeuronPrecursors",                          
            "ChromaffinCells",                           
            "SympatheticNeurons",                        
            "ParasympatheticNeurons"
          ),
          "Glia" = list(
            "SchwannCellPrecursors",                     
            "SchwannCells"
          )
        ),
        "Leukocytes" = list(
          "MyeloidCells" = list(
            "MonocytesMPOpos",                           
            "Monocytes",                                 
            "MonocyteDerivedCells",                      
            "MacrophagesCX3CR1pos",                      
            "MacrophagesTIMD4pos",                       
            "MacrophagesLYVE1pos",                       
            "MacrophagesATF3pos",                        
            "DendriticCellsType1",                       
            "DendriticCellsMature",                      
            "PlasmacytoidDendriticCells",                
            "MastCells",                                 
            "Megakaryocytes"
          ),
          "LymphoidCells" = list(
            "TCellsCD4pos",                              
            "TCellsCD8pos",                              
            "TregsCD4pos",                               
            "ProBCells",                                 
            "BCells",                                    
            "BCellsMS4A1pos",                            
            "NaturalKillerCells",                        
            "InnateLymphoidCells"
          )
        )
      )
      
      
      

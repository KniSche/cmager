#!/usr/bin/env Rscript


#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
######################################################
# ### Read in metadata
# Meta data
  meta.data = readRDS(file.path(script_dir, "resources", "metadata_grains.rds"))
 
  # Classification from the coarse down
  # This is the classification hierarchy - focusing on applying labels based on 
  # The coarse_grain annotation
      cluster_grouping = lapply(
        levels(meta.data[,"coarse_grain"]),
        function(k){
          c(
            as.character(unique(meta.data[meta.data[,"coarse_grain"] == k,"coarse_grain"])),
            as.character(unique(meta.data[meta.data[,"coarse_grain"] == k,"mid_grain"])),
            as.character(unique(meta.data[meta.data[,"coarse_grain"] == k,"fine_grain"]))
          )
        }
      )
      names(cluster_grouping) = levels(meta.data[,"coarse_grain"])
 

# Load in colour palettes
  coarse.col <- read.csv(file.path(script_dir, "resources", "coarse_colors.txt"), sep="\t", head=F)
  coarse.col = rbind(coarse.col, c("unknown", "black"))
  rownames(coarse.col) = coarse.col[,1]
  
  mid.col <- read.csv(file.path(script_dir, "resources", "mid_colors.txt"), sep="\t", head=F)
  mid.col = rbind(mid.col, c("unknown", "black"))
  rownames(mid.col) = mid.col[,1]
  
 
  fine.col.split <- read.csv(file.path(script_dir, "resources", "fine_colors_split.txt"), sep="\t", head=F)
  fine.col.split = rbind(fine.col.split, c("unknown", "black"))
  rownames(fine.col.split) = fine.col.split[,1]
  
######################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  



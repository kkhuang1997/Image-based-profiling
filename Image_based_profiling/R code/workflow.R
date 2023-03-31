source("./feature_extraction.R")
library("progress")
library("stringr")

## feature extraction workflow
root = "E:/suzhou_data/CRC DATA/extracted_features/images_jpg/"
feature_per_folder = vector(mode = "list", length = length(list.files(root)))
pb <- progress_bar$new(total = length(list.files(root)))
for(i in seq_along(list.files(root))){
  pb$tick()
  file = str_split(list.files(root)[[i]], fixed("_"))
  img = Load_img(file.path(root, list.files(root)[[i]]))
  
  z_stack = 1
  max_proj = img
  imgSeed = Seg_func(max_proj)
  # display(imgSeed, all = TRUE)
  # display(max_proj, all = TRUE)
  features = Fea_extrct_jpg(seg = imgSeed , ref = max_proj, 
                        org_line = file[[1]][1],
                        image_field = list.files(root)[[i]],
                        capture_date = file[[1]][2]
  )
  
  ## combine features by stack
  features = do.call(what = rbind, args = features)
  feature_names = colnames(features)
  feature_per_folder[[i]] = features
}
## save features
saveRDS(feature_per_folder, "../rds data/features_for_1960_jpg.rds")

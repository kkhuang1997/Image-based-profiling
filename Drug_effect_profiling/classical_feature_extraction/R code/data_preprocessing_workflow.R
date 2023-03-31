source("./utils.R")
library("progress")
library("MASS")
library("LMGene")
library("onlinePCA")
library("RColorBrewer")
library(paletteer)
library(ggbiplot)
library(dplyr)

## feature extraction workflow
root = "../train_B_red/"
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
saveRDS(feature_per_folder, "../RDS data/features_red.rds")

## metadata
org1_rep1_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org1_rep1_20220924\\field_1\\BF")
org1_rep1_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org1_rep1_20220924\\field_2\\BF")
org1_rep1_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org1_rep1_20220924\\field_3\\BF")
org1_rep1_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org1_rep1_20220924\\field_4\\BF")

org1_rep2_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org1_rep2_20220924\\field_1\\BF")
org1_rep2_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org1_rep2_20220924\\field_2\\BF")
org1_rep2_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org1_rep2_20220924\\field_3\\BF")
org1_rep2_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org1_rep2_20220924\\field_4\\BF")

org2_rep1_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org2_rep1_20220927\\field_1\\BF")
org2_rep1_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org2_rep1_20220927\\field_2\\BF")
org2_rep1_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org2_rep1_20220927\\field_3\\BF")
org2_rep1_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org2_rep1_20220927\\field_4\\BF")

org2_rep2_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org2_rep2_20220927\\field_1\\BF")
org2_rep2_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org2_rep2_20220927\\field_2\\BF")
org2_rep2_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org2_rep2_20220927\\field_3\\BF")
org2_rep2_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org2_rep2_20220927\\field_4\\BF")


org3_rep1_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org3_rep1_20220929\\field_1\\BF")
org3_rep1_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org3_rep1_20220929\\field_2\\BF")
org3_rep1_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org3_rep1_20220929\\field_3\\BF")
org3_rep1_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org3_rep1_20220929\\field_4\\BF")


org3_rep2_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org3_rep2_20220929\\field_1\\BF")
org3_rep2_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org3_rep2_20220929\\field_2\\BF")
org3_rep2_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org3_rep2_20220929\\field_3\\BF")
org3_rep2_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org3_rep2_20220929\\field_4\\BF")

org4_rep1_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org4_rep1_20220930\\field_1\\BF")
org4_rep1_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org4_rep1_20220930\\field_2\\BF")
org4_rep1_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org4_rep1_20220930\\field_3\\BF")
org4_rep1_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org4_rep1_20220930\\field_4\\BF")

org4_rep2_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org4_rep2_20220930\\field_1\\BF")
org4_rep2_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org4_rep2_20220930\\field_2\\BF")
org4_rep2_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org4_rep2_20220930\\field_3\\BF")
org4_rep2_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org4_rep2_20220930\\field_4\\BF")

org5_rep1_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org5_rep1_20221212\\field_1\\BF")
org5_rep1_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org5_rep1_20221212\\field_2\\BF")
org5_rep1_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org5_rep1_20221212\\field_3\\BF")
org5_rep1_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org5_rep1_20221212\\field_4\\BF")

org5_rep2_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org5_rep2_20221212\\field_1\\BF")
org5_rep2_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org5_rep2_20221212\\field_2\\BF")
org5_rep2_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org5_rep2_20221212\\field_3\\BF")
org5_rep2_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org5_rep2_20221212\\field_4\\BF")

org6_rep1_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org6_rep1_20221214\\field_1\\BF")
org6_rep1_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org6_rep1_20221214\\field_2\\BF")
org6_rep1_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org6_rep1_20221214\\field_3\\BF")
org6_rep1_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org6_rep1_20221214\\field_4\\BF")

org6_rep2_fov1 = Meta_extraction("E:\\image_based_profiling_tif\\org6_rep2_20221214\\field_1\\BF")
org6_rep2_fov2 = Meta_extraction("E:\\image_based_profiling_tif\\org6_rep2_20221214\\field_2\\BF")
org6_rep2_fov3 = Meta_extraction("E:\\image_based_profiling_tif\\org6_rep2_20221214\\field_3\\BF")
org6_rep2_fov4 = Meta_extraction("E:\\image_based_profiling_tif\\org6_rep2_20221214\\field_4\\BF")

meta_file_list = list(org1_rep1_fov1,
                      org1_rep1_fov2,
                      org1_rep1_fov3,
                      org1_rep1_fov4,
                      org1_rep2_fov1,
                      org1_rep2_fov2,
                      org1_rep2_fov3,
                      org1_rep2_fov4,
                      org2_rep1_fov1,
                      org2_rep1_fov2,
                      org2_rep1_fov3,
                      org2_rep1_fov4,
                      org2_rep2_fov1,
                      org2_rep2_fov2,
                      org2_rep2_fov3,
                      org2_rep2_fov4,
                      org3_rep1_fov1,
                      org3_rep1_fov2,
                      org3_rep1_fov3,
                      org3_rep1_fov4,
                      org3_rep2_fov1,
                      org3_rep2_fov2,
                      org3_rep2_fov3,
                      org3_rep2_fov4,
                      org4_rep1_fov1,
                      org4_rep1_fov2,
                      org4_rep1_fov3,
                      org4_rep1_fov4,
                      org4_rep2_fov1,
                      org4_rep2_fov2,
                      org4_rep2_fov3,
                      org4_rep2_fov4,
                      org5_rep1_fov1,
                      org5_rep1_fov2,
                      org5_rep1_fov3,
                      org5_rep1_fov4,
                      org5_rep2_fov1,
                      org5_rep2_fov2,
                      org5_rep2_fov3,
                      org5_rep2_fov4,
                      org6_rep1_fov1,
                      org6_rep1_fov2,
                      org6_rep1_fov3,
                      org6_rep1_fov4,
                      org6_rep2_fov1,
                      org6_rep2_fov2,
                      org6_rep2_fov3,
                      org6_rep2_fov4)

library(data.table)
meta_file = rbindlist(meta_file_list, use.names = TRUE)
write.csv(meta_file, "../table/meta_data.csv")


## combine original features and metadata
orig_features_red = readRDS("../RDS data/features_red.rds")
library(readxl)
meta_drugs <- read_excel("../table/meta_drugs.xlsx")
orig_features_red = orig_features_red %>%
  do.call(what = rbind, args = .) %>%
  as_tibble(orig_features_red) %>% 
  mutate(., "Org_size" = .["x.0.s.area"]) 
orig_features_meta = orig_features_red %>%
  left_join(meta_drugs, by = "Field") %>%
  subset(., select = - Line.x) %>%
  rename(., Line = Line.y, Size = Org_size)

## convert to feather format
library(feather)
orig_features_meta = as.matrix(orig_features_meta) ## size.x.0.s.area to size otherwise cannot convert to feather format
orig_features_meta = as.data.frame(orig_features_meta)
write_feather(orig_features_meta, "../RDS data/features_red_meta_drugs.feather")


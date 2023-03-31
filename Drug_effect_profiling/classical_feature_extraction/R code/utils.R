library("EBImage")
library("magrittr")
library("tibble")
library("ggplot2")
library("genefilter")
library("GGally")
library("MaxContrastProjection")
library("Rfast")
library("stringr")

## Load z-stack
Load_img <- function(filepath){
  img = readImage(filepath)
  return(img)
}

## divide slices into stack
Divide <- function(x, n )
{
  i <- ceiling(length(x)/n)
  split(x,x%/%i+1)
}


## segmentation
Seg_func <- function(input){
  threshold = otsu(input)
  img_th = combine( mapply(function(frame, th) frame > th, getFrames(max_proj), threshold, SIMPLIFY=FALSE) )
  ## morphology manipulation
  imgOpened = EBImage::opening(img_th, kern = makeBrush(11, shape = "disc"))
  imgClosed = EBImage::closing(imgOpened, kern = makeBrush(11,shape = "disc"))
  imgSeed = bwlabel(imgClosed)
  return(imgSeed)
}


# Calculate features for every image stack
Fea_extrct <- function(seg, ref, org_line, image_field, capture_date){
  num_stack = length(z_stack)
  features = vector(mode = "list", length = num_stack)
  for(fld in seq_len(num_stack)){
    img = Image(data = ref[,,fld])
    seg_map = Image(data = seg[,,fld])
    mask = matrix(
      as.integer(seg[,,fld] > 0), 
      nrow = nrow(seg), 
      ncol = ncol(seg)) 
    
    mask_erode = erode(mask, kern = makeBrush(11, "disc"))
    dmap = distmap(mask_erode, "euclidean")
    wshed = watershed(dmap, tolerance = 15)
    mask_labeled = propagate(x = mask, seeds = wshed, mask = mask)
    
    if(sum(mask_labeled != 0) == 0) {
      features[[fld]] = matrix(0, nrow = 0, ncol = 1572)
      # features_clumps[[fld]] = matrix(0, nrow = 0, ncol = 1572)
      # features_noseg[[fld]] = matrix(0, nrow = 0, ncol = 118)
      next
    }
    
    ## feature extraction
    features[[fld]] = computeFeatures(
      x = mask_labeled, ref = normalize(img), 
      haralick.scales = c(1, 2, 4, 8, 16, 24, 32, 48, 64, 86, 96, 128, 256, 512)
    )
    features[[fld]] = cbind(features[[fld]], "Line" = org_line , 
                            "Field" = image_field,
                            "Date" = capture_date
                            )
  }
  return(features)
}


# Calculate features for every jpg image
Fea_extrct_jpg <- function(seg, ref, org_line, image_field, capture_date){
  num_stack = length(z_stack)
  features = vector(mode = "list", length = num_stack)
  for(fld in seq_len(num_stack)){
    img = Image(data = ref)
    seg_map = Image(data = seg)
    mask = matrix(
      as.integer(seg > 0), 
      nrow = nrow(seg), 
      ncol = ncol(seg))  
    
    mask_erode = erode(mask, kern = makeBrush(5, "disc"))
    dmap = distmap(mask_erode, "euclidean")
    wshed = watershed(dmap, tolerance = 15)
    mask_labeled = propagate(x = mask, seeds = wshed, mask = mask)
    
    if(sum(mask_labeled != 0) == 0) {
      features[[fld]] = matrix(0, nrow = 0, ncol = 404)
      # features_clumps[[fld]] = matrix(0, nrow = 0, ncol = 1572)
      # features_noseg[[fld]] = matrix(0, nrow = 0, ncol = 118)
      next
    }
    
    ## feature extraction
    features[[fld]] = computeFeatures(
      x = mask_labeled, ref = normalize(img), 
      haralick.scales = c(1, 2, 4, 8, 16, 24, 32, 48, 64, 86, 96, 128, 256, 512)
    )
    features[[fld]] = cbind(features[[fld]], "Line" = org_line , 
                            "Field" = image_field,
                            "Date" = capture_date
    )
  }
  return(features)
}


## metadata extraction
Meta_extraction <- function(root_dir){
  # root_dir = "E:\\image_based_profiling_tif\\org1_rep1_20220924\\field_1\\BF"
  
  meta_file = str_split(root_dir,fixed("\\"))[[1]]
  Line = meta_file[3] %>%
    str_sub(., start = 1, end = 4)
  Replicate = meta_file[3] %>%
    str_sub(., start = 9, end = 9)
  Fov = meta_file[4] %>%
    str_sub(., start = 7, end = 7)
  
  meta_df = data.frame(Org_field = "", 
                       Org_line = "", 
                       Org_replicate = "", 
                       Org_fov = ""
                       )
  
  for( n in seq_along(list.files(root_dir))) {
    image_name = list.files(root_dir)[[n]] %>%
      str_remove_all(., pattern = "_拍照") %>%
      str_replace_all(., pattern = "tif", replacement = "jpg") 
    
    meta_df = rbind(meta_df, 
                      c(image_name, Line, Replicate, Fov))
    
    }
  return(meta_df)
}


## define a plot theme
theme_vignette <- function(base_size=14, base_family="Arial") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5), 
            plot.subtitle = element_text(hjust = 0.5),
            plot.caption = element_text(hjust = 0),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(1, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")))
}

theme_vignette_color <- function(...){
  library(scales)
  discrete_scale(
    aesthetics = c("colour", "color"), scale_name = "theme_vignette",
    palette = manual_pal(values = unname(colorScale)), ...)
}

theme_vignette_fill <- function(...){
  library(scales)
  discrete_scale(
    aesthetics = c("fill"), scale_name = "theme_vignette",
    palette = manual_pal(values = unname(colorScale)), ...)
}


# This function gets angles between row vectors
get_angles = function(row_vectors) {
  norm_vectors = t(apply(
    row_vectors, 1, function(x) x / sqrt(sum(x**2, na.rm = TRUE))))
  dotprod = norm_vectors %*% t(norm_vectors)
  diag(dotprod) = 1
  # return(acos(dotprod) * 180/pi)
  return(dotprod)
}

get_target_pathway = function(profiles) {
  target = c()
  pathway = c()
  for (idx in seq_along(rownames(profiles))) {
    if(str_sub(rownames(profiles)[idx], start = 9) %in% layout$Drug){
      pos = which(layout$Drug == str_sub(rownames(profiles)[idx], start = 9))
      target = append(target, layout[pos ,2])
      pathway = append(pathway, layout[pos, 3])
    }
  }
  result = rbind(target, pathway)
  return(result)
}


colorScale = stats::setNames(
  object = c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231"), 
  nm = c('CK33P02', 'CK38P01', 'CK40P01', 'CK44P01', 'CK51P06'))

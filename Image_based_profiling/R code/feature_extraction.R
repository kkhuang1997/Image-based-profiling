library("EBImage")
library("magrittr")
library("tibble")
library("ggplot2")
library("genefilter")
library("GGally")
library("MaxContrastProjection")
library("Rfast")

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
  img_th = combine( mapply(function(frame, th) frame < th, getFrames(max_proj), threshold, SIMPLIFY=FALSE) )
  ## morphology manipulation
  imgOpened = EBImage::opening(img_th, kern = makeBrush(11, shape = "disc"))
  imgClosed = EBImage::closing(imgOpened, kern = makeBrush(11,shape = "disc"))
  imgSeed = bwlabel(imgClosed)
  return(imgSeed)
}


# Calculate features for every image
Fea_extrct <- function(seg, ref, org_line, image_field, capture_date){
  num_stack = length(z_stack)
  features = vector(mode = "list", length = num_stack)
  for(fld in seq_len(num_stack)){
    img = Image(data = ref[,,fld])
    seg_map = Image(data = seg[,,fld])
    mask = matrix(
      as.integer(seg[,,fld] > 0), 
      nrow = nrow(seg), 
      ncol = ncol(seg))  ## set foreground 1
    
    # Small erosion to make the watershedding of touching objects better
    mask_erode = erode(mask, kern = makeBrush(11, "disc"))
    
    # Find distance to background and perform primitive watershedding. The 
    # tolerance of 15 was chosen empirically and is bound to change as the 
    # segmentation gets better
    dmap = distmap(mask_erode, "euclidean")
    wshed = watershed(dmap, tolerance = 15)
    
    # Undo the erosion by voronoi propagating the watershed labels into the 
    # original mask
    mask_labeled = propagate(x = mask, seeds = wshed, mask = mask)
    
    # This is a hacky solution to avoid errors due to missing segmentation.
    # WARNING: If anybody changes the number of features to compute (haralick scales), 
    # then the size of the array must be modified.
    if(sum(mask_labeled != 0) == 0) {
      features[[fld]] = matrix(0, nrow = 0, ncol = 1572)
      # features_clumps[[fld]] = matrix(0, nrow = 0, ncol = 1572)
      # features_noseg[[fld]] = matrix(0, nrow = 0, ncol = 118)
      next
    }
    
    ## processed segmap
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


Fea_extrct_jpg <- function(seg, ref, org_line, image_field, capture_date){
  num_stack = length(z_stack)
  features = vector(mode = "list", length = num_stack)
  for(fld in seq_len(num_stack)){
    img = Image(data = ref)
    seg_map = Image(data = seg)
    mask = matrix(
      as.integer(seg > 0), 
      nrow = nrow(seg), 
      ncol = ncol(seg))  ## set foreground 1
    
    # Small erosion to make the watershedding of touching objects better
    mask_erode = erode(mask, kern = makeBrush(11, "disc"))
    
    # Find distance to background and perform primitive watershedding. The 
    # tolerance of 15 was chosen empirically and is bound to change as the 
    # segmentation gets better
    dmap = distmap(mask_erode, "euclidean")
    wshed = watershed(dmap, tolerance = 15)
    
    # Undo the erosion by voronoi propagating the watershed labels into the 
    # original mask
    mask_labeled = propagate(x = mask, seeds = wshed, mask = mask)
    
    # This is a hacky solution to avoid errors due to missing segmentation.
    # WARNING: If anybody changes the number of features to compute (haralick scales), 
    # then the size of the array must be modified.
    if(sum(mask_labeled != 0) == 0) {
      features[[fld]] = matrix(0, nrow = 0, ncol = 1572)
      # features_clumps[[fld]] = matrix(0, nrow = 0, ncol = 1572)
      # features_noseg[[fld]] = matrix(0, nrow = 0, ncol = 118)
      next
    }
    
    ## processed segmap
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


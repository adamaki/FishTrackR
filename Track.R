# Lumpfish discrimination learning project
# Adam Brooker
# 24th June 2019


library(EBImage)
library(imager) 

#workingdir <- 'G:/Projects/Cleaner fish learning/Data/Test1'
workingdir <- ifelse(Sys.info()['user'] == 'Laptop', 'G:/Projects/Cleaner fish learning/Data/Test1', '/Users/adambrooker/R Projects/FishTrackR/ShortTest') # change to location of data

setwd(workingdir)

inputfile <- 'Week3Test'

# load image set, modify and save output
files <- list.files(path = workingdir, pattern = inputfile, all.files = FALSE, recursive = FALSE)

# select green channel
for(i in 1:length(files)){
  image <- readImage(files[i]) # load image file
  image <- channel(image, 'green') # Use green channel for max contrast against green tank
  image <- image*(1/max(image)) # increase contrast to max
  image <- (rotate(image, 22)[626:1730, 491:1600]) #change values to rotate and crop image
  image <- gblur(image, sigma = 3) # Gaussian smoothing filter
  writeImage(image, gsub('.jpg', '_green.jpg', files[i]))
}

# select blue channel
for(i in 1:length(files)){
  image <- readImage(files[i]) # load image file
  image <- channel(image, 'blue') # Use green channel for max contrast against green tank
  image <- image*(1/max(image)) # increase contrast to max
  image <- (rotate(image, 22)[626:1730, 491:1600]) #change values to rotate and crop image
  image <- gblur(image, sigma = 3) # Gaussian smoothing filter
  writeImage(image, gsub('.jpg', '_blue.jpg', files[i]))
}

# Load modified files as image stack
modfiles <- list.files(path = workingdir, pattern = '_blue', all.files = FALSE, recursive = FALSE)
image_stack <- readImage(modfiles, all = T)

stack_mean <- as.Image(rowMeans(image_stack, dims = 2)) # create mean image of stack

# subject background to leave fish
thresh_stack <- stack_mean

for(j in 1:length(modfiles)){
subimg <- stack_mean - image_stack[,,j]
subimg <- subimg*(1/max(subimg))
subimg <- subimg > 0.6
thresh_stack <- combine(thresh_stack, subimg)
}

thresh_stack <- thresh_stack[,,-c(1)] # remove first image, which is stack_mean

# segment images in z stack
for(k in 1:length(modfiles)){
  thresh_stack[,,k] <- bwlabel(thresh_stack[,,k])
}

display(colorLabels(thresh_stack), method = 'raster', all = T)

# adaptive thresholding techniques for detecting objects

# Laplacian high-pass filter for detecting edges
fhi <- matrix(1, nrow = 3, ncol = 3)
fhi[2, 2] <- -7.5
image_fhi <- filter2(image, fhi)
display(image_fhi)

#otsu filter
display(subimg < otsu(subimg))

# adaptive thresholding
disc <- makeBrush(121, 'disc')
disc <- disc/sum(disc)
offset <- 0.00001
image_bg <- filter2(subimg, disc)
image_th <- subimg < image_bg + offset
display(image_th, all = T)

# adaptive thresholding
image_th <- thresh(subimg, w = 30, h = 30, offset = 0.1)
display(1-image_th)

nmask <- thresh(subimg, w = 15, h = 15, offset = 0.0001)
nmask <- 1-nmask
nmask <- erode(nmask, makeBrush(9, shape = 'disc'))

nmask <- watershed(distmap(nmask), 1.5, 1)

ft <- computeFeatures.shape(nmask) ## need these for edge and perimeter
mf <- computeFeatures.moment(nmask, image) ## need these for intensity and size

for (i in seq_along(ft)) { 
  ft[[i]] <- cbind(ft[[i]], mf[[i]]) 
}

ft <- cbind(ft, mf)

nucimg <- rmObjects(nmask, 
                lapply(ft[,'s.area'],
                function(x)
                {which(x < 150 | x > 1000)}
                #{which(x[,'s.area'] < 150 | x[,'s.area'] > 10000)}
                ))


display(colorLabels(bwlabel(nmask)))
display(nmask)






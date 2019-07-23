# Lumpfish discrimination learning project
# Adam Brooker
# 24th June 2019


library(EBImage)
#library(imager) 
library(dplyr)

# functions (run once before using script)----------------------------------------------------------------------

# 1. read integer from command line input
readinteger <- function(message)
{ 
  num <- readline(prompt = message)
  if(!grepl("^[0-9]+$", num))
  {
    return(readinteger())
  }
  
  #return(as.integer(num))
  num <<- as.integer(num)
}



#workingdir <- 'G:/Projects/Cleaner fish learning/Data/Test1'
workingdir <- ifelse(Sys.info()['user'] == 'Laptop', 'G:/Projects/Cleaner fish learning/Data/Test2/subset2', '/Users/adambrooker/R Projects/FishTrackR/ShortTest') # change to location of data

setwd(workingdir)

inputfile <- 'Week3Test'
rotangle <- 22
xrange <- c(626,1730)
yrange <- c(491, 1600)

# load image set, modify and save output
files <- list.files(path = workingdir, pattern = inputfile, all.files = FALSE, recursive = FALSE)

# select green channel
#for(i in 1:length(files)){
#  image <- readImage(files[i]) # load image file
#  image <- channel(image, 'green') # Use green channel for max contrast against green tank
#  image <- (rotate(image, 22)[626:1730, 491:1600]) #change values to rotate and crop image
#  image <- image*(1/max(image)) # increase contrast to max  
#  image <- gblur(image, sigma = 3) # Gaussian smoothing filter
#  writeImage(image, gsub('.jpg', '_green.jpg', files[i]))
#}

# select blue channel
#for(i in 1:length(files)){
#  image <- readImage(files[i]) # load image file
#  image <- channel(image, 'blue') # Use blue channel for max contrast against green tank
#  image <- (rotate(image, 22)[626:1730, 491:1600]) #change values to rotate and crop image
#  image <- image*(1/max(image)) # increase contrast to max
#  image <- gblur(image, sigma = 3) # Gaussian smoothing filter
#  writeImage(image, gsub('.jpg', '_blue.jpg', files[i]))
#}

image_stack <- readImage(files, all = T)
mod_stack <- (rotate(image_stack, rotangle)[xrange[[1]]:xrange[[2]], yrange[[1]]:yrange[[2]],,])
blue_stack <- channel(mod_stack, 'blue')
blue_stack <- blue_stack*(1/max(blue_stack)) # increase contrast to max
blue_stack <- gblur(blue_stack, sigma = 3) # Gaussian smoothing filter

# Load modified files as image stack
#modfiles <- list.files(path = workingdir, pattern = '_blue', all.files = FALSE, recursive = FALSE)
#image_stack <- readImage(modfiles, all = T)
stack_mean <- as.Image(rowMeans(blue_stack, dims = 2)) # create mean image of stack

# subtract background to leave fish
thresh_stack <- stack_mean # seed thresholded image stack

for(j in 1:dim(blue_stack)[[3]]){
subimg <- blue_stack[,,j] - stack_mean
subimg <- subimg*(1/max(subimg))
subimg <- subimg > 0.75
thresh_stack <- EBImage::combine(thresh_stack, subimg)
}

thresh_stack <- thresh_stack[,,-c(1)] # remove first image, which is stack_mean

display(thresh_stack, method = 'raster', all = T)  

# segment images in z stack
thresh_stack <- bwlabel(thresh_stack)

# Remove small objects (noise) and large objects (people) from images
for(k in 1:dim(thresh_stack)[[3]]){
  sf <- computeFeatures.shape(thresh_stack[,,k], blue_stack[,,k]) # calculate shape features
  thresh_stack[,,k] <- rmObjects(thresh_stack[,,k], which(sf[,'s.area'] < 40 | sf[,'s.area'] > 2000))
}

# fill holes in fish and dilate to join gaps
thresh_stack <- fillHull(thresh_stack)
thresh_stack <- dilate(thresh_stack, makeBrush(11, shape = 'disc'))
thresh_stack <- bwlabel(thresh_stack)

# create coords list file
coords <- data.frame(frame = numeric(), fishx = numeric(), fishy = numeric(), errors = character())

# find centre of each fish object & write coords to file for each frame
for(m in 1:dim(thresh_stack)[[3]]){
  mf <- computeFeatures.moment(thresh_stack[,,m], blue_stack[,,m]) # calculate moment features  
  coords <- add_row(coords, frame = m, fishx = ifelse(is.matrix(get('mf')), round(mf[1,1]), NA), 
                    fishy = ifelse(is.matrix(get('mf')), round(mf[1,2]), NA), 
                    errors = ifelse(is.matrix(get('mf')) == F, 'No fish', ifelse(nrow(mf) > 1, 'Noise', 'None')))
}

display(thresh_stack, method = 'raster', all = T)   

# Run through errors and fix
for(p in 1:dim(thresh_stack)[[3]]){
  
 if(coords[p,4] == 'Noise'){
   mf <- computeFeatures.moment(thresh_stack[,,p], blue_stack[,,p]) # calculate moment features 
   obj <- data.frame(object = rownames(mf), fish_x = round(mf[,'m.cx']), fish_y = round(mf[,'m.cy']))
   print(obj)
   overlay <- paintObjects(thresh_stack[,,p], mod_stack[,,,p], col = c('#ff00ff', 'blue'), thick = T, closed = T, opac = c(1, 0))
   display(overlay) # display noisy image with colour labels
   readinteger('Enter oject number to remove: ') # input object number to remove using interactive image display
   thresh_stack[,,p] <- rmObjects(thresh_stack[,,p], num)
 } else {
   if(coords[p, 4] == 'No fish'){
     display(mod_stack[,,,p])
     Sys.sleep(0)
     readinteger('Enter x coordinate of fish: ')
     fish <- num
     readinteger('Enter y coordinate of fish: ')
     fish <- c(fish, num)
     thresh_stack[fish[[1]], fish[[2]], p] <- 1
     thresh_stack[,,p] <- dilate(thresh_stack[,,p], makeBrush(15, shape = 'disc'))
   }
 } 
}


# segment images in z stack
thresh_stack <- bwlabel(thresh_stack)

# find centre of each fish object & write coords to file for each frame
coords <- data.frame(frame = numeric(), fishx = numeric(), fishy = numeric(), errors = character())

for(m in 1:dim(thresh_stack)[[3]]){
  mf <- computeFeatures.moment(thresh_stack[,,m], blue_stack[,,m]) # calculate moment features  
  coords <- add_row(coords, frame = m, fishx = ifelse(is.matrix(get('mf')), round(mf[1,1]), NA), 
                    fishy = ifelse(is.matrix(get('mf')), round(mf[1,2]), NA), 
                    errors = ifelse(is.matrix(get('mf')) == F, 'No fish', ifelse(nrow(mf) > 1, 'Noise', 'None')))
}



# save segmented images as series
dir.create('Tracked')
setwd(paste0(workingdir, '/Tracked'))

for(n in 1:dim(thresh_stack)[[3]]){
  overlay <- paintObjects(thresh_stack[,,n], mod_stack[,,,n], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T)
  writeImage(overlay, gsub('.jpg', '_tracked.png', files[n]))
}


# Manual removal of noise and remarking of fish
frame <- 1 # change to frame number to remark
while(frame > 0){
  mf <- computeFeatures.moment(thresh_stack[,,frame], blue_stack[,,frame]) # calculate moment features  
  thresh_stack[,,frame] <- rmObjects(thresh_stack[,,frame], seq(1, nrow(mf), 1))
  display(mod_stack[,,,frame])
  Sys.sleep(5)
  readinteger('Enter x coordinate of fish: ')
  Sys.sleep(0.1)
  fish <- num
  readinteger('Enter y coordinate of fish: ')
  Sys.sleep(0.1)
  fish <- c(fish, num)
  thresh_stack[fish[[1]], fish[[2]], frame] <- 1
  thresh_stack[,,frame] <- dilate(thresh_stack[,,frame], makeBrush(15, shape = 'disc'))
  coords[frame,2] <- fish[[1]]
  coords[frame,3] <- fish[[2]]
  coords[frame,4] <- 'None'
  setwd(paste0(workingdir, '/Tracked'))
  overlay <- paintObjects(thresh_stack[,,frame], mod_stack[,,,frame], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T)
  writeImage(overlay, gsub('.jpg', '_tracked.png', files[frame]))
  frame <- 0
}


write.csv(coords, paste0(inputfile, '.csv'), row.names = F)


# Load segmented files as image stack
files <- list.files(path = paste0(workingdir, '/Tracked'), pattern = '_tracked.png', all.files = FALSE, recursive = FALSE)
thresh_stack <- readImage(files, all = T)
thresh_stack <- thresh_stack[,,1,] == 1 & thresh_stack[,,2,] == 0 & thresh_stack[,,3,] == 1 # threshold pixels with fuschia colour
colorMode(thresh_stack) <- Grayscale
thresh_stack <- fillHull(thresh_stack)
thresh_stack <- bwlabel(thresh_stack)









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






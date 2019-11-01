# Lumpfish discrimination learning project
# Adam Brooker
# 24th June 2019


library(EBImage)
#library(imager) 
library(dplyr)

# functions (run once before using script)----------------------------------------------------------------------

# read integer from command line input
readinteger <- function(message = 'Not a number')
{ 
  num <- readline(prompt = message)
  if(!grepl("^[0-9]+$", num))
  {
    return(readinteger())
  }
  
  #return(as.integer(num))
  num <<- as.integer(num)
}



# 1. Set working directory and input variables----------------------------
workingdir <- ifelse(Sys.info()['user'] == 'Laptop', 'G:/Projects/Cleaner fish learning/Data/Test3/101-300', '/Users/adambrooker/R Projects/FishTrackR/ShortTest') # change to location of data
setwd(workingdir)

modout <- bwlabel(channel(readImage(file.choose()), 'gray')) # choose modeloutling.png image

inputfile <- 'week3test1small'
rotangle <- 22 # image rotation angle to translate image to cartesian grid
xrange <- c(313,865) # x-axis crop dimensions
yrange <- c(246, 800) # y-axis crop dimensions
centre <- c(268, 264) # coords for centre of tank in cropped and rotated image
inrad <- 40 # radius of drain mask in pixels
outrad <- 240 # radius of tank mask in pixels
cal1 <- c(397, 164) # location of 1st calibration marker
cal2 <- c(573, 567) # location of 2nd calibration marker
caldist <- 100 # real distance between calibration markers
modpiv.l <- c(182, 261) # location of left model pivot axis in cropped and rotated image
modpiv.r <- c(385, 281) # location of right model pivot axis in cropped and rotated image

# calculated variables
cfactor <- 100/round(sqrt(abs(cal1[[1]]-cal2[[1]])^2+abs(cal1[[2]]-cal2[[2]])^2)) # calculate real distance conversion factor

# 2. load image set, modify images for analysis and create mean of image stack---------------------------------
files <- list.files(path = workingdir, pattern = inputfile, all.files = FALSE, recursive = FALSE)

mod_stack <- readImage(files[1])
mod_stack <- (rotate(mod_stack, rotangle)[xrange[[1]]:xrange[[2]], yrange[[1]]:yrange[[2]],])
#mod_stack <- channel(mod_stack, 'blue')

for(i in 2:length(files)){
  image <- readImage(files[i]) # load image file
  mod_stack <- EBImage::combine(mod_stack, rotate(image, rotangle)[xrange[[1]]:xrange[[2]], yrange[[1]]:yrange[[2]],])
  #mod_stack <- EBImage::combine(mod_stack, channel(image, 'blue'))
}

# create blue stack for tracking fish
blue_stack <- channel(mod_stack, 'blue')
blue_stack <- blue_stack*(1/max(blue_stack)) # increase contrast to max
blue_stack <- gblur(blue_stack, sigma = 3) # Gaussian smoothing filter
blue_stack_mean <- as.Image(rowMeans(blue_stack, dims = 2)) # create mean image of stack

# create rgb stack for tracking model

# extract red channel and modify (not working. trying green channel)
red_stack <- channel(mod_stack, 'red')
red_stack <- 1-red_stack # negative
red_stack <- red_stack*(1/max(red_stack)) # increase contrast to max
red_stack <- gblur(red_stack, sigma = 1) # Gaussian smoothing filter

green_stack <- channel(mod_stack, 'green')
green_stack <- 1-green_stack
green_stack <- green_stack*(1/max(green_stack)) # increase contrast to max
green_stack <- gblur(green_stack, sigma = 1) # Gaussian smoothing filter

rgb_stack <- red_stack*green_stack*blue_stack
#display(rgb_stack*(1/max(rgb_stack)))



# 3. Subtract mean image from stack and threshold image stack-----------------------------------------

# subtract background to leave fish
thresh_stack <- blue_stack_mean # seed thresholded image stack

for(j in 1:dim(blue_stack)[[3]]){
subimg <- blue_stack[,,j] - blue_stack_mean
subimg <- subimg*(1/max(subimg))
subimg <- subimg > 0.75
thresh_stack <- EBImage::combine(thresh_stack, subimg)
}

thresh_stack <- thresh_stack[,,-c(1)] # remove first image, which is stack_mean

display(thresh_stack, method = 'raster', all = T)  

# segment images in z stack
thresh_stack <- bwlabel(thresh_stack)

# 4. Create mask and remove object noise--------------------------------------------------

# create mask of tank area
tmask <- floodFill(mod_stack[,,,1], c(1, 1, 1, 1), col = 'black', tolerance = 255)
tmask <- drawCircle(tmask, centre[[1]], centre[[2]], radius = outrad, fill = T, col = 'white')
tmask <- drawCircle(tmask, centre[[1]], centre[[2]], radius = inrad, fill = T, col = 'black')
tmask <- channel(tmask, 'red')
tmask <- bwlabel(tmask)

# Remove small objects (noise), large objects (people) and objects outside mask from images
for(k in 1:dim(thresh_stack)[[3]]){
  sf <- computeFeatures.shape(thresh_stack[,,k], blue_stack[,,k]) # calculate shape features
  mf <- computeFeatures.moment(thresh_stack[,,k], blue_stack[,,k]) # calculate moment features 
  if(is.matrix(get('mf'))) {
    mf <- as.data.frame(mf)
    mf$tmask <- ifelse(tmask[matrix(data = c(round(mf[,'m.cx']), round(mf[,'m.cy'])), nrow(mf))] == 1, 1, 0) # test whether object centre coords are in mask
    }
  thresh_stack[,,k] <- rmObjects(thresh_stack[,,k], which(sf[,'s.area'] < 40 | sf[,'s.area'] > 2000)) # remove small and large objects
  thresh_stack[,,k] <- rmObjects(thresh_stack[,,k], which(mf[,'mask'] == 0)) # remove objects outside mask
}

# fill holes in fish and dilate to join gaps
thresh_stack <- fillHull(thresh_stack)
thresh_stack <- dilate(thresh_stack, makeBrush(11, shape = 'disc'))
thresh_stack <- bwlabel(thresh_stack)

# 5. Create list of fish coordinates and mark errors-----------------------------------------

# create coords list file
coords <- data.frame(frame = numeric(), fishpx = numeric(), fishpy = numeric(), errors = character())

# find centre of each fish object & write coords to file for each frame
for(m in 1:dim(thresh_stack)[[3]]){
  mf <- computeFeatures.moment(thresh_stack[,,m], blue_stack[,,m]) # calculate moment features  
  coords <- add_row(coords, frame = m, fishpx = ifelse(is.matrix(get('mf')), round(mf[1,1]), NA), 
                    fishpy = ifelse(is.matrix(get('mf')), round(mf[1,2]), NA), 
                    errors = ifelse(is.matrix(get('mf')) == F, 'No fish', ifelse(nrow(mf) > 1, 'Noise', 'None')))
}

coords$fishrx <- round(coords$fishpx * cfactor, 2)
coords$fishry <- round(coords$fishpy * cfactor, 2)
coords <- coords[,c(1, 2, 3, 5, 6, 4)]

display(thresh_stack, method = 'raster', all = T)   


# 6. save segmented images as series---------------------------------------------------
dir.create('Tracked')
setwd(paste0(workingdir, '/Tracked'))

for(n in 1:dim(thresh_stack)[[3]]){
  overlay <- paintObjects(thresh_stack[,,n], mod_stack[,,,n], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T)
  overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
  overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
  writeImage(overlay, gsub('.jpg', '_tracked.png', files[n]))
}


# 7. Run through errors and fix, saving new images (view images in external viewer as R viewer doesn't update in loops)----------------------------------
for(p in 1:dim(thresh_stack)[[3]]){
  
  if(coords[p,6] == 'Noise'){
    thresh_stack <- bwlabel(thresh_stack)
    mf <- computeFeatures.moment(thresh_stack[,,p], blue_stack[,,p]) # calculate moment features 
    obj <- data.frame(object = rownames(mf), fish_x = round(mf[,'m.cx']), fish_y = round(mf[,'m.cy']))
    print(obj)
    print(paste0('Frame: ', coords[p, 1]))
    #overlay <- paintObjects(thresh_stack[,,p], mod_stack[,,,p], col = c('#ff00ff', 'blue'), thick = T, closed = T, opac = c(1, 0))
    #display(overlay) # display noisy image with colour labels
    readinteger('Enter oject number to remove: ') # input object number to remove using interactive image display
    thresh_stack[,,p] <- rmObjects(thresh_stack[,,p], num)
    overlay <- paintObjects(thresh_stack[,,p], mod_stack[,,,p], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T)
    overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
    overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
    writeImage(overlay, gsub('.jpg', '_tracked.png', files[p]))
  } else {
    if(coords[p, 6] == 'No fish'){
      #display(mod_stack[,,,p])
      print(paste0('Frame: ', coords[p, 1]))
      readinteger('Enter x coordinate of fish: ')
      fish <- num
      readinteger('Enter y coordinate of fish: ')
      fish <- c(fish, num)
      thresh_stack[fish[[1]], fish[[2]], p] <- 1
      thresh_stack[,,p] <- dilate(thresh_stack[,,p], makeBrush(15, shape = 'disc'))
      overlay <- paintObjects(thresh_stack[,,p], mod_stack[,,,p], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T)
      overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
      overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
      writeImage(overlay, gsub('.jpg', '_tracked.png', files[p]))
    }
  } 
}


# 8. Manual removal of noise and remarking of fish for any remaining errors---------------------------------------------------
frame <- 95 # change to frame number to remark
thresh_stack <- bwlabel(thresh_stack)
while(frame > 0){
  mf <- computeFeatures.moment(thresh_stack[,,frame], blue_stack[,,frame]) # calculate moment features  
  if(is.matrix(get('mf'))) {thresh_stack[,,frame] <- rmObjects(thresh_stack[,,frame], seq(1, nrow(mf), 1))}
  display(mod_stack[,,,frame])
  #Sys.sleep(5)
  readinteger('Enter x coordinate of fish: ')
  #Sys.sleep(0.1)
  fish <- num
  readinteger('Enter y coordinate of fish: ')
  #Sys.sleep(0.1)
  fish <- c(fish, num)
  thresh_stack[fish[[1]], fish[[2]], frame] <- 1
  thresh_stack[,,frame] <- dilate(thresh_stack[,,frame], makeBrush(15, shape = 'disc'))
  coords[frame,2] <- fish[[1]]
  coords[frame,3] <- fish[[2]]
  coords[frame,6] <- 'None'
  setwd(paste0(workingdir, '/Tracked'))
  overlay <- paintObjects(thresh_stack[,,frame], mod_stack[,,,frame], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T)
  overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
  overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
  writeImage(overlay, gsub('.jpg', '_tracked.png', files[frame]))
  display(overlay)
  frame <- 0
}


# 9. model tracking code----------------------------------------------------------------------------------

model_stack <- rgb_stack

# adaptive threshold filter
disc <- makeBrush(31, 'disc')
disc <- disc / sum(disc)
offset <- 0.001

# create mask of model radii
mmask <- floodFill(mod_stack[,,,1], c(1, 1, 1, 1), col = 'black', tolerance = 255)
mmask <- drawCircle(mmask, modpiv.l[[1]], modpiv.l[[2]], radius = 22/cfactor, fill = T, col = 'white')
mmask <- drawCircle(mmask, modpiv.r[[1]], modpiv.r[[2]], radius = 22/cfactor, fill = T, col = 'white')
mmask <- drawCircle(mmask, centre[[1]], centre[[2]], radius = inrad, fill = T, col = 'black')
mmask <- channel(mmask, 'red')

#test_stack <- test
#mf.all <- data.frame()
#sf.all <- data.frame()

for(s in 1:dim(model_stack)[[3]]){
  
  test <- model_stack[,,s]
  
  test <- test > filter2(test, disc) + offset
  
  for(b in 1:10){
    test <- opening(test, makeBrush(5, shape = 'disc'))
  }
  
  test <- test * mmask
  
  test <- bwlabel(test)
  sf <- computeFeatures.shape(test, rgb_stack[,,1])
  test <- rmObjects(test, which(sf[,'s.area'] < 1000 | sf[,'s.area'] > 2000)) # keep objects 1000-2000 pixels area
  
  test <- watershed(distmap(test), 5)
  #display(colorLabels(test, normalize = T))
  
  sf <- computeFeatures.shape(test, rgb_stack[,,1])
  test <- rmObjects(test, which(sf[,'s.area'] < 1000 | sf[,'s.area'] > 2000))  # keep objects 1000-2000 pixels area
  
  mf <- computeFeatures.moment(test, rgb_stack[,,1])
  test <- rmObjects(test, which(mf[,'m.eccentricity'] < 0.96))  # remove objects that are not nearly a straight line
  
  mf <- computeFeatures.moment(test, rgb_stack[,,1])
  #test <- rmObjects(test, which( max(abs(mf[,'m.cx']-modpiv.l[[1]])+abs(mf[,'m.cy']-modpiv.l[[2]])) & mf[,'m.cx']<dim(test)[[1]]/2 | max(abs(mf[,'m.cx']-modpiv.r[[1]])+abs(mf[,'m.cy']-modpiv.r[[2]])) & mf[,'m.cx']>dim(test)[[1]]/2 ))
  while(nrow(mf)>2){
    test <- rmObjects(test, which.max( abs(mf[,'m.cx']-modpiv.l[[1]])+abs(mf[,'m.cy']-modpiv.l[[2]]) + abs(mf[,'m.cx']-modpiv.r[[1]])+abs(mf[,'m.cy']-modpiv.r[[2]]) ))
    mf <- computeFeatures.moment(test, rgb_stack[,,1])
  }
  
  mf <- computeFeatures.moment(test, rgb_stack[,,1])
  
  #display(paintObjects(mods, mod_stack[,,,1], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T))
  
  modang.l <- mf[which.min(mf[,'m.cx']), 'm.theta'] # select left model angle by min x coord
  modang.r <- mf[which.max(mf[,'m.cx']), 'm.theta'] # select right model angle by max x coord
  
  mod.l <- rotate(modout, modang.l*180/pi) # rotate left model by detected angle
  mod.r <- rotate(modout, 180+(modang.r*180/pi)) # rotate right model by detected angle
  
  mod.l <- translate(mod.l, c(modpiv.l[[1]] - (dim(mod.l)[[1]]/2), modpiv.l[[2]] - (dim(mod.l)[[2]]/2)))[0:dim(red_stack)[[1]], 0:dim(red_stack)[[2]]] # translate left model to correct position and crop
  mod.r <- translate(mod.r, c(modpiv.r[[1]] - (dim(mod.r)[[1]]/2), modpiv.r[[2]] - (dim(mod.r)[[2]]/2)))[0:dim(red_stack)[[1]], 0:dim(red_stack)[[2]]] # translate left model to correct position and crop
  
  #test_stack <- EBImage::combine(test_stack, test)
  #mf.all <- rbind(mf.all, mf)
  #sf.all <- rbind(sf.all, sf)
  
  model_stack[,,s] <- bwlabel(mod.l + mod.r) # combine left and right model outlines
  
  
}

#test_stack <- test_stack[,,-c(1)] 



# 10. Recalculate coordinates after all errors fixed and models tracked-------------------------------
thresh_stack <- bwlabel(thresh_stack)
display(thresh_stack, method = 'raster', all = T)

# create coords list file
coords <- data.frame(frame = numeric(), fishpx = numeric(), fishpy = numeric(), distmod.pl = numeric(), distmod.pr = numeric(), errors = character())

# find centre of each fish object & write coords to file for each frame
for(m in 1:dim(thresh_stack)[[3]]){
  #mf <- computeFeatures.moment(thresh_stack[,,m], blue_stack[,,m]) # calculate moment features  

# extract left model
mf <- computeFeatures.moment(model_stack[,,m], rgb_stack[,,m])
mod.l <- rmObjects(model_stack[,,m], which(mf[,'m.cx'] > dim(model_stack)[[1]]/2))
mod.r <- rmObjects(model_stack[,,m], which(mf[,'m.cx'] < dim(model_stack)[[1]]/2))

# get centre of fish coords
mf <- computeFeatures.moment(thresh_stack[,,m], rgb_stack[,,m])
fishpx <- c(round(mf[,'m.cx']), round(mf[,'m.cy']))

# get min fish distance to left model
modpx.l <- as.data.frame(which(mod.l == 1, arr.ind = T))
modpx.l$dist <- round(sqrt(abs(fishpx[[1]]-modpx.l[,1])^2+abs(fishpx[[2]]-modpx.l[,2])^2))

# get min fish distance to right model
modpx.r <- as.data.frame(which(mod.r == 1, arr.ind = T))
modpx.r$dist <- round(sqrt(abs(fishpx[[1]]-modpx.r[,1])^2+abs(fishpx[[2]]-modpx.r[,2])^2))

coords <- add_row(coords, frame = m, fishpx = ifelse(is.matrix(get('mf')), round(mf[1,1]), NA), 
                fishpy = ifelse(is.matrix(get('mf')), round(mf[1,2]), NA),
                distmod.pl = min(modpx.l$dist),
                distmod.pr = min(modpx.r$dist),
                errors = ifelse(is.matrix(get('mf')) == F, 'No fish', ifelse(nrow(mf) > 1, 'Noise', 'None')))

}


# convert pixel coordinates to mm coordinates
coords$fishrx <- round(coords$fishpx * cfactor, 2)
coords$fishry <- round(coords$fishpy * cfactor, 2)
coords$distmod.rl <- round(coords$distmod.pl * cfactor, 2)
coords$distmod.rr <- round(coords$distmod.pr * cfactor, 2)
coords <- coords[,c(1, 2, 3, 5, 6, 7, 8, 4)]


write.csv(coords, paste0(inputfile, '.csv'), row.names = F)




# Load segmented files as image stack (If need to reload old dataset)-----------------------------------
files <- list.files(path = paste0(workingdir, '/Tracked'), pattern = '_tracked.png', all.files = FALSE, recursive = FALSE)
setwd(paste0(workingdir, '/Tracked'))
thresh_stack <- readImage(files, all = T)
thresh_stack <- thresh_stack[,,1,] == 1 & thresh_stack[,,2,] == 0 & thresh_stack[,,3,] == 1 # threshold pixels with fuschia colour
colorMode(thresh_stack) <- Grayscale
thresh_stack <- fillHull(thresh_stack)
thresh_stack <- bwlabel(thresh_stack)







# Old code not required for lumpfish tank videos---------------------------------------------------

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






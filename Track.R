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
workingdir <- '/Users/adambrooker/R Projects/FishTrackR/Data/Control/Week 1/Training/Trial 1/Test 1' # change to location of data
setwd(workingdir)

modout <- bwlabel(channel(readImage(file.choose()), 'gray')) # choose modeloutling.png image

inputfile <- 'GOPR0001'
files <- list.files(path = workingdir, pattern = inputfile, all.files = FALSE, recursive = FALSE)
start <- 504 # start frame number
end <- 603 # end frame number
rotangle <- 21 # image rotation angle to translate image to cartesian grid
xrange <- c(201,497) # x-axis crop dimensions
yrange <- c(151, 457) # y-axis crop dimensions

centre <- c(152, 139) # coords for centre of tank in cropped and rotated image
inrad <- 27 # radius of drain mask in pixels
outrad <- 163 # radius of tank mask in pixels
cal1 <- c(6, 150) # location of 1st calibration marker in any image
cal2 <- c(293, 147) # location of 2nd calibration marker in any image
caldist <- 100 # real distance between calibration markers in cm
modpiv.l <- c(83, 147) # location of left model pivot axis in cropped and rotated image
modpiv.r <- c(220, 142) # location of right model pivot axis in cropped and rotated image
curve.fitting <- T # Turn curve fitting for model rotation on or off
#dir.l <- F # direction of left model (switch to F if wrong direction in output)
#dir.r <- T # direction of right model (switch to F if wrong direction in output)

# image testing to refine rotating and cropping
test.img <- readImage(files[[start]])
display(rotate(test.img, rotangle)) # display rotated image to find crop coordinates
test.img <- (rotate(test.img, rotangle)[xrange[[1]]:xrange[[2]], yrange[[1]]:yrange[[2]],])
display(test.img)
test.img <- drawCircle(test.img, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
test.img <- drawCircle(test.img, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
display(test.img)
test.img <- drawCircle(test.img, modpiv.l[[1]], modpiv.l[[2]], radius = 2, col = 'light blue', fill = T)
test.img <- drawCircle(test.img, modpiv.r[[1]], modpiv.r[[2]], radius = 2, col = 'light blue', fill = T)
display(test.img)

rm(test.img)

# calculated variables
cfactor <- caldist/round(sqrt(abs(cal1[[1]]-cal2[[1]])^2+abs(cal1[[2]]-cal2[[2]])^2)) # calculate real distance conversion factor
modout <- resize(modout, w = round(dim(modout)[[1]]*(0.2288/cfactor))) > 0 # resize model image


# 2. load image set, modify images for analysis and create mean of image stack---------------------------------

system.time({

files <- list.files(path = workingdir, pattern = inputfile, all.files = FALSE, recursive = FALSE)

mod_stack <- readImage(files[start])
mod_stack <- (rotate(mod_stack, rotangle, bg.col = 0.5)[xrange[[1]]:xrange[[2]], yrange[[1]]:yrange[[2]],])
#mod_stack <- channel(mod_stack, 'blue')

for(i in (start+1):end){
  image <- readImage(files[i]) # load image file
  mod_stack <- EBImage::combine(mod_stack, rotate(image, rotangle, bg.col = 0.5)[xrange[[1]]:xrange[[2]], yrange[[1]]:yrange[[2]],])
  #mod_stack <- EBImage::combine(mod_stack, channel(image, 'blue'))
}

# create grey stack for tracking fish
gray_stack <- channel(mod_stack, 'gray')
#gray_stack <- channel(mod_stack, 'gray')
gray_stack <- gray_stack*(1/max(gray_stack)) # increase contrast to max
gray_stack <- gblur(gray_stack, sigma = 3) # Gaussian smoothing filter
gray_stack_mean <- as.Image(rowMeans(gray_stack, dims = 2)) # create mean image of stack

# create rgb stack for tracking model

# extract blue channel and modify
blue_stack <- channel(mod_stack, 'blue')
#blue_stack <- 1-blue_stack # negative
blue_stack <- blue_stack*(1/max(blue_stack)) # increase contrast to max
blue_stack <- gblur(blue_stack, sigma = 3) # Gaussian smoothing filter

# extract red channel and modify
red_stack <- channel(mod_stack, 'red')
red_stack <- 1-red_stack # negative
red_stack <- red_stack*(1/max(red_stack)) # increase contrast to max
red_stack <- gblur(red_stack, sigma = 3) # Gaussian smoothing filter

green_stack <- channel(mod_stack, 'green')
green_stack <- 1-green_stack
green_stack <- green_stack*(1/max(green_stack)) # increase contrast to max
green_stack <- gblur(green_stack, sigma = 3) # Gaussian smoothing filter

rgb_stack <- red_stack*green_stack*blue_stack
#display(rgb_stack*(1/max(rgb_stack)))

Sys.sleep(2)
#}) # end of system time measurement


# 3. model tracking code----------------------------------------------------------------------------------

#system.time({

model_stack <- rgb_stack
#model_stack <- red_stack

# adaptive threshold filter
disc <- makeBrush(31, 'disc')
disc <- disc / sum(disc)
offset <- 0.001

# model size thresholds based on image size
#mod.min <- round((dim(model_stack)[[1]] * dim(model_stack)[[2]]) * 0.00325) # constant is area of image that is model, i.e. 0.325%
mod.min <- round((dim(model_stack)[[1]] * dim(model_stack)[[2]]) * 0.004) # constant is area of image that is model, i.e. 0.325%
#mod.max <- round((dim(model_stack)[[1]] * dim(model_stack)[[2]]) * 0.0065)
mod.max <- round((dim(model_stack)[[1]] * dim(model_stack)[[2]]) * 0.015)

# create mask of model radii
mmask <- floodFill(mod_stack[,,,1], c(1, 1, 1, 1), col = 'black', tolerance = 255)
mmask <- drawCircle(mmask, modpiv.l[[1]], modpiv.l[[2]], radius = 22/cfactor, fill = T, col = 'white')
mmask <- drawCircle(mmask, modpiv.r[[1]], modpiv.r[[2]], radius = 22/cfactor, fill = T, col = 'white')
mmask <- drawCircle(mmask, centre[[1]], centre[[2]], radius = inrad, fill = T, col = 'black')
mmask <- channel(mmask, 'red')

modangles <- data.frame(ang.l = numeric(), ang.r = numeric())

for(s in 1:dim(model_stack)[[3]]){
  
  test <- model_stack[,,s]
  
  test <- test > filter2(test, disc) + offset
  
  for(b in 1:10){
    test <- opening(test, makeBrush(5, shape = 'disc'))
  }
  
  test <- test * mmask
  
  test <- bwlabel(test)
  #sf <- computeFeatures.shape(test, rgb_stack[,,1])
  #test <- rmObjects(test, which(sf[,'s.area'] < mod.min | sf[,'s.area'] > mod.max)) # keep model-sized objects
  
  test <- watershed(distmap(test), 5) # splits joined objects using watershed function
  #display(colorLabels(test, normalize = T))
  
  sf <- computeFeatures.shape(test, rgb_stack[,,1])
  test <- rmObjects(test, which(sf[,'s.area'] < mod.min | sf[,'s.area'] > mod.max))  # keep model-sized objects
  
  mf <- computeFeatures.moment(test, rgb_stack[,,1])
  test <- rmObjects(test, which(mf[,'m.eccentricity'] < 0.90))  # remove objects that are not nearly a straight line
  
  mf <- computeFeatures.moment(test, rgb_stack[,,1])
  #test <- rmObjects(test, which( max(abs(mf[,'m.cx']-modpiv.l[[1]])+abs(mf[,'m.cy']-modpiv.l[[2]])) & mf[,'m.cx']<dim(test)[[1]]/2 | max(abs(mf[,'m.cx']-modpiv.r[[1]])+abs(mf[,'m.cy']-modpiv.r[[2]])) & mf[,'m.cx']>dim(test)[[1]]/2 ))
  if(is.matrix(get('mf'))) {
  while(nrow(mf)>2){
    test <- rmObjects(test, which.max( abs(mf[,'m.cx']-modpiv.l[[1]])+abs(mf[,'m.cy']-modpiv.l[[2]]) + abs(mf[,'m.cx']-modpiv.r[[1]])+abs(mf[,'m.cy']-modpiv.r[[2]]) )) # remove extra objects furthest away from model pivots defined in setup
    mf <- computeFeatures.moment(test, rgb_stack[,,1])
  }
  }
  
  mf <- computeFeatures.moment(test, rgb_stack[,,1])
  
  modang.l <- mf[which.min(mf[,'m.cx']), 'm.theta'] # select left model angle by min x coord
  modang.r <- mf[which.max(mf[,'m.cx']), 'm.theta'] # select right model angle by max x coord
  
  modangles <- add_row(modangles, ang.l = modang.l*180/pi, ang.r = modang.r*180/pi) # add model angles to list of angles 
  
  # nested if/else to get model direction right
  #if(dir.l == T){
   # if(dir.r == T){
    #  modangles <- add_row(modangles, ang.l = modang.l*180/pi, ang.r = modang.r*180/pi) # add model angles to list of angles  
    #} else {
    #  modangles <- add_row(modangles, ang.l = modang.l*180/pi, ang.r = 180+(modang.r*180/pi)) # add model angles to list of angles  
    #}
    
  #} else {
  #  if(dir.r == T){
  #    modangles <- add_row(modangles, ang.l = 180+(modang.l*180/pi), ang.r = modang.r*180/pi) # add model angles to list of angles  
  #  } else {
  #    modangles <- add_row(modangles, ang.l = 180+(modang.l*180/pi), ang.r = 180+(modang.r*180/pi)) # add model angles to list of angles  
  #  }
  #}
  
  
} # end of model thresholding loop

# convert radians angles to degrees
#modangles$ang.l <- modangles$ang.l*180/pi
#modangles$ang.r <- 180+(modangles$ang.r*180/pi)

# polynomial curve fitting to smooth turn angles of models
if (curve.fitting == T){
  fit <- lm(modangles$ang.l ~ poly(as.numeric(rownames(modangles)), 10, raw = T))
  modangles$pred.l <- predict(fit, data.frame(x = seq(1, nrow(modangles), 1)))
  fit <- lm(modangles$ang.r ~ poly(as.numeric(rownames(modangles)), 10, raw = T))
  modangles$pred.r <- predict(fit, data.frame(x = seq(1, nrow(modangles), 1)))
} else {
  modangles$pred.l <- modangles$ang.l
  modangles$pred.r <- modangles$ang.r
}

for(t in 1:dim(model_stack)[[3]]){
  
  mod.l <- rotate(modout, modangles[t, 'pred.l']) # rotate left model by detected angle
  mod.r <- rotate(modout, modangles[t, 'pred.r']) # rotate right model by detected angle
  
  mod.l <- translate(mod.l, c(modpiv.l[[1]] - (dim(mod.l)[[1]]/2), modpiv.l[[2]] - (dim(mod.l)[[2]]/2)))[0:dim(red_stack)[[1]], 0:dim(red_stack)[[2]]] # translate left model to correct position and crop
  mod.r <- translate(mod.r, c(modpiv.r[[1]] - (dim(mod.r)[[1]]/2), modpiv.r[[2]] - (dim(mod.r)[[2]]/2)))[0:dim(red_stack)[[1]], 0:dim(red_stack)[[2]]] # translate left model to correct position and crop
  
  model_stack[,,t] <- bwlabel(mod.l + mod.r) # combine left and right model outlines
  
  #overlay <- paintObjects(thresh_stack[,,t], mod_stack[,,,t], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T) # paint fish outline in purple
  #overlay <- paintObjects(model_stack[,,t], overlay, col = c('light blue', 'light blue'), opac = c(1, 0), thick = T) # paint model outline in light blue
  #overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
  #overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
  #writeImage(overlay, gsub('.jpg', '_tracked.png', files[[t]]))
}

#test_stack <- test_stack[,,-c(1)] 
Sys.sleep(2)
#}) # end of system time measurement

# 4. Subtract mean image from stack, threshold image stack, mask outside tank and remove noise-----------------------------------------

#system.time({
  
thresh_stack <- gray_stack_mean # seed thresholded image stack

# create mask of tank area
tmask <- floodFill(mod_stack[,,,1], c(1, 1, 1, 1), col = 'black', tolerance = 255)
tmask <- drawCircle(tmask, centre[[1]], centre[[2]], radius = outrad, fill = T, col = 'white')
tmask <- drawCircle(tmask, centre[[1]], centre[[2]], radius = inrad, fill = T, col = 'black')
tmask <- channel(tmask, 'red')
tmask <- bwlabel(tmask)

for(j in 1:dim(gray_stack)[[3]]){
#subimg <- blue_stack[,,j] - blue_stack_mean
subimg <- gray_stack_mean - gray_stack[,,j]
subimg <- subimg*(1/max(subimg))
fmask <- dilate(model_stack[,,j], makeBrush(5, shape = 'disc')) # make mask for models
subimg <- 1-(1-subimg + fmask) # mask out models
subimg <- 1-(1-subimg + 1-tmask) # mask outside tank
subimg <- subimg > 0.75
subimg <- bwlabel(subimg)
sf <- computeFeatures.shape(subimg, gray_stack[,,j]) # calculate shape features
mf <- computeFeatures.moment(subimg, gray_stack[,,j]) # calculate moment features 
#if(is.matrix(get('mf'))) {
#  mf <- as.data.frame(mf)
#  mf$tmask <- ifelse(tmask[matrix(data = c(round(mf[,'m.cx']), round(mf[,'m.cy'])), nrow(mf))] == 1, 1, 0) # test whether object centre coords are in mask
#  subimg <- rmObjects(subimg, which(mf[,'tmask'] == 0)) # remove objects outside mask
#}
subimg <- rmObjects(subimg, which(sf[,'s.area'] < round((dim(thresh_stack)[[1]] * dim(thresh_stack)[[2]]) * 0.00005) | sf[,'s.area'] > round((dim(thresh_stack)[[1]] * dim(thresh_stack)[[2]]) * 0.005))) # remove objects less than 0.005% of image area and greater than 1% of image area 
thresh_stack <- EBImage::combine(thresh_stack, subimg)
}

thresh_stack <- thresh_stack[,,-c(1)] # remove first image, which is stack_mean

# segment images in z stack
#thresh_stack <- bwlabel(thresh_stack)

# fill holes in fish and dilate to join gaps
thresh_stack <- fillHull(thresh_stack)
thresh_stack <- dilate(thresh_stack, makeBrush(11, shape = 'disc'))
thresh_stack <- bwlabel(thresh_stack)

#display(thresh_stack, method = 'raster', all = T)  

Sys.sleep(2)
#}) # end of system time measurement

# 6. Create list of fish coordinates and mark errors-----------------------------------------

#system.time({
  
# create coords list file
coords <- data.frame(frame = numeric(), fishpx = numeric(), fishpy = numeric(), errors = character())

# find centre of each fish object & write coords to file for each frame
for(m in 1:dim(thresh_stack)[[3]]){
  mf <- computeFeatures.moment(thresh_stack[,,m], blue_stack[,,m]) # calculate moment features  
  coords <- add_row(coords, frame = m+(start-1), fishpx = ifelse(is.matrix(get('mf')), round(mf[1,1]), NA), 
                    fishpy = ifelse(is.matrix(get('mf')), round(mf[1,2]), NA), 
                    errors = ifelse(is.matrix(get('mf')) == F, 'No fish', ifelse(nrow(mf) > 1, 'Noise', 'None')))
}

coords$fishrx <- round(coords$fishpx * cfactor, 2)
coords$fishry <- round(coords$fishpy * cfactor, 2)
coords <- coords[,c(1, 2, 3, 5, 6, 4)]

rownames(coords) <- coords$frame # rename rows to frame number

display(thresh_stack, method = 'raster', all = T)   

Sys.sleep(2)
#}) # end of system time measurement

# 7. save segmented images as series---------------------------------------------------
#system.time({
  
dir.create('Tracked')
setwd(paste0(workingdir, '/Tracked'))

for(n in 1:dim(thresh_stack)[[3]]){
  overlay <- paintObjects(thresh_stack[,,n], mod_stack[,,,n], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T) # paint fish outline in purple
  overlay <- paintObjects(model_stack[,,n], overlay, col = c('light blue', 'light blue'), opac = c(1, 0), thick = T) # paint model outline in light blue
  overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
  overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
  writeImage(overlay, gsub('.jpg', '_tracked.png', files[n+(start-1)]))
}

}) # end of system time measurement



# 8. Run through errors and fix, saving new images (view images in external viewer as R viewer doesn't update in loops)----------------------------------
for(p in 1:dim(thresh_stack)[[3]]){
  
  if(coords[p,6] == 'Noise'){
    thresh_stack <- bwlabel(thresh_stack)
    mf <- computeFeatures.moment(thresh_stack[,,p], blue_stack[,,p]) # calculate moment features 
    obj <- data.frame(object = rownames(mf), fish_x = round(mf[,'m.cx']), fish_y = round(mf[,'m.cy']))
    print(obj)
    print(paste0('Frame: ', coords[p, 1]))
    #overlay <- paintObjects(thresh_stack[,,p], mod_stack[,,,p], col = c('#ff00ff', 'blue'), thick = T, closed = T, opac = c(1, 0))
    #display(overlay) # display noisy image with colour labels
    readinteger('Enter oject number to keep: ') # input object number to keep. Enter '0' if all objects are not fish.
    if(num == 0){
      thresh_stack[,,p] <- rmObjects(thresh_stack[,,p], seq(1, nrow(mf), 1)) # remove all objects if fish not detected
    } else{
      thresh_stack[,,p] <- rmObjects(thresh_stack[,,p], which(seq(1, nrow(mf), 1) != num)) # remove objects not fish
    }
    overlay <- paintObjects(thresh_stack[,,p], mod_stack[,,,p], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T)
    overlay <- paintObjects(model_stack[,,p], overlay, col = c('light blue', 'light blue'), opac = c(1, 0), thick = T) # paint model outline in light blue
    overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
    overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
    writeImage(overlay, gsub('.jpg', '_tracked.png', files[p+(start-1)]))
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
      overlay <- paintObjects(model_stack[,,p], overlay, col = c('light blue', 'light blue'), opac = c(1, 0), thick = T) # paint model outline in light blue
      overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
      overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
      writeImage(overlay, gsub('.jpg', '_tracked.png', files[p+(start-1)]))
    }
  } 
}


# 9. Manual removal of noise and remarking of fish for any remaining errors---------------------------------------------------
mrem <- function(fnum){

frame <- fnum-(start-1) # change to frame number to remark
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
  overlay <- paintObjects(model_stack[,,frame], overlay, col = c('light blue', 'light blue'), opac = c(1, 0), thick = T) # paint model outline in light blue
  overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
  overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
  writeImage(overlay, gsub('.jpg', '_tracked.png', files[fnum]))
  display(overlay)
  coords <<- coords
  thresh_stack <<- thresh_stack
  frame <- 0
}

}

# 10. Forcing fish coordinates if fish remains in the same place for many frames---------------------------------------------------
fixpos <- function(stf, enf){
  
  startframe <- stf-start+1
  endframe <- enf-start+1
  
  readinteger('Enter x coordinate of fish: ')
  fish <- num
  readinteger('Enter y coordinate of fish: ')
  fish <- c(fish, num)
  
  for(frame in startframe:endframe){
    thresh_stack <- bwlabel(thresh_stack)
    mf <- computeFeatures.moment(thresh_stack[,,frame], blue_stack[,,frame]) # calculate moment features  
    if(is.matrix(get('mf'))) {thresh_stack[,,frame] <- rmObjects(thresh_stack[,,frame], seq(1, nrow(mf), 1))}

    thresh_stack[fish[[1]], fish[[2]], frame] <- 1
    thresh_stack[,,frame] <- dilate(thresh_stack[,,frame], makeBrush(15, shape = 'disc'))
    coords[frame,2] <- fish[[1]]
    coords[frame,3] <- fish[[2]]
    coords[frame,6] <- 'None'
    
    setwd(paste0(workingdir, '/Tracked'))
    overlay <- paintObjects(thresh_stack[,,frame], mod_stack[,,,frame], col = c('#ff00ff', '#ff00ff'), opac = c(1, 0), thick = T)
    overlay <- paintObjects(model_stack[,,frame], overlay, col = c('light blue', 'light blue'), opac = c(1, 0), thick = T) # paint model outline in light blue
    overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = inrad, col = 'yellow', fill = F)
    overlay <- drawCircle(overlay, centre[[1]], centre[[2]], radius = outrad, col = 'yellow', fill = F)
    writeImage(overlay, gsub('.jpg', '_tracked.png', files[frame+start-1]))
  }
  
  coords <<- coords
  thresh_stack <<- thresh_stack

}

# 11. Recalculate coordinates after all errors fixed and models tracked-------------------------------
system.time({
thresh_stack <- bwlabel(thresh_stack)
model_stack <- bwlabel(model_stack)
#display(thresh_stack, method = 'raster', all = T)

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

coords <- add_row(coords, frame = m+(start-1), fishpx = ifelse(is.matrix(get('mf')), round(mf[1,1]), NA), 
                fishpy = ifelse(is.matrix(get('mf')), round(mf[1,2]), NA),
                distmod.pl = min(modpx.l$dist),
                distmod.pr = min(modpx.r$dist),
                errors = ifelse(is.matrix(get('mf')) == F, 'No fish', ifelse(nrow(mf) > 1, 'Noise', 'None')))

}


# convert pixel coordinates to mm coordinates
coords$fish.rx <- round(coords$fishpx * cfactor, 2)
coords$fish.ry <- round(coords$fishpy * cfactor, 2)
coords$distmod.rl <- round(coords$distmod.pl * cfactor, 2)
coords$distmod.rr <- round(coords$distmod.pr * cfactor, 2)
coords <- coords[,c(1, 2, 3, 4, 5, 7, 8, 9, 10, 6)]


write.csv(coords, paste0(inputfile, '_', start, '-', end, '.csv'), row.names = F)

}) # end of system time measurement




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






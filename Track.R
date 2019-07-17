# Lumpfish discrimination learning project
# Adam Brooker
# 24th June 2019


library(EBImage)
library(imager) 

workingdir <- 'G:/Projects/Cleaner fish learning/Data/Test1'
setwd(workingdir)

inputfile <- 'LumpfishBehaviourTank'


# load first image to get boundaries
files <- list.files(path = workingdir, pattern = inputfile, all.files = FALSE, recursive = FALSE)

image <- readImage(files[1]) # load image file
image <- channel(image, 'green') # Use green channel for max contrast against green tank
image <- image*(1/max(image)) # increase contrast to max

image <- gblur(image, sigma = 3) # Gaussian smoothing filter

# Laplacian high-pass filter for detecting edges
fhi <- matrix(1, nrow = 3, ncol = 3)
fhi[2, 2] <- -7.5
image_fhi <- filter2(image, fhi)
display(image_fhi)


#otsu filter
display(image < otsu(image))

# adaptive thresholding
disc <- makeBrush(121, 'disc')
disc <- disc/sum(disc)
offset <- 0.00001
image_bg <- filter2(image, disc)
image_th <- image < image_bg + offset
display(image_th, all = T)

# adaptive thresholding
image_th <- thresh(image, w = 15, h = 15, offset = 0.0001)
display(1-image_th)

nmask <- thresh(image, w = 15, h = 15, offset = 0.0001)
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






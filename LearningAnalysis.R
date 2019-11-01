# Lumpfish learning analysis
# Copyright Adam Brooker 31st July 2019


library(ggplot2)
library(colorRamps)

workingdir <- ifelse(Sys.info()['user'] == 'Laptop', 'G:/Projects/Cleaner fish learning/Data/Test2', '/Users/adambrooker/R Projects/FishTrackR/ShortTest') # change to location of data
setwd(workingdir)

file <- 'Week3Test-complete.csv'

clist <- read.csv(file)


pen.col <- 'black'
pen.size <- 1.4
#plot.col <- rev(heat.colors(2, alpha = 1))
plot.col <- matlab.like(20)  

ggplot(clist, aes(fishx, fishy)) +
  geom_hex(bins = 20, alpha = 0.6) + 
  scale_fill_gradientn(colours=plot.col, space = 'Lab', limits = c(0, 20), na.value = plot.col[length(plot.col)], name = 'No. frames') +
  #annotate('segment', x = locations.lookup['7CSW', 'xmin'], xend = locations.lookup['7CSE', 'xmax'], y = locations.lookup['7CSW', 'ymin'], yend = locations.lookup['7CSE', 'ymin'], colour = pen.col, size = pen.size) + # pen boundary
  #annotate('curve', x = locations.lookup['7WHNW', 'xmin']+1, xend = locations.lookup['7WHNW', 'xmax']-1, y = locations.lookup['7WHNW', 'ymin']+1, yend = locations.lookup['7WHNW', 'ymax']-1, colour = pen.col, size = pen.size, curvature = 1) + # hide boundary
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  scale_x_continuous('x (pixels)', limits = c(0, 1100)) + scale_y_continuous('y (pixels)', limits = c(0,1100))


ggplot(clist, aes(fishx, fishy)) +
  geom_path() + 
  scale_x_continuous(limits = c(0, 1100)) +
  scale_y_continuous(limits = c(0, 1100)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))

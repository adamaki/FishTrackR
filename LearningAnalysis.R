# Lumpfish learning analysis
# Copyright Adam Brooker 31st July 2019


library(ggplot2)
library(colorRamps)

workingdir <- 'G:/Projects/Cleaner fish learning/Data/test3/Tracked' # change to location of data
setwd(workingdir)

file <- 'test_w3_n1_exp_t1.csv'

clist <- read.csv(file)


pen.col <- 'black'
pen.size <- 1.4
#plot.col <- rev(heat.colors(2, alpha = 1))
plot.col <- matlab.like(20)  

ggplot(clist, aes(fish.rx, fish.ry)) +
  geom_hex(bins = 20, alpha = 0.6) + 
  scale_fill_gradientn(colours=plot.col, space = 'Lab', limits = c(0, 20), na.value = plot.col[length(plot.col)], name = 'No. frames') +
  #annotate('segment', x = locations.lookup['7CSW', 'xmin'], xend = locations.lookup['7CSE', 'xmax'], y = locations.lookup['7CSW', 'ymin'], yend = locations.lookup['7CSE', 'ymin'], colour = pen.col, size = pen.size) + # pen boundary
  #annotate('curve', x = locations.lookup['7WHNW', 'xmin']+1, xend = locations.lookup['7WHNW', 'xmax']-1, y = locations.lookup['7WHNW', 'ymin']+1, yend = locations.lookup['7WHNW', 'ymax']-1, colour = pen.col, size = pen.size, curvature = 1) + # hide boundary
  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
  scale_x_continuous('x (pixels)', limits = c(0, 120)) + scale_y_continuous('y (pixels)', limits = c(0,120))


ggplot(clist, aes(fish.rx, fish.ry)) +
  geom_path() + 
  scale_x_continuous(limits = c(0, 120)) +
  scale_y_continuous(limits = c(0, 120)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))

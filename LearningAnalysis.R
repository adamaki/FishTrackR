# Lumpfish learning analysis
# Copyright Adam Brooker 31st July 2019


library(ggplot2)
library(colorRamps)
library(zoo)

workingdir <- '/Users/adambrooker/R Projects/FishTrackR/Data/Control/Week 1/Training/Trial 1' # change to location of data
setwd(workingdir)

file <- 'GOPR0001.csv'

clist <- read.csv(file)
clist$errors <- NULL

# Calculations-----------------------------------------

clist$velocity <- c(NA, sqrt(abs(diff(clist$fish.rx))^2 + abs(diff(clist$fish.ry))^2))
clist$move <- c(NA, ifelse(rollapply(clist$velocity, width = 3, FUN = mean, na.rm = T, align = 'center') > 0.9, T, F), NA)

summary(clist$move)


# Plots-----------------------------------------------

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
  scale_x_continuous('x (pixels)', limits = c(0, 100)) + 
  scale_y_reverse('y (pixels)', limits = c(100, 0))


ggplot(clist, aes(fish.rx, fish.ry)) +
  geom_path() + 
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_reverse(limits = c(100, 0)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))

ggplot(clist) +
  geom_line(aes(frame, distmod.rl)) +
  geom_line(aes(frame, distmod.rr), col = 'red') +
  scale_y_continuous(limits = c(0, 70), name = 'Distance from models (mm)', expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), name = 'Time (s)') +
  theme_minimal()











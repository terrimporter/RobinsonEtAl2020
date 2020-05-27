#Chloe Robinson, Feb 28 2020

library(ggmap) # map
library("ggsn") # map scale bar
library(scales) # comma
library(ggplot2) #ggmap
library(ggrepel) #geom_label_repel

# read in sites and coord
a <- read.table("Sites.csv", header = TRUE, sep = ",",stringsAsFactors = FALSE)
head(a)

a$label <- paste (a$Site, a$Region, sep=" ")

# get map, the higher the zoom, the more detail on the map
map <- get_stamenmap(bbox = c(left = -80.70, bottom = 43.20, 
                              right = -80.40, top = 43.60), 
                     zoom = 11)

p <- ggmap(map) +
  geom_point(data=a, aes(x=Long, y=Lat), size=1) +
  geom_label_repel(data=a, aes(x=Long, y=Lat, label = a$Region), box.padding = unit (2,"lines"))+
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none") +
  scalebar(x.min = -80.70, x.max = -80.40, y.min = 43.20, y.max = 43.60, location = "bottomright", 
           dist = 5, dist_unit = "km", transform = TRUE, model="WGS84",
           st.bottom = FALSE, st.color = "black", st.size = 3, st.dist = 0.025)

p
ggsave("Map.pdf", p, width=8, height=8)

###############################################################################


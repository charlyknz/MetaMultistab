source(here("/Users/charlottekunze/Desktop/phD/Meta_Multistab/MetaMultistab/code/02_MetaMultistab_Analysis.R"))

library(maps)

# change longitude in numeric
metadata$long <- as.numeric(metadata$long)
str(metadata)


# create data for world coordinates using  
# map_data() function 
world_coordinates <- map_data("world") 

# create world map using ggplot() function 
ggplot() + 
  
# geom_map() function takes world coordinates as input to plot world map 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(long, lat, map_id = region), 
    color = "#293133", fill = "white", linewidth = 0.3)+
  
# own data
geom_point( data = metadata, aes(long, lat, color = system), size=2.5) + 
  labs(x = 'Longitude', y = 'Latitude', color = 'System')+

#theme
theme_dark(base_size = 15)+
  theme(legend.position="bottom")

ggsave(plot = last_plot(), file = here('output/Extended_Worldmap.png'), width = 8, height = 6)

#This is code for mapping the different stock regions in the database, and mapping how many stocks are in each region

#load packages
library(ggplot2)
theme_set(theme_bw())
library(sf)
library(writexl)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(rnaturalearthhires)
library(rnaturalearth)
library(rnaturalearthdata)

loc <- "D:/Github/ICM/"

# Bring in the shapefiles
#read in shapefile
nafo.sf <- read_sf(paste0(loc,"Data/Map/NAFO/NAFO_Divisions_SHP/NAFO_Divisions_2021_poly_clipped.shp"))
ices.sf <- read_sf(paste0(loc,"Data/Map/ICES/ICES_Areas_20160601_cut_dense_3857.shp"))
alas.sf <- read_sf(paste0(loc,"Data/Map/Alaska/PVG_Registration_Areas_2012_Present_GCS_WGS1984.shp"))
#alas.sf <- st_make_valid(alas.sf) # Messy shapefile so need to clean it up.
# Now I want to combine these into 'big' areas starting with NAFO I want to have NAFO 0,1, and most of 2 as one polygon.
# Then NAFO 2J, 3, and 4 will be one polygon
# and NAFO 5 and 6 will be one polygon.
nafo.0.2 <- st_union(nafo.sf[nafo.sf$SubArea %in% c(0,1) | nafo.sf$Division %in% c("2G","2H"),])
nafo.2.4 <- st_union(nafo.sf[nafo.sf$SubArea %in% c(3,4) | nafo.sf$Division %in% c("2J"),])
nafo.5.6 <- st_union(nafo.sf[nafo.sf$SubArea %in% c(5) | nafo.sf$Division %in% c("6A","6B","6C"),])
nafo.6.other <- st_union(nafo.sf[nafo.sf$Division %in% c("6D","6E","6F","6G","6H"),])
# Now give these some attributes
nafo.0.2 <- st_sf(id= "0-2H",name = "NAFO 0-2H",geometry = nafo.0.2)
nafo.2.4 <- st_sf(id= "2J-4",name = "NAFO 2J-4",geometry = nafo.2.4)
nafo.5.6 <- st_sf(id= "5-6C",name = "NAFO 5-6C",geometry = nafo.5.6)
nafo.6.other <- st_sf(id= "6D-H",name = "NAFO 6D-H",geometry = nafo.6.other)
# Make a new NAFO object
new.nafo.sf <- rbind(nafo.0.2,nafo.2.4,nafo.5.6,nafo.6.other)
new.nafo.sf$num <- c(0,7,15,0)
new.nafo.sf$num.cat <- c(0,"5-9","10-15",0)
new.nafo.sf <- new.nafo.sf[new.nafo.sf$num > 0,]
new.nafo.sf$num.cat <- factor(new.nafo.sf$num.cat,levels = c("1-4","5-9","10-15"))

# Now combine into 'big areas' in Alaska, just 2, the gulf of Alaska and the Berring Sea/Aleutian Islands
alas.goa <- st_union(st_make_valid(alas.sf[alas.sf$REGISTRA_1 %in% c("K","M","L","A","E","H"),]))
alas.bs.ai <- alas.sf$geometry[alas.sf$REGISTRA_1 %in% c("O"),]
# Give some attributes...
alas.goa <- st_sf(id= "GOA",name = "Gulf of Alaska",geometry = alas.goa)
alas.bs.ai <- st_sf(id= "BS_AI",name = "Bering Sea and Aleutian Islands",geometry = alas.bs.ai)
# Then make a new NAFO object and put in the numbers
new.alas.sf <- rbind(alas.goa, alas.bs.ai)
new.alas.sf$num <- c(2,9)
new.alas.sf$num.cat <- c("1-4","5-9")
new.alas.sf$num.cat <- factor(new.alas.sf$num.cat,levels = c("1-4","5-9","10-15"))
# Finally let's get ICES areas sorted, this will be a most complicated one by miles...
# First split up into areas, based on ICES Ecoregions Baltic, North Sea, Celtic Seas, Bay of Biscay, NE Arctic, Iceland, Faroes I think is enough
ices.baltic <- st_union(ices.sf[ices.sf$SubArea==3 & ices.sf$Division != "a",])
ices.ns <- st_union(ices.sf[ices.sf$SubArea==4 | 
                            (ices.sf$SubArea == 7 & ices.sf$Division %in% c("d","e")) | 
                            (ices.sf$SubArea == 3 & ices.sf$Division %in% c("a")),])
# This is missing part of the Norwegian Sea at the moment....
ices.ne.arctic <- st_union(ices.sf[ices.sf$SubArea %in% c(1,2),])
# Iceland
ices.ice <- st_union(ices.sf[(ices.sf$SubArea == 5 & ices.sf$Division %in% c("a")) |
                              (ices.sf$SubArea == 12 & ices.sf$Division %in% c("a")) & !ices.sf$SubDivisio %in% c(1,3),]) 
# Faroes
ices.far <- st_union(ices.sf[(ices.sf$SubArea == 5 & ices.sf$Division %in% c("b")& ices.sf$SubDivisio %in% c(1,2)),]) 
# Celtic Seas
ices.cel <- ices.sf[ices.sf$SubArea %in% c(6,7),]
# Drop 6.b.2, 7.c.1 and 7.k.1, and 7.j.1
ices.cel <- ices.cel[-c(13:16),]
ices.cel <- st_union(ices.cel)
# Bay of Biscay and that's where we care about...
ices.bob <- ices.sf[ices.sf$SubArea %in% c(8,9),]
# Get rid of 8d18e1 and 9b1
ices.bob <- st_union(ices.bob[-c(8:10),])
# OK, not give them attributes
ices.baltic <- st_sf(id= "Baltic",name = "Baltic Eco-region",geometry = ices.baltic)
ices.ns <- st_sf(id= "Baltic",name = "Baltic Eco-region",geometry = ices.ns)
ices.ne.arctic <- st_sf(id= "NE Arctic",name = "Norwegian and Barents Sea Eco-regions",geometry = ices.ne.arctic)
ices.ice <- st_sf(id= "Iceland",name = "Icelandic Eco-region",geometry = ices.ice)
ices.far <- st_sf(id= "Faroes",name = "Faroes Eco-region",geometry = ices.far)
ices.cel <- st_sf(id= "Celtic",name = "Celtic Seas Eco-region",geometry = ices.cel)
ices.bob <- st_sf(id= "Bob and IC",name = "Bay of Biscay & Iberian Coast Eco-region",geometry = ices.bob)
# Now combine and count stocks
new.ices.sf <- rbind(ices.baltic, ices.ns,ices.ne.arctic,ices.ice,ices.far,ices.cel,ices.bob)
new.ices.sf$num <- c(5,14,5,3,4,12,1)
new.ices.sf$num.cat <- c("5-9","10-15","5-9","1-4","1-4","10-15","1-4")
new.ices.sf$num.cat <- factor(new.ices.sf$num.cat,levels = c("1-4","5-9","10-15"))

# Merge all objects into one using the Robinson projection
rob.crs <- "+proj=robin"
new.ices.sf <- st_transform(new.ices.sf,crs = rob.crs)
new.nafo.sf <- st_transform(new.nafo.sf,crs = rob.crs)
new.alas.sf <- st_transform(new.alas.sf,crs = rob.crs)
# Merge into a super sf
all.sf <- rbind(new.ices.sf,new.nafo.sf,new.alas.sf)

# Colors for the map...
cols <- c("#666666", "#E6AB02", "#7570B3")



################################## Map ###################################################################################
# Download and filter world data
world <- ne_countries(scale = "medium", returnclass = "sf")
# Reproject to Robinson
world.rob <- st_transform(world, crs = rob.crs)

# Plot
p <- ggplot(data = world.rob) +
      geom_sf(fill = "#e5d8bd", color = "black", size = 0.2) +
      geom_sf(data=all.sf, aes(fill = num.cat),color="white", size=10) + 
      theme_minimal() +
      scale_fill_manual(name = "Number of stocks",values=cols)+
      coord_sf(xlim = c(-17000000, 5000000), ylim = c(0.000000, 8500000), crs = rob.crs) 
  
save_plot(paste0(loc,("Figures/Map_of_stocks_by_region.png")),p,base_height = 8)
save_plot(paste0(loc,("Figures/figure_1.tiff")),p,base_height = 8)
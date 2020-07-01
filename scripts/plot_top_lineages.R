# Thanh Le Viet - QIB - 2020
rm(list=ls())
library(pacman)

p_load(tidyverse,
       googlesheets4,
       ggpmisc,
       rgeos,
       maps,
       ggpubr,
       leaflet,
       data.table,
       glue,
       viridis)

#Master metatdata table on googlesheet
gs_id <- "secret"
sheet_url <- glue("https://docs.google.com/spreadsheets/d/{gs_id}/edit#gid=0")
master_table <- read_sheet(sheet_url)

#private metadata downloaded from CLIMB
cog_meta <- fread("cog_global_2020-06-12_metadata_private.csv")

cog_meta_norw <- grepl("NORW", cog_meta$sequence_name) %>% 
  cog_meta[.,] %>% 
  mutate(central_sample_id = str_extract(sequence_name, "NORW-E[A-Z0-9]{4}")) %>% 
  select(central_sample_id, cog_uk_lineage = uk_lineage)

selected_master_table <- master_table %>% 
  select(central_sample_id, adm2_private) 


norwich <- master_table %>% 
  select(central_sample_id, adm2_private, lineage = uk_lineage) %>% 
  left_join(., cog_meta_norw) %>% 
  select(central_sample_id, adm2_private, lineage = cog_uk_lineage)


lineage_table_origin <- norwich %>%
  filter(!is.na(lineage)) %>% 
  group_by(lineage) %>% 
  summarise(counts = n()) %>% 
  arrange(desc(counts))

unique_admin2_postcode <- unique(norwich$adm2_private) %>%
.[!is.na(.)] #Remove NA values

#Original data downloaded from https://www.opendoorlogistics.com/wp-content/uploads/Data/UK-postcode-boundaries-Jan-2015.zip
england <- readRDS("data/england_postcode.rds")

england@data$id = getSpPPolygonsIDSlots(england)

# Find postcode ids that have samples
found <- england@data$name %in% unique_admin2_postcode

#Sort postcode areas by their centroid, long: top-> bottem; lat:left->right
selected_postcodes <- england[found,]
selected_postcodes_map <-  map(selected_postcodes, plot=FALSE, fill = TRUE)
selected_postcodes_map_centroids <- maps:::apply.polygon(selected_postcodes_map, maps:::centroid.polygon)
selected_postcodes_map_centroids_matrix <- Reduce(rbind, selected_postcodes_map_centroids)
dimnames(selected_postcodes_map_centroids_matrix) <- list(gsub("[^,]*,", "", names(selected_postcodes_map_centroids)),
                                            c("long", "lat"))
selected_postcodes_map_centroids_matrix_sorted <- selected_postcodes_map_centroids_matrix %>% 
  as.data.frame() %>% 
  mutate(postcode = row.names(.)) %>% 
  mutate(postcode= gsub(":[0-9]","", postcode)) %>% 
  arrange(long, desc(lat))

ordered_postcodes_0 <- unique(selected_postcodes_map_centroids_matrix_sorted$postcode)

#Count lineages by postcode
np_selected_postcodes <- data_frame(adm2_private = selected_postcodes$name) %>%
  left_join(., norwich) %>% 
  group_by(adm2_private, lineage) %>% 
  summarise(counts = n())


ordered_postcodes <- selected_postcodes$name
ordered_lineage <- str_sort(norwich$lineage, numeric = TRUE) %>% unique()

# Categorise number of samples by postcode
cat_samples <- cut(np_selected_postcodes$counts, c(0,1,5,10,20,30,100))
levels_cat_samples <- levels(cat_samples)

np_selected_postcodes$samples <- cat_samples
np_selected_postcodes$lineage <- gsub("\\/.*","",np_selected_postcodes$lineage)


manual_color <- viridis::viridis(6)

true_heatmap <- ggplot(np_selected_postcodes) +
  geom_tile(aes(factor(adm2_private, levels = ordered_postcodes_0), factor(lineage, levels = ordered_lineage), fill = samples), color = "white") +
  scale_fill_manual(values = manual_color) +
  coord_fixed(ratio = 1) +
  theme_linedraw() +
  # ggtitle("Circulation of UK SARS-Cov-2 viral lineages in Norfolk and Suffolk (Update 2020-05-15)") +
  ylab("Postcode") +
  xlab("Lineage") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text=element_text(size=2),
        plot.margin = margin(0.5, 0, 0, 0,"cm"))

# Map
# base_map_data <- readRDS("aoi_fortified_sub.rds")


# Create a base map to represent the postcodes with missing value in grey
# base_map_data <- readRDS("data/aoi_fortified_latest.rds") # This is for NORW

#Create a bounding box for the set of postcodes having reported lineages
bbox_has_data <- rgeos::bbox2SP(bbox = bbox(selected_postcodes), proj4string = CRS(proj4string(england)))

# Create a bounding box based on coordinates of the area of interest
# bbox_has_data <- rgeos::bbox2SP(n = n, s = s, w = w, e = e, proj4string = CRS(proj4string(england)))

#Clip area of interest
aoi <- gIntersection(england,bbox_has_data, byid = T)
#Get polygon I
aoi_id <- getSpPPolygonsIDSlots(aoi)
#Create polygon dataframe
aoi_df <- data.frame(id = aoi_id, row.names = aoi_id)
#Create SpatialPolygonDataFrame
aoi_spdf <- SpatialPolygonsDataFrame(aoi, aoi_df)
aoi_spdf@data <- aoi_spdf@data %>% 
  mutate(id = gsub(" 1", "", id))


#Tidy basemap data for ggplot
base_map_data <- broom::tidy(aoi_spdf, region = "id") %>% arrange()

top_lineages <- lineage_table_origin %>% 
  filter(!is.na(lineage), lineage != "",counts >=10 ) %>% 
  top_n(n = 10) %>%
  pull(lineage)

top_plots <- list()

for (i in seq_along(top_lineages)) {
  lineage_vector <- top_lineages[i]
  sub_lineage <- dplyr::filter(np_selected_postcodes, lineage == lineage_vector)
  #Subset postcodes having data
  sub_lineage_data <- england[match(sub_lineage$adm2_private,england@data$name),]
  sub_lineage_data@data <- sub_lineage_data@data %>% 
    mutate(counts = sub_lineage$samples)
  
  lineage_levels <- unique(sub_lineage$samples) %>% as.character()
  lineage_manual_color_index <- match(lineage_levels, levels_cat_samples) %>% sort()
  lineage_manual_color <- manual_color[lineage_manual_color_index]
  #lazy variable
  has_data <- sub_lineage_data
  
  #Tidy data for ggplot
    tidy_has_data <- fortify(has_data, region="name") %>% 
    arrange()
  #Mean label for each postcode
  sub_lineage_data_label <- tidy_has_data %>%
    group_by(id) %>%
    summarise(long = mean(long), lat = mean(lat))
  # print(i)
  # Create map
  top_plots[[i]] <- ggplot(base_map_data) +
    geom_map(aes(long, lat, map_id = order), map = base_map_data, fill = "grey70", color = "#000000", size = 0.2, show.legend = FALSE) +
    geom_map(data = has_data@data, aes(fill = counts, map_id = name), map = tidy_has_data, color = "white", size = 0.2, show.legend = TRUE) +
    geom_text(aes(label = id, x = long, y = lat), data = sub_lineage_data_label, size = 2, color = "tomato2") +
    scale_fill_manual(values = lineage_manual_color) +
    ggtitle(glue("{lineage_vector}")) +
    theme(legend.position="right") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}


# gg table of lineage distribution
tbody.style = tbody_style(color = "black",
                          hjust=1, x=0.8,
                          size = 9)
lineage_table <- lineage_table_origin %>%
  select(`top 10` = lineage, count = counts)
g.table <- ggtexttable(lineage_table %>% top_n(10),
                       rows = NULL,
                       theme = ttheme(
                         tbody.style = tbody.style,
                       ))
  
top_panel <- ggarrange(true_heatmap, g.table,
                       nrow = 1, ncol = 2,
                       widths = c(4,1),
                       labels = "Circulation of UK SARS-Cov-2 viral lineages in Norfolk and Suffolk (Update 2020-06-17)",
                       vjust = 1.5,
                       hjust = -0.7,
                       font.label = list(size = 8)
                       )
#Could be improved with a for-loop
bottom_panel_1 <- ggarrange(plotlist = top_plots[1:4], 
                          ncol = 4, nrow = 1, labels = "", 
                          label.x = 2.5, label.y = 1, align = "h")

bottom_panel_2 <- ggarrange(plotlist = top_plots[5:8], 
                          ncol = 4, nrow = 1, labels = "", 
                          label.x = 2.5, label.y = 1, align = "h")

bottom_panel_3 <- ggarrange(plotlist = top_plots[9:12], 
                            ncol = 4, nrow = 1, labels = "", 
                            label.x = 2.5, label.y = 1, align = "h")

bottom_panel_4 <- ggarrange(plotlist = top_plots[13:16],
ncol = 4, nrow = 1, labels = "",
label.x = 2.5, label.y = 1, align = "h")

bottom_panel_5 <- ggarrange(plotlist = top_plots[17:20],
                            ncol = 4, nrow = 1, labels = "", 
                            label.x = 2.5, label.y = 1, align = "h")

# Create combined graphic object for saving to a graphic device
gobj <- ggarrange(
  top_panel,
  bottom_panel_1,
  bottom_panel_2,
  bottom_panel_3,
  bottom_panel_4,
  nrow = 5,
  ncol = 1
)

ggsave("top_cog_uk_lineages.png", gobj, width = 27, height = 38, units = "cm")

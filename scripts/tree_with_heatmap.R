library(ape)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggnewscale)
library(stringr)
library(dplyr)


metadata <- read.csv("../data/dtol/dtol_all_samples.taxonomic_classification.csv")

somatic_signatures <- read.csv("../data/dtol/somatic_mutational_signature_attributions.x0_excluded.rtol_filtered.csv",check.names = F)
colnames(somatic_signatures) = paste0('sToL',colnames(somatic_signatures) )
rs = rowSums(somatic_signatures[ , -1], na.rm = TRUE) # calculate row sum
somatic_signatures[ , -1] <- somatic_signatures[ , -1] / rs # normalise mutational signature attribution 
somatic_signatures$sToL8_2 <- somatic_signatures$sToL8 + somatic_signatures$sToL2 # sum two CpG signatures = sToL8 + sToL2 
somatic_signatures[ , -1][somatic_signatures[ , -1] < 0.035] <- 0 # Set values < 0.035 to 0 
somatic_signatures[,1] <- ifelse(grepl("\\.", somatic_signatures[,1]),
                                 sub(".*\\.", "", somatic_signatures[,1]),
                                 somatic_signatures[,1])
somatic_signatures$label=NA
somatic_signatures$label=sapply(somatic_signatures[,1], function(x){
  Species = metadata$Species[metadata$Sample==x]
  if (length(metadata$Species[metadata$Species==Species]) == 1){
    label = Species
  }
  else{
    label = paste0(Species," (", x, ")")
  } 
  label
})
somatic_signatures$label[somatic_signatures$label=="Hemaris fuciformis"]="Hemaris fuciformis (iHemFuc2)"



germline_signatures <- read.csv("../data/dtol/germline_mutational_signature_attributions.x0_excluded.csv",check.names = F)
colnames(germline_signatures) = paste0('gToL',colnames(germline_signatures))
germline_signatures[,1] <- ifelse(grepl("\\.", germline_signatures[,1]),
                                 sub(".*\\.", "", germline_signatures[,1]),
                                 germline_signatures[,1])
germline_signatures$label=NA
germline_signatures$label=sapply(germline_signatures[,1], function(x){
  Species = metadata$Species[metadata$Sample==x]
  if (length(metadata$Species[metadata$Species==Species]) == 1){
    label = Species
  }
  else{
    label = paste0(Species," (", x, ")")
  } 
  label
})

germline_signatures$label[germline_signatures$label=="Hemaris fuciformis"]="Hemaris fuciformis (iHemFuc2)"



txt  <- readLines("../data/dtol/dtol_all_samples.taxonomic_tree.nwk")
txt <- gsub('_', '^', txt, fixed = TRUE)
txt2 <- gsub(" +", "_", txt)   # turn spaces inside labels into underscores
tr   <- read.tree(text = txt2)
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")
tr$tip.label <- gsub("\\^(.*?)\\^", "(\\1)", tr$tip.label)
org_labels =  tr$tip.label
tr_full <- tr   # full 764-tip tree; each figure below prunes it to the samples it has data for
# 4 artefact-dominated samples (cumulative artefact attribution 0.95-0.97) are not shown in the somatic trees
artefact_samples <- somatic_signatures$label[somatic_signatures[,1] %in% c("gfFlaVelt1", "ihDrePlat2", "ilYpoCagn5", "ilYpoPade1")]

get_expr_labels <- function(labels){
  esc <- function(s) gsub("'", "\\\\'", s, perl = TRUE)  # escape single quotes
  
  has_paren = grepl("\\(", labels)
  
  # get species name
  species   <- sub("\\s*\\(.*$", "", labels) 
  
  # get sample name
  paren     <- sub("^[^\\(]*", "", labels)
  
  species_e <- esc(species)
  paren_e   <- esc(paren)
  
  ifelse(
    has_paren,
    paste0("italic('", species_e, "')~'", paren_e, "'"),  # italic species + plain "(...)"
    paste0("italic('", species_e, "')")                   # all italic when no "(...)"
  )
}

# tr$tip.label = get_expr_labels(org_labels)
# ------------------Fig3 taxa enriched signatures--------------
# keep only samples with somatic attributions, minus the 4 artefact-dominated samples (678 tips); tips without data would otherwise be drawn white (= zero)
tr <- keep.tip(tr_full, setdiff(intersect(tr_full$tip.label, somatic_signatures$label), artefact_samples))
org_labels <- tr$tip.label
pal <- c(
  Coleoptera = "#A3843F",
  Chordata = "#765FA6",
  Viridiplantae = "#408145",
  Vespidae = "#EF9581",
  Tenthredinidae = "#F0D3CC"
)


anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno$Phylum <- sapply(anno$Species,function(x){
  unique(metadata$Phylum[metadata$Species==x])
})
anno$Class <- sapply(anno$Species,function(x){
  unique(metadata$Class[metadata$Species==x])
})
anno$Order <- sapply(anno$Species,function(x){
  unique(metadata$Order[metadata$Species==x])
})
anno$Family <- sapply(anno$Species,function(x){
  unique(metadata$Family[metadata$Species==x])
})

anno <- anno %>%
  mutate(ColorGroup = case_when(
    Kingdom %in% names(pal) ~ Kingdom,
    Phylum %in% names(pal) ~ Phylum,
    Class  %in% names(pal) ~ Class,
    Order  %in% names(pal) ~ Order,
    Family  %in% names(pal) ~ Family,
    TRUE ~ NA_character_   # << no group = no colour
  )) %>% filter(!is.na(ColorGroup)) 


grp_list <- split(anno$label, anno$ColorGroup)

hilight_df <- lapply(names(grp_list), function(g) {
  tips <- grp_list[[g]]
  idx  <- which(tr$tip.label %in% tips)
  if (length(idx) >= 2) {
    node_id <- ape::getMRCA(tr, idx)
    if (!is.na(node_id)) data.frame(node = node_id, ColorGroup = g)
  }
}) %>% bind_rows()

p <- ggtree(tr, layout = "circular", size = 0.25,color = "#090954") %<+% anno +
  geom_hilight(data = hilight_df,
               aes(node = node, fill = ColorGroup),
               alpha = 1) +
  # geom_tiplab(aes(label = label),
  #             offset = 0.3,
  #             fontface = "italic",
  #             family = "Helvetica",
  #             size = 1) +
  scale_fill_manual(values = pal, name = "Taxonomic rank", 
                    guide = guide_legend(override.aes = list(alpha = 1))) #+
  # theme(legend.position = "right",
  #       text = element_text(family = "Helvetica"))

p

ggsave("../figs/ToL_Fig3_noheatmap.pdf", plot = p,
       width = 7, height = 7, units = "in") 



# published sToL8+sToL2, sToL4, sToL13, sToL24 
hm <- somatic_signatures %>%
  select(label, sToL8_2, sToL4, sToL13, sToL24) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% org_labels) %>%   
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
# custom palettes for each signature
pal_1 <- c("#F7FCF5","#E1F3DC","#BCE4B5","#8ED08B","#56B567","#2C944C","#05712F","#00441B")  # greens 
pal_3 <- c("#FFF5EB","#FEE3C8","#FDC692","#FDA057","#F67824","#E05206","#AD3803","#7F2704")  # oranges
pal_4 <- c("#FCFBFD","#ECEBF4","#D1D2E7","#AFAED4","#8D89C0","#705EAA","#572C92","#3F007D")  # purples

# layout parameters
base_offset <- 0.05
band_width  <- 0.05
gap         <- 0.3

grid_col  <- "grey85" 


p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"sToL8_2", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_2, name = "sToL8 + sToL2",na.value = "white",
                       limits  = c(0, 1.0),
                       breaks  = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       labels  = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 1))
p1

p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, hm[,"sToL4", drop = F],
               offset = base_offset + band_width + gap,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_1, name = "sToL4",na.value = "white",
                       limits  = c(0, 0.9),
                       breaks  = c(0, 0.2, 0.4, 0.6, 0.8),
                       labels  = c("0", "0.2", "0.4", "0.6", "0.8"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 2))
p2

p3 <- p2 + ggnewscale::new_scale_fill()
p3 <- gheatmap(p3, hm["sToL13", drop = FALSE],
               offset = base_offset + 2*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_4, name = "sToL13",na.value = "white",
                       limits  = c(0, 0.5),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 3))
p3

p4 <- p3 + ggnewscale::new_scale_fill()
p4 <- gheatmap(p4, hm["sToL24", drop = FALSE],
               offset = base_offset + 3*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_3, name = "sToL24",na.value = "white",
                       limits  = c(0, 0.8),
                       breaks  = c(0, 0.2, 0.4, 0.6, 0.8),
                       labels  = c("0", "0.2", "0.4", "0.6", "0.8"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 4))
p4
p4 <- p4 + theme(legend.text = element_text(size = 6)) # as plot_figure_3.R
ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig3_heatmap_scale_right.pdf", plot = p4,
       width = 7, height = 7, units = "in") 

p4 <- p4+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 6))
p4

ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig3_heatmap_scale_bottom.pdf", plot = p4,
       width = 7, height = 7, units = "in") 


# ------------------Fig4--------------
# keep only samples with somatic attributions, minus the 4 X6-dominated samples (678 tips)
tr <- keep.tip(tr_full, setdiff(intersect(tr_full$tip.label, somatic_signatures$label), x6_samples))
org_labels <- tr$tip.label
pal <- c(
  Actiniaria = "#32A2DB",
  Actinopteri = "#35529A",
  Echinodermata = "#9BD0F2",
  Mollusca = "#70C4BE"
)


anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno$Phylum <- sapply(anno$Species,function(x){
  unique(metadata$Phylum[metadata$Species==x])
})
anno$Class <- sapply(anno$Species,function(x){
  unique(metadata$Class[metadata$Species==x])
})
anno$Order <- sapply(anno$Species,function(x){
  unique(metadata$Order[metadata$Species==x])
})
anno$Family <- sapply(anno$Species,function(x){
  unique(metadata$Family[metadata$Species==x])
})

anno <- anno %>%
  mutate(ColorGroup = case_when(
    Kingdom %in% names(pal) ~ Kingdom,
    Phylum %in% names(pal) ~ Phylum,
    Class  %in% names(pal) ~ Class,
    Order  %in% names(pal) ~ Order,
    Family  %in% names(pal) ~ Family,
    TRUE ~ NA_character_   # << no group = no colour
  )) %>% filter(!is.na(ColorGroup)) 


grp_list <- split(anno$label, anno$ColorGroup)

hilight_df <- lapply(names(grp_list), function(g) {
  tips <- grp_list[[g]]
  idx  <- which(tr$tip.label %in% tips)
  if (length(idx) >= 2) {
    node_id <- ape::getMRCA(tr, idx)
    if (!is.na(node_id)) data.frame(node = node_id, ColorGroup = g)
  }
}) %>% bind_rows()

p <- ggtree(tr, layout = "circular", size = 0.25,color = "#090954") %<+% anno +
  geom_hilight(data = hilight_df,
               aes(node = node, fill = ColorGroup),
               alpha = 1) +
  # geom_tiplab(aes(label = label),
  #             offset = 0.3,
  #             fontface = "italic",
  #             family = "Helvetica",
  #             size = 1) +
  scale_fill_manual(values = pal, name = "Taxonomic rank", 
                    guide = guide_legend(override.aes = list(alpha = 1))) #+
# theme(legend.position = "right",
#       text = element_text(family = "Helvetica"))

p

ggsave("../figs/ToL_Fig4_noheatmap.pdf", plot = p,
       width = 7, height = 7, units = "in") 



# published sToL9/13/18 
hm <- somatic_signatures %>%
  select(label, sToL9, sToL15, sToL21) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% org_labels) %>%   
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
# custom palettes for each signature
pal_1 <- c("#F7FCF5","#E1F3DC","#BCE4B5","#8ED08B","#56B567","#2C944C","#05712F","#00441B")  # greens (manuscript)
pal_2 <- c("#F7FBFF","#DBE9F6","#BAD6EB","#89BEDC","#539ECD","#2B7BBA","#0B559F","#08306B")  # blues
pal_3 <- c("#FFF5EB","#FEE3C8","#FDC692","#FDA057","#F67824","#E05206","#AD3803","#7F2704")  # oranges
pal_4 <- c("#FCFBFD","#ECEBF4","#D1D2E7","#AFAED4","#8D89C0","#705EAA","#572C92","#3F007D")  # purples


# layout parameters
base_offset <- 0.05
band_width  <- 0.05
gap         <- 0.3

grid_col  <- "grey85" 


p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"sToL9", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_2, name = "sToL9",na.value = "white",
                       limits  = c(0, 0.5),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 1))
p1

p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, hm[,"sToL15", drop = F],
               offset = base_offset + band_width + gap,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_3, name = "sToL15",na.value = "white",
                       limits  = c(0, 0.5),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 2))
p2

p3 <- p2 + ggnewscale::new_scale_fill()
p3 <- gheatmap(p3, hm["sToL21", drop = FALSE],
               offset = base_offset + 2*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_4, name = "sToL21",na.value = "white",
                       limits  = c(0, 0.5),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 4))
p3
p3 <- p3 + theme(legend.text = element_text(size = 6)) # as plot_figure_4.R
ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig4_heatmap_scale_right.pdf", plot = p3,
       width = 7, height = 7, units = "in") 

p3 <- p3+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 6))
p3

ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig4_heatmap_scale_bottom.pdf", plot = p3,
       width = 7, height = 7, units = "in") 

# ------------------Fig5--------------
# keep only samples with germline attributions (685 tips)
tr <- keep.tip(tr_full, intersect(tr_full$tip.label, germline_signatures$label))
org_labels <- tr$tip.label
pal <- c(
  Coleoptera = "#A3843F",
  Chordata = "#765FA6",
  Viridiplantae = "#408145",
  Fungi = "#925C35"
)


anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno$Phylum <- sapply(anno$Species,function(x){
  unique(metadata$Phylum[metadata$Species==x])
})
anno$Class <- sapply(anno$Species,function(x){
  unique(metadata$Class[metadata$Species==x])
})
anno$Order <- sapply(anno$Species,function(x){
  unique(metadata$Order[metadata$Species==x])
})
anno$Family <- sapply(anno$Species,function(x){
  unique(metadata$Family[metadata$Species==x])
})

anno <- anno %>%
  mutate(ColorGroup = case_when(
    Kingdom %in% names(pal) ~ Kingdom,
    Phylum %in% names(pal) ~ Phylum,
    Class  %in% names(pal) ~ Class,
    Order  %in% names(pal) ~ Order,
    Family  %in% names(pal) ~ Family,
    TRUE ~ NA_character_   # << no group = no colour
  )) %>% filter(!is.na(ColorGroup)) 


grp_list <- split(anno$label, anno$ColorGroup)

hilight_df <- lapply(names(grp_list), function(g) {
  tips <- grp_list[[g]]
  idx  <- which(tr$tip.label %in% tips)
  if (length(idx) >= 2) {
    node_id <- ape::getMRCA(tr, idx)
    if (!is.na(node_id)) data.frame(node = node_id, ColorGroup = g)
  }
}) %>% bind_rows()

p <- ggtree(tr, layout = "circular", size = 0.25,color = "#090954") %<+% anno +
  geom_hilight(data = hilight_df,
               aes(node = node, fill = ColorGroup),
               alpha = 1) +
  # geom_tiplab(aes(label = label),
  #             offset = 0.3,
  #             fontface = "italic",
  #             family = "Helvetica",
  #             size = 1) +
  scale_fill_manual(values = pal, name = "Taxonomic rank", 
                    guide = guide_legend(override.aes = list(alpha = 1))) #+
# theme(legend.position = "right",
#       text = element_text(family = "Helvetica"))

p

ggsave("../figs/ToL_Fig5_noheatmap.pdf", plot = p,
       width = 7, height = 7, units = "in") 



hm <- germline_signatures %>%
  select(label, gToL1, gToL3, gToL4, gToL6) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% org_labels) %>%   
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
# custom palettes for each signature
pal_1 <- c("#F7FCF5","#E1F3DC","#BCE4B5","#8ED08B","#56B567","#2C944C","#05712F","#00441B")  # greens (manuscript)
pal_2 <- c("#F7FBFF","#DBE9F6","#BAD6EB","#89BEDC","#539ECD","#2B7BBA","#0B559F","#08306B")  # blues
pal_3 <- c("#FFF5EB","#FEE3C8","#FDC692","#FDA057","#F67824","#E05206","#AD3803","#7F2704")  # oranges
pal_4 <- c("#FCFBFD","#ECEBF4","#D1D2E7","#AFAED4","#8D89C0","#705EAA","#572C92","#3F007D")  # purples

# layout parameters
base_offset <- 0.05
band_width  <- 0.05
gap         <- 0.3

grid_col  <- "grey85" 


p <- p + guides(fill = "none", colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"gToL1", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_3, name = "gToL1",na.value = "white",
                       limits  = c(0, 1.0),
                       breaks  = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       labels  = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 1))
p1

p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, hm[,"gToL3", drop = F],
               offset = base_offset + band_width + gap,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_2, name = "gToL3",na.value = "white",
                       limits  = c(0, 0.9),
                       breaks  = c(0, 0.2, 0.4, 0.6, 0.8, 0.9),
                       labels  = c("0", "0.2", "0.4", "0.6", "0.8", "0.9"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 2))
p2

p3 <- p2 + ggnewscale::new_scale_fill()
p3 <- gheatmap(p3, hm["gToL4", drop = FALSE],
               offset = base_offset + 2*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_1, name = "gToL4",na.value = "white",
                       limits  = c(0, 0.7),
                       breaks  = c(0, 0.2, 0.4, 0.6, 0.7),
                       labels  = c("0", "0.2", "0.4", "0.6", "0.7"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 3))
p3

p4 <- p3 + ggnewscale::new_scale_fill()
p4 <- gheatmap(p4, hm["gToL6", drop = FALSE],
               offset = base_offset + 3*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 2,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_4, name = "gToL6",na.value = "white",
                       limits  = c(0, 0.4),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 4))
p4
p4 <- p4 + theme(legend.text = element_text(size = 6)) # as plot_figure_5.R
ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig5_heatmap_scale_right.pdf", plot = p4,
       width = 7, height = 7, units = "in") 

p4 <- p4+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(size = 6))
p4

ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig5_heatmap_scale_bottom.pdf", plot = p4,
       width = 7, height = 7, units = "in") 


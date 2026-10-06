library(ape)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggnewscale)
library(stringr)
library(dplyr)


# Manuscript typography: keep all plot text in Helvetica at 5-7 pt.
theme_set(theme(
  text = element_text(family = "Helvetica", size = 6),
  axis.text = element_text(family = "Helvetica", size = 5),
  axis.title = element_text(family = "Helvetica", size = 6),
  plot.title = element_text(family = "Helvetica", size = 7),
  legend.text = element_text(family = "Helvetica", size = 6),
  legend.title = element_text(family = "Helvetica", size = 7)
))



metadata <- read.csv("../data/dtol/dtol_all_samples.taxonomic_classification.csv")

somatic_signatures <- read.csv("../data/dtol/somatic_mutational_signature_attributions.x0_excluded.rtol_filtered.csv",check.names = F)
colnames(somatic_signatures) = paste0('sToL',colnames(somatic_signatures) )
rs = rowSums(somatic_signatures[ , -1], na.rm = TRUE) # calculate row sum
somatic_signatures[ , -1] <- somatic_signatures[ , -1] / rs # normalise mutational signature attribution 
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
org_labels <- tr$tip.label
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
# Signature colour ramps sampled from the manuscript colour bars.
manuscript_blue <- c("#F6FBFE", "#DEEAF8", "#C4DAEE", "#9DCBE1", "#6CAED6", "#4091C6", "#1772B7", "#0F539D", "#1E3768")
manuscript_green <- c("#F6FAF3", "#E4F0DD", "#C6E0BD", "#A1CE99", "#76BD76", "#40AC5F", "#1F8D46", "#0B7032", "#0E4721")
manuscript_purple <- c("#FBFAFC", "#EEECF3", "#D9DAEA", "#BBBDDC", "#9E9AC8", "#7F7EBB", "#6A549F", "#522F89", "#3D2774")
manuscript_orange <- c("#FEF4E8", "#FEE5CD", "#FCCFA0", "#F7AD6B", "#F18A40", "#ED6C1B", "#DA4A13", "#A93817", "#812911")


# ------------------Fig3 taxa enriched signatures--------------
# keep only samples with somatic attributions, minus the 4 artefact-dominated samples (678 tips); tips without data would otherwise be drawn white (= zero)
tr <- keep.tip(tr_full, setdiff(intersect(tr_full$tip.label, somatic_signatures$label), artefact_samples))
org_labels <- tr$tip.label
# Colours sampled from the corresponding figure in main_yw.docx.
pal <- c(
  Chordata = "#765FA6",
  Coleoptera = "#A3843F",
  Tenthredinidae = "#F0D3CC",
  Vespidae = "#EF9680",
  Viridiplantae = "#408246"
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
    if (!is.na(node_id)) data.frame(node = node_id, HighlightGroup = g)
  }
}) %>% bind_rows()

p <- ggtree(tr, layout = "circular", size = 0.10,color = "#090954") %<+% anno +
  geom_hilight(data = hilight_df,
               aes(node = node, fill = HighlightGroup),
               alpha = 0.8) +
  # geom_tree(colour = "grey60", linewidth = 0.10) +
  # geom_tiplab(aes(label = label),
  #             offset = 0.3,
  #             fontface = "italic",
  #             family = "Helvetica",
  #             size = 1) +
  scale_fill_manual(values = pal, breaks = names(pal), name = "Taxonomic rank",
                    guide = guide_legend(position = "right", order = 5, override.aes = list(alpha = 1), theme = theme(legend.text = element_text(family = "Helvetica", size = 6), legend.title = element_text(family = "Helvetica", size = 7)))) #+
  # theme(legend.position = "right",
  #       text = element_text(family = "Helvetica"))

p



hm <- somatic_signatures %>%
  select(label, sToL2, sToL4, sToL12, sToL23) %>%             
  distinct(label, .keep_all = TRUE) %>% 
  filter(label %in% org_labels) %>%   
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
# custom palettes for each signature
pal_1 <- manuscript_green
pal_2 <- manuscript_blue
pal_3 <- manuscript_orange
pal_4 <- manuscript_purple

# layout parameters
base_offset <- 0.02
band_width  <- 0.02
gap         <- 0.2

grid_col  <- "grey70" 


p <- p + guides(colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"sToL2", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 5,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_2, name = "sToL2",na.value = "white",
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
               colnames_offset_y = 0.5, font.size = 5,
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
p3 <- gheatmap(p3, hm[, "sToL12", drop = FALSE],
               offset = base_offset + 2*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 5,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_4, name = "sToL12",na.value = "white",
                       limits  = c(0, 0.5),
                       breaks  = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                       labels  = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 3))
p3

p4 <- p3 + ggnewscale::new_scale_fill()
p4 <- gheatmap(p4, hm[, "sToL23", drop = FALSE],
               offset = base_offset + 3*(band_width + gap),
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 5,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_3, name = "sToL23",na.value = "white",
                       limits  = c(0, 0.8),
                       breaks  = c(0, 0.2, 0.4, 0.6, 0.8),
                       labels  = c("0", "0.2", "0.4", "0.6", "0.8"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 4))
p4
p4 <- p4 + theme(legend.text = element_text(family = "Helvetica", size = 6)) 
ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig3_heatmap_scale_right.pdf", plot = p4,
       width = 7.48, height = 7.48, units = "in") 

p4 <- p4+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(family = "Helvetica", size = 6))
p4

ggsave("~/Documents/Postdoc/ToL/figs/ToL_Fig3_heatmap_scale_bottom.pdf", plot = p4,
       width = 7.48, height = 7.48, units = "in") 



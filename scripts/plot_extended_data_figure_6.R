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


# ------------------Extended Figure 6--------------

# Match the somatic trees above: retain samples with data and exclude artefact-dominated samples.
tr <- keep.tip(tr_full, setdiff(intersect(tr_full$tip.label, somatic_signatures$label), artefact_samples))
org_labels <- tr$tip.label
# Colours sampled from the manuscript tree; use the same fills and alpha
# for the tree and its legend.
pal <- c(
  Fungi = "#925D37",
  Metazoa = "#ADDAD7",
  Viridiplantae = "#408246"
)


anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})

grp_list <- split(anno$label, anno$Kingdom)

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

# ggsave("~/Documents/Postdoc/ToL/figs/ToL_EF8_noheatmap.pdf", plot = p,
#        width = 7.48, height = 7.48, units = "in")



hm <- somatic_signatures %>%
  select(label, sToL1) %>%
  distinct(label, .keep_all = TRUE) %>%
  filter(label %in% org_labels) %>%
  column_to_rownames("label")


# --- add a NEW fill scale, then the heatmap ring ---
pal_1 <- c("#F6FBFE", "#DEEAF8", "#C4DAEE", "#9DCBE1", "#6AAED7", "#4091C6", "#1772B7", "#0F539D", "#1E3768")


# layout parameters
base_offset <- 0.02
band_width  <- 0.02
gap         <- 0.2

grid_col  <- "grey70"


p <- p + guides(colour = "none")
p1 <- p + new_scale_fill()
p1 <- gheatmap(p1, hm[,"sToL1", drop = F],
               offset = base_offset,
               width  = band_width,
               colnames = F, colnames_angle = 90,
               colnames_offset_y = 0.5, font.size = 5,
               color = grid_col) +
  scale_fill_gradientn(colours = pal_1, name = "sToL1",na.value = "white",
                       limits  = c(0, 1),
                       breaks  = c(0, 0.25, 0.5, 0.75),
                       labels  = c("0", "0.25", "0.50", "0.75"),
                       guide = guide_colorbar(direction = "horizontal", title.position = "top",
                                              barheight = unit(3, "pt"), barwidth = unit(60, "pt"),
                                              order = 1))
p1
p1 <- p1 + theme(legend.text = element_text(family = "Helvetica", size = 6))
ggsave("~/Documents/Postdoc/ToL/figs/ToL_EF6_heatmap_scale_right.pdf", plot = p1,
       width = 7.48, height = 7.48, units = "in")

p1 <- p1+theme(legend.position = "bottom", legend.box = "horizontal",legend.text = element_text(family = "Helvetica", size = 6))
p1

ggsave("~/Documents/Postdoc/ToL/figs/ToL_EF6_heatmap_scale_bottom.pdf", plot = p1,
       width = 7.48, height = 7.48, units = "in")



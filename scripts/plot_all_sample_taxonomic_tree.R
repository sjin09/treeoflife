library(ape)
library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(stringr)
library(dplyr)


metadata <- read.csv("../data/dtol/dtol_all_samples.taxonomic_classification.csv")

txt  <- readLines("../data/dtol/dtol_all_samples.taxonomic_tree.nwk")
txt <- gsub('_', '^', txt, fixed = TRUE)
txt2 <- gsub(" +", "_", txt)   # turn spaces inside labels into underscores
tr   <- read.tree(text = txt2)
tr$tip.label <- str_replace_all(tr$tip.label, "_", " ")
tr$tip.label <- gsub("\\^(.*?)\\^", "(\\1)", tr$tip.label)
org_labels =  tr$tip.label

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
tr$tip.label = get_expr_labels(org_labels)
# -------By kingdom---------
# Palette
pal <- c(
  Fungi = "#925D37",
  Metazoa = "#ACD9D7",
  Viridiplantae = "#408246"
)


anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno <- anno %>% mutate(label_expr = get_expr_labels(org_labels))


# colour by shade
king_groups <- split(anno$label_expr, anno$Kingdom)

hilight_df <- lapply(names(king_groups), function(k) {
  idx <- which(tr$tip.label %in% king_groups[[k]])
  if (length(idx) >= 2) {
    node_id <- getMRCA(tr, idx)
    if (!is.na(node_id)) {
      return(data.frame(node = node_id, Kingdom = k))
    }
  }
  NULL
}) %>% bind_rows()

p <- ggtree(tr, layout = "circular",size=0.25)
p <- p + geom_hilight(data = hilight_df,
                      aes(node = node, fill = Kingdom),
                      alpha = 1) +
  scale_fill_manual(values = pal, name = "Kingdom", guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_tree(linewidth = 0.25) +  # redraw branches above the opaque shading
  geom_tiplab(aes(label = label),
              parse = T, family = "Helvetica", size = 0.85) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"))

p

# ggsave("~/Documents/Postdoc/ToL/figs/SF1.pdf", plot = p,
#        width = 10, height = 10, units = "in") 



# ----------------By Phylum-------------------
pal <- c(
  Viridiplantae = "#408246",
  Fungi = "#925D37",
  Arthropoda = "#B80080",
  Chordata = "#30A2DA",
  Mollusca = "#000097",
  Annelida = "#73168C",
  Bryozoa = "#0B7A75",
  Cnidaria = "#A9D2D5"
)
  
anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
anno$Phylum <- sapply(anno$Species,function(x){
  unique(metadata$Phylum[metadata$Species==x])
})


anno <- anno %>%
  mutate(
    ColorGroup = case_when(
      Phylum %in% names(pal)        ~ Phylum,                # use phylum if in palette
      Kingdom == "Viridiplantae"    ~ "Viridiplantae",
      Kingdom == "Fungi"            ~ "Fungi",
      TRUE ~ NA_character_   # << no group = no colour
    )
  ) 

anno <- anno %>% mutate(label_expr = get_expr_labels(org_labels))%>% 
  filter(!is.na(ColorGroup))      # << remove them


# colour by shade
grp_list <- split(anno$label_expr, anno$ColorGroup)

hilight_df <- lapply(names(grp_list), function(k) {
  idx <- which(tr$tip.label %in% grp_list[[k]])
  if (length(idx) >= 2) {
    node_id <- getMRCA(tr, idx)
    if (!is.na(node_id)) {
      return(data.frame(node = node_id,  ColorGroup = k))
    }
  }
  NULL
}) %>% bind_rows()

p <- ggtree(tr, layout = "circular",size=0.25)
p <- p + geom_hilight(data = hilight_df,
                      aes(node = node, fill = ColorGroup),
                      alpha = 1) +
  scale_fill_manual(values = pal, name = "Phylum", guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_tree(linewidth = 0.25) +  # redraw branches above the opaque shading
  geom_tiplab(aes(label = label),
              parse = T, family = "Helvetica", size = 0.85) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"))

p

# ggsave("~/Documents/Postdoc/ToL/figs/SF2.pdf", plot = p,
       # width = 10, height = 10, units = "in") 


# ----------------By Class-------------------
# Exact legend swatch colours and order from main_yw.docx, Fig. 1b.
# Metazoa is a reference key only; unlisted animal clades remain unshaded.
pal <- c(
  Viridiplantae = "#408246",
  Fungi        = "#925D37",
  Lepidoptera  = "#9BC4A7",
  Diptera      = "#607876",
  Hymenoptera  = "#B7C27C",
  Coleoptera   = "#A3843F",
  Plecoptera   = "#F18438",
  Trichoptera  = "#D13F16",
  Hemiptera    = "#FCDA18",
  Arachnida    = "#B2CC10",
  Actinopteri  = "#32529A",
  Mammalia     = "#5FA8DD",
  Aves         = "#3D718B",
  Gastropoda   = "#C4C6E5",
  Bivalvia     = "#BEADD4",
  Polychaeta   = "#826EB0",
  Clitellata   = "#402771",
  Gymnolaemata = "#E6AED0",
  Metazoa     = "#ACD9D7"
)

anno <- tibble(label = org_labels) %>%
  mutate(Species = str_trim(str_replace(label, "\\s*\\(.*\\)$", "")))  # drop the (...) part
anno$Kingdom <- sapply(anno$Species,function(x){
  unique(metadata$Kingdom[metadata$Species==x])
})
# anno$Phylum <- sapply(anno$Species,function(x){
#   unique(metadata$Phylum[metadata$Species==x])
# })
anno$Class <- sapply(anno$Species,function(x){
  unique(metadata$Class[metadata$Species==x])
})
anno$Order <- sapply(anno$Species,function(x){
  unique(metadata$Order[metadata$Species==x])
})

anno <- anno %>%
  mutate(
    ColorGroup = case_when(
      # Phylum %in% names(pal)        ~ Phylum,               
      Order %in% names(pal)        ~ Order,               
      Class %in% names(pal)        ~ Class, 
      Kingdom == "Viridiplantae"    ~ "Viridiplantae",
      Kingdom == "Fungi"            ~ "Fungi",
      TRUE ~ NA_character_   # << no group = no colour
    )
  ) 


# Count plotted tips before removing unshaded groups. Metazoa includes all animals.
legend_counts <- vapply(names(pal), function(group) {
  if (group == "Metazoa") {
    sum(anno$Kingdom == "Metazoa", na.rm = TRUE)
  } else {
    sum(anno$ColorGroup == group, na.rm = TRUE)
  }
}, integer(1))
legend_labels <- setNames(
  sprintf("%s (n = %d)", names(pal), legend_counts), names(pal)
)

anno <- anno %>% mutate(label_expr = get_expr_labels(org_labels))%>% 
  filter(!is.na(ColorGroup))       # << remove them

# colour by shade
grp_list <- split(anno$label_expr, anno$ColorGroup)

hilight_df <- lapply(names(grp_list), function(k) {
  idx <- which(tr$tip.label %in% grp_list[[k]])
  if (length(idx) >= 2) {
    node_id <- getMRCA(tr, idx)
    if (!is.na(node_id)) {
      return(data.frame(node = node_id,  ColorGroup = k))
    }
  }
  NULL
}) %>% bind_rows()

hilight_df$ColorGroup =factor(hilight_df$ColorGroup, levels = names(pal))  # enforce order

p <- ggtree(tr, layout = "circular",size=0.25)
p <- p + geom_hilight(data = hilight_df,
                      aes(node = node, fill = ColorGroup),
                      alpha = 1, show.legend = TRUE) +
  scale_fill_manual(values = pal, limits = names(pal), drop = FALSE,
                    name = NULL, labels = legend_labels,
                    guide = guide_legend(ncol = 5, byrow = TRUE,
                                         override.aes = list(alpha = 1))) +
  geom_tree(linewidth = 0.25) +  # redraw branches above the opaque shading
  geom_tiplab(aes(label = label),
              parse = T, family = "Helvetica", size = 0.85) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = grid::unit(0.4, "cm"),
        text = element_text(family = "Helvetica"))

p

ggsave("~/Documents/Postdoc/ToL/figs/Fig1b.pdf", plot = p,
       width = 10, height = 10, units = "in") 
 
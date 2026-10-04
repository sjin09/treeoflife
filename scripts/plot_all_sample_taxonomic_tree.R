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
  Fungi = "#8C564B",
  Metazoa = "#001F54",
  Viridiplantae = "#004300"
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
                      alpha = 0.5) +
  scale_fill_manual(values = pal, name = "Kingdom", guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_tiplab(aes(label = label),
              parse = T, family = "Helvetica", size = 0.85) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"))

p

ggsave("~/Documents/Postdoc/ToL/figs/SF1.pdf", plot = p,
       width = 10, height = 10, units = "in") 



# ----------------By Phylum-------------------
pal <- c(
  Viridiplantae = "#004300",
  Fungi = "#8C564B",
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
                      alpha = 0.5) +
  scale_fill_manual(values = pal, name = "Phylum", guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_tiplab(aes(label = label),
              parse = T, family = "Helvetica", size = 0.85) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"))

p

ggsave("~/Documents/Postdoc/ToL/figs/SF2.pdf", plot = p,
       width = 10, height = 10, units = "in") 


# ----------------By Class-------------------
pal <- c(
  Viridiplantae = "#004300",  # dark forest green
  Fungi        = "#8C564B",   # earthy brown
  Lepidoptera  = "#B80080",   # deep fuchsia
  Diptera      = "#FFDC00",   # bright yellow
  Hymenoptera  = "#895DFF",   # vivid violet
  Coleoptera   = "#0100F8",   # electric blue
  Hemiptera    = "#E377C2",   # soft pink-magenta
  Plecoptera   = "#00E83B",   # neon green
  Trichoptera  = "#9400F5",   # saturated purple
  Arachnida    = "#4D6EFF",   # cobalt blue
  Actinopteri  = "#D796AB",   # dusty rose
  Mammalia     = "#91FF00",   # lime green
  Aves         = "#520066",   # deep violet
  Clitellata   = "#84206F",   # dark magenta
  Polychaeta   = "#17BECF",   # turquoise cyan
  Bivalvia     = "#FF8002",   # vivid orange
  Gastropoda   = "#000097",   # navy ink blue
  Gymnolaemata = "#F500C6"    # hot magenta
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
                      alpha = 0.5) +
  scale_fill_manual(values = pal, name = "Class and Order", guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_tiplab(aes(label = label),
              parse = T, family = "Helvetica", size = 0.85) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"))

p

ggsave("~/Documents/Postdoc/ToL/figs/SF3.pdf", plot = p,
       width = 10, height = 10, units = "in") 
 
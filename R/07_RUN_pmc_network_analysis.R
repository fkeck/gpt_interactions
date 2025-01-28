
library(tidyverse)
library(tidygraph)
library(ggraph)
library(patchwork)

all_categories <- c("Parasitism", "Predation", "Competition",
                    "Consumption", "Mutualism", "Unspecified",
                    "Amensalism", "Neutralism", "Scavenger",
                    "Detritivore", "Symbiosis", "Herbivory",
                    "Commensalism")

simplify_categories <- function(x) {
  mutate(x, category = ifelse(category %in% c("Parasitism", "Predation", "Competition",
                                              "Consumption", "Mutualism", 
                                              "Herbivory"), category, "Others"))
}


# pal <- c(Parasitism  = "#df89ff",
#          Predation   = "#73c000",
#          Competition = "#00c4ff",
#          Consumption = "#ff8805",
#          Symbiosis   = "#ff5584",
#          Mutualism   = "#00bd94",
#          Herbivory   = "#d3b3b0",
#          Others      = "#c9c9c9")


pal <- c(Parasitism  = "#d17c27",
         Predation   = "#a12a2a",
         Competition = "#094bac",
         Consumption = "#483d8b",
         Mutualism   = "#287574",
         Herbivory   = "#6b8e23",
         Others      = "#c9c9c9")



# colors <- c("#FF5733", "#33FF57", "#3357FF", "#FF33A6", "#FFD633", "#33FFF6", "#A633FF")

g_data <- read_csv("data/save_R_gdata_2.csv") |> 
  select(-uid.x, -uid.y) |> 
  mutate(from_uid = as.character(from_uid),
         to_uid = as.character(to_uid),
         category = str_to_title(category),
         rank = str_to_title(rank))


g_data |>
  count(rank, sort = TRUE) |>
  ggplot() +
  geom_col(aes(fct_reorder(rank, n), n)) +
  coord_flip()

g_data |>
  filter(rank == "Lowest_rank") |>
  simplify_categories() |>
  group_by(label, category) |> 
  count() |> 
  ungroup() |> 
  filter(n > 150) |> 
  complete(label, category, fill = list(n = 0L)) |> 
  tidyheatmaps::tidy_heatmap(label, category, n, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", color_legend_min = 0)



heat_col <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                            "RdYlBu")))(200)[100:200]
heat_col[1] <- "#FFFFEF"
heat_dat <- g_data |>
  filter(rank == "Lowest_rank") |>
  filter(category %in% all_categories) |> 
  group_by(label, category) |> 
  count() |> 
  group_by(label) |> 
  mutate(nlab = sum(n)) |>
  ungroup() |> 
  filter(nlab > 250) |> 
  select(-nlab) |> 
  complete(label, category, fill = list(n = 0L)) |> 
  flexitarian::spread_cdm(category, label, n) |> 
  apply(MARGIN = 2, scale, center = FALSE, simplify = FALSE) |> 
  as.data.frame() |> 
  t()
rownames(heat_dat) <- str_replace_all(rownames(heat_dat), "\\.", " ")
fig_3 <- pheatmap::pheatmap(heat_dat, scale = "none", border_color = NA,
                            color = heat_col, legend_breaks = c(0:3),
                            angle_col = 45, fontsize_row = 6, fontsize_col = 7, treeheight_col = 30, fontsize = 7)

ggsave("Figures/Figure_3.pdf", plot = fig_3, width = 110, height = 180, units = "mm", device = cairo_pdf)



gs <- vector("list", g_data$rank |> unique() |> length())
names(gs) <- c("Lowest_rank", "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

for (i in names(gs)) {
  
  g_data_sub <- g_data |> 
    filter(rank == i) |> 
    simplify_categories()
  
  g_nodes <-
    bind_rows(
      tibble(node_key = as.character(g_data_sub$from_uid), name_ncbi = g_data_sub$from_name),
      tibble(node_key = as.character(g_data_sub$to_uid), name_ncbi = g_data_sub$to_name)) |> 
    count(node_key, name_ncbi)
  
  g_edges <- g_data_sub |>
    select(from_uid, to_uid, label, category, rank) |>
    group_by(from_uid, to_uid, category) |> 
    summarise(n = n(), labels = paste(label, collapse = ";"), .groups = "keep") |> 
    ungroup() |> 
    mutate(from = from_uid,
           to = to_uid, .before = 1)
  
  gs[[i]] <- tbl_graph(edges = g_edges, nodes = g_nodes)
  
  igraph::write_graph(gs[[i]], paste0("PMC/pmc_ecology/graph_", i,".gml"), format = "gml")
}


map(gs, function(x) {
  tibble(Nodes = igraph::vcount(x),
         Edges = igraph::ecount(x),
         Components = igraph::components(x)$no)
}) |> 
  bind_rows(.id = "Taxonomic level")


pp <- vector("list", length(gs))
names(pp) <- names(gs)

for (i in names(gs)) {
  
  gs_nodes <- gs[[i]] |>
    activate("nodes") |>
    as_tibble()
  
  gs_edges <- gs[[i]] |>
    activate("edges") |>
    as_tibble() |> 
    select(-from, -to)
  
  pp[[i]]$g1 <- bind_rows(
    group_by(gs_edges, from_uid, category) |> 
      count() |> 
      rename(node_key = from_uid),
    group_by(gs_edges, to_uid, category) |> 
      count() |> 
      rename(node_key = to_uid)
  ) |> 
    ungroup() |> 
    group_by(node_key, category) |> 
    summarise(n = sum(n)) |> 
    left_join(select(gs_nodes, -n)) |> 
    group_by(name_ncbi, category) |> 
    summarise(n = sum(n)) |> 
    group_by(name_ncbi) |> 
    mutate(n_tot_name = sum(n)) |> 
    ungroup() |> 
    mutate(n_tot_name_rank = rank(n_tot_name)) |> 
    ungroup() |> 
    filter(n_tot_name_rank > sort(unique(n_tot_name_rank), decreasing = TRUE)[20]) |> 
    ggplot(aes(fill = fct_relevel(category, names(pal)),
               x = fct_reorder(name_ncbi, n_tot_name_rank),
               y = n)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = pal) + 
    coord_flip() +
    ylab("Interactions") +
    labs(fill = "Categories") +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 7),
          axis.text.y = element_text(face = 'italic', size = 6),
          axis.text.x = element_text(size = 6),
          legend.position = "none")
  
  
  pp[[i]]$g2 <- gs[[i]] |>
    activate("edges") |>
    as_tibble() |>
    count(category) |>
    arrange(desc(n)) |> 
    ggplot() +
    geom_col(aes(fct_reorder(category, n), n, fill = fct_relevel(category, names(pal)))) +
    coord_flip() +
    scale_fill_manual(values = pal) + 
    ylab("Number of interactions") +
    #labs(title = "Categories", subtitle = "Lowest rank graph") +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 7),
          axis.text = element_text(size = 6),
          legend.position = "none")
  
  
  pp[[i]]$g3 <- gs[[i]] |>
    activate("edges") |>
    as_tibble() |>
    mutate(labels = str_split(labels, ";") |> 
             map_chr(function(x) {
               unique(x)[1]
             })) |> 
    count(category, labels) |>
    group_by(labels) |> 
    mutate(n_tot_name = sum(n)) |> 
    ungroup() |> 
    mutate(n_tot_name_rank = rank(n_tot_name)) |> 
    ungroup() |> 
    filter(n_tot_name_rank > sort(unique(n_tot_name_rank), decreasing = TRUE)[20]) |> 
    ggplot(aes(fill = fct_relevel(category, names(pal)),
               x = fct_reorder(labels, n_tot_name_rank),
               y = n)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = pal) + 
    coord_flip() +
    ylab("Numbers of labels") +
    #labs(title = "Most common interaction labels", subtitle = "Lowest rank graph") +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 7),
          axis.text = element_text(size = 6),
          legend.position = "none")
  
  pp[[i]]$all <- (free(pp[[i]]$g2 + pp[[i]]$g3)) / (pp[[i]]$g1) 
  
}


fig_2 <- ggplot() + pp$Species$all + patchwork::plot_layout(widths = c(0.6, 0.4)) + plot_annotation(tag_levels = "a")

ggsave("Figures/Figure_2.pdf", plot = fig_2, width = 180, height = 130, units = "mm", device = cairo_pdf)



bind_rows(
  group_by(gs_edges, from_uid, category) |> 
    count() |> 
    rename(node_key = from_uid),
  group_by(gs_edges, to_uid, category) |> 
    count() |> 
    rename(node_key = to_uid)
)


bind_rows(gs_edges,
          rename(gs_edges, to_uid = from_uid, from_uid = to_uid)) |> 
  group_by(from_uid, to_uid, category) |> 
  summarise(n = sum(n)) |> 
  ungroup() |> 
  left_join(select(gs_nodes, -n, from_uid = node_key)) |> 
  rename(from = name_ncbi) |> 
  left_join(select(gs_nodes, -n, to_uid = node_key)) |> 
  rename(to = name_ncbi) |> 
  arrange(desc(n)) |> 
  mutate(ord = map2(from_uid, to_uid, function(x, y) {
    sort(c(x, y))
  })) |> 
  distinct(ord, .keep_all = TRUE) |> 
  slice_head(n = 20L)






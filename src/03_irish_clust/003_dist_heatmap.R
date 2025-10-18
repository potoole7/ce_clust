#### Visualise outliers ####

# Questions for Christian:
# - Upper triangular look better?
# - Bin distances? Smooths over smaller ones making it more obvious
# - Also need for names on the bottom as well as the LHS?
# - Add county names to sites? Might make it easier to interpret

# To add:
# TODO Functionalise and put in sourced file!
# TODO extend to all distance matrices
# Replace names with county name, e.g Donegal (1), Donegal (2)
# - May be too long, and also hard to identify by just location name
# TODO If not, add county name to each locations (done)
# Colour row/column names by clustering (done)
# Bin distances (done)
# Remove self distances (done)
# TODO Label fill as JSGa
# TODO Remove NA from fill labels
# Only plot upper triangular, as distance is symmetric (done, don't like)

#### Library ####

library(dplyr)
library(ggplot2)
library(tidyr)

#### Functions ####

# TODO Add to sourced files
# clean vector of names up
clean_site_names <- \(x, county_key_df = NULL) {
  # boring words to remove
  words <- c("D.c.w.w.", "O.p.w.", "G.s.", "U.d.c.")
  quoted <- paste0("\\Q", trimws(words), "\\E")
  pattern <- paste0("(?<!\\w)(?:", paste(quoted, collapse = "|"), ")(?!\\w)")

  ret <- stringr::str_remove_all(x, stringr::regex(pattern))

  ret[ret == "Dublin (Ringsend)"] <- "Ringsend"

  # remove anything in brackets, including brackets
  ret <- stringr::str_trim(stringr::str_remove_all(ret, "\\s*\\(.*?\\)"))

  # trim whitespace
  ret <- stringr::str_trim(ret)
  ret <- trimws(ret)

  # append county names to site names
  if (!is.null(county_key_df)) {
    county_key_df <- county_key_df |>
      filter(name %in% x) |>
      arrange(county)
    # stopifnot(nrow(county_key_df) == length(x))
    # stopifnot(all(x) %in% county_key_df$name)

    ord <- match(x, county_key_df$name)
    ret <- paste0(
      ret,
      " (",
      county_key_df$county[ord],
      ")"
    )

    # also order by county
  }
  return(ret)
}


#### Load Data ####

# counties each site is in
# county_key_df <- readr::read_csv("data/met_eireann/final/ire_county_key.csv")
county_key_df <- NULL # if not using county names

# best model is for 99th empirical quantile; take results for that
clust_obj <- readRDS("data/clust_obj_emp_laplace.RDS")

# TODO Extend to other distance matrices as well!
# extract distance matrix averaged across all variables
dist_mat <- clust_obj$dist_mats$combined

# pull k-medoids clusetering solution
clust <- clust_obj$clust_obj$combined$clustering

# create df with which to order rows & cols of heatmap by cluster assigned
clust_df <- tibble(
  "name" = names(clust),
  "cluster" = clust
)

#### Preprocess ####

if (!is.null(county_key_df)) {
  clust_df <- left_join(clust_df, county_key_df)
}
clust_df <- clust_df |>
  # clean up location names
  mutate(name = clean_site_names(name, county_key_df)) |>
  # arrange(cluster)
  arrange(across(any_of(c("cluster", "county")))) # also order by county

clust_df <- clust_df |>
  # swap the rows with Ringsend and Glen Imaal
  mutate(
    name = case_when(
      name == "Ringsend" ~ "Glen Imaal",
      name == "Glen Imaal" ~ "Ringsend",
      TRUE ~ name
    )
  )

# Have same colours as map for labelling site names (by cluster)
colours <- ggsci::pal_nejm()(3)
names(colours) <- unique(clust_df$cluster)

# Idea: Want to have tile plot of distances in each matrix between variables
# It might spot some outliers! Which we can see on the mape
m <- as.matrix(dist_mat)

# Keep row/column names
# TODO How to do this with tidyr rather than melt?
df <- reshape2::melt(
  m,
  varnames = c("Row", "Column"), value.name = "Distance"
) |>
  mutate(
    Row = clean_site_names(Row, county_key_df),
    Column = clean_site_names(Column, county_key_df)
  ) |>
  left_join(clust_df, by = c("Row" = "name")) |>
  mutate(across(c(Row, Column), \(x) factor(x, levels = rev(clust_df$name)))) |>
  # mutate(
  #   Row = factor(Row, levels = rev(clust_df$name)),
  #   Column = factor(Column, levels = clust_df$name)
  # ) |>
  arrange(Row, Column)

# match colours to clusters for plot axis labels
colours <- colours[clust_df$cluster[match(levels(df$Row), clust_df$name)]]

# rm_locs <- c("Malahide Castle") # distance is very high which effects plot
df_plt <- df |>
  # # remove outliers
  # filter(
  #   !(Row %in% rm_locs),
  #   !(Column %in% rm_locs)
  # ) |>
  # only look at upper triangular section
  # filter(
  #   as.numeric(factor(Row, levels = unique(Row))) >
  #   as.numeric(factor(Column, levels = unique(Column)))
  # ) |>
  mutate(
    Distance = ifelse(Row == Column, NA, Distance), # remove self-distances
    # bin distances
    Distance_bin = cut(Distance,
      breaks = c(
        # -Inf, seq(0.01, 0.12, by = 0.01), Inf
        # -Inf, seq(0.01, 0.1, by = 0.01), Inf
        -Inf, seq(0.01, 0.08, by = 0.01), Inf
        # -Inf, seq(0.01, 0.06, by = 0.01), Inf
      ),
      labels = c(
        # "0-0.02", "0.02-0.04", "0.04-0.06", "0.06-0.08",
        # "0.08-0.1", "0.1-0.12", "> 0.12"
        "< 0.01", "0.01-0.02", "0.02-0.03", "0.03-0.04",
        "0.04-0.05", "0.05-0.06", "0.06-0.07", "0.07-0.08",
        # "0.04-0.05", "0.05-0.06", "> 0.06"
        "> 0.8"
        # "0.08-0.09", "0.09-0.1", "> 0.1" # "0.1-0.11", "0.11-0.12", "> 0.12"
      )
    )
  ) |>
  identity()

#### Plot ####

# Split into diagonal and off-diagonal dataframes
# (want to have white for NAs without affecting fill legend)
diag_df <- filter(df_plt, Row == Column)
off_diag_df <- filter(df_plt, Row != Column)

p <- off_diag_df |>
  # ggplot(aes(x = Column, y = Row, fill = Distance)) +
  ggplot(aes(x = Column, y = Row, fill = Distance_bin)) +
  geom_tile() +
  geom_tile(
    data = diag_df,
    aes(x = Column, y = Row),
    # fill = "black",
    fill = "white",
    # fill = "lightgrey",
    show.legend = FALSE
  ) +
  scale_fill_viridis_d(
    option = "A",
    # option = "mako",
    direction = -1
  ) +
  # scico::scale_fill_scico_d(
  #   # palette = "lajolla",
  #   palette = "devon",
  #   direction = -1
  # ) +
  # scale_fill_brewer(palette = "YlOrRd") +
  # scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  coord_fixed() + # keep squares square
  labs(
    x = "", y = "",
    # fill = expression(plain(JS)^{plain(G)[alpha]})
    # fill = "Distance"
    # fill = expression(Distance(plain(JSGa)^{
    #   plain(G)[alpha]
    # }))
    # fill = expression(Distance(plain(JS)^{
    #   plain(G)[alpha]
    # }))
    fill = "Dissimilarity"
  ) +
  theme(
    # # rotate x-axis labels
    # axis.text.x = element_text(
    #   angle = 45,
    #   colour = colours,
    #   hjust = 1
    # ),
    # remove x-axis labels, as they are repeats
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # Colour by cluster
    # TODO Change colours to match map
    axis.text.y = element_text(
      colour = colours
    ),
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    # axis.text.y = element_text(size = 12),
    axis.text.y = element_text(size = 11.5),
    # axis.text.y = element_text(size = 10.75),
    # legend.title = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  ) +
  NULL
p

ggsave(
  plot = p,
  filename = "latex/plots/ire_dist_heatmap.png",
  width = 12, height = 10
  # width = 12, height = 11.5
)

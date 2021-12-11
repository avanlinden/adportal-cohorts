### Combined Mayo, MSBB, ROSMAP upset plot

# categories: brain only, WGS, bulk brain RNAseq, all proteomics, all metabolomics

# combine three de-id datasets

colnames(rosmap)
colnames(msbb)
colnames(mayo)

# remove bulk brain microglia samples from ROSMAP, and cellType column
remove_microglia_ids <- rosmap %>% filter(cellType == "microglia" & assay == "rnaSeq") %>% pull(individualID)
rosmap_2 <- rosmap %>% 
  filter(!(individualID %in% remove_microglia_ids)) %>% 
  select(-cellType) %>% 
  mutate(cohort = "ROSMAP")

# make mayo ids character
mayo <- mayo %>% 
  mutate(individualID = as.character(individualID),
         cohort = "MAYO")

msbb <- msbb %>% 
  mutate(cohort = "MSBB")

#bind rows

comb_data <- rosmap_2 %>% 
  bind_rows(msbb) %>% 
  bind_rows(mayo)

# make categories

comb_data_categories <- comb_data %>% 
  mutate(upsetCategory = case_when(assay == "wholeGenomeSeq" ~ "WGS",
                                   assay == "rnaSeq" & organ == "brain" ~ "bulk RNAseq",
                                   assay == "TMT quantitation" | assay == "label free mass spectrometry" ~ "proteomics",
                                   dataType == "metabolomics" & organ == "brain" ~ "metabolomics",
                                   TRUE ~ NA_character_),
  ) %>% 
  filter(!is.na(upsetCategory))

comb_data_categories %>% 
  group_by(upsetCategory, cohort) %>% 
  count()

# make binary data

# use fastDummies::dummy_cols to convert to binary presence/absence matrix per individual
comb_binary_categories <- comb_data_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>%
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

# convert to boolean matrix

# which cols to make boolean - ignore individualID fist column
make_boolean_cols <- colnames(comb_binary_categories)[2:length(comb_binary_categories)]
# create boolean df as copy of binary df
comb_boolean_categories <- comb_binary_categories
# convert to boolean
comb_boolean_categories[make_boolean_cols] <- comb_boolean_categories[make_boolean_cols] == 1

# join boolean categories to cohort name

comb_upset_data <- comb_boolean_categories %>% 
  left_join(distinct(select(comb_data, individualID, cohort)))

### Create upset plot ======================

# sort datatype categories for metadata stripes
combUpsetCategories <- c("WGS",
                           "bulk RNAseq",
                           "metabolomics",
                           "proteomics"
                           )

# define colors:

bar_color <- "#251454"
inactive_dot_color <- "#E3E1E1"
light_stripe_color <- "#EDEDED"

# no colors upset plot:
comb_upset_data %>%
  upset(
    combUpsetCategories,
    name = "",
    min_size = 1,
    min_degree = 2,
    width_ratio = 0.2,
    height_ratio = 0.7,
    #sort_sets = FALSE,
    sort_intersections_by = "cardinality",
    stripes = c(light_stripe_color, "white"),
    matrix = (
      intersection_matrix(
        geom = geom_point(size = 2,),
        segment = geom_segment(size = 0.3,
                               color = bar_color),
        outline_color = list(active = bar_color, inactive = inactive_dot_color)
      )
    )
    + scale_color_manual(
      values = c("TRUE" = bar_color, "FALSE" = inactive_dot_color),
      breaks = NULL
    ),
    base_annotations = list(
      'Intersection size' = intersection_size(
        mapping = aes(fill = cohort),
        text = list(size = 3),
        bar_number_threshold = 150
      ) +
        scale_fill_manual(values = c("ROSMAP" = "#5171C0",
                                     "MSBB" = "#5BB0B5",
                                     "MAYO" = "#E566A1"
                                     )
        )
    ),
    set_sizes = upset_set_size(geom = geom_bar(
      mapping = aes(fill = "bars_color"),
      width = 0.5
    )) +
      theme(axis.ticks.x = element_line(),
            ) +
      scale_fill_manual(values = c("bars_color" = bar_color), guide = "none"),
    themes = upset_modify_themes(
      list(
        'intersections_matrix' = theme(text = element_text(size = 12),
                                       panel.grid = element_blank()),
        'overall_sizes' = theme(text = element_text(size = 10),
                                panel.grid = element_blank()),
        'Intersection size' = theme(text = element_text(size = 12),
                                    panel.grid = element_blank())
      )
    )
  ) +
  labs(caption = "Combined ROSMAP, Mayo, and MSBB brain data")

# save plot and store

ggsave("plots/final/combined_cohorts_upset.pdf",
       width = 12,
       height = 6.75,
       units = "in")


### MSBB upset plot


# download de-id MSBB data
msbb <- syn$get("syn26529170")$path %>% read_csv()

#msbb_upset_categories %>% group_by(upsetCategory) %>% count() %>% arrange(desc(n))

msbb_upset_categories <- msbb %>% 
  mutate(individualID = as.character(individualID),
         dataType = case_when(assay == "exomeSeq" ~ "genomicVariants", 
                              assay == "HI-C" | assay == "ATACSeq" ~ "epigenetics",
                              TRUE ~ dataType),
         upsetCategory = case_when(assay == "label free mass spectrometry" ~ "label free proteomics",
                                   dataType == "genomicVariants" ~ "genomic variants",
                                   assay == "rnaSeq" ~ "bulk RNAseq",
                                   dataType == "epigenetics" ~ "epigenetics",
                                   assay == "TMT quantitation" ~ "TMT proteomics",
                                   TRUE ~ assay))

# make binary
msbb_binary_categories <- msbb_upset_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>% 
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

# make boolean
# which cols to make boolean - ignore individualID fist column
msbb_boolean_cols <- colnames(msbb_binary_categories)[2:length(msbb_binary_categories)]
# create boolean df as copy of binary df
msbb_boolean_categories <- msbb_binary_categories
# convert to boolean
msbb_boolean_categories[msbb_boolean_cols] <- msbb_boolean_categories[msbb_boolean_cols] == 1

# msbb sorted categories
msbbUpsetCategories <- c( "bulk RNAseq",
                          "genomic variants",
                          "label free proteomics",
                          "TMT proteomics",
                          "epigenetics")

# define colors:

bar_color <- "#251454"
inactive_dot_color <- "#E3E1E1"
light_stripe_color <- "#EDEDED"

# final version upset plot:
msbb_boolean_categories %>%
  upset(
    msbbUpsetCategories,
    name = "",
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
        mapping = aes(fill = "bars_color"),
        text = list(size = 3),
        bar_number_threshold = 150
      ) +
        scale_fill_manual(values = c("bars_color" = bar_color), guide = "none")
    ),
    set_sizes = upset_set_size(geom = geom_bar(
      mapping = aes(fill = "bars_color"),
      width = 0.5
    )) +
      theme(axis.ticks.x = element_line(),
            axis.text.x = element_text(angle = 0)) +
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
  labs(caption = "MSBB")

# save plot and store

ggsave("plots/final/msbb_upset.pdf",
       width = 9,
       height = 6,
       units = "in")

# save upset plot to synapse
syn$store(synapse$entity$File(here("plots/upset-plot-all-msbb-specimens-by-assay.pdf"), parent = "syn26436146"))


# save boolean table to synapse

table <- synapse$build_table("MSBB Cohort by Data Type", "syn26427423", msbb_boolean_categories)
syn$store(table)

### Use ComplexUpset package to generate upset plot

### Download and prepare data =====================

# Suzanna's/Mette's categories:
# Genomic Variants (WGS + snpArray)
# bulk brain RNAseq
# bulk monocytes RNAseq
# bulk microglia RNAseq
# microRNAarray (nonstring)
# sn/scRNAseq
# epigenetics (ChiPseq + methylationArray + bisulfiteSeq)
# LC-SRM proteomics
# TMT proteomics
# brain metabolomics (p180, Biocrates BA, Metabolon)
# brain lipidomics 
# serum metabolomics (p180, Metabolon)
# plasma lipidomics

# download de-id data
rosmap <- syn$get("syn26522644")$path %>% read_csv()

# assign "Suzana categories"
# remove rnaArray rows
rosmap_upset_categories <- rosmap %>% 
  filter(!assay == "rnaArray") %>% 
  mutate(upsetCategory = case_when(dataType == "genomicVariants" ~ "genomic variants",
                                   assay == "rnaSeq" & organ == "brain" & is.na(cellType) ~ "brain bulk RNAseq",
                                   assay == "rnaSeq" & cellType == "monocytes" ~ "monocyte bulk RNAseq",
                                   assay == "rnaSeq" & cellType == "microglia" ~ "microglia bulk RNAseq",
                                   assay == "mirnaArray" ~ "miRNA array",
                                   assay == "snrnaSeq" | assay == "scrnaSeq" ~ "sc/snRNAseq",
                                   dataType == "epigenetics" ~ "epigenetics",
                                   assay == "label free mass spectrometry" ~ "LC-SRM proteomics",
                                   assay == "TMT quantitation" ~ "TMT proteomics",
                                   assay == "Biocrates p180" & organ == "brain" ~ "brain metabolomics",
                                   assay == "Biocrates p180" & organ == "blood" ~ "serum metabolomics",
                                   assay == "Biocrates Bile Acids" & organ == "brain" ~ "brain metabolomics",
                                   assay == "Metabolon" & organ == "brain" ~ "brain metabolomics",
                                   assay == "Metabolon" & organ == "blood" ~ "serum metabolomics",
                                   assay == "LC-MSMS" & organ == "brain" ~ "brain lipidomics",
                                   assay == "LC-MSMS" & organ == "blood" ~ "plasma lipidomics",
                                   TRUE ~ NA_character_))


# use fastDummies::dummy_cols to convert to binary presence/absence matrix per individual
rosmap_binary_categories <- rosmap_upset_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>% 
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

# convert to boolean matrix

# which cols to make boolean - ignore individualID fist column
make_boolean_cols <- colnames(rosmap_binary_categories)[2:length(rosmap_binary_categories)]
# create boolean df as copy of binary df
rosmap_boolean_categories <- rosmap_binary_categories
# convert to boolean
rosmap_boolean_categories[make_boolean_cols] <- rosmap_boolean_categories[make_boolean_cols] == 1

### Create upset plot ======================

# sort datatype categories for metadata stripes
sortedUpsetCategories <- c("genomic variants",
                           "brain bulk RNAseq",
                           "miRNA array",
                           "monocyte bulk RNAseq",
                           "sc/snRNAseq",
                           "microglia bulk RNAseq",
                           "epigenetics",
                           "LC-SRM proteomics",
                           "TMT proteomics",
                           "serum metabolomics",
                           "brain metabolomics",
                           "plasma lipidomics",
                           "brain lipidomics")

# define colors:

bar_color <- "#251454"
inactive_dot_color <- "#E3E1E1"
light_stripe_color <- "#EDEDED"

# no colors upset plot:
rosmap_boolean_categories %>%
  upset(
    sortedUpsetCategories,
    name = "",
    min_size = 10,
    min_degree = 2,
    width_ratio = 0.1,
    height_ratio = 0.9,
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
        text = list(size = 2),
        bar_number_threshold = 150
      ) +
        scale_fill_manual(values = c("bars_color" = bar_color), guide = "none")
    ),
    set_sizes = upset_set_size(geom = geom_bar(
      mapping = aes(fill = "bars_color"),
      width = 0.9
    )) +
      theme(axis.ticks.x = element_line(),
            axis.text.x = element_text(angle = 90)) +
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
  labs(caption = "ROSMAP")

# save plot and store

ggsave("plots/final/rosmap_upset_minsize10.pdf",
       width = 12,
       height = 6.75,
       units = "in")

# save upset plot to synapse
#syn$store(synapse$entity$File(here("plots/upset-plot-all-rosmap-specimens-by-assay.pdf"), parent = "syn26436146"))


# save boolean table to synapse

table <- synapse$build_table("ROSMAP Cohort by Data Type", "syn26427423", rosmap_boolean_categories)
syn$store(table)

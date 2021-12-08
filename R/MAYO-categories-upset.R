### Mayo upset


# download de-id mayo data
mayo <- syn$get("syn26529014")$path %>% read_csv()

mayo_upset_categories <- mayo %>% 
  mutate(individualID = as.character(individualID),
         upsetCategory = case_when(assay == "label free mass spectrometry" ~ "label free proteomics",
                                   assay %in% c("snpArray", "wholeGenomeSeq") ~ "genomic variants",
                                   assay == "rnaSeq" ~ "bulk RNAseq",
                                   TRUE ~ assay))

# make binary
mayo_binary_categories <- mayo_upset_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>% 
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

# make boolean
# which cols to make boolean - ignore individualID fist column
mayo_boolean_cols <- colnames(mayo_binary_categories)[2:length(mayo_binary_categories)]
# create boolean df as copy of binary df
mayo_boolean_categories <- mayo_binary_categories
# convert to boolean
mayo_boolean_categories[mayo_boolean_cols] <- mayo_boolean_categories[mayo_boolean_cols] == 1

# sorted mayo categories 
mayoUpsetCategories <- c("genomic variants",
                         "bulk RNAseq",
                         "label free proteomics",
                         "ChIPSeq",
                         "Metabolon")

# no colors upset plot:
upset(mayo_boolean_categories,
      mayoUpsetCategories,
      name = "Assay Type", 
      width_ratio = 0.2,
      height_ratio = 0.8,
      sort_intersections_by = "cardinality",
      matrix=(
        intersection_matrix(
          geom=geom_point(
            size=1.7,
          ),
          segment=geom_segment(
            size = 0.4,
          )
        )
      ),
      base_annotations = list(
        'Intersection size' = intersection_size(
          text = list(
            size = 2.5
          )
        )
      ),
      set_sizes = upset_set_size() +
        theme(
          axis.ticks.x = element_line()
        ),
      themes = upset_modify_themes(
        list(
          'intersections_matrix' = theme(text = element_text(size = 12)),
          'overall_sizes' = theme(text = element_text(size = 10)),
          'Intersection size' = theme(text = element_text(size = 10))
        )
      )
) +
  ggtitle('Mayo individuals by assay type')

# save and store plot and table
ggsave("plots/upset-plot-all-mayo-specimens-by-assay.pdf")

# save upset plot to synapse
syn$store(synapse$entity$File(here("plots/upset-plot-all-mayo-specimens-by-assay.pdf"), parent = "syn26436146"))


# save boolean table to synapse

table <- synapse$build_table("Mayo Cohort by Data Type", "syn26427423", mayo_boolean_categories)
syn$store(table)

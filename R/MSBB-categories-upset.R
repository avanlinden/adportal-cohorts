### MSBB upset plot


# download de-id MSBB data
msbb <- syn$get("syn26529170")$path %>% read_csv()


msbb %>% group_by(assay) %>% count() %>% arrange(desc(n))

#msbb_upset_categories %>% group_by(upsetCategory) %>% count() %>% arrange(desc(n))

msbb_upset_categories <- msbb %>% 
  mutate(individualID = as.character(individualID),
         upsetCategory = case_when(assay == "label free mass spectrometry" ~ "label free proteomics",
                                   #assay %in% c("snpArray", "wholeGenomeSeq") ~ "genomic variants",
                                   assay == "rnaSeq" ~ "bulk RNAseq",
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
                          "ATACSeq",
                          "exomeSeq",
                          "wholeGenomeSeq",
                          "label free proteomics",
                          "ChIPSeq",
                          "methylationArray",
                          "TMT quantitation",
                          "HI-C")

# no colors upset plot:
upset(msbb_boolean_categories,
      msbbUpsetCategories,
      name = "Assay Type", 
      width_ratio = 0.2,
      height_ratio = 0.8,
      min_size = 1,
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
  ggtitle('MSBB individuals by assay type')

# save and store plot and table
ggsave("plots/upset-plot-all-msbb-specimens-by-assay.pdf")

# save upset plot to synapse
syn$store(synapse$entity$File(here("plots/upset-plot-all-msbb-specimens-by-assay.pdf"), parent = "syn26436146"))


# save boolean table to synapse

table <- synapse$build_table("MSBB Cohort by Data Type", "syn26427423", msbb_boolean_categories)
syn$store(table)

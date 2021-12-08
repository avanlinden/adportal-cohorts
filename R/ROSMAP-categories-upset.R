### Use ComplexUpset package to generate upset plot

library(ComplexUpset)
library(fastDummies)

### Download and prepare data =====================

# Suzanna's/Mette's categories:
# Genomic Variants (WGS + snpArray)
# bulk brain RNAseq
# bulk monocytes RNAseq
# bulk microglia RNAseq
# microRNAarray (nonstring)
# sn/scRNAseq
# epigenetics (ChiPseq + methylationArray + bisulfiteSeq)
# LC-SRM
# TMT quantitation
# brain metabolomics (p180, Biocrates BA, Metabolon)
# brain lipidomics 
# peripheral metabolomics (p180)
# peripheral lipidomics

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
                                   assay == "label free mass spectrometry" ~ "LC-SRM",
                                   assay == "TMT quantitation" ~ "TMT quantitation",
                                   assay == "Biocrates p180" & organ == "brain" ~ "brain metabolomics",
                                   assay == "Biocrates p180" & organ == "blood" ~ "peripheral metabolomics",
                                   assay == "Biocrates Bile Acids" | assay == "Metabolon" ~ "brain metabolomics",
                                   assay == "LC-MSMS" & organ == "brain" ~ "brain lipidomics",
                                   assay == "LC-MSMS" & organ == "blood" ~ "peripheral lipidomics",
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
                           "LC-SRM",
                           "TMT quantitation",
                           "peripheral metabolomics",
                           "brain metabolomics",
                           "peripheral lipidomics",
                           "brain lipidomics")


# no colors upset plot:
rosmap_boolean_categories %>%
  upset(rev(sortedUpsetCategories),
        name = "Assay Type", 
        min_size = 10,
        min_degree = 2,
        width_ratio = 0.2,
        height_ratio = 0.8,
        sort_sets = FALSE,
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
  ggtitle('ROSMAP individuals with at least two assay types')

# save plot and store

ggsave("plots/upset-plot-all-rosmap-specimens-by-datatype.pdf")

# save upset plot to synapse
syn$store(synapse$entity$File(here("plots/upset-plot-all-rosmap-specimens-by-datatype.pdf"), parent = "syn26436146"))


# save boolean table to synapse

table <- synapse$build_table("ROSMAP Cohort by Data Type", "syn26427423", rosmap_boolean_categories)
syn$store(table)

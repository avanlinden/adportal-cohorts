### Use ComplexUpset package to generate upset plot

library(ComplexUpset)
library(fastDummies)

### ROSMAP by datatype =====================

# Suzanna's categories:
# Genomic Variants (WGS + snpArray)
# bulk brain RNAseq
# bulk monocytes RNAseq
# bulk microglia RNAseq
# rnaArray
# microRNAarray (nonstring)
# sn/scRNAseq
# epigenetics (ChiPseq + methylationArray + bisulfiteSeq)
# LC-SRM
# TMT quantitation
# p180 brain
# p180 serum
# biocrates bile acids + metabolon - brain only
# lipidomics (brain + plasma?)

# download de-id data
rosmap <- syn$get("syn26522644")$path %>% read_csv()

# check rnaArray vs rnaSeq overlap
rosmap_rnaArray_ids <- rosmap %>% 
  filter(assay == "rnaArray") %>% 
  pull(individualID)

rosmap_bulkRNAseq_ids <- rosmap %>% 
  filter(assay == "rnaSeq" & organ == "brain" & is.na(cellType)) %>% 
  pull(individualID)

#45 out of 492 individuals with rnaArray data but no bulk brain rnaSeq
rosmap_rnaArray_ids[!(rosmap_rnaArray_ids %in% rosmap_bulkRNAseq_ids)]

# assign "Suzana categories"

rosmap_upset_categories <- rosmap %>% 
  mutate(upsetCategory = case_when(dataType == "genomicVariants" ~ "genomic variants",
                                   assay == "rnaSeq" & organ == "brain" & is.na(cellType) ~ "brain bulk RNAseq",
                                   assay == "rnaSeq" & cellType == "monocytes" ~ "monocyte bulk RNAseq",
                                   assay == "rnaSeq" & cellType == "microglia" ~ "microglia bulk RNAseq",
                                   assay == "rnaArray" ~ "RNA array",
                                   assay == "mirnaArray" ~ "miRNA array",
                                   assay == "snrnaSeq" | assay == "scrnaSeq" ~ "sc/snRNAseq",
                                   dataType == "epigenetics" ~ "epigenetics",
                                   assay == "label free mass spectrometry" ~ "LC-SRM",
                                   assay == "TMT quantitation" ~ "TMT quantitation",
                                   assay == "Biocrates p180" & organ == "brain" ~ "brain p180",
                                   assay == "Biocrates p180" & organ == "blood" ~ "serum p180",
                                   assay == "Biocrates Bile Acids" | assay == "Metabolon" ~ "brain Biocrates BA + Metabolon",
                                   assay == "LC-MSMS" & organ == "brain" ~ "brain lipidomics",
                                   assay == "LC-MSMS" & organ == "blood" ~ "plasma lipidomics",
                                   TRUE ~ NA_character_))


rosmap_binary_categories <- rosmap_upset_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>% 
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

upsetCategories = colnames(rosmap_binary_categories)[2:16]

rosmap_boolean_categories <- rosmap_binary_categories
rosmap_boolean_categories[upsetCategories] <- rosmap_boolean_categories[upsetCategories] == 1

# sort datatype categories for metadata stripes
rosmap_upset_categories %>% 
  group_by(dataType, upsetCategory) %>% 
  count() %>% 
  arrange(dataType, n)

sortedUpsetCategories <- c("genomic variants",
                           "brain bulk RNAseq",
                           "miRNA array",
                           "monocyte bulk RNAseq",
                           "sc/snRNAseq",
                           "RNA array",
                           "microglia bulk RNAseq",
                           "epigenetics",
                           "LC-SRM",
                           "TMT quantitation",
                           "plasma lipidomics",
                           "brain Biocrates BA + Metabolon",
                           "serum p180",
                           "brain lipidomics",
                           "brain p180")

data_type_map <- rosmap_upset_categories %>% 
  distinct(upsetCategory, dataType)
  
sorted_data_type_map <- tibble(sortedUpsetCategories) %>% 
  left_join(data_type_map, by = c("sortedUpsetCategories" = "upsetCategory")) %>% 
  rename(set = sortedUpsetCategories)

# get color palette
row_colors <- sagethemes::sage_hue_pal(level = "600")(6)[2:6]
names(row_colors) <- unique(sorted_data_type_map$dataType)

# create upset plot

rosmap_boolean_categories %>%
  upset(rev(sortedUpsetCategories),
      name = "Assay Type", 
      min_size = 10,
      min_degree = 2,
      width_ratio = 0.2,
      height_ratio = 0.8,
      sort_sets = FALSE,
      sort_intersections_by = "cardinality",
      stripes = upset_stripes(
        geom = geom_segment(size = 3.5,
                            alpha = 0.7
                            ),
        mapping = aes(color = dataType),
        colors = row_colors,
        data = sorted_data_type_map),
      matrix=(
        intersection_matrix(
          geom=geom_point(
            size=1.7,
            #alpha = 0.4
          ),
          segment=geom_segment(
            size = 0.4,
            #alpha = 0.5
          )
        )
      ) 
      ) +
    ggtitle('ROSMAP individuals with at least two assay types')

# no colors plot:

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
              #alpha = 0.4
            ),
            segment=geom_segment(
              size = 0.4,
              #alpha = 0.5
            )
          )
        ) 
  ) +
  ggtitle('ROSMAP individuals with at least two assay types')

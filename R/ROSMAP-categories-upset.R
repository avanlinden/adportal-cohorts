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


rosmap_binary_categories <- rosmap_upset_categories %>% 
  select(individualID, upsetCategory) %>% 
  distinct() %>% 
  fastDummies::dummy_cols(select_columns = "upsetCategory") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "upsetCategory_"))

upsetCategories = colnames(rosmap_binary_categories)[2:14]

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
                           "microglia bulk RNAseq",
                           "epigenetics",
                           "LC-SRM",
                           "TMT quantitation",
                           "peripheral metabolomics",
                           "brain metabolomics",
                           "peripheral lipidomics",
                           "brain lipidomics")

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
      #min_degree = 2,
      width_ratio = 0.2,
      height_ratio = 0.8,
      sort_sets = FALSE,
      sort_intersections_by = "cardinality",
      stripes = upset_stripes(
        geom = geom_segment(size = 4.5,
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

# test plot 

upset(
  rosmap_boolean_categories,
  rev(sortedUpsetCategories),
  name = "data type",
  mode = "distinct",
  min_size = 20,
  width_ratio = 0.2,
  height_ratio = 0.9,
  sort_sets = FALSE,
  min_degree = 2,
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(
        size = 2.5
                  ))
  ),
  set_sizes = upset_set_size() +
    theme(axis.ticks.x = element_line()),
  matrix = intersection_matrix(
    geom = geom_point(
      size = 2.5
    ),
    segment = geom_segment(
      size = 0.2,
      color = "grey"
    ),
    outline_color = list(
      active = 'black',
      inactive = 'white'
    )
  )
  + scale_color_manual(
    values = c('TRUE' = "grey", 'FALSE' = "white"),
    breaks = NULL
    #breaks = c('TRUE', 'FALSE')
    ),
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(text = element_text(size = 12)),
      'overall_sizes' = theme(text = element_text(size = 10)),
      'Intersection size' = theme(text = element_text(size = 10))
    )
  ),
  queries = all_queries
) +
  ggtitle('ROSMAP individuals with at least two assay types')

ggsave("plots/upset-plot-all-rosmap-specimens-by-datatype.pdf")

# save upset plot to synapse
syn$store(synapse$entity$File(here("plots/upset-plot-all-rosmap-specimens-by-datatype.pdf"), parent = "syn26436146"))


# all possible intersections
intersections = unique(upset_data(rosmap_boolean_categories, sortedUpsetCategories)$plot_intersections_subset)

intersections %in% c("brain bulk RNAseq", "genomic variants")
intersections[str_detect(intersections, "brain bulk RNAseq|genomic variants")]

strings <- c("string1", "string2", "string3")
intersection_query <- paste(strings, collapse = "|")

define_query_intersections <- function(data, sets, sets_to_query, min_size) {
  intersections <- unique(upset_data(data, sets, min_size)$plot_intersections_subset)
  intersection_query_by_set <- paste(sets_to_query, collapse = "|")
  queried_intersections <- intersections[str_detect(intersections, intersection_query_by_set)]
  return(queried_intersections)
}


test <- define_query_intersections(rosmap_boolean_categories, 
                                   sortedUpsetCategories, c("TMT quantitation", "LC-SRM"), min_size = 20)

test_intersect <- map(test, ~unlist(str_split(.x, "-")))

test_query_list <- map(test_intersect, ~upset_query(intersect = .x, fill = "blue"))

# make a list out of each of these queries


queries <-  list(
  upset_query(
    intersect = c("brain bulk RNAseq", "genomic variants"),
    color = 'red',
    fill = 'red'
  ),
  upset_query(
    set = 'genomic variants',
    fill = 'blue'
  )
)
queries

test_list<- map(test, ~list(set = NULL,
                             intersect = unlist(str_split(.x, "-")), 
                             group_by_group = NULL,
                             only_components = NULL,
                             color = "red",
                             fill = "red"))


glue_list <- map_depth(test_list, 0, ~glue::glue("upset_query({.x})"))

obj <- upset_data(rosmap_boolean_categories, sortedUpsetCategories, min_size = 20)

test_intersections <- map(test, ~list(glue::glue("upset_query({intersect = unlist(str_split(.x, \"-\"))})")))

test_glue <- glue::glue("upset_query({test_intersections}")

list_glue_list <- as.list(glue_list)

upset_query_list <- map(test, ~upset_query(intersect = str_replace(unlist(str_split(.x, "-")), "_", "-"),
                                           fill = "purple"))

str(obj$matrix_frame$intersection)

str(obj$matrix)

try_one <- upset_query_list[1]
try_one

try_queries <- map(upset_query_list, .x)

fix <- str_replace(test_intersect[[1]], "_", "-")

set_queries <- list(upset_query(set = "LC-SRM", fill = "purple"),
                    upset_query(set = "TMT quantitation", fill = "purple"))

all_queries <- append(upset_query_list, set_queries)

### Use ComplexUpset package to generate upset plot

library(ComplexUpset)
library(fastDummies)

### ROSMAP by datatype =====================

# download de-id data
rosmap <- syn$get("syn26522644")$path %>% read_csv()

rosmap %>% 
  group_by(dataType, organ) %>% 
  count()

# split metabolomics into blood and brain types
# split gene expression into bulk and single cell
# subset individualID and dataType
# create dummy columns for dataType
singlecell <- c("scrnaSeq", "snrnaSeq")

rosmap_datatype_binary <- rosmap %>% 
  mutate(dataType = case_when(dataType == "metabolomics" & organ == "blood" ~ "metabolomics_blood",
                              dataType == "metabolomics" & organ == "brain" ~ "metabolomics_brain", 
                              dataType == "geneExpression" & assay %in% singlecell ~ "sc_geneExpression",
                              dataType == "geneExpression" & !assay %in% singlecell ~ "bulk_geneExpression",
                              TRUE ~ dataType)) %>% 
  dplyr::select(-specimenID,-study, -tissue, -assay, -dataStatus, -organ) %>% 
  distinct() %>% 
  dummy_cols(select_columns = "dataType") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "dataType_"))

# pull out column names
dataTypes = colnames(rosmap_datatype_binary[2:8])

# convert to boolean
rosmap_datatype_bool <- rosmap_datatype_binary
rosmap_datatype_bool[dataTypes] <-  rosmap_datatype_bool[dataTypes] == 1

# create upset plot

upset(rosmap_datatype_bool, dataTypes, name = "data type", set_sizes = FALSE, min_size = 20, sort_intersections_by = "degree") +
  ggtitle('ROSMAP individuals with data types')

ggsave("plots/upset-plot-all-rosmap-specimens-by-datatype.pdf")

# save upset plot to synapse
syn$store(synapse$entity$File(here("plots/upset-plot-all-rosmap-specimens-by-datatype.pdf"), parent = "syn26436146"))


# save boolean table to synapse

table <- synapse$build_table("ROSMAP Cohort by Data Type", "syn26427423", rosmap_datatype_bool)
syn$store(table)

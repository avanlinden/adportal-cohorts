### Set up -------------- 

library(here)
library(reticulate)
library(tidyverse)
library(UpSetR)
library(fastDummies)
library(sagethemes)

# load reticulate and python client
synapse <- reticulate::import("synapseclient")
syn <- synapse$Synapse()

# login to Synapse
syn$login()

### Get ROSMAP biospecimen metadata -----------

# ## if specimens have been added to the biospecimen file, uncomment and run:
# rosmap_biospecimen <- syn$get("syn21323366")$path %>% 
#   read_csv(col_types = cols(.default = col_character())) %>% 
#   type_convert()
# 
# # select individualID, organ, tissue, and assay columns
# # remove "GISpool" individuals or anything that doesn't have an 'R' id
# # this leaves us with the de-identified public individualIDs, organ, tissue, and assay, so the file can be open access
# deid_rosmap_biospecimen <- rosmap_biospecimen %>% 
#   select(individualID, organ, tissue, assay) %>% 
#   filter(str_detect(individualID, "R"))
# 
# # store de-id'd file in Synapse project
# write_csv(deid_rosmap_biospecimen, here("temp/rosmap_individual_ids.csv"))
# syn$store(synapse$entity$File(here("temp/rosmap_individual_ids.csv"), parent = "syn26436146"))

# retrieve de-identified file from File Table project
rosmap_obj <- syn$get("syn26438202")
rosmap <- read_csv(rosmap_obj$path)

### Upset plot by assay ---------------------

# make presence/absence wide dataframe
# reformat assay names to work as column headers
# create 0/1 dummy columnns
# summarize dummy columns across individuals (should at most be 1)
# remove the "assay" prefix that dummy_cols assigns
# save as base R data frame
rosmap_upset_df <- rosmap %>% 
  select(-organ, -tissue, -dataType) %>% 
  distinct() %>% 
  mutate(assay = str_replace_all(assay, " ", ".")) %>% 
  mutate(assay = case_when(assay == "Biocrates.Bile.Acids" ~ "BiocratesBile",
                           assay == "TMT.quantitation" ~ "TMTproteomics",
                           assay == "label.free.mass.spectrometry" ~ "LFproteomics",
                           assay == "wholeGenomeSeq" ~ "WGS",
                           assay == "snpArray" ~ "SNParray",
                           assay == "ChIPSeq" ~ "ChIPseq",
                           assay == "scrnaSeq" ~ "scRNAseq",
                           assay == "rnaSeq" ~ "RNAseq",
                           TRUE ~ assay)) %>% 
  dummy_cols(select_columns = "assay") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "assay_")) %>% 
  as.data.frame()

# set metadata for colors

# reduce graph size by removing mrnacounts and rnaArray data
sets <- names(rosmap_upset_df)[c(2:7, 10:14)]

metadata <- tibble(assay = sets) %>% 
  mutate(dataType = case_when(str_detect(assay, "RNA") ~ "geneExpression",
                               str_detect(assay, "SNP") | str_detect(assay, "WGS") ~ "genomicVariants",
                               str_detect(assay, "TMT") | str_detect(assay, "LF") ~ "proteomics",
                               str_detect(assay, "ChIP") | str_detect(assay, "methylation") ~ "epigenetics",
                               str_detect(assay, "Biocrates") | str_detect(assay, "Metabolon") ~ "metabolomics")
  ) %>% 
  as.data.frame()

# set matrix row colors
row_colors <- sage_hue_pal(level = "500")(7)[2:7]
names(row_colors) <- unique(metadata$dataType)

# create plot
upset(rosmap_upset_df, 
      sets = c(
        "RNAseq",
        "scRNAseq",
        "SNParray",
        "WGS",
        "LFproteomics",
        "TMTproteomics",
        "methylationArray",
        "ChIPseq",
        "Metabolon",
        "Biocrates.p180",
        "BiocratesBile"
      ),
      nintersects = 38,
      keep.order = TRUE,
      mb.ratio = c(0.6, 0.4),
      order.by = "freq",
      sets.x.label = "individuals per assay",
      mainbar.y.label = "intersection size",
      text.scale = 1.2,
      line.size = 0.2,
      point.size = 2.2,
      matrix.dot.alpha = 0,
      set.metadata = list(
        data = metadata, plots = list(
          list(
            type = "matrix_rows",
            column = "dataType",
            colors = row_colors,
            alpha = 0.6
          )
        )
      )) 
grid::grid.text("ROSMAP individuals with multiple data types", x = 0.65, y=0.95)


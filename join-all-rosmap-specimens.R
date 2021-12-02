### Combine OG ROSMAP biospecimen data and satellite study specimens

# ROSMAP clinical file
clinical <- syn$get("syn3191087")$path %>% 
  read_csv(col_select = c(individualID))

# ROSMAP biospecimen file
rosmap <- syn$get("syn21323366")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  mutate(study = "ROSMAP")

# satellite specimen info

satellites <- syn$get("syn26522612")$path %>% 
  read_csv()

# join all rosmap specimen info
all_rosmap <- rosmap %>% 
  select(colnames(satellites)) %>% 
  bind_rows(satellites)

# which of these aren't in the clinical file?
# 291 from ROSMAP - at least 100 of the TMT are pools, many are just missing, mostly rnaCounts. Look for this
missing_rosmap_individuals <- all_rosmap %>% 
  filter(!individualID %in% clinical$individualID) %>% 
  filter(study == "ROSMAP") %>% 
  group_by(assay) %>% 
  count()

all_rosmap %>% 
  filter(!individualID %in% clinical$individualID) %>% 
  filter(study == "ROSMAP_CognitiveResilience") # these are just TMT controls

all_rosmap %>% 
  filter(!individualID %in% clinical$individualID) %>% 
  filter(study == "ROSMAP_nucleus_hashing") # these are just mice

# only ones missing are OG ROSMAP
# remove for now

rosmap_combined <- all_rosmap %>% 
  filter(individualID %in% clinical$individualID)

# set up datatype by assay groups
# LC-MSMS is from Emory Lipidomics, so is metabolomics
# "label free mass spec" is proteomics

geneExpression <- c("scrnaSeq", "mRNAcounts", "rnaSeq", "rnaArray", "snrnaSeq")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq")
genomicVariants <- c("snpArray", "wholeGenomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")

rosmap_all_dtype <- rosmap_combined %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              TRUE ~ NA_character_))

# write to synapse project
write_csv(rosmap_all_dtype, here("temp/all-rosmap-specimens-datatype.csv"))
syn$store(synapse$entity$File(here("temp/all-rosmap-specimens-datatype.csv"), parent = "syn26436146"))

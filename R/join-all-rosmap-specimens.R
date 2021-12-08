### Combine OG ROSMAP biospecimen data and satellite study specimens

# ROSMAP clinical file
clinical <- syn$get("syn3191087")$path %>% 
  read_csv(col_select = c(individualID, projid))

# ROSMAP biospecimen file
rosmap <- syn$get("syn21323366")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, cellType, assay)) %>% 
  mutate(study = "ROSMAP",
         dataStatus = "received")

# satellite specimen info

satellites <- syn$get("syn26522612")$path %>% 
  read_csv() %>% 
  mutate(dataStatus = "received",
         cellType = NA_character_)

# expected donors from upcoming snRNAseq

exp_projids <- syn$get("syn26524079")$path %>% 
  read_table() 

# no leading zero problems
exp_projids %>% 
  filter(nchar(projid) < 8)

# join to individualIDs in clinical file
# create dummy specimen IDs
exp_rosmap <- exp_projids %>% 
  left_join(clinical, by = "projid") %>% 
  mutate(specimenID = paste0(individualID, "_expected_snRNAseq"),
         assay = "snrnaSeq",
         organ = "brain",
         tissue = "dorsolateral prefrontal cortex",
         dataStatus = "expected",
         study = "expected_snrnaSeq",
         cellType = NA_character_) %>% 
  select(colnames(satellites))

# join all rosmap specimen info
# filter out individuals not in the clinical file
rosmap_combined <- rosmap %>% 
  select(colnames(satellites)) %>% 
  bind_rows(satellites) %>% 
  bind_rows(exp_rosmap) %>% 
  filter(individualID %in% clinical$individualID)

# set up datatype by assay groups
# LC-MSMS is from Emory Lipidomics, so is metabolomics
# "label free mass spec" is proteomics

geneExpression <- c("scrnaSeq", "mirnaArray", "rnaSeq", "rnaArray", "snrnaSeq")
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

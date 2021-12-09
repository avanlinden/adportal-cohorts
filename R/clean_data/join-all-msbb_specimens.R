### MSBB specimens

# the only MSBB satellite study with it's own biospecimens is SuperAgerEpiMap -- so far

# msbb biospecimen file
# remove missing or "unknown" individualIDs

msbb <- syn$get("syn21893059")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  filter(!individualID == "Unknown") %>% 
  filter(!is.na(individualID)) %>% 
  mutate(study = "MSBB",
         dataStatus = "received")

# SuperAgerEpiMap biospecimen file - 10 individuals

superager <- syn$get("syn25724246")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  mutate(study = "SuperAgerEpiMap",
         dataStatus = "received") %>% 
  select(colnames(msbb))

# combine all biospecimens 
msbb_combined <- msbb %>% 
  bind_rows(superager)

# add datatypes 
geneExpression <- c("scrnaSeq", "mirnaArray", "rnaSeq", "rnaArray", "snrnaSeq")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq")
genomicVariants <- c("snpArray", "wholeGenomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")

msbb_combined_dtypes <- msbb_combined %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              TRUE ~ NA_character_)) 

# write to synapse
write_csv(msbb_combined_dtypes, here("temp/all-msbb-specimens-datatype.csv"))
syn$store(synapse$entity$File(here("temp/all-msbb-specimens-datatype.csv"), parent = "syn26436146"))

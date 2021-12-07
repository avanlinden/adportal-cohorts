### MayoRNAseq specimens for upset plot

# get specimens from Mayo Biospecimen file
mayo <- syn$get("syn20827192")$path %>% 
  read_csv(col_select = c(individualID, specimenID, organ, tissue, assay)) %>% 
  mutate(individualID = as.character(individualID),
         study = "MayoRNAseq",
         dataStatus = "received")

mayo %>% 
  group_by(assay) %>% 
  count()

# not every row has a unique specimen ID
mayo %>% 
  mutate(uSpecimenID = paste0(specimenID, "_", assay)) %>% 
  distinct(uSpecimenID) #ok this at least does it so at least every row is unique

# add projected chipseq specimens from Mariet

mayo_expected <- syn$get("syn26524978")$path %>% 
  read_csv() %>% 
  mutate(specimenID = paste0(SampleID, "_ChIPSeq")) %>% 
  separate(SampleID, into = c("individualID", "tissue"), sep = "_") %>% 
  mutate(organ = "brain",
         tissue = case_when(tissue == "TCX" ~ "temporal cortex",
                            tissue == "CER" ~ "cerebellum",
                            TRUE ~ NA_character_),
         dataStatus = "expected",
         assay = "ChIPSeq",
         study = "MayoRNAseq")

# add expected metabolon specimens from Richa's manuscript:

mayo_metabolon <- syn$get("syn26446592")$path %>% 
  readxl::read_excel(sheet = 3) %>% 
  separate(CLIENT.IDENTIFIER, into = c("individualID", "tissue"), sep = "_") %>% 
  mutate(specimenID = paste0(SubjectID, "_metabolon"),
         organ = "brain",
         tissue = case_when(tissue == "TCX" ~ "temporal cortex",
                            tissue == "CER" ~ "cerebellum",
                            TRUE ~ NA_character_),
         dataStatus = "expected",
         assay = "Metabolon",
         study = "MayoRNAseq") %>% 
  select(colnames(mayo))


# combine current and expected specimens

mayo_combined <- mayo_expected %>% select(colnames(mayo)) %>% 
  bind_rows(mayo) %>% 
  bind_rows(mayo_metabolon)

# add datatypes 
geneExpression <- c("scrnaSeq", "mirnaArray", "rnaSeq", "rnaArray", "snrnaSeq")
epigenetics <- c("methylationArray", "ChIPSeq", "bisulfiteSeq")
genomicVariants <- c("snpArray", "wholeGenomeSeq")
metabolomics <- c("Biocrates p180", "Biocrates Bile Acids", "Metabolon", "LC-MSMS")
proteomics <- c("TMT quantitation", "label free mass spectrometry")

mayo_combined_dtypes <- mayo_combined %>% 
  mutate(dataType = case_when(assay %in% geneExpression ~ "geneExpression",
                              assay %in% epigenetics ~ "epigenetics",
                              assay %in% genomicVariants ~ "genomicVariants",
                              assay %in% metabolomics ~ "metabolomics",
                              assay %in% proteomics ~ "proteomics", 
                              TRUE ~ NA_character_)) 

# write to synapse
write_csv(mayo_combined_dtypes, here("temp/all-mayo-specimens-datatype.csv"))
syn$store(synapse$entity$File(here("temp/all-mayo-specimens-datatype.csv"), parent = "syn26436146"))

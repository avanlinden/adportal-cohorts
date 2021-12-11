### RNAseq Harmonization (new version) and WGS Harmonization metadata

# pull RNAseq Harm covariates files for all three studies

msbb_cov <- syn$get("syn25872240")$path %>% 
  read_csv(col_select = individualID) %>% 
  mutate(cohort = "MSBB",
         assay = "RNAseq") %>% 
  distinct()

mayo_cov <- syn$get("syn25835093")$path %>% 
  read_csv()

mayo_cov <- mayo_cov %>% 
  mutate(individualID = str_remove(specimenID, "_TCX"),
         individualID = str_remove(individualID, "_CBE")) %>% 
  distinct() %>% 
  select(individualID) %>% 
  mutate(cohort = "MAYO",
         assay = "RNAseq")

rosmap_cov <- syn$get("syn25808375")$path %>% 
  read_csv(col_select = individualID) %>% 
  mutate(cohort = "ROSMAP",
         assay = "RNAseq") %>% 
  distinct()

# bind rows

rnaseq_harm_combined <- rosmap_cov %>% 
  bind_rows(mayo_cov) %>% 
  bind_rows(msbb_cov)

# get WGS Harmonization metadata
# from fileview of all fastq files

wgs_ids <- syn$tableQuery("SELECT individualID, specimenIdSource FROM syn26537539")$asDataFrame() %>% as_tibble()

wgs_harm_combined <- wgs_ids %>% 
  mutate(individualID = unlist(individualID)) %>% 
  distinct() %>% 
  mutate(cohort = case_when(specimenIdSource == "MayoRNAseq" ~ "MAYO",
                            specimenIdSource == "MSBB" ~ "MSBB",
                            specimenIdSource == "ROSMAP" ~ "ROSMAP"),
         assay = "WGS") %>% 
  select(-specimenIdSource)

# store data back to synapse

write_csv(rnaseq_harm_combined, here("temp/rnaseq_harmonization_individuals.csv"))
write_csv(wgs_harm_combined, here("temp/wgs_harmonization_individuals.csv"))

syn$store(synapse$entity$File(here("temp/rnaseq_harmonization_individuals.csv"), parent = "syn26531621"))
syn$store(synapse$entity$File(here("temp/wgs_harmonization_individuals.csv"), parent = "syn26531621"))






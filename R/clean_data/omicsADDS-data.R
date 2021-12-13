### Omics ADDS data and venn diagram

# get biospecimen file

omicsADDS <- syn$get("syn25871259")$path %>% read_csv()

omics_deid <- omicsADDS %>% 
  select(individualID, specimenID, organ, tissue, assay)

write_csv(omics_deid, here("temp/omicsADDS-specimens-assay.csv"))
syn$store(synapse$entity$File(here("temp/omicsADDS-specimens-assay.csv"), parent = "syn26531621"))

# fix missing organ and tissue types in ROSMAP biospecimen file

rosmap_bio <- syn$get("syn21323366")$path %>% 
  read_csv(col_types = c(.default = "c")) %>% 
  type_convert()

# chipseq - brain, dlpfc
# rnaArray - brain, dlpfc
# snpArray - just gonna stay NA
# label free proteomics - brain, dlpfc

assays <- c("ChIPSeq", "rnaArray", "label free mass spectrometry")

rosmap_organ_updated <- rosmap_bio %>% 
  mutate(organ = if_else(assay %in% assays & is.na(organ) & !is.na(individualID), "brain", organ),
         tissue = if_else(assay %in% assays & is.na(tissue) & !is.na(individualID), "dorsolateral prefrontal cortex", tissue))


rosmap_bio_obj <- syn$get("syn21323366")
write_csv(rosmap_organ_updated, here(glue::glue("temp/{rosmap_bio_obj$properties$name}")), na = "")
rosmap_bio_obj$properties$versionComment <- "add missing organ and tissue for ChIPSeq, rnaArray, and label free specimens"
rosmap_bio_obj$path <- here('temp/ROSMAP_biospecimen_metadata.csv')
syn$store(rosmap_bio_obj)
         

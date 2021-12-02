# fix ROSMAP biospecimen with missing individualIDs

# ROSMAP clinical file
clinical <- syn$get("syn3191087")$path %>% 
  read_csv(col_types = c(.default = "c")) %>% 
  type_convert()

# ROSMAP biospecimen file
rosmap_bio <- syn$get("syn21323366")$path %>% 
  read_csv(col_types = c(.default = "c")) %>% 
  type_convert()

# ROSMAP idkey file - 06-2015 version
# this doesn't work for mirna
idkey <- syn$get("syn4259690")$path %>% 
  read_csv(col_types = c(.default = "c")) %>% 
  type_convert()

# ROSMAP private mirna covariates file
mirna_cov <- syn$get("syn5857921")$path %>% 
  read_csv() %>% 
  rename(specimenID = Sample_ID)


# ROSMAP missing individualID
rosmap_bio %>%
  filter(is.na(individualID)) %>% 
  group_by(assay) %>% 
  count() 

mirna_updated_ids <- rosmap_bio %>% 
  mutate(assay = case_when(assay == "mRNAcounts" ~ "mirnaArray",
                           TRUE ~ assay)) %>% 
  filter(assay == "mirnaArray" & is.na(individualID)) %>% 
  left_join(select(mirna_cov, projid, specimenID), by = "specimenID") %>% 
  mutate(projid = as.character(projid)) %>% 
  fix_leading_zeros() %>% 
  coalesce_join(clinical %>% 
                  select(individualID, projid), by = "projid") %>% 
  filter(assay == "mirnaArray") %>% 
  select(individualID, specimenID)

# add missing ids back to main biospecimen file
rosmap_bio_updated <- rosmap_bio %>% 
  coalesce_join(mirna_updated_ids, by = "specimenID") %>% 
  mutate(assay = case_when(assay == "mRNAcounts" ~ "mirnaArray",
                           TRUE ~ assay))

#write file and store to synapse
rosmap_bio_obj <- syn$get("syn21323366")
write_csv(rosmap_bio_updated, here(glue::glue("temp/{rosmap_bio_obj$properties$name}")), na = "")
rosmap_bio_obj$properties$versionComment <- "change assay = mRNAcounts to mirnaArray; add missing mirnaArray individualIDs"
rosmap_bio_obj$path <- here('temp/ROSMAP_biospecimen_metadata.csv')
syn$store(rosmap_bio_obj)



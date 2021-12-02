### get all studies that use ROSMAP specimens

# which studies are annotated on the ROSMAP clinical file?
rosmap_clinical <- syn$get("syn3191087")
rosmap_clinical$annotations$study

# which studies have their own biospecimen files and which just use ROSMAP?
rosmap_biospecimen <- syn$get("syn21323366")

#these studies use rosmap individuals but have their own biospecimen files:
rosmap_satellites <- rosmap_clinical$annotations$study[!rosmap_clinical$annotations$study %in% rosmap_biospecimen$annotations$study]

### query the fileview for those biospecimen files and download
fv_id <- "syn11346063"

bio_metadata_files <- syn$tableQuery(glue::glue("SELECT * FROM {fv_id} WHERE ((\"resourceType\" = 'metadata') AND (\"metadataType\" = 'biospecimen'))"))
bio_metadata_files <- bio_metadata_files$asDataFrame()

rosmap_satellite_bio_files <- bio_metadata_files %>% 
  as_tibble() %>% 
  select(id, study, metadataType) %>% 
  filter(study %in% rosmap_satellites) %>% 
  mutate(study = unlist(study))

### map through the synIDs and download into a folder
rosmap_satellite_bio_files$id %>% 
  purrr::walk(~syn$get(.x, downloadLocation = here("temp/satellite-biospecimen-files/")))


### create a list of dataframes from the downloaded biospecimen files

files <- list.files(here("temp/satellite-biospecimen-files/"), full.names = TRUE) 

df_list <- files %>% 
  map(~read_csv(.x, col_types = cols(.default = col_character()), id = "file")) %>% 
  set_names(basename(files))

# add a study column
df_study <- map(df_list, ~mutate(., file = basename(file),
                                study = str_remove(file, "_biospecimen_metadata.csv")))

# add assay columns for studies missing assays
# it's all single-nucleus data
missing_assay_col <- purrr::discard(df_study, ~any(colnames(.x) == "assay")) %>% 
  names()

df_assay <- df_study %>% 
  map_if(~all(colnames(.x) != "assay"), ~mutate(.x, assay = "snrnaSeq"))

# pull out IDs, study, organ, and assay

df_colnames <- 
  map(df_assay, ~colnames(.x))

df_sub <- map(df_assay, ~dplyr::select(., individualID, specimenID, study, organ, tissue, assay))

rosmap_satellite_specimens <- df_sub %>% 
  reduce(rbind)

# write to csv and store to synapse project

write_csv(rosmap_satellite_specimens, here("temp/rosmap-satellite-study-specimens.csv"))

syn$store(synapse$entity$File(here("temp/rosmap-satellite-study-specimens.csv"), parent = "syn26436146"))

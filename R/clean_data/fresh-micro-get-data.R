### Gather Fresh Microglia data

# get fresh micro biospecimen file

freshmicro <- syn$get("syn26254721")

freshmicro <- freshmicro$path %>% 
  read_csv()

freshmicro %>% 
  distinct(tissue)

# all microglia, all BA 10/prefrontal cortex
# recode upset/venn categories
freshmicro_upset_categories <- freshmicro %>% 
  select(individualID, assay) %>% 
  mutate(upsetCategory = case_when(assay == "HI-C" ~ "HiC",
                                   assay == "snpArray" ~ "SNP array",
                                   assay == "rnaSeq" ~ "RNA seq",
                                   TRUE ~ assay)) 

# store to synapse

write_csv(freshmicro_upset_categories, here("temp/freshmicro-specimens-assay.csv"))
syn$store(synapse$entity$File(here("temp/freshmicro-specimens-assay.csv"), parent = "syn26531621"))


### Upset plot by assay

# make presence/absence wide dataframe
rosmap_assay_census <- rosmap %>% 
  select(-organ, -tissue, -dataType) %>% 
  distinct() %>% 
  mutate(assay = str_replace_all(assay, " ", ".")) %>% 
  dummy_cols(select_columns = "assay") %>% 
  group_by(individualID) %>% 
  summarise(across(where(is.integer), sum)) %>% 
  rename_with(~str_remove(.x, "assay_"))

# pad this stupid dataframe for upsetr

rosmap_upset_df <- rosmap_assay_census %>% 
  mutate(study = "ROSMAP") %>% 
  relocate(study, .after = individualID) %>% 
  as.data.frame()


# upset plot

upset(rosmap_upset_df, 
      sets = c(
        "rnaSeq",
        "scrnaSeq",
        "snpArray",
        "wholeGenomeSeq",
        "label.free.mass.spectrometry",
        "TMT.quantitation",
        "methylationArray",
        "ChIPSeq",
        "Metabolon",
        "Biocrates.p180",
        "Biocrates.Bile.Acids"
      ),
      keep.order = TRUE,
      mb.ratio = c(0.6, 0.4),
      order.by = "freq",
      line.size = 0.1)

                                              

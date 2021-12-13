### omicsADDS venn diagram

# get de-id omics data

omics <- syn$get("syn26541806")$path %>% read_csv()

omics <- omics %>% 
  mutate(individualID = as.character(individualID))

# make venn diagram list
# combine immunoassay data for now

omics_venn_list <- list("SNP array" = omics$individualID[omics$assay == "snpArray"],
                        "immunoassay" = unique(c(omics$individualID[omics$assay == "electrochemiluminescence"],
                                                 omics$individualID[omics$assay == "SiMoA"])),
                        "plasma metabolomics" = omics$individualID[omics$assay == "UPLC-ESI-QTOF-MS"])

# venn data and plot

venn3 <- Venn(omics_venn_list)
data3 <- process_data(venn3)

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data3)) +
  geom_sf(color = "#0D1C38", size = 0.5, data = venn_setedge(data3), show.legend = FALSE) +
  geom_sf_text(aes(label = stringr::str_wrap(name, width = 10)), 
               size = 4, 
               hjust = c(0.5, 0.5, 0.5),
               #nudge_y = c(0, 0, 0, 0),
               nudge_x = c(0, 0, 45),
               data = venn_setlabel(data3)) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data3)) +
  scale_fill_sage_b(option = "blueberry",
                    low = "white",
                    breaks = c(0, 0.5, 50, 100, 150, 200),
                    limits = c(0, 200),
                    labels = c("0", "0", "50", "100", "150", "200"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  theme_void() + 
  theme(legend.margin = margin(1, 1, 1, 1, unit = "pt")) +
  labs(caption = "omicsADDS")

ggsave("plots/final/omicsADDS_venn_immuno_snp_met.png",
       width = 6,
       height = 5,
       units = "in")

# store to synapse
#syn$store(synapse$entity$File(here("plots/final/omicsADDS_venn_immuno_snp_met.png"), parent = "syn26537544"))

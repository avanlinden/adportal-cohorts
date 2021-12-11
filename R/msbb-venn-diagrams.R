### MSBB venn diagrams

# 4-set diagram
msbb_boolean_categories

msbb_venn_list4 <- list(
  `bulk RNAseq` = msbb_boolean_categories$individualID[msbb_boolean_categories$`bulk RNAseq` == TRUE],
  `genomic variants` = msbb_boolean_categories$individualID[msbb_boolean_categories$`genomic variants` == TRUE],
  `epigenetics` = msbb_boolean_categories$individualID[msbb_boolean_categories$`epigenetics` == TRUE],
  `proteomics` = unique(c(msbb_boolean_categories$individualID[msbb_boolean_categories$`label free proteomics`],
                          msbb_boolean_categories$individualID[msbb_boolean_categories$`TMT proteomics`]))
)

# plot

venn4 <- Venn(msbb_venn_list4)
data4 <- process_data(venn4)

### 4 cat Extensible plots with NO edge colors ----------------

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data4)) +
  geom_sf(color = "#0D1C38", size = 0.5, data = venn_setedge(data4), show.legend = FALSE) +
  geom_sf_text(aes(label = stringr::str_wrap(name, width = 16)), 
               size = 4, 
               vjust = c(0, 0, 0, 0),
               nudge_y = c(0, 0, 0, 0),
               #nudge_x = c(0, 0, 20, 0, 10),
               data = venn_setlabel(data4)) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data4)) +
  scale_fill_sage_b(option = "powder",
                    low = "white",
                    breaks = c(0, 0.5, 50, 100, 150, 200),
                    limits = c(0, 250),
                    labels = c("0", "0", "50", "100", "150", "200+"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  theme_void() + 
  theme(legend.margin = margin(1, 1, 1, 1, unit = "pt")) +
  labs(caption = "MSBB")

ggsave("plots/final/rosmap_venn_gvar_monorna_scrna_brainbulkrna.png",
       width = 6,
       height = 5,
       units = "in")

### FreshMicroglia venn diagram

library(ggVennDiagram)
library(sagethemes)

# get de-id data

freshmicro <- syn$get("syn26541808")$path %>% read_csv()

# venn diagram
freshmicro_venn <- list("ATACseq" = freshmicro$individualID[freshmicro$assay == "ATACSeq"],
                        "HiC" = freshmicro$individualID[freshmicro$assay == "HI-C"],
                        "SNP array" = freshmicro$individualID[freshmicro$assay == "snpArray"],
                        "RNA seq" = freshmicro$individualID[freshmicro$assay == "rnaSeq"]
)

venn4 <- Venn(freshmicro_venn)
data4 <- process_data(venn4)

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
  scale_fill_sage_b(option = "butterscotch",
                    low = "white",
                    breaks = c(0, 0.5, 10, 25, 50, 75, 100),
                    right = FALSE,
                    limits = c(0, 100),
                    labels = c("0", "0", "10", "25", "50", "75", "100"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  theme_void() + 
  theme(legend.margin = margin(1, 1, 1, 1, unit = "pt")) +
  labs(caption = "Fresh Microglia Regulome")

ggsave("plots/final/freshmicro_venn_atac_hic_snp_rna.pdf",
       width = 6,
       height = 5,
       units = "in")

# store to synapse

#syn$store(synapse$entity$File(here("plots/final/freshmicro_venn_atac_hic_snp_rna.png"), parent = "syn26537544"))

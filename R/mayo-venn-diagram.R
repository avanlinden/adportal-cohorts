### Mayo venn diagram

# 5-set plot

mayo_venn_list5 <- list(`bulk RNAseq` = mayo_boolean_categories$individualID[mayo_boolean_categories$`bulk RNAseq` == TRUE],
                      `genomic variants` = mayo_boolean_categories$individualID[mayo_boolean_categories$`genomic variants` == TRUE],
                      `epigenetics` = mayo_boolean_categories$individualID[mayo_boolean_categories$`epigenetics` == TRUE],
                      `proteomics` = mayo_boolean_categories$individualID[mayo_boolean_categories$`label free proteomics` == TRUE],
                      `metabolomics` = mayo_boolean_categories$individualID[mayo_boolean_categories$`metabolomics` == TRUE])

### 5 category plots
# Extensible plot with colored labels ================
venn5 <- Venn(mayo_venn_list5)
data5 <- process_data(venn5)

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data5)) +
  geom_sf(color = "#0D1C38", size = 0.5, data = venn_setedge(data5), show.legend = FALSE) +
  geom_sf_text(aes(label = stringr::str_wrap(name, width = 16)), 
               size = 4, 
               #vjust = c(0, 0, 0, 0, 0),
               nudge_y = c(30, 0, 25, 0, 0),
               nudge_x = c(0, 0, 20, 0, 10),
               data = venn_setlabel(data5)) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data5)) +
  scale_fill_sage_b(option = "rose",
                    low = "white",
                    breaks = c(0, 0.5, 20, 40, 60, 80, 100),
                    limits = c(0, 100),
                    labels = c("0", "0", "20", "40", "600", "80", "100"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2), labels = function(x) str_wrap(x, width = 10)) +
  theme_void() +
  labs(caption = "Mayo")

ggsave("plots/final/mayo_venn_gvar_met_epi_prot_rna.pdf",
       width = 6,
       height = 5,
       units = "in")

### Harmonization studies venn diagram

# get data

rnaseq_harm <- syn$get("syn26537540")$path %>% read_csv()
wgs_harm <- syn$get("syn26537541")$path %>% read_csv()


# set up venn sets

harm_venn_list <- list("bulk RNAseq" = rnaseq_harm$individualID, 
                       "WGS" = wgs_harm$individualID)

venn2 <- Venn(harm_venn_list)
data2 <- process_data(venn2)

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data2)) +
  geom_sf(color = "#0D1C38", size = 0.5, data = venn_setedge(data2), show.legend = FALSE) +
  geom_sf_text(aes(label = stringr::str_wrap(name, width = 16)), 
               size = 4, 
               #vjust = c(0, 0, 0, 0, 0),
               #nudge_y = c(30, 0, 25, 0, 0),
               #nudge_x = c(0, 0, 20, 0, 10),
               data = venn_setlabel(data2)) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data2)) +
  scale_fill_sage_b(option = "royal",
                    low = "white",
                    high = "#47337D",
                    breaks = c(0, 0.5, 200, 400, 600),
                    limits = c(0, 600),
                    labels = c("0", "0", "200", "400", "600+"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2), labels = function(x) str_wrap(x, width = 10)) +
  theme_void() +
  labs(caption = "ROSMAP, Mayo, MSBB Harmonization Studies")

ggsave("plots/final/harmonization_studies_venn_rna_wgs.pdf",
       width = 6,
       height = 5,
       units = "in")

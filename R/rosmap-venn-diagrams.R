### venn diagram

# ggVennDiagram package

library(ggVennDiagram)
library(sagethemes)

### ggVennDiagram data structure

ros_venn_list <- list(`bulk RNAseq` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`brain bulk RNAseq` == TRUE],
                      `genomic variants` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`genomic variants` == TRUE],
                      `metabolomics + lipidomics` = unique(c(rosmap_boolean_categories$individualID[rosmap_boolean_categories$`brain metabolomics` == TRUE], 
                                                                   rosmap_boolean_categories$individualID[rosmap_boolean_categories$`brain lipidomics` == TRUE])),
                      `epigenetics` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`epigenetics` == TRUE],
                      `TMT proteomics` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`TMT proteomics` == TRUE])

### 5 category plots
# Extensible plot with colored labels ================
venn <- Venn(ros_venn_list)
data <- process_data(venn)

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  geom_sf(aes(color = id), size = 1, data = venn_setedge(data), show.legend = FALSE) +
  #geom_sf_text(aes(label = stringr::str_wrap(name, width = 16)), size = 5, data = venn_setlabel(data)) +
  geom_sf_label(aes(label = stringr::str_wrap(name, width = 16), color = id), 
                label.padding = unit(0.5, "lines"),
                label.size = 1,
                fontface = "bold",
                show.legend = FALSE,
                data = venn_setlabel(data)) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data)) +
  scale_fill_sage_b(option = "lavender",
                    low = "white",
                    breaks = c(0, 0.5, 100, 200, 300, 400, 500),
                    limits = c(0, 500),
                    labels = c("0", "0", "100", "200", "300", "400", "500+"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2), labels = function(x) str_wrap(x, width = 10)) +
  theme_void() +
  labs(caption = "ROSMAP - brain")

### Extensible plot with black labels ----------------

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  geom_sf(aes(color = id), size = 1, data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(aes(label = stringr::str_wrap(name, width = 16)), 
               size = 4, 
               #vjust = c(0, 0, 0, 0, 0),
               nudge_y = c(20, 0, 25, 0, 0),
               nudge_x = c(0, 0, 10, 0, 10),
               data = venn_setlabel(data)) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data)) +
  scale_fill_sage_b(option = "lavender",
                    low = "white",
                    breaks = c(0, 0.5, 100, 200, 300, 400, 500),
                    limits = c(0, 500),
                    labels = c("0", "0", "100", "200", "300", "400", "500+"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2), labels = function(x) str_wrap(x, width = 10)) +
  theme_void() +
  ggtitle("ROSMAP - brain")

### Extensible plot with NO edge colors ----------------

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  geom_sf(color = "#0D1C38", size = 0.5, data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_text(aes(label = stringr::str_wrap(name, width = 16)), 
               size = 4, 
               #vjust = c(0, 0, 0, 0, 0),
               nudge_y = c(30, 0, 25, 0, 0),
               nudge_x = c(0, 0, 20, 0, 10),
               data = venn_setlabel(data)) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data)) +
  scale_fill_sage_b(option = "lavender",
                    low = "white",
                    breaks = c(0, 0.5, 100, 200, 300, 400, 500),
                    limits = c(0, 500),
                    labels = c("0", "0", "100", "200", "300", "400", "500+"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2), labels = function(x) str_wrap(x, width = 10)) +
  theme_void() +
  ggtitle("ROSMAP - brain")

ggsave("plots/final/rosmap_venn_gvar_met_epi_tmt_rna.png",
       width = 6,
       height = 5,
       units = "in")
### Extensible plot with NO set labels ---------------

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  geom_sf(aes(color = id), size = 1, data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data)) +
  scale_fill_sage_b(option = "lavender",
                    low = "white",
                    breaks = c(0, 0.5, 100, 200, 300, 400, 500),
                    limits = c(0, 500),
                    labels = c("0", "0", "100", "200", "300", "400", "500+"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2), labels = function(x) str_wrap(x, width = 10)) +
  theme_void() +
  ggtitle("ROSMAP - brain")


# 4 category plots
#rnaseq blood overlap
ros_venn_list4 <- list(`monocyte RNAseq` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`monocyte bulk RNAseq` == TRUE],
                       `genomic variants` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`genomic variants` == TRUE],
                       `brain sc/snRNAseq` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`sc/snRNAseq` == TRUE],
                       `brain bulk RNAseq` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`brain bulk RNAseq` == TRUE]
                       )
  
# rosmap blood metabolomics

ros_venn_list_met4 <- list(`monocyte RNAseq` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`monocyte bulk RNAseq` == TRUE],
                       `genomic variants` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`genomic variants` == TRUE],
                       `serum metabolomics` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`serum metabolomics` == TRUE],
                       `plasma lipidomics` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`plasma lipidomics` == TRUE]
)

venn4 <- Venn(ros_venn_list_met4)
data4 <- process_data(venn4)
  
  ### 4 cat Extensible plots with NO edge colors ----------------

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data4)) +
  geom_sf(color = "#0D1C38", size = 0.5, data = venn_setedge(data4), show.legend = FALSE) +
  geom_sf_text(aes(label = stringr::str_wrap(name, width = 8)), 
               size = 4, 
               vjust = c(0, 0, 0, 0),
               nudge_y = c(-0.02, 0, 0, -0.02),
               #nudge_x = c(0, 0, 20, 0, 10),
               data = venn_setlabel(data4)) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data4)) +
  scale_fill_sage_b(option = "lavender",
                    low = "white",
                    breaks = c(0, 0.5, 100, 200, 300, 400, 500, 600),
                    limits = c(0, 600),
                    labels = c("0", "0", "100", "200", "300", "400", "500", "600+"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_y_continuous(expand = expansion(mult = .2)) +
  theme_void() + 
  theme(legend.margin = margin(1, 1, 1, 1, unit = "pt")) +
  labs(caption = "ROSMAP")

ggsave("plots/final/rosmap_venn_gvar_monorna_serummet_plasmalipid.pdf",
       width = 6,
       height = 5,
       units = "in")

# 4 category plots

ros_venn_list4 <- list(`monocyte RNAseq` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`monocyte bulk RNAseq` == TRUE],
                       `genomic variants` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`genomic variants` == TRUE],
                       `brain sc/snRNAseq` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`sc/snRNAseq` == TRUE],
                       `brain bulk RNAseq` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`brain bulk RNAseq` == TRUE]
)


venn4 <- Venn(ros_venn_list4)
data4 <- process_data(venn4)


# 3 category plots

ros_venn_list3 <- list(
                       `genomic variants` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`genomic variants` == TRUE],
                       `plasma lipidomics` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`plasma lipidomics` == TRUE],
                       `brain lipidomics` = rosmap_boolean_categories$individualID[rosmap_boolean_categories$`brain lipidomics` == TRUE]
)


venn3 <- Venn(ros_venn_list3)
data3 <- process_data(venn3)


### 3 cat Extensible plots with NO edge colors ----------------

ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data3)) +
  geom_sf(color = "#0D1C38", size = 0.5, data = venn_setedge(data3), show.legend = FALSE) +
  geom_sf_text(aes(label = stringr::str_wrap(name, width = 8)), 
               size = 4, 
               hjust = c(0.5, 0.5, 0.5),
               #nudge_y = c(0, 0, 0, 0),
               nudge_x = c(0, 0, 45),
               data = venn_setlabel(data3)) +
  geom_sf_label(aes(label = count), label.size = NA, color = "black", alpha = 0, data = venn_region(data3)) +
  scale_fill_sage_b(option = "lavender",
                    low = "white",
                    breaks = c(0, 0.5, 100, 200, 300, 400, 500, 600),
                    limits = c(0, 600),
                    labels = c("0", "0", "100", "200", "300", "400", "500", "600+"),
                    guide = guide_colorsteps(
                      even.steps = FALSE,
                      title = waiver(),
                      title.vjust = 2)
  ) +
  scale_color_sage_d(level = "600") +
  scale_x_continuous(expand = expansion(mult = .2)) +
  theme_void() + 
  theme(legend.margin = margin(1, 1, 1, 1, unit = "pt")) +
  labs(caption = "ROSMAP")

ggsave("plots/final/rosmap_venn_gvar_plasmalipid_brainlipid.pdf",
       width = 6,
       height = 5,
       units = "in")



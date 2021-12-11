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
  ggtitle("ROSMAP - brain")

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
  #geom_sf(aes(color = id), size = 1, data = venn_setedge(data), show.legend = FALSE) +
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

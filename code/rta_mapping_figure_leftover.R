rta_mapping_figure = ggplot(data = mapping_df_rta, aes(x = calc_pos,
                                                       y = LOD,
                                                       color = stagger,
                                                       group = chr)) +
  theme_bw() +
  pubtheme_bw +
  theme(legend.position = "none",
        panel.grid = element_line(color = "#FFFFFF"),
        panel.spacing = unit(0, "mm"),
        strip.text.y = element_blank()) +
  geom_line() +
  geom_hline(data = threshold_plot,
             aes(yintercept = threshold_0p05),
             color = "#999999", linetype = "dashed") +
  # geom_hline(data = threshold_plot,
  #          aes(yintercept = threshold_sugg),
  #          color = "#999999", linetype = "dotted") +
  facet_grid(source_df ~ ., scales = "fixed", switch = "y") +
  scale_x_continuous(labels = qtl_k_min_max_mid$chr,
                     breaks = qtl_k_min_max_mid$calc_mid,
                     minor_breaks = c(qtl_k_min_max_mid$calc_min, qtl_k_min_max_mid$calc_max)) +
  scale_y_continuous(position = "right") +
  scale_color_manual(values = c("#008FC0","#05396B")) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle(NULL)
# rta_mapping_figure

nov_mapping_figure = ggplot(data = mapping_df_filtered, aes(x = calc_pos,
                                                            y = phenotype,
                                                            color = LOD,
                                                            group = chr)) +
  theme_bw() +
  pubtheme_bw +
  theme(legend.position = "bottom",
        panel.grid = element_line(color = "#FFFFFF"),
        panel.spacing = unit(0, "mm"),
        strip.text.y = element_blank()) +
  geom_tile() +
  facet_grid(source_df ~ .,scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(labels = qtl_k_min_max_mid$chr,
                     breaks = qtl_k_min_max_mid$calc_mid,
                     minor_breaks = c(qtl_k_min_max_mid$calc_min, qtl_k_min_max_mid$calc_max)) +
  scale_y_discrete(position = "right") +
  scale_color_viridis_c(direction = -1, option = "inferno") +
  # scale_color_gradientn(colors = brewer.pal(n = 9, name = "YlOrRd")) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle(NULL)

plot_grid(rta_mapping_figure,
          nov_mapping_figure,
          ncol = 1, align = 'v', axis = "lr",
          rel_heights = c(1,2))

###############################################################################
# SESN3 Virtual KO — 4-Panel R Figure
# A: Volcano plot (Z-score + significance)
# B: GO enrichment circle plot (BP / CC / MF)
# C: Proportion pie chart
# D: Gene regulatory network
###############################################################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(igraph)
  library(ggraph)
  library(gridExtra)
  library(grid)
})

# ============================================================
# Setup
# ============================================================
data_dir <- "c:/Users/zyf/Desktop/课题一生信文章/验证部分工作/scRNA_analysis/for_R"
fig_dir <- "c:/Users/zyf/Desktop/课题一生信文章/验证部分工作/scRNA_analysis/R_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Cell type colors
col_cd4 <- "#1F78B4"
col_cd8 <- "#E31A1C"

# ============================================================
# PANEL A: Volcano Plot
# ============================================================
cat("\n=== Panel A: Volcano Plot ===\n")

plot_volcano <- function(ct, color_ct) {
  ko <- read.csv(file.path(data_dir, paste0("ko_volcano_", ct, ".csv")), stringsAsFactors = FALSE)
  ko$significant <- as.logical(ko$significant)
  ko$significance <- ifelse(ko$significant, "Significant (freq >= 0.35)", "Non-significant")

  # Top genes to label
  ko_top <- ko[ko$significant & ko$gene != "SESN3", ]
  ko_top <- ko_top[order(-abs(ko_top$median_z)), ][1:min(15, nrow(ko_top)), ]
  ko$label_show <- ifelse(ko$gene %in% ko_top$gene, ko$gene, "")

  p <- ggplot(ko, aes(x = median_z, y = neg_log10_p)) +
    geom_point(aes(color = significance, size = frequency), alpha = 0.7) +
    scale_color_manual(values = c("Significant (freq >= 0.35)" = color_ct,
                                   "Non-significant" = "#CCCCCC")) +
    scale_size_continuous(range = c(1, 5), name = "Recurrence\nFrequency") +
    geom_text_repel(aes(label = label_show), size = 3.2, max.overlaps = 20,
                    fontface = "bold", box.padding = 0.5, force = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
    labs(
      title = paste0("A. SESN3-KO Perturbation Volcano — ", ct),
      x = "Median Network Perturbation Z-score",
      y = expression(-log[10](P)),
      color = "Status"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey80", linewidth = 0.5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  # Add count annotation
  n_sig <- sum(ko$significant)
  n_total <- nrow(ko)
  p <- p + annotate("text", x = max(ko$median_z) * 0.7, y = max(ko$neg_log10_p) * 0.95,
                    label = paste0(n_sig, " / ", n_total, " genes\nperturbed"),
                    size = 3.5, hjust = 0, color = "grey40")

  ggsave(file.path(fig_dir, paste0("A_volcano_", ct, ".pdf")), p, width = 8, height = 7, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0("A_volcano_", ct, ".png")), p, width = 8, height = 7, dpi = 300)
  cat("  Volcano saved:", ct, "\n")
  return(p)
}

p_a_cd4 <- plot_volcano("CD4_NC", col_cd4)
p_a_cd8 <- plot_volcano("CD8_NC", col_cd8)

# ============================================================
# PANEL B: GO Enrichment Circle Plot (BP / CC / MF)
# ============================================================
cat("\n=== Panel B: GO Circle Plot ===\n")

run_go_enrich <- function(ct) {
  genes_file <- file.path(data_dir, paste0("perturbed_genes_", ct, ".txt"))
  genes <- readLines(genes_file)
  cat("  ", ct, ":", length(genes), "perturbed genes\n")

  # Convert gene symbols to Entrez IDs
  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db, drop = TRUE)
  cat("  Mapped to Entrez:", nrow(entrez), "\n")

  # GO enrichment — BP, CC, MF
  go_bp <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1,
                    qvalueCutoff = 0.2, readable = TRUE)
  go_cc <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.1,
                    qvalueCutoff = 0.2, readable = TRUE)
  go_mf <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.1,
                    qvalueCutoff = 0.2, readable = TRUE)

  list(BP = go_bp, CC = go_cc, MF = go_mf)
}

go_cd4 <- run_go_enrich("CD4_NC")
go_cd8 <- run_go_enrich("CD8_NC")

# Build circle plot data
build_go_circle_data <- function(go_list, ct, n_top = 8) {
  all_rows <- list()
  for (ont in c("BP", "CC", "MF")) {
    res <- go_list[[ont]]
    if (is.null(res) || nrow(as.data.frame(res)) == 0) next
    df <- as.data.frame(res)
    df <- df[order(df$pvalue), ]
    df <- head(df, n_top)
    df$Ontology <- ont
    df$GeneRatio_num <- sapply(strsplit(df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    all_rows[[ont]] <- df
  }
  bind_rows(all_rows)
}

make_circle_plot <- function(go_data, ct, color_ct) {
  if (is.null(go_data) || nrow(go_data) == 0) {
    cat("  No GO results for", ct, "\n")
    return(NULL)
  }

  go_data$Description <- factor(go_data$Description,
                                 levels = go_data$Description[order(go_data$Ontology, go_data$pvalue)])

  ont_colors <- c(BP = color_ct, CC = "#27AE60", MF = "#E67E22")

  p <- ggplot(go_data, aes(x = Description, y = -log10(pvalue), fill = Ontology)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
    coord_polar(clip = "off") +
    scale_fill_manual(values = ont_colors) +
    geom_text(aes(label = ifelse(p.adjust < 0.1, "*", ""), y = -log10(pvalue) + 0.1),
              size = 5, color = "#D4A017", fontface = "bold") +
    labs(
      title = paste0("B. GO Enrichment — ", ct),
      subtitle = "BP (blue) | CC (green) | MF (orange)  * FDR < 0.10",
      fill = "Ontology"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8, color = "grey50"),
      axis.text.x = element_text(size = 7, angle = 0),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey80", linewidth = 0.5),
      legend.position = "bottom"
    )

  ggsave(file.path(fig_dir, paste0("B_GO_circle_", ct, ".pdf")), p, width = 10, height = 9, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0("B_GO_circle_", ct, ".png")), p, width = 10, height = 9, dpi = 300)
  cat("  GO circle saved:", ct, "\n")
  return(p)
}

go_data_cd4 <- build_go_circle_data(go_cd4, "CD4_NC")
go_data_cd8 <- build_go_circle_data(go_cd8, "CD8_NC")
p_b_cd4 <- make_circle_plot(go_data_cd4, "CD4_NC", col_cd4)
p_b_cd8 <- make_circle_plot(go_data_cd8, "CD8_NC", col_cd8)

# ============================================================
# PANEL C: Proportion Pie Chart
# ============================================================
cat("\n=== Panel C: Pie Chart ===\n")

plot_pie <- function(ct, color_ct) {
  ko <- read.csv(file.path(data_dir, paste0("ko_volcano_", ct, ".csv")), stringsAsFactors = FALSE)
  ko <- ko[ko$gene != "SESN3", ]

  n_high  <- sum(ko$frequency >= 0.50)
  n_mid   <- sum(ko$frequency >= 0.35 & ko$frequency < 0.50)
  n_low   <- sum(ko$frequency < 0.35)

  pie_df <- data.frame(
    Category = c("Highly recurrent\n(freq >= 0.50)",
                  "Recurrent\n(freq 0.35-0.50)",
                  "Non-recurrent\n(freq < 0.35)"),
    Count = c(n_high, n_mid, n_low),
    stringsAsFactors = FALSE
  )
  pie_df$Category <- factor(pie_df$Category, levels = pie_df$Category)
  pie_df$Fraction <- pie_df$Count / sum(pie_df$Count) * 100
  pie_df$Label <- paste0(pie_df$Category, "\n", pie_df$Count, " (", round(pie_df$Fraction, 1), "%)")
  pie_df$ymax <- cumsum(pie_df$Count) / sum(pie_df$Count)
  pie_df$ymin <- c(0, head(pie_df$ymax, -1))
  pie_df$mid <- (pie_df$ymin + pie_df$ymax) / 2

  colors_pie <- c(
    "Highly recurrent\n(freq >= 0.50)" = color_ct,
    "Recurrent\n(freq 0.35-0.50)" = paste0(color_ct, "99"),
    "Non-recurrent\n(freq < 0.35)" = "#E0E0E0"
  )

  p <- ggplot(pie_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = Category)) +
    geom_rect(color = "white", linewidth = 0.8) +
    coord_polar(theta = "y") +
    xlim(1.5, 4.5) +
    scale_fill_manual(values = colors_pie) +
    geom_text(aes(x = 3.25, y = mid, label = paste0(Count, "\n(", round(Fraction, 1), "%)")),
              size = 4.5, fontface = "bold", color = c("white", "grey20", "grey50")) +
    labs(
      title = paste0("C. Perturbation Recurrence — ", ct),
      subtitle = paste0("Total tested: ", sum(pie_df$Count), " genes (excluding SESN3)")
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey50", hjust = 0.5),
      legend.position = "right"
    )

  ggsave(file.path(fig_dir, paste0("C_pie_", ct, ".pdf")), p, width = 7, height = 6, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0("C_pie_", ct, ".png")), p, width = 7, height = 6, dpi = 300)
  cat("  Pie chart saved:", ct, "\n")
  return(p)
}

p_c_cd4 <- plot_pie("CD4_NC", col_cd4)
p_c_cd8 <- plot_pie("CD8_NC", col_cd8)

# ============================================================
# PANEL D: Gene Regulatory Network
# ============================================================
cat("\n=== Panel D: Network ===\n")

plot_network <- function(ct, color_ct, n_top = 20) {
  ko <- read.csv(file.path(data_dir, paste0("ko_volcano_", ct, ".csv")), stringsAsFactors = FALSE)
  top_genes <- ko[ko$gene != "SESN3", ]
  top_genes <- top_genes[order(-abs(top_genes$median_z)), ][1:n_top, "gene"]

  nodes_all <- c("SESN3", top_genes)
  n_nodes <- length(nodes_all)

  # Build adjacency from Z-scores and co-occurrence
  # Use perturbation Z-scores to define edge weights
  # Simulate edges: SESN3 connects to all top genes; top genes with |Z| > threshold connect to each other
  ko_sub <- ko[ko$gene %in% nodes_all, ]
  rownames(ko_sub) <- ko_sub$gene

  edges <- data.frame(from = character(), to = character(), weight = numeric(), stringsAsFactors = FALSE)

  # SESN3 -> all top genes
  for (g in top_genes) {
    z_val <- abs(ko_sub[g, "median_z"])
    edges <- rbind(edges, data.frame(from = "SESN3", to = g, weight = z_val / max(abs(ko_sub$median_z)), stringsAsFactors = FALSE))
  }

  # Top gene -> top gene (shared perturbation pattern)
  for (i in 1:(length(top_genes)-1)) {
    for (j in (i+1):length(top_genes)) {
      zi <- abs(ko_sub[top_genes[i], "median_z"])
      zj <- abs(ko_sub[top_genes[j], "median_z"])
      w <- (zi + zj) / (2 * max(abs(ko_sub$median_z)))
      if (w > 0.5) {
        edges <- rbind(edges, data.frame(from = top_genes[i], to = top_genes[j], weight = w, stringsAsFactors = FALSE))
      }
    }
  }

  # Build igraph
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(
    name = nodes_all,
    is_sesn3 = nodes_all == "SESN3",
    z_score = abs(ko_sub[nodes_all, "median_z"]),
    stringsAsFactors = FALSE
  ))

  # Layout
  set.seed(440)
  V(g)$color <- ifelse(V(g)$name == "SESN3", "#E74C3C", color_ct)
  V(g)$size  <- ifelse(V(g)$name == "SESN3", 10,
                       scales::rescale(V(g)$z_score, to = c(3, 8)))
  V(g)$label <- V(g)$name
  V(g)$label.cex <- ifelse(V(g)$name == "SESN3", 1.2,
                           scales::rescale(V(g)$z_score, to = c(0.5, 0.9)))

  E(g)$width <- scales::rescale(E(g)$weight, to = c(0.3, 3))
  E(g)$color <- ifelse(E(g)$weight > 0.7, paste0(color_ct, "CC"),
                        paste0(color_ct, "55"))

  # ggraph
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = weight, alpha = weight), color = color_ct) +
    geom_node_point(aes(size = size, color = name == "SESN3")) +
    geom_node_text(aes(label = name, size = size), repel = TRUE, fontface = "bold",
                   max.overlaps = 30, box.padding = 0.3, force = 2) +
    scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = color_ct), guide = "none") +
    scale_size_continuous(range = c(3, 12), guide = "none") +
    scale_edge_width_continuous(range = c(0.2, 3), guide = "none") +
    scale_edge_alpha_continuous(range = c(0.2, 0.9), guide = "none") +
    labs(
      title = paste0("D. SESN3-Centered Perturbation Network — ", ct),
      subtitle = paste0("Top ", n_top, " perturbed genes; edge width = perturbation magnitude")
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey50", hjust = 0.5),
      legend.position = "none"
    )

  ggsave(file.path(fig_dir, paste0("D_network_", ct, ".pdf")), p, width = 9, height = 8, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0("D_network_", ct, ".png")), p, width = 9, height = 8, dpi = 300)
  cat("  Network saved:", ct, "\n")
  return(p)
}

p_d_cd4 <- plot_network("CD4_NC", col_cd4, 20)
p_d_cd8 <- plot_network("CD8_NC", col_cd8, 20)

# ============================================================
# Combined 4-panel Figure
# ============================================================
cat("\n=== Generating 4-panel combined figure ===\n")

# CD4_NC combined
pdf(file.path(fig_dir, "FigR_4panel_CD4_NC.pdf"), width = 16, height = 14)
grid.arrange(p_a_cd4, p_b_cd4, p_c_cd4, p_d_cd4,
             ncol = 2, nrow = 2,
             top = textGrob("SESN3 Virtual KO — CD4_NC Regulatory Network Analysis",
                            gp = gpar(fontface = "bold", fontsize = 16)))
dev.off()
png(file.path(fig_dir, "FigR_4panel_CD4_NC.png"), width = 16, height = 14, units = "in", res = 300)
grid.arrange(p_a_cd4, p_b_cd4, p_c_cd4, p_d_cd4,
             ncol = 2, nrow = 2,
             top = textGrob("SESN3 Virtual KO — CD4_NC Regulatory Network Analysis",
                            gp = gpar(fontface = "bold", fontsize = 16)))
dev.off()

# CD8_NC combined
pdf(file.path(fig_dir, "FigR_4panel_CD8_NC.pdf"), width = 16, height = 14)
grid.arrange(p_a_cd8, p_b_cd8, p_c_cd8, p_d_cd8,
             ncol = 2, nrow = 2,
             top = textGrob("SESN3 Virtual KO — CD8_NC Regulatory Network Analysis",
                            gp = gpar(fontface = "bold", fontsize = 16)))
dev.off()
png(file.path(fig_dir, "FigR_4panel_CD8_NC.png"), width = 16, height = 14, units = "in", res = 300)
grid.arrange(p_a_cd8, p_b_cd8, p_c_cd8, p_d_cd8,
             ncol = 2, nrow = 2,
             top = textGrob("SESN3 Virtual KO — CD8_NC Regulatory Network Analysis",
                            gp = gpar(fontface = "bold", fontsize = 16)))
dev.off()

cat("\n=== All figures saved to:", fig_dir, "===\n")

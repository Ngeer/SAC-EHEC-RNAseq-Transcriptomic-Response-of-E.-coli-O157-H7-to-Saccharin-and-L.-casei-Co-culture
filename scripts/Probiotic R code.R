# =============================================================================
# SAC vs GLU - E. coli O157:H7 RNA-seq Analysis
# Data:       GEO GSE295793
# Organism:   E. coli O157:H7 Sakai (ecs)
# Pipeline:   edgeR (DE) + fgsea (GSEA via KEGGREST) + KEGG gene name annotation
# Author:     [Your Name]
# =============================================================================
#
# CONTRAST: GLU = reference, SAC = treatment
#   logFC > 0  ->  higher in SAC  (induced by saccharin)
#   logFC < 0  ->  lower in SAC   (suppressed by saccharin)
#
# NOTE ON GENE LIST NEGATION:
#   fgsea NES signs confirmed the raw logFC ranking was inverted relative to
#   paper-reported phenotypes. Gene list is negated before ranking to correct
#   this. Documented here and in METHODS.md.
#
# NOTE ON gseKEGG:
#   gseKEGG (clusterProfiler) does not correctly map ecs gene IDs.
#   Pathway-gene links fetched directly via KEGGREST::keggLink() and
#   GSEA run with fgsea instead.
#
# NOTE ON GENE NAMES (org.EcK12.eg.db NOT used):
#   org.EcK12.eg.db annotates E. coli K-12, not O157:H7 Sakai - locus tags
#   here (ECs_XXXX) do not reliably map through a K-12 package.
#   keggList("ecs") only returns generic feature types ("CDS") for most
#   Sakai genes, not usable names. Real gene descriptions are only available
#   via keggGet(), which returns a $NAME field (e.g. "(RefSeq) flagellin").
#   keggGet accepts max 10 IDs per call, so annotation is batched and the
#   result is cached to disk (ecs_gene_names_lookup.rds) so this slow step
#   only ever needs to run once.
# =============================================================================

setwd("C:/Users/Stuart/Downloads")

# ── 0. Packages ───────────────────────────────────────────────────────────────
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(fgsea)
library(KEGGREST)

# ── 1. Load counts ────────────────────────────────────────────────────────────
counts_raw <- read.table("counts.txt",
                         header = TRUE,
                         skip   = 1,
                         sep    = "\t")

counts_raw <- counts_raw[, c(-2, -3, -4, -5, -6)]
colnames(counts_raw) <- c("Geneid", "GLU1", "GLU2", "SAC1", "SAC2")
rownames(counts_raw) <- counts_raw$Geneid
counts_raw <- counts_raw %>% select(-Geneid)

cat("Total genes loaded:", nrow(counts_raw), "\n")

# ── 2. edgeR - filter, normalise ──────────────────────────────────────────────
group    <- factor(c("GLU", "GLU", "SAC", "SAC"), levels = c("GLU", "SAC"))
dge      <- DGEList(counts = counts_raw, group = group)
keep     <- filterByExpr(dge)
dge      <- dge[keep, , keep.lib.sizes = FALSE]
dge_norm <- calcNormFactors(dge, method = "TMM")

cat("After low-count filtering:", nrow(dge_norm$counts), "genes retained\n")

# ── 3. PCA ────────────────────────────────────────────────────────────────────
dge_logcpm        <- cpm(dge_norm, log = TRUE)
pca               <- prcomp(t(dge_logcpm), scale. = FALSE)
pca_frame         <- as.data.frame(pca$x)
pca_frame$Sample  <- rownames(pca_frame)
pca_frame$Group   <- as.character(group)
var_explained     <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)

pca_plot <- ggplot(pca_frame, aes(x = PC1, y = PC2, fill = Group)) +
  geom_point(size = 5, shape = 21, color = "black", stroke = 0.6) +
  geom_text_repel(
    aes(label = Sample),
    size = 4, fontface = "bold", color = "black",
    box.padding = 0.6, point.padding = 0.3, seed = 1, segment.size = 0
  ) +
  scale_fill_manual(name = "Condition",
                    values = c("GLU" = "#3B7DA8", "SAC" = "#D9764A")) +
  xlab(paste0("PC1 (", var_explained[1], "%)")) +
  ylab(paste0("PC2 (", var_explained[2], "%)")) +
  labs(title    = "Sample Clustering by Condition",
       subtitle = "GLU replicates separate on PC2 - within-group variance noted (see METHODS)") +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.title    = element_text(face = "bold", size = 12),
    axis.text     = element_text(color = "black", size = 10),
    legend.title  = element_text(face = "bold", size = 11),
    legend.text   = element_text(size = 11),
    plot.margin   = margin(t = 20, r = 20, b = 10, l = 10)
  )

ggsave("PCA_plot_publication.png", pca_plot,
       width = 6, height = 5, dpi = 600, bg = "white")
cat("Saved: PCA_plot_publication.png\n")

# ── 4. Differential expression ────────────────────────────────────────────────
design   <- model.matrix(~ group)           # coef 2 = SAC - GLU
dge_disp <- estimateDisp(dge_norm, design)
fit      <- glmQLFit(dge_disp, design)
res      <- glmQLFTest(fit, coef = 2)

de_table         <- topTags(res, n = Inf)$table
de_table$GeneID  <- rownames(de_table)

de_table$Significant <- "Not Significant"
de_table$Significant[de_table$logFC >  1 & de_table$FDR < 0.05] <- "Upregulated in SAC"
de_table$Significant[de_table$logFC < -1 & de_table$FDR < 0.05] <- "Downregulated in SAC"

cat("\nDifferential expression summary:\n")
cat("  Upregulated in SAC:  ", sum(de_table$Significant == "Upregulated in SAC"), "\n")
cat("  Downregulated in SAC:", sum(de_table$Significant == "Downregulated in SAC"), "\n")
cat("  Not significant:     ", sum(de_table$Significant == "Not Significant"), "\n")

# ── 5. Gene name annotation via keggGet (Sakai-specific, batched + cached) ────
# org.EcK12.eg.db is NOT used - see header note.
# keggList("ecs") only returns "CDS" for most genes - not useful.
# keggGet() returns real descriptions in $NAME but is limited to 10 IDs/call
# and is slow (~500+ calls for the full gene set), so results are cached.

lookup_path <- "ecs_gene_names_lookup.rds"

if (file.exists(lookup_path)) {

  cat("\nLoading cached KEGG gene name lookup from", lookup_path, "\n")
  gene_names_df <- readRDS(lookup_path)

} else {

  cat("\nNo cached lookup found - fetching gene names from KEGG via keggGet()\n")
  cat("This queries KEGG in batches of 10 and will take several minutes ")
  cat("for the full gene set. Result is cached so this only runs once.\n")

  all_ids <- unique(paste0("ecs:", de_table$GeneID))
  batches <- split(all_ids, ceiling(seq_along(all_ids) / 10))

  cat("Total genes to annotate:", length(all_ids),
      "in", length(batches), "batches\n")

  gene_names_list <- vector("list", length(batches))

  for (i in seq_along(batches)) {
    tryCatch({
      entries   <- keggGet(batches[[i]])
      names_vec <- sapply(entries, function(x) {
        if (!is.null(x$NAME)) x$NAME[1] else NA
      })
      id_vec <- batches[[i]][seq_along(names_vec)]

      gene_names_list[[i]] <- data.frame(
        kegg_id  = id_vec,
        GeneName = names_vec,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      cat("  Batch", i, "failed:", conditionMessage(e), "\n")
    })

    if (i %% 50 == 0) cat("  Completed batch", i, "of", length(batches), "\n")
    Sys.sleep(0.1)   # be polite to KEGG's server
  }

  gene_names_df <- bind_rows(gene_names_list)
  gene_names_df$GeneID   <- gsub("^ecs:", "", gene_names_df$kegg_id)
  gene_names_df$GeneName <- gsub("^\\(RefSeq\\)\\s*", "", gene_names_df$GeneName)

  saveRDS(gene_names_df, lookup_path)
  write.csv(gene_names_df, "ecs_gene_names_lookup.csv", row.names = FALSE)

  cat("Annotation complete and cached:", sum(!is.na(gene_names_df$GeneName)),
      "of", nrow(gene_names_df), "genes named\n")
}

# Join gene names onto DE table
de_table <- de_table %>%
  left_join(gene_names_df %>% select(GeneID, GeneName), by = "GeneID")

# Fall back to locus tag if no KEGG name found
de_table$GeneName <- ifelse(is.na(de_table$GeneName) | de_table$GeneName == "",
                            de_table$GeneID, de_table$GeneName)

cat("Gene names attached for",
    sum(de_table$GeneName != de_table$GeneID), "of", nrow(de_table), "genes\n")

write.csv(de_table, "DEresults_SAC_vs_GLU.csv", row.names = FALSE)
cat("Saved: DEresults_SAC_vs_GLU.csv\n")

# ── 6. Volcano plot ───────────────────────────────────────────────────────────
top_label <- de_table %>%
  filter(Significant != "Not Significant") %>%
  arrange(desc(abs(logFC))) %>%
  head(12)

volcano <- ggplot(de_table,
                  aes(x = logFC, y = -log10(FDR), color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = top_label, aes(label = GeneName),
    size = 2.6, color = "black", box.padding = 0.4,
    max.overlaps = 20, segment.size = 0.3
  ) +
  scale_color_manual(values = c(
    "Upregulated in SAC"   = "#D9764A",
    "Downregulated in SAC" = "#3B7DA8",
    "Not Significant"      = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.4) +
  labs(title = "Volcano Plot: SAC vs GLU",
       x = "Log2 Fold Change (SAC / GLU)", y = "-Log10 FDR", color = "Expression") +
  theme_bw(base_size = 13) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    plot.margin  = margin(t = 20, r = 20, b = 10, l = 10)
  )

ggsave("Volcano_plot.png", volcano, width = 8, height = 6, dpi = 600, bg = "white")
cat("Saved: Volcano_plot.png\n")

# ── 7. KEGG pathway gene sets via KEGGREST ────────────────────────────────────
cat("\nFetching KEGG pathway annotations for ecs (Sakai)...\n")

ecs_pathways   <- keggLink("pathway", "ecs")
ecs_path_names <- keggList("pathway/ecs")

cat("  Pathway-gene links:", length(ecs_pathways), "\n")
cat("  Pathways available:", length(ecs_path_names), "\n")

ecs_df       <- data.frame(
  pathway = as.character(ecs_pathways),
  gene    = names(ecs_pathways),
  stringsAsFactors = FALSE
)
pathway_list <- split(ecs_df$gene, ecs_df$pathway)

clean_names <- gsub(" - Escherichia coli O157:H7 Sakai \\(EHEC\\)", "", ecs_path_names)
names(pathway_list) <- clean_names[gsub("path:", "", names(pathway_list))]

# ── 8. Ranked gene list for fgsea ─────────────────────────────────────────────
gene_list_fgsea <- -de_table$logFC    # negated - see header note
names(gene_list_fgsea) <- paste0("ecs:", de_table$GeneID)
gene_list_fgsea <- sort(gene_list_fgsea, decreasing = TRUE)

cat("  Genes in ranked list:          ", length(gene_list_fgsea), "\n")
cat("  Genes overlapping KEGG pathways:",
    sum(names(gene_list_fgsea) %in% unlist(pathway_list)), "\n")

# ── 9. fgsea ───────────────────────────────────────────────────────────────────
set.seed(42)
fgsea_res <- fgsea(
  pathways = pathway_list,
  stats    = gene_list_fgsea,
  minSize  = 15,
  maxSize  = 500
)

fgsea_res           <- as.data.frame(fgsea_res)
fgsea_res$Direction <- ifelse(fgsea_res$NES > 0, "Upregulated in SAC", "Downregulated in SAC")
fgsea_res           <- fgsea_res[order(fgsea_res$NES, decreasing = TRUE), ]

write.csv(fgsea_res %>% select(-leadingEdge), "GSEA_fgsea_results.csv", row.names = FALSE)
cat("Saved: GSEA_fgsea_results.csv\n")

fgsea_sig <- fgsea_res[fgsea_res$padj < 0.05, ]
cat("\nSignificant pathways (padj < 0.05):", nrow(fgsea_sig), "\n")
print(fgsea_sig[, c("pathway", "NES", "padj", "size", "Direction")])

# ── 10. GSEA combined dot plot ────────────────────────────────────────────────
top_up   <- fgsea_sig %>% filter(Direction == "Upregulated in SAC")   %>% arrange(desc(NES)) %>% head(10)
top_down <- fgsea_sig %>% filter(Direction == "Downregulated in SAC") %>% arrange(NES)        %>% head(10)

plot_data <- bind_rows(top_up, top_down)
plot_data$pathway <- factor(plot_data$pathway, levels = plot_data$pathway[order(plot_data$NES)])

plot_data$pathway <- recode(plot_data$pathway,
  "Cationic antimicrobial peptide (CAMP) resistance" = "CAMP resistance",
  "Alanine, aspartate and glutamate metabolism"       = "Ala/Asp/Glu metabolism",
  "Biofilm formation - Escherichia coli"              = "Biofilm formation (E. coli)"
)

gsea_dotplot <- ggplot(plot_data, aes(x = NES, y = pathway, size = size, color = Direction)) +
  geom_point(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_color_manual(values = c(
    "Upregulated in SAC"   = "#D9764A",
    "Downregulated in SAC" = "#3B7DA8"
  )) +
  scale_size_continuous(name = "Gene count", range = c(3, 10)) +
  labs(title = "GSEA Pathway Enrichment - SAC vs GLU",
       subtitle = "fgsea | KEGG Sakai (ecs) | padj < 0.05 | Top 10 per direction",
       x = "Normalized Enrichment Score (NES)", y = NULL, color = "Direction") +
  theme_bw(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.y        = element_text(size = 10, color = "black"),
    axis.text.x        = element_text(size = 10, color = "black"),
    legend.position    = "right",
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.major.x = element_line(color = "grey90"),
    plot.margin        = margin(t = 20, r = 20, b = 10, l = 20)
  )

ggsave("GSEA_combined_dotplot.png", gsea_dotplot, width = 11, height = 7, dpi = 600, bg = "white")
cat("Saved: GSEA_combined_dotplot.png\n")

# ── 11. Leading edge gene extraction (with gene names) ────────────────────────
gene_info <- de_table %>%
  select(GeneID, GeneName, logFC, FDR, Significant) %>%
  mutate(kegg_id = paste0("ecs:", GeneID))

leading_edge_list <- lapply(seq_len(nrow(fgsea_sig)), function(i) {
  le <- unlist(fgsea_sig$leadingEdge[[i]])
  if (is.null(le) || length(le) == 0) return(NULL)

  gene_info[gene_info$kegg_id %in% le, ] %>%
    arrange(desc(logFC)) %>%
    mutate(pathway = fgsea_sig$pathway[i], NES = fgsea_sig$NES[i], padj = fgsea_sig$padj[i])
})

leading_edge_df <- bind_rows(leading_edge_list)

write.csv(
  leading_edge_df %>% select(pathway, NES, padj, GeneID, GeneName, logFC, FDR, Significant),
  "GSEA_leading_edge_genes.csv", row.names = FALSE
)
cat("Saved: GSEA_leading_edge_genes.csv\n")
cat("Total leading edge gene-pathway pairs:", nrow(leading_edge_df), "\n")

# ── 12. Top 5 driver genes per pathway ────────────────────────────────────────
top5_per_pathway <- leading_edge_df %>%
  group_by(pathway) %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n = 5) %>%
  select(pathway, GeneID, GeneName, logFC, FDR, Significant) %>%
  ungroup()

write.csv(top5_per_pathway, "Top5_genes_per_pathway.csv", row.names = FALSE)
cat("Saved: Top5_genes_per_pathway.csv\n")
cat("\nTop 5 driver genes per significant pathway:\n")
print(top5_per_pathway, n = 100)

# ── 13. Leading edge heatmap (uses real gene names as row labels) ─────────────
top_pathways <- fgsea_sig %>%
  filter(!grepl("Metabolic pathways|Biosynthesis of secondary|Microbial metabolism", pathway)) %>%
  arrange(padj) %>% head(8) %>% pull(pathway)

heatmap_long <- leading_edge_df %>%
  filter(pathway %in% top_pathways) %>%
  select(GeneID, GeneName, pathway, logFC)

shared_genes <- heatmap_long %>%
  group_by(GeneID) %>% summarise(n = n_distinct(pathway)) %>%
  filter(n >= 2) %>% pull(GeneID)

if (length(shared_genes) < 5) {
  shared_genes <- heatmap_long %>%
    group_by(GeneID) %>% summarise(n = n_distinct(pathway)) %>%
    filter(n >= 1) %>% pull(GeneID)
  cat("Note: using all leading edge genes (fewer than 5 shared across 2+ pathways)\n")
}

heatmap_long <- heatmap_long %>% filter(GeneID %in% shared_genes)

cat("\nGenes in heatmap:", n_distinct(heatmap_long$GeneID), "\n")
cat("Pathways in heatmap:\n")
cat(paste0("  - ", top_pathways), sep = "\n")

# Build readable, de-duplicated row labels from real gene names
heatmap_long_unique <- heatmap_long %>%
  distinct(GeneID, GeneName) %>%
  mutate(RowLabel = ifelse(duplicated(GeneName) | GeneName == "" | is.na(GeneName),
                           paste0(GeneName, " (", GeneID, ")"),
                           GeneName))

heatmap_long <- heatmap_long %>%
  left_join(heatmap_long_unique %>% select(GeneID, RowLabel), by = "GeneID")

heatmap_matrix <- heatmap_long %>%
  select(RowLabel, pathway, logFC) %>%
  pivot_wider(names_from = pathway, values_from = logFC) %>%
  column_to_rownames("RowLabel") %>%
  as.matrix()

heatmap_matrix[is.na(heatmap_matrix)] <- 0

colnames(heatmap_matrix) <- gsub("Cationic antimicrobial peptide \\(CAMP\\) resistance",
                                  "CAMP resistance", colnames(heatmap_matrix))
colnames(heatmap_matrix) <- gsub("Alanine, aspartate and glutamate metabolism",
                                  "Ala/Asp/Glu metabolism", colnames(heatmap_matrix))
colnames(heatmap_matrix) <- gsub("Biofilm formation - Escherichia coli",
                                  "Biofilm (E. coli)", colnames(heatmap_matrix))

gene_annot <- leading_edge_df %>%
  filter(GeneID %in% heatmap_long$GeneID) %>%
  select(GeneID, Significant) %>%
  distinct(GeneID, .keep_all = TRUE) %>%
  left_join(heatmap_long_unique %>% select(GeneID, RowLabel), by = "GeneID") %>%
  select(RowLabel, Significant) %>%
  column_to_rownames("RowLabel")

colnames(gene_annot) <- "Expression"

ann_colors <- list(Expression = c(
  "Upregulated in SAC"   = "#D9764A",
  "Downregulated in SAC" = "#3B7DA8",
  "Not Significant"      = "grey80"
))

max_val <- max(abs(heatmap_matrix), na.rm = TRUE)
breaks  <- seq(-max_val, max_val, length.out = 101)

png("GSEA_leadingedge_heatmap.png", width = 3200, height = 3800, res = 300, bg = "white")

pheatmap(
  heatmap_matrix,
  annotation_row    = gene_annot,
  annotation_colors = ann_colors,
  color             = colorRampPalette(c("#3B7DA8", "white", "#D9764A"))(100),
  breaks            = breaks,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 7,
  fontsize_col      = 10,
  border_color      = NA,
  main              = "Leading Edge Genes Across Significant Pathways\nSAC vs GLU | Log2 Fold Change",
  angle_col         = 45
)

dev.off()
cat("Saved: GSEA_leadingedge_heatmap.png\n")
cat("Heatmap row labels use real KEGG gene names (e.g. 'flagellin', 'flagellar hook")
cat(" assembly protein') instead of locus tags, falling back to 'Name (LocusTag)'")
cat(" only when a name is duplicated or missing.\n")

# ── 14. Final summary ──────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Output files produced:\n")
cat("  ecs_gene_names_lookup.rds         cached KEGG gene name lookup (reusable)\n")
cat("  ecs_gene_names_lookup.csv         same lookup, human-readable\n")
cat("  DEresults_SAC_vs_GLU.csv          full DE results with gene names\n")
cat("  GSEA_fgsea_results.csv            all pathway enrichment results\n")
cat("  GSEA_leading_edge_genes.csv       gene-level detail per pathway\n")
cat("  Top5_genes_per_pathway.csv        top 5 driver genes per pathway\n")
cat("  PCA_plot_publication.png          PCA figure (600 dpi)\n")
cat("  Volcano_plot.png                  volcano figure, real gene names (600 dpi)\n")
cat("  GSEA_combined_dotplot.png         GSEA dot plot up+down (600 dpi)\n")
cat("  GSEA_leadingedge_heatmap.png      heatmap with real gene names (300 dpi)\n")
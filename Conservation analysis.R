############################################################
## STQ/STP Conservation & GO Analysis
## Figures: 2D, 2E, S2B–F
## Author: <Francesca Conte>
## Date: <2025-12-01>
############################################################
# -----------------------------------------------------
# LOAD PACKAGES
# -----------------------------------------------------
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(org.Hs.eg.db)
library(purrr)
library(stringr)
library(tidyr)

# -----------------------------------------------------
# SET YOUR WORKING DIRECTORY

# -----------------------------------------------------
# -----------------------------------------------------
# LOAD DATA
# -----------------------------------------------------
phospho    <- read.delim("Phospho_results.txt", stringsAsFactors = FALSE) #Supplementary dataset 1
peptools   <- read.delim("Peptools_results.txt", stringsAsFactors = FALSE) #Supplementary dataset 2

# -----------------------------------------------------
# CREATE REGULATION COLUMNS
# -----------------------------------------------------
treatments <- c("CPT", "ETO", "Xray", "HU", "APH", "Gem",
                "H2O2", "FA", "MMS", "UV", "AsO2")

fc_cutoff <- log2(1.5)
p_cutoff  <- 0.05

# Loop over treatments and create Reg columns
for (tr in treatments) {
  
  # Correct column names
  logFC_col <- paste0("mean.log2.FC.", tr, "_vs_UT")
  pval_col  <- paste0(tr, "_vs_UT_adj.P.Val")
  reg_col   <- paste0("Reg", tr)
  
  # Safety check (prevents silent errors if a column is missing)
  if (!all(c(logFC_col, pval_col) %in% colnames(phospho))) {
    warning("Missing columns for ", tr)
    next
  }
  
  phospho[[reg_col]] <- ifelse(
    phospho[[logFC_col]] >=  fc_cutoff & phospho[[pval_col]] <= p_cutoff,  1,
    ifelse(
      phospho[[logFC_col]] <= -fc_cutoff & phospho[[pval_col]] <= p_cutoff, -1,
      0
    )
  )
}

# -----------------------------------------------------
# MERGE PHOSPHO + PEPTOOLS
# -----------------------------------------------------
phospho <- phospho %>%
  mutate(
    ID_position = paste(Uniprot_ID, paste(Aminoacid, Position, sep = ""), sep = "_")
  )

Conservation_data <- inner_join(
  phospho,
  peptools,
  by = "Phosphosite"
)

# -----------------------------------------------------
# DEFINE REGULATED SITES
# -----------------------------------------------------

reg_cols <- paste0("Reg", treatments)

Conservation_data <- Conservation_data %>%
  mutate(
    Regulated_full = if_any(all_of(reg_cols), ~ .x %in% c(1, -1))
  )

# -----------------------------------------------------
# HELPER: BUILD VIOLIN DATASET
# -----------------------------------------------------
build_violin_df <- function(df, score_col) {
  
  df %>%
    filter(!is.na(.data[[score_col]]) & .data[[score_col]] > 0) %>%
    mutate(
      Group = case_when(
        STPXRK.site & Regulated_full          ~ "STPXRK_regulated",
        STPXRK.site & !Regulated_full         ~ "STPXRK_nonregulated",
        STP.site & !STPXRK.site & Regulated_full  ~ "STP_regulated",
        STP.site & !STPXRK.site & !Regulated_full ~ "STP_nonregulated",
        STQ.site & Regulated_full             ~ "STQ_regulated",
        STQ.site & !Regulated_full            ~ "STQ_nonregulated",
        TRUE                                  ~ NA_character_
      ),
      logValue = pmin(-log10(.data[[score_col]] + 1e-6), 6)
    ) %>%
    filter(!is.na(Group))
}

# -----------------------------------------------------
# HELPER: DRAW VIOLIN PLOT
# -----------------------------------------------------
draw_violin <- function(df, score_label, outfile) {
  
  group_order <- c(
    "STQ_regulated", "STQ_nonregulated",
    "STP_regulated", "STP_nonregulated",
    "STPXRK_regulated", "STPXRK_nonregulated"
  )
  
  sample_sizes <- df %>%
    count(Group) %>%
    mutate(
      label = paste0(Group, " (n = ", n, ")"),
      Group = factor(Group, levels = group_order)
    )
  
  df <- left_join(df, sample_sizes, by = "Group") %>%
    mutate(label = factor(label, levels = sample_sizes$label))
  
  comparisons <- combn(levels(df$label), 2, simplify = FALSE)
  
  p <- ggplot(df, aes(x = label, y = logValue)) +
    geom_violin(trim = FALSE, fill = "lightgray", alpha = 0.4) +
    geom_boxplot(width = 0.3, outlier.shape = NA, fill = "white") +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = comparisons,
      p.adjust.method = "BH",
      label = "p.signif"
    ) +
    theme_classic() +
    coord_cartesian(ylim = c(0,5)) +
    labs(
      x = "Group",
      y = paste0("-log10(", score_label, ")"),
      title = score_label
    )
  
  ggsave(outfile, p, width = 9, height = 6)
  
  p
}

# -----------------------------------------------------
# METAZOA VIOLIN (Figure S2B)
# -----------------------------------------------------
v_metazoa <- build_violin_df(
  Conservation_data,
  "RLC_score_metazoa"
)

draw_violin(
  v_metazoa,
  "Metazoa Conservation Score",
  "Violinplot_Metazoa_Conservation_STQ_STP_STPXRK_log10.pdf"
)

# -----------------------------------------------------
# MAMMALIA VIOLIN (Figure S2C)
# -----------------------------------------------------
v_mammalia <- build_violin_df(
  Conservation_data,
  "RLC_score_mammalia"
)

draw_violin(
  v_mammalia,
  "Mammalia Conservation Score",
  "Violinplot_Mammalia_Conservation_STQ_STP_STPXRK_log10.pdf"
)

# -----------------------------------------------------
# GO HELPER FUNCTIONS
# -----------------------------------------------------
run_go <- function(gene_list, universe) {
  
  ego <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    universe = universe,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  ego <- clusterProfiler::simplify(ego, cutoff = 0.6)
  
  as.data.frame(ego) %>%
    mutate(FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
}

# -----------------------------------------------------
# DEFINE UNIVERSES
# -----------------------------------------------------
universe_conservation      <- unique(Conservation_data$Gene_name)
universe_phosphoproteome  <- unique(phospho$Gene_name)

# -----------------------------------------------------
# DEFINE INDUCED SITES - EXCLUDE ASO2-SPECIFIC SITES FOR THIS ANALYSIS
# -----------------------------------------------------
induced_sites <- Conservation_data %>%
  filter(
    RegCPT == 1 | RegETO == 1 | RegUV == 1 | RegMMS == 1 | RegFA == 1 |
      RegH2O2 == 1 | RegHU == 1 | RegGem == 1 | RegAPH == 1 | RegXray == 1
  )


# -----------------------------------------------------
# GO: STQ INDUCED (Figure 2D)
# -----------------------------------------------------
only_STQ_genes <- induced_sites %>%
  filter(STQ.site) %>%
  pull(Gene_name) %>%
  unique()

go_STQ <- run_go(only_STQ_genes, universe_conservation)

write.table(
  go_STQ,
  "GO ALL induced STQ sites.txt",
  sep = "\t",
  row.names = FALSE
)

# -----------------------------------------------------
# GO BARPLOT FUNCTION
# -----------------------------------------------------
plot_go_bar <- function(df, n_terms, out_file) {
  
  top_terms <- df %>%
    group_by(ONTOLOGY) %>%
    slice_max(FoldEnrichment, n = n_terms) %>%
    ungroup() %>%
    arrange(ONTOLOGY, FoldEnrichment)
  
  p <- ggplot(top_terms,
              aes(FoldEnrichment,
                  reorder(Description, FoldEnrichment),
                  fill = ONTOLOGY)) +
    geom_col() +
    geom_text(aes(label = sprintf("q = %.2g", qvalue)),
              hjust = -0.1,
              size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      x = "Fold Enrichment",
      y = NULL,
      fill = "Ontology"
    ) +
    theme_classic()
  
  ggsave(out_file, p, width = 8, height = 5)
  p
}

plot_go_bar(go_STQ, 2,
            "STQ induced sites Top2_GO_terms_by_Ontology_barplot.pdf")

# -----------------------------------------------------
# GO: STP INDUCED (no CDK targets) (Figure 2D)
# -----------------------------------------------------

only_STP_genes <- induced_sites %>%
  filter(STP.site, STPXRK.site == FALSE) %>%
  pull(Gene_name) %>%
  unique()

go_STP <- run_go(only_STP_genes, universe_conservation)

write.table(
  go_STP,
  "GO ALL STP no CDK targets induced sites.txt",
  sep = "\t",
  row.names = FALSE
)

plot_go_bar(
  go_STP,
  n_terms = 6,
  out_file = "STP no CDK targets induced sites Top6_GO_terms_barplot.pdf"
)

#######################################################################

# -----------------------------------------------------
# FILTER FOR VALID CONSERVATION SCORES
# -----------------------------------------------------
Conservation_filtered <- Conservation_data %>%
  filter(!is.na(RLC_score_metazoa) & RLC_score_metazoa > 0)

# -----------------------------------------------------
# DEFINE INDUCED SITES (include all treatments, explicitly RegAsO2)
# -----------------------------------------------------
all_induced_sites <- Conservation_filtered %>%
  filter(
    RegCPT  == 1 | RegETO  == 1 | RegUV  == 1 | RegMMS  == 1 | 
      RegFA   == 1 | RegH2O2 == 1 | RegHU  == 1 | RegGem  == 1 |
      RegAPH  == 1 | RegXray == 1 | RegAsO2 == 1
  )

# -----------------------------------------------------
# CONSERVATION STRATIFICATION
# -----------------------------------------------------
get_genes_by_range <- function(df, min = -Inf, max = Inf) {
  df %>%
    filter(STQ.site,
           RLC_score_metazoa >= min,
           RLC_score_metazoa < max) %>%
    pull(Gene_name) %>%
    unique()
}

genes_low    <- get_genes_by_range(all_induced_sites, 0.65, 1)
genes_medium <- get_genes_by_range(all_induced_sites, 0.35, 0.65)
genes_high   <- get_genes_by_range(all_induced_sites, 0,    0.35)

GO_list <- list(
  "Low conservation"    = genes_low,
  "Medium conservation" = genes_medium,
  "High conservation"   = genes_high
)

# -----------------------------------------------------
# CALCULATE PERCENTAGES OF INDUCED STQ SITES BY CONSERVATION LEVEL (Figure S2D)
# -----------------------------------------------------

# Filter Conservation_filtered to keep only significantly regulated phosphosites
regulated_sites <- subset(Conservation_filtered, 
                          RegCPT == 1 | RegCPT == -1 | RegETO == 1 | RegETO == -1 |
                            RegUV  == 1 | RegUV  == -1 | RegAsO2 == 1 | RegAsO2 == -1 |
                            RegMMS == 1 | RegMMS == -1 | RegFA  == 1 | RegFA  == -1 |
                            RegH2O2 == 1 | RegH2O2 == -1 | RegHU == 1 | RegHU == -1 |
                            RegGem == 1 | RegGem == -1 | RegAPH == 1 | RegAPH == -1 |
                            RegXray == 1 | RegXray == -1
)

regulated_STQ_sites <- regulated_sites %>%
  filter(STQ.site == TRUE)

# Calculate percentages relative to total regulated STQ sites
percentages <- data.frame(
  Conservation = c("High conservation", "Medium conservation", "Low conservation"),
  Percentage = c(
    nrow(high_conservation_induced_STQ) / nrow(regulated_STQ_sites) * 100,
    nrow(medium_conservation_induced_STQ) / nrow(regulated_STQ_sites) * 100,
    nrow(low_conservation_induced_STQ) / nrow(regulated_STQ_sites) * 100
  )
)

# Order by descending percentage
percentages$Conservation <- factor(
  percentages$Conservation,
  levels = percentages$Conservation[order(percentages$Percentage, decreasing = TRUE)]
)

# Plot horizontal barplot
FP_reg_plot <- ggplot(percentages, aes(x = Conservation, y = Percentage, fill = Conservation)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  coord_flip() +
  geom_text(aes(label = sprintf("%.2f%%", Percentage)), hjust = -0.1, size = 4.5) +
  scale_fill_manual(values = c("lightgrey", "grey60", "grey40")) +
  labs(
    x = NULL,
    y = "Percentage of induced STQ sites\n(relative to total regulated STQ sites)",
    title = "Induced STQ sites by conservation level"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
  )

ggsave("Barplot induced STQ sites stratified by conservation (percentage calculated against all regulated STQ sites).pdf", 
       plot = FP_reg_plot, width = 10, height = 10)

# -----------------------------------------------------
# GO COMPARE CLUSTER (Figure S2E)
# simplify function not working for compareCluster in new clusterProfiler versions, bug needs to be fixed
# -----------------------------------------------------
compare_clusters <- compareCluster(
  GO_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  universe = universe_phosphoproteome,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

compare_clusters <- simplify(compare_clusters)

go_comp_df <- as.data.frame(compare_clusters) %>%
  mutate(FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio)) %>%
  filter(Count >= 4)

top_terms <- go_comp_df %>%
  group_by(Cluster) %>%
  slice_max(FoldEnrichment, n = 8) %>%
  ungroup()

# -----------------------------------------------------
# DOTPLOT
# -----------------------------------------------------
p <- ggplot(top_terms,
            aes(Description,
                Cluster,
                size = Count,
                color = p.adjust)) +
  geom_point() +
  scale_color_gradientn(
    colours = c("#0c2c84", "#41b6c4", "#c7e9b4"),
    trans = "log10"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "GO Terms",
    y = "Condition",
    color = "Adjusted p-value",
    size = "Gene count"
  )

ggsave("GO STQ induced sites stratified by Metazoa conservation.pdf",
       p, width = 8, height = 8)

write.table(
  go_comp_df,
  "GO STQ induced sites stratified by Metazoa conservation.txt",
  sep = "\t",
  row.names = FALSE
)


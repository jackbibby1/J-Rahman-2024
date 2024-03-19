# load in packages --------------------------------------------------------
library(oligo)
library(limma)
library(tidyverse)
library(plotting)
library(magrittr)
library(pd.huex.1.0.st.v2)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)

# set global parameters ---------------------------------------------------
setwd("dra_ptgs2/")

`%notin%` <- Negate(`%in%`)

pt_genes <- c("FOXP3", "IFNG", "IL10", "PTGIS", "PTGIR", "PTGS2", "PTGES2", "IL1R2")

annotate_data <- function(id_data) {
  
  annotated_data <- id_data %>%
    data.frame() %>%
    rownames_to_column("id") %>%
    full_join(annots, "id") %>%
    dplyr::select(id, SYMBOL, everything()) %>%
    rename(gene = SYMBOL)
  
  return(annotated_data)
  
}

# read in data ------------------------------------------------------------
files <- list.files(path = "raw_files", 
                    pattern = ".CEL", full.names = T)
files <- grep(pattern = "P1", x = files, value = T, invert = T)
files <- grep(pattern = "donor4", x = files, value = T, invert = T)
raw_data <- read.celfiles(files)

# create metadata ---------------------------------------------------------
meta <- data.frame(file_path = files) %>%
  mutate(donor = str_extract(files, pattern = "donor[0-3]"),
         antagonist = str_extract(files, "plus|minus"),
         stimulation = str_extract(files, "NA|3-46"),
         stimulation = ifelse(stimulation == "NA", "unstim", "stim"),
         stim_ant = paste(stimulation, antagonist, sep = "_"))

# annotate data -----------------------------------------------------------
pData(raw_data)$stimulation <- meta$stimulation
pData(raw_data)$donor <- meta$donor
pData(raw_data)$antagonist <- meta$antagonist
pData(raw_data)$sample <- sampleNames(raw_data)
pData(raw_data)$stim_ant <- meta$stim_ant

# data normalisation ------------------------------------------------------
norm_data <- rma(raw_data, target = "core")


# data filtering ----------------------------------------------------------
hist(rowMedians(exprs(norm_data)), 
     breaks = 50, 
     col = "gray90")
abline(v = 4, col = "red", lwd = 2)

# 4 seems like a reasonable cutoff

keep_genes <- apply(exprs(norm_data), 1, function(x) sum(x > 4) >= 3)
table(keep_genes)

norm_data <- norm_data[keep_genes, ]

# and now looking at data distribution after filtering 

hist(rowMedians(exprs(norm_data)), 
     breaks = 50, 
     col = "gray90")
abline(v = 4, col = "red", lwd = 2)

# data visualisation after processing -------------------------------------
par(mfrow = c(1, 2))
boxplot(raw_data, which = "all", las = 2, main = "No normalisation")
boxplot(norm_data, las = 2, main = "RMA normalisation")

# pca using top 5 pcs -----------------------------------------------------
df <- exprs(norm_data)

pcs <- paste0("PC", 2:5)
plots <- lapply(pcs, function(x) {
  draw_pca(df, grouping = meta$stim_ant, return_data = F, log = F, y_axis = x)
})
patchwork::wrap_plots(plots, ncol = 4) + patchwork::plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# plot nice pca -----------------------------------------------------------
pca_df <- prcomp(t(df))
stdev_val <- pca_df$sdev^2
stdev_val <- round(stdev_val/sum(stdev_val) * 100, 1)

data.frame(pc1 = pca_df$x[, "PC1"],
           pc3 = pca_df$x[, "PC3"],
           group = factor(meta$stim_ant, levels = c("unstim_minus", "unstim_plus", "stim_minus", "stim_plus"))) %>%
  ggplot(aes(-pc1, pc3)) +
  geom_point(shape = 21, size = 4, aes(fill = group), stroke = 0.3, alpha = 0.8) +
  scale_fill_manual(values = c("gray90", "#fc956d", "#6d93fc", "#69db63"),
                    name = "Condition",
                    labels = c("unstim_minus" = "Unstim -DRA", "unstim_plus" = "Unstim +DRA",
                               "stim_minus" = "Stim -DRA", "stim_plus" = "Stim +DRA")) +
  labs(x = paste0("PC1: ", stdev_val[1], "%"),
       y = paste0("PC3: ", stdev_val[3], "%")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1,
        legend.key = element_blank(),
        legend.key.height = unit(5, "mm"),
        legend.spacing.x = unit(0.1, "mm"),
        legend.spacing.y = unit(1, "mm"),
        legend.box.spacing = unit(0, "mm"))

ggsave("pca_pc3.png", dpi = 600, width = 3.7, height = 3.7)

# plot ptgs2 values -------------------------------------------------------
annots <- AnnotationDbi::select(huex10sttranscriptcluster.db::huex10sttranscriptcluster.db, 
                                keys = featureNames(norm_data), 
                                columns = "SYMBOL",
                                keytype = "PROBEID") %>%
  rename("id" = "PROBEID")

df %>%
  annotate_data() %>%
  filter(gene == "PTGS2") %>%
  pivot_longer(cols = X3.46_donor1_2hr.minus_01.CEL:NA_donor3_plus_13.CEL,
               names_to = "sample", values_to = "expression") %>%
  mutate(expression = 2^expression,
         group = factor(meta$stim_ant, levels = c("unstim_minus", "unstim_plus", "stim_minus", "stim_plus"))) %>%
  filter(group != "unstim_plus") %>%
  ggplot(aes(group, expression)) +
  geom_boxplot(aes(fill = group), alpha = 0.4, linewidth = 0.4) +
  geom_jitter(aes(fill = group), shape = 21, width = 0.2, size = 3.4, stroke = 0.4) +
  scale_fill_manual(values = c("#fc956d", "#6d93fc", "#69db63")) +
  scale_x_discrete(labels = c("unstim_minus" = "Unstim", "stim_minus" = "Stim", "stim_plus" = "Stim +\n DRA")) +
  labs(y = "Normalised expression") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("ptgs2_boxplot.png", dpi = 600, width = 2, height = 2.7)


# heatmap for pc3 loadings ------------------------------------------------
pc3_loadings <- pca_df$rotation %>%
  annotate_data() %>%
  select(gene, PC1, PC3) %>%
  filter(gene != "") %>%
  arrange(desc(abs(PC3))) %>%
  head(10)

keep_samples <- meta$stim_ant %in% c("unstim_minus", "stim_minus", "stim_plus")
heatmap_names <- c("Stim-DRA", "Stim+DRA", "Stim-DRA", "Stim+DRA", "Stim+DRA", "Stim-DRA",
                   "Unstim-DRA", "Unstim-DRA", "Unstim-DRA")

hm_col <- colorRamp2(breaks = c(-2, 0, 2), colors = c("white", "grey93", "#ff582e"))

hm_df <- df %>%
  annotate_data() %>%
  filter(gene %in% pc3_loadings$gene) %>%
  select(-id) %>%
  mutate(gene = make.names(gene, unique = T)) %>%
  column_to_rownames("gene")

pc3_hm <- hm_df[, keep_samples] %>%
  t() %>%
  scale() %>%
  t() %>%
  Heatmap(name = "z-score",
          border = T,
          row_dend_width = unit(3, "mm"),
          column_dend_height = unit(4, "mm"),
          column_title = NULL,
          col = hm_col,
          rect_gp = gpar(lwd = 0.2, col = "white"),
          column_labels = heatmap_names,
          column_km = 2,
          row_names_gp = gpar(fontsize = 10, fontface = "italic"),
          column_names_gp = gpar(fontsize = 11.5))

png("pc3_heatmap.png", width = 3.7, height = 3.8, res = 600, units = "in")
draw(pc3_hm)
dev.off()


# statistical analysis ----------------------------------------------------
groups <- model.matrix(~0 + meta$stim_ant) %>%
  set_colnames(c("stim_minus", "stim_plus", "unstim_minus", "unstim_plus"))

fit <- lmFit(norm_data, design = groups)

contrast_matrix <- makeContrasts(unstim_minus-stim_minus, 
                                 unstim_minus-unstim_plus,
                                 stim_minus-stim_plus,
                                 unstim_minus-stim_plus,
                                 levels = groups)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

results <- decideTests(fit2)
vennDiagram(results, 
            names = c("unstim - DRA vs \n stim - DRA",
                      "unstim - DRA vs \n unstim + DRA",
                      "stim - DRA vs \n stim + DRA",
                      "unstim - DRA vs \n stim + DRA"),
            circle.col = c("red", "blue", "gray", "green"))


# exporting the c5ar signature --------------------------------------------
unstim_minus_stim_plus_degs %>%
  filter(adj.P.Val < 0.01 & logFC > 1) %>%
  head(150) %>%
  write.csv("c5_enriched_signature.csv", row.names = F)

unstim_minus_stim_plus_degs %>%
  filter(adj.P.Val < 0.01 & logFC < -1) %>%
  head(150) %>%
  write.csv("c5_depleted_signature.csv", row.names = F)



# data for input into gsea ------------------------------------------------
norm_annotated %>%
  select(contains(c("gene", "3.46"))) %>%
  drop_na() %>%
  rename(NAME = gene) %>%
  mutate(DESCRIPTION = 1) %>%
  select(NAME, DESCRIPTION, everything()) %>%
  set_colnames(c("NAME", "DESCRIPTION", "minus_1", "plus_1", "minus_2", "plus_2", "plus_3", "minus_3")) %>%
  write.table(file = "gsea_in.txt", sep = "\t", quote = F, row.names = F)


# reading and plotting gsea output ----------------------------------------
gsea_output <- rbind(
  
  read.delim("stim_comparison.Gsea.1680712855192/gsea_report_for_minus_1680712855192.tsv") %>%
    janitor::clean_names() %>%
    select(name, nes, fdr_q_val) %>%
    mutate(enrichment = "minus"),
  read.delim("stim_comparison.Gsea.1680712855192/gsea_report_for_plus_1680712855192.tsv") %>%
    janitor::clean_names() %>%
    select(name, nes, fdr_q_val) %>%
    mutate(enrichment = "plus")
  
) %>%
  mutate(rank = percent_rank(-nes))

ggplot(gsea_output, aes(rank, abs(nes))) +
  geom_point(shape = 21, fill = "gray80", size = 2.5, stroke = 0.05, alpha = 0.9) +
  geom_point(data = subset(gsea_output, name == "WP_PROSTAGLANDIN_SIGNALING"),
             shape = 21, fill = "#ff6017", size = 4, stroke = 0.4) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(0.3, 3.5)) +
  labs(x = "Pathway rank", y = "Absolute NES") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA))

ggsave("gsea_output.png", dpi = 600, height = 2.6, width = 3)












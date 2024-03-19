# load in packages --------------------------------------------------------
library(tidyverse)
library(edgeR)
library(ComplexHeatmap)
library(magrittr)
library(janitor)
library(ggrepel)
library(SCPA)
library(msigdbr)
library(Seurat)
library(circlize)

# set global parameters ---------------------------------------------------
setwd("ptgis-r_sirna/")

# read in raw data --------------------------------------------------------
df <- read.delim("/Volumes/data/jubayer_sirna/featurecounts_output/counts.txt", skip = 1, sep = "\t") %>%
  clean_names() %>%
  select(-c(chr, start, end, strand, length)) %>%
  column_to_rownames("geneid")

df <- set_colnames(df, str_extract(string = colnames(df), pattern = "s[0-9]+"))

# define metadata ---------------------------------------------------------
meta <- read.csv("/Volumes/data/jubayer_sirna/sample_information.csv") %>%
  mutate(group = paste(activation, gene, population, sep = "_"))


# check some raw qc -------------------------------------------------------
par(mfrow = c(1, 2))
barplot(colSums(df)/1e6, las = 2, main = "Sequencing depth", ylab = "Counts (millions)")
boxplot(df, outline = F, las = 2, main = "Count distribution")
par(mfrow = c(1, 1))


# data filtering ----------------------------------------------------------
df <- select(df, meta$sample_id) # reorganise samples so they match the metadata
y <- DGEList(counts = df, group = meta$group)
keep_genes <- filterByExpr(y, group = meta$group)
y <- y[keep_genes, , keep.lib.sizes = F]

# data normalisation ------------------------------------------------------
y <- calcNormFactors(y)
norm_data <- cpm(y) %>%
  data.frame() %>%
  set_colnames(meta$group) %>%
  clean_names()

dir.create("expression_data", showWarnings = F)
write.csv(norm_data, "expression_data/tmm_cpm_normalised_data.csv")


# looking at general pca --------------------------------------------------
lapply(c(2, 3, 4), function(x) {
  
  plotting::draw_pca(df = norm_data + 1, 
                     log = T, y_axis = paste0("PC", x),
                     grouping = factor(meta$group, levels = c("unstim_scr_bulk", 
                                                              "stim_scr_dn", "stim_ptgis_dn", "stim_ptgir_dn",
                                                              "stim_scr_ifng", "stim_ptgir_ifng")), 
                     return_data = F) +
    scale_fill_manual(values = box_col) +
    ggrepel::geom_text_repel(label = meta$donor)
  
}) %>% patchwork::wrap_plots(guides = "collect")

# there's a big batch effect across donors


# comparing dn populations ------------------------------------------------
# subset to only keep dn populations
dn_meta <- filter(meta, population == "dn")
dn_df <- select(df, dn_meta$sample_id)

# define batches and groups to compare
donor_batch <- factor(dn_meta$donor)
group <- factor(dn_meta$group)
group <- relevel(group, ref = "stim_scr_dn")

# creatae edger object and filter
y <- DGEList(counts = dn_df, group = group)
keep_genes <- filterByExpr(y, group = group)
y <- y[keep_genes, , keep.lib.sizes = F]

# normalisation
y <- calcNormFactors(y)

# design matrix
design <- model.matrix(~0 + donor_batch + group)
rownames(design) <- colnames(y)
design

# estimate dispersion
y <- estimateDisp(y, design, robust = T)

# test
fit <- glmQLFit(y, design, robust = TRUE)
ptgir_degs <- glmQLFTest(fit, coef = 4) %>%
  topTags(n = Inf) %>%
  data.frame() %>%
  mutate(comparison = "DN: Scr vs PTGIR")

ptgis_degs <- glmQLFTest(fit, coef = 5) %>%
  topTags(n = Inf) %>%
  data.frame() %>%
  mutate(comparison = "DN: Scr vs PTGIS")

# comparing the ifng populations ------------------------------------------
# subset to only keep dn populations
ifng_meta <- filter(meta, population == "ifng")
ifng_df <- select(df, ifng_meta$sample_id)

# define batches and groups to compare
donor_batch <- factor(ifng_meta$donor)
group <- factor(ifng_meta$group)
group <- relevel(group, ref = "stim_scr_ifng")

# creatae edger object and filter
y <- DGEList(counts = ifng_df, group = group)
keep_genes <- filterByExpr(y, group = group)
y <- y[keep_genes, , keep.lib.sizes = F]

# normalisation
y <- calcNormFactors(y)

# design matrix
design <- model.matrix(~0 + donor_batch + group)
rownames(design) <- colnames(y)
design

# estimate dispersion
y <- estimateDisp(y, design, robust = T)

# test
fit <- glmQLFit(y, design, robust = TRUE)
ifng_degs <- glmQLFTest(fit, coef = 3) %>%
  topTags(n = Inf) %>%
  data.frame() %>%
  mutate(comparison = "IFNG: Scr vs PTGIR")


# volcano plots -----------------------------------------------------------
norm_data <- read.csv("ptgis-r_sirna/expression_data/tmm_cpm_normalised_data.csv")
ptgir_degs <- read.csv("ptgis-r_sirna/differential_expression/stim_scr_dn_vs_stim_ptgir_dn.csv") %>%
  rename("gene" = "X") %>%
  mutate(direction = case_when(FDR < 0.05 & logFC > 0.6 ~ "up",
                               FDR < 0.05 & logFC < -0.6 ~ "down",
                               FDR > 0.05 | abs(logFC) < 0.6 ~ "none"))

ptgis_degs <- read.csv("ptgis-r_sirna/differential_expression/stim_scr_dn_vs_stim_ptgis_dn.csv") %>%
  rename("gene" = "X") %>%
  mutate(direction = case_when(FDR < 0.05 & logFC > 0.6 ~ "up",
                               FDR < 0.05 & logFC < -0.6 ~ "down",
                               FDR > 0.05 | abs(logFC) < 0.6 ~ "none"))

degs <- list(ptgir_degs, ptgis_degs)

lapply(degs, function(x) {
  
  ggplot(x, aes(logFC, -log10(FDR))) +
    geom_point(data = subset(x, direction == "none"), shape = 21, size = 1, stroke = 0, fill = "gray80", alpha = 0.2) +
    geom_point(data = subset(x, direction == "up"), shape = 21, size = 2.5, fill = "#63a4ff", stroke = 0.1, alpha = 0.85) +
    geom_point(data = subset(x, direction == "down"), shape = 21, size = 2.5, fill = "#ff6363", stroke = 0.1, alpha = 0.85) +
    geom_point(data = subset(x, gene %in% c("PTGIS", "PTGIR")), 
               shape = 21, size = 3.2, color = "black", stroke = 0.3, fill = "#63d47e") +
    theme(panel.background = element_blank(),
          axis.line = element_line(),
          axis.title = element_text(size = 12))
  
}) %>% patchwork::wrap_plots()

ggsave("paper_figures/volcano_plot.png", dpi = 600, height = 2.8, width = 5.5)


# look at overlap of degs -------------------------------------------------
diff_ex <- list(ptgir_degs, ptgis_degs, ifng_degs)

de_results <- lapply(diff_ex, function(x) {
  
  up_genes <- subset(x, FDR < 0.05 & logFC > 0.6)
  down_genes <- subset(x, FDR < 0.05 & logFC < -0.6)
  not_sig <- subset(x, FDR > 0.05 | abs(logFC) < 0.6)
  
  data.frame(up = as.numeric(nrow(up_genes)),
             down = as.numeric(nrow(down_genes)),
             none = as.numeric(nrow(not_sig)))
  
})

degs <- lapply(diff_ex[1:2], function(x){
  
  filter(x, FDR < 0.05 & abs(logFC) > 0.6) %>%
    rownames()
  
})

names(degs) <- c("PTGIR", "PTGIS")

ggvenn::ggvenn(degs, 
               c("PTGIR", "PTGIS"), 
               set_name_size = 5, 
               fill_color = c("tomato", "cornflowerblue"), 
               stroke_size = 0.5)









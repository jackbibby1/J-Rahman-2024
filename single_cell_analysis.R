# load in packages --------------------------------------------------------
library(Seurat)
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(magrittr)
library(msigdbr)
library(SCPA)
library(data.table)
library(edgeR)
library(harmony)
library(circlize)
library(scPipeline)
`%notin%` <- Negate(`%in%`)

# set global parameters ---------------------------------------------------
# pt genes
c5_enriched <- read.csv("dra_ptgs2/c5_enriched_signature.csv")$gene
c5_enriched <- c5_enriched[!is.na(c5_enriched)]

c5_depleted <- read.csv("dra_ptgs2/c5_depleted_signature.csv")$gene
c5_depleted <- c5_depleted[!is.na(c5_depleted)]

c5_sig_scpa <- c(c5_enriched, c5_depleted)
c5_sig_scpa <- data.frame(Pathway = rep("C5", times = length(c5_sig_scpa)), 
                          Genes = c5_sig_scpa)

ptgir <- read.csv("~/My Drive/analysis_for_people/jubayer/ptgis-r_sirna/differential_expression/stim_scr_dn_vs_stim_ptgir_dn.csv", row.names = 1) %>%
  subset(FDR < 0.05 & abs(logFC) > 0.6) %>%
  rownames()

ptgis <- read.csv("~/My Drive/analysis_for_people/jubayer/ptgis-r_sirna/differential_expression/stim_scr_dn_vs_stim_ptgis_dn.csv", row.names = 1) %>%
  subset(FDR < 0.05 & abs(logFC) > 0.6) %>%
  rownames()

ptgir_ifng <- read.csv("~/My Drive/analysis_for_people/jubayer/ptgis-r_sirna/differential_expression/stim_scr_ifng_vs_stim_ptgir_ifng.csv", row.names = 1) %>%
  subset(FDR < 0.05 & abs(logFC) > 0.6) %>%
  rownames()

prost_sig <- data.frame(Pathway = "prost_sig",
                        Genes = c(ptgir_ifng, intersect(ptgir, ptgis)) %>% unique())

pathways <- c("hallmark", "kegg", "wiki")
cp_sets <- msigdbr("Homo sapiens") %>%
  filter(grepl(paste(pathways, collapse = "|"), gs_name, ignore.case = T)) %>%
  format_pathways()

gene_sets <- c(cp_sets, list(prost_sig), list(c5_sig_scpa))


# crohn's data ------------------------------------------------------------
files <- list.files("GSE157477_iel_crohns/raw_data", full.names = T)
df <- lapply(files, function(x) {
  
  df <- Read10X(x) %>%
    CreateSeuratObject() %>%
    AddMetaData(metadata = str_extract(x, pattern = "GSM.*$"), col.name = "sample_id")
  
})


# sort metadata -----------------------------------------------------------
meta <- read.csv("GSE157477_iel_crohns/metadata.csv")

# making sure that the samples are in the same order as the metadata
sapply(df, function(x) str_extract(x$sample_id, "^.*?(?=_)")[1]) == meta$accession

# annotate data -----------------------------------------------------------
add_data <- c("donor", "sample", "disease", "compartment")
for (c in add_data) {
  for (i in 1:length(df)) {
    
    df[[i]] <- AddMetaData(df[[i]], metadata = meta[, c][i], col.name = c)
    
  }
}

# add mito percentage -----------------------------------------------------
df <- lapply(df, function(x) {
  NormalizeData(x) %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent_mt")
})
df <- lapply(df, function(x) subset(x, percent_mt < 10))

# merging the data --------------------------------------------------------
df <- merge(x = df[[1]], y = df[2:11], merge.data = T)

# downstream processing ---------------------------------------------------
df <- df %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "sample_id")

df <- df %>%
  RunUMAP(dims = 1:30, reduction = "harmony") %>%
  FindNeighbors(dims = 1:30, reduction = "harmony") %>%
  FindClusters(res = 0.5)

df$group <- paste(df$disease, df$sample, sep = "_")

saveRDS(df, file = "GSE157477_iel_crohns/harmony_integrated.rds")

# annotation of data ------------------------------------------------------
df <- readRDS("colitis_analysis/GSE157477_iel_crohns/harmony_integrated.rds")
df <- subset(df, group %in% c("healthy_normal", "cd_inflamed"))

df$neat_group <- factor(ifelse(df$group == "healthy_normal", "Healthy", "Crohn's"), levels = c("Healthy", "Crohn's"))

df@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  mutate(disease = df$neat_group) %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(shape = 21, size = 0.7, alpha = 0.7, aes(fill = disease), stroke = 0.05) +
  scale_fill_manual(name = "Status",
                    values = c("gray90", "tomato")) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

ggsave("colitis_analysis/GSE157477_iel_crohns/h_c_umap.png", dpi = 600,
       height = 5, width = 5)

df$broad <- case_when(df$seurat_clusters %in% c(0, 2, 6, 9, 10) ~ "cd8",
                      df$seurat_clusters %in% c(1, 3, 5, 7) ~ "cd4",
                      df$seurat_clusters %in% 4 ~ "nk",
                      df$seurat_clusters %in% 8 ~ "b",
                      df$seurat_clusters %in% 11 ~ "mac/mono")

df$fine <- case_when(df$seurat_clusters == 0 ~ "CD8 TEMRA",
                     df$seurat_clusters == 1 ~ "CD4 Tcm",
                     df$seurat_clusters == 2 ~ "CD8",
                     df$seurat_clusters == 3 ~ "CD4 Tem",
                     df$seurat_clusters == 4 ~ "NK",
                     df$seurat_clusters == 5 ~ "CD4 Treg",
                     df$seurat_clusters == 6 ~ "CD8 Trm",
                     df$seurat_clusters == 7 ~ "CD4 Trm",
                     df$seurat_clusters == 8 ~ "B cells",
                     df$seurat_clusters == 9 ~ "CD8 CTL",
                     df$seurat_clusters == 10 ~ "CD8 proliferating",
                     df$seurat_clusters == 11 ~ "Mono/mac")


# umap --------------------------------------------------------------------
col_scheme <- RColorBrewer::brewer.pal(n = 12, "Paired")

cluster_df <- df@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  mutate(cluster = df$fine)

ggplot(cluster_df, aes(UMAP_1, UMAP_2)) +
  geom_point(shape = 21, size = 0.7, alpha = 0.7, aes(fill = cluster), stroke = 0.05) +
  scale_fill_manual(name = "Cluster",
                    values = col_scheme) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(5, "mm")) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3)))

ggsave("colitis_analysis/GSE157477_iel_crohns/cluster_umap.png", dpi = 600,
       height = 5, width = 5)


# pathway analysis --------------------------------------------------------
split_df <- SplitObject(df, split.by = "group")
clusters <- split_df$healthy_normal$fine %>% unique()

scpa_out <- list()
for (i in clusters) {
  
  print(i)
  
  p1 <- seurat_extract(split_df$healthy_normal, meta1 = "fine", value_meta1 = i)
  p2 <- seurat_extract(split_df$cd_inflamed, meta1 = "fine", value_meta1 = i)
  
  scpa_out[[as.character(i)]] <- compare_pathways(samples = list(p1, p2), 
                                                  pathways = gene_sets,
                                                  downsample = 500,
                                                  max_genes = 500,
                                                  parallel = T,
                                                  cores = 12) %>%
    mutate(cluster = i)
  
}

scpa_out <- bind_rows(scpa_out)
saveRDS(scpa_out, "colitis_analysis/GSE157477_iel_crohns/scpa_out.rds")

# arthritis data ----------------------------------------------------------

# read in data ------------------------------------------------------------
files <- list.files("E-MTAB-9492_arthritis/raw_data", pattern = "", full.names = T)
df <- lapply(files, function(x) Read10X(x) %>% 
               CreateSeuratObject() %>%
               NormalizeData() %>%
               PercentageFeatureSet(pattern = "^MT-", col.name = "percent_mt"))

# read metadata -----------------------------------------------------------
meta <- read.csv("E-MTAB-9492_arthritis/metdata.csv")
meta <- meta[!duplicated(meta$sample_id), ] %>%
  mutate(sample_id = str_extract(sample_id, "[^_]+"))
meta <- meta[c(7, 8, 1, 4, 2, 5, 3, 6), ] %>%
  mutate(filname = files)


# annotate data -----------------------------------------------------------
for (i in 1:length(df)) {
  
  df[[i]]$tissue <- meta$tissue[i]
  df[[i]]$sample_id <- meta$sample_id[i]
  df[[i]]$seq_method <- meta$seq_method[i]
  df[[i]]$cell_type <- meta$cell[i]
  df[[i]]$filename <- meta$filname[i]
  df[[i]]$sex <- meta$sex[i]
  
}

# merging the data --------------------------------------------------------
df <- merge(x = df[[1]], 
            y = df[2:length(df)],
            merge.data = T)

# filtering base on mito and counts ---------------------------------------
plot(log2(df$nCount_RNA), log2(df$nFeature_RNA), 
     pch = 21, bg = "#f0635920", lwd = 0.1,
     main = "Keep cells in top right quadrant")
abline(v = 10.55, h = 9.1)
df <- subset(df, percent_mt < 10 & nCount_RNA > 2^10.55 & nFeature_RNA > 2^9.1)

# downstream processing ---------------------------------------------------
df <- df %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

df <- RunHarmony(df, group.by.vars = "sample_id", plot_convergence = T)

df <- df %>%
  RunUMAP(dims = 1:20, reduction = "harmony") %>%
  FindNeighbors(dims = 1:20, reduction = "harmony") %>%
  FindClusters(resolution = 0.3)


# annotate data -----------------------------------------------------------
df$fine <- case_when(df$seurat_clusters == 0 ~ "Naive/Tcm",
                     df$seurat_clusters == 1 ~ "CD8 Tem",
                     df$seurat_clusters == 2 ~ "CD4 Tem",
                     df$seurat_clusters == 3 ~ "CD4 Galectin+",
                     df$seurat_clusters == 4 ~ "CD8 Tem 2",
                     df$seurat_clusters == 5 ~ "CD4 Treg",
                     df$seurat_clusters == 6 ~ "CD8 Trm",
                     df$seurat_clusters == 7 ~ "MAIT",
                     df$seurat_clusters == 8 ~ "TRBV7+",
                     df$seurat_clusters == 9 ~ "CD8 CTL",
                     df$seurat_clusters == 10 ~ "CD8 Trm 2",
                     df$seurat_clusters == 11 ~ "CD4 AQP+",
                     df$seurat_clusters == 12 ~ "CD3-")


# plot umap ---------------------------------------------------------------
df <- subset(df, tissue != "synovial membrane")

df@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  mutate(disease = df$tissue) %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(shape = 21, size = 0.7, alpha = 0.7, aes(fill = disease), stroke = 0.05) +
  scale_fill_manual(name = "Status",
                    values = c("gray90", "tomato")) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

ggsave("colitis_analysis/E-MTAB-9492_arthritis/h_c_umap.png", dpi = 600,
       height = 5, width = 5)

df@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  mutate(cluster = df$fine) %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(shape = 21, size = 0.7, alpha = 0.7, aes(fill = cluster), stroke = 0.05) +
  scale_fill_manual(values = c(col_scheme, "gray85")) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

ggsave("colitis_analysis/E-MTAB-9492_arthritis/umap_clusters.png", dpi = 600,
       height = 5, width = 5)

# pathway analysis --------------------------------------------------------
clusters <- unique(df$fine)
scpa_out <- list()
for (i in clusters) {
  
  print(i)
  
  p1 <- seurat_extract(df, 
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "tissue", value_meta2 = "blood")
  p2 <- seurat_extract(df, 
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "tissue", value_meta2 = "synovial fluid")
  
  scpa_out[[as.character(i)]] <- compare_pathways(samples = list(p1, p2), 
                                                  pathways = gene_sets,
                                                  downsample = 500,
                                                  max_genes = 500,
                                                  parallel = T,
                                                  cores = 12) %>%
    mutate(cluster = i)
  
}

scpa_out <- bind_rows(scpa_out)
saveRDS(scpa_out, "colitis_analysis/E-MTAB-9492_arthritis/scpa_out.rds")


# jia data ----------------------------------------------------------------

# read in data ------------------------------------------------------------
files <- list.files("GSE160097_scrna_jia_memory/data", full.names = T)
df <- lapply(files, function(x) Read10X_h5(x) %>%
               CreateSeuratObject() %>%
               PercentageFeatureSet(pattern = "^MT-", col.name = "percent_mt"))

# filter data -------------------------------------------------------------
df <- lapply(df, function(x) subset(x, percent_mt < 10 & percent_mt > 1.5))


# add metadata ------------------------------------------------------------
meta <- data.frame(tissue = str_extract(string = files, pattern = "PB|SF|blood"),
                   cell = str_extract(string = files, pattern = "Treg|CD4|CD8"),
                   donor = str_extract(string = files, pattern = "p[0-9]"))

for (i in 1:length(df)) {
  
  df[[i]]$tissue <- meta$tissue[i]
  df[[i]]$cell <- meta$cell[i]
  df[[i]]$donor <- meta$donor[i]
  
}

# merging the data --------------------------------------------------------
df <- merge(x = df[[1]], 
            y = df[2:length(df)])


# downstream processing and more filtering --------------------------------
df$donor_tissue <- paste(df$donor, df$tissue, sep = "_")

df <- df %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

df <- RunHarmony(df, group.by.vars = "donor_tissue", plot_convergence = T)

ElbowPlot(df, ndims = 50, reduction = "harmony")

df <- df %>%
  RunUMAP(dims = 1:30, reduction = "harmony") %>%
  FindNeighbors(dims = 1:30, reduction = "harmony") %>%
  FindClusters(resolution = 0.3)

df <- FindClusters(df, resolution = 0.5)

df <- subset(df, seurat_clusters != 14)

df <- FindVariableFeatures(df) %>%
  ScaleData() %>%
  RunPCA()
df <- RunHarmony(df, group.by.vars = "donor_tissue", plot_convergence = T)
ElbowPlot(df, ndims = 50, reduction = "harmony")

df <- df %>%
  RunUMAP(dims = 1:30, reduction = "harmony") %>%
  FindNeighbors(dims = 1:30, reduction = "harmony") %>%
  FindClusters(resolution = 0.5)

# data annotation ---------------------------------------------------------
df$fine <- case_when(df$seurat_clusters == 0 ~ "cd4_treg_1",
                     df$seurat_clusters == 1 ~ "cd8_tem_1",
                     df$seurat_clusters == 2 ~ "cd4_tcm",
                     df$seurat_clusters == 3 ~ "cd4_tem",
                     df$seurat_clusters == 4 ~ "cd8_tem_2",
                     df$seurat_clusters == 5 ~ "cd8_tcm",
                     df$seurat_clusters == 6 ~ "cd4_cd40lg",
                     df$seurat_clusters == 7 ~ "cd4_treg_2",
                     df$seurat_clusters == 8 ~ "cd4_aqp3",
                     df$seurat_clusters == 9 ~ "cd8_tcm_2",
                     df$seurat_clusters == 10 ~ "cd8_tc1",
                     df$seurat_clusters == 11 ~ "cd4_ifi",
                     df$seurat_clusters == 12 ~ "cd8_tem_3",
                     df$seurat_clusters == 13 ~ "?")

DimPlot(df, split.by = "tissue", group.by = "cell") + DimPlot(df, split.by = "tissue", label = T)
VlnPlot(df, features = c("IFNG", "IL10"), split.by = "tissue", pt.size = 0)


# umap --------------------------------------------------------------------
df@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  mutate(disease = df$tissue) %>%
  slice_sample(n = nrow(df)) %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(shape = 21, size = 0.7, alpha = 0.7, aes(fill = disease), stroke = 0.05) +
  scale_fill_manual(name = "Status",
                    values = c("gray90", "tomato")) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

ggsave("colitis_analysis/GSE160097_scrna_jia_memory/h_c_umap.png", dpi = 600,
       height = 5, width = 5)

df@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  mutate(cluster = df$fine) %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(shape = 21, size = 0.7, alpha = 0.7, aes(fill = cluster), stroke = 0.05) +
  scale_fill_manual(values = c(col_scheme, "gray85", "gray60")) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        plot.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")

ggsave("colitis_analysis/GSE160097_scrna_jia_memory/umap_clusters.png", dpi = 600,
       height = 5, width = 5)


# pathway analysis --------------------------------------------------------
clusters <- sort(unique(df$fine))
scpa_out <- list()
for (i in clusters) {
  
  print(i)
  
  p1 <- seurat_extract(df, 
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "tissue", value_meta2 = "blood")
  p2 <- seurat_extract(df, 
                       meta1 = "fine", value_meta1 = i,
                       meta2 = "tissue", value_meta2 = "SF")
  
  scpa_out[[as.character(i)]] <- compare_pathways(samples = list(p1, p2), 
                                                  pathways = gene_sets,
                                                  downsample = 500,
                                                  max_genes = 500,
                                                  parallel = T,
                                                  cores = 12) %>%
    mutate(cluster = i)
  
}

scpa_out <- bind_rows(scpa_out)
saveRDS(scpa_out, file = "colitis_analysis/GSE160097_scrna_jia_memory/scpa_out.rds")

# read in pathway analysis results ----------------------------------------
scpa_d1 <- readRDS("colitis_analysis/GSE157477_iel_crohns/scpa_out.rds")
scpa_d2 <- readRDS("colitis_analysis/E-MTAB-9492_arthritis/scpa_out.rds")
scpa_d3 <- readRDS("colitis_analysis/GSE160097_scrna_jia_memory/scpa_out.rds")

scpa_d1 <- filter(scpa_d1, grepl("CD4", x = cluster))
scpa_d2 <- filter(scpa_d2, grepl("CD4|Naive/Tcm", x = cluster))
scpa_d3 <- filter(scpa_d3, grepl("cd4", x = cluster))

all_scpa <- list(scpa_d1, scpa_d2, scpa_d3)
for (i in 1:3) {
  
  all_scpa[[i]] <- mutate(all_scpa[[i]], dataset = paste0("dataset", i))
  
}

# summary heatmap ---------------------------------------------------------

highlight_pathways <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "prost_sig", "C5")

hm_data <- all_scpa %>%
  bind_rows() %>%
  mutate(group = paste(dataset, cluster, sep = "_")) %>%
  select(Pathway, qval, group) %>%
  pivot_wider(names_from = "group", values_from = qval) %>%
  column_to_rownames("Pathway")

hm_col <- colorRamp2(breaks = c(0, 2, 8), colors = c("#0086FF", "white", "#FF7900"))

position <- which(rownames(hm_data) %in% highlight_pathways)

col_an <- HeatmapAnnotation(Pathway = anno_mark(at = position,
                                                labels = c("IFNG response", "PTGIR/S", "C5AR2"),
                                                labels_gp = gpar(fontsize = 9),
                                                link_height = unit(3.5, "mm"),
                                                padding = unit(3.5, "mm"),
                                                link_gp = gpar(lwd = 0.5), 
                                                labels_rot = 35))

dataset <- substr(colnames(hm_data), 1, 8)
dataset <- case_when(dataset == "dataset1" ~ "Crohn's",
                     dataset == "dataset2" ~ "PSA",
                     dataset == "dataset3" ~ "JIA")

row_an <- rowAnnotation(Disease = dataset,
                        col = list(Disease = c("Crohn's" = "gray90", "PSA" = "red", "JIA" = "cornflowerblue")),
                        gp = gpar(col = "white", lwd = 0.3),
                        simple_anno_size = unit(0.35, "cm"), 
                        annotation_name_side = "top")

png("~/My Drive/analysis_for_people/jubayer/paper_figures/heatmap.png", res = 600, units = "in", width = 5.8, height = 3.4)
Heatmap(t(hm_data),
        border = T,
        name = "Qval",
        top_annotation = col_an,
        show_column_names = F,
        left_annotation = row_an,
        col = hm_col,
        row_dend_width = unit(3.5, "mm"),
        show_column_dend = F,
        row_labels = c("Tcm", "Trm", "Tem", "Treg", "Naive/Tcm", "AQP+",
                       "Galectin+", "Tem", "Treg", "Tem", "CD40LG+", "IFI+",
                       "Tcm", "Tem", "Treg", "Treg"))
dev.off()

# plot pathway ranks ------------------------------------------------------
rank_df <- scpa_d3 %>%
  filter(cluster != "cd4_treg_2") %>%
  mutate(cluster = case_when(cluster == "cd4_aqp3" ~ "Tem",
                             cluster == "cd4_cd40lg" ~ "CD40LG+",
                             cluster == "cd4_ifi" ~ "IFI+",
                             cluster == "cd4_tcm" ~ "Tcm",
                             cluster == "cd4_tem" ~ "Tem 2",
                             cluster == "cd4_treg_1" ~ "Treg"),
         cluster = factor(cluster, levels = c("Tcm", "Tem", "Tem 2", "CD40LG+", "IFI+", "Treg")))

col_scheme <- RColorBrewer::brewer.pal(3, "Set2")

ggplot(rank_df, aes(reorder(Pathway, qval), qval)) +
  geom_point(shape = 21, fill = "gray80", size = 2.5, alpha = 0.7, stroke = 0) +
  geom_point(data = subset(rank_df, Pathway %in% highlight_pathways), 
             shape = 21, aes(fill = Pathway), size = 3.2, stroke = 0.3) +
  facet_wrap(~cluster, ncol = 3) +
  scale_fill_manual(values = c("HALLMARK_INTERFERON_GAMMA_RESPONSE" = col_scheme[1], "prost_sig" = col_scheme[2], "C5" = col_scheme[3]),
                    label = c("HALLMARK_INTERFERON_GAMMA_RESPONSE" = "IFNG response", 
                              "prost_sig" = "PTGIR/S", 
                              "C5" = "C5AR2")) +
  labs(y = "Qval", title = "Pathway changes in JIA") +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(4.5, "mm"),
        legend.text = element_text(size = 10),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(hjust = 0, size = 10))

ggsave("paper_figures/rank_plot.png", width = 5.5, height = 3, dpi = 600)















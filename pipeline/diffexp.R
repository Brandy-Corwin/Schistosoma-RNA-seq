library(tidyverse)
library(broom)
library(DESeq2, lib.loc="/data/users/corwinbm5021/.conda/envs/biol343_project/lib/R/library/")


# For star
star_counts_df <- read_tsv("counting/counts/star_counts.tsv", comment = "#") |>
             mutate(across(where(is.numeric), as.integer))
write.csv(star_counts_df, "count_matrices/star_counts_df.csv")

star_counts_summary <- star_counts_df |>
    select(Geneid, contains('star.bam')) |>
    rename_with(~str_remove(., "counting/dedup/star.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)
write.csv(star_counts_summary, "count_matrices/star_counts_summary.csv")

star_sample_summary <- star_counts_df |>
    select(Geneid, contains('star.bam')) |>
    rename_with(~str_remove(., "counting/dedup/star.bam:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)
write.csv(star_sample_summary, "count_matrices/star_sample_summary.csv")

star_genes_to_remove = star_sample_summary$Geneid

star_counts_filt <- star_counts_summary |>
    filter(!Geneid %in% star_genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)

star_counts_m <- star_counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(star_counts_m) <- star_counts_filt$Geneid
write.csv(star_counts_m, "count_matrices/star_counts_m.csv")

star_dists <- dist(t(star_counts_m))

star_dists_df <- as.matrix(star_dists) |>
    as_tibble(rownames = 'sample')

star_dist_plot <- star_dists_df |>
    pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
    ggplot(aes(x = sample, y = comp, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    NULL
ggsave("plots/star_dist_plot_before_diffexp.png")

star_pca_fit <- t(log10(star_counts_m + 1)) |> 
  prcomp(scale = TRUE)

star_pca_plot_before_deseq2 <- star_pca_fit |>
    augment(t(star_counts_m)) |>
    dplyr::rename(sample = .rownames) |>
    mutate(sample = str_remove(sample, '-[0-9]+')) |>
    ggplot(aes(.fittedPC1, .fittedPC2, color = sample)) + 
    geom_point(size = 5)
ggsave("plots/star_pca_before_deseq2.png")

star_metadata <- data.frame(sample_id = colnames(star_counts_m)) |>
    mutate(tissue = str_sub(sample_id, 1, 3),
           rep = str_sub(sample_id, 5))
rownames(star_metadata) <- star_metadata$sample_id
star_metadata <- select(star_metadata, -sample_id)

all(rownames(star_metadata) == colnames(star_counts_m))

star_dds <- DESeqDataSetFromMatrix(countData = star_counts_m,
                              colData = star_metadata,
                              design = ~ tissue)
star_dds <- DESeq(star_dds)

star_res <- results(star_dds)

star_volcano_data <- as_tibble(star_res, rownames = 'gene_id')

star_volcano_plot <- star_volcano_data |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2)
ggsave("plots/star_volcano_plot.png")

star_vsd <- varianceStabilizingTransformation(star_dds)
write.csv(assay(star_vsd), "count_matrices/star_vsd.csv")

star_pca_after_deseq2 <- plotPCA(star_vsd, intgroup = 'tissue')
ggsave("plots/star_pca_after_deseq2.png")



# For hisat
hisat_counts_df <- read_tsv("counting/counts/hisat_counts.tsv", comment = "#") |>
             mutate(across(where(is.numeric), as.integer))
write.csv(hisat_counts_df, "count_matrices/hisat_counts_df.csv")

hisat_counts_summary <- hisat_counts_df |>
    select(Geneid, contains('hisat.bam')) |>
    rename_with(~str_remove(., "counting/dedup/hisat.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)
write.csv(hisat_counts_summary, "count_matrices/hisat_counts_summary.csv")

hisat_sample_summary <- hisat_counts_df |>
    select(Geneid, contains('hisat.bam')) |>
    rename_with(~str_remove(., "counting/dedup/hisat.bam:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)
write.csv(hisat_sample_summary, "count_matrices/hisat_sample_summary.csv")

hisat_genes_to_remove = hisat_sample_summary$Geneid

hisat_counts_filt <- hisat_counts_summary |>
    filter(!Geneid %in% hisat_genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)

hisat_counts_m <- hisat_counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(hisat_counts_m) <- hisat_counts_filt$Geneid
write.csv(hisat_counts_m, "count_matrices/hisat_counts_m.csv")

hisat_dists <- dist(t(hisat_counts_m))

hisat_dists_df <- as.matrix(hisat_dists) |>
    as_tibble(rownames = 'sample')

hisat_dist_plot <- hisat_dists_df |>
    pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
    ggplot(aes(x = sample, y = comp, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    NULL
ggsave("plots/hisat_dist_plot_before_diffexp.png")

hisat_pca_fit <- t(log10(hisat_counts_m + 1)) |> 
  prcomp(scale = TRUE)

hisat_pca_plot_before_deseq2 <- hisat_pca_fit |>
    augment(t(hisat_counts_m)) |>
    dplyr::rename(sample = .rownames) |>
    mutate(sample = str_remove(sample, '-[0-9]+')) |>
    ggplot(aes(.fittedPC1, .fittedPC2, color = sample)) + 
    geom_point(size = 5)
ggsave("plots/hisat_pca_before_deseq2.png")

hisat_metadata <- data.frame(sample_id = colnames(hisat_counts_m)) |>
    mutate(tissue = str_sub(sample_id, 1, 3),
           rep = str_sub(sample_id, 5))
rownames(hisat_metadata) <- hisat_metadata$sample_id
hisat_metadata <- select(hisat_metadata, -sample_id)

all(rownames(hisat_metadata) == colnames(hisat_counts_m))

hisat_dds <- DESeqDataSetFromMatrix(countData = hisat_counts_m,
                              colData = hisat_metadata,
                              design = ~ tissue)
hisat_dds <- DESeq(hisat_dds)

hisat_res <- results(hisat_dds)

hisat_volcano_data <- as_tibble(hisat_res, rownames = 'gene_id')

hisat_volcano_plot <- hisat_volcano_data |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2)
ggsave("plots/hisat_volcano_plot.png")

hisat_vsd <- varianceStabilizingTransformation(hisat_dds)
write.csv(assay(hisat_vsd), "count_matrices/hisat_vsd.csv")

hisat_pca_after_deseq2 <- plotPCA(hisat_vsd, intgroup = 'tissue')
ggsave("plots/hisat_pca_after_deseq2.png")


library(tidyverse)
library(broom)
library(DESeq2, lib.loc="/data/users/corwinbm5021/.conda/envs/biol343_project/lib/R/library/")


# For star
counts_df <- read_tsv('/data/classes/2025/spring/biol443/course_files/rnaseq_data/counts.tsv', comment = '#') |>
             mutate(across(where(is.numeric), as.integer))

counts_summary <- counts_df |>
    select(Geneid, contains('.bam')) |>
    rename_with(~str_remove(., "dedup/star/.*:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)

sample_summary <- counts_df |>
    select(Geneid, contains('.bam')) |>
    rename_with(~str_remove(., "dedup/star/.*:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)

genes_to_remove = sample_summary$Geneid

counts_filt <- counts_summary |>
    filter(!Geneid %in% genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)

counts_m <- counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(counts_m) <- counts_filt$Geneid

metadata <- data.frame(sample_id = colnames(counts_m)) |>
    mutate(tissue = str_sub(sample_id, 1, 3),
           rep = str_sub(sample_id, 5))
rownames(metadata) <- metadata$sample_id
metadata <- select(metadata, -sample_id)

all(rownames(metadata) == colnames(counts_m))

dds <- DESeqDataSetFromMatrix(countData = counts_m,
                              colData = metadata,
                              design = ~ tissue)
dds <- DESeq(dds)

res <- results(dds)

volcano_data <- as_tibble(res, rownames = 'gene_id')

volcano_plot <- volcano_data |> 
    ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point() +
    geom_vline(xintercept = 2) +
    geom_vline(xintercept = -2)
ggsave("plots_reduced/volcano_plot.png")

vsd <- varianceStabilizingTransformation(dds)

pca_after_deseq2 <- plotPCA(vsd, intgroup = 'tissue')
ggsave("plots_reduced/pca_plot.png")


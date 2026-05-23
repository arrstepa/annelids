
# Libraries and data preparation

library(dplyr)
library(stringr)
library(pheatmap)
library(readr)
library(tidyr)
library(ggplot2)


segment_labels <- as.character(1:12)

species_meta <- list(
  Pele = list(
    annot = data.frame(
      Region = c("Pharynx","Pharynx","Midgut","Midgut","Midgut","Midgut",
                 "Midgut","Midgut","Midgut","Midgut","Midgut","Hindgut"),
      row.names = segment_labels
    ),
    colors = list(Region = c(
      Pharynx = "#00cb00", Midgut = "#bf8c00", Hindgut = "#d7144c"
    ))
  ),
  Pdum = list(
    annot = data.frame(
      Region = c("Pharynx","Pharynx","Esophagus","Midgut","Midgut","Midgut",
                 "Midgut","Midgut","Midgut","Midgut","Midgut","Hindgut"),
      row.names = segment_labels
    ),
    colors = list(Region = c(
      Pharynx   = "#00cb00", Esophagus = "#91ad01",
      Midgut    = "#bf8c00", Hindgut   = "#d7144c"
    ))
  ),
  Amar = list(
    annot = data.frame(
      Region = c("Pharynx","Pharynx","Esophagus","Esophagus",
                 "Stomach","Stomach","Midgut","Midgut",
                 "Midgut","Midgut","Midgut","Hindgut"),
      row.names = segment_labels
    ),
    colors = list(Region = c(
      Pharynx   = "#00cb00", Esophagus = "#91ad01",
      Stomach   = "#777916", Midgut    = "#bf8c00", Hindgut = "#d7144c"
    ))
  )
)

SPECIES_ORDER <- c("Pele", "Pdum", "Amar")

# Data load

pfam <- read.table("demodata/pfam_gene_clean.tsv", sep="\t", stringsAsFactors=FALSE,
                   col.names=c("geneID","pfam","domain"))

tpm_raw <- read.table("demodata/all_tpm.txt", sep="\t", header=TRUE)
colnames(tpm_raw)[1] <- "geneID"
colnames(tpm_raw)[2:13] <- segment_labels
tpm <- tpm_raw[, c("geneID", segment_labels)]


# Merge pfam domains and TPM tables


merged <- pfam %>%
  inner_join(tpm, by="geneID")

# Make different tables for species

pele_full <- merged %>% filter(str_starts(geneID, "Pele"))
pdum_full <- merged %>% filter(str_starts(geneID, "Pdum"))
amar_full <- merged %>% filter(str_starts(geneID, "Amar"))


# Find common domains for all species in order to compare correctly


get_domains <- function(df) unique(df$domain[!is.na(df$domain) & df$domain != ""])

domains_pele <- get_domains(pele_full)
domains_pdum <- get_domains(pdum_full)
domains_amar <- get_domains(amar_full)

common_domains <- Reduce(intersect, list(domains_pele, domains_pdum, domains_amar))


pele_common <- pele_full %>% filter(domain %in% common_domains)
pdum_common <- pdum_full %>% filter(domain %in% common_domains)
amar_common <- amar_full %>% filter(domain %in% common_domains)

species_full_list   <- list(Pele = pele_full,   Pdum = pdum_full,   Amar = amar_full)
species_common_list <- list(Pele = pele_common, Pdum = pdum_common, Amar = amar_common)

# Utilitary functions for heatmaps

heatmap_colors <- colorRampPalette(c("#2166AC","#F7F7F7","#D6604D"))(100)
heatmap_breaks <- seq(-3, 3, length.out=101)

make_heatmap <- function(mat, title, col_annot, annot_colors,
                         row_annot=NULL, show_rownames=FALSE,
                         fontsize_row=7, gaps_row=NULL, row_clust=TRUE) {
  if (is.null(mat) || nrow(mat) < 2) {
    message("Слишком мало строк для: ", title)
    return(invisible(NULL))
  }
  pheatmap(mat,
           color             = heatmap_colors,
           breaks            = heatmap_breaks,
           cluster_rows      = row_clust,
           cluster_cols      = FALSE,
           show_rownames     = show_rownames,
           show_colnames     = TRUE,
           fontsize_row      = fontsize_row,
           annotation_col    = col_annot,
           annotation_row    = row_annot,
           annotation_colors = annot_colors,
           main              = title,
           border_color      = NA,
           treeheight_row    = 30,
           gaps_row          = gaps_row,
           na_col            = "white")
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.nan(z) | is.na(z)] <- 0
  z
}

# Gene clustering via anterio-posterior expression profiles

run_gene_clustering <- function(species_df, species_name,
                                n_clusters = 8, min_tpm = 0.5) {
  
  df <- species_df %>%
    distinct(geneID, .keep_all=TRUE) %>%
    dplyr::select(geneID, all_of(segment_labels))
  
  max_tpm <- apply(df[, segment_labels], 1, max, na.rm=TRUE)
  df <- df[max_tpm >= min_tpm, ]
  
  mat <- as.matrix(df[, segment_labels])
  rownames(mat) <- df$geneID
  mat_z <- row_zscore(mat)
  
  set.seed(42)
  km <- kmeans(mat_z, centers=n_clusters, nstart=25, iter.max=200)
  
  profiles <- lapply(1:n_clusters, function(k) {
    idx <- km$cluster == k
    colMeans(mat_z[idx, , drop=FALSE], na.rm=TRUE)
  })
  peak_pos    <- sapply(profiles, which.max)
  cluster_map <- setNames(rank(peak_pos, ties.method="first"), 1:n_clusters)
  cluster_ordered <- cluster_map[km$cluster]
  names(cluster_ordered) <- rownames(mat_z)
  
  gene_peak  <- apply(mat_z, 1, which.max)
  gene_order <- order(cluster_ordered, gene_peak)
  mat_sorted    <- mat_z[gene_order, ]
  clusters_sorted <- cluster_ordered[gene_order]
  
  print(table(clusters_sorted))
  
  list(
    mat      = mat_sorted,
    clusters = clusters_sorted,
    profiles = profiles[order(peak_pos)],
    peak_pos = sort(peak_pos)
  )
}

# Main heatmap with gene clusters

plot_gene_cluster_heatmap <- function(res, species_name) {
  meta <- species_meta[[species_name]]
  cluster_levels <- sort(unique(res$clusters))
  
  row_annot <- data.frame(
    Cluster = factor(paste0("C", res$clusters), levels=paste0("C", cluster_levels)),
    row.names = rownames(res$mat)
  )
  
  cluster_colors <- setNames(
    colorRampPalette(c(
      "#d73027","#f46d43","#fdae61","#fee090",
      "#abd9e9","#74add1","#4575b4","#762a83"
    ))(length(cluster_levels)),
    paste0("C", cluster_levels)
  )
  
  annot_colors <- c(meta$colors, list(Cluster = cluster_colors))
  gaps <- cumsum(as.numeric(table(res$clusters)))[-length(cluster_levels)]
  
  make_heatmap(
    mat           = res$mat,
    title         = paste0(species_name, " - all transcripts (n=", nrow(res$mat), ")"),
    col_annot     = meta$annot,
    annot_colors  = annot_colors,
    row_annot     = row_annot,
    show_rownames = FALSE,
    row_clust     = FALSE,
    gaps_row      = gaps
  )
}


# Clusterisation and heatmaps

res_pele <- run_gene_clustering(pele_common, "Pele", n_clusters=11)
res_pdum <- run_gene_clustering(pdum_common, "Pdum", n_clusters=12)
res_amar <- run_gene_clustering(amar_common, "Amar", n_clusters=9)

plot_gene_cluster_heatmap(res_pele, "P. elegans")
plot_gene_cluster_heatmap(res_pdum, "P. dumerilii")
plot_gene_cluster_heatmap(res_amar, "A. marina")

# List of domains in each cluster

extract_cluster_annotation <- function(res, species_df, species_name) {
  cluster_df <- data.frame(
    geneID  = names(res$clusters),
    cluster = as.integer(res$clusters)
  )
  
  annotated <- cluster_df %>%
    left_join(
      species_df %>% distinct(geneID, domain, pfam),
      by="geneID"
    ) %>%
    arrange(cluster, geneID)
  
  domain_summary <- annotated %>%
    filter(!is.na(domain), domain != "") %>%
    group_by(cluster, domain) %>%
    summarise(n_genes=n(), .groups="drop") %>%
    arrange(cluster, desc(n_genes))
  
  write.table(annotated,
              file = paste0(
                "results/tables/cluster_genes_", species_name,".tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  write.table(domain_summary,
              file = paste0(
                "results/tables/cluster_domains_", species_name, ".tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  
  list(genes=annotated, domains=domain_summary)
}


annot_pele <- extract_cluster_annotation(res_pele, pele_full, "P. elegans")
annot_pdum <- extract_cluster_annotation(res_pdum, pdum_full, "P. dumerilii")
annot_amar <- extract_cluster_annotation(res_amar, amar_full, "A. marina")

# Enrichment barplots for each cluster

plot_cluster_enrichment <- function(annot, species_name, top_n=15) {
  df <- annot$domains %>%
    filter(!str_detect(domain, "DUF|UU[0-9]"),
           !is.na(domain), domain != "") %>%
    group_by(cluster) %>%
    slice_max(n_genes, n=top_n) %>%
    ungroup()
  
  ggplot(df, aes(x=reorder(domain, n_genes), y=n_genes, fill=factor(cluster))) +
    geom_col(show.legend=FALSE) +
    coord_flip() +
    facet_wrap(~paste0("Cluster ", cluster), scales="free_y", ncol=2) +
    scale_fill_manual(
      values = colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1")
      )(length(unique(df$cluster)))
    ) +
    labs(title=paste("Domain enrichment by cluster:", species_name),
         x="Pfam-domain", y="Number of genes") +
    theme_bw(base_size=10) +
    theme(
      strip.text = element_text(face="bold", size=10),
      axis.text.y = element_text(size=5),
      axis.text.x = element_text(size=7),
      axis.title  = element_text(size=9)
    )
}

print(plot_cluster_enrichment(annot_pele, "P.elegans"))
print(plot_cluster_enrichment(annot_pdum, "P.dumerilii"))
print(plot_cluster_enrichment(annot_amar, "A.marina"))

# CLuster profiles

plot_cluster_profiles <- function(res, species_name, n_clusters=8) {
  cluster_colors <- colorRampPalette(
    c("#d73027","#f46d43","#fdae61","#74add1","#4575b4","#762a83","#1b7837","#8c510a")
  )(n_clusters)
  
  profiles_df <- lapply(1:n_clusters, function(k) {
    idx <- res$clusters == k
    if (sum(idx) == 0) return(NULL)
    data.frame(
      Segment = factor(segment_labels, levels=segment_labels),
      mean_z  = colMeans(res$mat[idx, , drop=FALSE], na.rm=TRUE),
      Cluster = paste0("C", k, " (n=", sum(idx), ")")
    )
  }) %>% bind_rows()
  
  p_overview <- ggplot(profiles_df,
                       aes(x=Segment, y=mean_z, group=Cluster, color=Cluster)) +
    geom_hline(yintercept=0, linetype="dashed", color="grey60") +
    geom_line(linewidth=1.2) +
    geom_point(size=2.5) +
    facet_wrap(~Cluster, ncol=2) +
    scale_color_manual(values=setNames(cluster_colors, unique(profiles_df$Cluster))) +
    theme_bw(base_size=11) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
    labs(title=paste("Cluster profile", species_name),
         x="Anterior → Posterior", y="Mean z-score")
  print(p_overview)
  
  for (k in 1:n_clusters) {
    idx <- res$clusters == k
    if (sum(idx) < 2) next
    
    mat_k <- res$mat[idx, ]
    mean_profile <- colMeans(mat_k, na.rm=TRUE)
    
    df_long <- as.data.frame(mat_k) %>%
      tibble::rownames_to_column("geneID") %>%
      pivot_longer(-geneID, names_to="Segment", values_to="zscore") %>%
      mutate(Segment=factor(Segment, levels=segment_labels))
    
    mean_df <- data.frame(
      Segment=factor(segment_labels, levels=segment_labels),
      zscore=mean_profile
    )
    
    p_k <- ggplot(df_long, aes(x=Segment, y=zscore, group=geneID)) +
      geom_line(alpha=min(0.6, 8/sqrt(nrow(mat_k))),
                color="grey40", linewidth=0.4) +
      geom_hline(yintercept=0, linetype="dashed", color="grey60") +
      geom_line(data=mean_df, aes(group=1), color="#d73027",
                linewidth=2.5, inherit.aes=FALSE,
                mapping=aes(x=Segment, y=zscore)) +
      theme_bw(base_size=11) +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      labs(title=paste0(species_name, " — Cluster ", k,
                        " (n=", sum(idx), " генов)"),
           x="A → P", y="z-score")
    print(p_k)
  }
}

plot_cluster_profiles(res_pele, "P. elegans",   n_clusters = 11)
plot_cluster_profiles(res_pdum, "P. dumerilii", n_clusters = 12)
plot_cluster_profiles(res_amar, "A. marina",    n_clusters = 9)

# Heatmaps for DUF domains (domains of unknown function)

make_combined_duf_heatmap <- function(
    top_n     = 90,
    min_genes = 3
) {
  get_duf_domain_mat <- function(species_df) {
    species_df %>%
      filter(domain %in% common_domains) %>%
      filter(str_detect(domain, "DUF")) %>%
      distinct(geneID, domain, .keep_all = TRUE) %>%
      group_by(domain) %>%
      summarise(
        across(all_of(segment_labels), mean, na.rm = TRUE),
        n_genes = n(),
        .groups = "drop"
      ) %>%
      filter(n_genes >= min_genes)
  }
  
  mat_p <- get_duf_domain_mat(pele_full)
  mat_d <- get_duf_domain_mat(pdum_full)
  mat_a <- get_duf_domain_mat(amar_full)
  
  shared <- Reduce(intersect, list(mat_p$domain, mat_d$domain, mat_a$domain))
  
  to_mat <- function(df) {
    m  <- df %>% filter(domain %in% shared) %>% arrange(domain)
    mx <- as.matrix(m[, segment_labels])
    rownames(mx) <- m$domain
    mx
  }
  
  m_p <- to_mat(mat_p)
  m_d <- to_mat(mat_d)
  m_a <- to_mat(mat_a)
  
  if (length(shared) > top_n) {
    mean_var <- (apply(m_p, 1, var, na.rm = TRUE) +
                   apply(m_d, 1, var, na.rm = TRUE) +
                   apply(m_a, 1, var, na.rm = TRUE)) / 3
    top_dom <- names(sort(mean_var, decreasing = TRUE))[1:top_n]
    m_p <- m_p[top_dom, ]
    m_d <- m_d[top_dom, ]
    m_a <- m_a[top_dom, ]
  }
  
  z_p <- row_zscore(m_p)
  z_d <- row_zscore(m_d)
  z_a <- row_zscore(m_a)
  
  colnames(z_p) <- paste0("Pele_", segment_labels)
  colnames(z_d) <- paste0("Pdum_", segment_labels)
  colnames(z_a) <- paste0("Amar_", segment_labels)
  
  spacer1  <- matrix(NA, nrow = nrow(z_p), ncol = 1,
                     dimnames = list(rownames(z_p), " "))
  spacer2  <- matrix(NA, nrow = nrow(z_p), ncol = 1,
                     dimnames = list(rownames(z_p), "  "))
  combined <- cbind(z_p, spacer1, z_d, spacer2, z_a)
  
  col_annot <- data.frame(
    Species = c(rep("P.elegans", 12), NA, rep("P.dumerilii", 12), NA, rep("A.marina", 12)),
    Region  = c(as.character(species_meta$Pele$annot$Region), NA,
                as.character(species_meta$Pdum$annot$Region), NA,
                as.character(species_meta$Amar$annot$Region)),
    row.names = colnames(combined)
  )
  
  all_regions <- c(species_meta$Pele$colors$Region,
                   species_meta$Pdum$colors$Region,
                   species_meta$Amar$colors$Region)
  all_regions <- all_regions[!duplicated(names(all_regions))]
  
  annot_colors <- list(
    Species = c("P.elegans"   = "#762a83",
                "P.dumerilii" = "#c51b7d",
                "A.marina"    = "#1b7837"),
    Region  = all_regions
  )
  
  row_clust <- hclust(dist(z_p))
  
  pheatmap(combined,
           color             = heatmap_colors,
           breaks            = heatmap_breaks,
           cluster_rows      = row_clust,
           cluster_cols      = FALSE,
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 7,
           fontsize_col      = 8,
           annotation_col    = col_annot,
           annotation_colors = annot_colors,
           na_col            = "white",
           main              = paste0("Common DUF domains across species ( top ", top_n, ")"),
           border_color      = NA,
           treeheight_row    = 40)
}


make_combined_duf_heatmap(top_n = 70)



# 13. РЕГЕНЕРАЦИЯ — ОТБОР КАНДИДАТОВ


regen_domains <- unique(c(
  "bHLH", "HLH", "Paired_box", "GATA", "Fork_head", "Ets", "T-box",
  "bZIP", "HMG_box", "AT_hook", "WD40", "Myb_DNA-binding",
  "Wnt_inhibitor", "Frizzled", "Dishevelled", "LRP", "BCL9",
  "Notch", "EGF_CA", "EGF", "Delta", "Serrate", "Ankyrin",
  "TGFb", "BMP", "Activin_recp", "GS",
  "FGF", "FGFR", "Hedgehog", "HH_signal",
  "Proliferating", "PCNA", "MCM", "Cyclin", "CDK", "RB",
  "Chromo", "Bromodomain", "PHD", "SET", "Histone", "HAT", "Polycomb",
  "Fibronectin", "Collagen", "Laminin", "Integrin", "Cadherin", "FN3",
  "Death", "Bcl", "IAP", "Caspase",
  "zf-C2H2", "zf-H2C2_2", "zf-H2C2_5", "zf-C4",
  "COesterase", "Bestrophin", "HRM", "7tm_1", "Kelch_1", "LRR_8"
))

regen_pattern <- paste(regen_domains, collapse = "|")

get_regen_genes <- function(species_df, species_name) {
  regen_genes <- species_df %>%
    filter(str_detect(domain, regen_pattern)) %>%
    distinct(geneID) %>%
    group_by(geneID) %>%
    summarise(sources = "Domain", .groups = "drop")
  regen_genes
}

regen_pele <- get_regen_genes(pele_full, "Pele")
regen_pdum <- get_regen_genes(pdum_full, "Pdum")
regen_amar <- get_regen_genes(amar_full, "Amar")

write.table(regen_pele, "results/tables/regen_candidates_Pele.tsv", sep="\t", row.names=FALSE)
write.table(regen_pdum, "results/tables/regen_candidates_Pdum.tsv", sep="\t", row.names=FALSE)
write.table(regen_amar, "results/tables/regen_candidates_Amar.tsv", sep="\t", row.names=FALSE)

# Only potentially regenerative domains heatmap

make_regen_heatmap <- function(regen_genes, species_df, species_name,
                               meta, min_tpm = 0.5) {
  regen_tpm <- tpm %>%
    filter(geneID %in% regen_genes$geneID)
  
  max_tpm <- apply(regen_tpm[, segment_labels], 1, max, na.rm = TRUE)
  regen_tpm <- regen_tpm[max_tpm >= min_tpm, ]
  
  if (nrow(regen_tpm) < 2) {
    message("Мало генов для регенерационного хитмапа: ", species_name)
    return(invisible(NULL))
  }
  
  mat <- as.matrix(regen_tpm[, segment_labels])
  rownames(mat) <- regen_tpm$geneID
  mat_z <- row_zscore(mat)
  
  set.seed(42)
  n_cl_regen <- min(6, nrow(mat_z) - 1)
  km <- kmeans(mat_z, centers = n_cl_regen, nstart = 25)
  
  peak_pos <- sapply(1:n_cl_regen, function(k)
    which.max(colMeans(mat_z[km$cluster == k, , drop = FALSE], na.rm = TRUE))
  )
  
  cmap  <- setNames(rank(peak_pos, ties.method = "first"), 1:n_cl_regen)
  cl_ord <- cmap[km$cluster]
  names(cl_ord) <- rownames(mat_z)
  
  gene_peak  <- apply(mat_z, 1, which.max)
  gene_order <- order(cl_ord, gene_peak)
  
  mat_sorted <- mat_z[gene_order, ]
  cl_sorted  <- cl_ord[gene_order]
  
  row_info <- data.frame(
    Cluster = factor(paste0("C", cl_sorted)),
    row.names = rownames(mat_sorted)
  )
  
  cl_cols <- setNames(
    colorRampPalette(
      c("#d73027","#fdae61","#abd9e9","#4575b4","#762a83","#1b7837")
    )(n_cl_regen),
    paste0("C", 1:n_cl_regen)
  )
  
  annot_colors <- c(meta$colors, list(Cluster = cl_cols))
  gaps <- cumsum(as.numeric(table(cl_sorted)))[-n_cl_regen]
  
  make_heatmap(
    mat = mat_sorted,
    title = paste("Regeneration associated domains", species_name),
    col_annot = meta$annot,
    annot_colors = annot_colors,
    row_annot = row_info,
    show_rownames = FALSE,
    row_clust = FALSE,
    gaps_row = gaps
  )
  
  invisible(list(mat = mat_sorted, clusters = cl_sorted))
}

regen_res_pele <- make_regen_heatmap(regen_pele, pele_full, "P.elegans",  species_meta$Pele)
regen_res_pdum <- make_regen_heatmap(regen_pdum, pdum_full, "P.dumerilii",species_meta$Pdum)
regen_res_amar <- make_regen_heatmap(regen_amar, amar_full, "A.marina",   species_meta$Amar)

# Regeneration-associated domains enrichment by clusters

extract_regen_cluster_annotation <- function(regen_res, species_df, species_name) {
  cluster_df <- data.frame(
    geneID  = names(regen_res$clusters),
    cluster = as.integer(regen_res$clusters)
  )
  
  annotated <- cluster_df %>%
    left_join(
      species_df %>% distinct(geneID, domain, pfam),
      by = "geneID"
    ) %>%
    arrange(cluster, geneID)
  
  domain_summary <- annotated %>%
    filter(!is.na(domain), domain != "") %>%
    group_by(cluster, domain) %>%
    summarise(n_genes = n(), .groups = "drop") %>%
    arrange(cluster, desc(n_genes))
  
  write.table(annotated,
              file = paste0("results/tables/regen_cluster_genes_", species_name, ".tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(domain_summary,
              file = paste0("results/tables/regen_cluster_domains_", species_name, ".tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  list(genes = annotated, domains = domain_summary)
}

plot_regen_cluster_enrichment <- function(regen_annot, species_name, top_n = 12) {
  df <- regen_annot$domains %>%
    filter(!is.na(domain), domain != "") %>%
    group_by(cluster) %>%
    slice_max(n_genes, n = top_n) %>%
    ungroup()
  
  ggplot(df, aes(x = reorder(domain, n_genes), y = n_genes,
                 fill = factor(cluster))) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    facet_wrap(~paste0("Cluster ", cluster), scales = "free_y", ncol = 2) +
    scale_fill_manual(
      values = colorRampPalette(
        c("#d73027", "#fdae61", "#abd9e9", "#4575b4", "#762a83", "#1b7837")
      )(length(unique(df$cluster)))
    ) +
    labs(
      title = paste("Regeneration-associated domain enrichment by cluster:", species_name),
      x = "Pfam domain", y = "Number of transcripts"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.text  = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 7),
      axis.title  = element_text(size = 9)
    )
}

rannot_pele <- extract_regen_cluster_annotation(regen_res_pele, pele_full, "Pele")
rannot_pdum <- extract_regen_cluster_annotation(regen_res_pdum, pdum_full, "Pdum")
rannot_amar <- extract_regen_cluster_annotation(regen_res_amar, amar_full, "Amar")

print(plot_regen_cluster_enrichment(rannot_pele, "P.elegans"))
print(plot_regen_cluster_enrichment(rannot_pdum, "P.dumerilii"))
print(plot_regen_cluster_enrichment(rannot_amar, "A.marina"))

# Comparative heatmap for regeneration-associated domains

make_comparative_regen_heatmap <- function(top_n = 100, min_genes = 2,
                                           title_suffix = "") {
  get_regen_domain_mat <- function(species_df) {
    species_df %>%
      filter(domain %in% common_domains) %>%
      filter(str_detect(domain, regen_pattern)) %>%
      filter(!is.na(domain), domain != "") %>%
      distinct(geneID, domain, .keep_all = TRUE) %>%
      group_by(domain) %>%
      summarise(across(all_of(segment_labels), mean, na.rm = TRUE),
                n_genes = n(), .groups = "drop") %>%
      filter(n_genes >= min_genes)
  }
  
  mat_p <- get_regen_domain_mat(pele_full)
  mat_d <- get_regen_domain_mat(pdum_full)
  mat_a <- get_regen_domain_mat(amar_full)
  
  shared <- Reduce(intersect, list(mat_p$domain, mat_d$domain, mat_a$domain))
  
  to_mat <- function(df) {
    m  <- df %>% filter(domain %in% shared) %>% arrange(domain)
    mx <- as.matrix(m[, segment_labels])
    rownames(mx) <- m$domain
    mx
  }
  
  m_p <- to_mat(mat_p)
  m_d <- to_mat(mat_d)
  m_a <- to_mat(mat_a)
  
  if (length(shared) > top_n) {
    mean_var <- (apply(m_p, 1, var, na.rm = TRUE) +
                   apply(m_d, 1, var, na.rm = TRUE) +
                   apply(m_a, 1, var, na.rm = TRUE)) / 3
    top_dom <- names(sort(mean_var, decreasing = TRUE))[1:top_n]
    m_p <- m_p[top_dom, ]
    m_d <- m_d[top_dom, ]
    m_a <- m_a[top_dom, ]
  }
  
  z_p <- row_zscore(m_p)
  z_d <- row_zscore(m_d)
  z_a <- row_zscore(m_a)
  
  colnames(z_p) <- paste0("Pele_", segment_labels)
  colnames(z_d) <- paste0("Pdum_", segment_labels)
  colnames(z_a) <- paste0("Amar_", segment_labels)
  
  spacer1  <- matrix(NA, nrow = nrow(z_p), ncol = 1,
                     dimnames = list(rownames(z_p), " "))
  spacer2  <- matrix(NA, nrow = nrow(z_p), ncol = 1,
                     dimnames = list(rownames(z_p), "  "))
  combined <- cbind(z_p, spacer1, z_d, spacer2, z_a)
  
  col_annot <- data.frame(
    Species = c(rep("P.elegans", 12), NA, rep("P.dumerilii", 12), NA, rep("A.marina", 12)),
    Region  = c(as.character(species_meta$Pele$annot$Region), NA,
                as.character(species_meta$Pdum$annot$Region), NA,
                as.character(species_meta$Amar$annot$Region)),
    row.names = colnames(combined)
  )
  
  all_regions <- c(species_meta$Pele$colors$Region,
                   species_meta$Pdum$colors$Region,
                   species_meta$Amar$colors$Region)
  all_regions <- all_regions[!duplicated(names(all_regions))]
  
  annot_colors <- list(
    Species = c("P.elegans"   = "#762a83",
                "P.dumerilii" = "#c51b7d",
                "A.marina"    = "#1b7837"),
    Region  = all_regions
  )

  row_clust <- hclust(dist(z_p))
  
  pheatmap(combined,
           color             = heatmap_colors,
           breaks            = heatmap_breaks,
           cluster_rows      = row_clust,
           cluster_cols      = FALSE,
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 7,
           fontsize_col      = 8,
           annotation_col    = col_annot,
           annotation_colors = annot_colors,
           na_col            = "white",
           main = paste0("Regeneration-associated domains: all species", title_suffix),
           border_color      = NA,
           treeheight_row    = 40)
}


make_comparative_regen_heatmap(top_n = 100, min_genes = 2)


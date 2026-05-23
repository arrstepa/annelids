library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(rstatix)


# Load data

pfam <- read.table("demodata/pfam_gene_clean.tsv", sep = "\t", stringsAsFactors = FALSE)
colnames(pfam) <- c("transcript", "pfam", "domain")

tpm <- read.table("demodata/all_tpm.txt", sep = "\t", header = TRUE)
colnames(tpm)[1] <- "transcript"


# Expression matrix

expr_cols <- grep("^X", colnames(tpm), value = TRUE)
expr_matrix <- as.matrix(tpm[, expr_cols])

positions <- seq_len(ncol(expr_matrix))

# Anterio-posterior pattern metrics

# Shannon entropy of the expression profile, normalized to [0, 1]:
### High entropy value means expression spread evenly across segments (smooth gradient)
### Low entropy - concentrated in few segments (localized)

calc_entropy <- function(x) {
  if (sum(x) == 0) return(NA)
  p <- x / sum(x)
  H <- -sum(p * log2(p + 1e-12))
  H / log2(length(x))
}


# Expression centroid: the "center of mass" position along the AP axis

calc_centroid <- function(x) {
  if (sum(x) == 0) return(NA)
  sum(positions * x) / sum(x)
}

# Spread: weighted standard deviation of position around the centroid

calc_spread <- function(x) {
  if (sum(x) == 0) return(NA)
  mu <- calc_centroid(x)
  sqrt(sum(x * (positions - mu)^2) / sum(x))
}

# Polarization: relative enrichment at the two body ends compared to the middle.

calc_polarization <- function(x) {
  anterior  <- mean(x[1:2])
  posterior <- mean(x[(length(x) - 1):length(x)])
  middle_idx <- floor(length(x) / 2)
  middle <- mean(x[(middle_idx - 1):(middle_idx + 1)])
  ((anterior + posterior) / 2 - middle) / (mean(x) + 1e-12)
}

# AP index: log2 ratio of anterior to posterior expression.
### > 0: anterior bias, 
### < 0: posterior bias.

calc_AP_index <- function(x) {
  anterior  <- mean(x[1:2])
  posterior <- mean(x[(length(x) - 1):length(x)])
  if (anterior == 0 & posterior == 0) return(0)
  log2((anterior + 1e-6) / (posterior + 1e-6))
}

# Total variation: sum of absolute changes between neighboring segments.
### High value means sharp transitions; 
### Low value means smooth gradients.

calc_total_variation <- function(x) {
  sum(abs(diff(x)))
}

# Number of local maxima in the expression profile.

calc_num_peaks <- function(x) {
  n <- length(x)
  sum(x[2:(n - 1)] > x[1:(n - 2)] & x[2:(n - 1)] > x[3:n])
}



# Metrics for all transcripts

metrics <- tpm %>%
  mutate(
    max_tpm = max_tpm,
    species          = sub("_.*", "", transcript),
    variance         = apply(expr_matrix, 1, var),
    sd               = apply(expr_matrix, 1, sd),
    mean_expr        = rowMeans(expr_matrix),
    cv               = sd / (mean_expr + 1e-12),
    entropy          = apply(expr_matrix, 1, calc_entropy),
    centroid         = apply(expr_matrix, 1, calc_centroid),
    spread           = apply(expr_matrix, 1, calc_spread),
    polarization     = apply(expr_matrix, 1, calc_polarization),
    AP_index         = apply(expr_matrix, 1, calc_AP_index),
    total_variation  = apply(expr_matrix, 1, calc_total_variation),
    n_peaks          = apply(expr_matrix, 1, calc_num_peaks)
  )

# Drop very lowly expressed transcripts as noise
metrics <- metrics %>% filter(max_tpm > 1)

# Join metrics with domain annotation (one row per transcript-domain pair)

df_with_tpm <- inner_join(metrics, pfam, by = "transcript")
df_with_tpm$species <- sub("_.*", "", df_with_tpm$transcript)



species_levels <- c("Pele", "Pdum", "Amar")  # we put species here in order of regenerative potential

metrics$species <- factor(metrics$species, levels = species_levels)
df_with_tpm$species <- factor(df_with_tpm$species, levels = species_levels)


# Metrics distributions for all transcripts

kruskal.test(entropy ~ species, data = metrics)
pairwise.wilcox.test(metrics$entropy, metrics$species, p.adjust.method = "BH")

metrics %>% wilcox_effsize(entropy ~ species,
                 comparisons = list(c("Amar", "Pele"), c("Pdum", "Pele")))

metrics %>% wilcox_effsize(AP_index ~ species,
                 comparisons = list(c("Amar", "Pele"), c("Pdum", "Pele")))


# Regeneration-associated transcripts and their metrics

regulatory_keywords <- unique(c(
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

regulatory_pattern <- paste(regulatory_keywords, collapse = "|")

regulatory_ids <- df_with_tpm %>%
  filter(grepl(regulatory_pattern, domain, ignore.case = TRUE)) %>%
  pull(transcript) %>%
  unique()

metrics$reg_group <- ifelse(
  metrics$transcript %in% regulatory_ids,
  "Regulatory",
  "Background"
)

metrics$reg_group <- factor(metrics$reg_group, levels = c("Regulatory", "Background"))

print(table(metrics$species, metrics$reg_group))

# Regulatory-associated transcripts are less abundant and we must balance them:


set.seed(42)

balanced <- metrics %>%
  group_by(species) %>%
  group_modify(~{
    reg <- .x %>% filter(reg_group == "Regulatory")
    bg  <- .x %>% filter(reg_group == "Background")
    if (nrow(reg) == 0 || nrow(bg) == 0) return(reg[0, ])
    bg_sub <- bg %>% sample_n(min(nrow(reg), nrow(bg)))
    bind_rows(reg, bg_sub)
  }) %>%
  ungroup()

# Violin + boxplot of a metric, split by regulatory/background, faceted by species

plot_group_metric <- function(data, metric, group_col) {
  ggplot(data, aes(x = .data[[group_col]], y = .data[[metric]], fill = .data[[group_col]])) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    facet_wrap(~species) +
    labs(y = metric, x = NULL) +
    theme_bw()
}

print(plot_group_metric(balanced, "entropy", "reg_group"))
print(plot_group_metric(balanced, "total_variation", "reg_group"))
print(plot_group_metric(balanced, "AP_index", "reg_group"))
print(plot_group_metric(balanced, "polarization", "reg_group"))


core_metrics <- c("entropy", "total_variation", "AP_index", "polarization")

# Wilcoxon test

regulatory_stats <- lapply(core_metrics, function(metric) {
  balanced %>%
    group_by(species) %>%
    wilcox_test(as.formula(paste(metric, "~ reg_group"))) %>%
    mutate(metric = metric)
}) %>% bind_rows()

# Wilcoxon test + effective size correction - due to the fact that there are too many samples and p-value will be extremely small

regulatory_effsize <- lapply(core_metrics, function(metric) {
  balanced %>%
    group_by(species) %>%
    wilcox_effsize(as.formula(paste(metric, "~ reg_group"))) %>%
    mutate(metric = metric)
}) %>% bind_rows()


AP_INDEX_QUANTILE <- 0.90  # treshold for highly polarized
MIN_DOMAIN_COUNT <- 20    # in order to skip rare domains we filter them

# anterior-biased and posterior-biased transcripts within each species

df_with_tpm <- df_with_tpm %>%
  group_by(species) %>%
  mutate(ap_threshold = quantile(AP_index, AP_INDEX_QUANTILE, na.rm = TRUE),
    high_pol = AP_index > ap_threshold
  ) %>%
  ungroup()

# Enrichment test 

enrichment_for_species <- function(df_sp) {
  domains <- unique(df_sp$domain)
  results <- list()
  
  for (d in domains) {
    has_domain <- df_sp$domain == d
    if (sum(has_domain) < MIN_DOMAIN_COUNT) next
    
    tab <- table(has_domain, df_sp$high_pol)
    if (nrow(tab) < 2 || ncol(tab) < 2) next
    
    ft <- fisher.test(tab)
    results[[d]] <- data.frame(
      species    = df_sp$species[1],
      domain     = d,
      odds_ratio = unname(ft$estimate),
      p_value    = ft$p.value,
      conf_low   = ft$conf.int[1],
      conf_high  = ft$conf.int[2],
      n          = sum(has_domain)
    )
  }
  
  res <- do.call(rbind, results)
  if (!is.null(res)) res$padj <- p.adjust(res$p_value, method = "BH")
  res
}

enrichment <- df_with_tpm %>%
  group_split(species) %>%
  lapply(enrichment_for_species) %>%
  bind_rows()

enrichment$species <- factor(enrichment$species, levels = species_levels)

# Significant domains

sig <- enrichment %>% filter(padj < 0.05)

# top enriched or depleted domains
sig %>% arrange(desc(odds_ratio)) %>% head(15) %>% print()
sig %>% arrange(odds_ratio) %>% head(15) %>% print()

# Volcano plot per species: right side enriched among anterior-biased, left side - among posterior-biased

top <- sig %>% filter(abs(log2(odds_ratio)) > 1)

ggplot(sig, aes(x = log2(odds_ratio), y = -log10(padj))) +
  geom_point(alpha = 0.7) +
  geom_text_repel(data = top, aes(label = domain), size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(
    ~species,
    scales = "free_y",
    ncol = 1,
    labeller = as_labeller(c(
      Pele = "italic('P. elegans')",
      Pdum = "italic('P. dumerilii')",
      Amar = "italic('A. marina')"
    ), label_parsed)
  ) +
  labs(x = "log2(odds ratio)", y = "-log10(adjusted p)") +
  theme_bw(base_size = 14)


# Now we can also calculate fractions of anterior/posterior and bipolar expressed regeneration-associated transcriptes in all species

regen_tpm <- metrics %>% filter(reg_group == "Regulatory")

regen_tpm$species <- factor(regen_tpm$species, levels = c("Pele", "Pdum", "Amar"))

# Mean expression at body ends

regen_tpm$anterior <- apply(regen_tpm[, expr_cols], 1, function(x) mean(x[1:2]))

regen_tpm$posterior <- apply(regen_tpm[, expr_cols], 1, function(x) mean(x[(length(x)-1):length(x)]))

# Relative expression normalized by maximal TPM

regen_tpm$anterior_rel <- regen_tpm$anterior / (regen_tpm$max_tpm + 1e-12)

regen_tpm$posterior_rel <- regen_tpm$posterior / (regen_tpm$max_tpm + 1e-12)

# Spatial classification threshold

thr <- 0.5

regen_tpm$pattern <- "other"

regen_tpm$pattern[regen_tpm$anterior_rel > thr & regen_tpm$posterior_rel > thr] <- "bipolar"

regen_tpm$pattern[regen_tpm$anterior_rel > thr & regen_tpm$posterior_rel <= thr] <- "anterior_only"

regen_tpm$pattern[regen_tpm$anterior_rel <= thr & regen_tpm$posterior_rel > thr] <- "posterior_only"

# Fraction of spatial expression classes

plot_df <- as.data.frame(prop.table(table(regen_tpm$species, regen_tpm$pattern), margin = 1))

colnames(plot_df) <- c(
  "species",
  "pattern",
  "fraction"
)

plot_df$pattern <- factor(
  plot_df$pattern,
  levels = c(
    "anterior_only",
    "posterior_only",
    "bipolar",
    "other"
  )
)

# Stacked barplot

ggplot(
  plot_df,
  aes(
    x = species,
    y = fraction,
    fill = pattern
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  
  scale_x_discrete(
    labels = c(
      Pele = "P. elegans",
      Pdum = "P. dumerilii",
      Amar = "A. marina"
    )
  ) +
  
  scale_y_continuous(
    expand = c(0, 0)
  ) +
  
  labs(
    title = "Spatial classes of regulatory transcripts",
    x = "Species",
    y = "Fraction of transcripts",
    fill = "Expression pattern"
  ) +
  
  theme_bw(base_size = 13)



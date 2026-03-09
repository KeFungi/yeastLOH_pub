library(tidyverse)
library(ggplot2)

# --- 1. Get List of Strains from meta_table.csv ---
meta_file <- "tables/meta_table.csv"
meta_data <- read_csv(meta_file, show_col_types = FALSE)

# Filter for strains included in analyses
good_strains <- meta_data %>%
    filter(`included in analyses` == "YES") %>%
    filter(!str_detect(strain, "^CNTRL-")) %>%
    select(strain) %>%
    distinct()

cat("Identified", nrow(good_strains), "strains included in analyses.\n")

# --- 2. Load BED Files ---
# Define file paths
loh_bed_path <- "LOH_detect/LOH_minSNP-5_typed.bed"
rev_loh_bed_path <- "LOH_detect/revLOH_minSNP-5_typed.bed"

# Read LOH (SK1)
loh_data <- read_tsv(loh_bed_path, col_names = FALSE, show_col_types = FALSE) %>%
    select(chr = X1, start = X2, end = X3, strain = X4) %>%
    mutate(type = "SK1")

# Read revLOH (BY4741)
rev_loh_data <- read_tsv(rev_loh_bed_path, col_names = FALSE, show_col_types = FALSE) %>%
    select(chr = X1, start = X2, end = X3, strain = X4) %>%
    mutate(type = "BY4741")

# Combine and Filter
# Join with the good_strains list
all_loh_data <- bind_rows(loh_data, rev_loh_data) %>%
    inner_join(good_strains, by = "strain")

# --- 3. Load Chromosome Sizes for Continuous Plotting ---
chr_sizes <- read_tsv("ref/S288C.chr.fasta.gz.fai", col_names = FALSE, show_col_types = FALSE) %>%
    select(chr = X1, length = X2)

# Ensure chromosomes are ordered correctly
chr_order <- c(
    "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", "chrVIII",
    "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI", "Mito"
)

chr_sizes <- chr_sizes %>%
    filter(chr %in% chr_order) %>%
    mutate(chr = factor(chr, levels = chr_order)) %>%
    arrange(chr) %>%
    mutate(
        cum_start = cumsum(lag(length, default = 0)),
        cum_end = cum_start + length,
        label_pos = (cum_start + cum_end) / 2
    )

# --- 4. Transform Coordinates for Continuous Plotting ---
all_loh_data <- all_loh_data %>%
    filter(chr %in% chr_order) %>%
    left_join(select(chr_sizes, chr, cum_start), by = "chr") %>%
    mutate(
        start_cum = start + cum_start,
        end_cum = end + cum_start
    )

# --- 5. Plotting ---
# Y-axis: Strains (Alphabetical Order)
# X-axis: Continuous Genome Position

# Order strains alphabetically
strain_levels <-
    good_strains %>%
    filter(strain %in% all_loh_data$strain) %>%
    arrange(strain) %>%
    pull(strain)

all_loh_data <- all_loh_data %>%
    mutate(strain = factor(strain, levels = strain_levels))

# Create backbone data
backbone_data <- crossing(chr_sizes, strain = factor(strain_levels, levels = strain_levels))

# Create detailed plot
p <- ggplot() +
    # Draw chromosome backbone bars (all grey)
    geom_rect(
        data = backbone_data,
        aes(
            xmin = cum_start, xmax = cum_end,
            ymin = as.numeric(strain) - 0.45,
            ymax = as.numeric(strain) + 0.45
        ),
        fill = "grey90", color = NA
    ) +

    # Draw LOH regions (SK1 and BY4741)
    geom_rect(
        data = all_loh_data,
        aes(
            xmin = start_cum, xmax = end_cum,
            ymin = as.numeric(strain) - 0.45,
            ymax = as.numeric(strain) + 0.45,
            fill = type
        ), color = NA
    ) +

    # Custom Colors
    scale_fill_manual(values = c("SK1" = "#D55E00", "BY4741" = "#0072B2")) +

    # Scales and Theme
    scale_x_continuous(breaks = chr_sizes$label_pos, labels = chr_sizes$chr, expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(strain_levels), labels = strain_levels, expand = c(0, 0.5)) +
    theme_bw() +
    theme(
        panel.background = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_blank(),
        legend.position = "bottom"
    ) +
    # Add vertical lines at chromosome boundaries
    geom_vline(xintercept = chr_sizes$cum_end, color = "black", linetype = "dotted", linewidth = 0.3) +
    labs(
        title = NULL,
        x = "Chromosome",
        y = "Strain",
        fill = "Genotype"
    )

# --- 6. Save Plot ---
# Calculate height based on number of strains
plot_height <- max(10, nrow(good_strains) * 0.25) # Ensure minimum height and scale for many strains

output_file <- "plots/LOH_track_all_good_strains.pdf"
ggsave(output_file, plot = p, width = 16, height = plot_height, limitsize = FALSE)

cat("All Strains Plot saved to", output_file, "\n")

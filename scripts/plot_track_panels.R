library(tidyverse)
library(ggplot2)

# --- 1. Get List of Strains from top5_list.tsv ---
# Read strains and their treatments directly from the summary TSV
strain_list <- read_tsv("tables/top5_list.tsv", show_col_types = FALSE) %>%
    select(strain, treatment = top5_in) %>%
    filter(!is.na(strain) & !is.na(treatment)) %>%
    # Format treatments nicely for the plot labels (e.g. "Acidic, Caffeine" -> "Acidic &\n Caffeine")
    mutate(treatment = str_replace_all(treatment, ", ", " &\n"))

cat("Identified", nrow(strain_list), "unique strains in", length(unique(strain_list$treatment)), "combined treatments.\n")

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
# Join with the consolidated strain list (1:many relationship strain -> bed regions)
all_loh_data <- bind_rows(loh_data, rev_loh_data) %>%
    inner_join(strain_list, by = "strain")

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

# Add cumulative start offsets to LOH data
all_loh_data <- all_loh_data %>%
    filter(chr %in% chr_order) %>%
    left_join(select(chr_sizes, chr, cum_start), by = "chr") %>%
    mutate(
        start_cum = start + cum_start,
        end_cum = end + cum_start
    )


# --- 5. Plotting ---
# We want to plot tracks continuously.
# Y-axis: Strains (grouped by Treatment)
# X-axis: Continuous Genome Position

# Create Unique Ids for Strains per Treatment
# Now treatment is the combined string, so uniqueness is guaranteed by strain name alone within treatment.
# But we still use the helper to ensure factor order is by Treatment then Strain.

strain_ordering <- strain_list %>%
    arrange(treatment, strain) %>%
    mutate(unique_id = paste(treatment, strain, sep = "__")) %>%
    mutate(unique_id = factor(unique_id, levels = unique_id))

# Map unique_id back to data
all_loh_data <- all_loh_data %>%
    mutate(unique_id = factor(paste(treatment, strain, sep = "__"), levels = levels(strain_ordering$unique_id)))

# Create backbone data with unique_ids
backbone_data <- crossing(chr_sizes,
    unique_id = levels(strain_ordering$unique_id)
) %>%
    mutate(unique_id = factor(unique_id, levels = levels(strain_ordering$unique_id))) %>%
    left_join(strain_ordering, by = "unique_id") # add treatment info for faceting

# Helper for Y-axis labels
clean_labels <- function(x) {
    str_match(x, "__(.+$)")[, 2]
}

# Create detailed plot
p <- ggplot() +
    # Draw chromosome backbone bars (all grey)
    geom_rect(
        data = backbone_data,
        aes(
            xmin = cum_start, xmax = cum_end,
            ymin = as.numeric(unique_id) - 0.4,
            ymax = as.numeric(unique_id) + 0.4
        ),
        fill = "grey90", color = NA
    ) +

    # Draw LOH regions (SK1 and BY4741)
    geom_rect(
        data = all_loh_data,
        aes(
            xmin = start_cum, xmax = end_cum,
            ymin = as.numeric(unique_id) - 0.4,
            ymax = as.numeric(unique_id) + 0.4,
            fill = type
        ), color = NA
    ) +

    # Facet by Treatment (vertical layout)
    facet_grid(treatment ~ ., scales = "free_y", space = "free_y") +

    # Custom Colors
    scale_fill_manual(values = c("SK1" = "#D55E00", "BY4741" = "#0072B2")) +

    # Scales and Theme
    scale_x_continuous(breaks = chr_sizes$label_pos, labels = chr_sizes$chr, expand = c(0, 0)) +
    scale_y_continuous(
        breaks = seq_along(levels(strain_ordering$unique_id)),
        labels = clean_labels(levels(strain_ordering$unique_id)),
        expand = c(0, 0.5)
    ) +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = "grey95"),
        strip.text.x = element_text(size = 14, angle = 0),
        strip.text.y = element_text(size = 14, angle = 0),
        axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
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
output_file <- "plots/LOH_track_panels.pdf"
ggsave(output_file, plot = p, width = 16, height = 10, limitsize = FALSE)

cat("Combined Treatment Plot saved to", output_file, "\n")

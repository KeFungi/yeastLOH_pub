#### Library ####
library(ggpubr)
library(rstatix)
library(plyranges)
library(tidytext)
library(khroma)
library(smplot2)
library(ggpmisc)
library(cowplot)
library(pracma)
library(tidyverse)
library(qtl)
library(qtl2)
library(qtl2convert)
library(multcompView)
library(GGally)

#### helper functions ####
rc <-
  function(seq){
    seq %>%
      strsplit("") %>%
      .[[1]] %>%
      rev() %>%
      chartr("ATGC","TACG",.) %>%
      paste0(collapse = "")
  }

inlier_mean <-
  function(x, IQR=1.5){
    quantile1_x <- quantile(x, 0.25)
    quantile3_x <- quantile(x, 0.75)
    IQR_x <- quantile3_x - quantile1_x
    mean(x[x>(quantile1_x-IQR_x*IQR)&x<(quantile3_x+IQR_x*IQR)])
  }

inlier_sd <-
  function(x, IQR=1.5){
    quantile1_x <- quantile(x, 0.25)
    quantile3_x <- quantile(x, 0.75)
    IQR_x <- quantile3_x - quantile1_x
    sd(x[x>(quantile1_x-IQR_x*IQR)&x<(quantile3_x+IQR_x*IQR)])
  }


ggplot_peaks <-
  function(peaks, map, tick_height,
           gap, bgcolor, algbgcolor,
           lwd=2, col="slateblue", xlab=NULL, ylab="",
           xlim=NULL, ylim=NULL, xaxs="i", yaxs="i",
           main="", mgp.x=c(2.6, 0.5, 0), mgp.y=c(2.6, 0.5, 0),
           mgp=NULL, las=1, lend=1, ljoin=1,
           hlines=NULL, hlines.col="white", hlines.lwd=1, hlines.lty=1,
           vlines=NULL, vlines.col="white", vlines.lwd=1, vlines.lty=1,
           point_size=0, vbars=TRUE,
           xaxt=ifelse(onechr, "y", "n"),
           yaxt="y",
           legend.title = "",
           legend.position = "right",
           ...) {
    dots <- list(...)
    onechr <- (length(map)==1) # single chromosome

    # make chr into factor
    chrs <- names(map)
    peaks$chr <- factor(peaks$chr, chrs)
    peaks$lodcolumn <- factor(peaks$lodcolumn)

    # color
    if(is.null(peaks$col)) {
      peaks$col <- "all"
      legend.position <- "none"
    }

    # get map lengths by chr.
    mapl <- t(sapply(map[chrs], range))
    mapl <- data.frame(mapl)
    colnames(mapl) <- c("lo","hi")
    mapl$chr <- factor(chrs, chrs)

    if(is.null(xlab)) {
      if(onechr) {
        if(names(map) == " ") xlab <- "Position"
        else xlab <- paste("Chr", names(map), "position")
      }
      else xlab <- "Chromosome"
    }

    p <- ggplot2::ggplot(peaks) +
      ggplot2::aes(x = .data$pos, y = .data$lodcolumn, group = .data$chr) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab)

    # Add color and legend.
    p <- p +
      ggplot2::theme(legend.position = legend.position)

    if(main != "") {
      p <- p +
        ggplot2::ggtitle(main)
    }

    if(!onechr) {
      p <- p +
        ggplot2::facet_grid(~ .data$chr, scales = "free_x", space = "free_x") +
        ggplot2::theme(panel.spacing = grid::unit(gap / 10000, "npc"))
    }

    # set up horizontal axis to match data.
    p <- p

    # include axis labels?
    if(yaxt == "n") {
      p <- p +
        ggplot2::theme(
          axis.text.y  = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank())
    }
    # X axis
    if(xaxt == "n") {
      p <- p +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank())
    }

    # grid lines
    if((length(vlines)==1 && is.na(vlines)) | !onechr) { # if vlines==NA (or mult chr), skip lines
      p <- p +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank())
    }
    if((length(hlines)==1 && is.na(hlines))) { # if hlines==NA, skip lines
      p <- p +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor.y = ggplot2::element_blank())
    }

    # Add box for each chr.
    p <- p +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = "black",
                                             fill=NA))

    p
  }

#### data input ####
input_tb <-
  read_csv("tables/bartender_matrix.csv") %>%
  rowwise() %>%
  mutate(barcode=rc(barcode))

strain_tb <-
  read_csv("tables/project_id.csv") %>%
  group_by(Description) %>%
  slice(1) %>%
  select(strain=Description) %>%
  ungroup()

control_strains <-
  filter(strain_tb, str_detect(strain, "^CNTRL-")) %>%
  pull(strain)

bc_tb <-
  read_csv("pub_tables/barcode_map.csv", col_names = c("strain", "barcode")) %>%
  left_join(strain_tb, .) %>%
  mutate(`barcode length`=nchar(barcode))

meta_tb <-
  read_csv("tables/meta_table.csv")

good_strains <-
  meta_tb %>%
  filter(`included in analyses`=="YES") %>%
  pull(strain)

heterozygosity_tb <-
    read_tsv("tables/runs.GTfilter.selected.rsid.vcf.gz.het") %>%
  mutate(heterozygosity_SNPs=(N_SITES-`O(HOM)`)/N_SITES,
         heterozygosity_genome=(N_SITES-`O(HOM)`)/12071326) %>%
  select(strain=INDV, heterozygosity_SNPs, heterozygosity_genome)

fordLOH_tb <-
  read_tsv("LOH_detect/LOH_minSNP-5_typed.bed", col_names = FALSE) %>%
  mutate(type="SK1") %>%
  filter(X4 %in% good_strains) %>%
  rename(LOH_type=X5)

revLOH_tb <-
  read_tsv("LOH_detect/revLOH_minSNP-5_typed.bed", col_names = FALSE) %>%
  mutate(type="BY4741") %>%
  filter(X4 %in% good_strains) %>%
  rename(LOH_type=X5)

LOH_type_tb <-
  rbind(fordLOH_tb, revLOH_tb) %>%
  rename(chr=X1, start=X2, end=X3, strain=X4) %>%
  mutate(length=end-start) %>%
  filter(strain %in% good_strains)

LOH_sep_tb <-
  LOH_type_tb %>%
  mutate(strain=factor(strain, levels=good_strains)) %>%
  group_by(strain, type, LOH_type) %>%
  summarise(n_LOH=n(), LOH_length=sum(length)) %>%
  ungroup() %>%
  complete(strain, type, LOH_type, fill = list(n_LOH=0, LOH_length=0))

LOH_long_tb <-
  LOH_sep_tb %>%
  group_by(strain, type, LOH_type) %>%
  summarise(n_LOH=sum(n_LOH), LOH_length=sum(LOH_length)) %>%
  mutate(type="both") %>%
  rbind(LOH_sep_tb, .) %>%
  complete(strain, type, LOH_type, fill = list(n_LOH=0, LOH_length=0)) %>%
  arrange(strain, type, LOH_type)

LOH_wide_tb <-
  LOH_long_tb %>%
  select(-LOH_type) %>%
  group_by(strain, type) %>%
  summarise(n_LOH=sum(n_LOH), LOH_length=sum(LOH_length)) %>%
  ungroup() %>%
  pivot_wider(id_cols = "strain", names_from = "type", values_from = c("n_LOH", "LOH_length"), names_sep = "-")

LOH_wide_wide_tb <-
  LOH_long_tb %>%
  group_by(strain, type, LOH_type) %>%
  summarise(n_LOH=sum(n_LOH), LOH_length=sum(LOH_length)) %>%
  pivot_wider(id_cols = "strain", names_from = c("type", "LOH_type"), values_from = c("n_LOH", "LOH_length"), names_sep = "-")

bedgraph_1 <-
  read_bed_graph("LOH_detect/any.depth.bedgraph") %>%
  as_tibble() %>%
  mutate(type="any")

bedgraph_2 <-
  read_bed_graph("LOH_detect/forward.depth.bedgraph") %>%
  as_tibble() %>%
  mutate(type="SK1")

bedgraph_3 <-
  read_bed_graph("LOH_detect/reverse.depth.bedgraph") %>%
  as_tibble() %>%
  mutate(type="BY4741")

chr_len_tb <-
  read_tsv("ref/S288C.chr.fasta.gz.fai", col_names = FALSE)

genome_length <-
  sum(chr_len_tb$X2)

SNP_loc <-
  read_tsv("bwa_haplotypecaller_finalvcf/runs.diploid.vcf.tsv.gz") %>%
  select(chr=`#[1]CHROM`, ci_lo=`[2]POS`) %>%
  mutate(lodcolumn=0.5, pos=ci_lo, ci_hi=ci_lo+200, lodindex="depth", col="transparent", chr=str_remove(chr, "chr")) %>%
  mutate(chr=factor(chr, levels = unique(chr)))

cen_loc <-
  read_bed("ref/centromere.bed") %>%
  as_tibble() %>%
  unique() %>%
  select(chr=seqnames, ci_lo=start, ci_hi=end) %>%
  mutate(lodcolumn=-1, pos=ci_lo, lodindex="depth", col="transparent", chr=str_remove(chr, "chr")) %>%
  mutate(lodindex="depth", col="transparent", chr=str_remove(chr, "chr")) %>%
  mutate(chr=factor(chr, levels = chr))

exp_treatments <-
  c("SC", "YPD-30C", "YPD-37C", "Acidic", "Caffeine", "H2O2", "worm-20C")


exp_stat_treatments <-
  c("SC", "YPD-30C", "YPD-37C", "Acidic", "Caffeine", "H2O2", "worm-20C", "(average)", "(maximum s)", "(minimum s)")

#### calculate fitness ####
#### initial freq ####
good_bc_tb <-
  bc_tb %>%
  filter(strain %in% c(good_strains, "P3-2C")) %>%
  inner_join(input_tb) %>%
  select(-barcode,-`barcode length`) %>%
  relocate(strain)

barseq_long_tb <-
  good_bc_tb %>%
  pivot_longer(-1, values_to = "count") %>%
  group_by(name) %>%
  mutate(freq=count/sum(count, na.rm=TRUE)) %>%
  select(-count) %>%
  mutate(name=str_remove(name, "-after"))

media_pool_tb <-
  barseq_long_tb %>%
  filter(str_detect(name, "^LOH-pool-")) %>%
  group_by(strain) %>%
  summarise(
    mean_freq_withOutlier=mean(freq, na.rm=TRUE),
    mean_freq=inlier_mean(freq),
    median_freq=median(freq, na.rm=TRUE)
  ) %>%
  arrange(mean_freq)

worm_pool_tb <-
  barseq_long_tb %>%
  filter(str_detect(name, "^Worm-pool-")) %>%
  group_by(strain) %>%
  summarise(
    mean_freq_withOutlier=mean(freq, na.rm=TRUE),
    mean_freq=inlier_mean(freq),
    median_freq=median(freq, na.rm=TRUE)
  ) %>%
  arrange(mean_freq)

barseq_long_tb %>%
  filter(str_detect(name, "-pool-")) %>%
  rowwise() %>%
  mutate(type=case_when(
    str_detect(name, "^Worm-pool-") ~ "worm",
    str_detect(name, "^LOH-pool-") ~ "media",
    TRUE ~ NA)
  ) %>%
  arrange(freq) %>%
  mutate(treatment=str_split(name, "-[^-]+$", simplify = TRUE)[[1]]) %>%
  ggplot(aes(x=reorder_within(strain, freq, type, median), y=freq)) +
  geom_boxplot() +
  facet_wrap(~type, scales = "free_x") +
  theme(axis.text.x = element_blank()) +
  geom_hline(yintercept = 0.0001, color="red") +
  xlab("strain") +
  ylab("frequency in inoculum")

ggsave("plots/beforeFreq.pdf", width = 12, height = 7)

#### after freq and fitness ####
media_tb <-
  barseq_long_tb %>%
  filter(!str_detect(name, "worm")) %>%
  filter(!str_detect(name, "-pool-")) %>%
  rename(after_freq=freq) %>%
  right_join(rename(media_pool_tb, before_freq=mean_freq)) %>%
  filter(before_freq > 0.0001) %>%
  mutate(enrichment=after_freq/before_freq)

worm_tb <-
  barseq_long_tb %>%
  filter(str_detect(name, "worm")) %>%
  filter(!str_detect(name, "-pool-")) %>%
  filter(!str_detect(name, "37C")) %>%
  rename(after_freq=freq) %>%
  right_join(rename(worm_pool_tb, before_freq=mean_freq)) %>%
  filter(before_freq > 0.0001) %>%
  mutate(enrichment=after_freq/before_freq)

control_tb <-
  rbind(media_tb, worm_tb) %>%
  filter(strain %in% control_strains) %>%
  group_by(name) %>%
  summarise(
    control_enrichment=inlier_mean(enrichment)
  )

exp_tb <-
  rbind(media_tb, worm_tb) %>%
  left_join(control_tb) %>%
  ungroup() %>%
  mutate(s=0.1*log(enrichment/control_enrichment)) %>%
  mutate(treatment=str_split(name, "-[^-]+$", simplify = TRUE)[,1]) %>%
  mutate(treatment=factor(treatment, levels=exp_stat_treatments)) %>%
  filter(!is.na(treatment)) %>%
  mutate(replicate=str_match(name, "[^-]+$")[,1]) %>%
  select(strain, experiment=name, treatment, replicate, before_freq, after_freq, enrichment, control_enrichment, s)

exp_tb %>%
  left_join(select(LOH_wide_tb, strain, `n_LOH-both`)) %>%
  filter(!treatment=="YPD-30C-35cyc") %>%
  group_by(treatment, strain) %>%
  mutate(quantile1_s=quantile(s, 0.25),
         quantile3_s=quantile(s, 0.75),
         IQR_s=quantile3_s - quantile1_s) %>%
  ungroup() %>%
  filter(s>quantile1_s-IQR_s*1.5, s<quantile3_s+IQR_s*1.5) %>%
  group_by(treatment) %>%
  summarise(lm=list(lm(s ~ strain ))) %>%
  mutate(anova=map(lm, ~anova(.x))) %>%
  mutate(anova=map(anova, ~mutate(.x, percent=`Sum Sq`/sum(`Sum Sq`)))) %>%
  mutate(anova=map(anova, function(x){x["sum",c("Df", "Sum Sq")] <- c(sum(x[,"Df"]), sum(x[,"Sum Sq"])); return(x)})) %>%
  select(-lm) %>%
  mutate(anova=map(anova, ~rownames_to_column(as.data.frame(.x), "Factor"))) %>%
  unnest(anova) %>%
  write_csv("tables/anova.csv")

control_E_plot1 <-
  exp_tb %>%
  filter(strain %in% control_strains) %>%
  mutate(treatment=fct_relevel(treatment, "worm-20C", after=Inf)) %>%
  filter(treatment!="YPD-30C-35cyc") %>%
  ggplot(aes(x=treatment, y=log2(enrichment),  color=strain)) +
  geom_abline(slope = 0, intercept = 0, lty=2, alpha=0.5) +
  geom_boxplot(outliers = FALSE) +
  geom_point(size=0.5, position = position_dodge(width=0.75)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

control_E_plot2 <-
  exp_tb %>%
  filter(strain %in% control_strains) %>%
  ggplot(aes(x=treatment, y=log2(enrichment),  color=replicate)) +
  geom_abline(slope = 0, intercept = 0, lty=2, alpha=0.5) +
  geom_boxplot(outliers = FALSE) +
  geom_point(size=0.5, position = position_dodge(width=0.75)) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle=30, hjust = 1))

plot_grid(control_E_plot1, control_E_plot2, ncol = 1)

ggsave("plots/control_enrichement.pdf", width = 8, height = 10)

exp_tb %>%
  left_join(LOH_wide_tb) %>%
  filter(`n_LOH-both`>0) %>%
  group_by(treatment, strain) %>%
  mutate(mean_s=inlier_mean(s)) %>%
  ggplot(aes(x=reorder_within(strain, s, treatment, inlier_mean), y=s)) +
  geom_boxplot() +
  geom_point(aes(y=mean_s), color="red", shape='-') +
  scale_x_reordered() +
  geom_hline(yintercept = 0, lty=2, alpha=0.5) +
  facet_wrap(~treatment, ncol=2, scales = "free") +
  theme_linedraw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(minor_breaks = NULL) +
  xlab("strain") +
  ylab("selection coefficient (s)")

ggsave("plots/fitness_box.pdf", width = 18, height = 8)

fit_exp_tb <-
  exp_tb %>%
  ungroup() %>%
  group_by(strain, treatment) %>%
  summarise(
    mean_s_withOutlier=mean(s, na.rm=TRUE),
    mean_s=inlier_mean(s),
    median_s=median(s, na.rm=TRUE)
  )

fit_stat_tb <-
  fit_exp_tb %>%
  group_by(strain) %>%
  summarise("(average)"=mean(mean_s),
            "(maximum s)"=max(mean_s),
            "(minimum s)"=min(mean_s)
  ) %>%
  pivot_longer(-1, names_to = "treatment", values_to = "mean_s")

#### final fitness ####
fit_tb <-
  rbind(fit_exp_tb, fit_stat_tb) %>%
  mutate(treatment=factor(treatment, levels=exp_stat_treatments)) %>%
  filter(!is.na(treatment)) %>%
  group_by(treatment) %>%
  left_join(LOH_wide_tb) %>%
  left_join(LOH_wide_wide_tb) %>%
  ungroup() %>%
  mutate(treatment=fct_relevel(treatment, "(average)", "(maximum s)", "(minimum s)", after=Inf)
  )

write_csv(fit_tb, "tables/fit_tb.csv")

#### fitness output table ####
fitness_readable_tb <-
  exp_tb %>%
  left_join(select(LOH_wide_tb, strain, `n_LOH-both`)) %>%
  group_by(strain, treatment) %>%
  summarise(
    mean_s=round(inlier_mean(s),3),
    sd=round(inlier_sd(s), 3)
  ) %>%
  ungroup() %>%
  filter(treatment!="YPD-30C-35cyc") %>%
  mutate(mean_text=paste0(mean_s, "±", sd)) %>%
  select(strain, treatment, mean_text) %>%
  pivot_wider(names_from = treatment, values_from = mean_text)

write_csv(fitness_readable_tb, "tables/fitness_by_strain_treatment.csv")

sd_strain_tb <-
  exp_tb %>%
  filter(treatment %in% exp_treatments) %>%
  mutate(treatment=str_replace_all(treatment, "-", "_")) %>%
  group_by(strain, treatment) %>%
  summarise(
    sd=round(inlier_sd(s), 3)
  ) %>%
  ungroup()

sd_anova <- aov(sd ~ treatment, data = sd_strain_tb)
sd_tukey <- TukeyHSD(sd_anova)
sd_cld <-
  multcompLetters(sd_tukey$treatment[,"p adj"]) %>%
  .$Letters %>%
  enframe() %>%
  rename(treatment=name,
         letters=value)

sd_labels <-
  sd_strain_tb %>%
  group_by(treatment) %>%
  summarise(
    y = max(sd)
  ) %>%
  left_join(sd_cld)

sd_strain_tb %>%
  ggplot(aes(x=treatment, y=sd)) +
  geom_boxplot() +
  geom_text(
    data = sd_labels, # Use our new data frame
    aes(x = treatment, y = 0, label = letters),
    vjust = 1,
    fontface = "bold"
  ) +
  ylab("standard deviation") +
  theme_linedraw()

ggsave("plots/sd_box.pdf", width = 8, height = 6)

fit_tb %>%
  filter(!str_detect(treatment, "worm")) %>%
  filter(treatment %in% exp_treatments) %>%
  ggplot(aes(x=mean_s, y=mean_s_withOutlier, color=treatment)) +
  geom_point(alpha=0.5) +
  sm_statCorr(color = "red", corr_method = "pearson", size=0.5, legends=TRUE) +
  scale_color_muted() +
  coord_fixed()

fit_tb %>%
  filter(!str_detect(treatment, "worm")) %>%
  filter(treatment %in% exp_treatments) %>%
  ggplot(aes(x=mean_s, y=median_s, color=treatment)) +
  geom_point(alpha=0.5) +
  sm_statCorr(color = "red", corr_method = "pearson", size=0.5, legends=TRUE) +
  scale_color_muted() +
  coord_fixed()

fit_tb %>%
  filter(!str_detect(treatment, "worm")) %>%
  filter(treatment %in% exp_treatments) %>%
  ggplot(aes(x=mean_s_withOutlier, y=median_s, color=treatment)) +
  geom_point(alpha=0.5) +
  sm_statCorr(color = "red", corr_method = "pearson", size=0.5, legends=TRUE) +
  scale_color_muted() +
  coord_fixed()

#### fitness hist ####
fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  filter(`n_LOH-both`>0) %>%
  group_by(treatment) %>%
  mutate(
    positive_rate=round(sum(mean_s>0)/n(), 3),
    mean_dis_s = round(mean(mean_s),3),
    sd=round(sd(mean_s),3)
  ) %>%
  ggplot(aes(x=mean_s)) +
  geom_histogram(breaks = seq(-0.65, 0.20, 0.05), fill="grey80") +
  geom_text(aes(x=-0.7, y=Inf, label=paste0("\n         mean = ", mean_dis_s)),  hjust = 0, size = 4, check_overlap = TRUE) +
  geom_text(aes(x=-0.7, y=Inf, label=paste0("\n\n\npositive rate = ", positive_rate)), hjust = 0, size=4, check_overlap = TRUE) +
  geom_text(aes(x=-0.7, y=Inf, label=paste0("\n\n\n\n\n               sd = ", sd)), hjust = 0, size=4, check_overlap = TRUE) +
  geom_vline(xintercept = 0, color="red", alpha=0.8, size=0.4) +
  #facet_wrap(vars(treatment), scale="free_y", ncol = 4) +
  facet_wrap(~treatment, ncol=3) +
  theme_linedraw() +
  theme(panel.grid = element_line(color="grey60")) +
  ylab("count") +
  xlab("selection coefficient (s)")

ggsave("plots/fitness_hist.pdf", width = 10, height = 8)

 #### fitness single hist ####
fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  filter(`n_LOH-both`==1) %>%
  group_by(treatment) %>%
  mutate(
    positive_rate=round(sum(mean_s>0)/n(), 3),
    mean_dis_s = round(mean(mean_s),3),
    sd=round(sd(mean_s),3)
  ) %>%
  ggplot(aes(x=mean_s)) +
  geom_histogram(breaks = seq(-0.65, 0.20, 0.05), fill="grey80") +
  geom_text(aes(x=-0.7, y=Inf, label=paste0("\n         mean = ", mean_dis_s)),  hjust = 0, size = 4, check_overlap = TRUE) +
  geom_text(aes(x=-0.7, y=Inf, label=paste0("\n\n\npositive rate = ", positive_rate)), hjust = 0, size=4, check_overlap = TRUE) +
  geom_text(aes(x=-0.7, y=Inf, label=paste0("\n\n\n\n\n               sd = ", sd)), hjust = 0, size=4, check_overlap = TRUE) +
  geom_vline(xintercept = 0, color="red", alpha=0.8, size=0.4) +
  #facet_wrap(vars(treatment), scale="free_y", ncol = 4) +
  facet_wrap(~treatment, ncol=3) +
  theme_linedraw() +
  ylab("count") +
  xlab("selection coefficient (s)")

ggsave("plots/fitness_hist_single.pdf", width = 10, height = 8)

fit_test1 <-
  fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% exp_treatments) %>%
  filter(!str_detect(treatment, "worm")) %>%
  filter(`n_LOH-both`==1) %>%
  group_by(treatment) %>%
  mutate(
    positive_rate=round(sum(mean_s>0)/n(), 3),
    mean_dis_s = round(mean(mean_s),3),
    sd=round(sd(mean_s),3)
  ) %>%
  ungroup() %>%
  select(treatment, positive_rate, mean_dis_s, sd) %>%
  unique() %>%
  mutate(subset="n_LOH==1")

fit_test2 <-
  fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% exp_treatments) %>%
  filter(!str_detect(treatment, "worm")) %>%
  filter(`n_LOH-both`>0) %>%
  group_by(treatment) %>%
  mutate(
    positive_rate=round(sum(mean_s>0)/n(), 3),
    mean_dis_s = round(mean(mean_s),3),
    sd=round(sd(mean_s),3)
  ) %>%
  ungroup() %>%
  select(treatment, positive_rate, mean_dis_s, sd) %>%
  unique() %>%
  mutate(subset="n_LOH>=1")

rbind(fit_test1, fit_test2) %>%
  dplyr::rename(`mean s`=mean_dis_s, `positive rate`=positive_rate) %>%
  pivot_longer(c(-treatment, -subset)) %>%
  pivot_wider(names_from = subset) %>%
  ggplot(aes(x=`n_LOH==1`, y=`n_LOH>=1`)) +
  sm_statCorr(r2=TRUE) +
  geom_point(aes(color=treatment)) +
  facet_wrap(~name, scales = "free") +
  theme_linedraw()

ggsave("plots/fitness_hist_compare.pdf", width = 12, height = 6)

rbind(fit_test1, fit_test2) %>%
  dplyr::rename(`mean s`=mean_dis_s, `positive rate`=positive_rate) %>%
  pivot_longer(c(-treatment, -subset)) %>%
  pivot_wider(names_from = subset) %>%
  ggplot(aes(x=`n_LOH==1`, y=`n_LOH>=1`)) +
  sm_statCorr() +
  geom_point(aes(color=treatment)) +
  facet_wrap(~name, scales = "free") +
  theme_linedraw()

ggsave("plots/fitness_hist_compare2.pdf", width = 12, height = 6)

#### report highest fitness gain for each treatment ####
fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(!treatment=="YPD-30C-after-35cyc") %>%
  filter(`n_LOH-both`>0) %>%
  group_by(treatment) %>%
  summarise(max(mean_s))

#### init freq vs fitness ####
exp_tb %>%
  select(strain, treatment, before_freq) %>%
  distinct() %>%
  right_join(fit_tb) %>%
  filter(treatment %in% c(exp_treatments)) %>%
  ggplot(aes(x=before_freq, y=mean_s)) +
  geom_point(alpha=0.2) +
  theme_linedraw() +
  facet_wrap(vars(treatment), nrow=2, scales = "free_y") +
  xlab("Initial Frequency in Bar-seq Pool") +
  ylab("selection coefficient (s)")

ggsave("plots/initFreq_vs_s.pdf", width = 12, height = 8)

#### ChrSize vs lenLOH ####
# size cap

chrLen_vs_LOHlen_1 <-
  LOH_type_tb %>%
  left_join(select(chr_len_tb, chr=X1, chr_len=X2)) %>%
  select(-type) %>%
  rename(type=`LOH_type`) %>%
  mutate(type=paste0(type, "-LOH")) %>%
  ggplot(aes(x=chr_len, length, color=type)) +
  geom_jitter(alpha=0.3, width=10000, height = 0) +
  geom_abline(slope = 1) +
  xlim(0, 1560000) +
  ylim(0, 1560000) +
  coord_equal() +
  theme_linedraw() +
  xlab("Chromosome Length") +
  ylab("LOH Tract Length") +
  theme(plot.margin = margin(20, 20, 20, 20))

chrLen_vs_LOHlen_2 <-
LOH_type_tb %>%
  left_join(select(chr_len_tb, chr=X1, chr_len=X2)) %>%
  select(-type) %>%
  rename(type=`LOH_type`) %>%
  mutate(type=paste0(type, "-LOH")) %>%
  ggplot(aes(x=chr_len, length, color=type)) +
  geom_jitter(alpha=0.3, width=10000, height = 0) +
  sm_statCorr(corr_method = "pearson", size=0.5, legends=TRUE, lty=2) +
  theme_linedraw() +
  xlab("Chromosome Length") +
  ylab("LOH Tract Length") +
  theme(plot.margin = margin(20, 20, 20, 20))

plot_grid(chrLen_vs_LOHlen_1, chrLen_vs_LOHlen_2, ncol=1, labels ="AUTO", rel_heights = c(1, 0.85), rel_widths = c(1, 0.85))
ggsave("plots/chrLen_vs_LOHlen.pdf", width = 8, height = 12)

LOH_type_tb %>%
  group_by(chr, LOH_type) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  complete(chr, LOH_type, fill=list(n=0)) %>%
  pivot_wider(names_from = "LOH_type", values_from = "n") %>%
  mutate(`pooled`=I+T) %>%
  pivot_longer(2:4) %>%
  left_join(select(chr_len_tb, chr=X1, chr_len=X2)) %>%
  mutate(name=case_match(name, "I" ~ "I-LOH", "T" ~ "T-LOH", "pooled" ~ "pooled")) %>%
  rename(`LOH type`=name, `Number of LOH`=value) %>%
  ggplot(aes(x=chr_len, y=`Number of LOH`, color=`LOH type`), alpha) +
  #geom_smooth(method="lm", fill="grey90") +
  geom_point(alpha=0.5) +
  sm_statCorr(corr_method = "pearson", size=0.5, legends=TRUE, lty=2) +
  scale_color_muted() +#values = c("pooled"="black", "I-LOH"="red", "T-LOH"="blue")) +
  theme_linedraw() +
  xlab("Chromosome Length")

ggsave("plots/ChrLen_vs_nLOH.pdf", height = 6, width = 8)

LOH_loc_tb <-
  cen_loc %>%
  mutate(chr=paste0("chr", chr)) %>%
  select(chr, ci_lo, ci_hi) %>%
  left_join(select(chr_len_tb, chr=X1, chr_len=X2)) %>%
  right_join(LOH_type_tb) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    lcdis1 = abs(ci_lo-start),
    rcdis1 = abs(ci_hi-start),
    lcdis2 = abs(ci_lo-end),
    rcdis2 = abs(ci_hi-end)
  ) %>%
  mutate(
    ltdis1 = start,
    rtdis1 = chr_len - start,
    ltdis2 = end,
    rtdis2 = chr_len - end
  ) %>%
  mutate(arm=ifelse(min(lcdis1, lcdis2) < min(rcdis1, rcdis2), "L", "R")) %>%
  mutate(min_cdis = ifelse(arm=="L", min(lcdis1, lcdis2) , min(rcdis1, rcdis2)),
         min_tdis = ifelse(arm=="L", max(ltdis1, ltdis2) , max(rtdis1, rtdis2))) %>%
  mutate(L_arm=ci_lo, R_arm=chr_len-ci_hi) %>%
  mutate(rel_cdist=ifelse(arm=="L", min_cdis/L_arm, min_cdis/R_arm),
         rel_tdist=ifelse(arm=="L", min_tdis/L_arm, min_tdis/R_arm)) %>%
  mutate(LOH_type=paste0(LOH_type, "-LOH")) %>%
  mutate(comb_type=interaction(type, LOH_type, sep=":")) %>%
  mutate(comb_type=factor(comb_type, levels=c("BY4741:I-LOH", "SK1:I-LOH", "BY4741:T-LOH", "SK1:T-LOH"))) %>%
    filter(rel_tdist<1, rel_cdist<1) %>%
  mutate(rel_length=ifelse(arm=="L", length/L_arm, length/R_arm))

type_comparison <-
  combn(c("BY4741:I-LOH", "SK1:I-LOH", "BY4741:T-LOH", "SK1:T-LOH"), m = 2) %>%
  t() %>%
  as_tibble() %>%
  mutate(v=map2(V1, V2, ~c(.x, .y))) %>%
  pull(v)

LOH_loc_plot1 <-
  LOH_loc_tb %>%
  ggplot(aes(x=comb_type, y=min_cdis)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_linedraw() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  xlab("LOH type") +
  ylab("Distance to Centromere")

LOH_loc_plot2 <-
  LOH_loc_tb %>%
  ggplot(aes(x=comb_type, y=min_tdis)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_linedraw() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  xlab("LOH type") +
  ylab("Distance to Telomere")

LOH_loc_plot3 <-
  LOH_loc_tb %>%
  ggplot(aes(x=comb_type, y=rel_cdist)) +
  geom_boxplot() +
  stat_compare_means() +
  theme_linedraw() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  xlab("LOH type") +
  ylab("Relative Distance to Centromere\n(Distance/Chromosome Arm Length)")

LOH_loc_plot4 <-
  LOH_loc_tb %>%
  ggplot(aes(x=comb_type, y=rel_tdist)) +
  geom_boxplot() +
  stat_compare_means(method="anova") +
  theme_linedraw() +
  theme(axis.text.x=element_text(angle=45, hjust = 1)) +
  xlab("LOH type") +
  ylab("Relative Distance to Telomere\n(Distance/Chromosome Arm Length)")

plot_grid(LOH_loc_plot1, LOH_loc_plot3, LOH_loc_plot2, LOH_loc_plot4, labels="AUTO")
ggsave("plots/LOH_loc_box.pdf", width = 10, height = 10)

LOH_length_dens_plot1 <-
  LOH_loc_tb %>%
  ggplot(aes(x=log10(length), color=type)) +
  geom_density() +
  theme_linedraw() +
  xlab("log10(LOH Length)") +
  ylab("Density")

LOH_length_dens_plot2 <-
  LOH_loc_tb %>%
  ggplot(aes(x=log10(rel_length), color=type)) +
  geom_density() +
  theme_linedraw() +
  xlab("log10(LOH Length/Chromosome Arm Length)") +
  ylab("Density")

LOH_length_dens_plot3 <-
  LOH_loc_tb %>%
  ggplot(aes(x=log10(length), color=LOH_type)) +
  geom_density() +
  theme_linedraw() +
  xlab("log10(LOH Length)") +
  ylab("Density")

LOH_length_dens_plot4 <-
  LOH_loc_tb %>%
  ggplot(aes(x=log10(rel_length), color=LOH_type)) +
  geom_density() +
  theme_linedraw() +
  xlab("log10(LOH Length/Chromosome Arm Length)") +
  ylab("Density")

plot_grid(LOH_length_dens_plot1, LOH_length_dens_plot2, LOH_length_dens_plot3, LOH_length_dens_plot4, ncol=2, labels = "AUTO")
ggsave("plots/LOH_length_density.pdf", width = 8, height = 6)

#### LOH len vs fitness ####
fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  filter(`LOH_length-both`>100) %>%
  select(strain, treatment, mean_s, starts_with("LOH_length")) %>%
  pivot_longer(starts_with("LOH_length"), values_to = "LOH_length", names_to = "type") %>%
  mutate(type=str_remove(type, "^LOH_length-")) %>%
  filter(type=="both") %>%
  ggplot(aes(x=log10(`LOH_length`), y=mean_s)) +
  geom_smooth(method="lm", formula = y~x) +
  geom_point(alpha=0.2) +
  stat_poly_eq(formula = y~x, parse=T, use_label("eq", "P"), small.p=TRUE, label.y = "bottom") +
  theme_linedraw() +
  theme(strip.text = element_text(size=12)) +
  theme(panel.grid = element_line(color="grey80")) +
  facet_wrap(vars(treatment), nrow=2, scales = "free_y") +
  xlab("log10(LOH length)") +
  ylab("selection coefficient (s)")

ggsave("plots/fitness_vs_LOHlen.pdf", width = 14, height = 6)

fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  filter(`LOH_length-both`>100) %>%
  select(strain, treatment, mean_s, starts_with("LOH_length")) %>%
  pivot_longer(starts_with("LOH_length"), values_to = "LOH_length", names_to = "type") %>%
  mutate(type=str_remove(type, "^LOH_length-")) %>%
  filter(type!="both") %>%
  filter(LOH_length>0) %>%
  ggplot(aes(x=log10(`LOH_length`), y=mean_s, color=type)) +
  geom_smooth(method="lm") +
  geom_point(alpha=0.2) +
  stat_poly_eq(parse=T, use_label("eq", "P"), formula=y~x, small.p=TRUE, label.y = "bottom") +
  theme_linedraw() +
  facet_wrap(vars(treatment), nrow=2, scales = "free_y") +
  scale_color_muted() +
  xlab("log10(LOH length)") +
  ylab("selection coefficient (s)")

ggsave("plots/fitness_vs_LOHlen_by_parenttype.pdf", width = 14, height = 6)

#### nLOH vs fitness ####
fitness_vs_nLOH_1 <-
  fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  select(strain, treatment, mean_s, starts_with("n_LOH")) %>%
  pivot_longer(starts_with("n_LOH"), values_to = "nLOH", names_to = "type") %>%
  mutate(type=str_remove(type, "^n_LOH-")) %>%
  filter(type=="both") %>%
  ggplot(aes(x=`nLOH`, y=mean_s)) +
  geom_smooth(method="lm", formula = y~x) +
  geom_jitter(alpha=0.2, height = 0, width = 0.2) +
  stat_poly_eq(parse=T, use_label("eq", "P"), formula=y~x, small.p=TRUE, label.y = "bottom") +
  theme_linedraw() +
  theme(plot.margin = margin(0, 1.1, 0, 0, unit = "inch")) +
  scale_x_continuous(breaks=0:10) +
  theme(panel.grid.major.x = element_blank()) +
  facet_wrap(vars(treatment), nrow=2, scales = "free_y") +
  xlab("number of LOH") +
  ylab("selection coefficient (s)")

fitness_vs_nLOH_2 <-
  fit_tb %>%
  filter(!str_detect(strain, "^n_LOH-")) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  select(-starts_with("n_LOH")) %>%
  left_join(LOH_wide_wide_tb) %>%
  select(strain, treatment, mean_s, starts_with("n_LOH")) %>%
  pivot_longer(starts_with("n_LOH"), values_to = "n_LOH", names_to = "type") %>%
  filter(str_detect(type, "-both-")) %>%
  mutate(type=str_remove(type, "^n_LOH-both-")) %>%
  mutate(type=paste(type, "LOH", sep="-")) %>%
  ggplot(aes(x=n_LOH, y=mean_s, color=type)) +
  geom_smooth(method="lm") +
  geom_jitter(alpha=0.2, height = 0, width = 0.2) +
  stat_poly_eq(parse=T, use_label("eq", "P"), formula=y~x, small.p=TRUE, label.y = "bottom") +
  theme_linedraw() +
  #theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "inch")) +
  facet_wrap(vars(treatment), nrow=2, scales = "free_y") +
  scale_color_muted() +
  xlab("number of LOH") +
  ylab("selection coefficient (s)")

plot_grid(fitness_vs_nLOH_1, fitness_vs_nLOH_2, ncol=1, labels="AUTO")
ggsave("plots/fitness_vs_nLOH.pdf", width = 14, height = 12)

#### heterozygosity vs fitness ####
fit_tb %>%
  left_join(heterozygosity_tb) %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  select(strain, treatment, mean_s, heterozygosity_genome, heterozygosity_SNPs) %>%
  ggplot(aes(x=1-`heterozygosity_genome`, y=mean_s)) +
  geom_smooth(method="lm", formula = y~x) +
  geom_point(alpha=0.2) +
  stat_poly_eq(parse=T, use_label("eq", "P"), formula=y~x, small.p=TRUE, label.y = "bottom") +
  theme_linedraw() +
  theme(panel.grid.major.x = element_blank()) +
  facet_wrap(vars(treatment), nrow=2, scales = "free_y") +
  xlab("1 - heterozygosity") +
  ylab("selection coefficient (s)")

ggsave("plots/fitness_vs_het.pdf", width = 14, height = 6)

strain_tb %>%
  filter(strain %in% good_strains) %>%
  filter(!strain %in% control_strains) %>%
  left_join(heterozygosity_tb) %>%
  left_join(LOH_wide_tb) %>%
  select(heterozygosity=heterozygosity_genome, `total length of LOH`=`LOH_length-both`) %>%
  ggplot(aes(x=1 - heterozygosity, y=`total length of LOH`)) +
  sm_statCorr() +
  geom_point(alpha=0.3) +
  theme_linedraw()

ggsave("plots/lenLOH_vs_het.pdf", width = 8, height = 8)

het_outlier_test_tb <-
strain_tb %>%
  filter(strain %in% good_strains) %>%
  filter(!strain %in% control_strains) %>%
  left_join(heterozygosity_tb) %>%
  left_join(LOH_wide_tb) %>%
  select(strain, heterozygosity=heterozygosity_genome, `total length of LOH`=`LOH_length-both`)

outlier <-
  lm(`total length of LOH`~heterozygosity, data=het_outlier_test_tb) %>% .$residuals %>% sort() %>% .[1:8] %>% names() %>% as.numeric() %>% het_outlier_test_tb[.,] %>% pull(strain)

strain_tb %>%
  filter(strain %in% good_strains) %>%
  left_join(heterozygosity_tb) %>%
  right_join(fit_tb) %>%
  select(strain, heterozygosity=heterozygosity_genome, `total length of LOH`=`LOH_length-both`, mean_s, treatment) %>%
  arrange(strain %in% outlier) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  ggplot(aes(x=1 - heterozygosity, y=`total length of LOH`)) +
  sm_statCorr() +
  geom_point(aes(color=strain %in% outlier)) +
  theme_linedraw()

strain_tb %>%
  filter(strain %in% good_strains) %>%
  filter(!strain %in% control_strains) %>%
  left_join(heterozygosity_tb) %>%
  right_join(fit_tb) %>%
  select(strain, heterozygosity=heterozygosity_genome, `total length of LOH`=`LOH_length-both`, mean_s, treatment) %>%
  filter(!strain %in% control_strains) %>%
  arrange(strain %in% outlier) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  ggplot(aes(x=1 - heterozygosity, y=mean_s)) +
  sm_statCorr() +
  geom_point(aes(color=strain %in% outlier)) +
  facet_wrap(~treatment) +
  theme_linedraw()

strain_tb %>%
  filter(strain %in% good_strains) %>%
  filter(!strain %in% control_strains) %>%
  left_join(heterozygosity_tb) %>%
  right_join(fit_tb) %>%
  select(strain, heterozygosity=heterozygosity_genome, `total length of LOH`=`LOH_length-both`, mean_s, treatment) %>%
  filter(!strain %in% control_strains) %>%
  arrange(strain %in% outlier) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  ggplot(aes(x=log10(`total length of LOH`), y=mean_s)) +
  sm_statCorr() +
  geom_point(aes(color=strain %in% outlier)) +
  facet_wrap(~treatment) +
  theme_linedraw()

#### 0LOH vs control vs LOH ####
my_comparisons <- list( c("control", "0 LOH"), c("control", ">0 LOH"), c("0 LOH", ">0 LOH"))

fit_tb %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  mutate(`LOH state`=
           case_when(strain %in% control_strains ~ "control",
                     `n_LOH-both` == 0 ~ "0 LOH",
                     `n_LOH-both` > 0 ~ ">0 LOH",
           )
         ) %>%
  ggboxplot(x = "LOH state", y = "mean_s") +
  stat_compare_means(method="wilcox.test", comparisons = my_comparisons, label = "p.signif") # Add brackets


control_vs_0LOH <-
fit_tb %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  mutate(`LOH_state`=
           case_when(strain %in% control_strains ~ "control",
                     `n_LOH-both` == 0 ~ "0 LOH",
                     `n_LOH-both` > 0 ~ ">0 LOH"
           )) %>%
  mutate(LOH_state=factor(LOH_state, levels=c("control", "0 LOH", ">0 LOH")))

control_vs_0LOH_test <-
  control_vs_0LOH %>%
  group_by(treatment) %>%
  tukey_hsd(mean_s ~ LOH_state, p.adjust.method="BH") %>%
  add_significance("p.adj") %>%
  add_y_position()

control_vs_0LOH %>%
    ggplot(aes(x=LOH_state, y=mean_s)) +
  geom_boxplot(outliers = FALSE) +
    facet_wrap(vars(treatment), nrow=2, scales = "free_y") +
    stat_pvalue_manual(control_vs_0LOH_test, label = "p.adj.signif",
                       tip.length = 0.01) +
  geom_jitter(width = 0.1, alpha=0.3, height = 0) +
  theme_linedraw() +
  scale_color_muted()

ggsave("plots/control_vs_0LOH.pdf", height = 10, width = 10)

#### top strain fitness across environemts ####
rank_tb <-
  fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% exp_treatments) %>%
  filter(!str_detect(treatment, "worm")) %>%
  filter(`n_LOH-both`>0) %>%
  arrange(desc(mean_s)) %>%
  group_by(treatment) %>%
  slice_head(n=5) %>%
  group_by(strain) %>%
  summarise(top5_in = paste(treatment, collapse = ", "),
            top5_in_vector=list(c(treatment))
  )

rank_tb %>%
  select(strain, top5_in) %>%
  write_tsv("tables/top5_list.tsv")

exp_tb %>%
  filter(!treatment=="YPD-30C-35cyc") %>%
  filter(!str_detect(treatment, "worm")) %>%
  right_join(select(rank_tb, strain, top5_in)) %>%
  group_by(strain, treatment) %>%
  mutate(mean_s=inlier_mean(s)) %>%
  ungroup() %>%
  arrange(top5_in, mean_s) %>%
  mutate(strain=factor(strain, levels=unique(strain))) %>%
  ggplot() +
  geom_boxplot(aes(x=strain, y=s, color=top5_in)) +
  geom_abline(slope = 0, intercept = 0, lty=2, alpha=0.5) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  facet_wrap(~treatment, scale="free_y", ncol = 2, dir="v") +
  scale_color_muted(name="top 5 in") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("selection coefficient (s)")

ggsave("plots/barseq_tops_enrichment.pdf", height = 8, width = 10)


#### pleiotropy analysis ####
pleiotropy_tb <-
  fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% exp_treatments) %>%
  filter(!str_detect(treatment, "worm")) %>%
  filter(!str_detect(treatment, "H2O2")) %>%
  filter(`n_LOH-both`>0) %>%
  group_by(strain) %>%
  summarise(
    max_s=max(mean_s),
    min_s=min(mean_s),
    mean_mean_s=mean(mean_s),
    gmean_mean_s=exp(mean(log(mean_s+1)))-1,
    sd=sd(mean_s),
    mean_cost=mean(pmin(sort(mean_s, decreasing=TRUE)[-1], 0)),
    max_cost=min(pmin(sort(mean_s, decreasing=TRUE)[-1], 0)),
    tradeoff_meancost=-mean_cost*(max(mean_s)),
    tradeoff_maxcost=-max_cost*(max(mean_s))
  ) %>%
  mutate(
    range=max_s-min_s,
  ) %>%
  ungroup() %>%
  mutate(
    direction=case_when(
      max_s > 0 & min_s < 0 ~ "antagonistic pleiotropy",
      max_s > 0 & min_s > 0 ~ "universal benefit",
      max_s < 0 & min_s < 0 ~ "universal cost")
  ) %>%
  left_join(LOH_wide_tb)

pleiotropy_tb %>%
  group_by(direction) %>%
  tally()

pleiotropy_tb %>%
  mutate(max_gain=max_s, max_cost=-min_s, mean_cost=-mean_cost) %>%
  select(strain, max_gain, max_cost, mean_cost, tradeoff_meancost, sd, range, `LOH_length-both`) %>%
  rename("Maximum Gain"=max_gain,
         "Maximum Cost"=max_cost,
         "Average Cost"=mean_cost,
         "Trade-Off (Max-Gain*Average-Cost)"=tradeoff_meancost,
         "Standard Deviation"=sd,
         "Range (Max-Gain - Max-Cost)"=range
  ) %>%
  pivot_longer(c(-strain, -`LOH_length-both`)) %>%
  mutate(name=factor(name, levels = c(
    "Maximum Gain",
    "Maximum Cost",
    "Average Cost",
    "Trade-Off (Max-Gain*Average-Cost)",
    "Standard Deviation",
    "Range (Max-Gain - Max-Cost)"
  ))
  ) %>%
  filter(`LOH_length-both`>100) %>%
  ggplot(aes(x=log10(`LOH_length-both`), y=value)) +
  geom_point(alpha=0.3) +
  geom_smooth(method=lm, formula=y ~ exp(x)) +
  stat_fit_tidy(method = "lm", method.args = list(formula = y ~ exp(x)),
                aes(label = sprintf("y=%.3g^x+%.3g,\np=%.3g",
                                    after_stat(`exp(x)_estimate`),
                                    after_stat(`Intercept_estimate`),
                                    after_stat(`exp(x)_p.value`)
                )
                )) +
  facet_wrap(~name, scales = "free_y") +
  theme_linedraw() +
  theme(strip.text = element_text(size=12)) +
  theme(panel.grid = element_line(color="grey80")) +
  xlab("log10(LOH length)") +
  ylim(0, NA)

ggsave("plots/tradeoff.pdf", width = 11, height = 6)

pleiotropy_tb %>%
  mutate(direction=factor(direction, levels=c("antagonistic pleiotropy", "universal cost", "universal benefit"))) %>%
  ggplot(aes(x=min_s, y=max_s)) +
  geom_point(aes(color=direction)) +
  coord_fixed() +
  geom_hline(yintercept = 0, lty=1, size=0.4) +
  geom_vline(xintercept = 0, lty=1, size=0.4) +
  geom_abline(slope=1, intercept = 0, lty=2, size=0.4) +
  geom_polygon(data=tibble(x=c(-Inf, Inf, Inf, -Inf), y=c(-Inf, Inf, -Inf, -Inf)), aes(x=x, y=y), fill="grey20", alpha=0.2) +
  xlim(-0.4, 0.4) +
  ylim(-0.4, 0.4) +
  xlab("minimum s") +
  ylab("maximum s") +
  scale_color_muted(drop=FALSE)

ggsave("plots/tradeoff_quadrant.pdf", width = 12, height = 8)

#### QTL ####
cM_unit <- 30/1e6

pheno_tb <-
  fit_tb %>%
  select(strain, treatment, mean_s) %>%
  pivot_wider(names_from = treatment, values_from = mean_s) %>%
  column_to_rownames("strain")

geno <-
  read.table("bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.named.vcf.raw", header = TRUE)

map_tb <-
  read.table("bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.named.vcf.map", header = FALSE) %>%
  select(chr=V1, marker=V2, cM=V3, pos=V4) %>%
  mutate(chr=str_replace(chr, "23", "chrX")) %>%
  mutate(cM=pos*cM_unit)

geno_dt <-
  geno %>%
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE) %>%
  column_to_rownames("IID")

colnames(geno_dt) <- str_remove(colnames(geno_dt), "_[^_]+$")
geno_dt[geno_dt==0] <- "A"
geno_dt[geno_dt==1] <- "H"
geno_dt[geno_dt==2] <- "B"
geno_dt[is.na(geno_dt)] <- "-"

geno_dt <-
  geno_dt[,map_tb$marker]

head1 <-
  c(names(pheno_tb), colnames(geno_dt))

head2 <-
  c(rep("", length(names(pheno_tb))), map_tb$chr)

head3 <-
  c(rep("", length(names(pheno_tb))), map_tb$cM)

dt <-
  left_join(rownames_to_column(pheno_tb, "ID"), rownames_to_column(geno_dt, "ID")) %>%
  column_to_rownames("ID")

cross_file <-
  tempfile()

write_lines(c(paste0(head1, collapse = ","), paste0(head2, collapse = ","), paste0(head3, collapse = ",")), cross_file)
write_csv(dt, cross_file, append = TRUE)

cross <-
  read.cross(format="csv", file=cross_file, convertXdata=FALSE)
cross2 <- convert2cross2(cross)
map <- insert_pseudomarkers(cross2$gmap, step=0.1)
pr <- calc_genoprob(cross2, map, cores=4)
out <- scan1(pr, cross2$pheno)

# QTL plots
par(mar=c(5.1, 4.1, 1.1, 1.1))
pdf(file = "plots/QTL_all.pdf", width = 11, height = 6)
ymx <- maxlod(out) # overall maximum LOD score
plot(out, map, lodcolumn=1, col=alpha(1, 0.3), ylim=c(0, ymx*1.02))
plot(out, map, lodcolumn=2, col=alpha(2, 0.3), add=TRUE)
plot(out, map, lodcolumn=3, col=alpha(3, 0.3), add=TRUE)
plot(out, map, lodcolumn=4, col=alpha(4, 0.3), add=TRUE)
plot(out, map, lodcolumn=5, col=alpha(5, 0.3), add=TRUE)
plot(out, map, lodcolumn=7, col=alpha(6, 0.3), add=TRUE)
plot(out, map, lodcolumn=8, col=alpha(7, 0.3), add=TRUE)
legend("topleft", lwd=2, col=1:7, colnames(out)[c(1, 2, 3, 4, 5, 7, 8)], pch=15, cex=1)
dev.off()

find_peaks(out, map, threshold=4) %>%
  mutate(pos_bp=pos/cM_unit)

par(mar=c(5.1, 4.1, 1.1, 1.1))
pdf(file = "plots/QTL_1.pdf", width = 8, height = 6)
plot(out, map, lodcolumn=1, col=alpha(1, 0.3))
legend("topleft", lwd=2, col=1, colnames(out)[c(1)], pch=15, cex=1)
dev.off()

par(mar=c(5.1, 4.1, 1.1, 1.1))
pdf(file = "plots/QTL_2.pdf", width = 8, height = 6)
plot(out, map, lodcolumn=2, col=alpha(1, 0.3))
legend("topleft", lwd=2, col=1, colnames(out)[c(2)], pch=15, cex=1)
dev.off()

par(mar=c(5.1, 4.1, 1.1, 1.1))
pdf(file = "plots/QTL_3.pdf", width = 8, height = 6)
plot(out, map, lodcolumn=3, col=alpha(1, 0.3))
legend("topleft", lwd=2, col=1, colnames(out)[c(3)], pch=15, cex=1)
dev.off()

par(mar=c(5.1, 4.1, 1.1, 1.1))
pdf(file = "plots/QTL_4.pdf", width = 8, height = 6)
plot(out, map, lodcolumn=4, col=alpha(1, 0.3))
legend("topleft", lwd=2, col=1, colnames(out)[c(4)], pch=15, cex=1)
dev.off()

par(mar=c(5.1, 4.1, 1.1, 1.1))
pdf(file = "plots/QTL_5.pdf", width = 8, height = 6)
plot(out, map, lodcolumn=5, col=alpha(1, 0.3))
legend("topleft", lwd=2, col=1, colnames(out)[c(5)], pch=15, cex=1)
dev.off()

par(mar=c(5.1, 4.1, 1.1, 1.1))
pdf(file = "plots/QTL_7.pdf", width = 8, height = 6)
plot(out, map, lodcolumn=7, col=alpha(1, 0.3))
legend("topleft", lwd=2, col=1, colnames(out)[c(7)], pch=15, cex=1)
dev.off()

par(mar=c(5.1, 4.1, 1.1, 1.1))
pdf(file = "plots/QTL_8.pdf", width = 8, height = 6)
plot(out, map, lodcolumn=8, col=alpha(1, 0.3))
legend("topleft", lwd=2, col=1, colnames(out)[c(8)], pch=15, cex=1)
dev.off()

#### LOH stat and map on chromosome ####
LOH_stat_plot1 <-
  LOH_long_tb %>%
  filter(!strain %in% control_strains) %>%
  mutate(strain=fct_drop(strain)) %>%
  complete(strain, type, fill=list(LOH_length=0)) %>%
  group_by(strain, type) %>%
  summarise(n_LOH=sum(n_LOH), LOH_length=sum(LOH_length)) %>%
  rename(`Cumulative LOH lenghth`=LOH_length) %>%
  mutate(
    bin = cut(`Cumulative LOH lenghth`, breaks = c(0, 1e-10, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), include.lowest = TRUE)
  ) %>%
  ungroup() %>%
  count(bin, type) %>%
  complete(bin, type, fill=list(n=0)) %>%
  mutate(bin=case_match(bin, "[0,1e-10]" ~ "(0,0]", "(1e-10,10]" ~ "(0,10]", .default=bin)) %>%
  mutate(type=factor(type, levels=c("BY4741", "both", "SK1"))) %>%
  ggplot(aes(x=bin, y=n, fill=type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  scale_fill_muted(name="LOH type") +
  theme_linedraw() +
  theme(legend.position = "none", plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "inch"),
        axis.text.x = element_text(angle=30, hjust = 1)) +
  xlab("Cumulative LOH Length in Strains") +
  ylab("count")

LOH_stat_plot2 <-
  rbind(LOH_type_tb, mutate(LOH_type_tb, type="both")) %>%
  filter(!strain %in% control_strains) %>%
  mutate(
    bin = cut(length, breaks = c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7))
  ) %>%
  count(bin, type) %>%
  complete(bin, type, fill=list(n=0)) %>%
  group_by(type) %>%
  mutate(sum=sum(n)) %>%
  mutate(frequency=n/sum) %>%
  mutate(type=factor(type, levels=c("BY4741", "both", "SK1"))) %>%
  ggplot(aes(x = bin, y = frequency, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "inch",),
        legend.position = "none"
  ) +
  xlab("Individual LOH Length") +
  scale_fill_muted(name="LOH type")

LOH_stat_plot3 <-
  LOH_long_tb %>%
  filter(!strain %in% control_strains) %>%
  ungroup() %>%
  mutate(strain=fct_drop(strain)) %>%
  group_by(type, strain) %>%
  summarise(n_LOH=sum(n_LOH), LOH_length=sum(LOH_length)) %>%
  complete(fill=list(n_LOH=0, LOH_length=0)) %>%
  group_by(n_LOH, type) %>%
  tally() %>%
  ungroup() %>%
  mutate(type=ifelse(type=="both", "pooled", type)) %>%
  mutate(type=factor(type, levels=c("BY4741", "pooled", "SK1"))) %>%
  ggplot(aes(x=n_LOH, y=n, fill=type)) +
  geom_col(position="dodge2", width=0.5) +
  scale_fill_muted(name="LOH type") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  theme_linedraw() +
  theme(legend.position = "none", plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "inch")) +
  xlab("Number of LOH in Strains") +
  ylab("count")

LOH_stat_legend1 <-
  get_legend(LOH_stat_plot3+theme_linedraw())

LOH_stat_row1 <-
  plot_grid(LOH_stat_plot1, LOH_stat_plot2, LOH_stat_plot3, LOH_stat_legend1, nrow=1, rel_widths = c(1, 1, 1, 0.3), labels = c("A", "B", "C", ""))

LOH_stat_plot4 <-
  rbind(bedgraph_2, bedgraph_3) %>%
  select(chr=seqnames, ci_lo=start, ci_hi=end, pos=start, lodcolumn=score, type) %>%
  mutate(lodindex="depth", col="transparent", chr=str_remove(chr, "chr")) %>%
  ggplot_peaks(map, gap=25, bgcolor="gray90", altbgcolor="gray85", tick_height = 0.8) +
  geom_rect(aes(xmin = ci_lo, xmax=ci_hi, ymin = 0, ymax=lodcolumn, fill=type), alpha=0.3) +
    geom_rect(data=mutate(LOH_type_tb, chr=str_remove(chr, "chr")), aes(xmin = (start+end)/2-5000, xmax = (start+end)/2+5000, ymax = 0, ymin=-0.5), alpha=0.3, inherit.aes = FALSE) +
  geom_rect(data=cen_loc, aes(xmin = ci_lo-15000, xmax=ci_hi+15000, ymax = -0.5, ymin=-1.0), alpha=1, fill="red", inherit.aes = FALSE) +
  geom_rect(data=SNP_loc, aes(xmin = ci_lo, xmax=ci_hi, ymax = -1.0, ymin=-1.5), alpha=0.5, inherit.aes = FALSE, fill="blue") +
  scale_y_discrete() +
  theme_bw(base_size = 11, base_family = "",
           base_line_size = 11/22, base_rect_size = 11/22) %+replace%
    theme(axis.text = element_text(colour = "black", size = rel(0.8)),
          axis.ticks = element_line(colour = "black", linewidth = rel(0.5)),
          panel.border = element_rect(fill = NA, colour = "black",
                                      linewidth = rel(1)), panel.grid = element_line(colour = "black"),
          panel.grid.major = element_line(linewidth = rel(0.1)),
          panel.grid.minor = element_line(linewidth = rel(0.05)),
          strip.background = element_rect(color = "black"),
          strip.text = element_text(colour = "black", size = rel(0.8),
                                    margin = margin(0.8 * 11/2/2, 0.8 * 11/2/2,
                                                    0.8 * 11/2/2, 0.8 * 11/2/2)), complete = TRUE) +
  theme(legend.position = "none", axis.ticks.x=element_blank(), axis.text.x=element_blank(), panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
  ylab("Depth (number of strains)")

LOH_stat_legend2 <-
  get_legend(
    ggplot(tibble(`LOH type`=rep(c("BY4741", "SK1"), 3), line=rep(factor(c("LOH event", "centromere", "SNP"), levels=c("LOH event", "centromere", "SNP")), 2))) +
      geom_line(aes(x=1, y=1, color=line)) +
      geom_raster(aes(fill=`LOH type`, x=1, y=1), alpha=0.3) +
      scale_color_manual(name="", values = c("centromere"="red", "SNP"="blue", "LOH event"="black")) +
      guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 2)
      ) +
      theme_linedraw()
  )

LOH_stat_row2 <-
  plot_grid(LOH_stat_plot4, LOH_stat_legend2, nrow=1, rel_widths = c(1, 0.1), labels=c("D",""))

LOH_stat_plot <-
  plot_grid(LOH_stat_row1, LOH_stat_row2, ncol = 1, rel_heights = c(1, 0.8))

ggsave("plots/LOH_stat.pdf", width = 11, height = 8)

#### LOH stats ####

SK1_cov <-
  bedgraph_2 %>%
  mutate(len=end-start) %>%
  pull(len) %>%
  sum() %>%
  `/`(genome_length)

BY4741_cov <-
  bedgraph_3 %>%
  mutate(len=end-start) %>%
  pull(len) %>%
  sum() %>%
  `/`(genome_length)

any_cov <-
  bedgraph_1 %>%
  mutate(len=end-start) %>%
  pull(len) %>%
  sum() %>%
  `/`(genome_length)

ind_LOH_stats1 <-
  LOH_type_tb %>%
  filter(!strain %in% control_strains) %>%
  group_by(type) %>%
  summarise(mean=mean(length),
            median=median(length),
            sd=sd(length),
            min=min(length),
            max=max(length))

ind_LOH_stats2 <-
  LOH_type_tb %>%
  filter(!strain %in% control_strains) %>%
  summarise(mean=mean(length),
            median=median(length),
            sd=sd(length),
            min=min(length),
            max=max(length)) %>%
  mutate(type="both")

ind_LOH_stats <-
  rbind(ind_LOH_stats1, ind_LOH_stats2)

ind_LOH_stats %>%
  knitr::kable()

LOH_long_tb %>%
  filter(!strain %in% control_strains) %>%
  group_by(type) %>%
  filter(!strain %in% control_strains) %>%
  summarise(across(c(-strain, -LOH_type), .fns = list(mean=mean, median=median, sd=sd, min=min, max=max))) %>%
  knitr::kable()

LOH_long_tb %>%
  filter(!strain %in% control_strains) %>%
  group_by(strain) %>%
  filter(sum(n_LOH)!=0) %>%
  group_by(type) %>%
  filter(!strain %in% control_strains) %>%
  summarise(across(c(-strain, -LOH_type), .fns = list(mean=mean, median=median, sd=sd, min=min, max=max))) %>%
  knitr::kable()

LOH_type_tb %>%
  filter(!strain %in% control_strains) %>%
  group_by(type, LOH_type) %>%
  summarise(n=n(), length_median=median(length), length_75=quantile(length, 0.75), length_25=quantile(length, 0.25)) %>%
  ungroup() %>%
  mutate(n_p=n/sum(n))

LOH_type_tb %>%
  filter(!strain %in% control_strains) %>%
  group_by(type) %>%
  summarise(n=n(), length_median=median(length), length_75=quantile(length, 0.75), length_25=quantile(length, 0.25)) %>%
  ungroup() %>%
  mutate(n_p=n/sum(n))

LOH_type_tb %>%
  filter(!strain %in% control_strains) %>%
  group_by(LOH_type) %>%
  summarise(n=n(), length_median=median(length), length_75=quantile(length, 0.75), length_25=quantile(length, 0.25)) %>%
  ungroup() %>%
  mutate(n_p=n/sum(n))

type_plot_1 <-
  LOH_sep_tb %>%
  filter(!strain %in% control_strains) %>%
  group_by(strain, type) %>%
  summarise(n=sum(n_LOH)) %>%
  ggplot(aes(x=type, y=n)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(alpha=0.3, width = 0.05, height = 0.2) +
  stat_compare_means() +
  theme_linedraw() +
  ylab("Number of LOH in a Strain") +
  xlab("LOH type")


type_plot_2 <-
  LOH_sep_tb %>%
  filter(!strain %in% control_strains) %>%
  group_by(strain, LOH_type) %>%
  summarise(n=sum(n_LOH)) %>%
  mutate(LOH_type=paste0(LOH_type, "-LOH")) %>%
  ggplot(aes(x=LOH_type, y=n)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(alpha=0.3, width = 0.05, height = 0.2) +
  stat_compare_means() +
  theme_linedraw() +
  ylab("Number of LOH in a Strain") +
  xlab("LOH type")

plot_grid(type_plot_1, type_plot_2, ncol=1, labels = "AUTO")
ggsave("plots/nLOH_strain_box.pdf", height = 10, width = 10)

LOH_wide_tb %>%
  filter(!strain %in% control_strains) %>%
  mutate(n=n()) %>%
  filter(`n_LOH-BY4741`==0) %>%
  reframe(n()/unique(n))

LOH_wide_tb %>%
  filter(!strain %in% control_strains) %>%
  mutate(n=n()) %>%
  filter(`LOH_length-SK1`>1000) %>%
  reframe(n()/unique(n))

LOH_wide_tb %>%
  filter(!strain %in% control_strains) %>%
  mutate(n=n()) %>%
  filter(`LOH_length-SK1`>1000) %>%
  reframe(n()/unique(n))

LOH_type_tb %>%
  mutate(bin=cut(length,breaks = c(0, 1000, 1000000, Inf))) %>%
  group_by(bin) %>%
  tally() %>%
  ungroup() %>%
  mutate(n_p=n/sum(n))

#### LOH correlation ####
LOH_wide_tb %>%
  filter(!strain %in% control_strains) %>%
  group_by(`n_LOH-SK1`, `n_LOH-BY4741`) %>%
  summarise(n=n()) %>%
  ggplot(aes(x= `n_LOH-SK1`, y=`n_LOH-BY4741`, fill=n)) +
  geom_raster() +
  coord_fixed() +
  theme_linedraw()

my_comparisons <- list( c("0 LOH", "1 LOH"), c("0 LOH", ">1 LOH"), c("1 LOH", ">1 LOH"))

LOH_wide_tb %>%
  filter(!strain %in% control_strains) %>%
  mutate(SK1_LOH=case_when(`n_LOH-SK1`==0 ~ "0 LOH", `n_LOH-SK1`==1 ~ "1 LOH", `n_LOH-SK1`>1 ~ ">1 LOH")) %>%
  mutate(SK1_LOH=factor(SK1_LOH, levels=c("0 LOH", "1 LOH",　">1 LOH"))) %>%
  ggplot(aes(x=SK1_LOH, y=`n_LOH-BY4741`)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(height = 0.2, width = 0.2, alpha=0.3) +
  stat_compare_means(method="wilcox.test", comparisons = my_comparisons, label = "p.signif") +
  theme_linedraw() +
  xlab("number of SK1-type LOH") +
  ylab("number of BY4741-type LOH")

ggsave("plots/SK1-type_LOH_vs_BY4741-type_LOH.pdf", width = 5, height = 7)


type_length_category <-
  LOH_type_tb %>%
  left_join(select(chr_len_tb, chr=X1, chr_len=X2)) %>%
  mutate(comb_type=interaction(type, LOH_type, sep=":"))

type_length_anova <- aov(length ~ comb_type, data = type_length_category)
type_length_tukey <- TukeyHSD(type_length_anova)
type_length_cld <-
  multcompLetters(type_length_tukey$comb_type[,"p adj"]) %>%
  .$Letters %>%
  enframe() %>%
  rename(comb_type=name,
         letters=value)

type_length_labels <-
  type_length_category %>%
  group_by(comb_type) %>%
  summarise(
    y = max(length)
  ) %>%
  left_join(type_length_cld) %>%
  mutate(comb_type=paste0(comb_type, "-LOH")) %>%
  mutate(comb_type=factor(comb_type, levels=c("BY4741:I-LOH", "SK1:I-LOH", "BY4741:T-LOH", "SK1:T-LOH")))

type_length_category %>%
  mutate(comb_type=paste0(comb_type, "-LOH")) %>%
  mutate(comb_type=factor(comb_type, levels=c("BY4741:I-LOH", "SK1:I-LOH", "BY4741:T-LOH", "SK1:T-LOH"))) %>%
  ggplot(aes(x=comb_type, y=log10(length))) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(alpha=0.3, width = 0.1, height = 0) +
  geom_text(
    data = type_length_labels, # Use our new data frame
    aes(x = comb_type, y = 0, label = letters),
    vjust = 1,
    fontface = "bold"
  ) +
  xlab("LOH Type") +
  ylab("LOH Tract Length (log10)") +
  theme_linedraw()

ggsave("plots/LOH_tract_type_vs_length_box.pdf", height = 8, width = 6)

# LOH distance
LOH_dist_tb <-
  crossing(filter(LOH_type_tb, type=="BY4741"), filter(LOH_type_tb, type=="SK1"), .name_repair="universal") %>%
  filter(strain...4==strain...11, chr...1 == chr...8, start...2 !=start...9) %>%
  rowwise() %>%
  mutate(dis1=abs(start...2-end...10), dis2=abs(end...3-start...9),
         min_dis=min(dis1, dis2)) %>%
  group_by(chr...1, start...2, end...3) %>%
  summarise(min_min_dis=min(min_dis)) %>%
  ungroup() %>%
  mutate(bin=cut(min_min_dis, c(0, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000), include.lowest = TRUE)) %>%
  count(bin)

sum(pull(LOH_dist_tb[1:4,],n))/sum(LOH_type_tb$type=="BY4741")

#### barseq stats ####
barseq_stats <-
  good_bc_tb %>%
  pivot_longer(-strain) %>%
  group_by(name) %>%
  summarise(depth=sum(value, na.rm = TRUE)) %>%
  mutate(name=str_remove(name, "-[^\\-]+$")) %>%
  group_by(name) %>%
  summarise(`number of replicate`=n(), `mean read`=mean(depth), `median read`=median(depth), `read sd`=sd(depth), `min read`=min(depth), `max read`=max(depth)) %>%
  mutate(name=str_remove(name, "-after")) %>%
  filter(!name %in% c("YPD-30C-35cyc", "worm-37C")) %>%
  mutate(name=case_match(name, "LOH-pool" ~ "inoculum (others)", "Worm-pool" ~ "inoculum (worm)", .default = name)) %>%
  arrange(str_detect(name, "^inoculum"))

write_csv(barseq_stats, "tables/barseq_stats.csv")

tibble(treatment=c("YPD", "BL", "H2O2", "HT", "CR", "NaCl", "EtOH")) %>%
  mutate(path=paste0("Vijayan_MA_lines/pgen.1011692.s014_", treatment, ".csv")) %>%
  mutate(df=map(path, ~read_csv(.x, skip = 1) %>% filter(Number_of_SNPs_supporting_LOH>=5) %>%  group_by(LOH_type) %>% tally())) %>%
  select(-path) %>%
  unnest(df) %>%
    pivot_wider(names_from = "LOH_type", values_from = "n") %>%
    mutate(sum=Interstitial+Terminal) %>%
    mutate(Interstitial_percent=Interstitial/sum,
           Terminal_percent=Terminal/sum)

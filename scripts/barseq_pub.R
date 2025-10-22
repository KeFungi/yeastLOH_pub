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

#### functions ####
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

fordLOH_tb <-
  read_tsv("LOH_detect/LOH_minSNP-5.bed", col_names = FALSE) %>%
  mutate(type="SK1")

revLOH_tb <-
  read_tsv("LOH_detect/revLOH_minSNP-5.bed", col_names = FALSE) %>%
  mutate(type="BY4741")

LOH_type_tb <-
  rbind(fordLOH_tb, revLOH_tb) %>%
  rename(chr=X1, start=X2, end=X3, strain=X4) %>%
  mutate(length=end-start) %>%
  filter(strain %in% good_strains)

LOH_sep_tb <-
  LOH_type_tb %>%
  mutate(strain=factor(strain, levels=good_strains)) %>%
  group_by(strain, type) %>%
  summarise(n_LOH=n(), LOH_length=sum(length)) %>%
  ungroup() %>%
  complete(strain, type, fill = list(n_LOH=0, LOH_length=0))

LOH_long_tb <-
  LOH_sep_tb %>%
  group_by(strain) %>%
  summarise(n_LOH=sum(n_LOH), LOH_length=sum(LOH_length)) %>%
  mutate(type="both") %>%
  rbind(LOH_sep_tb, .) %>%
  arrange(strain, type)

LOH_wide_tb <-
  LOH_long_tb %>%
  pivot_wider(id_cols = "strain", names_from = "type", values_from = c("n_LOH", "LOH_length"), names_sep = "-")

LOH_long_tb %>%
  filter(n_LOH>0) %>%
  ggplot(aes(x=n_LOH, y=log10(LOH_length))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ log(x)) +
  stat_poly_eq(method = "lm", formula = y ~ log(x), small.p = TRUE, use_label(c("eq", "P")), label.y = "bottom")

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
  c("Acidic", "Caffeine", "H2O2", "SC",
    "YPD-30C", "YPD-37C", "worm-20C")

exp_stat_treatments <-
  c("Acidic", "Caffeine", "H2O2", "SC",
    "YPD-30C", "YPD-37C", "worm-20C", "(average)", "(maximum s)", "(minimum s)")

LOH_dist_tb <-
  crossing(filter(LOH_type_tb, type=="BY4741"), filter(LOH_type_tb, type=="SK1"), .name_repair="universal") %>%
  filter(strain...4==strain...10, chr...1 == chr...7, start...2 !=start...8) %>%
  rowwise() %>%
  mutate(dis1=abs(start...2-end...9), dis2=abs(end...3-start...8),
         min_dis=min(dis1, dis2)) %>%
  group_by(chr...1, start...2, end...3) %>%
  summarise(min_min_dis=min(min_dis)) %>%
  ungroup() %>%
  mutate(bin=cut(min_min_dis, c(0, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000), include.lowest = TRUE)) %>%
  count(bin)

sum(pull(LOH_dist_tb[1:4,],n))/sum(LOH_type_tb$type=="BY4741")

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

media_pool_tb %>%
  ggplot(aes(x=mean_freq_withOutlier, y=mean_freq)) +
  geom_point(alpha=0.5) +
  sm_statCorr(color = "red", corr_method = "pearson", size=0.5, legends=TRUE) +
  scale_color_muted() +
  coord_fixed()

ggsave("plots/freq_mean_compare_1.png", width = 12, height = 7)

media_pool_tb %>%
  ggplot(aes(x=mean_freq_withOutlier, y=median_freq)) +
  geom_point(alpha=0.5) +
  sm_statCorr(color = "red", corr_method = "pearson", size=0.5, legends=TRUE) +
  scale_color_muted() +
  coord_fixed()

ggsave("plots/freq_mean_compare_2.png", width = 12, height = 7)

media_pool_tb %>%
  ggplot(aes(x=mean_freq, y=median_freq)) +
  geom_point(alpha=0.5) +
  sm_statCorr(color = "red", corr_method = "pearson", size=0.5, legends=TRUE) +
  scale_color_muted() +
  coord_fixed()

ggsave("plots/freq_mean_compare_3.png", width = 12, height = 7)

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
  mutate(treatment=fct_relevel(treatment, "worm-20C", after=Inf)) %>%
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
  filter(!str_detect(treatment, "worm")) %>%
  filter(!str_detect(treatment, "H2O2")) %>%
  group_by(strain) %>%
  summarise("(average)"=mean(mean_s),
            "(maximum s)"=max(mean_s),
            "(minimum s)"=min(mean_s)
  ) %>%
  pivot_longer(-1, names_to = "treatment", values_to = "mean_s")

fit_tb <-
  rbind(fit_exp_tb, fit_stat_tb) %>%
  group_by(treatment) %>%
  left_join(LOH_wide_tb) %>%
  ungroup() %>%
  mutate(treatment=fct_relevel(treatment, "(average)", "(maximum s)", "(minimum s)", after=Inf)
  )

write_csv(fit_tb, "tables/fit_tb.csv")

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

exp_tb %>%
  group_by(strain, treatment) %>%
  summarise(
    mean_s=round(inlier_mean(s),3),
    sd=round(inlier_sd(s), 3)
  ) %>%
  ungroup() %>%
  ggplot(aes(x=treatment, y=sd)) +
  geom_boxplot() +
  ylab("standard deviation")

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

fit_tb %>%
  pivot_wider(id_cols = "strain", names_from = "treatment", values_from = "mean_s") %>%
  ggplot(aes(x=`YPD-30C`, y=`YPD-30C-35cyc`)) +
  geom_point() +
  sm_statCorr(color = "red", corr_method = "pearson", size=0.5) +
  scale_color_muted() +
  coord_fixed()

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
  ylab("count") +
  xlab("selection coefficient (s)")

ggsave("plots/fitness_hist.pdf", width = 10, height = 8)

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

fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(!treatment=="YPD-30C-after-35cyc") %>%
  filter(`n_LOH-both`>0) %>%
  group_by(treatment) %>%
  summarise(max(mean_s))

fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  filter(`n_LOH-both`>0) %>%
  filter(`LOH_length-both`>100) %>%
  select(strain, treatment, mean_s, starts_with("LOH_length")) %>%
  pivot_longer(starts_with("LOH_length"), values_to = "LOH_length", names_to = "type") %>%
  mutate(type=str_remove(type, "^LOH_length-")) %>%
  filter(type=="both") %>%
  ggplot(aes(x=log10(`LOH_length`), y=mean_s)) +
  geom_smooth(method="lm", formula = y~x) +
  geom_jitter(alpha=0.2, width=0.1, height=0.01) +
  stat_poly_eq(parse=T, use_label("eq", "P"), formula=y~x, small.p=TRUE, label.y = "bottom") +
  theme_linedraw() +
  facet_wrap(vars(treatment), nrow=2, scales = "free_y") +
  xlab("log10(LOH length)") +
  ylab("selection coefficient (s)")

ggsave("plots/fitness_vs_LOHlen.pdf", width = 14, height = 6)

fit_tb %>%
  filter(!str_detect(strain, "^CNTRL-")) %>%
  filter(treatment %in% c(exp_treatments, "(average)")) %>%
  filter(`n_LOH-both`>0) %>%
  select(strain, treatment, mean_s, starts_with("LOH_length")) %>%
  pivot_longer(starts_with("LOH_length"), values_to = "LOH_length", names_to = "type") %>%
  filter(LOH_length>0) %>%
  mutate(type=str_remove(type, "^LOH_length-")) %>%
  filter(type!="both") %>%
  ggplot(aes(x=log10(`LOH_length`), y=mean_s, color=type)) +
  geom_smooth(method="lm") +
  geom_jitter(alpha=0.2, width=0.1, height=0.01) +
  stat_poly_eq(parse=T, use_label("eq", "P"), formula=y~x, small.p=TRUE, label.y = "bottom") +
  theme_linedraw() +
  facet_wrap(vars(treatment), nrow=2, scales = "free_y") +
  scale_color_muted() +
  xlab("log10(LOH length)") +
  ylab("selection coefficient (s)")

ggsave("plots/fitness_vs_LOHlen_by_type.pdf", width = 14, height = 6)

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
  facet_wrap(~treatment, scale="free_y", ncol = 2) +
  scale_color_muted(name="top 5 in") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("selection coefficient (s)")

ggsave("plots/barseq_tops_enrichment.pdf", height = 8, width = 10)

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
  ggplot(aes(x=log10(`LOH_length-both`), y=value)) +
  geom_point() +
  geom_smooth(method=lm, formula=y ~ exp(x)) +
  stat_fit_tidy(method = "lm", method.args = list(formula = y ~ exp(x)),
                aes(label = sprintf("y=%.3g^x+%.3g, p=%.3g",
                                    after_stat(`exp(x)_estimate`),
                                    after_stat(`Intercept_estimate`),
                                    after_stat(`exp(x)_p.value`)
                )
                )) +
  facet_wrap(~name, scales = "free_y") +
  theme_linedraw() +
  xlab("log10(LOH length)") +
  ylim(0, NA)

ggsave("plots/tradeoff.pdf", width = 14, height = 8)

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

ind_LOH_stats1 <-
  LOH_type_tb %>%
  group_by(type) %>%
  summarise(mean=mean(length),
            median=median(length),
            sd=sd(length),
            min=min(length),
            max=max(length))

ind_LOH_stats2 <-
  LOH_type_tb %>%
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
  group_by(type) %>%
  filter(!strain %in% control_strains) %>%
  summarise(across(-strain, .fns = list(mean=mean, median=median, sd=sd, min=min, max=max))) %>%
  knitr::kable()

LOH_long_tb %>%
  group_by(strain) %>%
  filter(sum(n_LOH)!=0) %>%
  group_by(type) %>%
  filter(!strain %in% control_strains) %>%
  summarise(across(-strain, .fns = list(mean=mean, median=median, sd=sd, min=min, max=max))) %>%
  knitr::kable()

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

cM_unit <- 30/1e6

pheno_tb <-
  fit_tb %>%
  select(strain, treatment, mean_s) %>%
  pivot_wider(names_from = treatment, values_from = mean_s) %>%
  column_to_rownames("strain")

geno <-
  read.table("bwa_haplotypecaller_finalvcf/runs.diploid.vcf.raw", header = TRUE)

map_tb <-
  read.table("bwa_haplotypecaller_finalvcf/runs.diploid.vcf.map", header = FALSE) %>%
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
operm <- scan1perm(pr, cross2$pheno, n_perm=1000)
summary(operm)

LOH_stat_plot1 <-
  LOH_long_tb %>%
  filter(!strain %in% control_strains) %>%
  mutate(strain=fct_drop(strain)) %>%
  complete(strain, type, fill=list(LOH_length=0)) %>%
  rename(`Cumulative LOH lenghth`=LOH_length) %>%
  group_by(strain, type) %>%
  mutate(
    bin = cut(`Cumulative LOH lenghth`, breaks = c(0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), include.lowest = TRUE)
  ) %>%
  ungroup() %>%
  count(bin, type) %>%
  complete(bin, type, fill=list(n=0)) %>%
  mutate(type=factor(type, levels=c("BY4741", "both", "SK1"))) %>%
  ggplot(aes(x=bin, y=n, fill=type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.5) +
  scale_fill_muted(name="LOH type") +
  theme_linedraw() +
  theme(legend.position = "none", plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "inch"),
        axis.text.x = element_text(angle=30, hjust = 1)) +
  xlab("Cumulative LOH Length") +
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
  complete(n_LOH, type, fill=list(n=0)) %>%
  group_by(n_LOH, type) %>%
  tally() %>%
  ungroup() %>%
  mutate(type=ifelse(type=="both", "combined", type)) %>%
  mutate(type=factor(type, levels=c("BY4741", "combined", "SK1"))) %>%
  ggplot(aes(x=n_LOH, y=n, fill=type)) +
  geom_col(position="dodge2", width=0.5) +
  scale_fill_muted(name="LOH type") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  theme_linedraw() +
  theme(legend.position = "none", plot.margin = unit(c(0.1, 0.1, 0.5, 0.1), "inch")) +
  xlab("Number of LOH") +
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
  geom_rect(data=cen_loc, aes(xmin = ci_lo-10000, xmax=ci_hi+10000, ymin = 0, ymax=-0.5), alpha=1, color="red") +
  geom_rect(data=SNP_loc, aes(xmin = ci_lo, xmax=ci_hi, ymin = -1, ymax=-0.5), alpha=1, inherit.aes = FALSE) +
  scale_y_discrete() +
  theme(legend.position = "none") +
  ylab("Depth (number of strains)")

LOH_stat_legend2 <-
  get_legend(
    ggplot(tibble(`LOH type`=c("BY4741", "SK1"), line=c("centromere", "SNP"))) +
      geom_line(aes(x=1, y=1, color=line)) +
      geom_raster(aes(fill=`LOH type`, x=1, y=1), alpha=0.3) +
      scale_color_manual(name="", values = c("centromere"="red", "SNP"="black")) +
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

dir.create("top5_list")
for (treatment in exp_treatments){
  rank_tb %>%
    rowwise() %>%
    mutate(pick=treatment %in% top5_in_vector) %>%
    filter(pick) %>%
    pull(strain) %>%
    write(paste0("tables/top5_",treatment, ".txt"))

  dir.create(paste0("top5_list/", treatment))
  from <-
    scan(paste0("tables/top5_",treatment, ".txt"), what=character()) %>%
    paste0("LOH_detect/LOH_minSNP-5_", ., ".bed")
  to=paste0("top5_list/", treatment)
  file.copy(from, to)
  from <-
    scan(paste0("tables/top5_",treatment, ".txt"), what=character()) %>%
    paste0("LOH_detect/revLOH_minSNP-5_", ., ".bed")
  to=paste0("top5_list/", treatment)
  file.copy(from, to)
}

top5_bed_tb <-
  tibble(treatment=exp_treatments) %>%
  mutate(fordbed=map(treatment, ~as_tibble(read_bed_graph(paste0("top5_list/LOH.", .x, ".bedgraph")))),
         revbed=map(treatment, ~as_tibble(read_bed_graph(paste0("top5_list/revLOH.", .x, ".bedgraph")))),
         bed=map2(fordbed, revbed, ~rbind(mutate(.x, type="SK1"), mutate(.y, type="BY4741")))) %>%
  select(treatment, bed) %>%
  unnest(bed) %>%
  select(chr=seqnames, ci_lo=start, ci_hi=end, pos=start, lodcolumn=score, type, treatment) %>%
  mutate(lodindex="depth", col="transparent", chr=str_remove(chr, "chr"))

top5_overlap_plots <-
map(exp_treatments, .f=function(x){
  top5_bed_tb %>%
    filter(treatment == x) %>%
    ggplot_peaks(map, gap=25, bgcolor="gray90", altbgcolor="gray85", tick_height = 0.8) +
    geom_rect(aes(xmin = ci_lo, xmax=ci_hi, ymin = 0, ymax=lodcolumn, fill=type), alpha=0.3) +
    scale_y_discrete() +
    ggtitle(x) +
    scale_fill_muted()
})

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

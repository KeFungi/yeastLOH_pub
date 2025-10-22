library(tidyverse)

strain_tb <-
  read_csv("tables/project_id.csv") %>%
  group_by(Description) %>%
  slice(1) %>%
  select(strain=Description) %>%
  ungroup()

het_tb <-
  read_tsv("tables/genomes.final.vcf.gz.het") %>%
  mutate(`homozygous rate`=`O(HOM)`/N_SITES) %>%
  select(strain=INDV, `homozygous rate`)

het_tb %>%
  ggplot(aes(x=`homozygous rate`)) +
  geom_histogram(breaks=seq(0,1,0.05)) +
  geom_vline(xintercept = 0.95, color="red") +
  theme_linedraw()

ggsave("plots/strain_homo_GT_rate.pdf", width = 8, height = 6)

missing_tb <-
  read_tsv("tables/genomes.final.vcf.gz.imiss") %>%
  select(strain=INDV, `missing genotype rate`=F_MISS)

missing_tb %>%
  ggplot(aes(x=`missing genotype rate`)) +
  geom_histogram(breaks=seq(0,0.17,0.01)) +
  geom_vline(xintercept = 0.05, color="red") +
  theme_linedraw()

ggsave("plots/strain_missing_GT_rate.pdf", width = 8, height = 6)

depth_tb <-
  read_tsv("tables/genomes.final.vcf.gz.idepth") %>%
  select(strain=INDV, `sequencing depth`=MEAN_DEPTH)

depth_tb %>%
  ggplot(aes(x=`sequencing depth`)) +
  geom_histogram()

barcode_tb <-
  read_csv("pub_tables/barcode_map.csv", col_names = c("strain", "barcode"))

meta_tb <-
  strain_tb %>%
  left_join(depth_tb) %>%
  left_join(missing_tb) %>%
  left_join(het_tb) %>%
  mutate(`included in analyses`= case_when(
    `missing genotype rate` > 0.05 ~ "NO, missing genotype rate > 5%",
    `homozygous rate` > 0.5 ~ "NO, haploid",
    strain %in% scan("pub_tables/dup_barcode.txt", what=character()) ~ "NO, duplicated barcode",
    strain %in% scan("pub_tables/unidentifiable_barcode.txt", what=character()) ~ "NO, unidentifiable barcode",
    strain == "P3-2C" ~ "NO, barcode underrepresented (<0.01%)",
    TRUE ~ "YES")
  )

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
  filter(strain %in% c(good_strains, "P3-2C"))

loh_pos_tb <-
  LOH_type_tb %>%
  mutate(pos_text=paste0(chr, ":", start+1, "-", end),
         length=end-start
  ) %>%
  group_by(strain, type) %>%
  summarise(n_LOH=n(), LOH_length=sum(length), concat_pos_text=paste0(pos_text, collapse = "; ")) %>%
  mutate(type_pos_text=paste0("(", type, ") ", concat_pos_text, ".")) %>%
  arrange(type) %>%
  group_by(strain) %>%
  summarise(`number of LOH`=sum(n_LOH),
            `total LOH length`=sum(LOH_length),
            `LOH position`=paste0(type_pos_text, collapse = " ")
  )

meta_readable_tb <-
  meta_tb %>%
  left_join(loh_pos_tb) %>%
  mutate(`number of LOH`=ifelse(is.na(`number of LOH`) & `included in analyses`=="YES", 0, `number of LOH`),
         `total LOH length`=ifelse(is.na(`total LOH length`) & `included in analyses`=="YES", 0, `total LOH length`)
  ) %>%
  mutate(`included in analyses`=fct_relevel(`included in analyses`, "YES")) %>%
  arrange(`included in analyses`) %>%
  relocate(-`included in analyses`) %>%
  left_join(barcode_tb) %>%
  mutate(barcode=ifelse(`included in analyses` %in% c("YES", "NO, barcode underrepresented (<0.01%)"), barcode, NA)) %>%
  relocate(strain, barcode)


write_csv(meta_readable_tb, "tables/meta_table.csv")

LOH_type_tb %>%
  filter(type=="SK1") %>%
  select(chr, start, end, strain) %>%
  write_tsv("tables/SK1_bedgraph_input.bed", col_names = FALSE)

LOH_type_tb %>%
  filter(type=="BY4741") %>%
  select(chr, start, end, strain) %>%
  write_tsv("tables/BY4741_bedgraph_input.bed", col_names = FALSE)


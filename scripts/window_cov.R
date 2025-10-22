library(tidyverse)
library(RcppRoll)

files <-
  list.files("bwa_mapping_depth", full.names = TRUE) %>%
  str_subset("\\.depth\\.txt")

cov_tb <- tibble()

for (f in files){
  id <- str_match(f, "bwa_mapping_depth\\/(.+)\\.genome\\.info\\.bam\\.depth\\.txt")[,2]
  tbi <-
    read_tsv(f, col_types="cii", col_names=FALSE) %>%
    group_by(X1) %>%
    mutate(!!id:=RcppRoll::roll_mean(X3, 5000L, by=5000L, fill=NA, na.rm=TRUE)) %>%
    ungroup() %>%
    filter(!is.na(.data[[id]]))
  if (nrow(cov_tb) == 0) {
    cov_tb <- select(tbi, chr=X1, position=X2, id)
  } else {
    cov_tb <- cbind(cov_tb, select(tbi, id))
  }
}

cov_tb %>%
  write_csv("tables/window_depth.csv")

cov_std <-
  read_csv("tables/window_depth.csv")

cov_std_tb <-
  cov_std %>%
  group_by(chr) %>%
  slice_tail(n=-10) %>%
  slice_head(n=-10) %>%
  ungroup() %>%
  mutate(across(c(-chr, -position), ~.x/mean(.x, na.rm=TRUE))) %>%
  rowwise() %>%
  mutate(across(c(-chr, -position), ~.x/mean(c_across(c(-chr, -position)), na.rm=TRUE))) %>%
  ungroup()

cov_std_chr_tb <-
  cov_std_tb %>%
  group_by(chr) %>%
  summarise(across(c(-position), ~mean(.x, na.rm=TRUE)))

cov_std_chr_tb %>%
  pivot_longer(-chr) %>%
  ggplot(aes(x=value)) +
  geom_histogram(breaks=seq(0.5, 1.6, 0.05)) +
  geom_vline(xintercept = c(0.75, 1.25)) +
  facet_wrap(facets = ~chr) +
  xlab("standarlized depth (cutoff=0.75, 1.25)")

cov_std_chr_tb %>%
  pivot_longer(-chr) %>%
  filter(value<0.75|value>1.25) %>%
  mutate(aneuploidy=
           case_when(value<0.75 ~ "<1",
                     value>1.25 ~ ">1"
           )
         ) %>%
  arrange(chr) %>%
  group_by(name) %>%
  select(-value, -chr) %>% unique() %>%
  rename(strain=name) %>%
  write_csv("tables/aneuploidy.txt")

plots <- list()
for (i in 1L:(dim(cov_std_tb)[2] %/% 10)+1) {
  i_s <- i*10-7
  i_e <- i*10+2
  if (i_e>dim(cov_std_tb)[2]) {i_e<-dim(cov_std_tb)[2]}

  plots[[i]] <-
    cov_std_tb %>%
    ungroup() %>%
    select(c(c(1, 2), c(i_s:i_e))) %>%
    mutate(loc=row_number()) %>%
    pivot_longer(c(-chr, -position, -loc)) %>%
    filter(chr=="chrIII") %>%
    ggplot(aes(y=value, x=loc, color=chr)) +
    geom_point() +
    geom_hline(yintercept = 1, lty=1) +
    geom_hline(yintercept = 1.5, lty=2) +
    geom_hline(yintercept = 0.5, lty=2) +
    facet_wrap(facets = "name", ncol=1)
}

+
  geom_histogram(breaks=seq(0.5, 1.6, 0.05)) +
  geom_vline(xintercept = c(0.6, 1.35))

cov_std_tb %>%
  pivot_longer(-chr, -position) %>%
  filter(value<0.6|value>1.35) %>% arrange(chr) %>%
  mutate(aneuploidy=
           case_when(value<0.6 ~ "<1",
                     value>1.35 ~ ">1"
                     )
         ) %>%
  group_by(name) %>%
  summarise(category=list(unique(aneuploidy))) %>% view()

            pull(name) %>% unique() %>%
  write_lines("tables/aneuploidy.txt")

cov_std_tb %>%
  pivot_longer(-chr) %>%
  ungroup() %>%
  filter(value<0.6|value>1.35) %>%
  mutate(n=ifelse(value<0.6, 0.5, 1.5)) %>%
  select(strain=name, n) ->
  cov_lab

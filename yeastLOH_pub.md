# reference
```bash
mkdir -p ref tables thinned_reads bowtie_barcode_mapping bwa_mapping bwa_haplotypecaller_gvcf bwa_haplotypecaller_log bwa_haplotypecaller_finalvcf

wget "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-5-1_20240529.tgz"
tar -zxf S288C_reference_genome_R64-5-1_20240529.tgz
seqkit grep -n -r -p ".*\[chromosome=.+\].*" S288C_reference_genome_R64-5-1_20240529/S288C_reference_sequence_R64-5-1_20240529.fsa.gz | seqkit replace -p ".*\[chromosome=(.+)\].*" -r "chr\$1" | bgzip > ref/S288C.chr.fasta.gz

samtools faidx ref/S288C.chr.fasta.gz
picard CreateSequenceDictionary -R ref/S288C.chr.fasta.gz -O ref/S288C.chr.dict

gzip -c -d S288C_reference_genome_R64-5-1_20240529/saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz | awk -F $'\t' -v OFS=$'\t' '$3 == "centromere" {print $1, $4-1, $5, $3}' > ref/centromere.bed
```

# copy data
```bash
mkdir cutadapt

cp -r /nfs/turbo/lsa-tyjames/mycology/next_gen_seq_data/Yeast_LOH/20240722_LOH_wholegenome1/11314-YK/ .
cp -r /nfs/turbo/lsa-tyjames/mycology/next_gen_seq_data/Yeast_LOH/20240821_LOH_wholegenome2/11790-YK/ .

csvtk concat 11314-YK/DemuxStats_11314-YK.csv 11790-YK/DemuxStats_11790-YK.csv  | csvtk mutate2 -e '$Project + "/" + $Sample_ID' -n path | csvtk sort -k Description:N -k '# Reads:nr' | csvtk uniq -f Description > tables/project_id.csv
```

# cutadapt
```bash
barpath=( $(csvtk cut -f 6 -U tables/project_id.csv) )
names=( $(csvtk cut -f 3 -U tables/project_id.csv) )

for i in ${!names[@]}; do
	f1=$(find . | grep "${barpath[$i]}" | grep "_R1_")
	f2=$(find . | grep "${barpath[$i]}" | grep "_R2_")
    o1="cutadapt/${names[$i]}_R1.trim.fastq.gz"
    o2="cutadapt/${names[$i]}_R2.trim.fastq.gz"
    sbatch -J cutadapt -A tyjames1 --mem=2g -o cutadapt/${names[$i]}.log --wrap="cutadapt -m 100 -q 20 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o $o1 -p $o2 $f1 $f2"
done
```

# bwa
```bash
bwa index ref/S288C.chr.fasta.gz

names=( $(csvtk cut -f 3 -U tables/project_id.csv) )
for i in ${!names[@]}; do
    f1="cutadapt/${names[$i]}_R1.trim.fastq.gz"
    f2="cutadapt/${names[$i]}_R2.trim.fastq.gz"
    sbatch -J bwa -A tyjames0 -t 24:00:00 --mem=2g -o bwa_mapping/${names[$i]}.genome.log --wrap="bwa mem ref/S288C.chr.fasta.gz $f1 $f2 > bwa_mapping/${names[$i]}.genome.sam"
done

names=( $(csvtk cut -f 3 -U tables/project_id.csv) )
for i in ${!names[@]}; do
  samfile=bwa_mapping/${names[$i]}.genome.sam
  bamfile=bwa_mapping/${names[$i]}.genome.sorted.bam
  logfile=bwa_mapping/${names[$i]}.sort.log
  if [[ ! -f ${bamfile}.bai ]]; then
      sbatch -J sort_bam -A tyjames0 -t 4:00:00 --mem=2g -o $logfile --wrap="samtools view -F 4 -h $samfile | samtools sort -O BAM -o $bamfile -; samtools index $bamfile"
  fi
done

names=( $(csvtk cut -f 3 -U tables/project_id.csv) )
for i in ${!names[@]}; do
  bamfile=bwa_mapping/${names[$i]}.genome.sorted.bam
  newbamfile=bwa_mapping/${names[$i]}.genome.info.bam
  logifle=bwa_mapping/${names[$i]}.info.log
  if [[ ! -f bwa_mapping/${names[$i]}.genome.info.bai ]]; then
    sbatch -J info_bam -A tyjames0 --mem=2g -o $logifle --wrap="picard AddOrReplaceReadGroups I=${bamfile} O=${newbamfile} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${names[$i]} TMP_DIR=/scratch/tyjames_root/tyjames0/yihongke/; picard BuildBamIndex INPUT=${newbamfile}"
  fi
done

names=( $(csvtk cut -f 3 -U tables/project_id.csv) )
mkdir -p bwa_mapping_depth
for i in ${!names[@]}; do
  bamfile=bwa_mapping/${names[$i]}.genome.info.bam
  logfile=bwa_mapping_depth/${names[$i]}.depth.log
  depthfile=bwa_mapping_depth/${names[$i]}.genome.info.bam.depth.txt
  if [[ -f ${bamfile%\.bam}.bai ]]; then
      sbatch -J depth_bam -A tyjames0 --mem=2g -o $logfile --wrap="samtools depth -aa $bamfile > $depthfile"
  fi
done

names=( $(csvtk cut -f 3 -U tables/project_id.csv) )
for i in ${!names[@]}; do
  inbam=bwa_mapping/${names[$i]}.genome.info.bam
  outvcf=bwa_haplotypecaller_gvcf/${names[$i]}.g.vcf.gz
  logfile=bwa_haplotypecaller_log/${names[$i]}_haplotypecaller.log
  if [[ ! -f bwa_haplotypecaller_gvcf/${names[$i]}.g.vcf.gz.tbi ]]; then
  sbatch -J hapcaller -A tyjames0 -t 24:00:00 --mem=8g -c 1 -o ${logfile} \
    --wrap="gatk --java-options '-Xmx8g' HaplotypeCaller \
    --tmp-dir /scratch/tyjames_root/tyjames0/yihongke/ -R ref/S288C.chr.fasta.gz \
    --native-pair-hmm-threads 1 \
    -I $inbam -O $outvcf -ploidy 2 -ERC GVCF"
  fi
done

runvflag=$(awk -F, '{print $3}' tables/project_id.csv | awk 'NR>1 {print "--variant bwa_haplotypecaller_gvcf/" $1 ".g.vcf.gz "}' | xargs echo)

sbatch -J combine_gvcf -t 48:00:00 -A tyjames0 -c 4 --mem=8g --wrap="gatk CombineGVCFs -R ref/S288C.chr.fasta.gz $runvflag -O bwa_haplotypecaller_gvcf/runs.g.vcf.gz"
sbatch -J genotype_gvcf -t 48:00:00 -A tyjames0 -c 4 --mem=8g --wrap="gatk GenotypeGVCFs -R ref/S288C.chr.fasta.gz -V bwa_haplotypecaller_gvcf/runs.g.vcf.gz -O bwa_haplotypecaller_gvcf/runs.raw.vcf.gz"

```

# find barcode
```bash
bowtie2-build ref/SceIKnockIn_20bp.fasta ref/SceIKnockIn_20bp
bowtie2-build ref/SceIKnockIn_19bp.fasta ref/SceIKnockIn_19bp
bowtie2-build ref/SceIKnockIn_18bp.fasta ref/SceIKnockIn_18bp
bowtie2-build ref/SceIKnockIn_17bp.fasta ref/SceIKnockIn_17bp

names=( $(csvtk cut -f 3 -U tables/project_id.csv) )
for i in ${!names[@]}; do
    f1="cutadapt/${names[$i]}_R1.trim.fastq.gz"
    f2="cutadapt/${names[$i]}_R2.trim.fastq.gz"
    sbatch -J barcod_bowtie -A tyjames0 --mem=2g -o bowtie_barcode_mapping/${names[$i]}.SceI_20bp.log --wrap="bowtie2 --no-unal --n-ceil L,0,0.2 --np 0 -x ref/SceIKnockIn_20bp -1 $f1 -2 $f2 > bowtie_barcode_mapping/${names[$i]}.SceI_20bp.sam"
    sbatch -J barcod_bowtie -A tyjames0 --mem=2g -o bowtie_barcode_mapping/${names[$i]}.SceI_19bp.log --wrap="bowtie2 --no-unal --n-ceil L,0,0.2 --np 0 -x ref/SceIKnockIn_19bp -1 $f1 -2 $f2 > bowtie_barcode_mapping/${names[$i]}.SceI_19bp.sam"
    sbatch -J barcod_bowtie -A tyjames0 --mem=2g -o bowtie_barcode_mapping/${names[$i]}.SceI_18bp.log --wrap="bowtie2 --no-unal --n-ceil L,0,0.2 --np 0 -x ref/SceIKnockIn_18bp -1 $f1 -2 $f2 > bowtie_barcode_mapping/${names[$i]}.SceI_18bp.sam"
    sbatch -J barcod_bowtie -A tyjames0 --mem=2g -o bowtie_barcode_mapping/${names[$i]}.SceI_18bp.log --wrap="bowtie2 --no-unal --n-ceil L,0,0.2 --np 0 -x ref/SceIKnockIn_17bp -1 $f1 -2 $f2 > bowtie_barcode_mapping/${names[$i]}.SceI_17bp.sam"
done

names=( $(csvtk cut -f 3 -U tables/project_id.csv) )
for i in ${!names[@]}; do
    samtools sort bowtie_barcode_mapping/${names[$i]}.SceI_20bp.sam > bowtie_barcode_mapping/${names[$i]}.SceI_20bp.bam
    samtools sort bowtie_barcode_mapping/${names[$i]}.SceI_19bp.sam > bowtie_barcode_mapping/${names[$i]}.SceI_19bp.bam
    samtools sort bowtie_barcode_mapping/${names[$i]}.SceI_18bp.sam > bowtie_barcode_mapping/${names[$i]}.SceI_18bp.bam
    samtools sort bowtie_barcode_mapping/${names[$i]}.SceI_17bp.sam > bowtie_barcode_mapping/${names[$i]}.SceI_17bp.bam
    bcftools mpileup -f ref/SceIKnockIn_20bp.fasta bowtie_barcode_mapping/${names[$i]}.SceI_20bp.bam | bcftools call -c -M -O z --ploidy 1 > bowtie_barcode_mapping/${names[$i]}.SceI_20bp.vcf.gz
    bcftools mpileup -f ref/SceIKnockIn_19bp.fasta bowtie_barcode_mapping/${names[$i]}.SceI_19bp.bam | bcftools call -c -M -O z --ploidy 1 > bowtie_barcode_mapping/${names[$i]}.SceI_19bp.vcf.gz
    bcftools mpileup -f ref/SceIKnockIn_18bp.fasta bowtie_barcode_mapping/${names[$i]}.SceI_18bp.bam | bcftools call -c -M -O z --ploidy 1 > bowtie_barcode_mapping/${names[$i]}.SceI_18bp.vcf.gz
    bcftools mpileup -f ref/SceIKnockIn_17bp.fasta bowtie_barcode_mapping/${names[$i]}.SceI_17bp.bam | bcftools call -c -M -O z --ploidy 1 > bowtie_barcode_mapping/${names[$i]}.SceI_17bp.vcf.gz
    bcftools index bowtie_barcode_mapping/${names[$i]}.SceI_20bp.vcf.gz
    bcftools index bowtie_barcode_mapping/${names[$i]}.SceI_19bp.vcf.gz
    bcftools index bowtie_barcode_mapping/${names[$i]}.SceI_18bp.vcf.gz
    bcftools index bowtie_barcode_mapping/${names[$i]}.SceI_17bp.vcf.gz
done

> tables/barcode_map.csv
> tables/barcode_20bp_map.csv
> tables/barcode_19bp_map.csv
> tables/barcode_18bp_map.csv
> tables/barcode_17bp_map.csv

names=( $(csvtk cut -f 3 -U tables/project_id.csv) )
for i in ${!names[@]}; do
	bar=$(bcftools view -H -r SceIKnockIn_20bp:2386-2405 bowtie_barcode_mapping/${names[$i]}.SceI_20bp.vcf.gz | awk '{print $5}' | tr -d '\n')
    echo "${names[$i]},${bar}" | tee /dev/tty >> tables/barcode_20bp_map.csv
    bar=$(bcftools view -H -r SceIKnockIn_19bp:2386-2404 bowtie_barcode_mapping/${names[$i]}.SceI_19bp.vcf.gz | awk '{print $5}' | tr -d '\n')
    echo "${names[$i]},${bar}" | tee /dev/tty >> tables/barcode_19bp_map.csv
    bar=$(bcftools view -H -r SceIKnockIn_18bp:2386-2403 bowtie_barcode_mapping/${names[$i]}.SceI_18bp.vcf.gz | awk '{print $5}' | tr -d '\n')
    echo "${names[$i]},${bar}" | tee /dev/tty >> tables/barcode_18bp_map.csv
    bar=$(bcftools view -H -r SceIKnockIn_17bp:2386-2402 bowtie_barcode_mapping/${names[$i]}.SceI_17bp.vcf.gz | awk '{print $5}' | tr -d '\n')
    echo "${names[$i]},${bar}" | tee /dev/tty >> tables/barcode_17bp_map.csv
done

grep -h -o -P '\w+-\w+,[A|T|C|G]+$' tables/barcode_20bp_map.csv tables/barcode_19bp_map.csv tables/barcode_18bp_map.csv tables/barcode_17bp_map.csv |
sort > tables/barcode_map.csv

awk -F, '{print $1}' tables/barcode_map.csv | uniq -c | awk '{if ($1>1) print}' # multiple possible barcodes
cat tables/barcode_map.csv | awk -F, '{print $1}' | sort | comm -13 - <(awk -F, 'NR>1{print $3}' tables/project_id.csv | sort) # missing barcodes
cat tables/barcode_map.csv | awk -F, '{print $1}' | sort | comm -3 - <(awk -F, 'NR>1{print $3}' tables/project_id.csv | sort) # missing or duplicated strains
```

manually edit tables/barcode_map.csv to pub_tables/barcode_map.csv


```R
library(tidyverse)
read_csv("pub_tables/barcode_map.csv", col_names = c("strain", "barcode")) %>% mutate(bar_length=nchar(barcode)) %>% write_csv("tables/barcode_map_parsed.csv") #parse barcode file
(read_csv("tables/barcode_map_parsed.csv") %>% mutate(bar_length=as.factor(bar_length)) %>% ggplot(aes(x=bar_length)) + geom_bar() + theme_light()) %>% ggsave("plots/barcode_hist.png", ., width=3, height=2) #barcode length histogram
read_csv("tables/barcode_map_parsed.csv") %>% filter(!is.na(barcode)) %>% group_by(barcode) %>% arrange(barcode) %>% mutate(n=n()) %>% filter(n>1) %>% print(n=100) #duplicated barcodes

read_csv("tables/barcode_map_parsed.csv") %>% filter(!is.na(barcode)) %>% group_by(barcode) %>% arrange(barcode) %>% mutate(n=n()) %>% filter(n>1) %>% pull(strain) %>% unique() %>% sort()
```

# A tibble: 31 × 4
# Groups:   barcode [15]
   strain barcode              bar_length     n
   <chr>  <chr>                     <dbl> <int>
 1 P3-3A  AAAGACAAAATACAATAAGA         20     2
 2 P3-3D  AAAGACAAAATACAATAAGA         20     2
 3 P3-3H  AACCGTTTCCCATAAGAGGC         20     2
 4 P3-4A  AACCGTTTCCCATAAGAGGC         20     2
 5 P1-12H ACAGGACTGAGAACGAACGA         20     2
 6 P1-7H  ACAGGACTGAGAACGAACGA         20     2
 7 P2-4D  ATGTGTATTTCTGGAAAAAA         20     2
 8 P3-4D  ATGTGTATTTCTGGAAAAAA         20     2
 9 P2-10G CAAACAAGGAGAAGTAAAAC         20     2
10 P2-9E  CAAACAAGGAGAAGTAAAAC         20     2
11 P1-3H  CAATAGCCCAAGGGATAATA         20     2
12 P2-11C CAATAGCCCAAGGGATAATA         20     2
13 P1-1H  CCACAGGGCTCTAGTCCAAA         20     2
14 P1-9E  CCACAGGGCTCTAGTCCAAA         20     2
15 P3-1H  CCTGATGGTCGCAAATGTTC         20     2
16 P3-4G  CCTGATGGTCGCAAATGTTC         20     2
17 P1-1D  CTACCCCCGTGACGAGTACC         20     2
18 P1-7C  CTACCCCCGTGACGAGTACC         20     2
19 P2-6G  GATCACAGAAAAATTTAAAT         20     2
20 P2-8G  GATCACAGAAAAATTTAAAT         20     2
21 P1-7D  GCACTGAGGCTCTACGGCGT         20     3
22 P1-7E  GCACTGAGGCTCTACGGCGT         20     3
23 P1-8C  GCACTGAGGCTCTACGGCGT         20     3
24 P3-7H  GGCAGTTGGCGGCGAAGAAG         20     2
25 P3-8H  GGCAGTTGGCGGCGAAGAAG         20     2
26 P1-11C GGGAAGACAAAACGAAGGAG         20     2
27 P1-6C  GGGAAGACAAAACGAAGGAG         20     2
28 P2-10H GGGCAGTAAAGGGGGGGAGT         20     2
29 P2-11H GGGCAGTAAAGGGGGGGAGT         20     2
30 P3-1C  GTGCTTGGCCTAATGAACAA         20     2
31 P3-5A  GTGCTTGGCCTAATGAACAA         20     2

# reference: SK1
```bash
fasterq-dump SRR1569638 -O ref/
fasterq-dump SRR1569870 -O ref/

f1="ref/SRR1569638_1.fastq"
f2="ref/SRR1569638_2.fastq"
o1="ref/SRR1569638_1.trim.fastq"
o2="ref/SRR1569638_2.trim.fastq"
sbatch -J SK1_cutadapt -A tyjames0 --mem=2g -o ref/SK1.cutadapt.log --wrap="cutadapt -m 100 -q 20 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o $o1 -p $o2 $f1 $f2"

f1="ref/SRR1569870_1.fastq"
f2="ref/SRR1569870_2.fastq"
o1="ref/SRR1569870_1.trim.fastq"
o2="ref/SRR1569870_2.trim.fastq"
sbatch -J BY4741_cutadapt -A tyjames0 --mem=2g -o ref/BY4741.cutadapt.log --wrap="cutadapt -m 100 -q 20 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o $o1 -p $o2 $f1 $f2"

f1="ref/SRR1569638_1.trim.fastq"
f2="ref/SRR1569638_2.trim.fastq"
logfile="ref/SK1.bwa.log"
sbatch -J SK1_bwa -A tyjames0 --mem=2g -o $logfile --wrap="bwa mem ref/S288C.chr.fasta.gz $f1 $f2 > ref/SK1.genome.sam"

f1="ref/SRR1569870_1.trim.fastq"
f2="ref/SRR1569870_2.trim.fastq"
logfile="ref/BY4741.bwa.log"
sbatch -J BY4741_bwa -A tyjames0 --mem=2g -o $logfile --wrap="bwa mem ref/S288C.chr.fasta.gz $f1 $f2 > ref/BY4741.genome.sam"


samfile=ref/SK1.genome.sam
bamfile=ref/SK1.genome.sorted.bam
logfile=ref/SK1.genome.sort.log
sbatch -J SK1_sort_bam -A tyjames0 --mem=2g -o $logfile --wrap="samtools view -F 4 -h $samfile | samtools sort -O BAM -o $bamfile -; samtools index $bamfile"

samfile=ref/BY4741.genome.sam
bamfile=ref/BY4741.genome.sorted.bam
logfile=ref/BY4741.genome.sort.log
sbatch -J BY4741_sort_bam -A tyjames0 --mem=2g -o $logfile --wrap="samtools view -F 4 -h $samfile | samtools sort -O BAM -o $bamfile -; samtools index $bamfile"


bamfile=ref/SK1.genome.sorted.bam
newbamfile=ref/SK1.genome.info.bam
logfile=ref/SK1.genome.info.log
sbatch -J SK1_info_bam -A tyjames0 --mem=2g -o $logfile --wrap="picard AddOrReplaceReadGroups I=${bamfile} O=${newbamfile} RGLB=ref RGPL=illumina RGPU=unit1 RGSM=SK1 TMP_DIR=/scratch/tyjames_root/tyjames0/yihongke/; picard BuildBamIndex INPUT=${newbamfile}"

bamfile=ref/BY4741.genome.sorted.bam
newbamfile=ref/BY4741.genome.info.bam
logfile=ref/BY4741.genome.info.log
sbatch -J BY4741_info_bam -A tyjames0 --mem=2g -o $logfile --wrap="picard AddOrReplaceReadGroups I=${bamfile} O=${newbamfile} RGLB=ref RGPL=illumina RGPU=unit1 RGSM=BY4741 TMP_DIR=/scratch/tyjames_root/tyjames0/yihongke/; picard BuildBamIndex INPUT=${newbamfile}"


inbam=ref/SK1.genome.info.bam
outvcf=ref/SK1.g.vcf.gz
logfile=ref/SK1.haplotypecaller.log
sbatch -J SK1_hapcaller -A tyjames0 -t 48:00:00 --mem=8g -c 1 -o ${logfile} --wrap="gatk --java-options '-Xmx8g' HaplotypeCaller --tmp-dir /scratch/tyjames_root/tyjames0/yihongke/ -R ref/S288C.chr.fasta.gz --native-pair-hmm-threads 1 -I $inbam -O $outvcf -ploidy 2 -ERC GVCF"

inbam=ref/BY4741.genome.info.bam
outvcf=ref/BY4741.g.vcf.gz
logfile=ref/BY4741.haplotypecaller.log
sbatch -J BY4741_hapcaller -A tyjames0 -t 48:00:00 --mem=8g -c 1 -o ${logfile} --wrap="gatk --java-options '-Xmx8g' HaplotypeCaller --tmp-dir /scratch/tyjames_root/tyjames0/yihongke/ -R ref/S288C.chr.fasta.gz --native-pair-hmm-threads 1 -I $inbam -O $outvcf -ploidy 2 -ERC GVCF"


sbatch -J combine_gvcf -t 48:00:00 -A tyjames0 -c 4 --mem=8g --wrap="gatk CombineGVCFs -R ref/S288C.chr.fasta.gz --variant ref/BY4741.g.vcf.gz --variant ref/SK1.g.vcf.gz -O ref/parents.g.vcf.gz"

gatk GenotypeGVCFs -R ref/S288C.chr.fasta.gz -V ref/parents.g.vcf.gz -O ref/parents.raw.vcf.gz

vcftools --gzvcf bwa_haplotypecaller_gvcf/runs.hardfilter.vcf.gz --depth

gatk VariantFiltration -R ref/S288C.chr.fasta.gz -V ref/parents.raw.vcf.gz -O ref/parents.hardfilter.vcf.gz --filter-name hard_filters --filter-expression 'BaseQRankSum<-5.0 || BaseQRankSum>5.0 || MQ<59.0 || QD<5.0 || FS>5.0 || ReadPosRankSum<-5.0 || ReadPosRankSum>5.0 || SOR>5.0'

gatk SelectVariants -R ref/S288C.chr.fasta.gz -V ref/parents.hardfilter.vcf.gz -O ref/parents.hardfilter.selected.vcf.gz --exclude-filtered --select-type-to-include SNP --restrict-alleles-to BIALLELIC --exclude-non-variants

gatk VariantFiltration -R ref/S288C.chr.fasta.gz -V ref/parents.hardfilter.selected.vcf.gz -O ref/parents.DPfilter.vcf.gz --genotype-filter-name 'isHetFilter' --genotype-filter-expression 'isHet == 1' --genotype-filter-name 'lowGT' --genotype-filter-expression 'DP < 20.0 || GQ <99'

gatk SelectVariants -R ref/S288C.chr.fasta.gz -V ref/parents.DPfilter.vcf.gz -O ref/parents.DPfilter.selected.vcf.gz --exclude-filtered --select-type-to-include SNP --restrict-alleles-to BIALLELIC --exclude-non-variants --set-filtered-gt-to-nocall

# make parent rsid list for filtering experiment genomes
echo -e "#CHROM\tPOS\tALT\tID" > ref/parents_rsid.tsv
bcftools query -f '%CHROM\t%POS\t%ALT[\t%GT]\n' ref/parents.DPfilter.selected.vcf.gz | awk '($4=="0/0" || $4=="0|0") && ($5=="1/1" || $5=="1|1") {print $1 "\t" $2 "\t" $3 "\trs" NR}' >> ref/parents_rsid.tsv
bgzip -f ref/parents_rsid.tsv
tabix -s 1 -b 2 -e 2 ref/parents_rsid.tsv.gz
gzip -c -d ref/parents_rsid.tsv.gz | awk '{print $4}' > ref/parents_rsid.list

cp ref/parents.raw.vcf.gz* bwa_haplotypecaller_finalvcf
cp ref/parents.DPfilter.selected.vcf.gz* bwa_haplotypecaller_finalvcf

```

# filter variants
```bash
# quality filter - flag

vcftools --gzvcf bwa_haplotypecaller_gvcf/runs.hardfilter.vcf.gz --depth --out bwa_haplotypecaller_gvcf/runs.raw.vcf.gz

gatk VariantFiltration -R ref/S288C.chr.fasta.gz -V bwa_haplotypecaller_gvcf/runs.raw.vcf.gz -O bwa_haplotypecaller_gvcf/runs.hardfilter.vcf.gz --filter-name hard_filters --filter-expression 'BaseQRankSum<-5.0 || BaseQRankSum>5.0 || MQ<59.0 || QD<5.0 || FS>5.0 || ReadPosRankSum<-5.0 || ReadPosRankSum>5.0 || SOR>5.0' --genotype-filter-name 'lowGT' --genotype-filter-expression 'DP < 20.0 || GQ <99'


# quality filter - exclude
gatk SelectVariants -R ref/S288C.chr.fasta.gz -V bwa_haplotypecaller_gvcf/runs.hardfilter.vcf.gz -O bwa_haplotypecaller_gvcf/runs.hardfilter.selected.vcf.gz --exclude-filtered --select-type-to-include SNP --restrict-alleles-to BIALLELIC --set-filtered-gt-to-nocall

# AF filter - flag
gatk VariantFiltration -R ref/S288C.chr.fasta.gz -V bwa_haplotypecaller_gvcf/runs.hardfilter.selected.vcf.gz -O bwa_haplotypecaller_gvcf/runs.GTfilter.vcf.gz --genotype-filter-name 'lowGT' --genotype-filter-expression 'DP < 20.0 || GQ <99' --filter-name AF_filters --filter-expression 'AN<386.0 || AF<0.45 || AF> 0.55' 

# AF filter - exclude
gatk SelectVariants -R ref/S288C.chr.fasta.gz -V bwa_haplotypecaller_gvcf/runs.GTfilter.vcf.gz -O bwa_haplotypecaller_gvcf/runs.GTfilter.selected.vcf.gz --exclude-filtered --select-type-to-include SNP --set-filtered-gt-to-nocall --restrict-alleles-to BIALLELIC --exclude-non-variants

# name variants by position
bcftools annotate -a ref/parents_rsid.tsv.gz -Oz -c CHROM,POS,ALT,ID -o bwa_haplotypecaller_gvcf/runs.GTfilter.selected.rsid.vcf.gz bwa_haplotypecaller_gvcf/runs.GTfilter.selected.vcf.gz
gatk IndexFeatureFile -I bwa_haplotypecaller_gvcf/runs.GTfilter.selected.rsid.vcf.gz


# filter by detected SK1-BY4741 variants
gatk SelectVariants -R ref/S288C.chr.fasta.gz -V bwa_haplotypecaller_gvcf/runs.GTfilter.selected.rsid.vcf.gz -O bwa_haplotypecaller_gvcf/genomes.final.vcf.gz --exclude-filtered --select-type-to-include SNP --set-filtered-gt-to-nocall --restrict-alleles-to BIALLELIC --keep-ids ref/parents_rsid.list 

cp bwa_haplotypecaller_gvcf/genomes.final.vcf.gz* bwa_haplotypecaller_finalvcf
cp bwa_haplotypecaller_gvcf/runs.raw.vcf.gz* bwa_haplotypecaller_finalvcf
cp bwa_haplotypecaller_gvcf/runs.GTfilter.selected.vcf.gz* bwa_haplotypecaller_finalvcf

bcftools index -f bwa_haplotypecaller_finalvcf/genomes.final.vcf.gz
vcftools --gzvcf bwa_haplotypecaller_finalvcf/genomes.final.vcf.gz --missing-indv --out tables/genomes.final.vcf.gz
vcftools --gzvcf bwa_haplotypecaller_finalvcf/genomes.final.vcf.gz --het --out tables/genomes.final.vcf.gz
vcftools --gzvcf bwa_haplotypecaller_finalvcf/genomes.final.vcf.gz --depth --out tables/genomes.final.vcf.gz

# exclude haploid samples
gatk SelectVariants  --exclude-sample-name P3-4H --exclude-sample-name CNTRL-1F --exclude-sample-name P1-11F --exclude-sample-name P1-7B --exclude-sample-name P2-2D -R ref/S288C.chr.fasta.gz -V bwa_haplotypecaller_finalvcf/genomes.final.vcf.gz -O bwa_haplotypecaller_finalvcf/runs.diploid.vcf.gz

bcftools query -H -f '%CHROM\t%POS\t[%GT\t]\n' bwa_haplotypecaller_finalvcf/runs.diploid.vcf.gz | gzip -c > bwa_haplotypecaller_finalvcf/runs.diploid.vcf.tsv.gz
```

# Appendix: strains excluded
duplicated barcode (pub_table/dup_barcode.txt): 31
P1-1D
P1-1H
P1-3H
P1-6C
P1-7C
P1-7D
P1-7E
P1-7H
P1-8C
P1-9E
P1-11C
P1-12H
P2-4D
P2-6G
P2-8G
P2-9E
P2-10G
P2-10H
P2-11C
P2-11H
P3-1C
P3-1H
P3-3A
P3-3D
P3-3H
P3-4A
P3-4D
P3-4G
P3-5A
P3-7H
P3-8H

mixed barcode (pub_table/unidentifiable_barcode.txt): 8
NFS-1F
NFS-1H
NFS-2E
NFS-8C
P1-11G
P2-1F
P2-8B
P3-6H

no insertion (pub_table/unidentifiable_barcode.txt): 1
P3-2D 

haploid (pub_table/haploid.txt): 2
P3-4H
CNTRL-1F

low barcode <0.0001 : 1

miss > 0.05:
P1-11F
P1-7B
P2-2D

# seperate bed file
```
awk '{
    # Use the value in the 4th column as the filename
    filename = "LOH_minSNP-5_" $4 ".bed"
    # Append the line to the corresponding file
    print $0 > filename
}' LOH_minSNP-5.bed

awk '{
    # Use the value in the 4th column as the filename
    filename = "revLOH_minSNP-5_" $4 ".bed"
    # Append the line to the corresponding file
    print $0 > filename
}' revLOH_minSNP-5.bed

```

# barseq - bartender
```bash
mkdir -p bartender
exp=( $(csvtk cut -f 1 -U pub_tables/barseq_fastq_path.csv) )
barpath=( $(csvtk cut -f 3 -U pub_tables/barseq_fastq_path.csv) )

for f in ${barpath[@]}; do
  sbatch -A tyjames0 --wrap="gzip -d -k $f"
done

barp='ATGTC[13-21]TAAAC'
for i in ${!exp[@]}; do
  inpath=${barpath[i]}
  outpath="bartender/${exp[i]}.extracted"
  sbatch -A tyjames0 --mem=8g --wrap="bartender_extractor_com -f ${inpath%\.gz} -o ${outpath} -q : -p $barp -m 2"
done

for i in ${!exp[@]}; do
  inpath="bartender/${exp[i]}.extracted_barcode.txt"
  outpath="bartender/${exp[i]}.clustered"
  sbatch -A tyjames0 --mem=8g --wrap="bartender_single_com -f ${inpath} -o ${outpath}"
done

python ./scripts/combine_barseq_exps.py
```


# QTL
```bash
bcftools annotate --set-id '%CHROM\_%POS' -Oz -o bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.named.vcf.gz bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.vcf.gz
plink --recode --vcf bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.named.vcf.gz --out bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.named.vcf --allow-extra-chr
plink --vcf bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.named.vcf.gz --recode A --out bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.named.vcf --allow-extra-chr

```

# LOH detect
```bash
python scripts/LOH_detect.py
```

# coverage
.bed file from meta_data.R

```bash
bedtools genomecov -bg -i <(bedtools sort -i tables/SK1_bedgraph_input.bed -faidx ref/S288C.chr.fasta.gz.fai) -g <(awk -v OFS='\t' '{print $1, $2}' ref/S288C.chr.fasta.gz.fai) > LOH_detect/forward.depth.bedgraph

bedtools genomecov -bg -i <(bedtools sort -i tables/BY4741_bedgraph_input.bed -faidx ref/S288C.chr.fasta.gz.fai) -g <(awk -v OFS='\t' '{print $1, $2}' ref/S288C.chr.fasta.gz.fai) > LOH_detect/reverse.depth.bedgraph

bedtools genomecov -bg -i <(bedtools sort -i <(cat tables/SK1_bedgraph_input.bed tables/BY4741_bedgraph_input.bed) -faidx ref/S288C.chr.fasta.gz.fai) -g <(awk -v OFS='\t' '{print $1, $2}' ref/S288C.chr.fasta.gz.fai) > LOH_detect/any.depth.bedgraph
```

# top-5 overlap
```bash
for i in "Acidic" "Caffeine" "H2O2" "SC" "YPD-30C" "YPD-37C" "worm-20C"; do
  bedtools genomecov -bg -i <(bedtools sort -i <(cat top5_list/$i/LOH*)) -g <(awk -v OFS='\t' '{print $1, $2}' ref/S288C.chr.fasta.gz.fai) > top5_list/LOH.$i.bedgraph
  done

for i in "Acidic" "Caffeine" "H2O2" "SC" "YPD-30C" "YPD-37C" "worm-20C"; do
  bedtools genomecov -bg -i <(bedtools sort -i <(cat top5_list/$i/revLOH*)) -g <(awk -v OFS='\t' '{print $1, $2}' ref/S288C.chr.fasta.gz.fai) > top5_list/revLOH.$i.bedgraph
  done
```

file to archive


ref/repeat_region.bed
ref/centromere.bed
ref/repeat_region.intervals
/nfs/turbo/lsa-tyjames/mycology/next_gen_seq_data/Yeast_LOH/20240722_LOH_wholegenome1/11314-YK/ .
/nfs/turbo/lsa-tyjames/mycology/next_gen_seq_data/Yeast_LOH/20240821_LOH_wholegenome2/11790-YK/ .
tables/project_id.csv


index ref/S288C.chr.fasta.gz
bwa_mapping/${names[$i]}.genome.info.bam

bwa_haplotypecaller_gvcf/runs.raw.vcf.gz
pub_tables/barcode_map.csv
tables/barcode_map_parsed.csv

ref/parents.raw.vcf.gz

bwa_haplotypecaller_finalvcf/parents.raw.vcf.gz*
bwa_haplotypecaller_finalvcf/parents.DPfilter.selected.vcf.gz*

bwa_haplotypecaller_finalvcf

pub_tables/barseq_fastq_path.csv


./scripts/combine_barseq_exps.py

bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.named.vcf.gz bwa_haplotypecaller_gvcf/runs.GTfilter.selected.vcf.gz

bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.named.vcf.gz bwa_haplotypecaller_finalvcf/runs.GTfilter.selected.vcf.gz

tables/SK1_bedgraph_input.bed 

LOH_detect/forward.depth.bedgraph

LOH_detect/

top5_list

ref/

bwa_mapping_depth/

ref/parents_rsid.list 

tables/genomes.final.vcf.gz het hom miss

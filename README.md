# Variant-calling
```
wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR333/ERR3335404/P7741_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR333/ERR3335404/P7741_R2.fastq.gz
```

```
sickle pe -f P7741_R1.fastq.gz -r P7741_R2.fastq.gz -t sanger -q 20-l 20 -g -o R1.trim.fastq.gz -p R2.trim.fastq.gz -s S_trim.fastq.gz
```

```
bwa mem ref/Ayg99.fna R1.trim.fastq.gz R2.trim.fastq.gz > output.sam
```
```
samtools view -S -b output.sam > output.bam
```
```
samtools sort -o output.sorted.bam output.bam
```
```
bwa mem -t 8 ref/Ayg99.fna R1.trim.fastq.gz R2.trim.fastq.gz | samtools sort -o outpud.sorted.bam -
```

```
samtools flagstats outpud.sorted.bam 
```

```
samtools flagstats outpud.sorted.bam > map-stats.txt
```
```
bcftools mpileup -O b -o raw.bcf -f ref/Ayg99.fna --threads 8 -q 20 -Q 30 outpud.sorted.bam
```
```
bcftools call --ploidy 1 -m -v -o variants.raw.vcf raw.bcf
```
#variyant miktarı
```
grep -v -c '^#' variants.raw.vcf 
```
32953

```
bcftools view -v snps variants.raw.vcf | grep -v -c '^#'
```
31295

```
bcftools view -v indels variants.raw.vcf | grep -v -c '^#'
```
1658
#Variyant position ları txt ye kaydetme

```
bcftools query -f '%POS\n' variants.raw.vcf > pos.txt
```


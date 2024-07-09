# Variant-calling

## Download Data
```bash
wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR333/ERR3335404/P7741_R1.fastq.gz
wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR333/ERR3335404/P7741_R2.fastq.gz
```

To download a selected dataset from websites containing biological data (such as ENA, NCBI, etc.) using the 'wget' command.

````
ls
````
P7741_R1.fastq.gz P7741_R2.fastq.gz

To list the downloaded files using the 'ls' command.

## Extract and view
````
gunzip P7741_R1.fastq.gz
````
To extract a compressed file using 'gunzip'.

````
head -n 4 P7741_R1.fastq 
````
The head command with the -n parameter is used to view the first sequence data in the dataset.


    @M01941:8:000000000-BRBPM:1:1101:16459:1430 1:N:0:23
    
    TCCTCGAGCTCGGTGGGCTCGAGGATCCGCTGGGCCAGCGGCAACAGATGCGGATGGTGCTCCGCGAGGACGCTTCCCGCGCTGGCCGTTGCGGTCACCGCCCCGATCGTCAGCCGTCCCGAGGTGTAGGCGGGATGCCCCAGCAGAAGGCGAAGGATCTCACCGCCTGCATAGCCGCTGGCACCCGCTACGGCCACCCTGATTGCATCGGTCGCGTGGGCCATTTCGAGGATTTTGCATGGTTATGCAAT
    
    '+'
    
    1>AA@AAAAFA?E000AGG00A00GFGFEGE?CFHCCE??A//?EFHHHF??/>@AGEFHHFGGCCCCCG?EGGGFFCCCCGGCCGGCCCC?CGCCFGCCC@C??CCGHHGGHHGGGGB?.-?CEFGGEF@@@@FFFFBE-AA--9FFFF?@?@FFFFFFFFFB?@BEFBFFFFFF@@??BBFEF???@@BF@?@FFFFE-FFBF/FF?BBB@@<@;-?AA-BFFB---@@FFFF/FBFFFFF/FFFF/B


The second line contains the DNA sequence, and the fourth line contains characters representing quality scores.

## Quality control

````
fastqc *.fastq*
````
The use of between two asterisks('*.fastq *') allows for quality control of all fastq files.

You can view the generated HTML file to access detailed information.

## Trim and Filter


If adapter contamination is detected using FastQC, it is preferable to use 'trimmomatic' at this stage. If no adapters are found, using 'sickl' is more practical.
### With trimmomatic 
````
trimmomatic PE -threads 4 P7741_R1.fastq.gz P7741_R2.fastq.gz \
	> R1.trim.fastq.gz R1.untrim.fastq.gz \
	> R2.trim.fastq.gz R2.untrim.fastq.gz \
	> SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:TruSeq3-PE.fa:2:40:15 
````
### With sickle

````
sickle pe -f P7741_R1.fastq.gz -r P7741_R2.fastq.gz -t sanger -q 20-l 20 -g -o R1.trim.fastq.gz -p R2.trim.fastq.gz -s S_trim.fastq.gz
````

## Aligmnet
The reference genome can be downloaded using the previously downloaded dataset path. After downloading, to organize the reference genome more neatly, a directory named ref has been created and it has been moved there.

BWA is used to align short reads. For long reads, consider looking into different alternatives.

```
bwa mem ref/Ayg99.fna R1.trim.fastq.gz R2.trim.fastq.gz > output.sam
```
Forward and reverse reads are aligned to the reference genome, and output is obtained in SAM format.
```
samtools view -S -b output.sam > output.bam
```
SAM files are large in size, so compressing them into BAM format helps reduce file size, making subsequent operations faster and more efficient on the computer.

```
samtools sort -o output.sorted.bam output.bam
```
It is used to sort aligned data.


We can handle these three lines of code in the following single line.
```
bwa mem -t 8 ref/Ayg99.fna R1.trim.fastq.gz R2.trim.fastq.gz | samtools sort -o outpud.sorted.bam -
```
## Statistical summary
```
samtools flagstats outpud.sorted.bam 
```
The flagstat command provides a statistical summary of alignment data. This command is used to evaluate and analyze the quality of alignment data after alignment.In particular, it provides important information about alignment quality, read pairing status, and post-alignment filtering operations.

    567208 + 0 in total (QC-passed reads + QC-failed reads)
    	538788 + 0 primary
    	0 + 0 secondary
    	28420 + 0 supplementary
    	0 + 0 duplicates
    	0 + 0 primary duplicates
    	519453 + 0 mapped (91.58% : N/A)
    	491033 + 0 primary mapped (91.14% : N/A)
    	538788 + 0 paired in sequencing
    	269394 + 0 read1
    	269394 + 0 read2
    	459160 + 0 properly paired (85.22% : N/A)
    	487012 + 0 with itself and mate mapped
    	4021 + 0 singletons (0.75% : N/A)
    	700 + 0 with mate mapped to a different chr
    	144 + 0 with mate mapped to a different chr (mapQ>=5)


The percentage of reads correctly aligned to the reference is 91.58%.
There are no duplicates found.In case of duplicates, you can use the samtools 'rmdup' command to create a new BAM file.

#### Save the statistics as a text file
```
samtools flagstats outpud.sorted.bam > map-stats.txt
```


## Variant detection

```
bcftools mpileup -O b -o raw.bcf -f ref/Ayg99.fna --threads 8 -q 20 -Q 30 outpud.sorted.bam
```
It is a tool used to make variant calls from alignment files. It compares with the reference genome to identify variations, recording possible variants for each nucleotide. It produces output in BCF format.

```
bcftools call --ploidy 1 -m -v -o variants.raw.vcf raw.bcf
```
The 'call' command converts BCF format to VCF format. It skips positions that match the reference and only calls variant positions.


#### Total variants
```
grep -v -c '^#' variants.raw.vcf 
```
32953

#### Number of SNP's
```
bcftools view -v snps variants.raw.vcf | grep -v -c '^#'
```
31295

#### Number of indels
```
bcftools view -v indels variants.raw.vcf | grep -v -c '^#'
```
1658

```
#### Saving variant positions to a text file
bcftools query -f '%POS\n' variants.raw.vcf > pos.txt
```


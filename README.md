Installation for each software is included in the last section.

## Step 0: Set up your environment
Recomendation: set the reference sequences, output folders, and subfolders before starting
```
mkdir ref_seq
mkdir output_folder
mkdir fastqc_out fastp_out hisat2 consensus exonerate
```
## Step 0: Quality Control
This step checks the quality of sequencing using `fastqc`. 
```

```
## Step 1: Read Preprocessing
This step trims adapters, filters low-quality reads, and generates reports using `fastp v. 0.24.0`. 
```
./fastp -i your_folder/specie.read_1.fastq.gz \
  -I your_folder/specie.read_2.fastq.gz \
  -o output_folder/fastp_out/specie.R1.fq.gz \
  -O output_folder/fastp_out/specie.R2.fq.gz \
  --unpaired1 output_folder/fastp_out/specie.unpaired_R1.fastq.gz \
  --unpaired2 output_folder/fastp_out/specie.unpaired_R2.fastq.gz \
  -w 32 &
```
We can concatenate the unpaired reads to maximize the reads that match with the section we are interested. 
```
#remember to be on the same folder as your files
cat PORE.unpaired_R1.fastq.gz PORE.unpaired_R2.fastq.gz > PORE.up.fastq.gz
```

## Step 2: Alignment and mapping reads to a target sequence
The target sequence could be a chromosome or a specific region that your genes of interes are located. We need the `.fasta` and `.gff3` files. They can be downloaded directly from NCBI and saved into your designated folder (e.g. `ref_seqs`).

1. Index your target sequence (e.g. a chromosome), this step would generate a series of `.ht2` files to be used as index for `HISAT2`.
```
hisat2-build ref_chromosome.fa chr_index
```
2. Align the reads to `HISAT2`.
```
conda activate hisat2
# Paired reads command
hisat2 -p 24 -x chr_index -1 output_folder/fastp_out/specie.R1.fq.gz -2 output_folder/fastp_out/specie.R2.fq.gz -S output_folder/hisat2/specie.chr7.sam

# Unpaired reads command
hisat2 -p 24 -x chr_index -U output_folder/fastp_out/specie.up.fastq.gz -S output_folder/hisat2/specie.chr7.up.sam
conda deactivate
```
3. Mapping the reads to the reference sequence
```
#Locate yourself on output_folder/hisat2/
conda activate pomacea
# Paired reads
samtools view -bS specie.chr7.sam > specie.chr7.bam
rm specie.chr7.sam
samtools view -H specie.chr7.bam
samtools sort specie.chr7.bam -o specie.chr7_sorted.bam
samtools index specie.chr7_sorted.bam

# Unpaired reads
samtools view -bS specie.chr7.up.sam > specie.chr7.up.bam
rm specie.chr7.up.sam
samtools view -H specie.chr7.up.bam
samtools sort specie.chr7.up.bam -o specie.chr7_sorted.up.bam
samtools index specie.chr7_sorted.up.bam

# Merge paired and unpaired 

samtools merge specie.all.chr7_sorted.bam specie.chr7_sorted.bam specie.chr7_sorted.up.bam
rm specie.chr7.bam specie.chr7_sorted.bam specie.chr7.up.bam specie.chr7_sorted.up.bam 
rm specie.chr7_sorted.bam.bai specie.chr7_sorted.up.bam.bai

samtools index specie.all.chr7_sorted.bam

### Check quality for indexed files
samtools flagstat specie.all.chr7_sorted.bam --threads 16 -O tsv 

```
4. Check the coverage and depth

```
#Locate yourself at output_folder/hisat2
#conda activate pomacea

#Calculate the coverage 
bedtools genomecov -ibam specie.all.chr7_sorted.bam -bga > specie_coverage.tsv

#intersect the coverage with the (expected) position for your genes based on your reference specie `.bed` file containing only your target gene. 
bedtools intersect -a specie_coverage.tsv -b pomacea.pv2.gene.feature.bed -wa -wb > specie.coverage_per_gene.tsv
```
You can check `specie_coverage.tsv` and `specie.coverage_per_gene.tsv` with the `coverage.R` script provided in this repository. The results would suggest if they region where your gene of interest is located have good coverage or not. 

5. Obtain the reads with bcftools 
```
### Generar un archivo máscara con las regiones que no coberturan las variantes
awk '$4 == 0 {print $1"\t"$2"\t"$3}' specie_coverage.tsv > specie.mask.bed

### Run bcftools
bcftools mpileup -Ou -f /home/bioupg02/pomacea_pr/ref_seq/LG7.fasta specie.all.chr7_sorted.bam | bcftools call -mv -Oz -o specie.chr7.variants.vcf.gz
bcftools index specie.chr7.variants.vcf.gz

bcftools consensus -f /home/bioupg02/pomacea_pr/ref_seq/LG7.fasta specie.chr7.variants.vcf.gz -m specie.mask.bed > /home/bioupg02/pomacea_pr/pomacea_out/consensus/specie.LG7_consensus.fa
### Warning: this tool skip variants, in case it overlaps with another variant. Check if this happen in your position of interest.
### For specie it doesn't happen but a possible solution is to normalize the variant (reduces it to two records?) 
### bcftools norm -m -both -f ref_chromosome.fasta specie.chr7.variants.vcf.gz -Oz -o specie.chr7.variants.norm.vcf.gz
### bcftools index specie.chr7.variants.norm.vcf.gz
### Then repeat consensus generation, and later indexing

samtools faidx /home/bioupg02/pomacea_pr/pomacea_out/consensus/specie.chr_consensus.fa
conda deactivate
```
## Step3: Extract the CDS of target genes
```
conda activate exonerate
exonerate --model protein2genome ref_proteins.fa ../consensus/specie.chr_consensus.fa --showvulgar no --showalignment yes > specie.exo.out
exonerate --model protein2genome ref_proteins.fa ../consensus/specie.chr_consensus.fa --showvulgar no --showalignment yes > specie.exo.out
conda deactivate
```
We need to organise and filter the data in order it'easier to take a decision, you can use `exo2tables.py` and get a `.tsv` summary you can inspect in any text reader or excel. 
```
conda activate pandas
# It would ask the file you would like to run to get the information. 
python3 exo2tables.py
conda deactivate
```

## Step4: Analizing your results 


## Step5: Which further analysis can I do with this results?

You get to the end of the pipeline!



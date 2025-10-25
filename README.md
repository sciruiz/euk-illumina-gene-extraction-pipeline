
# Gene Extraction Pipeline from Illumina Reads

This pipeline describes how to extract and validate genes of interest from Illumina paired-end reads of Mollusks (but might be useful for other Eukarya as well!). 

## Tools used in this pipeline

| Step            | Tool(s) with version      | 
| --------------- | ----------------- | 
| QC              | FastQC            | 
| Trimming        | fastp             | 
| Alignment       | HISAT2, SAMtools  | 
| Coverage        | BEDtools          | 
| Variant Calling | BCFtools          | 
| Gene Extraction | Exonerate, python | 
---

## Set Up the Environment

**Recommendation:**
Before starting, create folders for references, outputs, and intermediate results.

```bash
mkdir ref_seq
mkdir output_folder
mkdir fastqc_out fastp_out hisat2 consensus exonerate
```

---

## Step 1: Quality Control (FastQC) and Read Preprocessing (fastp)

This step checks the raw read quality using **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**.

```bash
# Example command
fastqc reads/specie.read_1.fastq.gz reads/specie.read_2.fastq.gz -o output_folder/fastqc_out/
```

---

Then, we can trim adapters, filters low-quality reads, and generates quality reports using **[FastP v0.24.0](https://github.com/OpenGene/fastp)**. Alternatively you can use Trimmomatic or other software. 

```bash
./fastp -i reads/specie.read_1.fastq.gz \
  -I reads/specie.read_2.fastq.gz \
  -o output_folder/fastp_out/specie.R1.fq.gz \
  -O output_folder/fastp_out/specie.R2.fq.gz \
  --unpaired1 output_folder/fastp_out/specie.unpaired_R1.fastq.gz \
  --unpaired2 output_folder/fastp_out/specie.unpaired_R2.fastq.gz \
  -w 32 
```

To include unpaired reads in later analysis (maximize usable reads):

```bash
cat specie.unpaired_R1.fastq.gz specie.unpaired_R2.fastq.gz > specie.up.fastq.gz
```

---

## Step 2: Alignment and Mapping to Target Sequence (HISAT2)

This step aligns the cleaned reads to your target genome section or chromosome using **[HISAT2](https://daehwankimlab.github.io/hisat2/)**.

### 2.0 Genome Retrieval 

You can download the region or full chromosome from NCBI using `efetch`. 
1. Go to NCBI Datasets, and search the species of your interest.
2. Scroll down, and see if there is an available genome at NCBI. In case not, maybe you can search in alternative websites such as EMBL.
3. Go to BLAST the reference genome.
4. By default, it will show the `blastn`, add your query sequence (e.g. mRNA, CDS, partial rRNA, etc.) and then go to BLAST. In case you have a protein sequence, select the option `tblastn`, add your query sequence and then go to BLAST.
5. If your query sequence hits a match, then you can see the range (e.g. 234432111 to 234433221) in the accesion (e.g. NC_037596.1) that it's more accurate for your analysis.

Now you can download your information to the folder of your interest (e.g. ref_seq). 
```
 efetch -db nucleotide -id NC_037596.1 -format fasta > species.fasta
```
Also you can download a portion of the sequence

⚠️ Only if you are sure that the sequence is on that region. 
```
efetch -db nucleotide -id NC_037596.1 -seq_start 2350000 -seq_stop 2720000 -format fasta > species_region.fasta
```
    
### 2.1 Index the Reference Sequence

Prepare `.ht2` index files for HISAT2:

```bash
hisat2-build ref_chromosome.fa chr_index
```

### 2.2 Align Reads

```bash
conda activate hisat2

# Paired reads
hisat2 -p 24 -x chr_index \
  -1 output_folder/fastp_out/specie.R1.fq.gz \
  -2 output_folder/fastp_out/specie.R2.fq.gz \
  -S output_folder/hisat2/specie.chr7.sam

# Unpaired reads
hisat2 -p 24 -x chr_index \
  -U output_folder/fastp_out/specie.up.fastq.gz \
  -S output_folder/hisat2/specie.chr7.up.sam

conda deactivate
```

---

### 2.3 Convert and Sort Alignments

Use **[SAMtools](http://www.htslib.org/)** to convert `.sam` → `.bam`, sort, and index:

```bash
# Move to output_folder/hisat2/
conda activate gene_extract

# Paired reads
samtools view -bS specie.chr7.sam > specie.chr7.bam
rm specie.chr7.sam
samtools sort specie.chr7.bam -o specie.chr7_sorted.bam
samtools index specie.chr7_sorted.bam

# Unpaired reads
samtools view -bS specie.chr7.up.sam > specie.chr7.up.bam
rm specie.chr7.up.sam
samtools sort specie.chr7.up.bam -o specie.chr7_sorted.up.bam
samtools index specie.chr7_sorted.up.bam

# Merge paired + unpaired BAMs
samtools merge specie.all.chr7_sorted.bam specie.chr7_sorted.bam specie.chr7_sorted.up.bam
samtools index specie.all.chr7_sorted.bam

# Clean up temporary files
rm specie.chr7.bam specie.chr7_sorted.bam specie.chr7.up.bam specie.chr7_sorted.up.bam 
rm specie.chr7_sorted.bam.bai specie.chr7_sorted.up.bam.bai

# Check mapping quality
samtools flagstat specie.all.chr7_sorted.bam --threads 16 -O tsv
```

---

### 2.4 Check Coverage and Depth

Use **[BEDTools](https://bedtools.readthedocs.io/)** to evaluate coverage across your target genes.

```bash
# Extract target gene CDS positions
grep -Ff target.genes.txt ref_seq/chr.gff3 | awk '$3 == "CDS"' > chr.target_cds.gff3
```

```bash
# Calculate coverage
bedtools genomecov -ibam specie.all.chr7_sorted.bam -bga > specie_coverage.tsv

# Convert CDS regions from GFF3 to BED
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $9, ".", $7}' ref_seq/chr.target_cds.gff3 > ref_seq/chr.target_cds.bed

# Intersect coverage with target genes
bedtools intersect -a specie_coverage.tsv -b ref_seq/chr.target_cds.bed -wa -wb > specie.coverage_per_gene.tsv
```

Visualize coverage results with the provided **`coverage.R`** script.
These metrics help assess sequencing depth and coverage robustness.

---

### 2.5 Generate Consensus Sequence (BCFtools)

Use **[BCFtools](http://samtools.github.io/bcftools/)** to call variants and build consensus.

```bash
# Mask uncovered regions
awk '$4 == 0 {print $1"\t"$2"\t"$3}' specie_coverage.tsv > specie.mask.bed

# Variant calling and consensus
bcftools mpileup -Ou -f ref_seq/chr.fasta specie.all.chr7_sorted.bam | bcftools call -mv -Oz -o specie.chr7.variants.vcf.gz
bcftools index specie.chr7.variants.vcf.gz

bcftools consensus -f ref_seq/chr.fasta specie.chr7.variants.vcf.gz -m specie.mask.bed > output_folder/consensus/specie.chr_consensus.fa
samtools faidx output_folder/consensus/specie.chr_consensus.fa
conda deactivate
```

**Note:**
If overlapping variants occur, normalize before consensus:

```bash
bcftools norm -m -both -f ref_seq/chr.fasta specie.chr7.variants.vcf.gz -Oz -o specie.chr7.variants.norm.vcf.gz
bcftools index specie.chr7.variants.norm.vcf.gz
```

---

## Step 3: Extract CDS of Target Genes (Exonerate)

Use **[Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)** to align reference protein sequences to your consensus sequence.

```bash
conda activate exonerate
exonerate --model protein2genome ref_proteins.fa ../consensus/specie.chr_consensus.fa \
  --showvulgar no --showalignment yes > specie.exo.out
conda deactivate
```

Process Exonerate results with the provided `exo2tables.py` script to generate a readable summary (`.tsv`).

```bash
conda activate pandas
python3 exo2tables.py
conda deactivate
```

Check:

1. **Raw score** – higher is better
2. **Position** – matches expected genomic region

---

Once you identify your alignment of interest, extract the FASTA sequence simply pasting your `Exonerate` alignmentt with this [Google Colab Notebook](https://colab.research.google.com/drive/18Gp7TfJn1hbEwWwZNy6iVgCyaSjYvaKe?usp=sharing) Each will print or download your FASTA sequence (or file), ready for downstream analysis.

---

## Step 4: Validate Gene Identity

Paste your sequence into:

* [ORFfinder (NCBI)](https://www.ncbi.nlm.nih.gov/orffinder)
* [BLAST Suite (NCBI)](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

These tools help confirm your gene identity.

---

## What’s Next?

If your genes seems to be what you were expecting... Congratulations! You can proceed to do a phylogenetic analysis, continue characterising the structure of the protein or another related analysis.

---

# Appendix: Software Installation

It is highly recommended to use **Conda** (or **Mamba**, for faster installations**) to manage environments and avoid dependency issues.
You can install [Miniconda](https://docs.anaconda.com/miniconda/) or [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) if not already available.

---

## Create a Base Bioinformatics Environment

You can either install all tools in one environment or separate them as suggested*.
*Update 2025: There is a conflict between pandas and exonerate, the second option is recommended. 

```bash
conda create -n gene_extract

conda activate gene_extract

conda install bioconda::fastqc
conda install bioconda::fastp
conda install bioconda::samtools
conda install bioconda::bedtools
conda install bioconda::bcftools
conda install bioconda::entrez-direct

conda deactivate
```

---

## HISAT2 — Read Alignment Environment

```bash
conda create -n hisat2 
conda install bioconda::hisat2
conda deactivate
```

---

## Exonerate — Gene Alignment and Extraction Environment

```bash
conda create -n exonerate 
conda install bioconda::exonerate
conda deactivate
```

---

## Python (Pandas & Utilities) — Post-processing
** In case it is not installed, usually python suites include pandas. 
```bash
conda create -n pandas
conda install pandas
conda deactivate
```

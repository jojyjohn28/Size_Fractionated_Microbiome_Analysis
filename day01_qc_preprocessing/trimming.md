#### Trimming & Filtering — Trimmomatic and Cutadapt (Day 1)

Trimming removes adapter contamination and low-quality bases so reads are suitable for downstream taxonomic profiling and genome-resolved workflows.

Clean input = more accurate classification, mapping, and abundance estimates.

---

#### Option A — Trimmomatic (classic, widely used)

#### Install

```bash
conda install -c bioconda trimmomatic
conda activate trimmomatic
```

or you can load the required module using

```bash
module load trimmomatic/0.39
 # if it is laredy installed on your HPC
```

##### Example (paired-end)

```bash
trimmomatic PE -threads 16 \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  sample_R1.trimmed.fastq.gz sample_R1.unpaired.fastq.gz \
  sample_R2.trimmed.fastq.gz sample_R2.unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 \
  TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:50
```

**What these parameters do**

✦ ILLUMINACLIP: removes adapter contamination using an adapter FASTA

✦ LEADING / TRAILING: trims low-quality bases at read ends

✦ SLIDINGWINDOW: trims when average quality drops within a sliding window

✦ MINLEN: removes very short reads after trimming

#### Option B — Cutadapt (explicit adapter control)

Cutadapt is useful when you want:

✦ precise adapter sequence control

✦ poly-A/T trimming

✦ complex read structures or special libraries

#### Install

```bash
conda install -c bioconda cutadapt
conda activate cutadapt
```

or you can load the required module using

```bash
module load cutadapt/4.9
 # if it is laredy installed on your HPC
```

##### Example (paired-end)

```bash
cutadapt \
  -j 16 \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -q 15,15 \
  --minimum-length 50 \
  -o sample_R1.trimmed.fastq.gz \
  -p sample_R2.trimmed.fastq.gz \
  sample_R1.fastq.gz sample_R2.fastq.gz
```

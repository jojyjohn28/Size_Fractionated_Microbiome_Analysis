#### FastQC — Raw Read Quality Assessment (Day 1)

FastQC provides a quick visual QC report for raw shotgun metagenomic reads.  
This is the first checkpoint before trimming, profiling, or assembly.

## Why run FastQC?

Raw reads often contain:

- Adapter contamination
- Low-quality tails (especially toward the 3′ end)
- Overrepresented sequences
- GC bias or unusual distributions

Identifying these issues early prevents false positives and poor mapping in downstream tools (Kaiju, MetaPhlAn, mOTUs).

## Installation

Using conda:

```bash
conda install -c bioconda fastqc
```

#### Run FastQC (single sample)

```bash
fastqc sample_R1.fastq.gz sample_R2.fastq.gz
```

For batch/loop see fastqc_loop.sh in folder scripts/ and make it Make executable:

```bash
chmod +x scripts/fastqc_loop.sh
```

#### What to check in the report

📌 Per base sequence quality: quality typically drops at read ends

📌 Adapter content: indicates trimming is needed

📌 Overrepresented sequences: often adapters/primers/low complexity

📌 Sequence length distribution: should match expected read length

📌 Per sequence GC content: major shifts can indicate contamination or mixed libraries

✅ Raw reads can look messy — that’s normal.
**The goal is to confirm what needs trimming, not to expect perfect reads.**
➡️ Next: Trimming and filtering

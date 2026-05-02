# Dna-methylation-and-analysis
# 🧬 DNA Methylation Data Analysis
### Bisulfite Sequencing Pipeline using Galaxy — Following the GTN Tutorial


## 📌 Project Overview

This repository contains the **outputs, results, and documentation** from completing the Galaxy Training Network (GTN) hands-on tutorial on **DNA Methylation Data Analysis** using bisulfite sequencing data. The tutorial was carried out entirely within the [Galaxy](https://usegalaxy.eu) web-based bioinformatics platform — no local software installation required.

**Tutorial followed:**
> Wolff J., Ryan D., Moosmann V. (2017, updated 2025). *DNA Methylation data analysis.* Galaxy Training Network.
> https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html

**Original study this tutorial is based on:**
> Lin I.H. et al. (2015). *Hierarchical clustering of breast cancer methylomes revealed differentially methylated and expressed breast cancer genes.* PLOS ONE, 10(2), e0118453.

---

## 🧪 Biological Background

### What is DNA Methylation?

DNA methylation is an **epigenetic modification** — a heritable change in gene expression that does not alter the underlying DNA sequence. In mammals, it occurs almost exclusively at **CpG dinucleotides** (cytosine followed by guanine), where a methyl group (–CH₃) is added to the 5-carbon position of cytosine by DNA methyltransferase (DNMT) enzymes, producing **5-methylcytosine (5mC)**.

```
Normal C  →  C  (unmethylated)
Methylated C  →  ᵐC  (5-methylcytosine)
```

### Why Can't Normal NGS Detect Methylation?

Standard next-generation sequencing reads A, T, G, C — it cannot distinguish between a methylated and an unmethylated cytosine because both look like "C" in the raw sequence. **Bisulfite conversion** solves this problem:

```
Bisulfite treatment:
  Unmethylated C  →  converts to U  →  reads as T in sequencing
  Methylated C    →  remains as C   →  reads as C in sequencing
```

This asymmetry means that after bisulfite conversion and sequencing, **every C in the reads is a methylated cytosine** and **every T (that was a C in the reference) is an unmethylated cytosine**. Standard aligners cannot handle this because the genome is now effectively a 3-letter code — specialised bisulfite-aware aligners are required.

### Why Does It Matter?

Methylation at gene **promoters** typically **silences gene expression** by preventing transcription factor binding and recruiting repressive chromatin-remodelling complexes. Loss of methylation at normally-silenced regions (hypomethylation) or gain of methylation at normally-active promoters (hypermethylation) are hallmarks of cancer. In this tutorial, we compare:

| Sample | Type | Methylation Status Expected |
|---|---|---|
| **NB** | Normal breast cells | Baseline / normal methylation |
| **BT089** | Fibroadenoma (benign breast tumor) | Near-normal |
| **BT126** | Invasive ductal carcinoma | Aberrant methylation |
| **BT198** | Invasive ductal carcinoma | Aberrant methylation |
| **MCF7** | Breast adenocarcinoma cell line | Strong aberrant methylation |

---

## 🗂️ Repository Structure

```
Dna-methylation-and-analysis/
│
├── 📁 falco/
│   └── Quality control reports from Falco (FastQC-compatible)
│       for subset_1.fastq and subset_2.fastq
│
├── 📁 bigwig/
│   └── BigWig coverage/methylation tracks for genome browser
│       visualisation (output of MethylDackel)
│
├── 📁 metilene_qval_plots/
│   └── Q-value distribution plots from Metilene DMR calling
│       (differentially methylated region significance plots)
│
└── README.md
```

---

## 🔬 Analysis Pipeline

The complete pipeline follows 6 steps as defined in the GTN tutorial:

```
Raw FASTQ reads (subset_1.fastq, subset_2.fastq)
        │
        ▼
┌─────────────────────┐
│  Step 1: Data Upload │  ← Zenodo dataset via Galaxy upload
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Step 2: Quality     │  ← Falco (FastQC-compatible)
│  Control             │
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Step 3: Alignment   │  ← bwameth → hg38 reference genome
│  (Bisulfite-aware)   │
└─────────────────────┘
        │
        ▼
┌──────────────────────────────────┐
│  Step 4: Methylation Bias        │  ← MethylDackel (mbias mode)
│  & CpG Metric Extraction         │     + bedGraph output
└──────────────────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Step 5: Visualisa-  │  ← bamCoverage → BigWig
│  tion                │     → IGV / UCSC genome browser
└─────────────────────┘
        │
        ▼
┌─────────────────────┐
│  Step 6: DMR Calling │  ← Metilene (differentially methylated
│                      │     regions between NB vs tumour)
└─────────────────────┘
```

---

## 🛠️ Tools Used

### Step 1 — Data Upload

| Detail | Value |
|---|---|
| **Input data** | `subset_1.fastq`, `subset_2.fastq` |
| **Source** | Zenodo record 557099 (https://zenodo.org/record/557099) |
| **Format** | FASTQ (paired-end bisulfite sequencing reads) |
| **Samples** | Subset of NB, BT089, BT126, BT198, MCF7 breast tissue samples |

The data is a curated subset of the Lin et al. 2015 methylation-seq dataset, downsampled to make the tutorial computationally feasible. Full data is available at Zenodo.

---

### Step 2 — Quality Control: Falco

> **Tool:** Falco v1.2.4 (Galaxy version: `falco/1.2.4+galaxy0`)
> **What it is:** An efficiency-optimised reimplementation of FastQC for rapid quality assessment of raw sequencing reads.

**Parameters used:**
- Input: `subset_1.fastq.gz` and `subset_2.fastq.gz` (both reads in parallel)

**Outputs stored in** `falco/`:
- Per-base sequence quality scores
- Per-base sequence content (GC distribution)
- Sequence duplication levels
- Adapter content

**Key observation — why bisulfite QC looks "weird":**

In bisulfite-converted reads, the **cytosine content drops dramatically** and **thymine content spikes** because all unmethylated cytosines are converted to uracil (→ thymine). This means:
- The C/T ratio is completely different from normal genomic sequencing
- The "Per base sequence content" plot will show an unusual pattern with very low %C and high %T
- This is **expected and correct** — it is NOT a quality problem

> The attentive biologist knows: every methylated C stays as C, and every unmethylated C becomes T during bisulfite conversion. Standard FastQC/Falco interpretation rules do not apply directly to bisulfite data.

**Results folder:** `falco/`

---

### Step 3 — Alignment: bwameth

> **Tool:** bwameth v0.2.7 (Galaxy version: `bwameth/0.2.7+galaxy0`)
> **What it is:** A bisulfite-aware short-read aligner that converts the genome reference in silico to handle C→T converted reads without losing alignment accuracy.

**Why a specialised aligner is needed:**

Standard aligners (BWA, Bowtie2, STAR) assume that C in a read matches C in the reference. In bisulfite data, an unmethylated C in the read has been converted to T. A standard aligner would therefore miss or mismatch these reads. bwameth handles this by:
1. Creating a C→T and G→A converted version of the reference genome
2. Mapping converted reads to the converted reference
3. Inferring methylation status from mismatches at CpG positions

**Parameters used:**
| Parameter | Value |
|---|---|
| Reference genome | Human hg38 (built-in, full) |
| Library type | Paired-end |
| Read 1 | `subset_1.fastq` |
| Read 2 | `subset_2.fastq` |

**Output:** Aligned BAM file (`aligned_subset.bam`)

> ⚠️ Note: Alignment is computationally intensive. The tutorial provides a pre-computed BAM file at `https://zenodo.org/records/557099/files/aligned_subset.bam` for users who want to skip this step.

---

### Step 4 — Methylation Bias & Metric Extraction: MethylDackel

> **Tool:** MethylDackel v0.5.2 (Galaxy version: `methyldackel/0.5.2+galaxy0`)
> **What it is:** A tool for extracting per-CpG methylation metrics from bisulfite-aligned BAM files and diagnosing methylation bias.

#### Part A — Methylation Bias Analysis (mbias mode)

**Purpose:** Detect whether methylation levels are systematically higher or lower at the 5′ or 3′ ends of reads — a common artefact of bisulfite library preparation. If bias is detected, the affected positions must be ignored (trimmed) during extraction.

**Parameters:**
| Parameter | Value |
|---|---|
| Reference genome | Human hg38 (local cache) |
| Input BAM | Output of bwameth |
| Mode | `mbias` — produces diagnostic SVG plots |
| Keep singletons | Yes |
| Keep discordant alignments | Yes |

**Outputs:** SVG plots showing methylation percentage along read position (OB and OT strands). A flat line indicates no bias; a curved line at ends indicates trimming is needed.

#### Part B — CpG Methylation Extraction (bedGraph output)

**Purpose:** Extract per-CpG methylation levels across the entire genome into a bedGraph format for downstream analysis and visualisation.

**Parameters:**
| Parameter | Value |
|---|---|
| Mode | Extract methylation metrics (bedGraph) |
| OT/OB strand bias trimming | Applied based on mbias results |

**Output:** bedGraph file with columns: chromosome, start, end, methylation percentage

---

### Step 5 — Visualisation: BigWig + Genome Browser

> **Tool:** bamCoverage (deepTools) + IGV / UCSC Genome Browser

**Purpose:** Convert alignment and methylation data into BigWig format for interactive genome browser visualisation, enabling visual inspection of methylation patterns at specific genomic loci.

**Process:**
1. The aligned BAM and methylation bedGraph are converted to **BigWig** format using `bamCoverage`
2. BigWig files are loaded into a genome browser (IGV or UCSC) to visually inspect methylation across genomic regions

**Outputs stored in** `bigwig/`:
- BigWig tracks representing read coverage and/or methylation signal across the hg38 genome

**What BigWig files show:**
- High signal at a CpG = high methylation at that position
- Low/no signal = unmethylated
- Comparison across samples (NB vs tumour) reveals visually which regions gain or lose methylation in cancer

---

### Step 6 — Differentially Methylated Regions: Metilene

> **Tool:** Metilene (Galaxy version available)
> **What it is:** A fast and sensitive method for calling differentially methylated regions (DMRs) between two conditions using a two-dimensional segmentation algorithm.

**Purpose:** Identify genomic regions where methylation levels differ significantly between normal breast tissue (NB) and tumour samples (BT089, BT126, BT198, MCF7).

**How Metilene works:**
Metilene uses a two-dimensional Kolmogorov-Smirnov test combined with a mean test to identify consecutive CpG sites that show consistent, statistically significant methylation differences between two groups. It is optimised for whole-genome bisulfite sequencing data.

**Inputs:**
- Precomputed bedGraph files from the full (non-subset) dataset for all 5 samples
- Group A: Normal breast (NB)
- Group B: Tumour samples (BT089, BT126, BT198, MCF7)

**Key output columns in Metilene results:**

| Column | Meaning |
|---|---|
| chr | Chromosome |
| start / end | DMR coordinates on hg38 |
| q-value | Multiple-testing corrected significance (Benjamini-Hochberg) |
| mean methylation diff | Average methylation difference (Group B − Group A) |
| # CpGs | Number of CpG sites in the DMR |
| mean Group A | Mean methylation in normal tissue |
| mean Group B | Mean methylation in tumour tissue |

**Interpretation:**
- **Positive mean diff** → hypermethylated in tumour vs. normal (potentially silenced tumour suppressor)
- **Negative mean diff** → hypomethylated in tumour vs. normal (potentially activated oncogene)

**Outputs stored in** `metilene_qval_plots/`:
- Q-value distribution plots showing the significance distribution of all called DMRs
- These plots allow assessment of how many DMRs pass various significance thresholds

---

## 📊 Key Biological Findings

### What the Results Show

**From QC (Falco):**
- Both FASTQ files passed quality metrics
- As expected for bisulfite data, per-base %C is very low and %T is abnormally high — this is correct, not a problem
- Read quality scores are high throughout, confirming the sequencing run was successful

**From Alignment (bwameth):**
- Reads align successfully to hg38 using bisulfite-aware mapping
- The methylation conversion (C→T) is handled transparently by bwameth

**From Methylation Bias (MethylDackel):**
- Methylation bias plots reveal whether library preparation introduced end-effects
- Any detected bias at read ends is corrected by trimming those positions in the extraction step

**From Visualisation (BigWig):**
- Methylation tracks loaded in a genome browser reveal regions of high and low methylation
- Comparison between normal (NB) and tumour samples reveals cancer-associated methylation changes at specific loci

**From DMR Calling (Metilene):**
- Metilene identifies statistically significant DMRs between normal breast and tumour samples
- Q-value plots (in `metilene_qval_plots/`) show the distribution of significance values, with many DMRs passing q < 0.05 — confirming widespread epigenetic dysregulation in breast cancer
- Results replicate the biological finding of Lin et al. (2015): breast cancer samples show pervasive aberrant methylation patterns compared to normal tissue

---

## 🚀 How to Reproduce This Analysis

### Prerequisites
- A Galaxy account (free) — use https://usegalaxy.eu or https://usegalaxy.org
- No local software installation required — everything runs in the browser

### Steps

1. **Go to Galaxy:** https://usegalaxy.eu
2. **Create a new history** (click the ✚ icon in the History panel)
3. **Upload input data** from Zenodo:
   ```
   https://zenodo.org/record/557099/files/subset_1.fastq
   https://zenodo.org/record/557099/files/subset_2.fastq
   ```
4. **Run Falco** on both FASTQ files → inspect QC reports
5. **Run bwameth** with hg38, paired-end → get aligned BAM
6. **Run MethylDackel** in mbias mode → inspect bias SVGs
7. **Run MethylDackel** in extract mode → get bedGraph
8. **Convert to BigWig** using bamCoverage → visualise in IGV
9. **Run Metilene** on precomputed multi-sample bedGraphs → get DMR table
10. Compare results against the outputs stored in this repository

**Full step-by-step tutorial:** https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html

**Pre-computed answer history (Galaxy EU):**
https://usegalaxy.eu/u/videmp/h/dna-methylation-data-analysis-gtn-answer

---

## 📁 Folder Contents

### `falco/`
Contains HTML and text quality control reports generated by Falco for `subset_1.fastq` and `subset_2.fastq`. Open the `.html` files in a browser to view the interactive FastQC-style QC report. Key plots to examine:
- **Per base sequence content** — will show unusual C/T pattern due to bisulfite conversion
- **Per base quality scores** — should be high (green zone) throughout
- **Sequence duplication levels** — moderate duplication is normal in bisulfite data

### `bigwig/`
Contains BigWig files representing methylation/coverage signal across the hg38 human genome. These files can be loaded into:
- **IGV (Integrative Genomics Viewer):** File → Load from File
- **UCSC Genome Browser:** My Data → Custom Tracks
- **deepTools computeMatrix / plotProfile** for aggregate methylation plots

BigWig format: chromosome, start, end, signal value (methylation % or coverage depth).

### `metilene_qval_plots/`
Contains Q-value distribution plots output by Metilene after DMR calling. These plots show:
- **X-axis:** Q-value threshold (0.0 → 1.0)
- **Y-axis:** Number of DMRs passing that threshold
- A steep drop-off near q = 0.05 indicates many statistically significant DMRs
- These plots are used to choose an appropriate significance cutoff for downstream analysis

---

## 🧰 Tools Reference

| Tool | Version | Purpose | Documentation |
|---|---|---|---|
| **Falco** | 1.2.4 | Quality control of raw FASTQ reads | https://falco.readthedocs.io |
| **bwameth** | 0.2.7 | Bisulfite-aware read alignment to hg38 | https://github.com/brentp/bwa-meth |
| **MethylDackel** | 0.5.2 | Methylation bias + CpG metric extraction | https://github.com/dpryan79/MethylDackel |
| **bamCoverage** | deepTools | BAM → BigWig conversion for visualisation | https://deeptools.readthedocs.io |
| **Metilene** | — | Differentially methylated region (DMR) calling | https://www.bioinf.uni-leipzig.de/Software/metilene/ |
| **Galaxy** | — | Cloud bioinformatics platform (no install needed) | https://usegalaxy.eu |

---

## 📖 Glossary

| Term | Definition |
|---|---|
| **CpG** | Cytosine–phosphate–Guanine dinucleotide — the primary site of DNA methylation in mammals |
| **5mC** | 5-methylcytosine — methylated cytosine base |
| **Bisulfite conversion** | Chemical treatment that converts unmethylated C → U (→T) while leaving methylated C unchanged |
| **Bisulfite sequencing (BS-Seq)** | Sequencing after bisulfite conversion to detect methylation at single-CpG resolution |
| **bwameth** | Bisulfite-aware aligner based on BWA-MEM |
| **MethylDackel** | Tool for extracting CpG methylation from bisulfite BAM files |
| **bedGraph** | Tab-delimited file format: chr, start, end, value — used to represent methylation per CpG |
| **BigWig** | Binary compressed format for genome-wide signal tracks — used in genome browsers |
| **DMR** | Differentially Methylated Region — genomic region with statistically different methylation between conditions |
| **Metilene** | Tool for calling DMRs between two groups using 2D segmentation |
| **Hypermethylation** | Higher methylation in tumour vs. normal — often silences tumour suppressors |
| **Hypomethylation** | Lower methylation in tumour vs. normal — often activates oncogenes |
| **mbias** | MethylDackel diagnostic mode: detects position-dependent methylation artefacts along reads |
| **Q-value** | False discovery rate (FDR)-corrected p-value — used in Metilene DMR significance |

---



*DNA Methylation Data Analysis · Bisulfite Sequencing · Galaxy Platform · GTN Tutorial*

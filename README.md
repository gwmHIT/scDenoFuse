# scDenoFuse

Gene fusion detection from full-length, paired-end single-cell RNA-seq data.

## Installation

Use Linux or WSL with Conda. Download or clone this repository and enter its root directory:

```bash
cd scDenoFuse
conda create -n scdenofuse --override-channels \
  -c conda-forge -c bioconda --strict-channel-priority \
  python=3.9 star=2.7.11b samtools bedtools r-base numpy pysam pyensembl
conda activate scdenofuse
python scDenoFuse.py --help
```

All internal modules are included. No separate scFusion installation or additional third-party R packages are required. `requirements.txt` lists the Python dependencies; the Conda command also installs the external tools.

## Usage

Use matching **GRCh38** genome FASTA, gene annotation GTF and mappability data. Prepare the two reference directories once:

```bash
python scDenoFuse.py BuildSTARIndex -g genome.fa -a annotation.gtf -s STAR_index -t 16
python scDenoFuse.py Index -g genome.fa -a annotation.gtf -d reference
```

The mappability input must be an uncompressed, tab-separated file with four columns: chromosome, start, end and score. Supply this file separately for the same reference assembly; bigWig input is not supported directly.

Name paired FASTQs `1_1.fastq.gz`, `1_2.fastq.gz`, `2_1.fastq.gz`, `2_2.fastq.gz`, etc., using consecutive positive cell numbers. Uncompressed `.fastq` files are also supported.

Run the complete pipeline:

```bash
python scDenoFuse.py -i fastq/ -s STAR_index/ -d reference/ \
  -m mappability.txt -o output/ -t 16
```

`Run` may be specified explicitly. Threads default to 16. Use a new empty output directory or an output directory created by this version. After successful result validation, the full pipeline removes `work/STARMapping/`; add `--keep-star` to retain BAM files for IGV or reuse. Failed runs and runs without a valid statistical threshold retain the mapping files.

Alternatively, run individual steps in order:

```bash
python scDenoFuse.py ReadMapping -i fastq/ -s STAR_index/ -o output/ -t 16
python scDenoFuse.py ReadProcessing -d reference/ -m mappability.txt -o output/ -t 16
python scDenoFuse.py FusionCandidate -d reference/ -o output/ -t 16
python scDenoFuse.py FusionReport -d reference/ -o output/ -t 16
python scDenoFuse.py Clean -o output/
```

`FusionReport` fits the model, applies filtering and exports results. To reapply filtering without refitting:

```bash
python scDenoFuse.py Filter -d reference/ -o output/
```

`Clean` removes only `work/STARMapping/`, preserving external symlink targets. Candidate generation, reporting and filtering remain available after cleanup; repeating ReadProcessing requires the mapping files.

## Output

```text
output/
├── Result.abridged.txt    # Final fusion calls
├── Result.full.txt        # Calls with per-cell evidence
├── work/                  # Intermediate files
└── logs/                  # Pipeline log, thresholds and run records
```

Both result files are tab-separated, with one row per fusion candidate and no fusion identifier column. Different breakpoints for the same gene pair remain separate rows.

| Column | Description |
|---|---|
| `FusionName` | Partner gene names joined by `--` |
| `Position1`, `Position2` | Corresponding breakpoints as `chromosome:position`, using 1-based coordinates |
| `strands` | Annotated gene strands, e.g. `+/-` |
| `directions` | Retained transcript portions: `u` = upstream/5′, `d` = downstream/3′; unknown values are `NA` |
| `CellSupport` | Number of cells with split-read support |
| `JunctionReadCount` | Total split-read support under the pipeline's deduplication rules |
| `SpanningFragCount` | Total spanning-fragment support across all cells |
| `Pvalue` | Model P-value, not FDR |

For a known `d/u` connection, partners, positions and strands are swapped together to display `u/d`.

`Result.full.txt` additionally contains:

| Column | Description |
|---|---|
| `SupportingCells` | Cell identifiers with split-read support |
| `SplitReadCountsByCell` | Split-read counts corresponding to `SupportingCells` |
| `SpanningCells` | Cell identifiers with spanning support, including cells without split-read support |
| `SpanningReadCountsByCell` | Spanning-fragment counts corresponding to `SpanningCells` |

Lists use JSON array notation, e.g. `[1, 5]` and `[3, 7]`. Each count list follows the order of its corresponding cell list.

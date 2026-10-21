<h1 align="center">vcf2curator</h1>

<p align="center">
  <a href="https://github.com/ClinSeq/vcf2curator/actions/workflows/tests.yml">
    <img src="https://github.com/ClinSeq/vcf2curator/actions/workflows/tests.yml/badge.svg" alt="Tests">
  </a>
</p>

A Python CLI that transforms raw bioinformatics pipeline output into clean,
schema-validated tables ready for ingestion into the **Curator** web app.

It consolidates a set of standalone scripts into one tool with a consistent
interface: a small set of plain subcommands with explicit arguments, a shared
Curator output writer, and consistent logging and error handling.

## Requirements

This project uses [`uv`](https://docs.astral.sh/uv/) exclusively for environment
and dependency management (no `pip`/`conda`/`poetry`). Python `>=3.11`.

```bash
# Install uv (once)
curl -LsSf https://astral.sh/uv/install.sh | sh
```

## Install

### As a standalone CLI (recommended)

Install `vcf2curator` as a real executable on your `PATH`, so it runs like any
other command — no `uv run` prefix:

```bash
# From the project root
uv tool install .

# ...or install editable, so local source edits are picked up without reinstalling
uv tool install --editable .
```

This places the executable at `~/.local/bin/vcf2curator`. Make sure
`~/.local/bin` is on your `PATH` (the uv installer adds it to your shell profile;
open a new shell or `source ~/.profile` if the command isn't found yet):

```bash
vcf2curator --help
```

Update or remove the installed tool later with:

```bash
uv tool install --editable . --reinstall   # refresh after pyproject/deps change
uv tool uninstall vcf2curator              # remove
```

### Run from the project environment (no install)

For a reproducible run pinned to `uv.lock` (e.g. in CI), skip the install and use
`uv run`, which executes inside the project's managed environment:

```bash
uv sync                        # create/sync .venv from uv.lock
uv run vcf2curator --help
```

## Usage

Examples below use the standalone `vcf2curator` command; prefix any of them with
`uv run` to run from the project environment instead.

```bash
# Convert a VEP-annotated somatic VCF into a Curator table
vcf2curator snvs \
    --input  sample.somatic.vep.vcf \
    --output sample.curator.tsv \
    --vcftype somatic \
    --oncokb  allAnnotatedVariants.txt \
    --cgc     cancer_gene_census.tsv \
    --hotspot-snv   hotspots_snv.csv \
    --hotspot-indel hotspots_indel.csv \
    --curation-ann  curation.csv

# Germline call set, whole-genome filters
vcf2curator snvs -i sample.germline.vep.vcf -o out.tsv \
    --vcftype germline --oncokb allAnnotatedVariants.txt --wgs

# Convert a GRIDSS SV VCF into a Curator table (panel data), with an evidence BAM
vcf2curator svs \
    --input  sample.gridss.vcf \
    --output sample.sv.curator.tsv \
    --vcftype somatic \
    --sdid P-0001-T-01 --project PROJ \
    --genes genes.bed --exons exons.bed --cgc cancer_gene_census.tsv \
    --targets capture_targets.bed \
    --curation-ann curation.csv \
    --bam tumor.bam --evidence-bam sample.sv.evidence.bam

# Combine one case's QC metrics (tumor + normal, panel) into a QC table
vcf2curator qc -o P-0001.qc.txt \
    --tumor-name  P-0001-T-01 \
    --tumor-hsmetrics  T.picard-hsmetrics.txt --tumor-markdup  T.picard-markdup.metrics.txt \
    --tumor-insertsize T.picard-insertsize.txt --tumor-contest T.contest.txt \
    --tumor-msings     T.MSI_Analysis.txt \
    --normal-name P-0001-N-01 \
    --normal-hsmetrics  N.picard-hsmetrics.txt --normal-markdup N.picard-markdup.metrics.txt \
    --normal-insertsize N.picard-insertsize.txt --normal-contest N.contest.txt

# Check a VCF is VEP-annotated with the required CSQ fields, without writing output
vcf2curator validate --input sample.vep.vcf

# See what would be written without touching disk
vcf2curator snvs ... --dry-run
```

Global flags: `--verbose/-v`, `--quiet/-q`. Distinct exit codes are used for
failure modes (input validation `3`, missing annotation resource `4`, schema
error `5`).

### `snvs` — SNV/indel VCF → Curator table

Reads a VEP-annotated somatic or germline VCF and emits the Curator table.
Annotation resources are supplied as arguments (`--oncokb`, `--cgc`,
`--hotspot-snv`, `--hotspot-indel`, `--curation-ann`); only `--oncokb` is
required. Filtering thresholds are configurable and default to the legacy
behaviour:

| Flag | Default | Effect |
|---|---|---|
| `--impact` (repeatable) | `HIGH`, `MODERATE` | IMPACT values treated as actionable |
| `--keep-filter` (repeatable) | somatic `PASS`,`LowQual` / germline `PASS` | FILTER values kept |
| `--include-splice / --no-include-splice` | on | keep splice-region variants |
| `--include-hotspot / --no-include-hotspot` | on | keep hotspot variants (somatic) |
| `--wgs` | off | apply whole-genome filters |
| `--wgs-min-alt` | `5` | min tumor ALT reads (somatic WGS) |

### `svs` — GRIDSS structural variant VCF → Curator table

Reads a GRIDSS VCF and writes one row per event; a translocation (TRA) gets one
row per breakend. Only PASS records are used. Single breakends and events
touching a contig outside 1–22/X/Y (with or without `chr`) are skipped.

- **SV type** (DEL/DUP/INS/INV/TRA) is derived from each record's breakend ALT,
  using the same rules as StructuralVariantAnnotation's `simpleEventType`, so no
  R step is needed.
- **Annotations:** `--genes` (BED with gene symbol in column 4) gives
  GENE_A/GENE_B, `--cgc` gives CGC_ANN, `--curation-ann` gives CURATE (matched
  on `--project` and `--vcftype`). `--exons` is a BED file.
- **CURATOR:** with `--targets` (panel BED), an event is YES when a breakend is
  within `--target-padding` bp of a target; with `--wgs`, when either gene is in
  CGC. Germline events below `--germline-min-support` reads, and non-TRA events
  that overlap no exon, are set to NO. Give exactly one of `--targets`/`--wgs`.
- **`--bam`** (optional) drops events whose supporting reads (INFO/BPNAMES) all
  start at the same place (SAME_START_READS). Add `--evidence-bam` to also write
  those reads to a sorted, indexed BAM for IGV review.

| Flag | Default | Effect |
|---|---|---|
| `--target-padding` | `150` | bp around a target that still counts as on target |
| `--germline-min-support` | `40` | germline: min supporting reads for CURATOR=YES |
| `--same-start-fraction` | `0.6` | SAME_START_READS: max share of reads at one position |
| `--min-distinct-reads` | `3` | SAME_START_READS: min distinct read positions |

### `qc` — QC metric files → QC overview table

Combines one case's QC files into a tab-separated table: a header, the tumor
row and, if given, the normal row. Each sample gets its own `--tumor-*` /
`--normal-*` options. Sample names come only from `--tumor-name` /
`--normal-name`, never from file names. Normal options are all-or-nothing;
leave them out for tumor-only cases.

| Design | Metrics option | Also required | Optional | Output columns |
|---|---|---|---|---|
| Panel | `--*-hsmetrics` | `--*-markdup`, `--*-insertsize`, `--*-contest` | `--*-msings` | SAMP, MEAN_TARGET_COVERAGE, FOLD_ENRICHMENT, dedupped_on_bait_rate, FOLD_80_BASE_PENALTY, READ_PAIRS_EXAMINED, PERCENT_DUPLICATION, contamination_%, MEDIAN_INSERT_SIZE, msing_score |
| WGS | `--*-wgsmetrics` | `--*-insertsize`, `--*-contest` | — | SAMP, MEAN_COVERAGE, FOLD_80_BASE_PENALTY, contamination_%, MEDIAN_INSERT_SIZE |

The design follows the metrics option given, and tumor and normal must use the
same one. Values are copied as written in the input files;
`dedupped_on_bait_rate` is `ON_BAIT_BASES / PF_BASES_ALIGNED`. Missing or
uncomputed values (Picard `?`, no mSINGS file) are written as `NA`. Only the
first data row of each Picard file is used.

## Development

Run these before committing:

```bash
uv run ruff check .
uv run ruff format --check .
uv run mypy
uv run pytest
```

## Project layout

```
vcf2curator/
├── cli.py               # builds the Typer app, registers commands
├── commands/            # one module per subcommand (snvs.py, svs.py, qc.py, validate.py, _shared.py)
├── utils/               # shared helpers (gridss.py, evidence.py, picard.py)
├── curator/             # schema.py + shared, schema-validated writer.py
├── annotations/         # OncoKB / CGC / hotspot / curation loaders, BED intervals
├── logging_config.py
└── exceptions.py
tests/                   # pytest suite + small synthetic fixtures/
```

Adding a new command = add `commands/<name>.py` holding the Typer function and
its logic (using the shared `CuratorWriter` for output and `commands/_shared.py`
for `--verbose/--quiet` and exit codes), then register it in `cli.py` with
`app.command()(<name>)`.


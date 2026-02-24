# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ICAMS (In-depth Characterization and Analysis of Mutational Signatures) is an R package for analyzing and visualizing mutational signatures from variant call files (VCFs). The package supports multiple mutation types (SBS, DBS, and indels) at various resolution levels and handles both count-based and density-based representations.

**Key publications**: Boot et al., Genome Research 2018 & 2020

**Current version**: 3.0.15 (see DESCRIPTION). Development branch: `4.0.0-branch`. Default remote branch: `v3.0.11-branch`.

## Development Commands

### Building and Testing

```r
# Load package for interactive development
devtools::load_all()

# Run all tests
devtools::test()

# Generate documentation from roxygen2 comments
devtools::document()

# Full package check (equivalent to R CMD check)
devtools::check()

# Build package
devtools::build()

# Install locally
devtools::install()
```

### Running Specific Tests

```r
# Run a single test file
testthat::test_file("tests/testthat/test_VCFsToCatalogs.R")

# Run tests matching a pattern
devtools::test(filter = "SBS96")
```

### CI/CD

GitHub Actions (`.github/workflows/R-CMD-check.yaml`) runs R CMD check on macOS, Windows, and Ubuntu (devel, release, oldrel-1). Build arguments: `--no-manual --compact-vignettes=gs+qpdf`. Currently triggered on pushes to `main`, `master`, `v3.0.10-branch`.

### Linting

`.lintr` config: UTF-8 encoding, 120-character line length, no assignment or object name linting enforced.

## Coding Style

- Exported functions use UpperCamelCase (`CanonicalizeID`, `ReadVCFs`); newer internal helpers use snake_case (`justify_indel`, `seg_simple`)
- Two-space indents, `<-` for assignment
- Roxygen2 with markdown enabled (`Roxygen: list(markdown = TRUE)` in DESCRIPTION)

## Architecture

### Core Data Structure: Catalogs

Catalogs are S3 objects extending matrices with specialized classes and attributes:

```r
# Class hierarchy
matrix -> SBS96Catalog / SBS192Catalog / SBS1536Catalog
          DBS78Catalog / DBS136Catalog / DBS144Catalog
          IndelCatalog / ID166Catalog
          COMPOSITECatalog (multiple types combined)

# Critical attributes
catalog.type    # "counts", "density", "counts.signature", "density.signature"
ref.genome      # "GRCh37", "GRCh38", "GRCm38"
region          # "genome", "exome", "transcript", "unknown"
abundance       # Named numeric vector of k-mer counts (for density)
class           # One or more catalog type classes
```

Always use `as.catalog()` to ensure proper attribute assignment. Never manually create catalogs without proper attributes.

### Module Organization

The R/ directory contains ~18K lines across 42 files. Major modules:

- **shiny_related_functions.R** (~2,600 lines): Interactive visualization (Shiny app)
- **utility_functions.R** (~2,260 lines): Catalog transformation, collapsing, manipulation
- **readVCF_etc.R** (~2,170 lines): VCF reading functions (ReadVCF, ReadStrelka*, ReadMutect*)
- **plot.R** (~2,060 lines): All visualization with S3 method dispatch
- **VCF_to_catalog_functions.R** (~1,660 lines): Core VCF-to-catalog pipeline orchestration
- **ID_functions.R** (~880 lines): Indel classification (microhomology, repeat detection)
- **infer_catalog_format.R** (~690 lines): Auto-detect catalog formats from files
- **other_catalog_formats.R** (~570 lines): SigProfiler/COSMIC format support
- **chromosome_name_functions.R** (~490 lines): Standardize chromosome naming ("chr1" vs "1")
- **sequence_context_functions.R** (~460 lines): Extract flanking sequence context from VCFs
- **strandbias_functions.R** (~400 lines): Transcriptional strand bias analysis

Indel classification has been modularized into several files:
- **categorize_1_justified_indel.R**: Koh indel categorization dispatch
- **justify_id_vcf.R**, **justify_indel.R**, **justify_indels_in_id_vcf_with_contexts.R**: Indel justification (canonical positioning)
- **annot_vcf_to_83_catalog.R**, **annot_vcf_to_89_catalog.R**, **annot_vcf_to_476_catalog.R**: Create catalogs for each classification scheme
- **gen_COSMIC_83_string.R**, **gen_koh_89_string.R**, **gen_koh_476_string.R**: Category string generation

### Rcpp

`src/segment_simple.cpp` provides a C++ segmentation implementation, exported as `segment_simple_cpp()`. Linked via `LinkingTo: Rcpp` in DESCRIPTION and `useDynLib(ICAMS, .registration = TRUE)` in NAMESPACE.

### VCF Processing Pipeline

```
ReadVCF() / ReadStrelkaXXXVCF() / ReadMutectVCF()   [readVCF_etc.R]
    |
MakeDataFrameFromVCF() -> Remove duplicates, filter chromosomes
    |
Annotate VCF:
  - AddSeqContext() -> Extract flanking sequences
  - AddTranscript() -> Add gene/strand info for transcriptional bias
  - AddRunInformation() -> Add repeat/microhomology info for indels
    |
CreateOneColXXXMatrix() -> Create matrix with canonical row order
    |
as.catalog() -> Add attributes
    |
cbind() -> Combine multiple samples
```

Functions return lists containing:
- `catalog`: The catalog matrix
- `discarded.variants`: Data frame explaining filtered variants (optional)
- `annotated.vcf`: Original VCF with added annotations (optional)

### S3 Method Dispatch Pattern

All catalog operations use S3 methods for type-specific behavior:

```r
# Generic function checks class and dispatches
PlotCatalog <- function(catalog, ...) {
  if (class(catalog)[1] %in% catalog.classes) {
    UseMethod("PlotCatalog")
  } else {
    # Convert to catalog first
  }
}

# Type-specific methods
PlotCatalog.SBS96Catalog <- function(catalog, ...) { ... }
PlotCatalog.DBS78Catalog <- function(catalog, ...) { ... }
```

Methods implemented: `PlotCatalog`, `PlotCatalogToPdf`, `WriteCatalog`, `[` (subsetting), `cbind`

### Exported Functions

53 exported functions plus 41 S3 methods. Key function groups:
- **VCF readers**: `ReadVCFs`, `ReadAndSplitVCFs`, `ReadAndSplitMutectVCFs`, `ReadAndSplitStrelkaSBSVCFs`, `SimpleReadVCF`
- **VCF-to-catalog pipelines**: `VCFsToCatalogs`, `VCFsToSBSCatalogs`, `VCFsToDBSCatalogs`, `VCFsToIDCatalogs`
- **Convenience wrappers**: `MutectVCFFilesToCatalog`, `StrelkaSBSVCFFilesToCatalog`, `StrelkaIDVCFFilesToCatalog`, plus `*ToPdf` and `*ToZipFile` variants
- **Catalog I/O**: `ReadCatalog`, `WriteCatalog`, `as.catalog`
- **Transformations**: `TransformCatalog`, `Collapse*` functions
- **Plotting**: `PlotCatalog`, `PlotCatalogToPdf`, `PlotTransBiasGeneExp`
- **VCF annotation**: `AnnotateSBSVCF`, `AnnotateDBSVCF`, `AnnotateIDVCF`
- **Indel functions**: `justify_indel`, `justify_id_vcf`, `categorize_1_justified_indel`, `annot_vcf_to_83_catalog`, `annot_vcf_to_89_catalog`, `annot_vcf_to_476_catalog`
- **Low-level indel**: `Canonicalize1Del`, `FindMaxRepeatDel`, `FindDelMH`
- **Utilities**: `revc`, `seg_simple`, `segment_simple_cpp`, `IsICAMSCatalog`, `GetFreebayesVAF`, `GetMutectVAF`, `GetStrelkaVAF`, `GetPCAWGConsensusVAF`

## Indel Classification System

ICAMS supports three indel classification schemes:

### COSMIC 83-Category System

Classification output format: `{DEL|INS}:{base|repeats|MH}:{length}:{count}`

Key functions: `FindMaxRepeatDel()` -> `FindDelMH()` -> `Canonicalize1Del()` / `Canonicalize1INS()` -> `CanonicalizeID()` (vectorized wrapper)

Examples:
- `DEL:T:1:2` = 1bp deletion of T in 2 tandem repeats
- `DEL:repeats:3:2` = 3bp deletion in 2 repeats
- `DEL:MH:5:3` = 5bp deletion with 3bp microhomology
- `INS:A:1:0` = 1bp insertion of A, no repeats

### Koh Classification Systems

Two additional, more granular schemes based on Koh et al.:

1. **Koh 89 categories** (`gen_koh_89_string.R`): Medium-resolution
2. **Koh 476 categories** (`gen_koh_476_string.R`): High-resolution

These consider additional factors: preceding/following bases, repeat unit count (R) with finer binning, and different treatment for insertions vs deletions.

**Indel justification**: Before classification, indels must be "justified" (canonically positioned) using `justify_indel()` and related functions. This ensures consistent classification for indels that can be represented in multiple ways within repeat sequences.

## Reference Genome Management

Three supported genomes:
- **GRCh37**: `BSgenome.Hsapiens.1000genomes.hs37d5`
- **GRCh38**: `BSgenome.Hsapiens.UCSC.hg38`
- **GRCm38**: `BSgenome.Mmusculus.UCSC.mm10`

BSgenome packages are Suggests (not Imports). Tests skip gracefully if genomes not installed:

```r
skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
```

## Key Dependencies

- **fastrc**: Provides `fast_rc` for fast reverse complement (replaces older custom `revc` implementation). Installed from GitHub via `Remotes: steverozen/fastrc`.
- **Rcpp**: C++ segmentation via `segment_simple.cpp`
- **Bioconductor**: BSgenome, Biostrings, GenomicRanges, IRanges, GenomeInfoDb
- **data.table**: High-performance data manipulation

## Abundance and Density Calculations

**Density** = counts / abundance (mutations per megabase of context)

Abundance varies by k-mer size, reference genome, region (genome/exome/transcript), and strand context. Pre-computed abundances stored in `sysdata.rda`. Use `TransformCatalog()` for conversions.

## Data Generation Pipeline

To regenerate internal data (only needed when updating reference genomes or transcript annotations):

1. **K-mer abundance files**: CSV files in `data-raw/new_masked_abundance/{GRCh37,GRCh38,GRCm38}/`
2. **Run scripts in order** (in `data-raw/code/`):
   - `load_abundance_from_files.R` -> Load k-mer counts
   - `create_catalogs.R` -> Generate catalog row orders
   - `create_order_for_DBS136_plotting.R` -> DBS136 plotting order
   - `create_ranges.R` -> Transcript ranges from GENCODE GTF
   - `create_gene_expression_data.R` -> Gene expression datasets
   - `create_ICAMS_SigPro_ID.R` -> ID format conversion matrices
   - `create_catalogs_COSMIC.R` -> COSMIC signature headers
3. **Save**: `save_global_variables.R` -> Creates `sysdata.rda`

## Multi-Format Support

ICAMS reads/writes multiple catalog formats:

- **ICAMS native**: CSV with mutation type row names
- **SigProfiler**: TSV format with different row ordering for ID catalogs
- **COSMIC**: CSV format from COSMIC signature database

`ReadCatalog()` auto-detects format. Conversion between ICAMS and SigProfiler ID formats uses `ICAMS.to.SigPro.ID` and `SigPro.to.ICAMS.ID` matrices.

## Parallel Processing

Many functions accept `num.of.cores` parameter. Implementation: `parallel::mclapply()` on Unix-like systems. Windows automatically falls back to sequential (num.of.cores forced to 1).

## Testing

- **76 test files** in `tests/testthat/` (testthat edition 3)
- Test data in `tests/testthat/testdata/` (VCFs, catalogs, expected outputs)
- Regression tests compare against saved `.csv` files
- Tests skip gracefully when BSgenome packages not installed

## Common Pitfalls

1. **Missing attributes**: Always use `as.catalog()` when creating catalogs programmatically
2. **Chromosome naming**: Use `CheckAndFixChrNames()` to handle "chr1" vs "1" inconsistencies
3. **BSgenome dependency**: Tests gracefully skip if reference genome not installed
4. **Discarded variants**: Check `discarded.variants` element of return lists to diagnose missing mutations
5. **Row order**: Each catalog type has a canonical row order; use `CatalogRowOrder()` and `CheckAndReorderRownames()`
6. **Abundance matching**: When transforming catalogs, ensure abundance matches region and reference genome

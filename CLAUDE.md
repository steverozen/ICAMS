# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ICAMS (In-depth Characterization and Analysis of Mutational Signatures) is an R package for analyzing and visualizing mutational signatures from variant call files (VCFs). The package supports multiple mutation types (SBS, DBS, and indels) at various resolution levels and handles both count-based and density-based representations.

**Key publications**: Boot et al., Genome Research 2018 & 2020

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

GitHub Actions runs R CMD check on macOS, Windows, and Ubuntu (devel, release, oldrel-1). Build arguments: `--no-manual --compact-vignettes=gs+qpdf`

## Architecture

### Core Data Structure: Catalogs

Catalogs are S3 objects extending matrices with specialized classes and attributes:

```r
# Class hierarchy
matrix → SBS96Catalog / SBS192Catalog / SBS1536Catalog
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

The R/ directory contains ~15K lines organized by functionality:

- **VCF_to_catalog_functions.R** (3,386 lines): Core VCF parsing and catalog creation
- **plot.R** (2,059 lines): All visualization with S3 method dispatch
- **shiny_related_functions.R** (2,301 lines): Interactive visualization
- **utility_functions.R** (2,261 lines): Catalog transformation, collapsing, manipulation
- **ID_functions.R** (803 lines): Indel classification (microhomology, repeat detection)
- **infer_catalog_format.R** (600 lines): Auto-detect catalog formats
- **other_catalog_formats.R** (493 lines): SigProfiler/COSMIC format support
- **chromosome_name_functions.R** (489 lines): Standardize chromosome naming
- **sequence_context_functions.R** (462 lines): Extract sequence context from VCFs
- **strandbias_functions.R** (404 lines): Transcriptional strand bias analysis
- **read_write_catalog.R** (244 lines): Catalog I/O

### VCF Processing Pipeline

```
ReadVCF() / ReadStrelkaXXXVCF() / ReadMutectVCF()
    ↓
MakeDataFrameFromVCF() → Remove duplicates, filter chromosomes
    ↓
Annotate VCF:
  - AddSeqContext() → Extract flanking sequences
  - AddTranscript() → Add gene/strand info for transcriptional bias
  - AddRunInformation() → Add repeat/microhomology info for indels
    ↓
CreateOneColXXXMatrix() → Create matrix with canonical row order
    ↓
as.catalog() → Add attributes
    ↓
cbind() → Combine multiple samples
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
# etc.
```

Methods implemented: `PlotCatalog`, `PlotCatalogToPdf`, `WriteCatalog`, `[` (subsetting), `cbind`

## Indel Classification System

The ID (insertion/deletion) classification is algorithmically complex:

1. **`FindMaxRepeatDel()`**: Count tandem repeat units in deletion context
2. **`FindDelMH()`**: Find microhomology at deletion boundaries
3. **`Canonicalize1Del()`**: Classify deletion type
   - 1bp deletions: Track deleted base + repeat count
   - Multi-bp in repeats: Track length + repeat count
   - Microhomology deletions: Track length + MH length
4. **`Canonicalize1INS()`**: Classify insertion type
   - 1bp insertions: Track inserted base + repeat count
   - Multi-bp insertions: Track length + repeat count
5. **`CanonicalizeID()`**: Vectorized wrapper for full VCF

Classification output format: `{DEL|INS}:{base|repeats|MH}:{length}:{count}`

Examples:
- `DEL:T:1:2` = 1bp deletion of T in 2 tandem repeats
- `DEL:repeats:3:2` = 3bp deletion in 2 repeats
- `DEL:MH:5:3` = 5bp deletion with 3bp microhomology
- `INS:A:1:0` = 1bp insertion of A, no repeats

## Reference Genome Management

Three supported genomes:
- **GRCh37**: `BSgenome.Hsapiens.1000genomes.hs37d5`
- **GRCh38**: `BSgenome.Hsapiens.UCSC.hg38`
- **GRCm38**: `BSgenome.Mmusculus.UCSC.mm10`

BSgenome packages are Suggests (not Imports) to reduce installation burden. Tests skip gracefully if genomes not installed:

```r
skip_if("" == system.file(package = "BSgenome.Hsapiens.1000genomes.hs37d5"))
```

Users must install separately: `BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")`

## Abundance and Density Calculations

**Density** = counts / abundance (mutations per megabase of context)

Abundance varies by:
- K-mer size (2bp, 3bp, 4bp, 5bp)
- Reference genome (GRCh37, GRCh38, GRCm38)
- Region (genome, exome, transcript)
- Strand context (stranded for transcribed regions)

Pre-computed abundances stored in `sysdata.rda` (internal package data). Use `TransformCatalog()` to convert between count-based and density-based representations or to normalize across different regions.

## Data Generation Pipeline

To regenerate internal data (only needed when updating reference genomes or transcript annotations):

1. **K-mer abundance files**: CSV files in `data-raw/new_masked_abundance/{GRCh37,GRCh38,GRCm38}/`
2. **Run scripts in order**:
   - `load_abundance_from_files.R` → Load k-mer counts
   - `create_catalogs.R` → Generate catalog row orders
   - `create_order_for_DBS136_plotting.R` → DBS136 plotting order
   - `create_ranges.R` → Transcript ranges from GENCODE GTF
   - `create_gene_expression_data.R` → Gene expression datasets
   - `create_ICAMS_SigPro_ID.R` → ID format conversion matrices
   - `create_catalogs_COSMIC.R` → COSMIC signature headers
3. **Save**: `save_global_variables.R` → Creates `sysdata.rda`

The pipeline is documented in `data-raw/code/save_global_variables.R`

## Multi-Format Support

ICAMS reads/writes multiple catalog formats:

- **ICAMS native**: CSV with mutation type row names
- **SigProfiler**: TSV format with different row ordering for ID catalogs
- **COSMIC**: CSV format from COSMIC signature database

`ReadCatalog()` auto-detects format. Use `ConvertCatalogToSigProfilerFormat()` for exports.

Conversion between ICAMS and SigProfiler ID formats uses matrices:
- `ICAMS.to.SigPro.ID` (83x83 sparse matrix)
- `SigPro.to.ICAMS.ID` (83x83 sparse matrix)

## Parallel Processing

Many functions accept `num.of.cores` parameter:

```r
ReadVCFs(..., num.of.cores = 4)
VCFsToCatalogs(..., num.of.cores = 4)
```

Implementation: `parallel::mclapply()` on Unix-like systems (Linux, macOS)

**Windows limitation**: No fork support, automatically falls back to sequential processing (num.of.cores forced to 1)

## Testing Conventions

- **63 test files** in `tests/testthat/`
- Each catalog type has dedicated tests for plotting, I/O, transformations
- Test data: `tests/testthat/testdata/` (VCFs, catalogs, expected outputs)
- Regression tests compare against saved `.csv` files
- Use `expect_equal()` with tolerance for numeric comparisons
- Test parallel execution separately from sequential

## Documentation

- **Roxygen2**: All exported functions documented with `@title`, `@param`, `@return`, `@export`, `@examples`
- Generate with: `devtools::document()`
- Reference manual: See `data-raw/ICAMS_3.0.9.pdf` (version-specific)
- Function count: 88 exported functions

## Version and Branch Strategy

- Current version: **3.0.11** (see DESCRIPTION)
- Active branch: **v3.0.11-branch**
- CRAN releases use version-tagged branches
- Main development typically on version branches, not master

## Common Pitfalls

1. **Missing attributes**: Always use `as.catalog()` when creating catalogs programmatically
2. **Chromosome naming**: Use `CheckAndFixChrNames()` to handle "chr1" vs "1" inconsistencies
3. **BSgenome dependency**: Tests gracefully skip if reference genome not installed
4. **Discarded variants**: Check `discarded.variants` element of return lists to diagnose missing mutations
5. **Row order**: Each catalog type has a canonical row order; use `CatalogRowOrder()` and `CheckAndReorderRownames()`
6. **Abundance matching**: When transforming catalogs, ensure abundance matches region and reference genome

## Citation

If adding features that should be cited, follow format in README:

> Rozen SG, Jiang NH, Boot A, Liu M, Wu Y, Huang MN, Chang JG (2025). ICAMS: In-depth Characterization and Analysis of Mutational Signatures. R package version 3.0.11, https://CRAN.R-project.org/package=ICAMS.

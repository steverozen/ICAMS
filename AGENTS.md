# Repository Guidelines

## Project Structure & Module Organization
The package follows standard R layout. Core functions reside in `R/`, with roxygen comments feeding `man/` for documentation and `NAMESPACE` for exports. Unit tests live under `tests/testthat/`, while `inst/extdata/` and `inst/example.csv` hold example catalogs and helper spreadsheets referenced in tests and vignettes. Data stored with the package is bundled in `data/`, and scripts that regenerate those objects sit in `data-raw/`. The pkgdown site builds into `docs/`; update it only after passing package checks.

## Build, Test, and Development Commands
Use `Rscript -e "devtools::load_all()"` to refresh the package environment during interactive work. Run `Rscript -e "devtools::test()"` before commits to execute the `testthat` suite. For a release check, call `R CMD build .` followed by `R CMD check ICAMS_<version>.tar.gz`. Update the documentation site with `Rscript -e "pkgdown::build_site()"` once checks succeed.

## Coding Style & Naming Conventions
Code is written in base R with two-space indents and `<-` for assignments. Exported functions typically use UpperCamelCase (`CanonicalizeID`), while internal helpers may use snake_case when scoped to a single file. Keep arguments descriptive and default-heavy to mirror existing APIs. Document every exported function with roxygen2 blocks, placing examples inside `\dontrun{}` when they rely on large external files. Prefer vectorised operations over loops where practical and reuse helpers from `utility_functions.R` before adding new ones.

## Testing Guidelines
Tests use `testthat` snapshot and expectation helpers; add new files as `tests/testthat/test-<topic>.R`. Reuse sample VCFs or catalogs from `inst/extdata/` to keep runtimes short. When tests create temporary artifacts, write them beneath `tests/tmp/` and clean them up within the test. Aim to cover new branches and error paths; quick coverage checks can be run via `Rscript -e "covr::report()"` when available.

## Commit & Pull Request Guidelines
Recent history favors concise, lower-case summaries (for example, “minor potential bug”). Follow that style, keeping the first line under 72 characters and expanding rationale in the body if needed. Pull requests should link related issues, outline interface changes, and list any new data files. Include the commands run (tests, build, pkgdown) and attach screenshots when modifying plot output to help reviewers verify visual changes.

## Security & Configuration Tips
Genome reference data can be large and sensitive; never commit raw patient VCFs. Store credentials and private paths outside the repo and reference them via environment variables. When generating new datasets in `data-raw/`, document provenance and checksum expectations so downstream users can reproduce the build.

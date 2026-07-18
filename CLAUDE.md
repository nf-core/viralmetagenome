# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**nf-core/viralmetagenome** is a Nextflow-based bioinformatics pipeline for reconstructing consensus viral genomes and identifying intra-host variants from metagenomic or enriched sequencing data. It's an official nf-core pipeline following the nf-core template 4.0.2, built with Nextflow ≥25.10.4.

### Key Characteristics
- **Type**: Bioinformatics analysis pipeline (Nextflow)
- **Primary Language**: Groovy (Nextflow DSL2)
- **Supported Execution**: Docker, Singularity, Conda, AWS
- **Testing Framework**: nf-test
- **CI/CD**: GitHub Actions

## Architecture & Structure

### High-Level Pipeline Design

The pipeline is organized as a directed acyclic graph (DAG) of reusable modules and subworkflows:

1. **Read QC & Preprocessing** → **Metagenomic Classification** → **Assembly** → **Polishing** → **Variant Calling** → **Consensus Generation** → **Validation & Annotation** → **Reporting**

Optional branches include host removal, contig scaffolding, BLAST reference identification, taxonomic pre-clustering, and iterative refinement.

### Directory Structure

- **`main.nf`**: Entry point that initializes the pipeline, runs validation, and calls the main workflow
- **`workflows/viralmetagenome.nf`**: Main analysis workflow orchestrating all processing steps (25KB, highly modular)
- **`subworkflows/`**:
  - `local/`: Custom pipeline-specific workflows (e.g., preprocessing chains, assembly variants, polishing loops)
  - `nf-core/`: Reusable nf-core subworkflows (imported from nf-core modules)
- **`modules/`**:
  - `local/`: Custom tools not in nf-core (14 custom modules: BLAST filtering, clustering extraction, variant conversion, SSPACE scaffolding, etc.)
  - `nf-core/`: ~100 curated modules from nf-core (managed via `modules.json`)
- **`conf/`**: Nextflow configuration files (profiles, institution-specific settings)
- **`tests/`**: nf-test test definitions and snapshots (.nf.test files)
- **`docs/`**: User-facing documentation
- **`bin/`**: Utility scripts (Python/Bash)
- **`assets/`**: Static resources (logos, example configs)

### Configuration System

- **`nextflow.config`**: ~300+ parameters controlling every aspect (inputs, tool selection, thresholds, resource limits)
- **`nextflow_schema.json`**: Schema defining parameters and their valid values (used for validation and UI generation)
- **`.nf-core.yml`**: Metadata including linting rules and template version

### Module & Subworkflow Organization

Modules follow nf-core conventions:
- Each module in `modules/nf-core/<tool>/<subtool>/main.nf` wraps a single bioinformatics tool
- Modules declare `input`, `output`, and `script` blocks in DSL2 syntax
- Local modules in `modules/local/<task>/` handle pipeline-specific logic

Subworkflows compose modules into logical units:
- Local subworkflows chain multiple modules with conditional logic
- Enable parameter-driven tool selection (e.g., `trim_tool: 'fastp' | 'trimmomatic'`, `assemblers: 'spades,megahit'`)
- Handle branching for optional steps

## Development Workflow

### Setting Up

```bash
# Clone and navigate
git clone https://github.com/nf-core/viralmetagenome.git
cd viralmetagenome

# Fetch the designated development branch (replace with actual branch name if different)
git fetch origin
git checkout claude/claude-md-docs-ga3uu4  # Or your assigned feature branch
```

### Running the Pipeline

```bash
# Test run with minimal data (recommended for development)
nextflow run . -profile test,docker -resume

# With custom parameters
nextflow run . -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --skip_preprocessing false
```

### Testing

```bash
# Run all nf-tests
nf-test test

# Run a specific test
nf-test test tests/default.nf.test

# Run tests with a specific profile
nf-test test --profile docker

# Update snapshots after intentional changes
nf-test test --update-snapshot
```

Available test profiles in `tests/`:
- `default.nf.test` — Full pipeline test with typical parameters
- `test_minimal.nf.test` — Minimal input (fast smoke test)
- `test_full.nf.test` — Full-size dataset (AWS CI only)
- `test_umi.nf.test` — UMI-based deduplication flow
- `test_group.nf.test` — Sample grouping scenarios
- `test_fail_db.nf.test` — Database failure handling
- `test_fail_mapped.nf.test` — Read mapping failure handling

### Linting & Validation

```bash
# Check Nextflow syntax and nf-core standards
nf-core lint .

# Validate pipeline schema and parameters
nextflow run . --help  # Shows all parameters with descriptions

# Pre-commit hooks (if configured)
pre-commit run --all-files
```

### Commits & Branching

- Develop on the designated feature branch (`claude/claude-md-docs-ga3uu4`)
- Write clear, descriptive commit messages
- Push to the same branch: `git push -u origin claude/claude-md-docs-ga3uu4`
- Do NOT force-push unless explicitly authorizing it
- Use `git status` before commits to verify staging

## Key Conventions

### Parameter Management

All pipeline options are defined in `nextflow.config` and validated against `nextflow_schema.json`. Key parameter types:

- **Tool Selection**: `trim_tool`, `assemblers`, `cluster_method`, `read_classifiers` — comma-separated lists enabling multiple tools in sequence or branching
- **Skip Flags**: `skip_*` parameters disable entire workflows (e.g., `skip_preprocessing`, `skip_assembly`, `skip_polishing`)
- **Database URLs**: Remote URLs (S3, HTTP) or local paths for Kraken2, Kaiju, BLAST references
- **Thresholds**: `identity_threshold`, `perc_reads_contig`, `min_contig_size`, etc. — control filtering and clustering

### Module Naming & Patterns

- Local modules use `SNAKE_CASE` for workflow step names (e.g., `EXTRACT_CLUSTER`, `IVAR_VARIANTS_TO_VCF`)
- nf-core modules maintain their canonical names (e.g., `FASTQC`, `BOWTIE2_ALIGN`)
- Input channels: use descriptive names indicating data type (e.g., `reads_merged`, `contigs_clustered`)

### Channel Operations

Channels carry tuples `[meta, data]` where:
- `meta` is a Groovy Map containing sample ID, file paths, group info, etc.
- `data` is the actual file(s) or values

Key operations:
- `.groupTuple()` — Combine outputs from parallel samples
- `.branch()` — Conditional routing based on parameters or metadata
- `.mix()` — Merge channels
- `.flatten()` — Expand nested structures

### Workflow Structure

- Main workflow in `workflows/viralmetagenome.nf` takes samplesheet and orchestrates 20+ processing steps
- Subworkflows handle complex branching (e.g., separate assembly methods, optional polishing iterations)
- Each step is conditional on skip flags and tool selection parameters
- Use `emit:` to expose outputs for downstream use

## Common Tasks

### Adding a New Tool or Step

1. **Check if nf-core module exists** — Browse [nf-core/modules](https://github.com/nf-core/modules) for existing tools
2. **Import and include** in `workflows/viralmetagenome.nf`:
   ```groovy
   include { TOOL_NAME } from '../modules/nf-core/tool/subtool'
   ```
3. **Define parameters** in `nextflow.config` for tool-specific options
4. **Add schema entries** in `nextflow_schema.json` for validation
5. **Integrate into subworkflow** or main workflow with proper channel handling
6. **Add tests** by updating relevant `.nf.test` file with new output checks

### Adding a Custom Module

1. Create `modules/local/<task>/` with `main.nf`, `meta.yml`, and `tests/` subdirectory
2. Follow DSL2 structure:
   ```groovy
   process MODULE_NAME {
       input:
       tuple val(meta), path(input_file)
       output:
       tuple val(meta), path("*.out"), emit: output
       script:
       // Implementation
   }
   ```
3. Document in `meta.yml` (input/output descriptions, authors)
4. Add to test snapshots if used in pipeline

### Modifying Parameters

1. Edit `nextflow.config` with new parameter and default value
2. Add entry to `nextflow_schema.json` with type, description, and validation rules
3. Update affected workflow/subworkflow to use the parameter
4. Run `nf-core lint` to validate schema
5. Update tests if parameter affects behavior

### Debugging Pipeline Issues

- **Check logs**: `nextflow log <run-name> -f stdout | tail -100`
- **Inspect work directory**: `cat work/<hash>/.command.log` for individual process output
- **Run with `-resume`** to skip completed tasks and restart from failure
- **Enable debug mode**: `nextflow run . -profile test,docker -resume --debug`
- **Test a single module**: Use nf-test or write a minimal workflow

## Important Files & Their Roles

| File | Purpose |
|------|---------|
| `main.nf` | Pipeline entry point; handles parameter validation and initialization |
| `workflows/viralmetagenome.nf` | Main orchestration; defines DAG of all processing steps |
| `nextflow.config` | Parameter defaults and profile definitions |
| `nextflow_schema.json` | Parameter schema for validation and UI generation |
| `modules.json` | Registry of nf-core modules and their versions |
| `tests/*.nf.test` | Test definitions and assertions |
| `.github/workflows/*.yml` | CI/CD pipelines (linting, nf-test, AWS full tests) |

## Testing Strategy

The pipeline uses **nf-test** for end-to-end validation:

- **Snapshot-based testing**: Each test compares actual outputs to stored snapshots
- **Profiles**: Tests run against the `test` profile with small datasets
- **Multiple scenarios**: Tests cover happy paths, edge cases (UMI, grouping), and failure modes
- **CI integration**: GitHub Actions run all tests on PR; full-size tests on AWS

When modifying workflows or adding steps:
1. Run affected tests to ensure they still pass
2. Update snapshots if changes are intentional: `nf-test test --update-snapshot`
3. Add new tests for significant new functionality

## Working with nf-core Standards

This is an official nf-core pipeline, so:

- Follow [nf-core module guidelines](https://nf-co.re/developers/module_guidelines)
- Use `nf-core lint` to catch style and structure issues
- Keep local modules minimal; prefer nf-core modules when available
- Maintain backward compatibility with existing parameters when possible
- Document all parameters in `nextflow_schema.json`
- Use consistent naming: PascalCase for processes, snake_case for parameters and files

## Resources

- **Pipeline Docs**: https://nf-co.re/viralmetagenome
- **nf-core Docs**: https://nf-co.re/docs
- **Nextflow Docs**: https://nextflow.io/docs
- **GitHub Issues**: https://github.com/nf-core/viralmetagenome/issues
- **Slack**: #viralmetagenome channel in nf-core Slack

# nf-core/viralmetagenome: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [date]

### `Added`

- Add argument to control the maximum mpileup depth in `custom_mpileup.py` script ([#176](https://github.com/nf-core/viralmetagenome/pull/176)) (by @Joon-Klaps)
- Removing redundant `samtools_sort` after `BAM_DEDUPLICATE` ([#177](https://github.com/nf-core/viralmetagenome/pull/177)) (by @Joon-Klaps)
- Removing `seqkit replace` and move logic to `blast_filter.py` ([#178](https://github.com/nf-core/viralmetagenome/pull/178)) (by @Joon-Klaps)
- Giving Kaiju more memory ([#179](https://github.com/nf-core/viralmetagenome/pull/179)) (by @Joon-Klaps)
- Add option for sample merging based on group column in samplesheet ([#180](https://github.com/nf-core/viralmetagenome/pull/180)) (by @Joon-Klaps)
- Creating testdatasets using nf-core repo ([#183](https://github.com/nf-core/viralmetagenome/pull/183)) (by @Joon-Klaps)
- Add option to annotate snps with Snpeff for mapping constraint route ([#186](https://github.com/Joon-Klaps/viralgenie/pull/186)) (by @Joon-Klaps)
- Add `nf-tests` for `test` profile ([#189](https://github.com/nf-core/viralmetagenome/pull/189)) (by @Joon-Klaps)
- Update docs ([#200](https://github.com/nf-core/viralmetagenome/pull/200)) (by @Joon-Klaps)
- Template update for nf-core/tools v3.3.2 ([#202](https://github.com/nf-core/viralmetagenome/pull/202)) (by @Joon-Klaps)
- Add support for both string as integer inputs in samplesheets ([#230](https://github.com/nf-core/viralmetagenome/pull/230)) (by @Joon-Klaps)

### `Fixed`

- Fix local module `ivar_variants_to_vcf` to handle empty tsv files ([#197](https://github.com/Joon-Klaps/viralmetagenome/pull/197/)) (by @Joon-Klaps)
- Migrate lib dir functions to utils_nfcore_viralmetagenome_pipeline ([#194](https://github.com/nf-core/viralmetagenome/pull/194)) (by @Joon-Klaps)
- Fix inconsistent dependency versions across modules ([#208](https://github.com/nf-core/viralmetagenome/pull/208)) (by @Joon-Klaps)
- Fix conda issue unrecognized arguments: --mkdir ([#210](https://github.com/nf-core/viralmetagenome/pull/210)) (by @Joon-Klaps)
- Fix writing no sequence for select_reference.py to the first reference of the multifasta ([#214](https://github.com/nf-core/viralmetagenome/pull/214)) (by @Joon-Klaps)
- Fix main language detection to ignore generated files ([#224](https://github.com/nf-core/viralmetagenome/pull/224)) (by @Joon-Klaps)

### `Dependencies`

### `Deprecated`

- Refactor `params.skip_annotation` to `params.skip_consensus_annotation`([#181](https://github.com/nf-core/viralmetagenome/pull/181)) (by @Joon-Klaps)
- Deprecate `params.skip_nocov_to_reference` ([#212](https://github.com/nf-core/viralmetagenome/pull/212)) (by @Joon-Klaps)
- Deprecate `BWAMEM` as mapping tool ([#212](https://github.com/nf-core/viralmetagenome/pull/212)) (by @Joon-Klaps)

## v0.1.2 - 2025-02-28

Second release of thenf-core/viralmetagenome pipeline. Focusing on user experience and bug fixes.

### `Added`

- Set default umitools dedup strategy to cluster ([#126](https://github.com/nf-core/viralmetagenome/pull/126)) (by @Joon-Klaps)
- Include both krakenreport &nodes.dmp in taxonomy filtering ([#128](https://github.com/nf-core/viralmetagenome/pull/128)) (by @Joon-Klaps)
- Add Sspace indiv to each assembler seperatly ([#132](https://github.com/nf-core/viralmetagenome/pull/132)) (by @Joon-Klaps)
- Add read & contig decomplexification using prinseq++ ([#133](https://github.com/nf-core/viralmetagenome/pull/133)) (by @Joon-Klaps)
- Add option to filter contig clusters based on cumulative read coverage ([#138](https://github.com/nf-core/viralmetagenome/pull/138)) (by @Joon-Klaps)
- Reffurbish mqc implementation ([#139](https://github.com/nf-core/viralmetagenome/pull/139)) (by @Joon-Klaps)
- Adding mash-screen output to result table ([#140](https://github.com/nf-core/viralmetagenome/pull/140)) (by @Joon-Klaps)
- Add logic to allow samples with no reference hits to be analysed ([#141](https://github.com/nf-core/viralmetagenome/pull/141)) (by @Joon-Klaps)
- Add visualisation for hybrid scaffold ([#143](https://github.com/nf-core/viralmetagenome/pull/143)) (by @Joon-Klaps)
- Add new module to inculde custom mpileup-vcf file for intra-host analyses ([#151](https://github.com/nf-core/viralmetagenome/pull/151)) (by @Joon-Klaps)
- Update docs ([#150](https://github.com/nf-core/viralmetagenome/pull/150)) (by @Joon-Klaps)
- Make custom-mpileup.py postion 1 index based and not 0 index to follow bcftools ([#153](https://github.com/nf-core/viralmetagenome/pull/153)) (by @Joon-Klaps)
- Update docs for more streamlined docs & figures ([#154](https://github.com/nf-core/viralmetagenome/pull/154)) (by @Joon-Klaps)
- Add column in custom mpileup - Shannon entropy ([#156](https://github.com/nf-core/viralmetagenome/pull/156)) (by @Joon-Klaps)
- Constrain -> Constraint & further python script debugging ([#161](https://github.com/nf-core/viralmetagenome/pull/161)) (by @Joon-Klaps)
- include empty samples in multiqc sample overview ([#162](https://github.com/nf-core/viralmetagenome/pull/162)) (by @Joon-Klaps)
- Include samtools stats pre dedup & post dedup in overview tables ([#163](https://github.com/nf-core/viralmetagenome/pull/163)) (by @Joon-Klaps)
- adding prokka for gene detection & annotation ([#165](https://github.com/nf-core/viralmetagenome/pull/165)) (by @Joon-Klaps)

### `Fixed`

- OOM with longer contigs for nocov_to_reference, uses more RAM now ([#125](https://github.com/nf-core/viralmetagenome/pull/125)) (by @Joon-Klaps)
- fixing null output from global prefix ([#147](https://github.com/nf-core/viralmetagenome/pull/147)) (by @Joon-Klaps)
- Fix empty filtered clusters ([#148](https://github.com/nf-core/viralmetagenome/pull/148)) (by @Joon-Klaps)
- Fixing missing columns from general stats & add general stats sample filtering ([#149](https://github.com/nf-core/viralmetagenome/pull/149)) (by @Joon-Klaps)
- process.shell template fix ([#157](https://github.com/nf-core/viralmetagenome/pull/157)) (by @Joon-Klaps) - see also [nf-core/tools #3416](https://github.com/nf-core/tools/pull/3416)

### `Parameters`

- New parameter mmseqs_cluster_mode default to 0 ([#130](https://github.com/nf-core/viralmetagenome/pull/130)) (by @Joon-Klaps) **DEPRECATED**
- Refactor module arguments to pipeline arguments ([#166](https://github.com/nf-core/viralmetagenome/pull/166)) (by @Joon-Klaps)

## v0.1.1 - 2024-05-08

Initial release of nf-core/viralmetagenome, created with the [nf-core](https://nf-co.re/) template.

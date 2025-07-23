# nf-core/viralmetagenome: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.3dev - [date]

### Enhancement

- Add argument to control the maximum mpileup depth in `custom_mpileup.py` script ([#176](https://github.com/nf-core/viralmetagenome/pull/176))
- Removing redundant `samtools_sort` after `BAM_DEDUPLICATE` ([#177](https://github.com/nf-core/viralmetagenome/pull/177))
- Removing `seqkit replace` and move logic to `blast_filter.py` ([#178](https://github.com/nf-core/viralmetagenome/pull/178))
- Giving Kaiju more memory ([#179](https://github.com/nf-core/viralmetagenome/pull/179))
- Add option for sample merging based on group coulumn in samplesheet ([#180](https://github.com/nf-core/viralmetagenome/pull/180))
- Creating testdatasets using nf-core repo ([#183](https://github.com/nf-core/viralmetagenome/pull/183))
- Add option to annotate snps with Snpeff for mapping constraint route ([#186](https://github.com/Joon-Klaps/viralgenie/pull/186))
- Template update for nf-core/tools v3.3.2 ([#202](https://github.com/nf-core/viralmetagenome/pull/202))

### `Fixed`

- Fix local module `ivar_variants_to_vcf` to handle empty tsv files ([#197](https://github.com/Joon-Klaps/viralmetagenome/pull/197/))
- Migrate lib dir functions to utils_nfcore_viralmetagenome_pipeline ([#194](https://github.com/nf-core/viralmetagenome/pull/194))

### `Parameters`

## v0.1.2 - 2025-02-28

Second release of the viralmetagenome pipeline. Focusing on user experience and bug fixes.

### `Enhancement`

- Set default umitools dedup strategy to cluster ([#126](https://github.com/nf-core/viralmetagenome/pull/126))
- Include both krakenreport &nodes.dmp in taxonomy filtering ([#128](https://github.com/nf-core/viralmetagenome/pull/128))
- Add Sspace indiv to each assembler seperatly ([#132](https://github.com/nf-core/viralmetagenome/pull/132))
- Add read & contig decomplexification using prinseq++ ([#133](https://github.com/nf-core/viralmetagenome/pull/133))
- Add option to filter contig clusters based on cumulative read coverage ([#138](https://github.com/nf-core/viralmetagenome/pull/138))
- Reffurbish mqc implementation ([#139](https://github.com/nf-core/viralmetagenome/pull/139))
- Adding mash-screen output to result table ([#140](https://github.com/nf-core/viralmetagenome/pull/140))
- Add logic to allow samples with no reference hits to be analysed ([#141](https://github.com/nf-core/viralmetagenome/pull/141))
- Add visualisation for hybrid scaffold ([#143](https://github.com/nf-core/viralmetagenome/pull/143))
- Add new module to inculde custom mpileup-vcf file for intra-host analyses ([#151](https://github.com/nf-core/viralmetagenome/pull/151))
- Update docs ([#150](https://github.com/nf-core/viralmetagenome/pull/150))
- Make custom-mpileup.py postion 1 index based and not 0 index to follow bcftools ([#153](https://github.com/nf-core/viralmetagenome/pull/153))
- Update docs for more streamlined docs & figures ([#154](https://github.com/nf-core/viralmetagenome/pull/154))
- Add column in custom mpileup - Shannon entropy ([#156](https://github.com/nf-core/viralmetagenome/pull/156))
- Constrain -> Constraint & further python script debugging ([#161](https://github.com/nf-core/viralmetagenome/pull/161))
- include empty samples in multiqc sample overview ([#162](https://github.com/nf-core/viralmetagenome/pull/162))
- Include samtools stats pre dedup & post dedup in overview tables ([#163](https://github.com/nf-core/viralmetagenome/pull/163))
- adding prokka for gene detection & annotation ([#165](https://github.com/nf-core/viralmetagenome/pull/165))

### `Fixed`

- OOM with longer contigs for nocov_to_reference, uses more RAM now ([#125](https://github.com/nf-core/viralmetagenome/pull/125))
- fixing null output from global prefix ([#147](https://github.com/nf-core/viralmetagenome/pull/147))
- Fix empty filtered clusters ([#148](https://github.com/nf-core/viralmetagenome/pull/148))
- Fixing missing columns from general stats & add general stats sample filtering ([#149](https://github.com/nf-core/viralmetagenome/pull/149))
- process.shell template fix ([#157](https://github.com/nf-core/viralmetagenome/pull/157)) - see also [nf-core/tools #3416](https://github.com/nf-core/tools/pull/3416)

### `Parameters`

- New parameter mmseqs_cluster_mode default to 0 ([#130](https://github.com/nf-core/viralmetagenome/pull/130)) **DEPRECATED**
- Refactor module arguments to pipeline arguments ([#166](https://github.com/nf-core/viralmetagenome/pull/166))

## v0.1.1 - 2024-05-08

Initial release of nf-core/viralmetagenome, created with the [nf-core](https://nf-co.re/) template.

# WES-filtering

Filtering and QC of the somatic SNV/CNV/SV calls of a
WES-snakemake run. Every call is labelled with a FILTER value
(PASS or the reason it failed), and the filtered tables plus a per-sample `metadata.csv`
are written for downstream analyses.

The code is project-agnostic. It sits at `src/WES-filtering` of a project organized as
`data/ results/ src/ metadata/`. Each run is a dated experiment directory under `results/`
that holds its own `config.yaml` and outputs.

## Setup

```bash
cd <project>
git submodule add https://github.com/dmkv1/WES-filtering.git src/WES-filtering
```

The project commit then pins the WES-filtering version. To update:
`git submodule update --remote src/WES-filtering`.

R packages: data.table, dplyr, tidyr, openxlsx, qs2, yaml, ggplot2, ggh4x,
ggbeeswarm, scales.

WES-snakemake must run with `keep_strand_bias_calls: true` (strand-bias rescue on) and write
`results/combined/combined_{snvs,cnvs,svs,qc}.tsv` and `results/metadata/samples.tsv`.

## Running an experiment

```bash
src/WES-filtering/new_experiment [--topic TOPIC] [--results-dir DIR] [--date YYYY-MM-DD]
cd results/<today>_WES-filtering[-TOPIC]
./runall --sample-annotations --load --plot 2>&1 | tee runall.log
```

`new_experiment` copies `template/` to `results/<date>_WES-filtering[-TOPIC]/`. It sets
`results_dir` to `--results-dir`, or else to the chronologically last
`results/*/WES-snakemake/results`. The experiment runs on the defaults as is. Edit its
`config.yaml` for project-specific settings, and note what changed in its `README.md`.

`runall` flags (combinable, always executed in this order):
- `--sample-annotations`: rebuild `<metadata_dir>/sample_annotations.csv`, keeping edits.
- `--load`: write `metadata.csv` and `combined_{snvs,cnvs,svs}.qs2`.
- `--plot`: diagnostic plots (needs `--load` outputs).

`--load` also writes the WES-filtering version (`git describe`) to
`WES-filtering.version`, and the experiment name to `<metadata_dir>/metadata.ver`.

## Config

An experiment's `config.yaml` is layered over `template/config.yaml`, which holds the
defaults and documents every key (`config.R`). Nested maps (`snv_filter_params`, ...) merge
key by key. `sv_driver_loci` is never inherited. Relative paths resolve against the
experiment directory.

Defaults:
- `metadata_dir: ../../metadata`.
- `group1: sample_type`, no `group2`, no `disease`, no `patient_annotations`.
- `cytoband`: the bundled `resources/cytoBand.hg38.txt.gz` (UCSC hg38, downloaded
  2026-09-23).
- On: strand-bias rescue (AF >= 0.2, no recurrence in other patients) and the HLA class I
  high-AF flag.
- Off, because they need project gene lists: germline rescue, clustered-call rescue, the
  hotspot recurrence exemption, and the gene-list export.
- No SV driver loci. The template has commented DLBCL (BCL2/BCL6/MYC) and MCL (CCND1-3)
  entries.

Project-level settings usually changed:
- `disease`
- `group1` / `group2`
- `patient_annotations: patients.csv`
- `snv_genelist` with `snv_genelist_diagnosis`, `hotspot_list` and
  `clustered_rescue_genes`, which enable the rescues
- `sv_driver_loci`

## Outputs

In the experiment directory:
- `analysis-SNVs/combined_snvs.qs2`, `SNV_counts.csv`, `filters_SNV.png`,
  `af_SNV_<Disease>.png` (`af_SNV.png` without `disease`), `filter_SNVs.log`
- `analysis-CNVs/combined_cnvs.qs2`, `cn0_rel_depth_distribution.png`, `cnv_tracks.png`,
  `filter_CNVs.log`
- `analysis-SV/combined_svs.qs2`, `filter_SVs.log`

`<metadata_dir>/metadata.csv` has one row per tumor sample and is regenerated on every
`--load`. Its columns are `ID, sID, sample_type, gender, capture_kit, tumor_fraction,
Disease, calling_mode, group1, [group2,]`, then the patient_annotations and
sample_annotations columns, then `QC`. `QC` is `excluded` when `exclude_reason` is set in
sample_annotations, else `high_SNV_count` (> 2000 PASS SNVs), else `PASS`. Downstream
analyses keep `QC == "PASS"`.

Hand-edited inputs in `metadata_dir`:
- `patient_annotations` (csv keyed by `ID`). Columns other than `notes` are joined into
  `metadata.csv`. A `Disease` column overrides `disease`. A column repeating a
  `samples.tsv` column must agree with it.
- `sample_annotations` (csv keyed by `sID`; `exclude_reason`, `sample_notes`, more
  allowed). Rebuilt by `--sample-annotations`; edits are kept.

## Filters

SNVs (`analysis-SNVs/filter_SNVs.R`), in order:
1. Rescues of calls the pipeline keeps with a Mutect2 label. Each rescued call becomes PASS
   and goes through the filters below:
   - strand-bias rescue (`strandbias_rescue.R`)
   - germline rescue (tumor-only samples)
   - clustered-call rescue (aSHM)
   Unrescued calls keep their Mutect2 label.
2. Gene blocklist (`flags_genes.csv`: `FLAGS_large_gene`, `FLAGS_polymorphic`,
   `FLAGS_repetitive`); HLA class I calls at AF > 0.9 (`FLAGS_high_AF`).
3. `low_alt_reads` (AD_ALT <= 3), `low_AF` (< 0.05), `tumor_only_low_AF` (< 0.2),
   `germline_ctrl` (AF_ctrl >= 0.05 and AF_ctrl/AF_tumor >= 0.25), `germline_gnomad`
   (gnomAD AF > 0.05).
4. Per-sample QC: > 2000 PASS SNVs makes the sample `high_SNV_count`.
5. Consequence allowlist (`SNV_consequence_terms.tsv`): `low_priority_consequence`. These
   calls are quality-clean and count towards TMB.

CNVs (`analysis-CNVs/filter_CNVs.R`):
- Segments at the sample's expected copy number (ploidy- and sex-aware) are dropped.
- `artifact_zero_cov`: cn = 0 while coverage remains.
- `low_probe_small_segment`.
- `FLAGS_all_genes`.

SVs (`analysis-SV/filter_SVs.R`):
- `FLAGS_all_genes`.
- `acmg_benign`.
- `ig_tcr_locus`.
- `recurrent_artifact`.
- `single_caller_TRA` (DELLY-only translocations).

Ig/TCR-partnered translocations and `sv_driver_loci` are exempt from the last two.

The header of each script documents every label.

## Layout

```
config.R                 read_wes_config(): template defaults + experiment config
load_WES_results.R       --load entry point; sources the three filters
build_sample_annotations.R
analysis-SNVs/           filter_SNVs.R, strandbias_rescue.R, snv_consequence_helpers.R, plot_snv_filters.R
analysis-CNVs/           filter_CNVs.R, plot_cn0_rel_depth.R, plot_cnv_tracks.R
analysis-SV/             filter_SVs.R
template/                runall, config.yaml (defaults), flags_genes.csv, SNV_consequence_terms.tsv, README.md
resources/               cytoBand.hg38.txt.gz
new_experiment           creates results/<date>_WES-filtering[-topic]/ from template/
```

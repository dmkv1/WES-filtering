# WES-filtering

WES-filtering of the WES-snakemake run in `results_dir` (`config.yaml`), with the
src/WES-filtering version in `WES-filtering.version`.

```bash
./runall --sample-annotations --load --plot 2>&1 | tee runall.log
```

`metadata_dir/metadata.csv` is regenerated from this experiment (`metadata.ver`).

## Config

Settings changed from the template defaults, and why.

## Results

QC-PASS / excluded sample counts, PASS SNV/CNV/SV counts, rescue counts, observations.

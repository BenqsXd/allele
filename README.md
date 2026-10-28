# DNA Variant Analyzer

Live: https://benqsxd.github.io/allele/

Offline, single-page analyzer for consumer raw DNA files (MyHeritage, 23andMe, AncestryDNA, FTDNA).

## Use
Open `index.html` in a browser (double-click), choose the raw `.csv`/`.txt` file. Nothing leaves your computer.
The report can be downloaded as a standalone HTML file or printed / saved as PDF.

## Files
- `index.html` – page shell and styles
- `app.js` – parsing, scoring, comparison and the interface (Overview, Categories, Gene scan, Traits, Compare)
- `db.js` – curated markers: effect alleles, weights, evidence level, literature references
- `clinvar.js` – generated from NCBI ClinVar; used for the full gene scan
- `build_clinvar.py` – regenerates `clinvar.js` (`python build_clinvar.py --refresh` to pull the latest ClinVar)
- `test.js` – self-check (`node test.js`)

`data/variant_summary.txt.gz` (~450 MB) is only needed to rebuild `clinvar.js` and can be deleted.

## Scoring
Each marker: +1 best genotype, −1 worst, heterozygote in between (additive / dominant / recessive model).
Weight = effect size (1–3) × evidence (strong 1, moderate 0.75, limited 0.5).
Category % = 50 + 50 × weighted mean → 100% best possible combination, 50% balanced, 0% worst possible.

Not a medical test. Low-pass / imputed data makes rare-variant calls unreliable.

# Allele — DNA Variant Analyzer

Version 1.0.2 · [Download](https://github.com/BenqsXd/allele/releases/latest)

Live: https://benqsxd.github.io/allele/

Offline, single-page analyzer for consumer raw DNA files (MyHeritage, 23andMe, AncestryDNA, FTDNA).

## Use
Open `index.html` in a browser (double-click), choose the raw `.csv`/`.txt` file. Nothing leaves your computer.
The report can be downloaded as a standalone HTML file or printed / saved as PDF.

## Read check (CRAM)
MyHeritage also offers the sequencing reads as a `.cram` file. Load the person's raw DNA file, then use **Check against reads** (or drop the `.cram` on the page). The app scans the reads in the browser (about 10–15 s for a 1.5 GB file) and shows, for every ClinVar hit and curated marker, whether the reads contain the variant. The CRAM must be aligned to GRCh37, as MyHeritage's are. Reference bases for the checked sites are built in, so no reference genome download is needed.

## Privacy
- DNA files are read with the browser's File API and analysed in the page. No code sends them anywhere.
- A Content-Security-Policy (`connect-src 'none'`, `form-action 'none'`) makes the browser block any upload, fetch, beacon or socket from the page, so this holds even if the code changed by mistake.
- Fonts and scripts load only from this site; there are no third-party requests, analytics or cookies.
- Once the page has loaded, it works without a connection. Links to dbSNP, ClinVar, PubMed etc. open only when clicked and carry only the rsID.

## Files
- `index.html` – page shell and styles
- `app.js` – parsing, scoring, comparison and the interface (Overview, Categories, Gene scan, Traits, Compare)
- `db.js` – curated markers: effect alleles, weights, evidence level, literature references
- `clinvar.js` – generated from NCBI ClinVar; used for the full gene scan
- `build_clinvar.py` – regenerates `clinvar.js` (`python build_clinvar.py --refresh` to pull the latest ClinVar)
- `vendor/cram.js` – CRAM reader ([@gmod/cram](https://github.com/GMOD/cram-js), MIT), loaded only for a read check
- `fonts/` – Geist and Geist Mono (SIL Open Font License, `fonts/OFL.txt`)
- `test.js` – self-check (`node test.js`)

`data/variant_summary.txt.gz` (~450 MB) is only needed to rebuild `clinvar.js` and can be deleted.

## Scoring
Each marker: +1 best genotype, −1 worst, heterozygote in between (additive / dominant / recessive model).
Weight = effect size (1–3) × evidence (strong 1, moderate 0.75, limited 0.5).
Category % = 50 + 50 × weighted mean → 100% best possible combination, 50% balanced, 0% worst possible.

Not a medical test. Low-pass / imputed data makes rare-variant calls unreliable.

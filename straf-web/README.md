# straf-web

Pure-browser TypeScript reimplementation of [STRAF](https://straf.fr).
All computation runs client-side; uploaded genotype files never leave the browser.

## Status

Proof-of-concept covering the most-used parts of the original Shiny app:

- TSV parsing (diploid + haploid, missing as `0`)
- Allele frequency table (per population)
- Forensic / popgen indices: **N, Nall, GD (Hexp), PIC, PM, PD, Hobs, PE, TPI**
- Combined PM, PD, PE across loci
- Hardy–Weinberg test (asymptotic χ²; exact MCMC test from Genepop is not yet ported)
- PCA on the individual × allele dosage matrix (Jacobi eigendecomp on the smaller Gram)
- Plotly interactive plots (basic dist, ~1 MB) and TSV downloads

Not yet ported:

- Pairwise Fst between populations (hierfstat::pairwise.WCfst)
- Linkage disequilibrium tests (Genepop MCMC)
- MDS + tree
- Reference-frequency comparison (STRidER)
- File conversions (Genepop, Familias, Arlequin, FASTA)
- Haploid haplotype statistics (haplotype diversity, pairwise distance)

## Run

```sh
npm install
npm run dev      # http://localhost:5173
npm run build    # static bundle in dist/
npm test         # numerical parity tests
```

## Bundle size

Targets a static deployable site. Production build (gzipped, expected):

- App code: ~10 KB
- Plotly basic dist: ~280 KB

If the Plotly bundle is still too large for the use case, swap to uPlot or
hand-rolled SVG (the only plot types used are bar charts and 2D scatter).

## Architecture

```
src/
  parser.ts             TSV → Genotypes (allele strings preserved verbatim)
  stats/
    freq.ts             allele frequency table + per-locus frequency maps
    forensic.ts         N, Nall, GD, PIC, PM, PD, Hobs, PE, TPI
    hwe.ts              χ² HWE test + chi-square CDF (Lanczos log-Γ)
    pca.ts              dosage matrix + Jacobi eigendecomposition
  ui/
    table.ts            HTML table rendering + TSV download
    plots.ts            Plotly wrappers
  main.ts               UI orchestration (no framework)
  style.css
```

## Numerical reference

`src/stats/forensic.test.ts` runs a small parity check against the example
dataset shipped with the original app. Spot-check values against R/STRAF for
your own data before relying on this for any forensic reporting.

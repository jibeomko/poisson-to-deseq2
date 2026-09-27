# poisson-to-deseq2

DESeq2 is easy to run, but I realized I couldn't really explain what happens between the count matrix and the `padj` column. So I went through it one step at a time: why a Poisson model isn't enough, how the dispersion gets estimated and shrunk, where the standard error comes from, and which genes `padj` is actually adjusted over. For each step I redid the numbers by hand or by simulation and compared them with what DESeq2 gives.

The notes are in Korean. I followed a Korean DESeq2 study text chapter by chapter (the text itself isn't in this repo) and wrote down the places where the text and the package don't agree. Steps 01 to 06 also have a notebook or a script, if you just want to run the numbers.

## Notes and code

| # | Note | What it's about | Code |
|---|---|---|---|
| 00 | [Overview](notes/00_overview.md) | What DESeq2 compares, and the steps from counts to padj | |
| 01 | [Poisson and overdispersion](notes/01_poisson_simulation.md) | Why a Poisson model can't describe the spread between biological replicates | [notebook (Python)](01_poisson_simulation.ipynb) |
| 02 | [Negative binomial and size factors](notes/02_negative_binomial.md) | What the dispersion α and the size factors each take care of | [notebook (R)](02_negative_binomial.ipynb) |
| 03 | [Dispersion estimation](notes/03_dispersion_estimation.md) | Estimating α from three replicates: likelihood, Cox-Reid, trend and MAP | [script](03_dispersion_estimation.R) |
| 04 | [The negative binomial GLM](notes/04_glm_condition_batch.md) | How the design matrix and IRLS give the log2 fold change and its SE | [script](04_glm_condition_batch.R) |
| 05 | [Wald test and LRT](notes/05_wald_vs_lrt.md) | Why the same fold change can get very different p-values, and what the LRT tests | [script](05_wald_vs_lrt.R) |
| 06 | [Multiple testing](notes/06_multiple_testing.md) | What padj is adjusted over, and how to test "down under starvation, back up once glucose returns" | [script](06_multiple_testing.R) |
| 07 | [LFC shrinkage and QC](notes/07_lfc_shrinkage_and_qc.md) | What `lfcShrink()` changes, and where the NAs in a results table come from | |
| 08 | [One gene, end to end](notes/08_one_gene_end_to_end.md) | Everything above for a single gene, by hand and inside DESeq2 | |

The notes do the explaining. The notebooks and scripts only reproduce the main numbers, and they stop with an error if anything comes out different from what the note says. [notes/README.md](notes/README.md) has a short guide to the notes.

## One gene, all the way through

All the notes follow the same gene: three control replicates and three glucose-starved ones, with counts 100, 130, 90 and 200, 250, 180 and size factors fixed at 1. Note 08 takes it through every step.

| Step | Value |
|---|---|
| Dispersion, plain negative binomial maximum likelihood | 0.014786 |
| Dispersion, Cox-Reid adjusted | 0.025385 |
| Dispersion, MAP with the study text's teaching prior | 0.053147 |
| log2 fold change | 0.977280 |
| Standard error (log2 scale) | 0.289057 |
| Wald statistic | 3.3809 |
| p-value | 0.00072244 |

![Three objective functions for the example gene's dispersion](figures/03_gene_a_objectives.png)

One thing to keep in mind: the prior in the MAP step (centred at 0.08, SD 0.7 on the log scale) is just a value the study text picks for teaching. If you put the same counts through `DESeq()` together with 2,000 simulated background genes (size factors still fixed at 1), the trend and the prior width come from those genes instead, and the final dispersion ends up around 0.096 (p ≈ 0.0097). If you fix the dispersion at the MAP value above (0.053147342 before rounding) and rerun only the Wald test, you get the hand-computed numbers back: SE 0.289057, Wald statistic 3.38092, p-value 0.000722.

## How I checked things

Everything ran on R 4.5.2 with DESeq2 1.50.2 (Bioconductor 3.22), plus NumPy and SciPy for a few simulations in note 01. Every code block in the notes was actually run, and the output under it is exactly what it printed. Whenever a note says DESeq2 does something a certain way, I checked it against the package source with `deparse()` and `args()`, and read the C++ parts in the 1.50.2 source tarball. I didn't install apeglm, ashr, glmGamPoi, stageR or tximport, so anything about those comes from their source code and help pages, and the notes say so.

Later I rebuilt the environment from `environment.yml` and reran the notebooks, the scripts and all the notes except 07. The results were the same, apart from floating-point noise around 1e-15, the sizes of the files note 08 saves, and small rendering differences in some figures.

A few things I didn't expect (each note has the full list at the end):

- Dispersion outliers are flagged when the log gene-wise estimate is more than 2 × √`varLogDispEsts` above the log trend. That width is a robust spread of the residuals, not the prior SD.
- For `lfcShrink(type = "normal")`, `lfcSE` is not a posterior SD. That description only holds for apeglm and ashr.
- A p-value becomes NA for all-zero genes, for Cook's outliers, and for genes whose observation weights make the coefficients impossible to estimate. A gene whose coefficients didn't converge still gets an ordinary-looking p-value.

## Running the code

```bash
conda env create -f environment.yml
conda activate poisson-to-deseq2
Rscript 03_dispersion_estimation.R     # same for 04, 05 and 06
jupyter lab                            # 01 uses the Python kernel, 02 the R kernel
```

Each script takes a few seconds and ends with "All checks passed." Activate the environment before starting Jupyter, because the R kernel starts whichever `R` comes first on your path. In this environment, loading DESeq2 or ggplot2 prints a warning that a package was built under R 4.5.3. It doesn't change any of the results.

Inside a note, the code blocks are meant to be run from top to bottom in one R session unless the note says otherwise. The figure code writes to `figures/` and should be run from the `notes/` folder.

## Layout

```
.
├── README.md
├── environment.yml
├── 01_poisson_simulation.ipynb   Python notebook
├── 02_negative_binomial.ipynb    R notebook
├── 03_dispersion_estimation.R
├── 04_glm_condition_batch.R
├── 05_wald_vs_lrt.R
├── 06_multiple_testing.R
├── notes/     study notes 00–08 (Korean) and their index
├── figures/   figures used in the notes
└── data/      empty for now
```

## References

- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 2014;15:550. https://doi.org/10.1186/s13059-014-0550-8
- DESeq2 vignette: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B* 1995;57:289–300.
- Bourgon R, Gentleman R, Huber W. Independent filtering increases detection power for high-throughput experiments. *PNAS* 2010;107:9546–9551. https://doi.org/10.1073/pnas.0914005107
- Van den Berge K, Soneson C, Robinson MD, Clement L. stageR: a general stage-wise method for controlling the gene-level false discovery rate in differential expression and differential transcript usage. *Genome Biology* 2017;18:151. https://doi.org/10.1186/s13059-017-1277-0
- Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 2019;35:2084–2092. https://doi.org/10.1093/bioinformatics/bty895

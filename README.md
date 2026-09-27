# poisson-to-deseq2

Working through what DESeq2 does between a count matrix and a column of adjusted p-values, one step at a time: why a Poisson model is not enough, how the dispersion is estimated and shrunk, where the standard error comes from, and which genes the adjusted p-value is adjusted over. At each step the numbers are recomputed by hand or by simulation and checked against DESeq2's own output and source code.

The notes are in Korean. They follow a Korean DESeq2 study text chapter by chapter (the text itself is not included here) and record where that text and the actual package behaviour differ. Steps 01 to 06 also have a notebook or script that reproduces the note's main numbers on its own.

## Notes and code

| # | Note | Question | Code |
|---|---|---|---|
| 00 | [Overview](notes/00_overview.md) | What does DESeq2 compare, and which steps lead from counts to padj? | |
| 01 | [Poisson and overdispersion](notes/01_poisson_simulation.md) | Why doesn't a Poisson model describe the spread between biological replicates? | [notebook (Python)](01_poisson_simulation.ipynb) |
| 02 | [Negative binomial and size factors](notes/02_negative_binomial.md) | What do the dispersion α and the size factor each account for? | [notebook (R)](02_negative_binomial.ipynb) |
| 03 | [Dispersion estimation](notes/03_dispersion_estimation.md) | How is α estimated from three replicates, and what do the Cox-Reid adjustment and MAP shrinkage add? | [script](03_dispersion_estimation.R) |
| 04 | [The negative binomial GLM](notes/04_glm_condition_batch.md) | How do the design matrix and IRLS produce the log2 fold change and its standard error? | [script](04_glm_condition_batch.R) |
| 05 | [Wald test and LRT](notes/05_wald_vs_lrt.md) | Why can the same fold change get very different p-values, and what does the LRT test instead? | [script](05_wald_vs_lrt.R) |
| 06 | [Multiple testing](notes/06_multiple_testing.md) | Which genes is padj adjusted over, and how do you test "down under starvation, back up once glucose returns"? | [script](06_multiple_testing.R) |
| 07 | [LFC shrinkage and QC](notes/07_lfc_shrinkage_and_qc.md) | What does `lfcShrink()` change, and where do the NAs in a results table come from? | |
| 08 | [One gene, end to end](notes/08_one_gene_end_to_end.md) | All of the above for a single gene, by hand and inside DESeq2 | |

The notes explain each step. The notebooks and scripts only reproduce the numbers, and each stops with an error if a result differs from the one printed in its note. [notes/README.md](notes/README.md) explains how to read the notes.

## One gene, all the way through

The running example is a gene with three control and three glucose-starved replicates: counts 100, 130, 90 and 200, 250, 180, with size factors fixed at 1. Note 08 carries it through every step.

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

The prior in the MAP step, centred at 0.08 with a standard deviation of 0.7 on the log scale, is a value the study text fixes for teaching. When the same counts go through `DESeq()` alongside 2,000 simulated background genes (size factors still fixed at 1), the dispersion trend and the prior width are estimated from those genes instead, and the final dispersion comes out at about 0.096 (p ≈ 0.0097). Fixing the dispersion at the MAP value above (0.053147342 before rounding) and rerunning only the Wald test gives back the hand-computed standard error 0.289057, Wald statistic 3.38092, and p-value 0.000722.

## How the checks were done

Everything was run with R 4.5.2 and DESeq2 1.50.2 (Bioconductor 3.22). A few simulations in note 01 use Python with NumPy and SciPy. Every code block in the notes was executed, and the output under it is what the run printed. Statements about how DESeq2 behaves were checked against the installed package's source with `deparse()` and `args()`, and the C++ routines against the 1.50.2 source tarball. apeglm, ashr, glmGamPoi, stageR, and tximport were not installed, so the claims about them rest on source code and help pages only; the notes mark these as not run.

The notebooks, the scripts, and notes 00–06 and 08 were run again in a fresh environment built from `environment.yml`. They gave the same results, apart from floating-point noise around 1e-15, the sizes of files that note 08 writes, and small rendering differences in some figure files.

Each note ends with a table of the places where the study text and DESeq2 1.50.2 disagree, or where the text leaves out an implementation detail. Three examples:

- Dispersion outliers are flagged when the log gene-wise estimate exceeds the log trend by more than 2 × √`varLogDispEsts`, a robust spread of the residuals, not the prior SD.
- For `lfcShrink(type = "normal")`, `lfcSE` is not a posterior SD. That description holds only for apeglm and ashr.
- A p-value becomes NA for all-zero genes, for Cook's outliers, and for genes whose observation weights make the coefficients impossible to estimate. A gene whose coefficients did not converge keeps an ordinary-looking p-value.

## Running the code

```bash
conda env create -f environment.yml
conda activate poisson-to-deseq2
Rscript 03_dispersion_estimation.R     # likewise 04, 05 and 06
jupyter lab                            # 01 runs on the Python kernel, 02 on the R kernel
```

Each script finishes in under ten seconds and ends with "All checks passed." Activate the environment before starting Jupyter, because the R kernel starts whichever `R` is first on the path. In this environment, loading DESeq2 or ggplot2 prints a warning that a package was built under R 4.5.3. It does not change any result.

Inside a note, code blocks are meant to be run top to bottom in one R session unless the note says otherwise. Figure code writes to `figures/` and should be run from the `notes/` directory.

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

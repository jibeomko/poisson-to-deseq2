# DESeq2 study notes

When you run DESeq2, the results table gives you a log2 fold change, lfcSE, pvalue and padj. I wanted to know how these values are computed from the raw counts and which statistical steps they go through. So I followed a personal study text (*DESeq2 Study Notes: From Counts to Conclusions*) chapter by chapter and computed each step myself; this folder is the record. The explanations and numbers in the study text were checked by running them again in an installed DESeq2.
I wrote with a reader in mind who knows means, variances and p-values but has only heard the names of the likelihood and the GLM.

These are English translations of the Korean notes in the [parent folder](../README.md). The code blocks are the same code: Korean comments, Korean labels inside strings and the matching lines of printed output were translated, and no number was changed. Run the code from the `notes/` folder, as in the originals, so that the figure paths (`../figures/`) resolve.

## Reading order

If this is your first time, start with 00_overview.md and read in order. Each note opens with one question, shows the numbers first, and gives the formulas after that. Source-code digging and detailed checks are folded into "Going deeper" at the end of each note, so you can skip them on a first read.

All the notes use the same example: gene A from chapter 16 of the study text, counted in cells grown three times each in ordinary medium (Ctrl) and in medium without glucose (Starvation). All size factors are set to 1. The example names differ from the study text's, and quotations from the study text were changed to match them.

| Condition | count | Mean |
|---|---|---:|
| Ctrl | 100, 130, 90 | 106.667 |
| Starvation | 200, 250, 180 | 210 |

Along the notes, these six numbers change like this.

| Step | Gene A | Note |
|---|---|---|
| dispersion: plain NB maximum likelihood → Cox-Reid adjustment | 0.014786 → 0.025385 | [03](03_dispersion_estimation.md) |
| final dispersion: MAP with the teaching prior set by the study text | 0.053147 | [03](03_dispersion_estimation.md) |
| log2 fold change and standard error | 0.977280, 0.289057 | [04](04_glm_condition_batch.md) |
| Wald statistic and p-value | 3.3809, 0.00072244 | [05](05_wald_vs_lrt.md) |

[08](08_one_gene_end_to_end.md) redoes this calculation from beginning to end in one go and compares it with the same counts put into an actual `DESeq()` run. It also has an example that adds a third condition, starved cells given glucose again (Starvation+Glucose, called Glucose for short in the notes).

## List of notes

| Note | What I wanted to know | Study text | Code |
|---|---|---|---|
| [00. What does DESeq2 ask of six counts?](00_overview.md) | What DESeq2 compares, and which steps lie between the counts and padj | Reading guide, chapter 1, appendices C and D | |
| [01. Why do counts differ so much between replicates of the same condition?](01_poisson_simulation.md) | Why a Poisson model cannot describe the spread between biological replicates | 2.1–2.2 | [notebook (Python)](../../01_poisson_simulation.ipynb) |
| [02. How much do counts spread, and how is the difference in sequencing depth corrected?](02_negative_binomial.md) | What the negative binomial dispersion α and the size factors each do | 2.3–2.4, chapter 3 | [notebook (R)](../../02_negative_binomial.ipynb) |
| [03. How is the dispersion α set from three replicates?](03_dispersion_estimation.md) | Setting α with the likelihood, the Cox-Reid adjustment, the trend and the prior | Chapters 5–7 | [R script](../../03_dispersion_estimation.R) |
| [04. How do the LFC and SE come out of a model with condition and batch?](04_glm_condition_batch.md) | The design matrix, IRLS, standard errors, and the covariance of a comparison (contrast) | Chapters 4 and 8 | [R script](../../04_glm_condition_batch.R) |
| [05. Why do the p-values differ when the fold change is the same?](05_wald_vs_lrt.md) | What the Wald test and the LRT each test | Chapters 9, 10 and 12 | [R script](../../05_wald_vs_lrt.R) |
| [06. How should p-values be read when thousands of genes are tested at once?](06_multiple_testing.md) | BH adjustment, independent filtering, and adjusting several comparisons and rescue claims | Chapters 11 and 14 | [R script](../../06_multiple_testing.R) |
| [07. Can a large fold change from small counts be trusted?](07_lfc_shrinkage_and_qc.md) | LFC shrinkage, NAs in the results table, QC plots | Chapter 13 | |
| [08. How do six counts become one p-value?](08_one_gene_end_to_end.md) | Computing gene A to the end by hand and comparing with actual DESeq2 | Chapters 15 and 16 | |

The notebooks and scripts in the code column reproduce only the main numbers of the notes; the explanations are in the notes. When run, they compare their results with the numbers written in the notes and stop with an error if any of them differ.

## How things were checked

Everything ran on R 4.5.2 and DESeq2 1.50.2 (Bioconductor 3.22), and some of the simulations in 01 were done in Python (numpy, scipy). Every code block in the notes was actually run, and the output right below it is what came out at the time (translated as described above). Within a note the code was run from the top in one session, and where a note splits sessions it says so. Run from the `notes/` folder, the figure code saves PNGs to `../figures/`. The `<bytecode: 0x...>` addresses printed with functions differ from run to run, so you can ignore them.

Later I built a fresh environment from the repository's `environment.yml` and reran the notes except 07 and the code files. Apart from floating-point differences around 1e-15, the sizes of the files 08 saves and rendering differences in some figures, the results were the same.

Statements about how DESeq2 behaves were checked by pulling the installed package source out with `deparse()` or `args()`. The parts written in C++, such as `fitDisp` and `fitBeta`, exist only as compiled files in the installed package, so I read them in the 1.50.2 source tarball. apeglm, ashr, glmGamPoi, stageR and tximport were not installed, so what the notes say about them was written from their source and help pages only, and the notes mark it as not confirmed by running.

## Where the results differ from the study text

Most of the study text's explanations and numbers matched the results of running them, and all the numbers of chapter 16 are reproduced too. Here are the differences worth remembering.

- To match the Z, p and confidence interval of chapter 16 to the last digit, the MAP dispersion has to be converged fully by tightening the tolerance of `optimize()`. With the study text's code as it is, it stops at 0.053146542 ([08](08_one_gene_end_to_end.md)).
- The prior of chapter 16 (centre 0.08, SD 0.7) is a value the study text set for explanation. Putting the same counts into an actual `DESeq()` estimates the trend and the prior width from the other genes, so the final dispersion comes out not at 0.053 but at 0.093–0.098 depending on the background genes ([08](08_one_gene_end_to_end.md)).
- The study text does not give the criterion for dispersion outliers. In fact a gene is an outlier when its log gene-wise value is more than `2 × sqrt(varLogDispEsts)` above the log trend, and this width is a robust measure of the spread of the residuals, not the prior SD ([03](03_dispersion_estimation.md)).
- Section 6.5 of the study text has a gene that barely wobbles within conditions. The study text says only that this gene's α "can be small" under `~condition`, but running it shows it does not just get small; it hits DESeq2's lower limit of 1e-8 ([03](03_dispersion_estimation.md)).
- The statement "the mean is first estimated with an initial dispersion" holds only for some designs. For designs that are solved by group means alone, like `~condition`, it uses group means that do not depend on α ([00](00_overview.md)).
- The `lfcSE` returned by `lfcShrink(type = "normal")` is not a posterior SD. The posterior SD description fits only apeglm and ashr ([07](07_lfc_shrinkage_and_qc.md)).
- The "fitting problems" of appendix B 17 split into two cases. A β that does not converge does not make p NA; only rows whose coefficients cannot be estimated because of weights get NA even for p ([07](07_lfc_shrinkage_and_qc.md)).

Each note's "Going deeper" lists every place where it differs from the study text.

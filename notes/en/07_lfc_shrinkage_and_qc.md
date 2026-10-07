# 07. Can a large fold change from small counts be trusted? LFC shrinkage and checking the results table

When a gene with single-digit counts shows a fold change of more than 50-fold, I hesitate to trust that number as it is. In this note I pull such LFCs toward 0 with `lfcShrink()`, sort out where the NAs in the results table come from, and read the QC plots separately from the test. It is easier if you have seen the dispersion shrinkage of [03](03_dispersion_estimation.md), the Wald test of [05](05_wald_vs_lrt.md) and the padj of [06](06_multiple_testing.md) first.

> Study text chapter 13

## 1. A 56-fold difference from 3 vs 3: can it be trusted?

Let me start with gene A, the shared example. Its counts (the number of reads assigned to the gene in one sample) are 100, 130, 90 in Ctrl and 200, 250, 180 in Starvation. The log2 fold change (LFC) DESeq2 gives is 0.977; the LFC is the log2 of the ratio of two condition means, so this means about a 1.97-fold increase.

The SE (standard error) of this estimate is 0.289. The SE says how much the estimate would move if the same experiment were repeated. The 95% interval made as estimate ± 1.96 × SE is [0.41, 1.54], which as a fold is between 1.3-fold and 2.9-fold. With counts in the hundreds, the estimate is fairly solid.

Now a gene with far fewer counts. With DESeq2's built-in simulation function I made data with 3000 genes and 3 samples of condition A vs 3 of B, and picked gene1349.

```r
suppressPackageStartupMessages(library(DESeq2)); options(width = 110)
set.seed(7)
dds <- makeExampleDESeqDataSet(n = 3000, m = 6, betaSD = 1)   # simulated data: 3000 genes, 3 A vs 3 B
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds)
counts(dds)["gene1349", ]
signif(as.data.frame(res["gene1349", c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]), 4)
sel <- which(res$baseMean < 5 & abs(res$log2FoldChange) > 3 & res$lfcSE > 1)   # low expression but large LFC
oneZero <- rowSums(counts(dds)[sel, 1:3]) == 0 | rowSums(counts(dds)[sel, 4:6]) == 0
c(n_sel = length(sel), one_group_all_zero = sum(oneZero))
ok <- !is.na(res$log2FoldChange)                                               # overall distribution: most are small
cat("genes with LFC:", sum(ok), "| share with |LFC| < 1:", round(mean(abs(res$log2FoldChange[ok]) < 1), 3), "\n")
```

```
sample1 sample2 sample3 sample4 sample5 sample6 
     10       9      11       0       0       0 
         baseMean log2FoldChange lfcSE    pvalue     padj
gene1349    4.886         -5.808 1.699 0.0006278 0.008435
             n_sel one_group_all_zero 
                81                 38 
genes with LFC: 2984 | share with |LFC| < 1: 0.539
```

A has 10, 9, 11 while B has 0, 0, 0. DESeq2 gives an LFC of −5.81, that is, B 56 times lower than A. p = 0.00063, and padj (the p-value adjusted for testing many genes together, [06](06_multiple_testing.md)) is also significant at 0.0084. But the SE is as large as 1.70. The 95% interval is [−9.14, −2.48], which as a fold runs from about a 5.6-fold decrease to about a 560-fold decrease.

In fact, what this data says is roughly "B is much lower than A". If all three B samples are 0, a B mean of 0.1 or of 0.01 is about equally plausible. Yet the LFC has to write exactly that difference as a single number.

Looking at the MLE makes the situation clearer. The likelihood, with the observations held fixed, says how plausible a candidate parameter makes those observations. The MLE (maximum likelihood estimate) is the parameter value with the largest likelihood (covered in detail in [03](03_dispersion_estimation.md)). If B is all 0, the likelihood keeps increasing as the LFC is pushed more negative, so there is no finite MLE. −5.81 is just where DESeq2 stopped because it puts a lower bound of 0.5 (`minmu`) on the expected count (the mean count the model expects in that sample) when it computes the coefficients. Lowering the bound to 0.1 gives −8.13, and to 0.01 gives −11.45 ("When one group is all 0" under Going deeper).

Nor is a group of all 0s the only problem. With small counts, chance fluctuations show up strongly in ratios. For example, a Poisson count with mean 10 has a standard deviation of about 3.2, around 30% of the mean, but with mean 1000 the standard deviation is about 32, 3% of the mean ([01](01_poisson_simulation.md)).

The latter part of the output above describes the simulated data as a whole. baseMean is, for each gene, the mean of the normalized counts of the six samples. A normalized count is the count divided by the size factor, and the size factor is the multiplier that corrects for sequencing depth differing between samples. 81 genes had baseMean below 5 together with |LFC| > 3 and SE > 1, and 38 of them had one group all 0.

## 2. How does lfcShrink pull the LFC toward 0?

The last line of the output in section 1 shows that of the 2984 genes with a computed LFC, 54% have |LFC| below 1. More than half the genes change by less than 2-fold. It is then better to read the extreme LFCs of genes with little information as pulled toward this "overall tendency" and to leave genes with lots of information as they are. Pulling estimates that carry little information toward the overall tendency like this is called shrinkage.

The overall tendency is expressed as a prior, a distribution for where a parameter is likely to be before the data are seen. `lfcShrink(type = "normal")` takes the prior of the LFC to be a normal distribution centred on 0, N(0, σ²), and sets σ² from the distribution of the MLE LFCs of all genes. Estimating the prior from the data as a whole like this is called empirical Bayes. What `type = "normal"` returns is the most plausible value when the likelihood and the prior are considered together, the MAP (maximum a posteriori).

The dispersion had a similar shrinkage ([03](03_dispersion_estimation.md)). The dispersion α expresses how much more the counts spread than the Poisson in variance = μ + αμ² ([02](02_negative_binomial.md)). The two shrinkages differ in the direction they pull and in where they are used.

| | dispersion shrinkage | LFC shrinkage |
|---|---|---|
| What is pulled | the dispersion α | the LFC |
| Where to (prior centre) | the trend: a curve showing where the dispersion roughly lies as a function of mean expression (not 0) | 0 (normal, apeglm) or a mixture near 0 (ashr) |
| Why | To reduce the instability of per-gene dispersion estimates | To shrink the noisy effect sizes of low-expression, high-variance genes |
| When | Automatically inside `DESeq()` (`estimateDispersionsMAP`) | `lfcShrink()` called separately after `DESeq()` finishes |
| Where the result goes | The whole subsequent GLM fit and test (reflected in the test's SE) | Ranking, plotting and reporting of effect sizes. In the default call only two columns change, `log2FoldChange` and `lfcSE` |

The GLM (generalized linear model) in the table is the model that writes the mean of the counts on the log scale as a sum of condition, batch and so on ([04](04_glm_condition_batch.md)).

The expression that `type = "normal"` maximizes looks like this.

$$
\ell(\beta) - \frac{\beta^2}{2\sigma^2}
$$

β is the LFC being sought (in log2 units), and ℓ(β) is the log likelihood, which says how well β explains this gene's counts. The MLE is the value that maximizes ℓ(β) alone. The β²/(2σ²) behind it is a penalty that grows as β moves away from 0, and it is stronger the smaller σ² is.

For a gene with many counts, ℓ(β) is sharp near the MLE. Moving even slightly off the MLE drops the likelihood a lot, so the penalty can hardly move β. Conversely, for a gene with few counts ℓ(β) is flat, so many β are about equally plausible; the penalty wins and pulls it a long way toward 0.

Approximating the likelihood by a normal distribution, this effect can be written in one line.

$$
\hat\beta_{\text{shrunken}} \approx \hat\beta_{\text{MLE}} \times \frac{\sigma^2}{\sigma^2 + \text{SE}^2}
$$

β̂_MLE and SE are the LFC and its standard error given by `results()`. σ²/(σ² + SE²) is a factor between 0 and 1: close to 1 if the SE is much smaller than σ, and close to 0 if the SE is large.

This formula is only an approximation to aid understanding; DESeq2 finds the MAP directly with a penalized IRLS (the iteratively reweighted least squares calculation that finds the coefficients; ridge-penalized IRLS). Let me borrow σ² = 1.591, which DESeq2 set for the simulated data above, and plug in gene A. SE² = 0.289² = 0.0836, so the factor is 1.591 / (1.591 + 0.0836) = 0.950, and 0.977 × 0.950 = 0.928. Computing it directly with DESeq2 gives the same value. In the code below, the first value of `betaPriorVar`, 1e6, means there is effectively no prior on the intercept.

```r
suppressPackageStartupMessages(library(DESeq2))
cts <- matrix(c(100L, 130L, 90L, 200L, 250L, 180L), nrow = 1, dimnames = list("geneA", paste0("s", 1:6)))
cd  <- data.frame(condition = factor(rep(c("Ctrl", "Starvation"), each = 3)))
dA  <- DESeqDataSetFromMatrix(cts, cd, ~ condition)
sizeFactors(dA) <- rep(1, 6)
dispersions(dA) <- 0.053147                                  # fix the final α
mle <- results(nbinomWaldTest(dA, quiet = TRUE))
s2  <- 1.591035                                              # prior variance estimated by lfcShrink on the section 1 data (log2)
map <- results(nbinomWaldTest(dA, betaPrior = TRUE, betaPriorVar = c(1e6, s2),
                                modelMatrixType = "standard", quiet = TRUE))
se  <- mle$lfcSE
cat(sprintf("MLE  LFC %.6f  SE %.6f\n", mle$log2FoldChange, se))
cat(sprintf("MAP  LFC %.6f  (prior N(0, %.3f))\n", map$log2FoldChange, s2))
cat(sprintf("approx s2/(s2+SE^2) = %.4f  ->  %.4f\n", s2 / (s2 + se^2), mle$log2FoldChange * s2 / (s2 + se^2)))
```

```
MLE  LFC 0.977280  SE 0.289057
MAP  LFC 0.928519  (prior N(0, 1.591))
approx s2/(s2+SE^2) = 0.9501  ->  0.9285
```

Gene A shrank only a little, from 0.977 to 0.929 (about 1.97-fold → 1.90-fold). Let me do the same calculation for gene1349 of section 1 (run continuing from the code of section 1).

```r
shr <- lfcShrink(dds, coef = "condition_B_vs_A", type = "normal", quiet = TRUE)
s2 <- priorInfo(shr)$betaPriorVar[["conditionB"]]      # prior variance σ² estimated from all genes (log2 units)
s2
g <- c("gene1349", "gene2613", "gene2055")
mle <- res[g, "log2FoldChange"]; se <- res[g, "lfcSE"]
round(data.frame(baseMean = res[g, "baseMean"], MLE = mle, SE = se, factor = s2 / (s2 + se^2),
                 approx = mle * s2 / (s2 + se^2), shrunken = shr[g, "log2FoldChange"], row.names = g), 3)
```

```
[1] 1.591035
         baseMean    MLE    SE factor approx shrunken
gene1349    4.886 -5.808 1.699  0.355 -2.064   -1.969
gene2613    3.950  5.457 3.754  0.101  0.554    0.352
gene2055 1205.552  2.063 0.313  0.942  1.943    1.943
```

The `factor` column is the factor, `approx` the value of the approximation, and `shrunken` the value DESeq2 gave. gene1349 has a factor of 0.355, so −5.81 shrank a lot, to −1.97 (about a 3.9-fold decrease). gene2613, with an SE of 3.75, was squashed from 5.46 to 0.35, while gene2055, with counts above 1000, barely moved, from 2.06 to 1.94. The approximation is nearly exact for genes with lots of information, and for genes with little information it gets only the direction and the rough size right.

An MA plot puts the mean count on the x axis and the LFC on the y axis. Putting before and after shrinkage side by side shows the effect at a glance.

![MA plot before and after LFC shrinkage](../../figures/07_ma_before_after_shrink.png)

The left panel is the MLE of `results()`, the right the values of `lfcShrink(type = "normal")`. The points spread widely at the left end, where the mean count is small, gather near 0 on the right, while the side with large means is almost unchanged.

<details>
<summary>Code for the figure</summary>

```r
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
set.seed(7)
dds <- makeExampleDESeqDataSet(n = 3000, m = 6, betaSD = 1)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds)
shr <- lfcShrink(dds, coef = "condition_B_vs_A", type = "normal", quiet = TRUE)
panel <- c("MLE: results()", "Shrunken: lfcShrink(type = \"normal\")")
df <- rbind(data.frame(panel = panel[1], gene = rownames(res), mean = res$baseMean, lfc = res$log2FoldChange),
            data.frame(panel = panel[2], gene = rownames(shr), mean = shr$baseMean, lfc = shr$log2FoldChange))
df <- df[df$mean > 0, ]                      # all-zero genes have nothing to plot
df$panel <- factor(df$panel, panel)
lab <- df[df$gene %in% c("gene1349", "gene2613", "gene2055"), ]
lab$hj <- ifelse(lab$gene == "gene2055", 1.15, -0.15)   # label at the right edge goes on the left
p <- ggplot(df, aes(mean, lfc)) +
  geom_hline(yintercept = 0, colour = "#52514e", linewidth = 0.4) +
  geom_point(colour = "grey70", size = 0.6, alpha = 0.6) +
  geom_point(data = lab, colour = "#eb6834", size = 2) +
  geom_text(data = lab, aes(label = gene, hjust = hj), size = 3, vjust = -0.5) +
  scale_x_log10(breaks = 10^(-1:3), labels = c("0.1", "1", "10", "100", "1000")) + facet_wrap(~ panel) +
  labs(x = "Mean of normalized counts (log scale)", y = "log2 fold change (B vs A)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"))
out <- "../figures/07_ma_before_after_shrink.png"   # run from notes/
ggsave(out, p, width = 7, height = 4.5, dpi = 150, bg = "white")
cat("genes plotted per panel:", sum(df$panel == panel[1]), "| LFC range MLE:", round(range(res$log2FoldChange, na.rm = TRUE), 2),
    "| shrunken:", round(range(shr$log2FoldChange, na.rm = TRUE), 2), "\n")
cat("file exists:", file.exists(out), file.size(out), "bytes\n")
```

```
genes plotted per panel: 2984 | LFC range MLE: -6.04 6.02 | shrunken: -3.33 2.99 
file exists: TRUE 182810 bytes
```

</details>

In numbers: collecting only genes with |MLE LFC| > 0.5 and taking the median of |shrunken| / |MLE| per baseMean bin gave 0.270 for baseMean 5 or less, 0.616 for 5–20, 0.817 for 20–100, 0.895 for 100–1000 and 0.941 above 1000. The higher the expression, the closer to 1, that is, the less it was pulled. This calculation is in experiment 1. In this note, "experiment N" always refers to a collapsible block under Going deeper.

`lfcShrink()` has three methods (`type`). All three use the posterior, the distribution that combines the likelihood and the prior. The mode is the peak of that distribution, and the mean its average. The three differ in the shape of the prior and in the value they return.

| `type=` | Shape of the prior | Value returned | This environment |
|---|---|---|---|
| `"apeglm"` (default) | A heavy-tailed distribution (help page: adaptive Student's t) | posterior mode (MAP) | Not installed, cannot run |
| `"ashr"` | A mixture of several normal distributions centred on 0 | posterior mean (labelled "MMSE", meaning the estimate that minimizes the mean squared error) | Not installed, cannot run |
| `"normal"` | Normal distribution N(0, σ²) | posterior mode (MAP) | Built into DESeq2. Every run in this note |

So writing "lfcShrink results are all MAP" is wrong for ashr. And a heavy-tailed prior like apeglm's squashes large effects less.

apeglm and ashr are not installed in this environment. But the default of `type` is `"apeglm"`, so omitting `type` stops with an installation error. The study text's example in section 13.2, `lfcShrink(dds, coef="condition_Starvation_vs_Ctrl", type="apeglm")`, has the right coefficient name but stops for the same reason. So all the shrinkage in this note was run with `type = "normal"`, and what is said about apeglm and ashr was checked in the source and help pages (Going deeper).

## 3. What changes in the results table after shrinking?

Only two columns change, `log2FoldChange` and `lfcSE`; the p-values and padj stay the same (run continuing from the code of sections 1–2).

```r
c(pvalue_same = identical(res$pvalue, shr$pvalue), padj_same = identical(res$padj, shr$padj))
z <- shr["gene1349", "log2FoldChange"] / shr["gene1349", "lfcSE"]        # a calculation that should not be done
round(c(stat_MLE = res["gene1349", "stat"], z_from_shrunken = z), 3)
signif(c(p_saved = res["gene1349", "pvalue"], p_from_z = 2 * pnorm(-abs(z))), 3)
```

```
pvalue_same   padj_same 
       TRUE        TRUE 
       stat_MLE z_from_shrunken 
         -3.419          -3.142 
 p_saved p_from_z 
0.000628 0.001680
```

`pvalue` and `padj` are identical, down to the values. The test is still the Wald test done with the MLE and its SE ([05](05_wald_vs_lrt.md)), which compares (estimate − 0) / SE with the standard normal. The shrunken LFC is not for testing but for ranking, plotting and reporting effect sizes.

What, then, if the shrunken LFC is divided by the new `lfcSE` to make a new Z? At first I too thought that was fine. But gene1349's Wald statistic is −3.419, while the value made that way is −3.142, and the p from it is 0.00168, different from the stored p of 0.000628. This new number is neither an MLE test nor Bayesian inference, just a mix of the two.

This happens because the meaning of `lfcSE` changed. The column name is the same, but what it holds differs depending on the function that made the object.

| Function that made the object | Meaning of `lfcSE` | gene1349 |
|---|---|---|
| `results()` | Wald SE of the MLE | 1.70 |
| `lfcShrink(type = "apeglm")`, `"ashr"` | posterior SD (the standard deviation of the posterior distribution) | Not installed, not run |
| `lfcShrink(type = "normal")` | SE of the penalized estimate. Neither a posterior SD nor a Wald SE | 0.63 |

The help page says the `lfcSE` of normal is a posterior SD, but the implementation differs. Computing the posterior SD directly with a Laplace approximation gives 0.94 for gene1349, and the 0.63 that normal returns is smaller (experiment 1). The column description under normal also stays "standard error". apeglm and ashr rewrite this description as "posterior SD" in the source.

So the same calculation, "estimate ± 1.96 × lfcSE", builds different intervals.

| Bar in a figure or table | Calculation | Kind | gene1349 |
|---|---|---|---|
| Wald 95% interval | LFC ± 1.96·lfcSE of `results()` | Nominal (frequentist) interval of the MLE | [−9.14, −2.48] |
| apeglm posterior interval | `interval` of the apeglm fit object obtained with `returnList=TRUE` | Credible interval (an interval from the posterior) | `"interval"` is in the source's `ape.cols`. Not installed, not run |
| normal's MAP ± 1.96·lfcSE | `lfcShrink(type = "normal")` | Neither | [−3.20, −0.74] |
| (Reference) MAP ± 1.96·Laplace SD | Computed directly | Approximate posterior interval | [−3.81, −0.12] |
| padj | BH (Benjamini–Hochberg) adjustment | Not an interval but an FDR level | 0.0084 |

gene2613, with an SE of 3.75, has a Wald interval of [−1.90, 12.81], which includes 0, while normal's MAP ± 1.96·lfcSE is [−0.42, 1.12], narrow enough to make the uncertainty look smaller than it is. So the caption of a figure should say which of these the points and bars are. The points of `plotMA(res)` are MLEs, and those of `plotMA(shr)` shrunken LFCs.

Passing `svalue = TRUE` to apeglm or ashr gives s-values instead of p-values. The s-value is "the average probability that the sign of the LFC is wrong, among genes whose s-value is equal to or smaller than this gene's". padj, on the other hand, was "the smallest level that controls the expected proportion (FDR) of genes with no real difference (β = 0) mixed into the list when this gene is called a discovery" ([06](06_multiple_testing.md)). Both are list-level values, but the s-value deals with getting the sign wrong and padj with genes of no difference getting mixed in.

There are exceptions to the p-values staying the same, too. `lfcThreshold` is the argument that makes the test "|LFC| is larger than this value" instead of "the LFC is not 0"; giving normal `lfcThreshold > 0` changes `stat`, `pvalue` and `padj` as well, into the threshold test of the shrunken fit (experiment 2). Also, apeglm and ashr drop the `stat` column in the default call. The columns for each type are summarized in the lfcShrink source block under Going deeper.

## 4. How to shrink when comparing Starvation and Glucose directly?

Take the three groups of the shared example, Ctrl, Starvation and Starvation+Glucose (Glucose for short). With Ctrl as the reference there are two coefficients, `condition_Starvation_vs_Ctrl` and `condition_Glucose_vs_Ctrl`. "Glucose vs Starvation" is not a coefficient but the difference of two coefficients, a contrast. A contrast is a vector that sets which two conditions (or which combination of coefficients) to compare ([04](04_glm_condition_batch.md)).

Does shrinking the two coefficients separately and subtracting them give the shrunken LFC of Glucose vs Starvation, then? I checked with three-group simulated data (A=Ctrl, B=Starvation, C=Glucose, 3 each, 2000 genes) (experiment 2).

For the MLE the subtraction is exact. The maximum difference between (C − A) − (B − A) and `results(contrast = c("condition", "C", "B"))` is 8.9e-16. The exception is 14 genes where B and C are all 0; for these `results()` overwrites the contrast's LFC with 0 and its p with 1.

But the shrunken values differ. The difference between the two shrunken coefficients subtracted and `lfcShrink(contrast = ...)` was at most 0.69, with a median of 0.05. For example, the low-expression gene158 (baseMean 4.7) has MLE 4.27, subtraction 0.66 and contrast shrink 1.35.

The difference comes from the prior. The design matrix here is the table that records in numbers which condition, batch or pair each sample belongs to. The `coef` path sets a separate prior variance for each coefficient in the ordinary design matrix with a reference level (Starvation 0.80, Glucose 1.695). The `contrast` path refits with a different form of design matrix that gives every level the same prior variance (1.195). Re-setting Starvation as the reference (`relevel`) and calling with `coef` gives yet another set of prior variances (0.8, 1.089), differing from the contrast path by up to 1.24 and from the subtraction by up to 1.28. apeglm accepts only `coef`, and each call re-estimates the posterior of that one coefficient only.

In the end the three methods use different priors, so none of them alone is "the right answer". The study text's conclusion, that subtracting two shrunken coefficients is not the same as the proper posterior shrinkage of the contrast, is right. But for normal the reason is which coefficients the prior was put on, because even for the same comparison the prior differs depending on how the coefficients were defined (the reference level, the form of the design matrix). When reporting, state which path was used and record the prior variances with `priorInfo()`.

## 5. What do the NAs in the results table mean?

Back to gene A. Suppose the first Ctrl sample were 1000 instead of 100, a situation where an accident such as sample contamination made just one sample jump. The dispersion α was fixed at 0.053147.

```r
suppressPackageStartupMessages(library(DESeq2))
fitA <- function(ctrl1) {
  cts <- matrix(as.integer(c(ctrl1, 130, 90, 200, 250, 180)), nrow = 1, dimnames = list("geneA", paste0("s", 1:6)))
  d <- DESeqDataSetFromMatrix(cts, data.frame(condition = factor(rep(c("Ctrl", "Starvation"), each = 3))), ~ condition)
  sizeFactors(d) <- rep(1, 6); dispersions(d) <- 0.053147
  nbinomWaldTest(d, quiet = TRUE)
}
for (x in c(100, 1000)) {
  d <- fitA(x)
  r <- tryCatch(results(d), error = function(e) conditionMessage(e))
  cat("Ctrl1 =", x, "| cooks:", round(assays(d)[["cooks"]][1, ], 3), "| maxCooks:", round(mcols(d)$maxCooks, 3),
      "| cutoff:", round(qf(0.99, 2, 4), 2), "\n")
  if (is.character(r)) cat("  results error:", r, "\n") else cat("  LFC", round(r$log2FoldChange, 4), "pvalue", signif(r$pvalue, 4), "\n")
  cat("  robust alpha:", DESeq2:::robustMethodOfMomentsDisp(d, model.matrix(~condition, colData(d))), " H:", round(assays(d)[["H"]][1,], 4), "\n")
}
```

```
Ctrl1 = 100 | cooks: 0.03 0.363 0.185 0.019 0.304 0.171 | maxCooks: 0.363 | cutoff: 18 
  LFC 0.9773 pvalue 0.0007224 
  robust alpha: 0.04  H: 0.3333 0.3333 0.3333 0.3333 0.3333 0.3333 
Ctrl1 = 1000 | cooks: 18.801 4.088 5.355 0.019 0.304 0.171 | maxCooks: 18.801 | cutoff: 18 
  LFC -0.9535 pvalue NA 
  robust alpha: 0.04  H: 0.3333 0.3333 0.3333 0.3333 0.3333 0.3333
```

In the original data the maximum of `cooks` is 0.363. Changing the first sample to 1000 makes that sample's value jump to 18.80 and `pvalue` becomes NA. The LFC remains, at −0.95, a value with its direction flipped because the Ctrl mean jumped to 406.7.

`cooks` is Cook's distance, which expresses how much one sample's count drags the coefficient estimates. DESeq2 computes it like this.

$$
\text{Cook}_{ij} = \frac{(K_{ij} - \mu_{ij})^2}{V_{ij}} \cdot \frac{1}{p} \cdot \frac{h_{jj}}{(1 - h_{jj})^2}, \qquad V_{ij} = \mu_{ij} + \tilde\alpha_i\, \mu_{ij}^2
$$

K_ij is the count of gene i in sample j, and μ_ij the expected count of that sample. V_ij is the variance of that count, and the α̃_i in it is a dispersion computed separately for this calculation only, so that it is less swayed by jumping values (lower bound 0.04; gene A's is 0.04 too). The first factor, (K − μ)²/V, measures how far the count is from its expectation. p is not a p-value but the number of coefficients, so 2 for two groups. h_jj is the leverage, how strongly that sample is positioned to pull its own fitted value. With 3 per group it is 1/3, and h/(1 − h)² = 0.75.

I plugged in the case with 1000 directly. μ = (1000 + 130 + 90) / 3 = 406.67, V = 406.67 + 0.04 × 406.67² = 7021.8, and (1000 − 406.67)² / 7021.8 = 50.14; multiplying by 1/2 and 0.75 gives 18.80, the same as the output.

The cutoff is the 99% quantile of the F distribution, F₀.₉₉(p, n − p), where n is the number of samples. With n = 6 and p = 2, F₀.₉₉(2, 4) = 18.0. 18.80 > 18, so DESeq2 sets this gene's p to NA and keeps the LFC. So "an LFC but only p NA" is the mark of a Cook's outlier.

NAs in the results table appear in different columns depending on the cause.

| What the results table shows | Cause | What to check | Where checked |
|---|---|---|---|
| `baseMean = 0`, LFC, lfcSE, stat, p and padj all NA | Count 0 in every sample (`allZero`) | The raw counts | Experiment 4 |
| An LFC, but p and padj NA | A count outlier whose Cook's distance exceeds the cutoff | `assays(dds)[["cooks"]]`, `mcols(dds)$maxCooks` | Experiment 3 (a), (e) |
| A p, but padj NA | Independent filtering: baseMean below `filterThreshold`, so left out of the adjustment ([06](06_multiple_testing.md)) | `metadata(res)$filterThreshold` | Experiment 4, Problem 17 |
| baseMean above 0, yet everything from LFC to padj NA | The gene's coefficients cannot be estimated because of weights (`weightsFail`, handled like all-zero) | `mcols(dds)$weightsFail` | "When even p is NA" under Going deeper |
| (Not NA) a convergence warning | Both IRLS and its fallback, the optim calculation, failed to converge (`betaConv == FALSE`) | `mcols(dds)$betaConv`, `betaIter`, full rank of the design matrix | Experiment 4 (0 non-converged). Even with deliberately created non-convergence, p is not NA ("When even p is NA" under Going deeper) |
| (Not NA) LFC 0, stat 0, p 1 | Both groups compared with `contrast=` are all 0 (`cleanContrast`) | The counts of the two groups | Experiment 2 (gene357) |

Counting in the data of experiment 1 (3000 genes), there were 16 all-zero genes with everything NA and 579 genes with only padj NA. All 579 had a baseMean below the `filterThreshold` of 4.4313. There were no Cook's outliers (maximum 17.54 < 18) and no convergence failures (experiment 4).

But the Cook's filter does not always operate. DESeq2 makes the Cook's call only for samples in a design cell with 3 or more samples. A cell here is a set of samples sharing the same row of the design matrix, that is, samples with the same combination of condition and batch. I put 5000 into one gene of simulated data and checked by design (experiment 3).

| Design | What happens |
|---|---|
| 3 vs 3 | Cook's 36.9 > 18, so p is NA. With the filter off (`cooksCutoff=FALSE`), p = 0.00177 |
| 2 vs 2 | No cell with 3 or more, so `maxCooks` is all NA. With no filter, p = 0.00451 comes out as it is |
| 3 vs 2 | Only the samples of the group of 3 are judged. A large Cook's in the group of 2 is ignored |
| `~ pair + condition` (paired) | Every sample has a different row of the design matrix, so 1 per cell. No Cook's filter at all |
| 7 or more in a cell (e.g. 7 vs 7) | The jumping count is replaced with the trimmed mean (the mean without the top and bottom 20%) × size factor and the gene refitted (5000 → 5). The original count stays in `counts(dds)`. If every sample is replaceable, the Cook's filter of later `results()` calls is turned off |

If the shared example obtained the three conditions from the same donors and analysed them with `~ pair + condition`, the Cook's filter would not operate. In that case, whether one sample drives the conclusion has to be checked directly with `plotCounts()`.

Finally, an NA is not a biological conclusion. The Cook's outlier gene with p NA actually showed a large difference, LFC −7.9 (experiment 3). A gene with only padj NA was just left out of the adjustment because its counts are low. So don't read NA as "no difference in expression" or "not expressed". A rule that counts NA as "not a discovery" can be used in meta-analyses or sensitivity comparisons, but then record the rule and the cause of the NA (the column pattern in the table above) together.

## 6. What do VST and PCA show, and what don't they do?

Gene A's Ctrl mean is 106.7 and its Starvation mean 210. Under the negative binomial (NB) the variance is μ + αμ², so with α = 0.053147 the variances are 106.667 + 0.053147 × 106.667² ≈ 711 and 210 + 0.053147 × 210² ≈ 2554. The mean is about 2 times, but the variance 3.6 times. Like this, the variance of counts depends on the mean. VST (variance stabilizing transformation) and rlog (regularized log) are transformations that reduce this mean–variance dependence and the noise of low counts. The transformed values are used for exploratory plots such as PCA, clustering and heatmaps.

But the transformed values do not go into the test. DESeq2's NB GLM puts the mean–variance relationship of the raw counts and the size factors directly into the model, while the transformed values are on a log scale with that relationship already squashed, so the information the test needs is gone. Feeding transformed values into `DESeqDataSetFromMatrix()` is rejected because they are not integers, but that is only a symptom. Turning them back into integers, as with `round(2^vst)`, does not fix the problem.

The `blind` argument is often misunderstood too. `blind = TRUE` switches the design to `~1`, which uses no condition information at all, when estimating the dispersion trend, and `blind = FALSE` uses the current design as it is. But neither removes batch or condition effects from the transformed matrix.

Experiment 5 shows this. In 8 samples with condition (A, B) and batch (b1, b2) laid out in balance, I made 600 genes 2 times higher in b2, ran `DESeq()` with `~ batch + condition`, and drew a PCA from the `vst(blind = FALSE)` values. PCA (principal component analysis) summarizes the data with a few axes that show the largest differences between samples.

![PCA of VST data before and after removing batch for plotting](../../figures/07_pca_vst_batch.png)

On the left, even with `vst(blind = FALSE)`, PC2 splits by batch (batch means on PC2 ±5.58). The right panel is `limma::removeBatchEffect()` applied to a copy of the same matrix. This function estimates the batch effect in the transformed values with a linear model and subtracts it, which brought the batch means on PC2 to 0.

<details>
<summary>Code for the figure</summary>

```r
suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })
set.seed(4)                                                  # same data as experiment 5
dds <- makeExampleDESeqDataSet(n = 2000, m = 8, betaSD = 0.5)
dds$batch <- factor(rep(c("b1", "b2"), 4))
bf <- rep(1, nrow(dds)); bf[1:600] <- 2
cts <- counts(dds); j <- dds$batch == "b2"
cts[, j] <- matrix(rpois(sum(j) * nrow(dds), lambda = cts[, j] * bf), nrow = nrow(dds))
counts(dds) <- cts; design(dds) <- ~ batch + condition
dds <- DESeq(dds, quiet = TRUE)
vsd <- vst(dds, blind = FALSE)
vsdR <- vsd                                                  # copy for plotting only
assay(vsdR) <- limma::removeBatchEffect(assay(vsd), batch = vsd$batch,
                                        design = model.matrix(~ condition, colData(vsd)))
getpc <- function(v, title) {
  pd <- plotPCA(v, intgroup = c("condition", "batch"), returnData = TRUE)
  pv <- round(100 * attr(pd, "percentVar"))
  data.frame(pd[, c("PC1", "PC2", "condition", "batch")], panel = sprintf("%s\nPC1 %d%%, PC2 %d%%", title, pv[1], pv[2]))
}
df <- rbind(getpc(vsd, "vst(blind = FALSE)"), getpc(vsdR, "same + limma::removeBatchEffect (plot only)"))
df$panel <- factor(df$panel, unique(df$panel))
p <- ggplot(df, aes(PC1, PC2, colour = condition, shape = batch)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 3) + facet_wrap(~ panel) +
  scale_colour_manual(values = c(A = "#2a78d6", B = "#eb6834")) + scale_shape_manual(values = c(b1 = 16, b2 = 17)) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey95"))
out <- "../figures/07_pca_vst_batch.png"            # run from notes/
ggsave(out, p, width = 7, height = 4.5, dpi = 150, bg = "white")
print(aggregate(cbind(PC1, PC2) ~ panel + batch, df, function(x) round(mean(x), 2)))
cat("file exists:", file.exists(out), file.size(out), "bytes\n")
```

```
using ntop=500 top features by variance
using ntop=500 top features by variance
                                                          panel batch   PC1
1                          vst(blind = FALSE)\nPC1 24%, PC2 23%    b1  1.65
2 same + limma::removeBatchEffect (plot only)\nPC1 31%, PC2 16%    b1  0.00
3                          vst(blind = FALSE)\nPC1 24%, PC2 23%    b2 -1.65
4 same + limma::removeBatchEffect (plot only)\nPC1 31%, PC2 16%    b2  0.00
    PC2
1  5.58
2  0.00
3 -5.58
4  0.00
file exists: TRUE 29403 bytes
```

</details>

The batch-removed matrix on the right is for plotting only. The test is still done in `DESeq()` with the raw counts and the `~ batch + condition` design. If the design is changed so that batch and condition overlap completely, `DESeq()` refuses to fit with "full model matrix is less than full rank". In such a design, even if the PCA splits cleanly, the data cannot tell whether it is a condition effect or a batch effect.

There is also a small difference in the choice of function. The study text's example, `varianceStabilizingTransformation(dds, blind = FALSE)`, uses the trend `DESeq()` already fitted as it is. `vst(dds, blind = FALSE)`, on the other hand, picks 1000 genes (`nsub`) with baseMean > 5 and refits the trend, so its values differ slightly (up to 0.289). `rlog(blind = FALSE)` did not remove batch either (PC2 ±6.0).

The QC plots can be read like this.

- In a PCA (`plotPCA`, by default the top 500 genes by variance), look at batch, donor, RNA quality (quality variables in colData such as RIN) and unusual samples, not just condition.
- In a dispersion plot (`plotDispEsts`), look at whether the trend fits well and whether some genes lie far above it.
- In an MA plot (`plotMA`), tell apart the large MLEs of low-expression genes and how shrinkage squashes them (the MA plot in section 2).
- For genes of interest, take out the normalized counts per sample with `plotCounts(dds, gene, intgroup, returnData = TRUE)` and see whether one sample drives the conclusion. The returned `count` column is normalized count + 0.5 with the default (`transform = TRUE`).

A PCA not splitting as expected does not mean DESeq2 went wrong. Conversely, a clean split does not mean confounding in the design has been resolved. Data quality, experimental design and model diagnostics have to be looked at together.

## Summary

- The MLE LFC of genes with few counts jumps around a lot. gene1349 was −5.81 (56-fold), but its 95% interval ran from about a 5.6-fold to a 560-fold decrease.
- `lfcShrink()` uses a prior centred on 0 to pull LFCs with little information a lot and LFCs with lots of information a little. The factor is roughly σ²/(σ² + SE²), so gene A went 0.977 → 0.929 and gene1349 −5.81 → −1.97.
- In the default call `pvalue` and `padj` stay the same. Don't rebuild p by dividing the shrunken LFC by `lfcSE`; `lfcSE` also means different things depending on the function that made it (Wald SE, posterior SD, SE of a penalized estimate), so note what the bars in a figure are.
- The difference of two shrunken coefficients, the contrast path and the relevel path use different priors, so their values differ too. Record which path was used.
- Read the cause of an NA from the column pattern. Everything NA is all-zero, only p NA is a Cook's outlier, only padj NA is independent filtering. The Cook's filter operates only when a cell has 3 or more samples.
- Shrinkage (stabilizing estimates), testing and visualization (VST/rlog, PCA) do different jobs. Testing is done with the Wald test or the LRT (likelihood ratio test) ([05](05_wald_vs_lrt.md)), and `blind = FALSE` does not remove batch either.

## Exercises

**Problem 17.** You found a gene that has a raw p but whose padj alone is NA. Which procedure should be checked first?

<details>
<summary>Solution</summary>

Check independent filtering first. `results()` sets padj to NA for genes whose baseMean is below `filterThreshold` (this happens inside `pvalueAdjustment()`).

In experiment 4, 579 genes had a p but only padj NA. `metadata(res)$filterThreshold` was 4.4313, and the maximum baseMean of the 579 was 4.4307, all below the threshold. 0 genes at or above the threshold had only padj NA. Extracting again with `results(dds, independentFiltering = FALSE)` brings such NAs to 0.

If even p is NA, the cause is different. If baseMean is 0, it is all-zero. If baseMean is above 0, check whether it is a Cook's outlier with `mcols(dds)$maxCooks` above the cutoff F₀.₉₉(p, n − p) (experiment 3 (a): baseMean 783, LFC −7.9, p NA).

Of the study text's answer in appendix B, the independent filtering, all-zero and count outlier parts match reality. The "fitting problems" the study text also lists split into two cases. When IRLS did not converge (`betaConv = FALSE`), p does not become NA. On the other hand, rows whose coefficients cannot be estimated because of weights are flagged `weightsFail`, handled like all-zero, and get NA even for p. The data of this note have 0 non-converged genes, so I deliberately created both cases to check ("When even p is NA" under Going deeper). A solution of the same problem with other data (116 genes, filterThreshold 1.59) is in [06_multiple_testing.md](06_multiple_testing.md).

</details>

**Problem 18.** Can Z be recomputed from the shrunken LFC and lfcSE returned by apeglm to replace the default Wald p-value? Is the point estimate of ashr a MAP?

<details>
<summary>Solution</summary>

No to both.

First, Z must not be rebuilt. Doing the same calculation with the runnable `type = "normal"`, gene1349's Wald stat is −3.419 while `LFC_shr / SE_shr` is −3.142. 2Φ(−|Z|) from the latter is 0.00168, but the p stored by `results()` is 0.000628. apeglm's `lfcSE` puts in `fit$sd` (the posterior SD) in the source and changes the description to "posterior SD", so it is not a Wald SE. normal's `lfcSE` is the SE of a penalized estimate (a sandwich SE), not even a posterior SD (experiment 1), but the conclusion is the same. Dividing the two numbers does not give a defined test statistic.

And the point estimate of ashr is not a MAP but the posterior mean (`res$log2FoldChange <- fit$result$PosteriorMean`, label "MMSE"). It is apeglm and normal that return the posterior mode (MAP).

The conclusion is the same as appendix B of the study text. But the reason given, "the lfcSE of shrinkage is a posterior SD", holds only for apeglm and ashr, not for normal.

</details>

## Going deeper

<details>
<summary>Correspondence with the study text, and the reproduction environment</summary>

The sections of the body correspond to chapter 13 of the study text as follows: sections 1–2 ↔ 13.1–13.2, section 3 ↔ 13.3, section 4 ↔ the contrast part of 13.2, section 5 ↔ 13.5, and section 6 ↔ 13.4 and 13.6. Dispersion shrinkage is covered in [03_dispersion_estimation.md](03_dispersion_estimation.md), and multiple testing and independent filtering in [06_multiple_testing.md](06_multiple_testing.md).

The code was run with R 4.5.2 and DESeq2 1.50.2 on 2026-09-26. The `r` blocks of the experiments below were all run as they are, and the block right below each is its output. Experiments 4 and 6 were run continuing the R session of experiment 1 (reusing `dds`, `res`, `shr`, `d` and `pick`). The simulated-data code of sections 1–3 (gene1349) and the minmu block were run continuously in one session, and the gene A code and the figure code were each run separately (the figure code from the `notes/` folder).

The installation status and function signatures were checked like this.

```r
suppressPackageStartupMessages(library(DESeq2))
cat("apeglm:", requireNamespace("apeglm", quietly = TRUE), " ashr:", requireNamespace("ashr", quietly = TRUE), "\n")
args(lfcShrink)
f <- c(formals(DESeq)["minReplicatesForReplace"], formals(replaceOutliers)[c("trim", "minReplicates")],
       formals(DESeq2:::estimateBetaPriorVar)[c("betaPriorMethod", "upperQuantile")])
sapply(f, deparse)
```

```
apeglm: FALSE  ashr: FALSE
function (dds, coef, contrast, res, type = c("apeglm", "ashr",
    "normal"), lfcThreshold = 0, svalue = FALSE, returnList = FALSE,
    format = c("DataFrame", "GRanges", "GRangesList"), saveCols = NULL,
    apeAdapt = TRUE, apeMethod = "nbinomCR", parallel = FALSE,
    BPPARAM = bpparam(), quiet = FALSE, ...)
NULL
        minReplicatesForReplace                            trim
                            "7"                           "0.2"
                  minReplicates                 betaPriorMethod
                            "7" "c(\"weighted\", \"quantile\")"
                  upperQuantile
                         "0.05"
```

With no apeglm or ashr, every run in this note is `type="normal"`. What is said about apeglm/ashr was checked only in the `lfcShrink` source and help pages.

</details>

<details>
<summary>Calling the study text's example code as it is: apeglm and ashr stop at the installation check</summary>

I ran the study text's `lfcShrink(dds, coef="condition_Starvation_vs_Ctrl", type="apeglm")` from section 13.2 on synthetic data made to produce the same coefficient name (Ctrl reference, Starvation, Glucose, 3 each). I also tried a call with `type` omitted, a `contrast` call and `type="ashr"`.

```r
suppressPackageStartupMessages(library(DESeq2)); options(width = 120)
cat("apeglm:", requireNamespace("apeglm", quietly = TRUE), " ashr:", requireNamespace("ashr", quietly = TRUE), "\n")
set.seed(1)
dds <- makeExampleDESeqDataSet(n = 500, m = 9)
dds$condition <- factor(rep(c("Ctrl", "Starvation", "Glucose"), each = 3), levels = c("Ctrl", "Starvation", "Glucose"))
dds <- DESeq(dds, quiet = TRUE)
resultsNames(dds)
try1 <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) e)
  cat(sprintf("%-36s -> %s\n", label, if (inherits(r, "error")) paste("Error:", conditionMessage(r)) else paste("ran,", nrow(r), "rows")))
}
try1('coef, type="apeglm" [textbook 13.2]', lfcShrink(dds, coef = "condition_Starvation_vs_Ctrl", type = "apeglm"))
try1('coef, type omitted (default)',        lfcShrink(dds, coef = "condition_Starvation_vs_Ctrl"))
try1('contrast, type="apeglm"',             lfcShrink(dds, contrast = c("condition", "Glucose", "Starvation"), type = "apeglm"))
try1('coef, type="ashr"',                   lfcShrink(dds, coef = "condition_Starvation_vs_Ctrl", type = "ashr"))
try1('contrast, type="ashr"',               lfcShrink(dds, contrast = c("condition", "Glucose", "Starvation"), type = "ashr"))
try1('coef, type="normal"',                 lfcShrink(dds, coef = "condition_Starvation_vs_Ctrl", type = "normal", quiet = TRUE))
```

```
apeglm: FALSE  ashr: FALSE
[1] "Intercept"                    "condition_Starvation_vs_Ctrl" "condition_Glucose_vs_Ctrl"
coef, type="apeglm" [textbook 13.2]  -> Error: type='apeglm' requires installing the Bioconductor package 'apeglm'
coef, type omitted (default)         -> Error: type='apeglm' requires installing the Bioconductor package 'apeglm'
contrast, type="apeglm"              -> Error: type='apeglm' requires installing the Bioconductor package 'apeglm'
coef, type="ashr"                    -> Error: type='ashr' requires installing the CRAN package 'ashr'
contrast, type="ashr"                -> Error: type='ashr' requires installing the CRAN package 'ashr'
coef, type="normal"                  -> ran, 500 rows
```

Three things can be read from the results. First, the study text's coefficient name `condition_Starvation_vs_Ctrl` really is in `resultsNames()`. The study text's example call stops not because of the name but because apeglm is not installed. Next, the default of `type` is `"apeglm"` (`match.arg`, line 12), so omitting `type` gives the same error. To run shrinkage in this environment, `type="normal"` must be written out. Finally, `contrast` + apeglm is rejected in the source as "only for use with 'coef'" (lines 170–171), but the installation check (lines 167–168) comes first, so in this environment only the installation message appears. So "apeglm takes only coef" was confirmed from the source, not observed by running it. ashr stops at the installation check on lines 283–284 for both `coef` and `contrast`.

The study text's section 15.1, "depending on whether they are installed, optional parts are skipped with a message", describes `DESeq2_workbook.R` (that file is not in the repository, so I could not run it). `lfcShrink()` itself does not skip; it `stop()`s. To skip, the caller has to branch first with `requireNamespace()`.

</details>

<details>
<summary>Where does the MLE stop when one group is all 0? (minmu)</summary>

Run continuing from the code of section 1. `minmu` is, as the help page says, the "lower bound on the estimated count while fitting the GLM", and the default of `nbinomWaldTest()` is 0.5.

```r
for (mm in c(0.5, 0.1, 0.01)) {                      # minmu: lower bound on the expected count during the GLM fit (default 0.5)
  r <- results(nbinomWaldTest(dds, minmu = mm, quiet = TRUE))
  cat("minmu", mm, "| gene1349 LFC", round(r["gene1349", "log2FoldChange"], 3), "SE", round(r["gene1349", "lfcSE"], 3),
      "| gene2055 LFC", round(r["gene2055", "log2FoldChange"], 4), "\n")
}
```

```
minmu 0.5 | gene1349 LFC -5.808 SE 1.699 | gene2055 LFC 2.0626 
minmu 0.1 | gene1349 LFC -8.13 SE 2.904 | gene2055 LFC 2.0626 
minmu 0.01 | gene1349 LFC -11.452 SE 8.418 | gene2055 LFC 2.0626
```

The LFC and SE of gene1349 (B is 0, 0, 0) keep growing as the bound is lowered. There is no finite MLE, so the calculation stops leaning on the bound. gene2055, with many counts in both groups, is 2.0626 regardless of the bound.

</details>

<details>
<summary>Formula details: the normal prior, lfcSE, apeglm/ashr, the s-value, the Cook's cutoff, VST</summary>

The symbols of this note follow appendix C. But Cook's distance is written $\mathrm{Cook}_{ij}$ so as not to clash with the LRT statistic $D$ of appendix C. The number of samples $n$ and the number of coefficients $p$ are the symbols of section 4.5 of the study text (a $p$ different from a p-value; the DESeq2 source writes the number of samples as `m`).

First, the LFC MAP under the normal prior. Putting $\beta_k \sim N(0, \sigma_k^2)$ on the log2 coefficients $\beta$ of gene $i$, the posterior log objective is

$$
\ell(\beta) - \sum_k \frac{\beta_k^2}{2\sigma_k^2}
$$

and the $\hat\beta_{MAP}$ that maximizes it is returned. The implementation is the ridge-penalized IRLS of `nbinomWaldTest(betaPrior=TRUE, betaPriorVar=...)`. The prior variances $\sigma_k^2$ are on the log2 scale, but the IRLS runs on natural-log coefficients $b = \beta \ln 2$, so the penalty on $b$ is $\Lambda = \mathrm{diag}\big(1/(\sigma_k^2 \ln^2 2)\big)$ (the recomputation in experiment 1 used this scale and matched the returned values to the fourth decimal place). $\sigma_k^2$ is estimated by `estimateBetaPriorVar()`. With the default `betaPriorMethod="weighted"`, it drops allZero genes and $|\hat\beta| \ge 10$ and uses the normal variance matched to the weighted upper quantile (`upperQuantile=0.05`) of the MLE $\hat\beta$, with weights $1/\{1/\bar q_i + \alpha_{tr}(\bar q_i)\}$ (`1/(1/baseMean + dispFit)`). The intercept gets $10^6$, so effectively no prior. In the expanded model matrix (the `contrast=` path), `averagePriorsOverLevels()` averages over levels, so every level has the same variance (the 1.195 of experiment 2). Genes with lots of information have a large curvature of $\ell$, so the penalty is relatively small, while for low-expression genes $\ell$ is flat and they are pulled a long way toward 0.

The returned `lfcSE` is of sandwich form, built from $W$ at the MAP (multiplied by $1/\ln 2$ to convert to the $\log_2$ scale).

$$
\mathrm{lfcSE}_{normal} = \frac{1}{\ln 2}\sqrt{\left[(X^TWX+\Lambda)^{-1}\,X^TWX\,(X^TWX+\Lambda)^{-1}\right]_{kk}}
\quad\ne\quad
\frac{1}{\ln 2}\sqrt{\left[(X^TWX+\Lambda)^{-1}\right]_{kk}} \;(\text{Laplace posterior SD})
$$

`?lfcShrink` describes apeglm as an "adaptive Student's t prior shrinkage estimator", and the title of its reference is "Heavy-tailed prior distributions ..." (Zhu et al. 2018). It finds the posterior mode under this prior (`fit$map`). Thanks to the heavy tail, large effects are squashed less. Whether the prior is exactly Cauchy (a t with 1 degree of freedom) I could not check, because apeglm is not installed.

ashr puts a mixture prior $\beta \sim \sum_k \pi_k N(0, \tau_k^2)$ on $\hat\beta \mid \beta \sim N(\beta, \widehat{SE}^2)$, estimates $\pi_k$ from the data and returns the posterior mean. `lfcShrink` calls it as `ashr::ash(..., method = "shrink")`, and in the ashr source (GitHub `R/ash.R`, not installed locally) `method == "shrink"` means `pointmass = FALSE; prior = "uniform"`, so there is no point mass at 0.

Let me also define the lfsr and the s-value. Given a posterior, the local false sign rate is $\mathrm{lfsr}_i = \min\{P(\beta_i \ge 0 \mid K), P(\beta_i \le 0 \mid K)\}$ (the definition of Stephens 2017; if there is probability mass at 0, the equality signs matter). The s-value is defined as the cumulative mean after sorting the lfsr in ascending order, $s_{(r)} = \frac{1}{r}\sum_{u \le r} \mathrm{lfsr}_{(u)}$.

Cook's distance, as implemented in DESeq2, is computed with the number of coefficients $p$, the diagonal of the hat matrix $h_{jj}$ and a robust method-of-moments dispersion $\tilde\alpha_i$.

$$
V_{ij} = \mu_{ij} + \tilde\alpha_i \mu_{ij}^2, \qquad
\mathrm{Cook}_{ij} = \frac{(K_{ij} - \mu_{ij})^2}{V_{ij}} \cdot \frac{1}{p} \cdot \frac{h_{jj}}{(1-h_{jj})^2}
$$

The default cutoff is $F_{0.99}(p,\; n-p)$ (`qf(0.99, p, m - p)` in the source). For $n=6, p=2$ it is 18, for $n=5, p=2$ 30.82, and for $n=14, p=2$ 6.93.

Finally, VST (parametric trend). For the trend $\alpha_{tr}(x) = a_1/x + a_0$, `getVarianceStabilizedData` applies

$$
\mathrm{vst}(x) = \log_2\!\left( \frac{1 + a_1 + 2 a_0 x + 2\sqrt{a_0 x (1 + a_1 + a_0 x)}}{4 a_0} \right)
$$

to the observed normalized counts $x = K_{ij}/s_j$ (`extraPois` = $a_1$, `asymptDisp` = $a_0$). The $q_{ij}$ of appendix C is the expected value $\mu_{ij}/s_j$, so it is not used here.

</details>

<details>
<summary>A detailed comparison of the three types</summary>

| `type=` | prior | Point estimate returned | `lfcSE` returned | `coef` / `contrast` | Package needed |
|---|---|---|---|---|---|
| `apeglm` | adaptive Student's t (heavy-tailed; help-page wording) | posterior **mode** (`fit$map`) | posterior SD (`fit$sd`) | `coef` only | apeglm (not installed) |
| `ashr` | adaptive mixture of normals | posterior **mean** (`PosteriorMean`, labelled "MMSE") | posterior SD (`PosteriorSD`) | With `res` both are ignored; otherwise `results()` is built from `coef`/`contrast` | ashr (not installed) |
| `normal` | $N(0,\sigma_\beta^2)$ | posterior **mode** (ridge-penalized IRLS) | Sandwich SE of the penalized MAP estimate (not a posterior SD, experiment 1) | `coef` (standard model matrix) or `contrast` (expanded model matrix) | Built into DESeq2 |

Remember that of the three, only ashr returns the posterior mean. The study text's example in section 13.2 (Ctrl reference, `condition_Starvation_vs_Ctrl`, Glucose − Starvation) was run in experiment 2 as A=Ctrl, B=Starvation, C=Glucose.

</details>

<details>
<summary>What I checked in the lfcShrink source</summary>

The excerpts below are taken from the output of `deparse(DESeq2::function)`, with line numbers given (omissions marked `...`). They are not code to run.

`lfcShrink()` is used after running `DESeq()` with `betaPrior=FALSE` (the default). It stops if `resultsNames` is empty, with one exception: with `type="apeglm"` + `apeAdapt=FALSE` it does not stop and builds the model matrix column names from the design. An object run with `betaPrior=TRUE` stops regardless of type.

```text
# excerpt (not run)
# deparse(DESeq2::lfcShrink) lines 15–33
    if (length(resultsNames(dds)) == 0) {
        if (type != "apeglm" | (type == "apeglm" & apeAdapt)) {
            stop("first run DESeq() before running lfcShrink()")
        }
    }
    ...
    betaPrior <- attr(dds, "betaPrior")
    if (!is.null(betaPrior) && betaPrior) {
        stop("lfcShrink() should be used downstream of DESeq() with betaPrior=FALSE (the default)")
    }
    ...
    if (type == "apeglm" & !apeAdapt & length(resultsNames(dds)) ==
        0) {
        resultsNamesDDS <- colnames(model.matrix(design(dds),
            data = colData(dds)))
    }
```

The return values and labels by type are as follows. normal takes two columns (five if `lfcThreshold>0`) from the result of rerunning `nbinomWaldTest(betaPrior=TRUE)`, apeglm uses `map`/`sd`, and ashr `PosteriorMean`/`PosteriorSD`.

```text
# excerpt (not run)
# normal: lines 129–134, 153–162
        betaPriorVar <- estimateBetaPriorVar(dds, modelMatrix = modelMatrix)
        stopifnot(length(betaPriorVar) > 0)
        if (!parallel) {
            dds.shr <- nbinomWaldTest(dds, betaPrior = TRUE,
                betaPriorVar = betaPriorVar, modelMatrix = modelMatrix,
                modelMatrixType = modelMatrixType, quiet = TRUE)
    ...
        if (lfcThreshold > 0) {
            change.cols <- c("log2FoldChange", "lfcSE", "stat",
                "pvalue", "padj")
        }
        else {
            change.cols <- c("log2FoldChange", "lfcSE")
        }
        for (column in change.cols) {
            res[[column]] <- res.shr[[column]]
        }
# apeglm: lines 259–263
        res$log2FoldChange <- log2(exp(1)) * fit$map[, coefNum]
        res$lfcSE <- log2(exp(1)) * fit$sd[, coefNum]
        mcols(res)$description[2] <- sub("MLE", "MAP", mcols(res)$description[2])
        mcols(res)$description[3] <- sub("standard error", "posterior SD",
            mcols(res)$description[3])
# ashr: lines 290–296
        fit <- ashr::ash(betahat, sebetahat, mixcompdist = "normal",
            method = "shrink", ...)
        res$log2FoldChange <- fit$result$PosteriorMean
        res$lfcSE <- fit$result$PosteriorSD
        mcols(res)$description[2] <- sub("MLE", "MMSE", mcols(res)$description[2])
        mcols(res)$description[3] <- sub("standard error", "posterior SD",
            mcols(res)$description[3])
```

The restrictions on `coef`/`contrast`:

```text
# excerpt (not run)
# normal: lines 51–53, 108–113
            if (type == "normal" & is.numeric(contrast)) {
                stop("for type='normal' and numeric contrast, user must provide 'res' object")
            }
    ...
        if (missing(contrast)) {
            modelMatrixType <- "standard"
        }
        else {
            modelMatrixType <- "expanded"
            expMM <- makeExpandedModelMatrix(dds)
# apeglm: lines 170–172
        if (!missing(contrast)) {
            stop("type='apeglm' shrinkage only for use with 'coef'")
        }
```

The s-value path and the column layout differ between apeglm and ashr.

```text
# excerpt (not run)
# apeglm: lines 217–222, 264–280
        if (lfcThreshold > 0) {
            message(paste0("computing FSOS 'false sign or small' s-values (T=",
                round(lfcThreshold, 3), ")"))
            svalue <- TRUE
            apeT <- log(2) * lfcThreshold
        }
    ...
        if (svalue) {
            coefAlphaSpaces <- gsub("_", " ", coefAlpha)
            res <- res[, 1:3]
            if (lfcThreshold > 0) {
                res$svalue <- as.numeric(apeglm::svalue(fit$thresh))
                description <- paste0("FSOS s-value (T=", round(lfcThreshold,
                  3), "): ", coefAlphaSpaces)
            }
            else {
                res$svalue <- as.numeric(fit$svalue)
                description <- paste0("s-value: ", coefAlphaSpaces)
            }
            mcols(res)[4, ] <- DataFrame(type = "results", description = description)
        }
        else {
            res <- res[, c(1:3, 5:6)]
        }
# ashr: lines 298–318 (the lfcThreshold == 0 path and the > 0 path)
        if (lfcThreshold == 0) {
            if (svalue) {
                res <- res[, 1:3]
                res$svalue <- fit$result$svalue
            }
            else {
                res <- res[, c(1:3, 5:6)]
            }
        }
        else {
            message(paste0("computing FSOS 'false sign or small' s-values (T=",
                round(lfcThreshold, 3), ")"))
            svalue <- TRUE
    ...
            res$svalue <- apeglm_svalue(lfsr)
        }
# common: lines 320–323
    if (svalue) {
        mcols(res)[4, ] <- DataFrame(type = "results", description = paste("s-value:",
            coefAlphaSpaces))
    }
```

As a table:

| type | `lfcThreshold=0`, `svalue=FALSE` (default) | `svalue=TRUE` | `lfcThreshold>0` |
|---|---|---|---|
| normal | 6 columns, only `log2FoldChange`/`lfcSE` replaced (experiment 1) | Error `object 'coefAlphaSpaces' not found` (experiment 1) | `stat`/`pvalue`/`padj` too replaced by the threshold test of the shrunken fit (experiment 2) |
| apeglm | 5 columns `baseMean, log2FoldChange, lfcSE, pvalue, padj` (`stat` removed) | `pvalue`/`padj` removed, `svalue` added | Forced to `svalue=TRUE`, so `pvalue`/`padj` removed and the FSOS s-value added |
| ashr | 5 columns (`stat` removed) | `pvalue`/`padj` removed, `svalue` added | No columns removed. Keeps the MLE `stat`/`pvalue`/`padj` and appends the FSOS `svalue` |

ashr + `lfcThreshold>0` keeps the p-values, unlike the help page's description of `svalue` ("should p-values and adjusted p-values be replaced with s-values"). Also, on this path the common block appears to overwrite the description of the 4th column (`stat`) with "s-value: ...". But without ashr I could confirm neither by running it.

</details>

<details>
<summary>The results() source: the Cook's cutoff, NA handling, betaConv</summary>

```text
# excerpt (not run)
# deparse(DESeq2::results) lines 213–218, 225–246
    m <- nrow(attr(object, "dispModelMatrix"))
    p <- ncol(attr(object, "dispModelMatrix"))
    defaultCutoff <- qf(0.99, p, m - p)
    if (missing(cooksCutoff)) {
        cooksCutoff <- defaultCutoff
    }
    ...
        cooksOutlier <- mcols(object)$maxCooks > cooksCutoff
        if (any(cooksOutlier, na.rm = TRUE) & is(design(object),
            "formula")) {
            designVars <- all.vars(design(object))
            if (length(designVars) == 1) {
                var <- colData(object)[[designVars]]
                if (is(var, "factor") && nlevels(var) == 2) {
                  dontFilter <- logical(sum(cooksOutlier, na.rm = TRUE))
                  for (i in seq_along(dontFilter)) {
                    ii <- which(cooksOutlier)[i]
                    outCount <- counts(object)[ii, which.max(assays(object)[["cooks"]][ii,
                      ])]
                    if (sum(counts(object)[ii, ] > outCount) >=
                      3) {
                      dontFilter[i] <- TRUE
                    }
                  }
                  cooksOutlier[which(cooksOutlier)][dontFilter] <- FALSE
                }
            }
        }
        res$pvalue[which(cooksOutlier)] <- NA
```

`res$pvalue[which(cooksOutlier)] <- NA`, so a Cook's outlier gets NA only for p, and LFC/lfcSE/stat remain. There is an exception when the design has only a single 2-level factor: if 3 or more counts, across both groups, are larger than the count of the sample with the largest Cook's (`outCount`), it is not filtered. The common case is a value that jumps within its own group while the other group is all larger (experiment 3-(a2): A's 400 against B's 1000, 1100, 950). So it does not mean "only low outliers are let off". `betaConv` does not appear in the `results` source. Non-convergence does not create NA; it is reported only by the message of `nbinomWaldTest` ("rows did not converge in beta, labelled in mcols(object)\$betaConv").

`betaConv` is not just the result of IRLS either. Rows where IRLS did not converge or was unstable are refitted with optim, and `betaConv` is overwritten with the optim result. So `betaConv == FALSE` means that both IRLS and the optim fallback failed to converge.

```text
# excerpt (not run)
# deparse(DESeq2:::fitNbinomGLMs) lines 136, 143–145, 152–159
    betaConv <- betaRes$iter < maxit
    ...
    rowsForOptim <- if (useOptim) {
        which(!betaConv | !rowStable | !rowVarPositive)
    }
    ...
    if (length(rowsForOptim) > 0) {
        resOptim <- fitNbinomGLMsOptim(object, modelMatrix, lambda,
            rowsForOptim, rowStable, normalizationFactors, alpha_hat,
            weights, useWeights, betaMatrix, betaSE, betaConv,
            beta_mat, mu, logLike, minmu = minmu)
        betaMatrix <- resOptim$betaMatrix
        betaSE <- resOptim$betaSE
        betaConv <- resOptim$betaConv
```

</details>

<details>
<summary>When even p is NA: non-convergence of β and weights (appendix B 17)</summary>

The answer to appendix B 17 lists all-zero, count outliers and "fitting problems" as the causes when even p is NA. The data of this note have no non-converged genes at all, so I split the "fitting problems" into two cases and created them deliberately. (a) is non-convergence of β: IRLS is cut off at 2 iterations, both with the optim refit on and off. (b) is the case where coefficients cannot be estimated because of weights: setting the weights of all group B samples of gene1 to 0 leaves no information for estimating the mean of B.

The two blocks below were run in a fresh R session, separate from the other experiments of this note.

```r
suppressPackageStartupMessages(library(DESeq2))
set.seed(1)
dds <- makeExampleDESeqDataSet(n = 1000, m = 6)
dds <- estimateDispersions(estimateSizeFactors(dds), quiet = TRUE)

# (a) non-convergence of β: cut IRLS at 2 iterations (optim fallback on / off)
ok   <- nbinomWaldTest(dds, quiet = TRUE)
opt  <- nbinomWaldTest(dds, maxit = 2, quiet = TRUE)
bad  <- nbinomWaldTest(dds, maxit = 2, useOptim = FALSE, quiet = TRUE)
ref  <- results(ok)
for (d in list(opt, bad)) {
  r  <- results(d); nc <- which(!mcols(d)$betaConv)
  cat("betaConv FALSE:", length(nc), "| of these pvalue NA:", sum(is.na(r$pvalue[nc])),
      "| baseMean>0:", sum(r$baseMean[nc] > 0),
      "| max|log2FC - default fit|:", signif(max(abs(r$log2FoldChange[nc] - ref$log2FoldChange[nc])), 3),
      "| pvalue NA total:", sum(is.na(r$pvalue)), "\n")
}
cat("default fit  betaConv FALSE:", sum(!mcols(ok)$betaConv, na.rm = TRUE),
    "| pvalue NA:", sum(is.na(ref$pvalue)), "| allZero:", sum(mcols(ok)$allZero), "\n")

# (b) a row whose design is degenerated by weights: all group B weights of gene1 set to 0
w <- matrix(1, nrow(dds), ncol(dds), dimnames = dimnames(dds)); w[1, dds$condition == "B"] <- 0
assays(dds)[["weights"]] <- w
ddsW <- DESeq(dds, quiet = TRUE)
print(cbind(counts(ddsW)[1, , drop = FALSE], weightsFail = mcols(ddsW)$weightsFail[1],
            allZero = mcols(ddsW)$allZero[1]))
print(results(ddsW)[1, ])
```

```
betaConv FALSE: 3 | of these pvalue NA: 0 | baseMean>0: 3 | max|log2FC - default fit|: 4.68e-08 | pvalue NA total: 8 
betaConv FALSE: 992 | of these pvalue NA: 0 | baseMean>0: 992 | max|log2FC - default fit|: 11.7 | pvalue NA total: 8 
default fit  betaConv FALSE: 0 | pvalue NA: 8 | allZero: 8 
Warning message:
In getAndCheckWeights(object, modelMatrix, weightThreshold = weightThreshold) :
  for 1 row(s), the weights as supplied won't allow parameter estimation, producing a
  degenerate design matrix. These rows have been flagged in mcols(dds)$weightsFail
  and treated as if the row contained all zeros (mcols(dds)$allZero set to TRUE).
  If you are blocking for donors/organisms, consider design = ~0+donor+condition.
      sample1 sample2 sample3 sample4 sample5 sample6 weightsFail allZero
gene1      14       1       1       4       3      14           1       1
log2 fold change (MLE): condition B vs A 
Wald test p-value: condition B vs A 
DataFrame with 1 row and 6 columns
       baseMean log2FoldChange     lfcSE      stat    pvalue      padj
      <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
gene1   2.51101             NA        NA        NA        NA        NA
```

Starting with (a), the default fit has 0 non-converged genes. With the optim refit on, 3 remain where both IRLS and optim failed, and the estimates of these 3 match the default fit within $5\times10^{-8}$. With optim off, 992 are non-converged and their log2FC is off by up to 11.7. But in both cases the p of the non-converged rows is not NA. The 8 with p NA are all all-zero rows. Non-convergence attaches a normal-looking p to an estimate that may be wrong, so it does not show up as NA. `mcols(dds)$betaConv` has to be checked directly.

In (b), `getAndCheckWeights()` flags rows whose design degenerated because of weights as `weightsFail` and sets `allZero` to TRUE (lines 33–34 in the source search below). So gene1 has counts that are not 0 and a baseMean of 2.51, yet everything from log2FoldChange to padj is NA. Of the study text's "fitting problems", only this case, where the coefficients cannot be estimated at all, makes p NA.

As a table:

| Cause | baseMean | log2FoldChange · lfcSE · stat | pvalue | Column to check | Study text appendix B 17 |
|---|---|---|---|---|---|
| all-zero | 0 | NA | NA | `allZero` | Matches |
| Cook's outlier | > 0 | Values present | NA | `maxCooks` | Matches (experiment 3-(a)) |
| Design degenerated by weights | > 0 | NA | NA | `weightsFail`, `allZero` | Matches (the cannot-estimate part of "fitting problems") |
| Non-convergence of β | > 0 | Values present (not trustworthy) | Value present | `betaConv` | Does not match (p is not NA) |

But the excerpt above alone cannot say there are no other places that put in NA. So I `deparse()`d `results()` and every function that produces p, and searched for all lines containing `NA`, `betaConv` or an assignment to p. The line numbers are those of the default `deparse()`.

```r
suppressPackageStartupMessages(library(DESeq2))
ns <- asNamespace("DESeq2")
show_lines <- function(f, pat) {
  d <- deparse(get(f, envir = ns))
  i <- grep(pat, d)
  i <- i[!grepl("assays\\(", d[i])]                                  # skip the NA rows of the mu/H/cooks assays
  i <- sort(c(i, i[grepl("DataFrameWithNARows\\(.*, *$", d[i])] + 1))  # include the NA-row argument
  cat(sprintf("%s (%d lines): %s\n", f, length(d), toString(i)))
  cat(sprintf("  %3d: %s\n", i, trimws(d[i])), sep = "")
}
for (f in c("results", "getPvalue", "cleanContrast", "getContrast", "nbinomWaldTest", "nbinomLRT"))
  show_lines(f, "NA|[bB]etaConv|[pP]value.*<-")   # NA, betaConv, assignment to p
show_lines("getAndCheckWeights", "allZero.*<-|weightsFail =")
```

```
results (290 lines): 119, 160, 167, 186, 192, 196, 198, 200, 204, 208, 211, 246, 257
  119: pvalue <- getPvalue(object, test, name)
  160: newPvalue <- mapply(pfunc_lfc, T, LFC, SE)
  167: newPvalue <- mapply(pfunc_lfc, T, LFC, SE, df)
  186: newPvalue <- mapply(pfunc_lfc, lfcThreshold,
  192: newPvalue <- pmin(1, 2 * pfunc((abs(LFC) - T)/SE))
  196: pvalueAbove <- pfunc((T - LFC)/SE)
  198: pvalueBelow <- pfunc((LFC + T)/SE)
  200: newPvalue <- pmax(pvalueAbove, pvalueBelow)
  204: newPvalue <- pfunc((LFC - T)/SE)
  208: newPvalue <- pfunc((-T - LFC)/SE)
  211: res$pvalue <- newPvalue
  246: res$pvalue[which(cooksOutlier)] <- NA
  257: res$pvalue[nowZero] <- 1
getPvalue (15 lines): 
cleanContrast (202 lines): 44, 75, 91, 191, 195
   44: pvalue <- getPvalue(object, test, name)
   75: pvalue <- getPvalue(object, test, swapName)
   91: "lfcSE", "stat", "pvalue"), "description"] <- contrastDescriptions
  191: res$pvalue[contrastAllZero] <- 1
  195: pvalue <- getPvalue(object, test, name = NULL)
getContrast (58 lines): 45, 49, 53, 54
   45: contrastPvalue <- 2 * pt(abs(contrastStatistic), df = df,
   49: contrastPvalue <- 2 * pnorm(abs(contrastStatistic), lower.tail = FALSE)
   53: contrastResults <- buildDataFrameWithNARows(contrastList,
   54: mcols(object)$allZero)
nbinomWaldTest (191 lines): 138, 140, 143, 145, 146, 147, 149, 162, 164
  138: df <- ifelse(df > 0, df, NA)
  140: WaldPvalue <- 2 * pt(abs(WaldStatistic), df = df, lower.tail = FALSE)
  143: WaldPvalue <- 2 * pnorm(abs(WaldStatistic), lower.tail = FALSE)
  145: colnames(WaldPvalue) <- paste0("WaldPvalue_", modelMatrixNames)
  146: betaConv <- fit$betaConv
  147: if (any(!betaConv)) {
  149: message(paste(sum(!betaConv), "rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest"))
  162: list(betaConv = betaConv, betaIter = fit$betaIter, deviance = -2 *
  164: WaldResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
nbinomLRT (184 lines): 79, 120, 129, 131, 133, 135, 149, 151, 156, 159
   79: LRTPvalue <- pchisq(LRTStatistic, df = df, lower.tail = FALSE)
  120: LRTPvalue <- qlr$pval
  129: betaSE = array(NA, dim(fit_full$Beta), dimnames = list(rownames(fit_full$Beta),
  131: betaConv = rep(TRUE, nrow(objectNZ)), betaIter = rep(NA,
  133: reducedModel <- list(betaConv = rep(TRUE, nrow(objectNZ)))
  135: maxCooks <- rep(NA, nrow(objectNZ))
  149: if (any(!fullModel$betaConv)) {
  151: message(paste(sum(!fullModel$betaConv), "rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT"))
  156: fullBetaConv = fullModel$betaConv, reducedBetaConv = reducedModel$betaConv,
  159: LRTResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
getAndCheckWeights (48 lines): 33, 34
   33: mcols(object)$allZero[!weights.ok] <- TRUE
   34: weightsDF <- DataFrame(weightsFail = !weights.ok)
```

- In `results()` there are four places that put a value into p. Line 119 takes out the stored p, and lines 160–211 recompute p from the LFC and SE depending on `lfcThreshold` and `altHypothesis`. Line 246 puts NA into Cook's outliers, and line 257 puts 1 into rows whose baseMean became 0 after outlier replacement. Line 246 is the only one that puts in NA directly, and `betaConv` appears nowhere in the 290 lines. `getPvalue()` has no such line, and `cleanContrast()` only puts 1 into rows where both groups compared are 0 (line 191).
- Lines 53–54 of `getContrast()`, line 164 of `nbinomWaldTest()` and line 159 of `nbinomLRT()` decide the rows to fill with NA by `allZero` alone. Lines 33–34 of `getAndCheckWeights()` set `allZero` to TRUE for rows whose design degenerated because of weights, so the row in (b) also becomes NA via this path.
- The other p formulas (lines 45 and 49 of `getContrast()`, lines 140 and 143 of `nbinomWaldTest()`, line 79 of `nbinomLRT()`) only compute p from the statistic, so p becomes NA only when the statistic or the degrees of freedom are NA. Line 138 of `nbinomWaldTest()` (`df <- ifelse(df > 0, df, NA)`) is such a case, but it is inside the non-default `useT=TRUE` branch, and I did not run it. In (a), the rows with p NA numbered 8 in all three fits, the same as the number of all-zero rows.
- Lines 120–135 of `nbinomLRT()` are the glmGamPoi branch, which sets all `betaConv` to TRUE. The remaining `betaConv` lines print a message or store the `mcols` column. Non-convergence in the LRT was checked only in the source.

</details>

<details>
<summary>Does results(contrast=) refit the coefficients?</summary>

It does not. `getContrast()` calls `fitBeta(..., maxitSEXP = 0)`, so it leaves the coefficients as they are and computes only $c^\top\hat b$ and $\sqrt{c^\top\hat\Sigma c}$ (the excerpt and the reproduction of the covariance are in the contrast part of [note 04](04_glm_condition_batch.md)). `cleanContrast()` overwrites the LFC and stat with 0 and the p with 1 for genes where both groups compared are 0 but which are not all-zero overall (the run excerpt is in the implementation check of [06](06_multiple_testing.md)). So the MLE contrast equals the difference of the coefficients exactly, the only exception being genes where both groups compared are 0 (experiment 2).

</details>

<details>
<summary>The replacement conditions of DESeq() and the definition of maxCooks (source)</summary>

```text
# excerpt (not run)
# deparse(DESeq2::DESeq) lines 139–143
    sufficientReps <- any(nOrMoreInCell(attr(object, "modelMatrix"),
        minReplicatesForReplace))
    if (sufficientReps) {
        object <- refitWithoutOutliers(object, test = test, betaPrior = betaPrior,
            full = full, reduced = reduced, quiet = quiet, minReplicatesForReplace = minReplicatesForReplace,
# deparse(DESeq2:::nOrMoreInCell): are there n or more samples sharing the same model matrix row (design cell)?
function (modelMatrix, n)
{
    row_hash <- apply(modelMatrix, 1, paste0, collapse = "_")
    hash_table <- table(row_hash)
    numEqual <- as.vector(unname(hash_table[row_hash]))
    numEqual >= n
}
# deparse(DESeq2::replaceOutliers) lines 8–10, 19–21, 27–28, 33–35
    if (minReplicates < 3) {
        stop("at least 3 replicates are necessary in order to indentify a sample as a count outlier")
    }
    ...
    if (missing(cooksCutoff)) {
        cooksCutoff <- qf(0.99, p, m - p)
    }
    ...
    trimBaseMean <- apply(counts(object, normalized = TRUE),
        1, mean, trim = trim)
    ...
    else {
        as.integer(outer(trimBaseMean, sizeFactors(object), "*"))
    }
# deparse(DESeq2:::recordMaxCooks)
function (design, colData, modelMatrix, cooks, numRow)
{
    samplesForCooks <- nOrMoreInCell(modelMatrix, n = 3)
    p <- ncol(modelMatrix)
    m <- nrow(modelMatrix)
    maxCooks <- if ((m > p) & any(samplesForCooks)) {
        apply(cooks[, samplesForCooks, drop = FALSE], 1, max)
    }
    else {
        rep(NA, numRow)
    }
    maxCooks
}
# deparse(DESeq2:::refitWithoutOutliers) lines 48–64
        if (all(object$replaceable)) {
            mcols(object)$maxCooks <- NA
        }
        else {
            replaceCooks <- assays(object)[["cooks"]]
            replaceCooks[, object$replaceable] <- 0
            mcols(object)$maxCooks <- recordMaxCooks(design(object),
                colData(object), attr(object, "dispModelMatrix"),
                replaceCooks, nrow(object))
        }
    }
    if (nrefit > 0) {
        assays(object, withDimnames = FALSE)[["replaceCounts"]] <- counts(object)
        assays(object, withDimnames = FALSE)[["replaceCooks"]] <- assays(object)[["cooks"]]
        counts(object) <- assays(object)[["originalCounts"]]
        assays(object, withDimnames = FALSE)[["cooks"]] <- cooks
        assays(object)[["originalCounts"]] <- NULL
```

Three rules come out of this. First, only samples in a cell with 3 or more go into the `maxCooks` calculation. If there are no such samples (2 vs 2, or a paired design with 1 per cell), `maxCooks` is all NA and there is no Cook's filter at all. If it is mixed, like 3 vs 2, only the samples of the group of 3 are judged, and a large Cook's in the group of 2 is ignored (experiment 3-(e)). Second, if any cell has 7 or more, the outlier counts of the samples in that cell are replaced with the trimmed mean × size factor and refitted, and if every sample is replaceable, the Cook's filter of later `results()` calls is turned off too. Third, a paired design like `~pair+condition` has a different model matrix row for every sample, so every cell has size 1 and it falls under the first rule. The study text's sentence, "don't describe it as if it always happens in every paired design", agrees with this rule.

</details>

<details>
<summary>The VST, vst, plotCounts and plotPCA source</summary>

```text
# excerpt (not run)
# deparse(DESeq2::varianceStabilizingTransformation) lines 17–25
    if (blind) {
        design(object) <- ~1
    }
    if (blind | is.null(attr(dispersionFunction(object), "fitType"))) {
        object <- estimateDispersionsGeneEst(object, quiet = TRUE)
        object <- estimateDispersionsFit(object, quiet = TRUE,
            fitType)
    }
    vsd <- getVarianceStabilizedData(object)
# deparse(DESeq2::vst) lines 15–17, 23–38: refit the trend on a subset, then call VST with blind = FALSE
        if (blind) {
            design(object) <- ~1
        }
    ...
    baseMean <- MatrixGenerics::rowMeans(counts(object, normalized = TRUE))
    if (sum(baseMean > 5) < nsub) {
        stop("less than 'nsub' rows with mean normalized count > 5, \n  it is recommended to use varianceStabilizingTransformation directly")
    }
    object.sub <- object[baseMean > 5, ]
    baseMean <- baseMean[baseMean > 5]
    o <- order(baseMean)
    idx <- o[round(seq(from = 1, to = length(o), length = nsub))]
    object.sub <- object.sub[idx, ]
    object.sub <- estimateDispersionsGeneEst(object.sub, quiet = TRUE)
    object.sub <- estimateDispersionsFit(object.sub, fitType = fitType,
        quiet = TRUE)
    suppressMessages({
        dispersionFunction(object) <- dispersionFunction(object.sub)
    })
    vsd <- varianceStabilizingTransformation(object, blind = FALSE)
# deparse(DESeq2::plotCounts) lines 15–18
    if (missing(pc)) {
        pc <- if (transform)
            0.5
        else 0
```

The defaults of the DESeqTransform method of `plotPCA` are `intgroup = "condition", ntop = 500, returnData = FALSE`.

</details>

<details>
<summary>Experiment 1. What happens when low-expression genes with a large MLE and SE are shrunk? (type="normal")</summary>

```r
suppressPackageStartupMessages(library(DESeq2)); options(width = 120)
set.seed(7)
dds <- makeExampleDESeqDataSet(n = 3000, m = 6, betaSD = 1)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds)
shr <- lfcShrink(dds, coef = "condition_B_vs_A", type = "normal", quiet = TRUE)
print(data.frame(results = mcols(res)$description, lfcShrink = mcols(shr)$description))
cat("pvalue identical:", identical(res$pvalue, shr$pvalue), " padj identical:", identical(res$padj, shr$padj),
    " colnames(shr):", colnames(shr), "\n")
print(priorInfo(shr)$betaPriorVar)
cat("svalue=TRUE + normal:", tryCatch(lfcShrink(dds, coef = 2, type = "normal", svalue = TRUE, quiet = TRUE),
                                      error = conditionMessage), "\n")
d <- data.frame(baseMean = res$baseMean, LFC_MLE = res$log2FoldChange, SE_MLE = res$lfcSE,
                LFC_shr = shr$log2FoldChange, SE_shr = shr$lfcSE, stat = res$stat,
                pvalue = res$pvalue, padj = res$padj, row.names = rownames(res))
d$Z_recomputed <- d$LFC_shr / d$SE_shr
sel <- which(d$baseMean < 5 & abs(d$LFC_MLE) > 3 & d$SE_MLE > 1)
hi  <- which(d$baseMean > 500 & abs(d$LFC_MLE) > 2)
pick <- c(sel[order(-abs(d$LFC_MLE[sel]))][1:4], hi[order(d$SE_MLE[hi])][1:2])
cat("n(baseMean<5 & |LFC_MLE|>3 & SE_MLE>1):", length(sel), "\n")
print(signif(d[pick, ], 4))
ok <- !is.na(d$LFC_MLE) & abs(d$LFC_MLE) > 0.5          # only genes with |LFC_MLE| > 0.5
print(round(tapply(abs(d$LFC_shr[ok]) / abs(d$LFC_MLE[ok]),
                   cut(d$baseMean[ok], c(0, 5, 20, 100, 1000, Inf)), median), 3))
cat("fold = 2^|LFC_shr| (first 3):", round(2^abs(d$LFC_shr[pick[1:3]]), 2),
    "| 2*pnorm(-|Z_recomputed|) (gene1349):", signif(2 * pnorm(-abs(d$Z_recomputed[pick[1]])), 3), "\n")
```

```
                                    results                                 lfcShrink
1 mean of normalized counts for all samples mean of normalized counts for all samples
2  log2 fold change (MLE): condition B vs A  log2 fold change (MAP): condition B vs A
3          standard error: condition B vs A          standard error: condition B vs A
4          Wald statistic: condition B vs A          Wald statistic: condition B vs A
5       Wald test p-value: condition B vs A       Wald test p-value: condition B vs A
6                      BH adjusted p-values                      BH adjusted p-values
pvalue identical: TRUE  padj identical: TRUE  colnames(shr): baseMean log2FoldChange lfcSE stat pvalue padj
   Intercept   conditionB
1.000000e+06 1.591035e+00
svalue=TRUE + normal: object 'coefAlphaSpaces' not found
n(baseMean<5 & |LFC_MLE|>3 & SE_MLE>1): 81
         baseMean LFC_MLE SE_MLE LFC_shr SE_shr   stat    pvalue      padj Z_recomputed
gene1349    4.886  -5.808 1.6990 -1.9690 0.6266 -3.419 6.278e-04 8.435e-03      -3.1420
gene1107    4.738  -5.764 1.7320 -1.8880 0.6244 -3.328 8.746e-04 1.096e-02      -3.0230
gene2003    4.434  -5.668 1.8080 -1.7170 0.6182 -3.134 1.722e-03 1.874e-02      -2.7780
gene2613    3.950   5.457 3.7540  0.3516 0.3931  1.454 1.461e-01        NA       0.8944
gene2055 1206.000   2.063 0.3131  1.9430 0.2949  6.587 4.487e-11 1.199e-08       6.5880
gene80    786.300   2.436 0.3431  2.2680 0.3193  7.099 1.254e-12 6.973e-10       7.1020
      (0,5]      (5,20]    (20,100] (100,1e+03] (1e+03,Inf]
      0.270       0.616       0.817       0.895       0.941
fold = 2^|LFC_shr| (first 3): 3.91 3.7 3.29 | 2*pnorm(-|Z_recomputed|) (gene1349): 0.00168
```

The labels show that only `log2FoldChange` changes from MLE → MAP, while the `lfcSE` label stays "standard error" under normal (apeglm/ashr rewrite it as "posterior SD" in the source). `pvalue` and `padj` are identical down to the values. Giving `svalue=TRUE` together with normal makes the common block of the source (lines 320–323) refer to the undefined `coefAlphaSpaces`, and it errors.

The genes picked for the table are the top 4 by $|\hat\beta_{MLE}|$ among the 81 with `baseMean < 5`, $|\hat\beta_{MLE}| > 3$ and $SE > 1$, and the 2 with the smallest SE among those with `baseMean > 500` and $|\hat\beta_{MLE}| > 2$. `stat` is the MLE Wald statistic, and `Z_recomputed = LFC_shr / SE_shr` is the calculation that should not be done, done on purpose.

The low-expression genes whose MLE was |LFC| 5.5–5.8 (44–56-fold, larger than the study text's "log2FC 5 = 32-fold") were squashed to LFC −1.7 to −2.0, about 3.3–3.9-fold. The SE also dropped from 1.7 to 0.6, but this 0.6 is a sandwich SE, so the posterior uncertainty is larger than that (below). gene2613, with an SE of 3.75 and almost no information, was squashed to 0.35, while the high-expression genes with SEs around 0.3 shrank by only 6–7%. The median of $|\hat\beta_{shr}|/|\hat\beta_{MLE}|$ per baseMean bin, for genes with $|\hat\beta_{MLE}| > 0.5$ only, ran from 0.270 to 0.941, closer to 1 the higher the expression. The last column leads into the answer to Problem 18. gene1349's Wald stat is −3.419, but dividing the shrunk values gives −3.142, and the p from that is 0.00168, different from the stored 0.000628.

Next I checked what kind of SE the `lfcSE` of normal is. I redid the MAP fit with the prior variances `lfcShrink` used and computed the two formulas with $W$ at that MAP.

```r
pv  <- priorInfo(shr)$betaPriorVar                 # prior variances used by lfcShrink (log2 scale)
fit <- nbinomWaldTest(dds, betaPrior = TRUE, betaPriorVar = pv, modelMatrixType = "standard", quiet = TRUE)
X <- model.matrix(design(dds), colData(dds))
L <- diag(1 / (pv * log(2)^2))                     # ridge penalty on the natural-log coefficients b
se <- t(sapply(rownames(d)[pick[c(1, 4, 5)]], function(g) {
  b  <- log(2) * unlist(mcols(fit)[g, c("Intercept", "condition_B_vs_A")])
  mu <- pmax(sizeFactors(dds) * exp(drop(X %*% b)), 0.5)
  I  <- crossprod(X, X * (mu / (1 + dispersions(dds)[rownames(dds) == g] * mu)))   # X'WX
  P  <- solve(I + L)                                                              # (X'WX + Λ)^-1
  c(lfcSE_shr = shr[g, "lfcSE"], sandwich = log2(exp(1)) * sqrt((P %*% I %*% P)[2, 2]),
    laplace_postSD = log2(exp(1)) * sqrt(P[2, 2]), SE_MLE = res[g, "lfcSE"])
}))
print(round(se, 4))
ci <- function(est, s) sprintf("[%.2f, %.2f]", est - 1.96 * s, est + 1.96 * s)
g3 <- rownames(se)
print(data.frame(Wald95_MLE = ci(res[g3, "log2FoldChange"], res[g3, "lfcSE"]),
                 MAP_1.96_lfcSE = ci(shr[g3, "log2FoldChange"], se[, "lfcSE_shr"]),
                 MAP_1.96_laplaceSD = ci(shr[g3, "log2FoldChange"], se[, "laplace_postSD"]), row.names = g3))
```

```
         lfcSE_shr sandwich laplace_postSD SE_MLE
gene1349    0.6266   0.6266         0.9413 1.6987
gene2613    0.3931   0.3931         1.1906 3.7540
gene2055    0.2949   0.2949         0.3038 0.3131
             Wald95_MLE MAP_1.96_lfcSE MAP_1.96_laplaceSD
gene1349 [-9.14, -2.48] [-3.20, -0.74]     [-3.81, -0.12]
gene2613 [-1.90, 12.81]  [-0.42, 1.12]      [-1.98, 2.69]
gene2055   [1.45, 2.68]   [1.36, 2.52]       [1.35, 2.54]
```

The returned `lfcSE` equals the sandwich formula to the fourth decimal place. It differs from the Laplace-approximation posterior SD (gene1349 0.63 vs 0.94, gene2613 0.39 vs 1.19), while for the high-expression gene2055 the prior has little influence and all three are similar. The second table of the output shows the same "±1.96 × lfcSE" becoming three different intervals. gene2613's Wald interval is [−1.90, 12.81], which includes 0, while normal's MAP ± 1.96·lfcSE is [−0.42, 1.12], narrow enough to make the uncertainty look small. This is why the caption must say which one the bars in a figure are.

</details>

<details>
<summary>Experiment 2. Does subtracting two shrunken coefficients give the shrinkage of the contrast? (A=Ctrl, B=Starvation, C=Glucose)</summary>

```r
suppressPackageStartupMessages(library(DESeq2)); options(width = 120)
set.seed(3)
dds <- makeExampleDESeqDataSet(n = 2000, m = 9, betaSD = 1)
dds$condition <- factor(rep(c("A", "B", "C"), each = 3))      # A=Ctrl, B=Starvation, C=Glucose
dds <- DESeq(dds, quiet = TRUE)
rBA <- results(dds, name = "condition_B_vs_A"); rCA <- results(dds, name = "condition_C_vs_A")
rCB <- results(dds, contrast = c("condition", "C", "B"))
dlt <- (rCA$log2FoldChange - rBA$log2FoldChange) - rCB$log2FoldChange
zeroBC <- rowSums(counts(dds)[, dds$condition != "A"]) == 0 & !mcols(dds)$allZero
cat("MLE: max|(C-A)-(B-A)-(C vs B)| =", signif(max(abs(dlt), na.rm = TRUE), 4),
    "| n(|diff|>1e-6):", sum(abs(dlt) > 1e-6, na.rm = TRUE), "| n(B,C all zero):", sum(zeroBC),
    "| max excluding them:", signif(max(abs(dlt[!zeroBC]), na.rm = TRUE), 3), "\n")
g0 <- which(zeroBC)[1]
cat(rownames(dds)[g0], "counts:", counts(dds)[g0, ], "| (C-A)-(B-A):", round(dlt[g0], 4),
    "| results(contrast) LFC/stat/pvalue:", unlist(as.data.frame(rCB[g0, c("log2FoldChange", "stat", "pvalue")])), "\n")
sB  <- lfcShrink(dds, coef = "condition_B_vs_A", type = "normal", quiet = TRUE)
sC  <- lfcShrink(dds, coef = "condition_C_vs_A", type = "normal", quiet = TRUE)
sCB <- lfcShrink(dds, contrast = c("condition", "C", "B"), type = "normal", quiet = TRUE)
cat("betaPriorVar coef path (standard MM):", signif(priorInfo(sB)$betaPriorVar, 4),
    "| same in sC:", identical(priorInfo(sB)$betaPriorVar, priorInfo(sC)$betaPriorVar), "\n")
cat("betaPriorVar contrast path (expanded MM):", signif(priorInfo(sCB)$betaPriorVar, 4), "\n")
fitS <- nbinomWaldTest(dds, betaPrior = TRUE, betaPriorVar = priorInfo(sB)$betaPriorVar,
                       modelMatrixType = "standard", quiet = TRUE)
cat("one joint MAP fit: max|sB - fit B| =", max(abs(sB$log2FoldChange - mcols(fitS)$condition_B_vs_A), na.rm = TRUE),
    "| max|sC - fit C| =", max(abs(sC$log2FoldChange - mcols(fitS)$condition_C_vs_A), na.rm = TRUE), "\n")
dd <- (sC$log2FoldChange - sB$log2FoldChange) - sCB$log2FoldChange
cat("shrunk: max|(sC-sB) - shrink(C vs B)| =", round(max(abs(dd), na.rm = TRUE), 4),
    "| excluding B,C all-zero:", round(max(abs(dd[!zeroBC]), na.rm = TRUE), 4),
    "| median:", round(median(abs(dd), na.rm = TRUE), 4), "\n")
tab <- data.frame(baseMean = rCB$baseMean, MLE_CvsB = rCB$log2FoldChange, sC_minus_sB = sC$log2FoldChange - sB$log2FoldChange,
                  shrink_CvsB = sCB$log2FoldChange, row.names = rownames(dds))
print(round(tab[order(-abs(dd))[1:2], ], 3))
dds2 <- dds; dds2$condition <- relevel(dds2$condition, "B"); dds2 <- nbinomWaldTest(dds2, quiet = TRUE)
sCB2 <- lfcShrink(dds2, coef = "condition_C_vs_B", type = "normal", quiet = TRUE)
cat("relevel(B) coef path betaPriorVar:", signif(priorInfo(sCB2)$betaPriorVar, 4),
    "| max|relevel coef - contrast shrink| =", round(max(abs(sCB2$log2FoldChange - sCB$log2FoldChange), na.rm = TRUE), 3),
    "| max|relevel coef - (sC-sB)| =", round(max(abs(sCB2$log2FoldChange - (sC$log2FoldChange - sB$log2FoldChange)), na.rm = TRUE), 3), "\n")
sT <- lfcShrink(dds, coef = "condition_B_vs_A", type = "normal", lfcThreshold = 1, quiet = TRUE)
cat("lfcThreshold=1 (normal): stat identical to results()?", identical(sT$stat, rBA$stat),
    "| pvalue identical?", identical(sT$pvalue, rBA$pvalue), "\n")
```

```
MLE: max|(C-A)-(B-A)-(C vs B)| = 0.03896 | n(|diff|>1e-6): 14 | n(B,C all zero): 14 | max excluding them: 8.88e-16
gene357 counts: 0 5 0 0 0 0 0 0 0 | (C-A)-(B-A): -0.039 | results(contrast) LFC/stat/pvalue: 0 0 1
betaPriorVar coef path (standard MM): 1e+06 0.8 1.695 | same in sC: TRUE
betaPriorVar contrast path (expanded MM): 1e+06 1.195 1.195 1.195
one joint MAP fit: max|sB - fit B| = 0 | max|sC - fit C| = 0
shrunk: max|(sC-sB) - shrink(C vs B)| = 0.6911 | excluding B,C all-zero: 0.6911 | median: 0.0504
         baseMean MLE_CvsB sC_minus_sB shrink_CvsB
gene158     4.746    4.266       0.658       1.349
gene1242    2.064    1.824      -0.094       0.418
relevel(B) coef path betaPriorVar: 1e+06 0.8 1.089 | max|relevel coef - contrast shrink| = 1.239 | max|relevel coef - (sC-sB)| = 1.278
lfcThreshold=1 (normal): stat identical to results()? FALSE | pvalue identical? FALSE
```

The MLE contrast is exactly additive. `getContrast()` calls `fitBeta` with `maxitSEXP = 0` and does not refit the coefficients, and leaving out the 14 genes where B and C are all 0, the difference is 8.9e-16. The 0.039 comes entirely from those 14. `cleanContrast()` overwrites the LFC with 0, stat with 0 and p with 1 for genes where both groups compared are 0 (gene357: counts 0 5 0 0 0 0 0 0 0), while the coefficient path's (C−A)−(B−A) stays at −0.039.

On the shrinkage side: with `type="normal"`, the two calls `coef=B` and `coef=C` do the same joint MAP fit, with the same standard model matrix and the same prior variances (0.80, 1.695), and just take out different columns (each differs from the joint fit by 0). So sC − sB is exactly the contrast of that joint posterior mode, $c^T\hat\beta_{MAP}$. The subtraction itself does not break. The difference of up to 0.69 arises because the priors differ: the `contrast=` path refits in the expanded model matrix with the same variance, 1.195, for every level. For the low-expression gene158, sC − sB is 0.66 and the contrast shrink 1.35. Leaving out the genes where B and C are all 0, the maximum difference is still 0.69. The `coef` path after relevel(B) has yet different prior variances (0.8, 1.089) and matches neither (maximum differences 1.24, 1.28). apeglm re-estimates the posterior for only one `coef` per call (help page: "re-estimates posterior LFCs for the coefficient specified by coef"), so the two calls are different posteriors, and their difference is not the mode of any single posterior (not run, since it is not installed). The study text's conclusion, "subtracting two shrunken coefficients is not the proper posterior shrinkage of the contrast", is right. But for normal the reason is not the nonlinearity of the posterior mode but which parametrization the prior was put on. Record which path was used with `priorInfo()`.

Giving `lfcThreshold=1` together with `type="normal"` changed `stat` and `pvalue` too. It is a new threshold test with the shrunk coefficients, an exception to "the pvalue is always the MLE test".

</details>

<details>
<summary>Experiment 3. Cook's distance: six cases that depend on the design and the number of samples</summary>

```r
suppressPackageStartupMessages(library(DESeq2))
fitWith <- function(dds, gene, sample, value, quiet = TRUE) {
  counts(dds)[gene, sample] <- as.integer(value); DESeq(dds, quiet = quiet) }
report <- function(dds, gene = "gene1", ...) {
  X <- attr(dds, "dispModelMatrix"); r <- results(dds, ...)[gene, ]
  cat(sprintf("  %s counts: %s | cooks: %s | maxCooks: %s | cutoff qf(.99,%d,%d)=%.2f\n  baseMean %.1f LFC %.3f pvalue %s\n",
      gene, paste(counts(dds)[gene, ], collapse = " "), paste(round(assays(dds)[["cooks"]][gene, ], 2), collapse = " "),
      round(mcols(dds)[gene, "maxCooks"], 2), ncol(X), nrow(X) - ncol(X), qf(0.99, ncol(X), nrow(X) - ncol(X)),
      r$baseMean, r$log2FoldChange, signif(r$pvalue, 3)))
}
cat("=== (a) 3 vs 3, gene1 sample 1 <- 5000\n")
set.seed(1); a <- fitWith(makeExampleDESeqDataSet(n = 1000, m = 6), "gene1", 1, 5000)
report(a); cat("  cooksCutoff=FALSE -> pvalue", signif(results(a, cooksCutoff = FALSE)["gene1", "pvalue"], 3), "\n")
cat("=== (a2) 3 vs 3, 2-group dontFilter exception\n")
set.seed(1); e <- makeExampleDESeqDataSet(n = 1000, m = 6)
counts(e)["gene1", ] <- c(400L, 10L, 12L, 1000L, 1100L, 950L)   # 400: a high value within A
counts(e)["gene2", ] <- c(5000L, 10L, 12L, 20L, 22L, 18L)       # control: no count larger than 5000
e <- DESeq(e, quiet = TRUE)
for (g in c("gene1", "gene2")) {
  report(e, g); flagged <- counts(e)[g, which.max(assays(e)[["cooks"]][g, ])]
  cat("  n(count > flagged count", flagged, "):", sum(counts(e)[g, ] > flagged),
      "| pvalue with cooksCutoff=FALSE:", signif(results(e, cooksCutoff = FALSE)[g, "pvalue"], 3), "\n")
}
cat("=== (b) 7 vs 7 (minReplicatesForReplace = 7)\n")
set.seed(1); b <- fitWith(makeExampleDESeqDataSet(n = 1000, m = 14), "gene1", 1, 5000, quiet = FALSE)
cat("  assays:", assayNames(b), "| all replaceable:", all(b$replaceable), "\n")
cat("  counts(b) gene1:    ", counts(b)["gene1", ], "\n")
cat("  replaceCounts gene1:", assays(b)[["replaceCounts"]]["gene1", ], "\n")
cat("  as.integer(mean(normalized, trim = 0.2) * sf[1]) =",
    as.integer(mean(counts(b, normalized = TRUE)["gene1", ], trim = 0.2) * sizeFactors(b)[1]), "\n")
report(b)
cat("=== (c) ~ pair + condition, 3 pairs x 2\n")
set.seed(1); p0 <- makeExampleDESeqDataSet(n = 500, m = 6)
p0$pair <- factor(rep(1:3, 2)); design(p0) <- ~ pair + condition
cc <- fitWith(p0, "gene1", 1, 5000)
cat("  nOrMoreInCell(X, 3):", DESeq2:::nOrMoreInCell(attr(cc, "modelMatrix"), 3),
    "| maxCooks all NA:", all(is.na(mcols(cc)$maxCooks)), "\n"); report(cc)
cat("=== (d) 2 vs 2\n")
set.seed(1); d4 <- fitWith(makeExampleDESeqDataSet(n = 500, m = 4), "gene1", 1, 5000)
cat("  maxCooks all NA:", all(is.na(mcols(d4)$maxCooks)), "\n"); report(d4)
cat("=== (e) 3 vs 2 (A A A B B)\n")
set.seed(1); e0 <- makeExampleDESeqDataSet(n = 500, m = 5)
for (s in c(1, 4)) {
  e5 <- fitWith(e0, "gene1", s, 5000)
  cat("  inject into sample", s, "| nOrMoreInCell(X, 3):", DESeq2:::nOrMoreInCell(attr(e5, "modelMatrix"), 3),
      "| maxCooks all NA:", all(is.na(mcols(e5)$maxCooks)), "\n"); report(e5)
}
```

```
=== (a) 3 vs 3, gene1 sample 1 <- 5000
  gene1 counts: 5000 1 1 4 3 14 | cooks: 36.93 9.22 9.21 0.37 0.6 1.97 | maxCooks: 36.93 | cutoff qf(.99,2,4)=18.00
  baseMean 783.8 LFC -7.875 pvalue NA
  cooksCutoff=FALSE -> pvalue 0.00177
=== (a2) 3 vs 3, 2-group dontFilter exception
  gene1 counts: 400 10 12 1000 1100 950 | cooks: 26.85 6.82 6.51 0.01 0.13 0.08 | maxCooks: 26.85 | cutoff qf(.99,2,4)=18.00
  baseMean 556.8 LFC 2.895 pvalue 0.0449
  n(count > flagged count 400 ): 3 | pvalue with cooksCutoff=FALSE: 0.0449
  gene2 counts: 5000 10 12 20 22 18 | cooks: 36.44 9.12 9.09 0 0.07 0.07 | maxCooks: 36.44 | cutoff qf(.99,2,4)=18.00
  baseMean 793.1 LFC -6.342 pvalue NA
  n(count > flagged count 5000 ): 0 | pvalue with cooksCutoff=FALSE: 0.00176
=== (b) 7 vs 7 (minReplicatesForReplace = 7)
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
-- replacing outliers and refitting for 1 genes
-- DESeq argument 'minReplicatesForReplace' = 7
-- original counts are preserved in counts(dds)
estimating dispersions
fitting model and testing
  assays: counts mu H cooks replaceCounts replaceCooks | all replaceable: TRUE
  counts(b) gene1:     5000 1 1 4 3 14 1 6 1 3 2 10 6 21
  replaceCounts gene1: 5 1 1 4 3 14 1 6 1 3 2 10 6 21
  as.integer(mean(normalized, trim = 0.2) * sf[1]) = 5
  gene1 counts: 5000 1 1 4 3 14 1 6 1 3 2 10 6 21 | cooks: 83.69 2.34 2.34 2.32 2.32 2.26 2.34 0.01 0.39 0.17 0.26 0.09 0.01 2.11 | maxCooks: NA | cutoff qf(.99,2,12)=6.93
  baseMean 5.3 LFC 0.753 pvalue 0.324
=== (c) ~ pair + condition, 3 pairs x 2
  nOrMoreInCell(X, 3): FALSE FALSE FALSE FALSE FALSE FALSE | maxCooks all NA: TRUE
  gene1 counts: 5000 1 1 9 3 3 | cooks: 0.23 0.07 0.05 0.23 0.06 0.05 | maxCooks: NA | cutoff qf(.99,4,2)=99.25
  baseMean 793.2 LFC -0.384 pvalue 0.855
=== (d) 2 vs 2
  maxCooks all NA: TRUE
  gene1 counts: 5000 1 1 9 | cooks: 0.22 0.22 0.13 0.13 | maxCooks: NA | cutoff qf(.99,2,2)=99.00
  baseMean 1173.3 LFC -8.866 pvalue 0.00451
=== (e) 3 vs 2 (A A A B B)
  inject into sample 1 | nOrMoreInCell(X, 3): TRUE TRUE TRUE FALSE FALSE | maxCooks all NA: FALSE
  gene1 counts: 5000 1 1 9 3 | cooks: 36.91 9.25 9.2 1.28 1.33 | maxCooks: 36.91 | cutoff qf(.99,2,3)=30.82
  baseMean 973.3 LFC -8.128 pvalue NA
  inject into sample 4 | nOrMoreInCell(X, 3): TRUE TRUE TRUE FALSE FALSE | maxCooks all NA: FALSE
  gene1 counts: 4 1 1 5000 3 | cooks: 0.7 0.22 0.14 24.7 24.71 | maxCooks: 0.7 | cutoff qf(.99,2,3)=30.82
  baseMean 976.3 LFC 10.291 pvalue 1.59e-05
```

(a) In 3 vs 3, Cook's 36.9 > 18, so only p is NA. The LFC of −7.9 remains, so "an LFC but p NA" is the trace of Cook's. With `cooksCutoff=FALSE`, p is 0.00177.

(a2) is the 2-group exception. gene1's 400 jumps within A, so it was flagged with Cook's 26.85 > 18, but the three values of B, 1000, 1100 and 950, are larger than 400, so it stayed via `dontFilter` (p 0.0449, the same as with `cooksCutoff=FALSE`). gene2, added as a control, has no count larger than 5000 and was filtered (p NA; 0.00176 with the filter off).

(b) In 7 vs 7, 5000 was replaced by 5, a value based on the trimmed mean, and refitted. The original count is in `counts(dds)` and the replacement in `assays(dds)[["replaceCounts"]]`. The `cooks` assay shows the value before replacement (83.69), but every sample is replaceable, so `maxCooks` was changed to NA and there is no Cook's filter afterwards.

(c) In the paired design, all six rows of the model matrix differ, so 1 per cell, and `maxCooks` is all NA. Even with 5000 in it, nothing is done. The pair coefficient absorbs that sample's count, so Cook's itself is small, 0.23, and the cutoff is high too, $F_{0.99}(4,2) = 99.25$.

(d) 2 vs 2 also has `maxCooks` NA, so p = 0.00451 is reported as it is, with no Cook's filter.

(e) In 3 vs 2, only the A cell (3 samples) has `nOrMoreInCell` TRUE. Putting 5000 into a sample of A gives `maxCooks` 36.91 > 30.82, and p becomes NA. Putting it into a sample of B leaves out the Cook's of the two B samples (24.7), `maxCooks` becomes A's maximum of 0.7, and p = 1.59e-05 is reported. In fact this 24.7 is below the cutoff of 30.82, so it would not have been filtered even if included. What this run shows is that in a design with mixed cells, `maxCooks` is computed only from the samples of cells with 3 or more.

</details>

<details>
<summary>Experiment 4. Counting each row of the NA table (continuing the session of experiment 1)</summary>

```r
naP <- is.na(res$pvalue); naQ <- is.na(res$padj)
g0 <- rownames(res)[which(res$baseMean == 0)[1]]
cat("allZero (baseMean==0):", sum(res$baseMean == 0), "|", g0, ":"); print(unlist(as.data.frame(res[g0, ])))
cat("pvalue NA & baseMean>0 (Cook's):", sum(naP & res$baseMean > 0),
    "| max(maxCooks):", round(max(mcols(dds)$maxCooks, na.rm = TRUE), 2), "< cutoff", qf(0.99, 2, 4), "\n")
r2 <- results(dds, independentFiltering = FALSE)
cat("pvalue ok & padj NA (independent filtering):", sum(!naP & naQ),
    "| same with independentFiltering=FALSE:", sum(!is.na(r2$pvalue) & is.na(r2$padj)), "\n")
thr <- metadata(res)$filterThreshold
cat("filterThreshold:", round(thr, 4), "| filterTheta:", round(metadata(res)$filterTheta, 3),
    "| max baseMean (padj-NA-only):", round(max(res$baseMean[!naP & naQ]), 4),
    "| padj-NA-only with baseMean >= threshold:", sum(!naP & naQ & res$baseMean >= thr), "\n")
cat("betaConv FALSE:", sum(!mcols(dds)$betaConv, na.rm = TRUE), "| betaConv NA:", sum(is.na(mcols(dds)$betaConv)),
    "| NA rows == allZero rows:", identical(is.na(mcols(dds)$betaConv), mcols(dds)$allZero), "\n")
```

```
allZero (baseMean==0): 16 | gene189 :      baseMean log2FoldChange          lfcSE           stat         pvalue           padj
             0             NA             NA             NA             NA             NA
pvalue NA & baseMean>0 (Cook's): 0 | max(maxCooks): 17.54 < cutoff 18
pvalue ok & padj NA (independent filtering): 579 | same with independentFiltering=FALSE: 0
filterThreshold: 4.4313 | filterTheta: 0.198 | max baseMean (padj-NA-only): 4.4307 | padj-NA-only with baseMean >= threshold: 0
betaConv FALSE: 0 | betaConv NA: 16 | NA rows == allZero rows: TRUE
```

The three kinds of NA leave different column patterns. In this data there were no Cook's outliers (maximum 17.54 < 18). There were no convergence failures either (the 16 NA in `betaConv` are the allZero genes), but in the source a failure does not turn into NA, so `mcols(dds)$betaConv` has to be checked directly.

</details>

<details>
<summary>Experiment 5. Does blind= remove batch?</summary>

In a balanced design (2 each of condition × batch, m=8), I put a 2-fold effect of batch b2 into 600 genes, ran `DESeq()` with `design = ~batch + condition`, and computed the transformation in several ways.

```r
suppressPackageStartupMessages(library(DESeq2))
set.seed(4)
dds <- makeExampleDESeqDataSet(n = 2000, m = 8, betaSD = 0.5)
dds$batch <- factor(rep(c("b1", "b2"), 4))              # 2 each of condition x batch (balanced)
bf <- rep(1, nrow(dds)); bf[1:600] <- 2                  # genes 1..600 are 2x in b2
cts <- counts(dds); j <- dds$batch == "b2"
cts[, j] <- matrix(rpois(sum(j) * nrow(dds), lambda = cts[, j] * bf), nrow = nrow(dds))
counts(dds) <- cts; design(dds) <- ~ batch + condition
dds <- DESeq(dds, quiet = TRUE)
pca <- function(v) { pd <- plotPCA(v, intgroup = c("batch", "condition"), returnData = TRUE)
  cat(sprintf("percentVar PC1/PC2: %.3f %.3f | PC1 mean by condition: A %.2f, B %.2f | PC2 mean by batch: b1 %.2f, b2 %.2f\n",
      attr(pd, "percentVar")[1], attr(pd, "percentVar")[2], tapply(pd$PC1, pd$condition, mean)[1],
      tapply(pd$PC1, pd$condition, mean)[2], tapply(pd$PC2, pd$batch, mean)[1], tapply(pd$PC2, pd$batch, mean)[2])) }
vsd_T <- vst(dds, blind = TRUE); vsd_F <- vst(dds, blind = FALSE)
vsd_V <- varianceStabilizingTransformation(dds, blind = FALSE)          # study text's example function: reuses the trend of DESeq()
t_rlog <- system.time(rld_F <- rlog(dds, blind = FALSE))["elapsed"]
cat("design(dds) still:", deparse(design(dds)), "\n")
cat("max|vst(T) - vst(F)| =", round(max(abs(assay(vsd_T) - assay(vsd_F))), 4),
    "| max|vst(F) - VST(F)| =", round(max(abs(assay(vsd_F) - assay(vsd_V))), 4),
    "| max|rlog(F) - VST(F)| =", round(max(abs(assay(rld_F) - assay(vsd_V))), 4), "| rlog elapsed (s):", t_rlog, "\n")
cat("vst  blind=TRUE   "); pca(vsd_T)
cat("vst  blind=FALSE  "); pca(vsd_F)
cat("VST  blind=FALSE  "); pca(vsd_V)
cat("rlog blind=FALSE  "); pca(rld_F)
cat("DESeqDataSetFromMatrix(assay(vsd_F), ...):",
    tryCatch(DESeqDataSetFromMatrix(assay(vsd_F), colData(dds), ~ condition), error = conditionMessage), "\n")
vsd_R <- vsd_F                                                       # copy for visualization only
assay(vsd_R) <- limma::removeBatchEffect(assay(vsd_F), batch = vsd_F$batch,
                                         design = model.matrix(~ condition, colData(vsd_F)))
cat("removeBatchEffect "); pca(vsd_R)
pcT <- plotCounts(dds, "gene1", intgroup = c("condition", "batch"), returnData = TRUE)
pcF <- plotCounts(dds, "gene1", intgroup = c("condition", "batch"), returnData = TRUE, transform = FALSE)
nc <- counts(dds, normalized = TRUE)["gene1", ]
cat("plotCounts count - normalized count: transform=TRUE", unique(round(pcT$count - nc, 6)),
    "| transform=FALSE", unique(round(pcF$count - nc, 6)), "\n")
cf <- dds; cf$batch <- factor(ifelse(cf$condition == "A", "b1", "b2"))   # batch = condition (complete confounding)
cat("confounded ~ batch + condition:", tryCatch({ DESeq(cf, quiet = TRUE); "fitted" }, error = conditionMessage), "\n")
```

```
design(dds) still: ~batch + condition
max|vst(T) - vst(F)| = 0.9555 | max|vst(F) - VST(F)| = 0.2894 | max|rlog(F) - VST(F)| = 7.165 | rlog elapsed (s): 0.11
vst  blind=TRUE   using ntop=500 top features by variance
percentVar PC1/PC2: 0.228 0.210 | PC1 mean by condition: A -7.01, B 7.01 | PC2 mean by batch: b1 -6.85, b2 6.85
vst  blind=FALSE  using ntop=500 top features by variance
percentVar PC1/PC2: 0.238 0.229 | PC1 mean by condition: A -5.62, B 5.62 | PC2 mean by batch: b1 5.58, b2 -5.58
VST  blind=FALSE  using ntop=500 top features by variance
percentVar PC1/PC2: 0.240 0.233 | PC1 mean by condition: A 4.94, B -4.94 | PC2 mean by batch: b1 4.91, b2 -4.91
rlog blind=FALSE  using ntop=500 top features by variance
percentVar PC1/PC2: 0.240 0.226 | PC1 mean by condition: A -6.17, B 6.17 | PC2 mean by batch: b1 6.02, b2 -6.02
DESeqDataSetFromMatrix(assay(vsd_F), ...): some values in assay are not integers
removeBatchEffect using ntop=500 top features by variance
percentVar PC1/PC2: 0.308 0.158 | PC1 mean by condition: A -6.02, B 6.02 | PC2 mean by batch: b1 0.00, b2 0.00
plotCounts count - normalized count: transform=TRUE 0.5 | transform=FALSE 0
confounded ~ batch + condition: full model matrix is less than full rank
```

Even with `blind=FALSE`, PC2 splits by batch at ±5.6, because `blind` only changes which design the dispersion trend is estimated with. The study text's example function `varianceStabilizingTransformation(blind=FALSE)` reuses the trend of `DESeq()`, while `vst(blind=FALSE)` refits the trend on a subset, so the values differ by up to 0.289. In both, PC2 splits by batch (±5.6, ±4.9). `rlog(blind=FALSE)` also split PC2 at ±6.0. rlog is a different transformation from VST, so the values themselves differ a lot (maximum difference 7.17). At this size (2000 × 8), rlog took just over 0.1 seconds, so no cost difference showed; I did not check its speed with many samples.

Only in the copy with `removeBatchEffect` applied did the batch means on PC2 become 0. This matrix is not integer and cannot be test input; more fundamentally, it has lost the mean–variance structure of the raw counts, so it is used only for plotting. The test is done with the raw counts and `~ batch + condition`. If batch and condition are changed to overlap completely, `DESeq()` rejects it with "full model matrix is less than full rank", and in such a design even a clean split in the PCA cannot tell a condition effect from a batch effect.

The `count` column of `plotCounts(..., returnData=TRUE)` is normalized count + 0.5 with `transform=TRUE` (the default), and the normalized count as it is with `transform=FALSE` (source: `pc <- if (transform) 0.5 else 0`). Don't forget this 0.5 when copying values into a table.

</details>

<details>
<summary>Experiment 6. How does the s-value differ from padj? (An illustrative calculation, not DESeq2 output)</summary>

Without apeglm/ashr, real s-values cannot be computed. Instead, from the `normal` result of experiment 1, I approximated the posterior as $N(\hat\beta_{shr}, \mathrm{lfcSE}^2)$ and reproduced only the definitions of the lfsr and s-value (continuing the session of experiment 1).

```r
lfsr <- pnorm(-abs(shr$log2FoldChange) / shr$lfcSE)   # approximate the posterior as N(LFC_shr, lfcSE^2) (illustration)
svalue <- function(l) { o <- order(l); out <- numeric(length(l)); out[o] <- cumsum(l[o]) / seq_along(o); out }
ok <- !is.na(lfsr); sv <- rep(NA_real_, length(lfsr)); sv[ok] <- svalue(lfsr[ok])
g <- c("gene1349", "gene2613", "gene3"); i <- match(g, rownames(res))
print(signif(data.frame(d[g, c("baseMean", "LFC_MLE", "LFC_shr", "SE_shr", "pvalue", "padj")],
                        lfsr = lfsr[i], svalue = sv[i]), 4))
cat("n(padj<0.05):", sum(res$padj < 0.05, na.rm = TRUE), "| n(svalue<0.05):", sum(sv < 0.05, na.rm = TRUE), "\n")
cat("Spearman cor(pvalue, svalue):", round(cor(res$pvalue, sv, method = "spearman", use = "complete.obs"), 3),
    "| padj<0.05 but svalue>=0.05:", sum(res$padj < 0.05 & sv >= 0.05, na.rm = TRUE),
    "| svalue<0.05 but (padj>=0.05 or NA):", sum(sv < 0.05 & (is.na(res$padj) | res$padj >= 0.05), na.rm = TRUE), "\n")
```

```
         baseMean LFC_MLE LFC_shr SE_shr    pvalue     padj      lfsr    svalue
gene1349    4.886  -5.808 -1.9690 0.6266 0.0006278 0.008435 0.0008396 0.0001459
gene2613    3.950   5.457  0.3516 0.3931 0.1461000       NA 0.1856000 0.0600100
gene3       4.533  -1.375 -0.6348 0.6295 0.3109000 0.557200 0.1566000 0.0495500
n(padj<0.05): 313 | n(svalue<0.05): 1479
Spearman cor(pvalue, svalue): 0.999 | padj<0.05 but svalue>=0.05: 0 | svalue<0.05 but (padj>=0.05 or NA): 1166
```

Cutting both at the same 0.05 gives different numbers, 313 vs 1479. But this example only shows the formulas; it is no evidence about how different the two quantities are in size. To begin with, the `lfcSE` used here is a sandwich SE, smaller than the posterior SD (experiment 1: gene1349 0.63 vs 0.94). The lfsr and s-values come out correspondingly smaller, so 1479 is inflated. Also, in this approximation the s-value ranking is almost the same as the p ranking (Spearman 0.999; 0 genes with padj < 0.05 but s-value ≥ 0.05). So 313 vs 1479 is a difference not of ranking but of the scale on which they are cut. The difference in meaning lies in the definitions. The s-value is "the average probability of a wrong sign in the list of genes with equal or smaller s-values", and padj is "the smallest level that controls the expected proportion (FDR) of true nulls in the list rejected up to this gene". Real apeglm/ashr s-values come from each method's posterior, so they differ from the numbers above.

</details>

<details>
<summary>Common misconceptions, with evidence</summary>

- Does lfcShrink "correct" the p-values too? Not in the default call. With normal, `pvalue` and `padj` are identical down to the values (experiment 1), and apeglm/ashr keep the MLE `pvalue` and `padj` as they are and just drop `stat` (source). The exceptions differ by type. normal + `lfcThreshold>0` replaces stat/pvalue/padj with the threshold test of the shrunken fit (experiment 2). apeglm removes pvalue/padj and adds s-values with `svalue=TRUE` or `lfcThreshold>0`. ashr removes pvalue/padj with `svalue=TRUE` (`lfcThreshold=0`), but with `lfcThreshold>0` it keeps the MLE pvalue/padj and appends the FSOS s-value (source).
- Is the shrunken LFC / lfcSE the Wald Z? No. It gives a different number (−3.142 vs −3.419). normal's lfcSE is the sandwich SE of the MAP, and apeglm/ashr's is the posterior SD, so neither is the MLE Wald SE.
- Is normal's lfcSE a posterior SD too? The help page says so, but the implemented value is a sandwich SE, and for low-expression genes it is much smaller than the Laplace posterior SD (0.63 vs 0.94, experiment 1).
- Are all lfcShrink results MAP? ashr is `PosteriorMean` (source), labelled "MMSE".
- Does shrinking two `coef`s and subtracting them equal the contrast shrink? They differ by up to 0.69, and for the low-expression gene158 it is 0.66 vs 1.35 (experiment 2). With normal the two calls are the same joint MAP, so the subtraction itself is right, but the `contrast` path has a different prior (expanded MM, 1.195). apeglm refuses the `contrast` argument altogether.
- Does the MLE too drift slightly between the coefficient difference and the contrast because of refitting? There is no refit (`maxit = 0`), and it is exactly additive (8.9e-16). The 14 that differ are genes where both groups compared are 0, which `results(contrast=)` overwrote with LFC 0, p 1 (experiment 2).
- Does Cook's outlier handling always happen? Only samples in cells with 3 or more go into the call. Without such samples (2 vs 2, or a paired design with 1 per cell), `maxCooks` is all NA, and with 3 vs 2 only the group of 3 is judged. If a cell has 7 or more, replacement happens, and if every sample is replaceable, the filter is turned off afterwards (experiment 3).
- In 2 groups, are only low outliers let off? The condition is "3 or more counts larger than the flagged count, across both groups". A value high within its own group also stays if the other group is larger (experiment 3-(a2)).
- If p is NA, is it "not expressed"? The Cook's outlier was a gene with baseMean 783. Tell the causes of NA apart by the column pattern (experiment 4).
- If padj is NA, is it "no difference"? Only padj NA means it was removed from the test family by independent filtering, and the Cook's gene with p NA too had LFC −7.9 (experiment 3-(a)). When counting NA as a non-discovery, record that rule and the cause.
- Does `blind=FALSE` give a PCA with batch removed? PC2 still splits by batch (experiment 5). The design is used only for dispersion estimation.
- Does `vst(blind=FALSE)` use the trend of `DESeq()` as it is? `vst()` refits the trend on a subset. The one that uses it as it is is `varianceStabilizingTransformation(blind=FALSE)`, and the two differ by up to 0.289 (experiment 5).
- Can VST values be fed into DESeq2 for testing? The NB GLM models the mean–variance relationship of the raw counts and the size factor offset, so transformed values cannot be its input. `DESeqDataSetFromMatrix` rejecting them as non-integer is only a symptom.
- If the PCA splits well, is the design fine too? A completely confounded batch is rejected by `DESeq()` with "full model matrix is less than full rank" (experiment 5). Whether it splits has to be looked at together with quality, design and model diagnostics.

</details>

<details>
<summary>Where the results differ from the study text</summary>

Only the places where the study text's explanation and the actual behaviour of DESeq2 differ, or where details not in the study text came out, are collected here. The rest of the study text's explanations matched the actual behaviour of DESeq2.

| Study text claim (chapter 13) | Result | Evidence |
|---|---|---|
| LFC shrinkage is a separate step after `DESeq()` | Matches (with an exception path) | `lfcShrink` lines 15–19: stops if `resultsNames` is empty, except for apeglm + `apeAdapt=FALSE`. Lines 23–26: always stops if `betaPrior=TRUE` |
| apeglm = heavy-tailed prior, posterior mode | Matches (help page and source) / details of the prior's form not confirmed | Help page "adaptive Student's t prior", reference title "Heavy-tailed prior distributions"; `fit$map`, `sub("MLE","MAP")`. Whether it is Cauchy not confirmed, apeglm not installed |
| ashr = posterior mean | Matches (source) | `fit$result$PosteriorMean`, label "MMSE"; cannot be run, not installed |
| apeglm supports only coefficients | Matches (source) / not observed by running | `stop("type='apeglm' shrinkage only for use with 'coef'")` (lines 170–171). In this environment, calling `contrast` + apeglm hits the installation check on lines 167–168 first, so only the installation message appears (the "Calling the study text's example code as it is" block) |
| The difference of two shrunken coefficients ≠ the posterior shrinkage of the contrast | Matches (conclusion) / reason supplemented | Experiment 2: up to 0.69, gene158 0.66 vs 1.35. normal's two coef calls are the same joint MAP (difference 0), so the subtraction is preserved, and the difference comes from the different prior of the contrast path (expanded MM 1.195 vs standard 0.80/1.695). apeglm is a different posterior per call (not run) |
| Re-setting Starvation as the reference and refitting makes it a coef | Matches (but under normal the values differ from the contrast path) | `condition_C_vs_B` appears after `relevel` + `nbinomWaldTest`; the prior variances (0.8, 1.089) differ, so up to 1.24 from the contrast path (experiment 2) |
| The `lfcSE` of `lfcShrink()` is a posterior SD | Matches for apeglm/ashr (source) / does not match for normal | apeglm `fit$sd`, ashr `PosteriorSD`, label "posterior SD". normal is a sandwich SE: gene1349 lfcSE 0.6266 = sandwich 0.6266 ≠ Laplace posterior SD 0.9413 (experiment 1), and its label stays "standard error". The help page says posterior SD for all three types |
| pvalue and padj remain | Matches (default) / exceptions by type | normal's default gives identical values (experiment 1); normal + `lfcThreshold>0` replaces them (experiment 2); apeglm removes them with `svalue=TRUE` or `lfcThreshold>0`; ashr removes them with `svalue=TRUE` (`lfcThreshold=0`) and keeps them with `lfcThreshold>0` (source) |
| The s-value and the BH padj mean different things | Matches (conceptual) | Help page "probability of false signs among the tests with equal or smaller s-value"; experiment 6 only reproduces the definition (rank correlation 0.999) |
| Example code `varianceStabilizingTransformation(dds, blind=FALSE)` | Matches (behaviour) / a caution | Reuses the trend of `DESeq()`. `vst(blind=FALSE)` refits the trend on a subset and differs by up to 0.289 (experiment 5) |
| Automatic replacement does not always happen in a paired design | Matches (more strongly: a paired design with 1 per cell has no Cook's filter at all) | Experiment 3-(c): `nOrMoreInCell(X,3)` all FALSE, `maxCooks` all NA |
| NA table: baseMean=0 → p, padj NA | Matches (LFC, lfcSE and stat NA too) | Experiment 4, gene189 |
| NA table: convergence warning → check β convergence and the model matrix | Matches (columns to check) / causes not confirmed | The `betaConv` column exists (experiment 4, 0 non-converged). Even deliberately created non-convergence with `maxit = 2` does not make p NA (the "When even p is NA" block). The causes the study text lists (lack of information, extreme counts, design) not confirmed, since no naturally occurring non-convergence case arose |
| Appendix B 17: if p is NA too, it is all-zero, a count outlier or a fitting problem | Partly matches | All-zero and Cook's in experiments 3 and 4. "Fitting problems" split in two. Rows whose design degenerated because of weights get NA even for p (matches), while non-convergence of β does not make p NA (does not match). Both were created deliberately and run in the "When even p is NA" block, and a search of the whole source confirmed that the only line in `results()` that puts NA into p is the Cook's rule |
| Appendix B 18: the lfcSE of shrinkage is a posterior SD | Partly does not match | normal is a sandwich SE (experiment 1). The conclusions (no recomputing Z, ashr = posterior mean) match |
| (A detail not in the study text) In 2 groups, no filtering if 3 or more counts are larger than the flagged count | Confirmed (source + run) | `dontFilter` in `results` lines 225–246; experiment 3-(a2): with 400 vs 1000·1100·950, maxCooks 26.85 > 18 yet p 0.0449 kept, while the control gene2 is NA |
| Don't judge DESeq2 or the design from whether the PCA splits alone | Matches (partly run) | Experiment 5: the balanced design gives PC1=condition, PC2=batch; complete confounding is rejected by `DESeq()`. The PCA of the confounded data itself was not run |
| `type="apeglm"` in the example code (`coef="condition_Starvation_vs_Ctrl"`) | The coefficient name matches / cannot be run in this environment | The "Calling the study text's example code as it is" block (run): `condition_Starvation_vs_Ctrl` is in `resultsNames`; `requireNamespace("apeglm")` is FALSE and the call stops with `Error: type='apeglm' requires installing the Bioconductor package 'apeglm'`. Omitting `type` (default apeglm) and the `contrast` call give the same error. `type="ashr"` gives `Error: type='ashr' requires installing the CRAN package 'ashr'` for both `coef` and `contrast`. Only `type="normal"` runs (line 500) |
| (15.1) Without apeglm/ashr, optional shrinkage is skipped with a message | Partly matches | Same block as above: normal runs without apeglm/ashr, and apeglm/ashr stop with the messages above. But `lfcShrink()` `stop()`s rather than skipping, so "skips" holds only when the workbook branches with `requireNamespace()`. `DESeq2_workbook.R` is not in the repository, so that branch was not confirmed |

</details>

In the next note, [08](08_one_gene_end_to_end.md), I follow by hand the calculation in which gene A's six counts become α, SE, Z and p, and compare it with the same counts put into actual DESeq2. Then I run the three-group paired practice code of chapter 15 of the study text on data with known true effects. The conditions under which the Cook's filter and the refit seen in this note switch on (design cell size) come up again in that chapter 15 practice.

---

← Previous: [06. Multiple testing](06_multiple_testing.md) · Next: [08. One gene, end to end](08_one_gene_end_to_end.md) →

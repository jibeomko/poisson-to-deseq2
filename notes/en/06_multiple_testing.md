# 06. How should p-values be read when thousands of genes are tested at once?

A p-value for one gene, as in [note 05](05_wald_vs_lrt.md), is not hard to read, but RNA-seq tests thousands of genes at once. How far, then, can a list picked at p < 0.05 be trusted, and which genes does DESeq2's `padj` adjust together? Starting from computing BH by hand, I follow it through to three comparisons and to the rescue claim "starvation lowered it and glucose brought it back". To give the punchline first: what padj means depends on which p-values were adjusted together as one set, so you have to decide first what you want to claim.

> Study text chapters 11 and 14 · Code: [06_multiple_testing.R](../../06_multiple_testing.R)

All code, including what is inside the collapsible blocks, was run from top to bottom in a single R session, and later blocks use objects from earlier blocks (`dds`, `res_list`, `lfc` and so on) as they are.

## 1. If genes are picked at p < 0.05, how many of them are fake?

The Wald test p-value of gene A (Ctrl 100, 130, 90 / Starvation 200, 250, 180) found in [note 05](05_wald_vs_lrt.md) was about 0.00072. A p-value is "the probability, assuming there is really no difference (the null hypothesis), of a result as extreme as the current one or more arising by chance". Looking at one gene alone, 0.00072 is fairly strong evidence.

But RNA-seq tests 20000 genes at once. Even if none of the 20000 changed, on average 20000 × 0.00072244 ≈ 14.4 genes would have p of 0.00072244 or less. With the threshold at p < 0.05, the average is 1000. So how many fakes are mixed into a list cannot be told from one gene's p-value.

I ran a small simulation with 1000 genes that did not change (null genes) and 100 that really changed (signal genes). The test statistic Z is "estimate ÷ SE", like the Wald Z of note 05. I made the Z of null genes follow the standard normal and the Z of signal genes a normal with mean 3, and drew a histogram of the two-sided p-values.

![p-value histogram of 1000 null and 100 signal genes](../../figures/06_pvalue_histogram.png)

*Figure 1.* The p-values of the null genes (grey) spread evenly between 0 and 1, about 50 per bin, while the signal genes (blue) pile up near 0. The p < 0.05 bin holds 85 signal and 48 null genes together.

Calling all 133 with p < 0.05 discoveries makes 48 of them, about 36%, fake. The p-values of null genes spread evenly (a uniform distribution), so about 5% of the null genes fall into any bin. So the more null genes there are, the more fakes in the p < 0.05 bin, in proportion. This is the multiple-testing problem.

<details>
<summary>Code for figure 1</summary>

The data are the same as the first repetition of the simulation below (200 repetitions).

```r
set.seed(20260926)
z <- c(rnorm(1000, 0), rnorm(100, 3)); is_sig <- rep(c(FALSE, TRUE), c(1000, 100))
p <- 2 * pnorm(-abs(z))
col_null <- "#b8b7b0"; col_sig <- "#2a78d6"; col_line <- "#eb6834"; ink <- "#52514e"
png("../figures/06_pvalue_histogram.png", 7, 4.5, "in", res = 150, bg = "white")
par(mar = c(4.2, 4.2, 2.6, 1), col.axis = ink, col.lab = ink, fg = ink, las = 1)
br <- seq(0, 1, 0.05)
h <- rbind(null = hist(p[!is_sig], br, plot = FALSE)$counts, signal = hist(p[is_sig], br, plot = FALSE)$counts)
barplot(h, space = 0, col = c(col_null, col_sig), border = "white", ylim = c(0, 140),
        xlab = "raw p-value (bins of width 0.05)", ylab = "number of genes")
axis(1, at = seq(0, 20, 4), labels = seq(0, 1, 0.2))
abline(h = 50, lty = 2, col = ink)
text(19.8, 66, "expected nulls per bin = 1000 x 0.05 = 50", adj = 1, cex = 0.8, col = ink)
text(1.3, h[1, 1] + h[2, 1] - 8, sprintf("p < 0.05: %d nulls + %d signals", h[1, 1], h[2, 1]), adj = 0, cex = 0.8, col = ink)
legend("topright", c("true null (1000 genes)", "true signal (100 genes)"), fill = c(col_null, col_sig),
       border = NA, bty = "n", cex = 0.85, inset = c(0, 0.08))
title("Simulated p-values: nulls are flat, signals pile up near 0", adj = 0, cex.main = 1, font.main = 1, col.main = "#0b0b0b")
invisible(dev.off())
cat("figure 1 written:", file.exists("../figures/06_pvalue_histogram.png"),
    "| p < 0.05: null", h[1, 1], "signal", h[2, 1], "\n")
#> figure 1 written: TRUE | p < 0.05: null 48 signal 85
```

</details>

### What do the FWER and the FDR each guard against?

With $R$ the number of genes discovered and $V$ the number of false discoveries among them (picked as discoveries though actually null), the p < 0.05 list of figure 1 has $R=133$ and $V=48$, so the proportion of false discoveries, the FDP, is 48/133 ≈ 0.36. Applying BH at the 0.05 level to the same data (section 2) gives $R=50$ and $V=2$, lowering the FDP to 0.04. This proportion $\mathrm{FDP}=V/\max(R,1)$ is a value from one experiment and differs from experiment to experiment.

$$\mathrm{FWER}=P(V\ge 1),\qquad \mathrm{FDR}=E[\mathrm{FDP}]=E\!\left[\frac{V}{\max(R,1)}\right].$$

The FWER (family-wise error rate) is the probability of even one false discovery. The Bonferroni and Holm adjustments guard against it. Bonferroni multiplies every p-value by the number of tests $m$, and Holm keeps the same FWER while reducing the multiplier by rank, $m, m-1, \dots$, so it is a little less strict. The FDR (false discovery rate) is the expected proportion of false discoveries in the list of discoveries; think of it as the average FDP if the same experiment were repeated many times. The BH (Benjamini–Hochberg) adjustment guards against this. The $\max(R,1)$ in the formula is a device to avoid dividing by 0 when there are no discoveries at all.

So FDR 5% means neither "each gene has a 5% chance of being wrong" nor a guarantee that "exactly 5 of 100 discoveries are wrong". It is an average over repetitions of the procedure. At first I read padj 0.05 as "the probability that this gene is wrong", but repeating the same simulation 200 times made the difference clear.

```r
set.seed(20260926); B <- 200; m0 <- 1000; m1 <- 100; alpha <- 0.05
fdp_bh <- fdp_bonf <- fdp_raw <- power_bh <- fwer_holm <- numeric(B)
for (b in seq_len(B)) {
  z <- c(rnorm(m0, 0), rnorm(m1, 3)); truth <- rep(c(FALSE, TRUE), c(m0, m1))   # 1000 null + 100 signal
  p <- 2 * pnorm(-abs(z))
  R  <- p.adjust(p, "BH") < alpha;         fdp_bh[b] <- sum(R & !truth) / max(sum(R), 1); power_bh[b] <- mean(R[truth])
  Rb <- p.adjust(p, "bonferroni") < alpha; fdp_bonf[b] <- sum(Rb & !truth) / max(sum(Rb), 1)
  Rh <- p.adjust(p, "holm") < alpha;       fwer_holm[b] <- any(Rh & !truth)
  Rr <- p < alpha;                         fdp_raw[b] <- sum(Rr & !truth) / max(sum(Rr), 1)
}
cat(sprintf("raw p < 0.05 : FDR = %.3f\n", mean(fdp_raw)))
cat(sprintf("BH           : FDR = %.4f  (pi0 * alpha = %.4f)  power = %.3f  P(V>=1) = %.2f\n",
            mean(fdp_bh), m0 / (m0 + m1) * alpha, mean(power_bh), mean(fdp_bh > 0)))
cat(sprintf("Bonferroni   : FDR = %.4f\n", mean(fdp_bonf)))
cat(sprintf("Holm         : P(V>=1) = %.3f\n", mean(fwer_holm)))
cat(sprintf("BH FDP over %d datasets: 5%% = %.3f  median = %.3f  95%% = %.3f  max = %.3f\n",
            B, quantile(fdp_bh, .05), median(fdp_bh), quantile(fdp_bh, .95), max(fdp_bh)))
#> raw p < 0.05 : FDR = 0.372
#> BH           : FDR = 0.0459  (pi0 * alpha = 0.0455)  power = 0.480  P(V>=1) = 0.88
#> Bonferroni   : FDR = 0.0035
#> Holm         : P(V>=1) = 0.045
#> BH FDP over 200 datasets: 5% = 0.000  median = 0.044  95% = 0.098  max = 0.130
```

The FDR of BH (average over 200) is below the target of 0.05, but the FDP of each data set swings from 0 to 0.13. A list picked at p < 0.05 with no adjustment has an FDR as high as 0.372. BH's 0.0459 is almost the same as the theoretical value $\pi_0\alpha$ = 0.0455. $\pi_0$ is the proportion of null genes in the total, here 1000/1100. Meanwhile it found 48% of the signal genes (power, the proportion of genes that really changed that were found).

But BH does not guard the FWER as well. The probability that the BH list contained at least one false discovery was 0.88. Holm holds that probability at 0.045, keeping FWER 0.05, and Bonferroni is strict enough to drop the FDR to 0.0035.

When looking for candidates among thousands of genes, the FDR, which accepts a few fakes and manages their proportion, is the usual choice. DESeq2's default adjustment is BH too.

## 2. How does BH change the p-values?

Let me start with the five p-values of section 11.2 of the study text, 0.001, 0.010, 0.030, 0.040 and 0.200. There are $m=5$ of them.

The idea of BH is this. Line up the p-values from smallest to largest and compare the $k$-th p-value with the reference line $(k/m)\,q$, where $q$ is the target FDR (e.g. 0.05). The line rises with the rank: rank 1 is compared with the strict threshold $q/m$, as in Bonferroni, but the threshold loosens further down. Find the last rank $k^*$ below the line, and accept ranks 1 through $k^*$ all as discoveries.

![sorted p-values with the BH line](../../figures/06_bh_line.png)

*Figure 2.* The 1100 p-values of figure 1 lined up from smallest, with only the first 65 drawn. The orange line is the BH reference line $(k/1100)\times0.05$. The last rank below the line is 50, so there are 50 discoveries, 2 of which are null genes.

The BH-adjusted p-value (padj) writes this procedure down as one number per gene. It is "the smallest $q$ at which this gene remains a discovery", so picking genes with padj of 0.05 or less gives the same list as running BH with $q=0.05$. As a formula:

$$\tilde p_{(i)}=\min\!\left\{1,\ \min_{k\ge i}\frac{m}{k}\,p_{(k)}\right\}$$

$p_{(k)}$ is the $k$-th p-value from the smallest, and $m$ is the number of p-values adjusted together. $\frac{m}{k}p_{(k)}$ is $p_{(k)}\le(k/m)\,q$ solved for $q$, so if this value is at most $q$, the $k$-th p-value lies below the reference line. $\min_{k\ge i}$ means take a smaller value from a later rank if there is one, because if a later rank is accepted, all earlier ranks are accepted too. The outer $\min\{1,\cdot\}$ caps padj at 1.

Computing it by hand, rank 3 of the study text's table is $5\times0.030/3=0.050$. Rank 4 after it is $5\times0.040/4=0.050$ and rank 5 is $5\times0.200/5=0.200$, so there is no smaller value, and the padj of rank 3 is 0.050.

```r
p <- c(0.001, 0.010, 0.030, 0.040, 0.200)
m <- length(p); o <- order(p); ps <- p[o]
raw  <- m * ps / seq_len(m)            # m p_(i) / i
step <- rev(cummin(rev(raw)))          # cumulative minimum from the back
padj_hand <- pmin(1, step)
print(data.frame(rank = 1:m, p = ps, m_p_over_i = raw, BH_hand = padj_hand, p.adjust_BH = p.adjust(p, "BH")[o],
                 bonferroni = p.adjust(p, "bonferroni")[o], holm = p.adjust(p, "holm")[o]))
p3 <- c(0.010, 0.012, 0.013)
cat("p =", p3, "| m p/i =", 3 * p3 / (1:3), "| BH =", p.adjust(p3, "BH"), "\n")
#>   rank     p m_p_over_i BH_hand p.adjust_BH bonferroni  holm
#> 1    1 0.001      0.005   0.005       0.005      0.005 0.005
#> 2    2 0.010      0.025   0.025       0.025      0.050 0.040
#> 3    3 0.030      0.050   0.050       0.050      0.150 0.090
#> 4    4 0.040      0.050   0.050       0.050      0.200 0.090
#> 5    5 0.200      0.200   0.200       0.200      1.000 0.200
#> p = 0.01 0.012 0.013 | m p/i = 0.03 0.018 0.013 | BH = 0.013 0.013 0.013
```

The hand-computed `BH_hand` equals the study text's table, 0.005, 0.025, 0.050, 0.050, 0.200, and also R's `p.adjust`.

In fact, in this table "taking from behind" changes no value at all. That ranks 3 and 4 are both 0.050 is because their $m p/i$ are equal to begin with. The example where taking from behind actually does something is the last line. With p = 0.010, 0.012, 0.013, $m p/i$ is 0.030, 0.018, 0.013, and BH brings the later 0.013 forward and makes all three 0.013. The table's Bonferroni is every p multiplied by $m$ (0.005, 0.050, 0.150, 0.200, 1.000), and Holm is in between.

<details>
<summary>Are the hand calculation and p.adjust identical down to the bit? With the p.adjust source</summary>

```r
cat("identical(hand, p.adjust):", identical(padj_hand, p.adjust(p, "BH")[o]),
    "| all.equal:", isTRUE(all.equal(padj_hand, p.adjust(p, "BH")[o])),
    "| max abs diff:", max(abs(padj_hand - p.adjust(p, "BH")[o])), "\n")
cat("did the step-up minimum change any m p/i value?", any(step != raw), "\n")
pa_src <- deparse(p.adjust)
cat(pa_src[grep("nna <- !is.na|p <- p\\[nna\\]|lp <- length", pa_src)], pa_src[grep("BH = \\{", pa_src) + 0:4],
    pa_src[grep("BY = \\{", pa_src) + 0:5], sep = "\n")
#> identical(hand, p.adjust): FALSE | all.equal: TRUE | max abs diff: 6.938894e-18
#> did the step-up minimum change any m p/i value? FALSE
#>     if (all(nna <- !is.na(p)))
#>     else p <- p[nna]
#>     lp <- length(p)
#>     }, BH = {
#>         i <- lp:1L
#>         o <- order(p, decreasing = TRUE)
#>         ro <- order(o)
#>         pmin(1, cummin(n/i * p[o]))[ro]
#>     }, BY = {
#>         i <- lp:1L
#>         o <- order(p, decreasing = TRUE)
#>         ro <- order(o)
#>         q <- sum(1/(1L:n))
#>         pmin(1, cummin(q * n/i * p[o]))[ro]
```

The hand calculation and `p.adjust` are the same within floating-point error. `identical` is FALSE (a difference of $6.9\times10^{-18}$ at rank 3), but `all.equal` is TRUE. The BH branch of `p.adjust` is `pmin(1, cummin(n/i * p[o]))[ro]`, which sweeps backwards from the largest p and takes the cumulative minimum. Only the order of computation differs; it is the same formula. The BY branch multiplies the same calculation once more by `q <- sum(1/(1L:n))` (`q * n/i * p[o]`). And `p.adjust` drops `NA` (`p <- p[nna]`) and sets `n` from the count that remains.

</details>

<details>
<summary>Code for figure 2</summary>

```r
set.seed(20260926)
z <- c(rnorm(1000, 0), rnorm(100, 3)); is_sig <- rep(c(FALSE, TRUE), c(1000, 100))
p <- 2 * pnorm(-abs(z)); o <- order(p); ps <- p[o]; sig_o <- is_sig[o]; m <- length(p); q <- 0.05; kk <- 1:65
kstar <- max(which(ps <= seq_len(m) / m * q))
png("../figures/06_bh_line.png", 7, 4.5, "in", res = 150, bg = "white")
par(mar = c(4.2, 5.6, 2.6, 1), col.axis = ink, col.lab = ink, fg = ink, las = 1)
plot(kk, ps[kk], pch = 16, cex = 0.8, col = ifelse(sig_o[kk], col_sig, col_null),
     xlab = "rank k (smallest p first)", ylab = "", bty = "l")
title(ylab = "sorted raw p-value", line = 4.2)
abline(0, q / m, col = col_line, lwd = 2)
abline(v = kstar + 0.5, lty = 3, col = ink)
text(66, 0.0012, "BH line: (k / m) x q\nm = 1100, q = 0.05", adj = c(1, 0), cex = 0.8, col = col_line)
text(kstar - 1, max(ps[kk]) * 0.75, sprintf("last rank under the line: k = %d\n-> %d discoveries (%d true nulls)",
     kstar, kstar, sum(!sig_o[1:kstar])), adj = 1, cex = 0.8, col = ink)
legend("topleft", c("true null", "true signal"), pch = 16, col = c(col_null, col_sig), bty = "n", cex = 0.85)
title("Benjamini-Hochberg: compare each sorted p-value with a rising line", adj = 0, cex.main = 1, font.main = 1, col.main = "#0b0b0b")
invisible(dev.off())
cat("figure 2 written:", file.exists("../figures/06_bh_line.png"), "| k* =", kstar,
    "| true nulls among them:", sum(!sig_o[1:kstar]), "| equals BH at 0.05:", kstar == sum(p.adjust(p, "BH") <= q), "\n")
#> figure 2 written: TRUE | k* = 50 | true nulls among them: 2 | equals BH at 0.05: TRUE
```

</details>

### Why does the same p-value get a different padj?

The BH reference line is set by the other p-values adjusted together. So even the very same p = 0.03 gets a different padj depending on which set it belongs to.

```r
cat("family A (0.001,0.01,0.03,0.04,0.2): padj of 0.03 =", p.adjust(c(0.001, 0.01, 0.03, 0.04, 0.2), "BH")[3], "\n")
cat("family B (0.03,0.5,0.6,0.7,0.9)    : padj of 0.03 =", p.adjust(c(0.03, 0.5, 0.6, 0.7, 0.9), "BH")[1], "\n")
set.seed(1)
cat("family C (0.03 + 1000 x U(0,1))    : padj of 0.03 =", p.adjust(c(0.03, runif(1000)), "BH")[1], "\n")
#> family A (0.001,0.01,0.03,0.04,0.2): padj of 0.03 = 0.05
#> family B (0.03,0.5,0.6,0.7,0.9)    : padj of 0.03 = 0.15
#> family C (0.03 + 1000 x U(0,1))    : padj of 0.03 = 0.9418149
```

The same 0.03 becomes 0.05 in set A, 0.15 in B, and as high as 0.94 in C, where it is mixed with 1000 null p-values.

The padj of gene A is the same. It is not determined by gene A's counts alone; a value comes out only once it is also settled which genes it was adjusted together with. The set of hypotheses adjusted together like this is called the test family.

BH's FDR guarantee also has a condition. When the p-values are independent or only positively dependent (a condition called PRDS), $\mathrm{FDR}\le\pi_0 q\le q$ holds. To make it hold under any dependence structure, the BY (Benjamini–Yekutieli) adjustment is used; BY multiplies the BH calculation once more by $\sum_{k=1}^m 1/k$, and the discoveries shrink accordingly (section 4).

## 3. Which genes does DESeq2's padj adjust together?

Now on to DESeq2. I made check data with three conditions, Ctrl, Starvation and Glucose, 4 samples each, and 2000 genes. Glucose here is the short name for starved cells given glucose again (Starvation+Glucose), not ordinary medium with extra glucose. I made it as a design without donor pairing, so the design is `~condition`.

I also planted the answer for each gene. The values in the table below are log2 fold changes (LFC, the log2 of the ratio of two condition means), and the effect size $|\beta|$ was drawn between 1 and 2.5. In gene5 and gene7 I also planted outliers by changing one count to an abnormally large value.

| Kind of gene | Number | Starvation vs Ctrl | Glucose vs Ctrl | Meaning |
|---|---|---|---|---|
| null | 1500 | 0 | 0 | No change at all |
| starvation_only | 150 | $\beta$ | $\beta$ | Starvation changes it and Glucose does not bring it back |
| rescued | 150 | $\beta$ | 0 | Glucose brings it back to the Ctrl level |
| partial | 100 | $\beta$ | $0.5\beta$ | Brought about halfway back |
| glucose_only | 100 | 0 | $\beta$ | Changes only under the Glucose condition |

<details>
<summary>Code that made the check data</summary>

I reused `trueIntercept` and `trueDisp` from `makeExampleDESeqDataSet` and generated counts the same way as that function ($\mu=2^{X\beta}$, `rnbinom(size=1/α)`).

```r
suppressPackageStartupMessages(library(DESeq2))
set.seed(2026); n <- 2000; m <- 12
base <- makeExampleDESeqDataSet(n = n, m = m)
cond <- factor(rep(c("Ctrl", "Starvation", "Glucose"), each = 4), levels = c("Ctrl", "Starvation", "Glucose"))
X <- model.matrix(~ cond)
cls <- rep(c("null", "starvation_only", "rescued", "partial", "glucose_only"), c(1500, 150, 150, 100, 100))
eff <- sample(c(-1, 1), n, TRUE) * runif(n, 1, 2.5)
bS <- ifelse(cls %in% c("starvation_only", "rescued", "partial"), eff, 0)
bG <- ifelse(cls == "starvation_only", bS, ifelse(cls == "partial", 0.5 * bS, ifelse(cls == "glucose_only", eff, 0)))
mu <- t(2^(X %*% t(cbind(mcols(base)$trueIntercept, bS, bG))))
cnt <- matrix(rnbinom(m * n, mu = mu, size = 1 / mcols(base)$trueDisp), ncol = m,
              dimnames = list(paste0("gene", 1:n), paste0("sample", 1:m)))
cnt[5, 1] <- 50000L; cnt[7, 6] <- 80000L
dds <- DESeqDataSetFromMatrix(cnt, DataFrame(condition = cond), ~ condition)
dds <- DESeq(dds, quiet = TRUE)
cat("resultsNames:", resultsNames(dds), "| all-zero genes:", sum(rowSums(counts(dds)) == 0), "\n")
#> converting counts to integer mode
#> resultsNames: Intercept condition_Starvation_vs_Ctrl condition_Glucose_vs_Ctrl | all-zero genes: 6
```

</details>

### One call to `results()` is one family

When `results()` is called once, DESeq2 collects only the `pvalue` column of that results table and applies BH to it. It does not combine them with the p-values of other calls (source in [Going deeper](#going-deeper)). But some genes are left out of this family, and those genes come out with `padj` of `NA`. There are two ways to be left out.

One is when `pvalue` itself is `NA`: genes with no p-value at all. Genes whose counts are 0 in every sample and genes judged by Cook's distance to have an outlier fall here. Cook's distance measures how strongly one sample shakes the coefficient estimates of that gene (the cutoff and a table of NA causes are in [note 07](07_lfc_shrinkage_and_qc.md), which covers section 13.5 of the study text). The other is when `pvalue` exists but only `padj` is `NA`: genes removed from the family by independent filtering.

The idea of independent filtering (IF) is simple. A gene with only a handful of reads can hardly give a small p-value even if there is a difference between conditions. Putting such genes into the family only makes $m$ larger, which flattens the slope $q/m$ of the BH reference line, so the other genes face a stricter threshold. So DESeq2 removes genes with a low baseMean from the family. baseMean is the mean of the normalized counts (count / size factor) over all samples, computed without looking at which sample belongs to which condition.

Where to cut is decided by DESeq2. It takes 50 quantiles of baseMean (values below which some percentage lies) as candidate cutoffs, counts the genes with `padj < alpha` for each candidate, and picks the lowest of the cutoffs whose count is close to the maximum. The default of `alpha` is 0.1.

```r
r_if   <- results(dds, contrast = c("condition", "Starvation", "Ctrl"))                # defaults: IF on, alpha = 0.1
r_if05 <- results(dds, contrast = c("condition", "Starvation", "Ctrl"), alpha = 0.05)
r_no   <- results(dds, contrast = c("condition", "Starvation", "Ctrl"), independentFiltering = FALSE)
show <- function(r) c(NA_pvalue = sum(is.na(r$pvalue)), NA_padj = sum(is.na(r$padj)), padj_0.05 = sum(r$padj < 0.05, na.rm = TRUE),
                      padj_0.1 = sum(r$padj < 0.1, na.rm = TRUE), baseMean_cutoff = round(unname(c(metadata(r)$filterThreshold, NA)[1]), 2))
print(rbind(r_if = show(r_if), r_if05 = show(r_if05), r_no = show(r_no)))
cat("stat, pvalue identical across the three calls:", identical(r_if$stat, r_no$stat),
    identical(r_if$pvalue, r_no$pvalue), identical(r_if05$pvalue, r_if$pvalue), "\n")
filt <- is.na(r_if$padj) & !is.na(r_if$pvalue); pna <- which(is.na(r_if$pvalue))
cook <- pna[which(mcols(dds)$maxCooks[pna] > qf(.99, 3, 9))]
cat("padj NA only:", sum(filt), "genes ; all below the baseMean cutoff:", all(r_if$baseMean[filt] < metadata(r_if)$filterThreshold), "\n")
cat("pvalue NA:", length(pna), "= all-zero genes", sum(r_if$baseMean[pna] == 0), "+ Cook's outliers", length(cook), "(", rownames(dds)[cook], ")\n")
#>        NA_pvalue NA_padj padj_0.05 padj_0.1 baseMean_cutoff
#> r_if           9     125       161      201            1.59
#> r_if05         9     472       163      204            5.41
#> r_no           9       9       158      198              NA
#> stat, pvalue identical across the three calls: TRUE TRUE TRUE
#> padj NA only: 116 genes ; all below the baseMean cutoff: TRUE
#> pvalue NA: 9 = all-zero genes 6 + Cook's outliers 3 ( gene5 gene7 gene1571 )
```

Across the three calls, `stat` and `pvalue` are exactly the same; only `padj` changes.

Looking into the 9 with `pvalue` `NA`, 6 are genes whose counts are 0 in every sample and the other 3 are Cook's outliers. The deliberately planted gene5 and gene7 were caught here. The 116 with a `pvalue` but `padj` `NA` all have a baseMean below the cutoff of 1.59. Turning IF off (`r_no`) fills in their `padj`, leaving only 9 `NA`.

Calling with `alpha=0.05` raises the cutoff to 5.41 and increases the `NA` to 472. Genes with padj < 0.05 number 158 without IF, 161 with the default (alpha 0.1), and 163 with `alpha=0.05`. `alpha` here is only the criterion IF uses for counting "how many padj are below this value" when it picks the cutoff, so changing it does not recompute the coefficients or the Wald p-values. `summary()` also uses this value as its default threshold. So if you intend to report at padj < 0.05, the right thing is to pass `alpha=0.05` to `results()` as well. Section 11.6 of the study text and `?results` recommend the same.

<details>
<summary>A closer look at IF: the padj of the passing genes, the cutoff grid, the actual FDP</summary>

```r
cat("chosen theta (quantile of baseMean):", round(metadata(r_if)$filterTheta, 3), "(alpha 0.1) ->",
    round(metadata(r_if05)$filterTheta, 3), "(alpha 0.05) ; max baseMean among filtered:", round(max(r_if$baseMean[filt]), 3), "\n")
cat("pvalue NA genes:", pna, "\n  baseMean:", round(r_if$baseMean[pna], 2), "\n  maxCooks:", round(mcols(dds)$maxCooks[pna], 2),
    "\n  cutoff qf(.99, p=3, m-p=9) =", round(qf(.99, 3, 9), 3), "\n")
pass <- !is.na(r_if$padj)
cat("padj(IF) == p.adjust(pvalue[pass], 'BH'):", isTRUE(all.equal(r_if$padj[pass], p.adjust(r_if$pvalue[pass], "BH"))),
    "| padj(noIF) == p.adjust(all non-NA p):",
    isTRUE(all.equal(r_no$padj[!is.na(r_no$pvalue)], p.adjust(r_no$pvalue[!is.na(r_no$pvalue)], "BH"))), "\n")
both <- pass & !is.na(r_no$padj); dd <- r_if$padj[both] - r_no$padj[both]
cat("kept genes:", sum(both), "; max(padj_IF - padj_noIF) =", signif(max(dd), 3), "; < 0:", sum(dd < 0),
    "; < -1e-12 (beyond float error):", sum(dd < -1e-12), "; mean diff =", round(mean(dd), 4), "\n")
fr <- metadata(r_if)$filterNumRej
cat("filterNumRej: theta grid length =", nrow(fr), "; numRej range =", range(fr$numRej), "\n")
tS <- bS != 0; fdp <- function(R) sum(R & !tS, na.rm = TRUE) / max(sum(R, na.rm = TRUE), 1)
cat(sprintf("true FDP (one dataset) at padj<0.1: IF=%.3f (R=%d)  noIF=%.3f (R=%d) | at padj<0.05: IF(alpha=.05)=%.3f (R=%d)  noIF=%.3f (R=%d)\n",
    fdp(r_if$padj < 0.1), sum(r_if$padj < 0.1, na.rm = TRUE), fdp(r_no$padj < 0.1), sum(r_no$padj < 0.1, na.rm = TRUE),
    fdp(r_if05$padj < 0.05), sum(r_if05$padj < 0.05, na.rm = TRUE), fdp(r_no$padj < 0.05), sum(r_no$padj < 0.05, na.rm = TRUE)))
summary(r_if05)
#> chosen theta (quantile of baseMean): 0.061 (alpha 0.1) -> 0.235 (alpha 0.05) ; max baseMean among filtered: 1.589
#> pvalue NA genes: 5 7 15 237 316 709 917 1191 1571
#>   baseMean: 3846.02 6455.97 0 0 0 0 0 0 3.64
#>   maxCooks: 33.27 33.27 NA NA NA NA NA NA 8.27
#>   cutoff qf(.99, p=3, m-p=9) = 6.992
#> padj(IF) == p.adjust(pvalue[pass], 'BH'): TRUE | padj(noIF) == p.adjust(all non-NA p): TRUE
#> kept genes: 1875 ; max(padj_IF - padj_noIF) = -4.79e-19 ; < 0: 1875 ; < -1e-12 (beyond float error): 1871 ; mean diff = -0.0144
#> filterNumRej: theta grid length = 50 ; numRej range = 30 205
#> true FDP (one dataset) at padj<0.1: IF=0.124 (R=201)  noIF=0.121 (R=198) | at padj<0.05: IF(alpha=.05)=0.067 (R=163)  noIF=0.057 (R=158)
#>
#> out of 1994 with nonzero total read count
#> adjusted p-value < 0.05
#> LFC > 0 (up)       : 102, 5.1%
#> LFC < 0 (down)     : 61, 3.1%
#> outliers [1]       : 3, 0.15%
#> low counts [2]     : 463, 23%
#> (mean count < 5)
#> [1] see 'cooksCutoff' argument of ?results
#> [2] see 'independentFiltering' argument of ?results
```

- There are 50 candidate cutoffs (`theta`), and the discoveries counted for each (`numRej`) ranged from 30 to 205. With the default `alpha=0.1` the 0.061 quantile was chosen, and with `alpha=0.05` the 0.235 quantile.
- Of the 9 with `pvalue` `NA`, the 6 all-zero genes have `baseMean=0` and `maxCooks=NA`. The other 3 have `maxCooks` of 33.27, 33.27 and 8.27, above the cutoff `qf(.99, 3, 9)` = 6.992 (3 coefficients, 12 samples).
- The padj of the passing genes is exactly `p.adjust(…, "BH")` applied to that subset only.
- In this data, for all 1875 passing genes, the padj with IF was less than or equal to the padj without IF. 1871 were smaller by more than the floating-point tolerance of $10^{-12}$, and the other 4 were equal within it. This is because the p-values of the cut low-expression genes were mostly large. It is why IF gains power, but it is not always so mathematically (counterexample below).
- Because this data has known answers, the actual FDP can be counted too. At padj < 0.1 it is 0.124 with IF (R=201) and 0.121 without (R=198); at padj < 0.05 it is 0.067 with IF (`alpha=0.05`, R=163) and 0.057 without (R=158). These are values from one data set, not the FDR itself, but they show the direction that IF increases the number of discoveries while the error rate stays at a similar level.
- That `summary(r_if05)` counts "adjusted p-value < 0.05" is because it reads `metadata(res)$alpha`.

</details>

<details>
<summary>Does IF always lower the padj of the remaining genes? (Counterexample)</summary>

If the genes being cut have small p-values, the ranks of the remaining genes drop by more than the family size does, and their padj actually grows. Below, the three genes with `filter` 1 are cut.

```r
DESeq2:::filtered_p(filter = c(1, 1, 1, 5, 5), test = c(0.001, 0.002, 0.003, 0.04, 0.5),
                    theta = c(0, 0.7), method = "BH")
#>         0%  70%
#> [1,] 0.005   NA
#> [2,] 0.005   NA
#> [3,] 0.005   NA
#> [4,] 0.050 0.08
#> [5,] 0.500 0.50
```

The fourth gene (p = 0.04) has padj 0.050 in the whole family, but after the three genes with small p are cut it becomes 0.08.

</details>

### How is IF different from "pre-selecting by p-value"?

Can genes with large p-values then be removed beforehand and BH applied? IF too ends up selecting only some genes and adjusting them, so I was confused about how the two differ. The difference is whether the criterion for removal is independent of the null genes' p-values. baseMean is computed without looking at the conditions, so for null genes it is nearly independent of the p-value. So even after cutting on baseMean, the p-values of the remaining null genes still spread evenly between 0 and 1. If instead only p < 0.5 is kept, the remaining null p-values pile up between 0 and 0.5, but BH treats them as if they were spread between 0 and 1.

I repeated a situation where every gene is null 1000 times. In that case any discovery is fake, so the FDR and the FWER are both the same, $P(V\ge1)$. BH without a filter gives 0.041, and BH after keeping only the top 50% by baseMean gives 0.046, staying near 0.05. But BH after keeping only p < 0.5 gives 0.091, almost double.

<details>
<summary>Simulation code and output comparing IF with selection by p-value</summary>

```r
set.seed(5); B <- 1000; m <- 2000
v <- matrix(NA, B, 3, dimnames = list(NULL, c("noFilter", "filter_indep_stat(baseMean top50%)", "filter_by_p(p<0.5)")))
for (b in seq_len(B)) {
  mu <- rlnorm(m, 3, 1.5); p <- runif(m)                                           # null: p ~ U(0,1), independent of baseMean
  v[b, 1] <- any(p.adjust(p, "BH") < .05)
  keep <- mu > quantile(mu, .5); v[b, 2] <- any(p.adjust(p[keep], "BH") < .05)   # select on an independent statistic (fixed cutoff)
  keep <- p < 0.5;               v[b, 3] <- any(p.adjust(p[keep], "BH") < .05)   # select on p
}
cat("global null (all 2000 genes null), B=1000, FDR = P(V>=1) at BH 0.05:\n")
for (k in 1:3) cat(sprintf("  %-38s %.3f  (MC SE %.3f)\n", colnames(v)[k], mean(v[, k]), sd(v[, k]) / sqrt(B)))
#> global null (all 2000 genes null), B=1000, FDR = P(V>=1) at BH 0.05:
#>   noFilter                               0.041  (MC SE 0.006)
#>   filter_indep_stat(baseMean top50%)     0.046  (MC SE 0.007)
#>   filter_by_p(p<0.5)                     0.091  (MC SE 0.009)
```

</details>

This difference is exactly what section 11.6 of the study text means by "adjusting after picking only genes with small raw p is not the same principle". But this simulation fixed the cutoff in advance. When, as in DESeq2, the cutoff is chosen by counting discoveries with the same p-values, the guarantee is approximate ([Going deeper](#going-deeper)).

## 4. With three comparisons, what is adjusted together?

With three conditions there are three pairs to compare. For example, if starvation halves the expression ($\Delta_S=-1$) and glucose doubles it again ($\Delta_R=+1$), the difference between Glucose and Ctrl is $\Delta_E=0$. Writing the log2 means of the conditions as $C$, $S$, $G$, the three comparisons of section 14.1 of the study text are written like this.

$$\Delta_S=S-C,\qquad \Delta_R=G-S,\qquad \Delta_E=G-C=\Delta_S+\Delta_R .$$

$\Delta_S$ is the starvation effect, $\Delta_R$ the response when glucose is given again in the starved state (R for rescue), and $\Delta_E$ the difference from Ctrl that remains even after glucose is given again.

In DESeq2, each comparison is given a contrast (a vector that sets which two conditions to compare) and `results()` is called separately. If the three conditions came from the same donors, fit with `~pair+condition`. As section 3 showed, one `results()` is one family, so calling it three times applies BH three times separately.

### What is the FDR of three separately adjusted lists combined?

I simulated 2000 times with 1000 genes that have no effect at all in any of the three comparisons. Applying BH separately in each comparison and merging the discoveries into one list, the FDR for one comparison alone was kept at 0.045, but the FDR of the merged list of all three was 0.128.

<details>
<summary>Simulation code and output for the merged list of three separately adjusted comparisons</summary>

```r
# global null: no effect in any of the three contrasts. Three contrasts sharing the group means C,S,G (m=1000 genes, 2000 runs)
set.seed(16); B <- 2000; mg <- 1000; se_g <- 0.3; se_d <- sqrt(2) * se_g
any_sep <- any_glob <- sep_k <- numeric(B)
for (b in seq_len(B)) {
  Cm <- rnorm(mg, 0, se_g); Sm <- rnorm(mg, 0, se_g); Gm <- rnorm(mg, 0, se_g)
  Pm <- 2 * pnorm(-abs(cbind(Sm - Cm, Gm - Sm, Gm - Cm) / se_d))
  sep  <- apply(Pm, 2, p.adjust, method = "BH") < .05
  any_sep[b] <- any(sep); sep_k[b] <- any(sep[, 1]); any_glob[b] <- any(p.adjust(as.vector(Pm), "BH") < .05)
}
cat(sprintf("global null, 3 contrasts: FDR of contrast S_C alone (per-contrast BH) = %.3f ; FDR of the combined list = %.3f ; global BH = %.3f\n",
            mean(sep_k), mean(any_sep), mean(any_glob)))
set.seed(17); any_ind <- numeric(B)                                     # comparison: the three comparisons independent
for (b in seq_len(B)) any_ind[b] <- any(apply(matrix(runif(mg * 3), mg), 2, p.adjust, method = "BH") < .05)
cat(sprintf("3 independent contrasts: FDR of the combined list = %.3f (theory 1 - 0.95^3 = %.3f) ; MC SE of 0.128 = %.4f\n",
            mean(any_ind), 1 - 0.95^3, sqrt(mean(any_sep) * (1 - mean(any_sep)) / B)))
#> global null, 3 contrasts: FDR of contrast S_C alone (per-contrast BH) = 0.045 ; FDR of the combined list = 0.128 ; global BH = 0.046
#> 3 independent contrasts: FDR of the combined list = 0.138 (theory 1 - 0.95^3 = 0.143) ; MC SE of 0.128 = 0.0075
```

If the three comparisons were independent, the theoretical value would be $1-0.95^3=0.143$ (the same simulation with independent p-values gives 0.138). That 0.128 is a little lower seems to be because the three comparisons share the same group means and move together (section 5). But the difference is around twice the simulation standard error (about 0.008), so it is not large.

</details>

Structurally, per-comparison BH guarantees only $E[V_k/\max(R_k,1)]$ within comparison $k$, and guarantees nothing about the merged list's $E[\sum_k V_k/\max(\sum_k R_k,1)]$. Planning the comparisons in advance does not control the error rate of the merged list either.

The solution is the global BH of section 11.4 of the study text: collect the original p-values (raw p) of the three comparisons and adjust them in one go. Calling a pair of one gene and one comparison a cell, the whole set of gene×comparison cells becomes one family. The `global BH` in the simulation above is this, and its FDR was 0.046. Running the study text's code as it is on the data of section 3:

```r
res_list <- list(
  S_C = results(dds, contrast = c("condition", "Starvation", "Ctrl"), independentFiltering = FALSE),
  G_S = results(dds, contrast = c("condition", "Glucose", "Starvation"), independentFiltering = FALSE),
  G_C = results(dds, contrast = c("condition", "Glucose", "Ctrl"),   independentFiltering = FALSE))
P <- do.call(cbind, lapply(res_list, function(x) x$pvalue)); rownames(P) <- rownames(dds)
cat("NA per contrast in P:", colSums(is.na(P)), "\n")
P[is.na(P)] <- 1                                                         # untestable cells conservatively treated as non-discoveries
Q_global <- matrix(p.adjust(as.vector(P), method = "BH"), nrow = nrow(P), dimnames = dimnames(P))
Q_BY     <- matrix(p.adjust(as.vector(P), method = "BY"), nrow = nrow(P), dimnames = dimnames(P))
Q_sep    <- sapply(res_list, function(x) x$padj)                          # per-comparison BH
Qw       <- matrix(p.adjust(as.vector(Q_sep), "BH"), nrow = nrow(P))      # wrong procedure: BH applied again to padj
cat("dim(Q_global):", dim(Q_global), "| BY factor sum(1/k) for m = 6000 cells:", round(sum(1 / (1:6000)), 2), "\n")
cat("per-contrast BH<0.05:", colSums(Q_sep < 0.05, na.rm = TRUE), " total =", sum(Q_sep < 0.05, na.rm = TRUE), "\n")
cat("global BH<0.05      :", colSums(Q_global < 0.05), " total =", sum(Q_global < 0.05), "\n")
cat("global BY<0.05      :", colSums(Q_BY < 0.05), " total =", sum(Q_BY < 0.05), "\n")
cat("BH re-applied to padj (NOT global BH) <0.05 total:", sum(Qw < 0.05, na.rm = TRUE), "\n")
#> NA per contrast in P: 9 9 9
#> dim(Q_global): 2000 3 | BY factor sum(1/k) for m = 6000 cells: 9.28
#> per-contrast BH<0.05: 158 106 96  total = 360
#> global BH<0.05      : 153 109 105  total = 367
#> global BY<0.05      : 92 62 57  total = 211
#> BH re-applied to padj (NOT global BH) <0.05 total: 162
```

Of 6000 cells (2000 genes × 3 comparisons), the number with padj < 0.05 is 360 for per-comparison BH, 367 for global BH, 211 for BY, and 162 for "BH applied again to padj".

The sum of per-comparison BH (360) and global BH (367) are different procedures, so the lists differ too. Nor is global BH always stricter, because cell by cell, cells where the global padj is larger are mixed with cells where it is smaller (collapsible block below). When a comparison with many small p-values joins the same family, the rank $k^*$ that BH accepts grows, so the effective threshold $(k^*/m)\,q$ rises, and the padj of another comparison can actually get smaller.

Applying BH again to already adjusted padj shrinks it to 162. These are values adjusted twice, and the study text does not regard this as a principled global adjustment. What is collected is the raw p. BY multiplies by a further $\sum 1/k\approx9.28$ over 6000 cells, so the discoveries dropped from 367 to 211, by about 43%. The three comparisons are correlated because of the shared comparison group (section 5), and under this correlation global BH held at 0.046 in the global-null simulation above. Of course this is a result seen under this correlation structure, not a general proof.

And what global BH controls is the FDR at the level of gene×comparison cells. A list that merges the cells into genes (genes discovered in at least one comparison: 234 for per-comparison BH, 236 for global BH) does not automatically get the same guarantee.

<details>
<summary>More checks of global BH: cell-level comparison, actual FDP, NA handling</summary>

Section 11.4 of the study text fixes the gene set with a common pre-filter and turns off per-comparison IF. Here I used all 2000 genes without a pre-filter and set `independentFiltering=FALSE` so that the three comparisons share the same gene set. With IF on, each `results()` picks its own cutoff, so the gene set differs between comparisons within one family. I did not turn off Cook's filtering. The 9 `NA` in `P` include the 3 Cook's, and `P[is.na(P)] <- 1` keeps them in the family.

```r
cat("union of genes: per-contrast =", sum(rowSums(Q_sep < 0.05, na.rm = TRUE) > 0), " global =", sum(rowSums(Q_global < 0.05) > 0), "\n")
ok <- !is.na(Q_sep)                                   # exclude pvalue-NA cells; compare with a 1e-12 tolerance
cat("Q_global > Q_sep + 1e-12 in", sum(Q_global[ok] > Q_sep[ok] + 1e-12), "cells; Q_global < Q_sep - 1e-12 in",
    sum(Q_global[ok] < Q_sep[ok] - 1e-12), "cells (of", sum(ok), ")\n")
truth <- cbind(S_C = bS != 0, G_S = (bG - bS) != 0, G_C = bG != 0)
cat(sprintf("true FDP of gene x contrast cells at 0.05 (one dataset): per-contrast BH = %.3f ; global BH = %.3f ; BH-of-padj = %.3f\n",
    sum(Q_sep < 0.05 & !truth, na.rm = TRUE) / sum(Q_sep < 0.05, na.rm = TRUE), sum(Q_global < 0.05 & !truth) / sum(Q_global < 0.05),
    sum(Qw < 0.05 & !truth, na.rm = TRUE) / sum(Qw < 0.05, na.rm = TRUE)))
p1 <- res_list$S_C$pvalue
cat("p.adjust with NA kept vs NA->1: n non-NA =", sum(!is.na(p1)), "; #padj<0.05:", sum(p.adjust(p1, "BH") < 0.05, na.rm = TRUE),
    "vs", sum(p.adjust(ifelse(is.na(p1), 1, p1), "BH") < 0.05), "\n")
#> union of genes: per-contrast = 234  global = 236
#> Q_global > Q_sep + 1e-12 in 2760 cells; Q_global < Q_sep - 1e-12 in 3073 cells (of 5973 )
#> true FDP of gene x contrast cells at 0.05 (one dataset): per-contrast BH = 0.064 ; global BH = 0.074 ; BH-of-padj = 0.006
#> p.adjust with NA kept vs NA->1: n non-NA = 1991 ; #padj<0.05: 158 vs 158
```

- Cell by cell, the global padj is larger in 2760 cells and smaller in 3073 (the 5973 cells with a `pvalue`, tolerance $10^{-12}$).
- The cell-level FDP counted against the answers is 0.064 for per-comparison BH, 0.074 for global BH and 0.006 for re-adjusting padj. These are values from one data set. As section 1 showed, the FDP of one experiment swings from 0 to 0.13, so these numbers say nothing about which procedure is better or about BH failing. The 0.006 of re-adjusting padj shows over-adjustment.
- `p.adjust` drops `NA` and sets the family size from the count that remains. `P[is.na(P)] <- 1` is the conservative choice of counting `NA` in the family size. Here there are only 9 `NA`, so the number of discoveries was the same (158 vs 158).

</details>

### What decides the family?

Which family to use is decided backwards from what you will report (section 11.3 of the study text).

| Claim to report | Test family to consider | Possible approach | This note |
|---|---|---|---|
| Report each comparison as a separate question | Genes per comparison | BH in each comparison, with the adjustment scope stated | Sections 2 and 3 |
| Report the individual differences of all comparisons as one list | Gene×comparison cells | Collect raw p, global BH/BY | Section 4 |
| Find genes with a condition effect and confirm the detailed comparisons | A hierarchical structure per gene | LRT screen + stage-wise method | Below in section 4 |
| Claim genes in which two changes hold at the same time | Genes for the joint hypothesis | Define a joint p, then adjust per gene | Section 5 |

That does not mean every analysis should use the largest family. Decide first the final claim the reader will interpret and the unit of the error rate. As chapter 11 of the study text sums it up, instead of "was it adjusted?", ask "which hypotheses were grouped together, and whose error rate was controlled?"

About the third row, just one point. Picking "genes with a condition effect" first with the LRT ([note 05](05_wald_vs_lrt.md); tests several coefficients at once) and then applying ordinary BH to the per-comparison p-values within the picked genes creates the same problem as "pre-selecting by p-value" in section 3, because the selection criterion (LRT) and the per-comparison p-values come from the same group means. In a simple simulation that draws estimates straight from a normal distribution without making counts (the normal-model simulation below), the cell-level FDR of this approach was 0.071, above 0.05, and global BH on the same p-values gave 0.043. The stage-wise method (stageR) that deals with this problem is written up in the collapsible block below.

<details>
<summary>Filtering by LRT first and then applying BH again: simulation and stageR</summary>

As the omnibus test I used a $\chi^2$ test with 2 degrees of freedom instead of the LRT. Unlike baseMean, the omnibus statistic is not independent of the per-comparison p-values of null cells.

```r
# 3-group normal model: 1600 fully null, 200 with only S shifted (Δ=1), 200 with only G shifted. cell = gene x {S-C, G-S, G-C}
set.seed(11); B <- 200; mg <- 2000; se_g <- 0.3; se_d <- sqrt(2) * se_g
shift <- rep(c("none", "S", "G"), c(1600, 200, 200))
tr <- cbind(S_C = shift == "S", G_S = shift != "none", G_C = shift == "G")   # cells with a true difference
res11 <- matrix(NA, B, 5)
for (b in seq_len(B)) {
  M <- cbind(C = rnorm(mg, 0, se_g), S = rnorm(mg, (shift == "S") * 1, se_g), G = rnorm(mg, (shift == "G") * 1, se_g))
  chi <- rowSums((M - rowMeans(M))^2) / se_g^2; p_omni <- pchisq(chi, 2, lower.tail = FALSE)   # omnibus (same role as the LRT)
  Pw <- 2 * pnorm(-abs(cbind(M[, 2] - M[, 1], M[, 3] - M[, 2], M[, 3] - M[, 1]) / se_d))
  sel <- p.adjust(p_omni, "BH") < .05                                            # screening
  Qs <- matrix(FALSE, mg, 3); Qs[sel, ] <- p.adjust(as.vector(Pw[sel, ]), "BH") < .05   # ordinary flat BH after selection
  Qg <- matrix(p.adjust(as.vector(Pw), "BH") < .05, mg)                          # global BH without selection
  fdpC <- function(Q) sum(Q & !tr) / max(sum(Q), 1)
  res11[b, ] <- c(fdpC(Qs), sum(Qs), fdpC(Qg), sum(Qg), mean(sel[shift == "none"]))
}
cat(sprintf("screen (omnibus BH<.05) then flat BH on pairwise p: cell FDR=%.3f  mean claims=%.1f\n", mean(res11[, 1]), mean(res11[, 2])))
cat(sprintf("global BH on all gene x contrast p (no screen)     : cell FDR=%.3f  mean claims=%.1f\n", mean(res11[, 3]), mean(res11[, 4])))
cat(sprintf("share of complete-null genes passing the screen    : %.4f\n", mean(res11[, 5])))
#> screen (omnibus BH<.05) then flat BH on pairwise p: cell FDR=0.071  mean claims=248.5
#> global BH on all gene x contrast p (no screen)     : cell FDR=0.043  mean claims=170.7
#> share of complete-null genes passing the screen    : 0.0032
```

Flat BH after selection gives more discoveries (an average of 248.5 vs 170.7), but its cell-level FDR is 0.071, above the nominal 0.05. The proportion of fully null genes passing the selection was 0.0032. As section 11.5 of the study text says, the reasoning "ANOVA was done first, so the post-hoc tests are safe" does not by itself control the hierarchical error rate over thousands of genes.

The stage-wise procedure of stageR (Van den Berge et al. 2017) has two stages. First, screening applies BH at α to the omnibus p-values and picks $R$ genes. Then confirmation tests the per-comparison p-values within the picked genes with an FWER method (Holm and the like) at the adjusted level $\alpha\cdot R/m$. The target is the overall gene-level FDR defined in the paper, a different quantity from the flat gene×contrast FDR. So, as section 11.5 of the study text says, the confirmation correction, the hypothesis structure and the target α used have to be stated in the report. stageR is not installed in this environment (`requireNamespace("stageR")` → FALSE), so I could not confirm stageR itself by running it.

</details>

## 5. How do you test "starvation lowered it and glucose brought it back"?

The single word rescue mixes four claims of different strength (section 14.2 of the study text). This section covers the first row, and section 6 the second and third.

| Claim | Information needed | Wording to avoid | This note |
|---|---|---|---|
| Opposite-direction response | Signs, effects and uncertainty of $\Delta_S$ and $\Delta_R$ | Jumping straight to complete rescue | Section 5 (joint p) |
| Moving closer to Ctrl | Size and uncertainty of the remaining difference | Judging from the sign flip alone | Section 6 |
| Practical equivalence | A prespecified margin $\varepsilon$ and an equivalence test | p>0.05 against Ctrl, so identical | Section 6 (`lessAbs`) |
| Mechanistic rescue | The results above plus further causal and functional validation | Settling the mechanism from transcript changes alone | Cannot be answered by statistical tests alone |

As chapter 14 of the study text sums it up, only after telling these four apart and defining the statistical endpoint of the analysis (the claim to be answered by a test) can the scope of the adjustment be decided.

### Wouldn't the intersection of two lists do?

The common method is this. Make a list of genes with padj < 0.05 that went down in Starvation vs Ctrl, and a list of genes with padj < 0.05 that went up in Glucose vs Starvation. Then call the intersection of the two lists "rescue genes". Section 14.4 of the study text says this can be used as an exploratory definition of candidates, but there is no guarantee that the FDR of the joint claim for the intersection is 5%.

I worked it out with the normal-model simulation. I split 2000 genes into four kinds, and only the 100 of A have both changes true ($\Delta_S=-1.2$, $\Delta_R=+1.2$). The 400 of B have only the starvation effect true, the 400 of C only the glucose response, and the 1100 of D have no change. The SE of one comparison is 0.424, and I repeated it 200 times.

| Method | FDR of the joint claim | Average number of claims |
|---|---|---|
| (a) Intersection of two one-sided BH lists | 0.132 | 60.5 |
| (a') Two two-sided padj < 0.05 + opposite signs | 0.077 | 42.1 |
| (b) BH on the joint p $=\max(p_S, p_R)$ | 0.021 | 15.9 |
| (c) BH on $\min(p_S, p_R)$ (wrong method) | 0.880 | |

Even though each list keeps FDR 5%, the FDR of the joint claim built from the intersection is 7.7% with the common two-sided + sign approach and 13% for the intersection of one-sided lists. The fake joint claims came mostly from B and C, where only one effect is true (summed over 200 repetitions: B 682, C 683, D 247). For example, a B gene has a real starvation effect, so it gets into the $\Delta_S$ list almost certainly. If one of the 5% fakes the $\Delta_R$ list allows lands on it, it stays in the intersection as it is.

<details>
<summary>Simulation code and output: intersection vs joint p</summary>

```r
set.seed(1); B <- 200; se_g <- 0.30                         # group-mean SE; contrast SE = sqrt(2)*se_g = 0.424
n_cls <- c(A = 100, B = 400, C = 400, D = 1100); cls2 <- rep(names(n_cls), n_cls)
dS_true <- c(A = -1.2, B = -1.2, C = 0, D = 0)[cls2]; dR_true <- c(A = 1.2, B = 0, C = 1.2, D = 0)[cls2]
joint_true <- cls2 == "A"; fdpA <- function(R) sum(R & !joint_true) / max(sum(R), 1)
out <- matrix(NA, B, 9); false_by_cls <- setNames(numeric(4), names(n_cls))
for (b in seq_len(B)) {
  Cm <- rnorm(length(cls2), 0, se_g); Sm <- rnorm(length(cls2), dS_true, se_g); Gm <- rnorm(length(cls2), dS_true + dR_true, se_g)
  dS <- Sm - Cm; dR <- Gm - Sm; se_d <- sqrt(2) * se_g                             # share the same Sm
  pS <- pnorm(dS / se_d); pR <- pnorm(dR / se_d, lower.tail = FALSE)               # correct one-sided p
  inter1 <- p.adjust(pS, "BH") < .05 & p.adjust(pR, "BH") < .05                     # (a) intersection of two lists
  p2S <- 2 * pnorm(-abs(dS / se_d)); p2R <- 2 * pnorm(-abs(dR / se_d))
  inter2 <- p.adjust(p2S, "BH") < .05 & p.adjust(p2R, "BH") < .05 & dS < 0 & dR > 0  # (a') two-sided BH + signs
  maxp <- p.adjust(pmax(pS, pR), "BH") < .05                                        # (b) BH across genes on the IUT joint p
  minp <- p.adjust(pmin(pS, pR), "BH") < .05                                        # (c) the wrong min
  Sm2 <- rnorm(length(cls2), dS_true, se_g); pR_ind <- pnorm((Gm - Sm2) / se_d, lower.tail = FALSE)   # independent copy of S
  inter_ind <- p.adjust(pS, "BH") < .05 & p.adjust(pR_ind, "BH") < .05
  out[b, ] <- c(fdpA(inter1), sum(inter1), fdpA(inter2), sum(inter2), fdpA(maxp), sum(maxp), fdpA(inter_ind),
                cor(dS[cls2 == "D"], dR[cls2 == "D"]), fdpA(minp))
  false_by_cls <- false_by_cls + table(factor(cls2[inter1 & !joint_true], levels = names(n_cls)))
}
cat(sprintf("  (a) two one-sided BH<.05 lists, intersect       : FDR=%.3f  mean #claims=%.1f\n", mean(out[, 1]), mean(out[, 2])))
cat(sprintf("  (a') two two-sided BH<.05 lists + opposite sign : FDR=%.3f  mean #claims=%.1f\n", mean(out[, 3]), mean(out[, 4])))
cat(sprintf("  (b) IUT p_joint=max(pS,pR), BH<.05 over genes   : FDR=%.3f  mean #claims=%.1f\n", mean(out[, 5]), mean(out[, 6])))
cat(sprintf("  (c) WRONG p_joint=min(pS,pR), BH<.05            : FDR=%.3f\n", mean(out[, 9])))
cat(sprintf("  (a) again but with independent S copies         : FDR=%.3f   (shared-S adds %.3f)\n", mean(out[, 7]), mean(out[, 1]) - mean(out[, 7])))
cat(sprintf("  cor(dS_hat, dR_hat) among null genes = %.3f\n", mean(out[, 8])))
cat("  false claims of (a) by class, summed over 200 reps:", paste(names(false_by_cls), false_by_cls, collapse = ", "), "\n")
#>   (a) two one-sided BH<.05 lists, intersect       : FDR=0.132  mean #claims=60.5
#>   (a') two two-sided BH<.05 lists + opposite sign : FDR=0.077  mean #claims=42.1
#>   (b) IUT p_joint=max(pS,pR), BH<.05 over genes   : FDR=0.021  mean #claims=15.9
#>   (c) WRONG p_joint=min(pS,pR), BH<.05            : FDR=0.880
#>   (a) again but with independent S copies         : FDR=0.096   (shared-S adds 0.036)
#>   cor(dS_hat, dR_hat) among null genes = -0.499
#>   false claims of (a) by class, summed over 200 reps: A 0, B 682, C 683, D 247
```

(a') is exactly the "both lists BH<0.05 + opposite signs" rule of section 14.4 of the study text, yet it is 7.7%. The negative correlation that comes from $\Delta_S$ and $\Delta_R$ sharing the same $S$ (below) adds 3.6 percentage points to the intersection FDR. Computing $\Delta_R$ from an independent copy of $S$ gives an intersection FDR of 0.096. (b) is conservative (0.021) because the joint null hypothesis is a composite null that bundles several cases.

</details>

### Why does the joint p-value use the larger one?

If the directions were set in advance, the joint hypothesis can be set up directly. The alternative to be claimed is $H_1:\ \Delta_S<0\ \text{and}\ \Delta_R>0$. Calling the one-sided p-values $p_S$ (for starvation having lowered it) and $p_R$ (for glucose having raised it), section 14.4 of the study text uses the joint p-value of the intersection–union test (IUT). The IUT is a method for testing a claim that several conditions all hold.

$$p_{joint}=\max(p_S,\ p_R).$$

Plugging in the numbers of Problem 20, $p_S=0.003$ and $p_R=0.08$, so $p_{joint}=\max(0.003, 0.08)=0.08$. Why the larger one rather than the smaller? The joint null hypothesis is "at least one of the two does not hold". For $\max(p_S,p_R)\le t$, the p-value of the side that does not hold (the side whose null is true) must also be at most $t$, and the probability of that is at most $t$. So the level of the joint test does not exceed $t$. This argument does not even need the assumption that the two comparisons are independent. All that is needed is that each p-value is a proper p-value, that is, $P(p\le t)\le t$ when its null hypothesis is true.

Conversely, choosing the smaller one (0.003) reads a gene where only one side is true as evidence of "both true". Imitating the situation where only one side is true ($\Delta_S=-1.2$ true, $\Delta_R=0$) a million times makes the difference clear. The probability that the max rule wrongly says "both true" was 0.0498, not above 0.05, while for the min rule it was 0.88.

<details>
<summary>The level of the max and min rules under a partial null (simulation)</summary>

```r
set.seed(4); N <- 1e6; se_g <- 0.3; se_d <- sqrt(2) * se_g
Cm <- rnorm(N, 0, se_g); Sm <- rnorm(N, -1.2, se_g); Gm <- rnorm(N, -1.2, se_g)   # Δ_S=-1.2 true, Δ_R=0 (partial null)
pS <- pnorm((Sm - Cm) / se_d); pR <- pnorm((Gm - Sm) / se_d, lower.tail = FALSE)
cat(sprintf("shared S      : P(max(pS,pR)<=.05)=%.4f  P(min(pS,pR)<=.05)=%.4f  cor(dS,dR)=%.3f\n",
            mean(pmax(pS, pR) <= .05), mean(pmin(pS, pR) <= .05), cor(Sm - Cm, Gm - Sm)))
Sm2 <- rnorm(N, -1.2, se_g); pR2 <- pnorm((Gm - Sm2) / se_d, lower.tail = FALSE)   # Δ_R from an independent copy of S
cat(sprintf("independent S : P(max(pS,pR)<=.05)=%.4f  P(min(pS,pR)<=.05)=%.4f  cor(dS,dR)=%.3f\n",
            mean(pmax(pS, pR2) <= .05), mean(pmin(pS, pR2) <= .05), cor(Sm - Cm, Gm - Sm2)))
pw <- pnorm(1.2 / se_d - qnorm(.95))                                               # P(pS <= .05): power of the Δ_S test
cat(sprintf("theory if independent: power(Δ_S)=%.4f ; P(max) = power*0.05 = %.4f ; P(min) = %.4f\n", pw, pw * .05, pw + .05 - pw * .05))
#> shared S      : P(max(pS,pR)<=.05)=0.0498  P(min(pS,pR)<=.05)=0.8817  cor(dS,dR)=-0.500
#> independent S : P(max(pS,pR)<=.05)=0.0438  P(min(pS,pR)<=.05)=0.8872  cor(dS,dR)=0.001
#> theory if independent: power(Δ_S)=0.8817 ; P(max) = power*0.05 = 0.0441 ; P(min) = 0.8876
```

If the two comparisons were independent, this probability for the max rule would be $\text{power}(\Delta_S)\times0.05=0.0441$ (simulation 0.0438). In the actual structure, the negative correlation of the two comparisons pushes it up to 0.0498, but 0.05 is kept. The closer the power of the $\Delta_S$ test is to 1, the closer the actual level of the IUT gets to α.

</details>

After computing $p_{joint}$ for each gene, BH is applied once across the genes ((b) in the table above). Section 14.4 of the study text adds some cautions here. The dependence conditions between genes that BH relies on and the normal approximation of the Wald test are still needed. If you proceeded by picking genes with the first test and then looking at the second result, that selection process has to be stated in the report. If you also search for the opposite direction (Starvation-up/Glucose-down), choosing the direction must be adjusted for too. In that case you can use $\min\{1,\ 2\min(p^{\downarrow\uparrow}_{joint},p^{\uparrow\downarrow}_{joint})\}$, the smaller of the two directions' joint p-values times 2, capped at 1; this is a Bonferroni adjustment over two directions, so it is conservative.

This joint test is built directly from the hypothesis, so DESeq2 does not do it automatically. All six `altHypothesis` options of `results()` are hypotheses about one coefficient (or one contrast).

### What happens when two comparisons share the Starvation mean

$\Delta_S=S-C$ and $\Delta_R=G-S$ contain the same $S$ with opposite signs. So if the Starvation mean happens to be estimated low, $\Delta_S$ moves more negative and $\Delta_R$ more positive, the direction that looks like rescue even with no effect at all. If the three group means are independent of each other, the following holds (section 14.5 of the study text).

$$\mathrm{Cov}(\hat\Delta_S,\hat\Delta_R)=\mathrm{Cov}(S-C,\ G-S)=-\mathrm{Var}(S).$$

If the variances of the three group means are all equal to $v$, then $\mathrm{Var}(\hat\Delta_S)=\mathrm{Var}(\hat\Delta_R)=2v$, so the correlation is $-v/2v=-\tfrac12$. Computing it with the actual null genes of the section 3 data:

```r
lfc <- sapply(res_list, function(x) x$log2FoldChange)
nul <- cls == "null" & !is.na(lfc[, 1])
cat("cor(Δ_S_hat, Δ_R_hat) among true-null genes =", round(cor(lfc[nul, "S_C"], lfc[nul, "G_S"]), 3),
    " (theory -0.5) ; cor(Δ_S_hat, Δ_E_hat) =", round(cor(lfc[nul, "S_C"], lfc[nul, "G_C"]), 3), "(theory +0.5)\n")
#> cor(Δ_S_hat, Δ_R_hat) among true-null genes = -0.565  (theory -0.5) ; cor(Δ_S_hat, Δ_E_hat) = 0.492 (theory +0.5)
```

Even for genes with no change at all, the correlation of $\hat\Delta_S$ and $\hat\Delta_R$ is −0.565, near the theoretical −0.5.

In a paired design where the three conditions came from the same donors (`~pair+condition`), the donor effect adds a term to the covariance formula, so the formula above is not used as it is. The principle that the shared comparison group creates a negative correlation remains, though. In a simulation with donor effects, the correlation for null genes was −0.528 (collapsible block below). So "the correlation between the starvation effect and the glucose response is negative" across genes is not, by itself, evidence of a rescue mechanism.

<details>
<summary>How does the covariance change in a paired design?</summary>

Let the donor variance be $\sigma_p^2$, the residual variance $\sigma^2$ and the number of donors $n$. The variance of one group mean is $\mathrm{Var}(\bar S)=(\sigma_p^2+\sigma^2)/n$, but the effect of the same donor cancels when taking differences, so $\mathrm{Cov}(\bar S-\bar C,\ \bar G-\bar S)=-\sigma^2/n$. So $-\mathrm{Var}(\bar S)$ is not substituted as it is. It is still negative.

```r
# paired design: 4 donors x three conditions, all genes null, a donor effect (log2 SD 0.5) per gene
set.seed(7); n2 <- 1000
pair  <- factor(rep(1:4, 3)); cond2 <- factor(rep(c("Ctrl", "Starvation", "Glucose"), each = 4), levels = c("Ctrl", "Starvation", "Glucose"))
base2 <- makeExampleDESeqDataSet(n = n2, m = 12)
mu2   <- 2^(mcols(base2)$trueIntercept + matrix(rnorm(n2 * 4, 0, 0.5), n2, 4)[, as.integer(pair)])
cnt2  <- matrix(rnbinom(n2 * 12, mu = mu2, size = 1 / mcols(base2)$trueDisp), n2)
dp <- DESeq(DESeqDataSetFromMatrix(cnt2, DataFrame(pair = pair, condition = cond2), ~ pair + condition), quiet = TRUE)
a <- results(dp, contrast = c("condition", "Starvation", "Ctrl")); b2 <- results(dp, contrast = c("condition", "Glucose", "Starvation"))
okp <- !is.na(a$pvalue)
cat("~pair+condition, all-null genes: cor(Δ_S_hat, Δ_R_hat) =", round(cor(a$log2FoldChange[okp], b2$log2FoldChange[okp]), 3), "\n")
#> converting counts to integer mode
#> ~pair+condition, all-null genes: cor(Δ_S_hat, Δ_R_hat) = -0.528
```

</details>

<details>
<summary>Arithmetic of the three comparisons: does Δ_E = Δ_S + Δ_R hold exactly? Do the SEs add too?</summary>

```r
se <- sapply(res_list, function(x) x$lfcSE); pv <- sapply(res_list, function(x) x$pvalue)
gap <- abs(lfc[, "G_C"] - (lfc[, "S_C"] + lfc[, "G_S"]))
cat("max|Δ_E - (Δ_S + Δ_R)| =", max(gap, na.rm = TRUE), "; genes with gap>1e-3:", sum(gap > 1e-3, na.rm = TRUE),
    "; median gap =", signif(median(gap, na.rm = TRUE), 3), "\n")
cat("results(G vs S) lfc vs coef difference: max|diff| =",
    signif(max(abs(lfc[, "G_S"] - (coef(dds)[, 3] - coef(dds)[, 2])), na.rm = TRUE), 3), "\n")
r_num <- results(dds, contrast = c(0, -1, 1), independentFiltering = FALSE)
cat("character c(condition,Glucose,Starvation) vs numeric c(0,-1,1): same lfc?", isTRUE(all.equal(r_num$log2FoldChange, lfc[, "G_S"])), "\n")
k <- which(gap > 1e-3)
cat("gap genes:", length(k), "; in each, some contrast has lfc==0 & pvalue==1 (contrastAllZero rule)?",
    all(apply(lfc[k, ] == 0 & pv[k, ] == 1, 1, any)), "\n")
cat("gene150: counts Ctrl/Starvation/Glucose =", tapply(counts(dds)[150, ], cond, sum), "; lfc S_C G_S G_C =", round(lfc[150, ], 3), "\n")
cat("SE not additive: max|SE_E - (SE_S+SE_R)| =", round(max(abs(se[, "G_C"] - (se[, "S_C"] + se[, "G_S"])), na.rm = TRUE), 3),
    "; median SE_E^2 / (SE_S^2 + SE_R^2) =", round(median(se[, "G_C"]^2 / (se[, "S_C"]^2 + se[, "G_S"]^2), na.rm = TRUE), 3), "\n")
set.seed(2); N <- 1e6; v3 <- c(C = 0.3, S = 0.5, G = 0.2)          # different group-mean variances (normal group-mean model)
Cm <- rnorm(N, 0, sqrt(v3["C"])); Sm <- rnorm(N, 0, sqrt(v3["S"])); Gm <- rnorm(N, 0, sqrt(v3["G"]))
cat(sprintf("Var(C,S,G)=(%.1f, %.1f, %.1f): Cov(S-C, G-S)=%.4f [-Var(S)] ; Cov(S-C, G-C)=%.4f [+Var(C)]\n",
            v3["C"], v3["S"], v3["G"], cov(Sm - Cm, Gm - Sm), cov(Sm - Cm, Gm - Cm)))
#> max|Δ_E - (Δ_S + Δ_R)| = 0.02133857 ; genes with gap>1e-3: 9 ; median gap = 1.39e-17
#> results(G vs S) lfc vs coef difference: max|diff| = 1.78e-15
#> character c(condition,Glucose,Starvation) vs numeric c(0,-1,1): same lfc? TRUE
#> gap genes: 9 ; in each, some contrast has lfc==0 & pvalue==1 (contrastAllZero rule)? TRUE
#> gene150: counts Ctrl/Starvation/Glucose = 0 0 19 ; lfc S_C G_S G_C = 0 4.69 4.71
#> SE not additive: max|SE_E - (SE_S+SE_R)| = 3.82 ; median SE_E^2 / (SE_S^2 + SE_R^2) = 0.5
#> Var(C,S,G)=(0.3, 0.5, 0.2): Cov(S-C, G-S)=-0.4994 [-Var(S)] ; Cov(S-C, G-C)=0.3004 [+Var(C)]
```

- $\hat\Delta_E=\hat\Delta_S+\hat\Delta_R$ holds exactly by coefficient arithmetic (the LFC of G vs S and the coefficient difference agree within $1.8\times10^{-15}$). The character contrast and the numeric contrast `c(0,-1,1)` give the same LFC.
- The exception is 9 genes. If the counts of both groups being compared are all 0, DESeq2 puts `lfc=0, p=1` into that comparison (the `contrastAllZero` rule; source under Going deeper). For example, gene150 is all 0 in Ctrl and Starvation, so S_C is fixed at 0, while G_S = 4.690 and G_C = 4.710.
- The SEs do not add. If the three groups have similar variances, $SE_E^2\approx\tfrac12(SE_S^2+SE_R^2)$ (median 0.5).
- The last line is a normal simulation in which the variances of the three group means were deliberately set to differ ($\mathrm{Var}(C)=0.3$, $\mathrm{Var}(S)=0.5$, $\mathrm{Var}(G)=0.2$, $N=10^6$). $\mathrm{Cov}(S-C,G-S)=-0.4994\approx-\mathrm{Var}(S)$ and $\mathrm{Cov}(S-C,G-C)=0.3004\approx\mathrm{Var}(C)$, so the covariance is set by the variance of the shared group. The $\mathrm{cor}(\hat\Delta_S,\hat\Delta_E)=0.492$ of null genes (theory +0.5) is for the same reason.

</details>

<details>
<summary>Applying the same joint test to the DESeq2 results of section 3</summary>

```r
zS <- res_list$S_C$stat; zR <- res_list$G_S$stat
pS <- pnorm(zS); pR <- pnorm(zR, lower.tail = FALSE)
joint_true  <- (bS < 0) & ((bG - bS) > 0)                          # Δ_S<0 & Δ_R>0 (Starvation-down among rescued+partial)
joint_true2 <- joint_true | ((bS > 0) & ((bG - bS) < 0))           # either of the two directions
cat("true joint genes: Starvation-down/Glucose-up =", sum(joint_true), "; either direction =", sum(joint_true2), "\n")
fdpJ <- function(R, tr = joint_true) c(FDP = round(sum(R & !tr, na.rm = TRUE) / max(sum(R, na.rm = TRUE), 1), 3),
                                        claims = sum(R, na.rm = TRUE), true = sum(R & tr, na.rm = TRUE))
inter1 <- p.adjust(pS, "BH") < .05 & p.adjust(pR, "BH") < .05
inter2 <- res_list$S_C$padj < .05 & res_list$G_S$padj < .05 & lfc[, "S_C"] < 0 & lfc[, "G_S"] > 0
p_joint <- pmax(pS, pR); iut <- p.adjust(p_joint, "BH") < .05
p_joint_rev <- pmax(pnorm(zS, lower.tail = FALSE), pnorm(zR))
both_dir <- p.adjust(pmin(1, 2 * pmin(p_joint, p_joint_rev)), "BH") < .05
cat("                                           FDP claims true\n")
cat("(a) one-sided BH lists intersect        :", fdpJ(inter1), "\n")
cat("(a') two-sided padj<.05 both + signs    :", fdpJ(inter2), "\n")
cat("(b) IUT max(pS,pR) then BH over genes   :", fdpJ(iut), "\n")
cat("(b') both directions, 2*min(pj) then BH :", fdpJ(both_dir, joint_true2), "  [truth = either direction]\n")
#> true joint genes: Starvation-down/Glucose-up = 137 ; either direction = 250
#>                                            FDP claims true
#> (a) one-sided BH lists intersect        : 0 20 20
#> (a') two-sided padj<.05 both + signs    : 0 22 22
#> (b) IUT max(pS,pR) then BH over genes   : 0 16 16
#> (b') both directions, 2*min(pj) then BH : 0.038 52 50   [truth = either direction]
```

In this data, with large effects ($|\beta|\ge1$) and only 4 replicates, none of the three procedures made fake claims. The point is not that the intersection is always wrong but that it has no guarantee. The conditions under which the guarantee breaks (many genes with only one effect true and middling power) are shown by the simulation above. The FDP and true of (a)–(b) were scored with the 137 Starvation-down/Glucose-up genes as the answers, and (b') with the 250 genes of either direction as the answers. (b') is the two-direction Bonferroni combination the study text describes.

</details>

## 6. How do you show "after glucose was given again, it matched Ctrl"?

### If the signs are opposite, did it move closer to Ctrl?

Consider a gene with $\Delta_S=-1$ and $\Delta_R=+2.5$. Starvation lowered it and glucose raised it, so the signs satisfy the rescue condition. But $\Delta_E=-1+2.5=+1.5$, so the Glucose condition is $2^{1.5}\approx2.8$ times higher than Ctrl. Going past Ctrl like this is called overshoot. The signs alone cannot say it moved closer to Ctrl.

Translating "moved closer to Ctrl" into a formula gives $|\Delta_E|<|\Delta_S|$. But counting this inequality on estimates, it holds about half the time by chance alone. In the section 3 data, the proportion satisfying this condition was 0.83 for rescued genes, which really went back, and 0.80 for partial. But it was also 0.51 for null genes, with no starvation effect at all, and 0.53 for starvation_only, which glucose did not bring back (0.14 for glucose_only).

<details>
<summary>Code and output counting |Δ_E| < |Δ_S| by kind</summary>

```r
closer <- abs(lfc[, "G_C"]) < abs(lfc[, "S_C"])                 # descriptive condition |Δ_E| < |Δ_S| of study text 14.2
cat("share with |Δ_E_hat| < |Δ_S_hat| by class:",
    paste(names(tapply(closer, cls, mean)), round(tapply(closer, cls, mean, na.rm = TRUE), 2), collapse = ", "), "\n")
#> share with |Δ_E_hat| < |Δ_S_hat| by class: glucose_only 0.14, null 0.51, partial 0.8, rescued 0.83, starvation_only 0.53
```

</details>

Comparisons between estimates carry their own uncertainty. Don't write "moved closer" without looking at both the size and the uncertainty of the remaining difference.

### Does "not significant" mean "the same"?

Suppose the p-value of Glucose vs Ctrl is 0.4. It could be because the two groups really are the same, or because there were too few samples to catch a difference. A large p-value alone cannot tell the two apart. To claim "the difference is small enough", set a margin in advance and run an equivalence test.

Section 14.3 of the study text sets the margin at 1.2-fold. In log2 that is $T=\log_2 1.2=0.263$, and anything between about 1/1.2 and 1.2-fold counts as "practically the same". This value is set for teaching and is not a criterion that fits every gene or assay. Here the hypotheses are flipped: the null hypothesis is $H_0: |\beta|\ge T$ (outside the margin) and the alternative $H_1: |\beta|<T$ (inside the margin). This is tested with two one-sided tests (TOST), and the larger of the two p-values is used.

$$p=\max\left\{P\!\left(Z<\frac{\hat\beta-T}{SE}\right),\ P\!\left(Z>\frac{\hat\beta+T}{SE}\right)\right\}$$

$\hat\beta$ is the estimated LFC (Glucose vs Ctrl), $SE$ its standard error (how much the estimate would move if the same experiment were repeated), and $Z$ a standard normal variable. The first term is the one-sided p-value for rejecting "$\beta\ge T$", and the second for rejecting "$\beta\le -T$". Both sides have to be rejected to say it is inside the margin, so the larger one is used. It is the same structure as the max rule of section 5.

Let me plug in gene A's SE of 0.289. Even if the estimated LFC were exactly 0, $p=P(Z<-0.263/0.289)\approx0.181$. With this precision, 1.2-fold equivalence cannot be shown. For $p<0.05$ at $\hat\beta=0$, we need $SE<T/z_{0.95}=0.263/1.645=0.16$, where $z_{0.95}=1.645$ is the value that cuts off the upper 5% of the standard normal. In DESeq2, `altHypothesis="lessAbs"` is this test.

```r
Tt <- log2(1.2)
r_eq <- results(dds, contrast = c("condition", "Glucose", "Ctrl"), lfcThreshold = Tt, altHypothesis = "lessAbs", alpha = 0.05)
gc <- res_list$G_C
cat("equivalent within 1.2x (padj<0.05):", sum(r_eq$padj < 0.05, na.rm = TRUE), "of", sum(!is.na(r_eq$padj)), "genes\n")
cat("two-sided p(G-C) > 0.3:", sum(gc$pvalue > 0.3, na.rm = TRUE), "genes ; of them equivalent at 1.2x:",
    sum(gc$pvalue > 0.3 & r_eq$padj < 0.05, na.rm = TRUE), "\n")
cat("lfcSE(G-C) over genes, 5/25/50/75/95%:", round(quantile(gc$lfcSE, c(.05, .25, .5, .75, .95), na.rm = TRUE), 2),
    "| needed even at b=0: SE < T/z_.95 =", round(Tt / qnorm(.95), 3), "-> genes:", sum(gc$lfcSE < Tt / qnorm(.95), na.rm = TRUE), "\n")
tost <- function(b, s) max(pnorm((b - Tt) / s), pnorm((b + Tt) / s, lower.tail = FALSE))
cat(sprintf("gene A SE, b = 0          : TOST p = %.3f\n", tost(0, 0.289057)))
cat(sprintf("Ex19: b = 0.30, SE = 0.36 : two-sided p = %.3f ; lessAbs(T=log2 1.2=%.3f) p = %.3f\n", 2 * pnorm(-0.30 / 0.36), Tt, tost(0.30, 0.36)))
#> equivalent within 1.2x (padj<0.05): 0 of 1991 genes
#> two-sided p(G-C) > 0.3: 1292 genes ; of them equivalent at 1.2x: 0
#> lfcSE(G-C) over genes, 5/25/50/75/95%: 0.36 0.49 0.69 1.03 2.17 | needed even at b=0: SE < T/z_.95 = 0.16 -> genes: 0
#> gene A SE, b = 0          : TOST p = 0.181
#> Ex19: b = 0.30, SE = 0.36 : two-sided p = 0.405 ; lessAbs(T=log2 1.2=0.263) p = 0.541
```

There are 1292 genes whose two-sided Glucose vs Ctrl p is above 0.3, yet 0 genes can be called equivalent within 1.2-fold.

The reason is the SE. In this data even the 5% quantile of the Glucose vs Ctrl `lfcSE` is 0.36 (median 0.69), so no gene is below 0.16. An equivalence claim is decided by the number of samples and the margin, not by "the p was large". The last line is the situation of Problem 19. With an estimated LFC of 0.30 and an SE of 0.36, the two-sided p is 0.405 but the equivalence p is 0.541.

Widening the margin to 2-fold ($T=1$) makes 89 come out equivalent. But 77 of them are null genes that starvation never changed in the first place (11 rescued, 1 starvation_only). So the equivalence Glucose ≈ Ctrl alone cannot pick out rescue genes. Rescue is the joint claim "$\Delta_S\ne0$ (or $<0$) and $|\Delta_E|<\varepsilon$", handled as in section 5 with the maximum of two p-values. Equivalence p-values also need an adjustment that fits the family of genes being reported.

<details>
<summary>More checks of lessAbs: comparison with a hand-computed TOST, and widening the margin</summary>

```r
cat("lessAbs: description(pvalue) =", mcols(r_eq)$description[colnames(r_eq) == "pvalue"], "; lfcThreshold in metadata =", metadata(r_eq)$lfcThreshold, "\n")
g <- which(!is.na(r_eq$pvalue))[1:3]
cat("manual TOST p = max(P(Z<(b-T)/SE), P(Z>(b+T)/SE)):",
    signif(pmax(pnorm((r_eq$log2FoldChange[g] - Tt) / r_eq$lfcSE[g]), pnorm((r_eq$log2FoldChange[g] + Tt) / r_eq$lfcSE[g], lower.tail = FALSE)), 4),
    " vs DESeq2:", signif(r_eq$pvalue[g], 4), "\n")
i <- which(gc$pvalue > 0.3 & gc$pvalue < 0.5 & cls == "rescued")[1:4]
cat("examples (rescued genes, p(G-C) in 0.3-0.5): lfc=", round(gc$log2FoldChange[i], 2), " SE=", round(gc$lfcSE[i], 2),
    " p_two=", round(gc$pvalue[i], 2), " p_lessAbs=", round(r_eq$pvalue[i], 2), "\n")
for (fc in c(1.5, 2)) {
  r <- results(dds, contrast = c("condition", "Glucose", "Ctrl"), lfcThreshold = log2(fc), altHypothesis = "lessAbs", alpha = 0.05)
  e <- which(r$padj < .05)
  cat(sprintf("lessAbs tolerance %.1fx (T=%.3f): padj<0.05 -> %d genes; by class: %s ; true |Δ_E|=|β_G| < T: %d\n", fc, log2(fc), length(e),
              paste(names(table(cls[e])), table(cls[e]), collapse = ", "), sum(abs(bG[e]) < log2(fc))))
}
cat("true |β_G| of the starvation_only call:", round(abs(bG[e][cls[e] == "starvation_only"]), 3), "\n")
#> lessAbs: description(pvalue) = Wald test p-value: condition Glucose vs Ctrl ; lfcThreshold in metadata = 0.2630344
#> manual TOST p = max(P(Z<(b-T)/SE), P(Z>(b+T)/SE)): 0.4965 0.4211 0.5701  vs DESeq2: 0.4965 0.4211 0.5701
#> examples (rescued genes, p(G-C) in 0.3-0.5): lfc= -0.63 -1.59 0.36 -0.86  SE= 0.65 2.13 0.45 0.97  p_two= 0.33 0.45 0.42 0.37  p_lessAbs= 0.71 0.73 0.59 0.73
#> lessAbs tolerance 1.5x (T=0.585): padj<0.05 -> 0 genes; by class:  ; true |Δ_E|=|β_G| < T: 0
#> lessAbs tolerance 2.0x (T=1.000): padj<0.05 -> 89 genes; by class: null 77, rescued 11, starvation_only 1 ; true |Δ_E|=|β_G| < T: 88
#> true |β_G| of the starvation_only call: 1.216
```

- DESeq2's lessAbs p-values equal the hand-computed TOST to the fourth decimal place. The description of `pvalue` in the results table still says "Wald test p-value".
- The four rescued genes with a two-sided p of 0.3–0.5 have SEs of 0.45–2.13, so their lessAbs p is 0.59–0.73.
- With a 1.5-fold margin there are 0, and with 2-fold 89. 88 of the 89 really are equivalent, with $|\Delta_E|=|\beta_G|<1$, and the only false one is 1 starvation_only gene with $|\beta_G|=1.216$.

</details>

### Can it be summarized with a single rescue ratio?

When you want to summarize how far it went back with a single number, section 14.6 of the study text gives this ratio.

$$R_{\log}=-\frac{\Delta_R}{\Delta_S}$$

On the log2 scale, 0 means no reversal, 1 means a return to the log mean of Ctrl, and above 1 means overshoot. For example, $\Delta_S=-2$ (a 4-fold decrease) and $\Delta_R=+1$ give $R_{\log}=0.5$. But on the count scale it is different. With Ctrl at 1, Starvation is $2^{-2}=0.25$ and Glucose $2^{-1}=0.5$, so the recovery fraction is $(0.5-0.25)/(1-0.25)=0.333$. So state which scale the ratio is on.

Computing it by kind of gene in the section 3 data, rescued genes with large effects have a median of 0.96, matching the true value of 1 well. But null genes, whose starvation effect is 0, also give a plausible-looking "half rescue" with a median of 0.48 (collapsible block below).

For null genes, $R_{\log}$ is the ratio of two estimates with mean 0. The denominator is near 0, so the value swings widely: the middle 50% of the values spread between −0.38 and 1.27 (the 25%–75% quantile range, IQR). Moreover, the negative correlation of section 5 tilts it toward positive values. The median of the ratio of two normal estimates with mean 0 is "(correlation of the two) × (ratio of their standard deviations)", here $\mathrm{corr}(-\hat\Delta_R,\hat\Delta_S)\cdot SE_R/SE_S=0.5\times1=0.5$. Even genes with true $R=1$ are unstable if $\Delta_S$ is small. In a normal simulation, $\Delta_S=-0.2$ gives an IQR of [−0.05, 1.43], while $\Delta_S=-2$ stabilizes it at [0.90, 1.11] (collapsible block below).

So section 14.6 of the study text says to report, instead of a single ratio, a minimum criterion for effect size, the uncertainty of the estimates and consistency across donors together.

<details>
<summary>Rescue ratio: values by kind of gene, a normal simulation, ratios of shrunken LFCs</summary>

```r
R_log <- -lfc[, "G_S"] / lfc[, "S_C"]
for (k in c("rescued", "partial", "starvation_only", "null")) {
  x <- R_log[cls == k & !is.na(R_log)]
  cat(sprintf("R_log %-12s true=%-4s median=%.2f  IQR=[%.2f, %.2f]  P(R<0)=%.2f\n", k,
              c(rescued = "1", partial = "0.5", starvation_only = "0", null = "0/0")[k], median(x), quantile(x, .25), quantile(x, .75), mean(x < 0)))
}
#> R_log rescued      true=1    median=0.96  IQR=[0.70, 1.24]  P(R<0)=0.07
#> R_log partial      true=0.5  median=0.47  IQR=[0.18, 0.80]  P(R<0)=0.15
#> R_log starvation_only true=0    median=0.03  IQR=[-0.31, 0.34]  P(R<0)=0.46
#> R_log null         true=0/0  median=0.48  IQR=[-0.38, 1.27]  P(R<0)=0.34
```

```r
set.seed(3); N <- 1e6; sg <- 0.3 / sqrt(2)                   # group-mean SE; contrast SE = 0.3, corr(Δ_S_hat, Δ_R_hat) = -0.5
rr <- function(dS, dR) { Cm <- rnorm(N, 0, sg); Sm <- rnorm(N, dS, sg); Gm <- rnorm(N, dS + dR, sg); -(Gm - Sm) / (Sm - Cm) }
for (d in c(-0.2, -2)) { R <- rr(d, -d)
  cat(sprintf("true R=1, Δ_S=%.1f, contrast SE 0.3 (shared S): median=%.2f  IQR=[%.2f, %.2f]  P(R<0)=%.2f\n",
              d, median(R), quantile(R, .25), quantile(R, .75), mean(R < 0))) }
cat(sprintf("null (Δ_S=Δ_R=0): median R=%.3f  [theory: corr(-Δ_R, Δ_S)*SE_R/SE_S = 0.5]\n", median(rr(0, 0))))
cat("count-scale recovery for Δ_S=-2, Δ_R=+1:", round((2^-1 - 2^-2) / (1 - 2^-2), 3), "vs R_log =", -1 / -2, "\n")
#> true R=1, Δ_S=-0.2, contrast SE 0.3 (shared S): median=0.67  IQR=[-0.05, 1.43]  P(R<0)=0.26
#> true R=1, Δ_S=-2.0, contrast SE 0.3 (shared S): median=1.00  IQR=[0.90, 1.11]  P(R<0)=0.00
#> null (Δ_S=Δ_R=0): median R=0.500  [theory: corr(-Δ_R, Δ_S)*SE_R/SE_S = 0.5]
#> count-scale recovery for Δ_S=-2, Δ_R=+1: 0.333 vs R_log = 0.5
```

Even with true $R=1$, $\Delta_S=-0.2$ gives a median of 0.67 and $P(R<0)=0.26$. For nulls, the median of the normal simulation is 0.500, matching the 0.48 of the DESeq2 null genes.

A ratio made by dividing two shrunken LFCs is not a proper joint posterior ratio either. The ratio of two posterior summaries generally differs from the posterior summary of the ratio ($E[a]/E[b]\ne E[a/b]$), and two values shrunk separately do not reflect the correlation from the shared $S$. apeglm and ashr are not in this environment, so I did not confirm this part by running it.

</details>

## Summary

- The p-values of null genes spread evenly between 0 and 1. So a list picked at p < 0.05 alone has fakes mixed in, in proportion to the number of null genes (FDR 0.37 in the simulation).
- The FDR is the average proportion of fakes in the list of discoveries over repetitions. BH controls it, but the actual proportion in one experiment swings. If there must be no false discovery at all, use the FWER (Holm).
- padj depends on the p-values adjusted together. DESeq2's default family is the genes that pass Cook's and independent filtering within one call to `results()`. If even `pvalue` is `NA` it is all-zero or a Cook's outlier; if only `padj` is `NA` it is independent filtering. If the final threshold is padj < 0.05, pass `alpha=0.05`.
- To report several comparisons as one list, collect the raw p-values and adjust them in one go (global BH). Don't adjust padj again.
- "Starvation lowered it and glucose raised it" is tested with $p_{joint}=\max(p_S,p_R)$. The intersection of two significant lists has no FDR guarantee, and $\Delta_S$ and $\Delta_R$ share the same Starvation mean, so they are negatively correlated.
- "Not significant" is not "the same". Equivalence is tested with a margin set in advance and `altHypothesis="lessAbs"`. An opposite-direction response, moving closer to Ctrl, practical equivalence and mechanistic rescue are different claims.

The next note, [07](07_lfc_shrinkage_and_qc.md), covers how LFC shrinkage for ranking and effect-size reporting (chapter 13) does not change the p-values and padj in the default call (`lfcThreshold=0`, `svalue=FALSE`), and the exceptions. [Note 08](08_one_gene_end_to_end.md) follows one gene from counts to padj from beginning to end.

## Exercises

The problems of appendix A of the study text that correspond to the content of this note (chapters 11 and 14).

**Problem 16.** Each of three comparisons was BH-adjusted. Is the FDR of the merged list of individual gene×contrast results automatically 5% too? For a global adjustment, are raw p or padj collected?

<details>
<summary>Solution</summary>

It is not guaranteed automatically. Per-comparison BH controls only $E[V_k/\max(R_k,1)]$ within each comparison, not the merged list's $E[\sum V_k/\max(\sum R_k,1)]$. In the global-null simulation of section 4, one comparison was 0.045 but the merged list was 0.128, and global BH on the collected raw p was 0.046. In the section 3 data too, the two procedures give different lists (per-comparison sum 360 cells, global 367 cells). The cell FDPs there, 0.064 and 0.074, are values from one data set and no basis for ranking them. Applying BH again to padj gives 162 cells (FDP 0.006), an over-adjustment with a different meaning. What is collected is the raw p, the same as the study text's answer.

</details>

**Problem 17.** You found a gene that has a raw p but whose padj alone is NA. Which procedure should be checked first?

<details>
<summary>Solution</summary>

Independent filtering. In section 3, all 116 such genes had `baseMean < metadata(res)$filterThreshold` = 1.59, and calling again with `independentFiltering=FALSE` filled in their padj (125 `NA` → 9). Conversely, the 9 with `pvalue` itself `NA` were 6 all-zero genes and 3 Cook's outliers (`maxCooks` 33.27, 33.27, 8.27 > `qf(.99, 3, 9)` = 6.992).

The study text's answer also lists "fitting problems" as a cause when even p is NA, which splits into two cases. Non-convergence of IRLS is only flagged with `mcols(dds)$betaConv=FALSE` and a message; it does not turn p into `NA`. The only line in `results()` that puts `NA` into `pvalue` is the Cook's rule (source check under Going deeper). On the other hand, rows whose coefficients cannot be estimated because of weights are flagged `weightsFail`, handled like all-zero, and get `NA` even for p. The results of deliberately creating and running both cases are in "When even p is NA" under Going deeper in [note 07](07_lfc_shrinkage_and_qc.md). So it only partly agrees with the study text's answer.

</details>

**Problem 19.** Starvation−Ctrl<0 and Glucose−Starvation>0 are both significant, and Glucose−Ctrl has p=0.4. Is this enough to claim complete normalization?

<details>
<summary>Solution</summary>

It is not enough. In numbers, with $\hat\Delta_E=0.30$ and SE 0.36, the two-sided p is 0.405. But the lessAbs p with a 1.2-fold margin ($T=0.263$) is $\max\{P(Z<(0.30-0.263)/0.36),\ P(Z>(0.30+0.263)/0.36)\}=0.541$, so it fails the equivalence test (the `Ex19` line in section 6). The sign flip holds under overshoot too, and in the section 3 data, 0 of the 1292 genes with a two-sided p > 0.3 passed equivalence. "Complete normalization" is a separate claim that needs a prespecified $\varepsilon$ and an equivalence test, a joint claim with the starvation effect, and a gene-family adjustment of those p-values. The same as the study text's answer.

</details>

**Problem 20.** The one-sided p-values of Starvation-down and Glucose-up are 0.003 and 0.08. What is the intersection–union joint p-value that both conditions hold at the same time? Why must the smaller one not be chosen?

<details>
<summary>Solution</summary>

$p_{joint}=\max(0.003, 0.08)=0.08$. The joint null hypothesis is "at least one of the two is null", so for $\max\le t$, the p of the side whose null is true must also be at most $t$, and the probability of that is at most $t$. With $\min$, for genes where only one side is true (genes with only a starvation effect), $P(\min\le0.05)=0.882$, so "at the same time" would be wrongly claimed almost always. Under the same conditions $\max$ is 0.0498, not above 0.05 (section 5, the shared-$S$ structure). If both directions are explored, the choice of direction is also adjusted for with $\min\{1,\ 2\min(p^{\downarrow\uparrow}_{joint}, p^{\uparrow\downarrow}_{joint})\}$. The same as the study text's answer.

</details>

## Going deeper

<details>
<summary>What I checked directly in the DESeq2 source</summary>

First, the multiple-testing defaults of `results()`.

```r
suppressPackageStartupMessages(library(DESeq2))
args(DESeq2::results)
#> function (object, contrast, name, lfcThreshold = 0, altHypothesis = c("greaterAbs",
#>     "greaterAbsUPSHOT", "lessAbs", "greater", "less", "greaterAbs2014"),
#>     listValues = c(1, -1), cooksCutoff, independentFiltering = TRUE,
#>     alpha = 0.1, filter, theta, pAdjustMethod = "BH", filterFun,
#>     format = c("DataFrame", "GRanges", "GRangesList"), saveCols = NULL,
#>     test = c("Wald", "LRT"), addMLE = FALSE, tidy = FALSE, parallel = FALSE,
#>     BPPARAM = bpparam(), minmu = 0.5)
#> NULL
```

The defaults are `independentFiltering=TRUE`, `alpha=0.1` and `pAdjustMethod="BH"`. All six options of `altHypothesis` are hypotheses about one coefficient (or one contrast); there is no joint hypothesis over two contrasts. `alpha` is passed on to the p-value adjustment step after its range is checked, and the adjustment is done by `pvalueAdjustment()` (or a user-supplied `filterFun`).

```r
d <- deparse(DESeq2::results)
cat(d[c(17:18, 270:276)], sep = "\n")
#>     stopifnot(length(alpha) == 1)
#>     stopifnot(alpha > 0 & alpha < 1)
#>     if (missing(filterFun)) {
#>         deseqRes <- pvalueAdjustment(deseqRes, independentFiltering,
#>             filter, theta, alpha, pAdjustMethod)
#>     }
#>     else {
#>         deseqRes <- filterFun(deseqRes, filter, alpha, pAdjustMethod)
#>     }
```

Next, the needed lines of `DESeq2:::pvalueAdjustment`, picked out by line number (fallback rules in between, such as `0.9 * maxFit`, are left out).

```r
pa <- deparse(DESeq2:::pvalueAdjustment)
cat(pa[c(3:22, 30:35, 46:48, 55:60)], sep = "\n")
#>     if (independentFiltering) {
#>         if (missing(filter)) {
#>             filter <- res$baseMean
#>         }
#>         if (missing(theta)) {
#>             lowerQuantile <- mean(filter == 0)
#>             if (lowerQuantile < 0.95)
#>                 upperQuantile <- 0.95
#>             else upperQuantile <- 1
#>             theta <- seq(lowerQuantile, upperQuantile, length = 50)
#>         }
#>         stopifnot(length(theta) > 1)
#>         stopifnot(length(filter) == nrow(res))
#>         filtPadj <- filtered_p(filter = filter, test = res$pvalue,
#>             theta = theta, method = pAdjustMethod)
#>         numRej <- MatrixGenerics::colSums(filtPadj < alpha, na.rm = TRUE)
#>         lo.fit <- lowess(numRej ~ theta, f = 1/5)
#>         if (max(numRej) <= 10) {
#>             j <- 1
#>         }
#>             maxFit <- max(lo.fit$y)
#>             rmse <- sqrt(mean(residual^2))
#>             thresh <- maxFit - rmse
#>             j <- if (any(numRej > thresh)) {
#>                 which(numRej > thresh)[1]
#>             }
#>         padj <- filtPadj[, j, drop = TRUE]
#>         cutoffs <- quantile(filter, theta)
#>         filterThreshold <- cutoffs[j]
#>         metadata(res)[["alpha"]] <- alpha
#>     }
#>     else {
#>         padj <- p.adjust(res$pvalue, method = pAdjustMethod)
#>     }
#>     res$padj <- padj
```

Four facts can be read from this.

1. The only thing adjusted is the `pvalue` column of this `res` object. There is no code that combines it with other `results()` calls.
2. The filter statistic is `baseMean`, and the `theta` grid has 50 points from the proportion of genes with baseMean 0 up to 0.95.
3. Within the padj calculation, `alpha` is used only as the criterion for picking the cutoff by counting "how many are significant" at each `theta` of the grid, `numRej <- colSums(filtPadj < alpha)`. It does not enter the p-value calculation. Otherwise its range is checked with `stopifnot`, it is stored as `metadata(res)$alpha` and becomes the default threshold of `summary()` (the "adjusted p-value < 0.05" of `summary(r_if05)` in section 3), and it is passed on to `filterFun` if one is given.
4. With `independentFiltering=FALSE` it is just `p.adjust(..., "BH")`.

`filtered_p` lives inside DESeq2, not in genefilter. IF works even though genefilter is not installed in this environment, and `?results` also says "filtered_p R code is now copied into DESeq2 package".

```r
cat("filtered_p in DESeq2 namespace:", exists("filtered_p", envir = asNamespace("DESeq2")),
    "| genefilter installed:", requireNamespace("genefilter", quietly = TRUE), "\n")
cat(deparse(DESeq2:::filtered_p)[6:17], sep = "\n")
#> filtered_p in DESeq2 namespace: TRUE | genefilter installed: FALSE
#>     cutoffs <- quantile(U1, theta)
#>     result <- matrix(NA_real_, length(U1), length(cutoffs))
#>     colnames(result) <- names(cutoffs)
#>     for (i in 1:length(cutoffs)) {
#>         use <- U1 >= cutoffs[i]
#>         if (any(use)) {
#>             if (is.function(test))
#>                 U2 <- test(data[use, ])
#>             else U2 <- test[use]
#>             result[use, i] <- p.adjust(U2, method)
#>         }
#>     }
```

For each `theta`, only genes with "baseMean ≥ quantile" are BH-adjusted, and the padj of cut genes stays `NA_real_`. The description of `alpha` in `?results` says the same: *"the significance cutoff used for optimizing the independent filtering (by default 0.1). If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value."*

I also looked for the paths by which `pvalue` becomes `NA` (Problem 17).

```r
d <- deparse(DESeq2::results); w <- deparse(DESeq2::nbinomWaldTest)
cat(grep("pvalue.*<- NA", d, value = TRUE), sep = "\n")                      # the line in results() that puts NA into pvalue
cat(grep("betaConv <- fit|did not converge|WithNARows\\(resultsList|df <- ifelse", w, value = TRUE), sep = "\n")
#>         res$pvalue[which(cooksOutlier)] <- NA
#>         df <- ifelse(df > 0, df, NA)
#>     betaConv <- fit$betaConv
#>             message(paste(sum(!betaConv), "rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest"))
#>     WaldResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
```

The only line in `results()` that puts `NA` into `pvalue` is `res$pvalue[which(cooksOutlier)] <- NA`. All-zero genes become `NA` rows when `nbinomWaldTest` builds the results, via `buildDataFrameWithNARows(..., mcols(object)$allZero)`. Non-convergence of IRLS is only recorded in `betaConv` with a message. `nbinomWaldTest` also has the line `df <- ifelse(df > 0, df, NA)`, but that handles the degrees of freedom when `useT=TRUE` uses the t distribution, and has nothing to do with the default path (`useT=FALSE`, the normal distribution).

`altHypothesis="lessAbs"` is the maximum of two one-sided p-values, that is, TOST. The same function also stops if `lfcThreshold` is 0 or the object was fitted with `betaPrior=TRUE`.

```r
d <- deparse(DESeq2::results)
cat(d[c(34:35, 141:142, 194:200)], sep = "\n")
#>     if (lfcThreshold == 0 & altHypothesis == "lessAbs") {
#>         stop("when testing altHypothesis='lessAbs', set the argument lfcThreshold to a positive value")
#>         if (altHypothesis == "lessAbs" & attr(object, "betaPrior")) {
#>             stop("testing altHypothesis='lessAbs' requires setting the DESeq() argument betaPrior=FALSE")
#>         else if (altHypothesis == "lessAbs") {
#>             newStatAbove <- pmax((T - LFC)/SE, 0)
#>             pvalueAbove <- pfunc((T - LFC)/SE)
#>             newStatBelow <- pmax((LFC + T)/SE, 0)
#>             pvalueBelow <- pfunc((LFC + T)/SE)
#>             newStat <- pmin(newStatAbove, newStatBelow)
#>             newPvalue <- pmax(pvalueAbove, pvalueBelow)
```

The code of section 14.3 of the study text passes `lfcThreshold=log2(1.2)`, so it satisfies the first condition. If `betaPrior` is omitted, `DESeq()` sets `if (missing(betaPrior)) betaPrior <- FALSE` in its body, so the second condition is also met on the default path ([note 04](04_glm_condition_batch.md)).

A contrast between two non-reference levels, such as Glucose vs Starvation, is not a refit. `DESeq2:::getContrast` passes the stored coefficients and calls `fitBeta(..., maxitSEXP = 0, ...)` to compute only $c^\top\hat b$ and $\sqrt{c^\top\Sigma c}$ (the covariance calculation of `getContrast` is in Going deeper C of [note 04](04_glm_condition_batch.md), and contrast extraction from an LRT object in [note 05](05_wald_vs_lrt.md)). But `cleanContrast` has a rule that sets `log2FoldChange=0, stat=0, pvalue=1` if the counts of both groups being compared are all 0, so for such genes $\hat\Delta_E=\hat\Delta_S+\hat\Delta_R$ breaks (the "Arithmetic of the three comparisons" block in section 5).

```r
g <- deparse(DESeq2:::getContrast)
cat(g[grep("fitBeta\\(", g):grep("maxitSEXP", g)], sep = "\n")
cat(deparse(DESeq2:::cleanContrast)[187:191], sep = "\n")
#>     betaRes <- fitBeta(ySEXP = countsMatrix, xSEXP = modelMatrix,
#>         nfSEXP = normalizationFactors, alpha_hatSEXP = alpha_hat,
#>         contrastSEXP = contrast, beta_matSEXP = beta_mat, lambdaSEXP = lambda,
#>         weightsSEXP = weights, useWeightsSEXP = useWeights, tolSEXP = 1e-08,
#>         maxitSEXP = 0, useQRSEXP = FALSE, minmuSEXP = minmu)
#>     contrastAllZero <- contrastAllZero & !mcols(object)$allZero
#>     if (sum(contrastAllZero) > 0) {
#>         res$log2FoldChange[contrastAllZero] <- 0
#>         res$stat[contrastAllZero] <- 0
#>         res$pvalue[contrastAllZero] <- 1
```

</details>

<details>
<summary>When does independent filtering keep the FDR?</summary>

If the filter statistic $U_1$ (baseMean) and the test statistic $U_2$ (the Wald p) are independent under the null hypothesis, then keeping only genes with $U_1\ge$ a cutoff fixed in advance and applying BH leaves the remaining null p-values still U(0,1), so the FDR is kept (Bourgon et al. 2010). DESeq2 picks the cutoff on a grid of baseMean quantiles (the `theta` argument of `results()`) as the lowest quantile whose "number with padj < alpha" is close to the maximum. This `theta` is a different symbol from the $\theta=\log\alpha$ of appendix C. Because this choice is made by looking at the same p-values, the guarantee for a fixed cutoff applies, strictly speaking, only approximately (picking from 50 grid points can add a little optimism). The simulation in section 3 checked only a fixed 50% cutoff.

</details>

<details>
<summary>Common misconceptions</summary>

| Misconception | In fact |
|---|---|
| Each gene with padj < 0.05 has a 5% chance of being wrong | The FDR is the expected proportion over the whole list. The FDP of one data set spread over 0–13% (section 1) |
| Running three `results()` makes DESeq2 adjust the whole thing on its own | `pvalueAdjustment` looks only at the `pvalue` column of one `res`. To combine them, collect the raw p |
| The comparisons were planned, so BH on each makes the whole 5% too | Under the global null, per comparison it was 0.045, but the merged list was 0.128 (section 4) |
| Global BH is always more conservative | Cell by cell, 2760 got larger and 3073 smaller. The BH reference line is set by the p distribution of the whole family |
| If padj is NA, the test failed | If `pvalue` exists and only `padj` is NA, it is independent filtering (low baseMean). If `pvalue` is NA too, it is Cook's or all-zero |
| Changing `alpha` changes the results (β, p) | `stat` and `pvalue` are identical. In the padj calculation `alpha` is the criterion for counting significant genes when picking the filter cutoff; otherwise it is stored as `metadata(res)$alpha` and used as the default threshold of `summary()` |
| When IF shrinks the family, the padj of the remaining genes always gets smaller | If the cut genes have small p, it actually grows (the counterexample in section 3, 0.050 → 0.08). In the section 3 data it got smaller because the p-values of the cut low-expression genes were large |
| Independent filtering is also "adjusting only a selected subset", so it is the same as selection by p | If the filter statistic is independent of p under the null hypothesis, the FDR is kept (0.046 vs 0.091) |
| If two comparisons each have BH < 0.05 with opposite signs, the rescue list has FDR 5% too | The FDR of the joint claim rose to 7.7–13%. Adjust the joint p $=\max(p_1,p_2)$ across genes |
| The joint p should be the stronger evidence, $\min(p_1,p_2)$ | Under a partial null, $P(\min\le0.05)=0.88$. For "at the same time", the weaker side has to pass |
| Glucose−Ctrl has p = 0.4, so it returned to Ctrl | Of the 1292 with a two-sided p > 0.3, 0 passed 1.2-fold equivalence. Equivalence is a separate test |
| Once Glucose ≈ Ctrl equivalence is confirmed, it is rescue | Of the 89 equivalent at a 2-fold margin, 77 were null genes that starvation never changed (section 6) |
| If the signs are opposite, it moved closer to Ctrl | Overshoot ($\Delta_S=-1,\ \Delta_R=+2.5 \Rightarrow \Delta_E=+1.5$) also satisfies both signs |
| If $\lvert\hat\Delta_E\rvert<\lvert\hat\Delta_S\rvert$, it moved closer to Ctrl | 51% of null genes and 53% of starvation_only genes also satisfied it (section 6). Comparisons between estimates carry their own uncertainty |
| A negative correlation between $\hat\Delta_S$ and $\hat\Delta_R$ across genes is evidence of a rescue mechanism | Even fully null genes give a correlation of −0.5 because of the shared $S$ (−0.53 in the paired design too) |

</details>

<details>
<summary>Where the results differ from the study text</summary>

| Study text claim | Result | Evidence |
|---|---|---|
| 11.5 The screening/confirmation structure of stageR and its gene-level FDR target | Not confirmed (by running) | stageR not installed. The procedure is described from the paper |
| 14.1 The third condition adds rescue on top of the treatment (here Starvation+Glucose); `~pair+condition` with donor pairing | Matches (definition of the design) | The main data of this note are unpaired `~condition`; a `~pair+condition` fit is run in the paired block of section 5 |
| 14.1 $\Delta_E=\Delta_S+\Delta_R$ | Matches (with exceptions) | 1.8e-15 by coefficient arithmetic; in 9 contrasts where both groups are all-zero, the `contrastAllZero` rule (lfc=0, p=1) breaks the additive relationship |
| 14.3 DESeq2 supports `altHypothesis="lessAbs"`, and the study text's code runs | Matches (with an addition) | Runs; p is the TOST `pmax(pvalueAbove, pvalueBelow)`. The source requires `lfcThreshold>0` and `betaPrior=FALSE` (not mentioned in the study text; met with the defaults) |
| 14.4 BH's dependence conditions and the Wald approximation are still needed; state the selection process | Not confirmed (a statement) | A statement about conditions, so not run separately |
| 14.6 Report an effect-size criterion, uncertainty and donor consistency together; the ratio of two shrunken LFCs is not a joint posterior ratio | Not confirmed (by running) | apeglm and ashr not installed. A statement based on the general fact $E[a]/E[b]\ne E[a/b]$ |
| Appendix B, answer 17: if p is NA too, it is all-zero, a count outlier or a fitting problem | Partly matches | All-zero and Cook's confirmed (section 3). Of the "fitting problems", non-convergence of β does not make p NA: `nbinomWaldTest` only flags `betaConv` and prints a message, and the only line in `results()` that puts NA into `pvalue` is `res$pvalue[which(cooksOutlier)] <- NA` (grep on `deparse`). Rows whose design degenerates because of weights get NA even for p. Both cases are run in the "When even p is NA" block of [07](07_lfc_shrinkage_and_qc.md) |

The rest of the study text's explanations, not in the table, all held up when run directly or worked out from the definitions.

</details>

---

← Previous: [05. The Wald test and the LRT](05_wald_vs_lrt.md) · Next: [07. LFC shrinkage and QC](07_lfc_shrinkage_and_qc.md) →

# 08. How do six counts become one p-value?

I looked at the dispersion of [03](03_dispersion_estimation.md), the GLM of [04](04_glm_condition_batch.md) and the Wald test of [05](05_wald_vs_lrt.md) one at a time, but never followed six counts all the way to a single p-value. So with gene A from chapter 16 of the study text I compute α → SE → Z → p by hand, and put the same counts into actual DESeq2 to compare what is the same and what differs. Then I run the three-group paired practice code of chapter 15 of the study text on data with known true effects, score the output, and go as far as testing the two-direction claim "it goes down and comes back".

> Study text chapters 15 and 16

## 1. If the mean doubles, can a conclusion be drawn right away?

Gene A of chapter 16 of the study text gave the counts below in six samples. A count is the number of reads assigned to the gene in one sample. Because sequencing depth differs between samples, a multiplier that corrects for it, the size factor, is used; in this example they are all 1, so the counts can be compared as they are.

| Condition | count | Mean |
|---|---|---:|
| Ctrl | 100, 130, 90 | 106.667 |
| Starvation | 200, 250, 180 | 210 |

The Starvation mean is about twice the Ctrl mean. log2 of the ratio of the two means is called the log2 fold change (LFC), here log2(210 / 106.667) = 0.977. An LFC of 1 means 2-fold, and −1 means half.

But the three Ctrl samples, under exactly the same condition, still wobble from 90 to 130. To judge whether a twofold difference could arise from such wobble alone, we need the SE (standard error), how much the LFC would move if the same experiment were repeated. And the SE depends on the dispersion α, the spread between replicates, so the calculation starts from α.

Chapter 16 of the study text splits this calculation into five steps, and this note follows the same order. The concepts were covered in earlier notes, so I only touch on them briefly.

| Step | What it asks | Result for gene A | Note that first covered it (study text chapter) |
|---|---|---|---|
| ① maximize the NB likelihood | How large is the spread α between replicates? | α = 0.014786 | [02](02_negative_binomial.md) (2.3), [03](03_dispersion_estimation.md) (chapter 5) |
| ② Cox-Reid adjustment | What if the loss from estimating the mean from the same data is made up? | α = 0.025385 | [03](03_dispersion_estimation.md) (chapter 6) |
| ③ MAP with a prior | What if the tendency of other genes is also taken into account? | α = 0.053147 | [03](03_dispersion_estimation.md) (chapter 7) |
| ④ weights and SE | At this α, how much does the LFC wobble? | SE = 0.289057 (log2) | [04](04_glm_condition_batch.md) (chapters 4 and 8) |
| ⑤ Wald test | How many SEs does the LFC lie from 0? | Z = 3.3809, p = 0.00072244 | [05](05_wald_vs_lrt.md) (chapter 9) |

This example is easy to solve by hand because the design has only two columns, the intercept and a Starvation indicator, and the size factors are all 1. In this case each group's mean estimate equals the arithmetic group mean (106.667, 210) whatever α is, so changing α leaves the LFC unchanged and changes only the SE. That lets the effect of α on the SE be seen on its own. But, as section 16.1 of the study text also adds, this does not mean general designs with offsets, pairs or batches also solve as arithmetic means, and section 5 shows an example where they actually do not.

Section 16.1 of the study text also notes that these numbers are not output from running DESeq2. So section 4 compares them with actual DESeq2. The chapter 15 practice in sections 5–6 also brings in one-sided p-values, joint p-values and BH ([06](06_multiple_testing.md), study text chapters 11 and 14), and size factors, the mu assay and outliers ([00](00_overview.md), [07](07_lfc_shrinkage_and_qc.md), study text chapter 13).

The sentences with which the study text sums up the two chapters are the backbone of this note.

> Chapter 15: "`dds` is not just a results table but a record in which the counts, design, normalization, dispersion and coefficient estimation are linked. It stores the inputs and outputs of each step together."
>
> Chapter 16: "Likelihood maximization, bias adjustment, prior shrinkage, coefficient uncertainty and hypothesis testing are connected but distinct operations. You should be able to trace this order in a single numerical example."

The code was run with `Rscript` in a working folder outside the repository, sections 2–4 and sections 5–6 each in one R session from top to bottom (including the code in the collapsible blocks, in the order it appears), so later blocks reuse objects from earlier blocks (`k`, `cr()`, `a_map`, `se_from_alpha()`, `dds` and so on).

## 2. What should the spread α between replicates be set to?

### A variance larger than the Poisson predicts

In a Poisson distribution the variance equals the mean. If gene A followed a Poisson, the Ctrl variance should be about 107 and the Starvation variance around 210. Let me compute the actual variances.

```r
options(digits = 8)
k <- c(100, 130, 90, 200, 250, 180)            # gene A: 3 Ctrl, 3 Starvation
grp <- rep(c("Ctrl", "Starvation"), each = 3)
m <- tapply(k, grp, mean); v <- tapply(k, grp, var)
tab1 <- rbind(mean = m, var = v, poisson_var = m, moment_alpha = (v - m) / m^2)
print(noquote(formatC(tab1, format = "f", digits = 4)), right = TRUE)
cat("LFC = log2(210 / 106.667) =", log2(m[["Starvation"]] / m[["Ctrl"]]), "\n")
```

```
                 Ctrl Starvation
mean         106.6667   210.0000
var          433.3333  1300.0000
poisson_var  106.6667   210.0000
moment_alpha   0.0287     0.0247
LFC = log2(210 / 106.667) = 0.97727992
```

Comparing the `var` row with the `poisson_var` row, the actual variance is 4–6 times the Poisson prediction.

When the spread between replicates of the same condition is larger than the Poisson predicts like this, it is called overdispersion. The negative binomial (NB) is a count distribution that allows overdispersion, with variance μ + αμ². μ is the expected count, the mean count the model expects in that sample, and α is the dispersion. α is not the variance itself but the size of the spread piled on top of the Poisson variance μ, so α = 0 is the same as the Poisson ([02](02_negative_binomial.md)).

Solving the variance formula for α gives α = (variance − μ) / μ². For Ctrl that is (433.3 − 106.7) / 106.7² = 0.0287, the value in the `moment_alpha` row of the output. Both groups are in the range 0.025–0.029.

### Finding α with the likelihood

This formula is good for getting a feel, but DESeq2 uses the likelihood. The likelihood, with the observed counts held fixed, says how plausible a candidate α makes those counts. The α with the largest value is the MLE (maximum likelihood estimate).

For example, setting α = 0.0148 and adding up the logs of the probabilities of the six counts gives −27.03, and α = 0.08 gives −28.97 (grid table below). The larger one, 0.0148, explains the data better. α is greater than 0 and changes over orders of magnitude, so it is searched as θ = log α, with the means μ_j fixed at the group means. As a formula:

$$\ell(\theta) = \sum_{j=1}^{6} \log \mathrm{NB}\!\left(K_j;\ \mu_j,\ \text{size} = e^{-\theta}\right)$$

- $K_j$ is the count of sample $j$: Ctrl 100, 130, 90 and Starvation 200, 250, 180.
- $\mu_j$ is the expected count of sample $j$: 106.667 for Ctrl and 210 for Starvation.
- $\text{size} = e^{-\theta} = 1/\alpha$ is the notation used by R's `dnbinom()`.
- $\log \mathrm{NB}(\cdot)$ is the log of the probability of the count $K_j$ at that α, and the sum of the six is $\ell$.

### The price of estimating the mean from the same data: the Cox-Reid adjustment

The moment calculation gave 0.025–0.029, but maximizing $\ell$ gives a smaller 0.0148. This is because the means μ_j were estimated from the same six counts. A group mean always sits right in the middle of that group's counts, so the differences between counts and mean look smaller than the real spread, and α is estimated too small. The sample variance `var()` makes up this loss by dividing by n − 1 instead of n, but maximizing $\ell$ has no such device. The Cox-Reid adjustment is exactly the correction term that reduces this underestimation.

$$\ell_{CR}(\theta) = \ell(\theta) - \tfrac12 \log\det\!\left(X^\top W X\right), \qquad W = \mathrm{diag}\!\left(\frac{\mu_j}{1 + \alpha\mu_j}\right)$$

- $X$ is the design matrix, a 6×2 table recording in numbers which condition each sample is in: the first column is all 1 (intercept), and the second is 1 for Starvation and 0 for Ctrl.
- $W$ is the matrix with each sample's weight $\mu_j/(1+\alpha\mu_j)$ on its diagonal. The weight is the amount of information the sample gives about the mean; it comes up again in section 3.
- $X^\top W X$ is the 2×2 matrix that gathers the information about the mean coefficients. In this design its determinant (det) is (sum of Ctrl weights) × (sum of Starvation weights).
- $-\tfrac12 \log\det(\cdot)$ is the correction term. The smaller α, the larger the weights and the information, and so the more is subtracted, which penalizes small α and raises the estimate.

### Borrowing the tendency of other genes: the prior and the MAP

With only three replicates, the estimate of α wobbles a lot. So DESeq2 reduces this wobble by consulting the α of other genes with similar mean expression. The curve showing where the dispersion roughly lies as a function of mean expression is called the trend, and the distribution for where α is likely to be before the data are seen is called the prior. DESeq2 estimates a prior centred on the trend from all genes (empirical Bayes) and finds the MAP (maximum a posteriori), the most plausible value when the likelihood and the prior are considered together. Pulling estimates that carry little information toward the overall tendency like this is called shrinkage ([03](03_dispersion_estimation.md)).

But the trend cannot be estimated from one gene. So section 16.3 of the study text sets a teaching prior directly: centre log 0.08 and standard deviation 0.7, and the study text states that "these values were not estimated from real data".

$$\log \text{post}(\theta) = \ell_{CR}(\theta) - \frac{(\theta - \log 0.08)^2}{2 \cdot 0.7^2}$$

- $\log 0.08$ is the centre of the prior, the trend value assumed to have been obtained from other genes.
- $0.7$ is the standard deviation of the prior on the log α scale. The larger it is, the weaker the prior.
- The second term is the penalty from the prior. It is 0 at α = 0.08 and grows quadratically as log α moves away from log 0.08.

### The three α values with the study text's code

```r
# the code of study text 16.6 (optimize() default tol)
X <- cbind(Intercept=1, Starvation=c(0,0,0,1,1,1))
mu <- rep(c(mean(k[1:3]), mean(k[4:6])), each=3)
ll <- function(theta) sum(dnbinom(k, mu=mu, size=exp(-theta), log=TRUE))
cr <- function(theta) {
  a <- exp(theta); w <- mu / (1 + a*mu)
  info <- crossprod(X, sweep(X, 1, w, `*`))
  ll(theta) - 0.5*as.numeric(determinant(info, logarithm=TRUE)$modulus)
}
post <- function(theta) cr(theta) - 0.5*((theta-log(0.08))/0.7)^2
bounds <- log(c(1e-5, 2))
a_mle <- exp(optimize(ll, bounds, maximum=TRUE)$maximum)
a_cr  <- exp(optimize(cr, bounds, maximum=TRUE)$maximum)
a_map <- exp(optimize(post, bounds, maximum=TRUE)$maximum)
cat("mu =", mu, "\n"); cat(sprintf("a_mle = %.6f\na_cr  = %.6f\na_map = %.6f\n", a_mle, a_cr, a_map))
```

```
mu = 106.66667 106.66667 106.66667 210 210 210 
a_mle = 0.014786
a_cr  = 0.025385
a_map = 0.053147
```

All three values agree with sections 16.2–16.3 of the study text to the sixth decimal place. But the calculations after this use the value converged to the end with a tighter tolerance for `optimize()` (α_MAP = 0.053147342), because the Z and p of section 3 respond to this difference even though the sixth decimal place is the same.

<details>
<summary>The tolerance of optimize() and the study text's numbers</summary>

The default `tol` of `optimize()` is `.Machine$double.eps^0.25` = 1.2e-4 on the θ scale, which is loose for reporting 6 digits. The block below recomputes `a_mle`, `a_cr` and `a_map` with `tol = 1e-12` and overwrites them.

```r
# optimize()'s default tol is .Machine$double.eps^0.25 = 1.2e-4 on the θ scale, loose for reporting 6 digits.
a_map_def <- a_map
a_mle <- exp(optimize(ll, bounds, maximum=TRUE, tol=1e-12)$maximum)
a_cr  <- exp(optimize(cr, bounds, maximum=TRUE, tol=1e-12)$maximum)
a_map <- exp(optimize(post, bounds, maximum=TRUE, tol=1e-12)$maximum)   # later calculations use this converged value
cat(sprintf("tol=1e-12: %.6f %.6f %.6f\n", a_mle, a_cr, a_map))
dpost <- function(a, h=1e-6) (post(log(a)+h) - post(log(a)-h)) / (2*h)  # d log post / dθ
cat(sprintf("a_map: default tol %.9f (gradient %.1e) | tol=1e-12 %.9f (gradient %.1e)\n", a_map_def, dpost(a_map_def), a_map, dpost(a_map)))
```

```
tol=1e-12: 0.014786 0.025385 0.053147
a_map: default tol 0.053146542 (gradient 4.5e-05) | tol=1e-12 0.053147342 (gradient 3.6e-09)
```

The α_MAP obtained with the study text's code as it is is 0.053146542, and the value converged with `tol = 1e-12` is 0.053147342. The slopes are 4.5e-05 vs 3.6e-09, so the former had not yet reached the maximum. Rounded to 6 digits both are 0.053147, so α alone cannot tell them apart, but the Z and p of 16.4 respond to this difference (see "Comparing with the study text's numbers digit by digit" in section 3). So every calculation after this block uses the converged value.

</details>

Plugging α_MAP into the three objectives shows the size of each term.

```r
a <- 0.053147342                                 # converged α_MAP (see the collapsible block above)
w <- mu / (1 + a*mu); info <- crossprod(X, X * w)
cat("w =", round(w, 4), "\n")
cat("det(XtWX) =", det(info), "  sum(w_Ctrl) * sum(w_Starvation) =", sum(w[1:3]) * sum(w[4:6]), "\n")
cat("ll =", ll(log(a)), "  -0.5*log det =", -0.5*log(det(info)), "  ll_CR =", cr(log(a)), "\n")
cat("prior penalty =", -0.5*((log(a) - log(0.08))/0.7)^2, "  log post =", post(log(a)), "\n")
```

```
w = 15.9943 15.9943 15.9943 17.2684 17.2684 17.2684 
det(XtWX) = 2485.7609   sum(w_Ctrl) * sum(w_Starvation) = 2485.7609 
ll = -28.182241   -0.5*log det = -3.909167   ll_CR = -32.091408 
prior penalty = -0.17066029   log post = -32.262068 
```

Subtracting the correction term 3.91 from $\ell$ = −28.18 gives $\ell_{CR}$ = −32.09, and subtracting the prior penalty 0.17 from that gives log post = −32.26. The penalty 0.17 is the difference between log 0.0531 and log 0.08, −0.409, squared and divided by 2 × 0.49. The weights 15.99 and 17.27 are used again when computing the SE in section 3.

### The three curves in one figure

![Gene A: three objective functions over alpha](../../figures/03_gene_a_objectives.png)

A figure of the three objectives each lowered so that its maximum is 0, taken as it is from [03](03_dispersion_estimation.md), where figure 3 of the study text was redrawn. Each curve has a different baseline, so heights must not be compared between curves.

The maximum moves to the right in the order MLE 0.0148 → CR 0.0254 → MAP 0.0531. The 0.0254 enlarged by the Cox-Reid adjustment came close to the mean of the two groups' moment values, 0.0267, but is not the same. The adjustment is an approximation, and the study text also writes in 16.2 that "this does not mean Cox-Reid is a formula that always increases α by the same proportion in every dataset" (the principle of the adjustment is in section 3 of [03](03_dispersion_estimation.md)). Attaching the prior makes α grow further to 0.0531, between the gene-wise value (0.0254, found from this gene's data alone) and the prior centre 0.08. At first I thought shrinkage would pull the value down, but it actually pulls toward the trend, so gene A, which was below the trend, was pulled up.

The region near the maximum is flat. For α between 0.01 and 0.03, $\ell$ is at most 0.36 below its maximum and $\ell_{CR}$ at most 0.54 below. Six counts cannot pin α down precisely. Also, only the curve with the prior is steep on the left, but not because the prior is asymmetric. The penalty is set only by the log distance from log 0.08: the left end of the grid, α = 0.001, is 4.38 log units from the centre and the right end, α = 0.5, 1.83 log units, so the penalties diverge, 19.59 vs 3.43.

<details>
<summary>Figure 3 in numbers: the three objectives over a grid of α</summary>

```r
grid <- c(0.001,0.002,0.005,0.01,0.014786,0.02,0.025385,0.03,0.05,0.053147,0.08,0.1,0.2,0.5)
th <- log(grid)
tab <- data.frame(alpha=grid, ll=sapply(th, ll), cr=sapply(th, cr), post=sapply(th, post))
tab$ll_rel <- tab$ll - ll(log(a_mle)); tab$cr_rel <- tab$cr - cr(log(a_cr)); tab$post_rel <- tab$post - post(log(a_map))
print(round(tab, 4))
cat("|log(alpha) - log(0.08)| at 0.001, 0.5:", abs(log(c(0.001, 0.5)) - log(0.08)), "\n")
cat("prior penalty at 0.001, 0.5:", 0.5*((log(c(0.001, 0.5)) - log(0.08))/0.7)^2, "\n")
cat("shifted ll, cr, post at alpha = 0.6:", round(c(ll(log(0.6)) - ll(log(a_mle)), cr(log(0.6)) - cr(log(a_cr)), post(log(0.6)) - post(log(a_map))), 2), "\n")
cat("alpha where shifted post = -12:", exp(uniroot(function(t) post(t) - post(log(a_map)) + 12, log(c(1e-3, a_map)))$root), "\n")
```

```
    alpha       ll       cr     post  ll_rel  cr_rel post_rel
1  0.0010 -29.5421 -35.5031 -55.0972 -2.5145 -3.7337 -22.8351
2  0.0020 -28.8493 -34.6843 -48.5699 -1.8218 -2.9149 -16.3078
3  0.0050 -27.7393 -33.2736 -41.1178 -0.7117 -1.5042  -8.8557
4  0.0100 -27.1330 -32.3114 -36.7237 -0.1055 -0.5420  -4.4616
5  0.0148 -27.0275 -31.9551 -34.8638  0.0000 -0.1857  -2.6017
6  0.0200 -27.0938 -31.8055 -33.7665 -0.0663 -0.0361  -1.5045
7  0.0254 -27.2403 -31.7694 -33.1139 -0.2128  0.0000  -0.8518
8  0.0300 -27.3912 -31.7868 -32.7684 -0.3637 -0.0174  -0.5064
9  0.0500 -28.0793 -32.0423 -32.2677 -1.0518 -0.2729  -0.0056
10 0.0531 -28.1822 -32.0914 -32.2621 -1.1547 -0.3220   0.0000
11 0.0800 -28.9677 -32.5078 -32.5078 -1.9402 -0.7383  -0.2457
12 0.1000 -29.4579 -32.7910 -32.8418 -2.4303 -1.0216  -0.5798
13 0.2000 -31.1914 -33.8648 -34.7215 -4.1639 -2.0954  -2.4595
14 0.5000 -33.8391 -35.6168 -39.0437 -6.8115 -3.8474  -6.7816
|log(alpha) - log(0.08)| at 0.001, 0.5: 4.3820266 1.8325815 
prior penalty at 0.001, 0.5: 19.594038 3.4268927 
shifted ll, cr, post at alpha = 0.6: -7.38 -4.23 -7.88 
alpha where shifted post = -12: 0.0033178508 
```

Rows 5, 7 and 10 are α_MLE, α_CR and α_MAP respectively (`*_rel` = 0 for that curve), and row 11 (α = 0.08) is the prior centre. As in the caption of figure 3 of the study text, each curve is shifted so that its own maximum is 0, so the `*_rel` columns must not be used to compare heights between curves. The maxima of the absolute columns (`ll`, `cr`, `post`) are −27.0275, −31.7694 and −32.2621, but these too are values of different objectives, so comparing sizes between curves means nothing. The absolute values are only useful for relationships within one row. For example, at α = 0.08 the prior penalty is 0, so `post` and `cr` are equal (−32.5078).

In the 0.01–0.03 range (rows 4–8), the lowest `ll_rel` is −0.3637 and the lowest `cr_rel` −0.5420. The last two lines are values picked to overlay on the study text's figure. At α = 0.6 the three curves are −7.38, −4.23 and −7.88, and the prior curve drops below −12 at α ≈ 0.0033. Comparing by eye with the right end and lower-left shape of figure 3 of the study text (p.35), they matched these values. The code for the figure is in [03](03_dispersion_estimation.md).

</details>

## 3. Once α is set, how do the SE and p-value come out?

### Weight: the information one sample gives

Once α is set, a weight can be computed for each sample. The weight is the amount of information one sample gives about its group mean. Under a Poisson (α = 0) the weight would be μ itself, 106.7 for a Ctrl sample and 210 for a Starvation sample, so more counts mean more information. But with α = 0.0531, Ctrl becomes 106.667 / (1 + 0.0531 × 106.667) = 15.99 and Starvation 17.27. The mean is twice as large, yet the information is almost the same.

$$w_j = \frac{\mu_j}{1 + \alpha\mu_j}$$

Even as μ grows, $w_j$ only approaches 1/α = 18.8 and grows no further. Under the NB, no matter how many counts there are, the relative spread between replicates (CV² = 1/μ + α) does not go below α.

### SE: the sum of the reciprocals of the per-group weight sums

Adding the weights of the three Ctrl samples gives 3 × 15.994 = 47.98, and Starvation 3 × 17.268 = 51.81. Adding the reciprocals of the two sums gives 1/47.98 + 1/51.81 = 0.02084 + 0.01930 = 0.04014, whose square root is 0.2004; divided by ln 2 = 0.693 it is 0.2891. This is gene A's SE.

This calculation comes from the GLM. The GLM (generalized linear model) is the model that writes the mean of the counts on the log scale as a sum of condition, batch and so on ([04](04_glm_condition_batch.md)). In this example log(mean) = $b_0$ + $b_1$ × (1 if Starvation), $b_1$ is the Starvation effect in natural-log units, and the LFC $\beta$ in log2 units is $b_1 / \ln 2$. In a two-group model the SE of $b_1$ becomes a simple formula.

$$\mathrm{SE}(\hat b_1) = \sqrt{\frac{1}{\sum_{\text{Ctrl}} w_j} + \frac{1}{\sum_{\text{Starvation}} w_j}}, \qquad \mathrm{SE}(\hat\beta) = \frac{\mathrm{SE}(\hat b_1)}{\ln 2}$$

For each group the reciprocal of its information (weight sum) is the variance of that group's mean (on the log scale), and the variance of the difference between the two groups is the sum of the two.

In a general design, instead of this formula, the inverse of the information matrix seen in section 2, $(X^\top W X)^{-1}$, is used. This is the covariance matrix of the coefficients, and the square roots of its diagonal elements are the SEs of the coefficients. In a two-group model the diagonal element at the Starvation coefficient equals the value under the square root in the formula above.

### Wald test: how many SEs from 0?

Dividing the LFC by the SE gives Z = 0.977280 / 0.289057 ≈ 3.381, meaning the LFC lies 3.4 SEs away from 0. The Wald test compares (estimate − 0) / SE with the standard normal like this ([05](05_wald_vs_lrt.md)). Here p ≈ 0.00072, and the 95% CI is 0.977 ± 1.96 × 0.289 = [0.411, 1.544]. The study text calls this interval a nominal CI, meaning one without adjustment for looking at many genes together.

$$Z = \frac{\hat\beta}{\mathrm{SE}(\hat\beta)}, \qquad p = 2\,\Phi(-|Z|), \qquad \text{95\% CI} = \hat\beta \pm 1.96\,\mathrm{SE}(\hat\beta)$$

$\Phi$ is the cumulative distribution function of the standard normal, so $2\Phi(-|Z|)$ is the probability of a value more extreme than $|Z|$ in either tail. 1.96 is the value that cuts off the middle 95% of the standard normal. Let me run the calculation so far as code.

```r
se_from_alpha <- function(a) {
  w <- mu / (1 + a*mu); info <- crossprod(X, sweep(X, 1, w, `*`)); cov_nat <- solve(info)
  list(w=w, info=info, cov=cov_nat, se_nat=sqrt(cov_nat[2,2]), se_log2=sqrt(cov_nat[2,2])/log(2))
}
r <- se_from_alpha(a_map)                      # a_map = 0.053147342 (converged value)
cat("W diag:", r$w, "\n"); print(r$cov)
cat("1/sum(w_Ctrl) + 1/sum(w_Starvation) =", 1/sum(r$w[1:3]) + 1/sum(r$w[4:6]), "\n")
cat("SE (natural log) =", r$se_nat, "  SE (log2) =", r$se_log2, "\n")
lfc <- log2(210/(320/3)); Z <- lfc / r$se_log2; p <- 2*pnorm(-abs(Z))
cat(sprintf("LFC=%.6f  SE_log2=%.7f  Z=%.7f  p=%.8f\n", lfc, r$se_log2, Z, p))
cat(sprintf("95%% CI (z=1.96): [%.6f, %.6f]\n", lfc - 1.96*r$se_log2, lfc + 1.96*r$se_log2))
```

```
W diag: 15.994282 15.994282 15.994282 17.268399 17.268399 17.268399
              Intercept   Starvation
Intercept   0.020840781 -0.020840781
Starvation -0.020840781  0.040143863
1/sum(w_Ctrl) + 1/sum(w_Starvation) = 0.040143863
SE (natural log) = 0.20035934   SE (log2) = 0.28905742
LFC=0.977280  SE_log2=0.2890574  Z=3.3809197  p=0.00072244
95% CI (z=1.96): [0.410727, 1.543832]
```

The bottom-right value of the covariance matrix, 0.040143863, equals the hand calculation on the line right below it. The last two lines also agree with the values of section 16.4 of the study text (LFC 0.977280, SE 0.289057, Z 3.380919, p 0.00072244, CI [0.410727, 1.543833]) at the rounding level.

Section 16.4 of the study text does not produce a padj here. padj (the adjusted p-value) is the p-value after a multiple-testing adjustment together with other genes, and the study text gives the reason: "a real BH padj can only be computed with the p-values of the other genes and the set of adjustment targets". A real padj comes out in section 4.

<details>
<summary>Comparing with the study text's numbers digit by digit</summary>

```r
print(r$info)
cat(sprintf("CI z=1.96: [%.6f, %.6f]   CI z=qnorm(0.975): [%.6f, %.6f]\n",
            lfc-1.96*r$se_log2, lfc+1.96*r$se_log2, lfc-qnorm(0.975)*r$se_log2, lfc+qnorm(0.975)*r$se_log2))
cat(sprintf("textbook fraction 0.977280/0.289057 = %.6f   2*0.977280 - 0.410727 = %.6f\n", 0.977280/0.289057, 2*0.977280 - 0.410727))
# comparison: a_map with the default tol, and the 6-digit rounded 0.053147
for (a in c(a_map_def, 0.053147)) { s <- se_from_alpha(a)$se_log2; z <- lfc/s
  cat(sprintf("alpha=%.9f: SE=%.7f Z=%.6f p=%.8f CI(1.96)=[%.6f, %.6f]\n", a, s, z, 2*pnorm(-abs(z)), lfc-1.96*s, lfc+1.96*s)) }
```

```
           Intercept Starvation
Intercept  99.788045  51.805198
Starvation 51.805198  51.805198
CI z=1.96: [0.410727, 1.543832]   CI z=qnorm(0.975): [0.410738, 1.543822]
textbook fraction 0.977280/0.289057 = 3.380925   2*0.977280 - 0.410727 = 1.543833
alpha=0.053146542: SE=0.2890555 Z=3.380942 p=0.00072238 CI(1.96)=[0.410731, 1.543829]
alpha=0.053147000: SE=0.2890566 Z=3.380929 p=0.00072241 CI(1.96)=[0.410729, 1.543831]
```

The table below collects the output of the two blocks above.

| Quantity | Study text | Recomputed (converged α_MAP = 0.053147342) | Reference: study text code, default tol (α = 0.053146542) | Verdict |
|---|---:|---:|---:|---|
| α_MLE / α_CR / α_MAP | 0.014786 / 0.025385 / 0.053147 | 0.014786 / 0.025385 / 0.053147 | Same (6 digits) | Matches |
| LFC | 0.977280 | 0.977280 | 0.977280 | Matches |
| SE (log2) | 0.289057 | 0.2890574 | 0.2890555 | Matches |
| Z | 3.380919 | 3.3809197 | 3.380942 | Matches (differs by 1 in the last digit, rounding level) |
| p | 0.00072244 | 0.00072244 | 0.00072238 | Matches |
| 95% CI | [0.410727, 1.543833] | z = 1.96: [0.410727, 1.543832]; qnorm(0.975): [0.410738, 1.543822] | [0.410731, 1.543829] (1.96) | Matches (the upper end differs by 1 in the last digit) |

The study text's Z, p and CI agree with values computed from the converged α_MAP. The difference I recorded in an earlier version of this note as a "mismatch in the 5th and 7th decimal places" was not the study text's fault; it was because this note used the α from `optimize()`'s default tol (0.053146542). Two minor points remain.

1. Dividing the fraction written in the study text's formula, 0.977280 / 0.289057, as it is gives 3.380925. The displayed SE is a rounded value, while Z was computed with the SE before rounding.
2. The CI seems to use z = 1.96, since qnorm(0.975) gives [0.410738, 1.543822], further off. The study text's upper end 1.543833 equals 2 × 0.977280 − 0.410727 = 1.543833, but whether the study text computed it this way cannot be known.

Even plugging in α rounded to 6 digits, 0.053147, gives Z = 3.380929, p = 0.00072241, off from the study text (last line). To report 6 digits, first make sure of convergence by tightening the tol of `optimize()`.

</details>

### What changes if only α changes?

Let me repeat the same calculation with four values of α (MLE, Cox-Reid, MAP, and 0, corresponding to the Poisson).

```r
res <- t(sapply(c(a_mle, a_cr, a_map, 0), function(a) { rr <- se_from_alpha(a); z <- lfc/rr$se_log2
  c(alpha=a, w_ctrl=rr$w[1], w_starvation=rr$w[4], SE_log2=rr$se_log2, Z=z, p=2*pnorm(-abs(z))) }))
print(format(as.data.frame(res), digits=6))
cat("p ratio Poisson/MAP:", res[4,"p"]/res[3,"p"], "  MAP/MLE:", res[3,"p"]/res[1,"p"], "\n")
```

```
      alpha   w_ctrl w_starvation   SE_log2       Z           p
1 0.0147860  41.3890      51.1564 0.1741401 5.61203 1.99965e-08
2 0.0253849  28.7688      33.1710 0.2122064 4.60533 4.11819e-06
3 0.0531473  15.9943      17.2684 0.2890574 3.38092 7.22437e-04
4 0.0000000 106.6667     210.0000 0.0990355 9.86797 5.73117e-23
p ratio Poisson/MAP: 7.9331117e-20   MAP/MLE: 36128.058
```

The LFC is the same, 0.977280, in all four cases; what changes is only the weight, SE and p. Wrongly assuming a Poisson (row 4) shrinks the SE to 0.099, making p about 10⁻¹⁹ times the MAP's, and using the MAP instead of the MLE makes p about 3.6 × 10⁴ times larger. The effect size stays the same and only the strength of evidence changes. This is what section 16.5 of the study text sets out to show.

For gene A, shrinkage enlarged α and so p grew too, but the opposite direction also happens. Genes whose gene-wise value is above the trend get a lower α (gene11 in section 5: 0.821 → 0.480). If the mean stays the same, the weights grow, the SE shrinks and p gets smaller too. So which way shrinkage moves p differs from gene to gene.

The study text places this example as one that singles out just the dispersion path of the chapter 12 question "the difference between putting in only two groups and putting in three groups together". That question itself was covered in section 6 of [05](05_wald_vs_lrt.md). The study text also adds "in real datasets normalization and β can change too, so this example should not be generalized into a numerical result for every situation", and this note has two such examples. In section 4 (A), estimating the size factors from the data changes the LFC from 0.977 to 0.874, and in the pair design of section 5, when α changes, the means and β move along with it.

## 4. When the same counts go into actual DESeq2, what is the same and what differs?

### Slipping gene A in among background genes

Gene A must not be put into DESeq2 as a single row. Section 15.3 of the study text gives the reason: "If you pick a small gene set arbitrarily and fit DESeq2 from scratch, the information for learning the size factors and the dispersion trend changes. Usually you fit with an appropriate full gene universe and extract only the genes of interest from the results." So I slipped gene A in, as a single row called `teach_gene`, among 2000 background genes made with `makeExampleDESeqDataSet(n = 2000, m = 6, betaSD = 1)` (seed 1). The condition names A/B of the example data were changed to Ctrl/Starvation.

It was run in three ways. (A) is the default `DESeq()`, which estimates even the size factors from the data, and (B) fixes the size factors at 1 to match the study text's means. (C) takes (B), changes only gene A's α to the converged α_MAP of section 2 (0.053147342), and redoes only the Wald test.

```r
suppressPackageStartupMessages(library(DESeq2)); options(digits = 7)
set.seed(1)
dds2 <- makeExampleDESeqDataSet(n = 2000, m = 6, betaSD = 1)        # 2000 background genes
levels(dds2$condition) <- c("Ctrl", "Starvation")
cts <- rbind(counts(dds2), teach_gene = c(100L,130L,90L,200L,250L,180L))   # add gene A as one row
dds2 <- DESeqDataSetFromMatrix(cts, colData(dds2)[, "condition", drop=FALSE], design = ~ condition)
ddsA <- DESeq(dds2, quiet = TRUE)                                              # (A) default DESeq()
ddsB <- dds2; sizeFactors(ddsB) <- rep(1, 6); ddsB <- DESeq(ddsB, quiet = TRUE) # (B) size factor = 1
ddsC <- ddsB                                                                   # (C) B with only α set to the study text's value
dispersions(ddsC)[rownames(ddsC) == "teach_gene"] <- a_map                    # a_map = 0.053147342
ddsC <- nbinomWaldTest(ddsC, betaPrior = FALSE)
one <- function(d) { m <- mcols(d)["teach_gene", ]; r <- results(d)["teach_gene", ]
  c(dispGeneEst = m$dispGeneEst, dispFit = m$dispFit, dispMAP = m$dispMAP, dispersion = m$dispersion,
    log2FC = r$log2FoldChange, lfcSE = r$lfcSE, stat = r$stat, pvalue = r$pvalue) }
print(signif(sapply(list(A = ddsA, B = ddsB, C = ddsC), one), 6))
```

```
found results columns, replacing these
                    A          B           C
dispGeneEst 0.0266060 0.02538760 0.025387600
dispFit     0.1340590 0.13326300 0.133263000
dispMAP     0.0969760 0.09583870 0.095838700
dispersion  0.0969760 0.09583870 0.053147300
log2FC      0.8741090 0.97728000 0.977280000
lfcSE       0.3799470 0.37787800 0.289057000
stat        2.3006100 2.58623000 3.380920000
pvalue      0.0214138 0.00970315 0.000722434
```

The `found results columns, replacing these` on the first line is the message printed when `nbinomWaldTest()` is called again in (C), meaning the existing Wald result columns were replaced with new values.

First, the gene-wise estimate is the same as the study text's. The `dispGeneEst` of (B), 0.025388, equals the study text's Cox-Reid value 0.025385, and the 3 × 10⁻⁶ difference is optimization error arising because the objective is flat near its maximum (the first output line of "Reconstructing DESeq2's prior by hand" below: a difference of 7 × 10⁻⁹ in $\ell_{CR}$ between the two values). (A) is slightly different, 0.026606, because its size factors are not 1, so the means changed.

But the MAP differs. The `dispMAP` of (B), 0.0958, is larger than the study text's 0.0531, because both the centre and the width of the prior differ from the study text's. We look at this just below.

Fixing α at the study text's value makes the rest the same. (C)'s lfcSE 0.289057, stat 3.38092 and p 0.000722434 equal the hand calculation of section 3, and the log2FC is 0.977280 in both (B) and (C). The statement of 16.5, that changing α does not change the means, holds in actual fitting too.

<details>
<summary>The full output of (A)(B)(C) and a comparison table with the study text</summary>

```r
show_gene <- function(d, tag) {
  cat(sprintf("[%s] sizeFactors: %s\n", tag, paste(round(sizeFactors(d), 4), collapse=" ")))
  print(as.data.frame(mcols(d)["teach_gene", c("baseMean","dispGeneEst","dispFit","dispMAP","dispersion","dispOutlier","betaConv")]))
  print(as.data.frame(results(d)["teach_gene", ]))
  cat(sprintf("[%s] mu assay: %s\n", tag, paste(round(assays(d)[["mu"]]["teach_gene", ], 4), collapse=" ")))
  cat(sprintf("[%s] dispPriorVar: %s\n", tag, attr(dispersionFunction(d), "dispPriorVar")))
}
show_gene(ddsA, "A"); show_gene(ddsB, "B"); show_gene(ddsC, "C")
rC <- results(ddsC)["teach_gene", ]
cat(sprintf("[C] lfcSE=%.7f stat=%.7f p=%.9f padj=%.8f\n", rC$lfcSE, rC$stat, rC$pvalue, rC$padj))
```

```
[A] sizeFactors: 0.9917 0.9675 1.0027 1.0501 1.0703 1.0501
           baseMean dispGeneEst   dispFit    dispMAP dispersion dispOutlier
teach_gene 153.4003  0.02660597 0.1340589 0.09697603 0.09697603       FALSE
           betaConv
teach_gene     TRUE
           baseMean log2FoldChange     lfcSE     stat     pvalue      padj
teach_gene 153.4003      0.8741087 0.3799469 2.300608 0.02141381 0.1177086
[A] mu assay: 107.401 104.7743 108.5855 208.4326 212.4497 208.4417
[A] dispPriorVar: 0.25
[B] sizeFactors: 1 1 1 1 1 1
           baseMean dispGeneEst  dispFit   dispMAP dispersion dispOutlier
teach_gene 158.3333  0.02538756 0.133263 0.0958387  0.0958387       FALSE
           betaConv
teach_gene     TRUE
           baseMean log2FoldChange    lfcSE     stat      pvalue       padj
teach_gene 158.3333      0.9772803 0.377878 2.586232 0.009703152 0.06112221
[B] mu assay: 106.6666 106.6666 106.6666 210 210 210
[B] dispPriorVar: 0.25
[C] sizeFactors: 1 1 1 1 1 1
           baseMean dispGeneEst  dispFit   dispMAP dispersion dispOutlier
teach_gene 158.3333  0.02538756 0.133263 0.0958387 0.05314734       FALSE
           betaConv
teach_gene     TRUE
           baseMean log2FoldChange     lfcSE     stat       pvalue        padj
teach_gene 158.3333      0.9772801 0.2890574 3.380921 0.0007224337 0.009247152
[C] mu assay: 106.6666 106.6666 106.6666 210 210 210
[C] dispPriorVar: 0.25
[C] lfcSE=0.2890574 stat=3.3809208 p=0.000722434 padj=0.00924715
```

The output above, collected:

| | (A) default DESeq() | (B) size factor = 1 | (C) B + α fixed (converged α_MAP) | Study text by hand |
|---|---:|---:|---:|---:|
| sizeFactors | 0.9917 0.9675 1.0027 1.0501 1.0703 1.0501 | 1 (×6) | 1 (×6) | 1 |
| baseMean | 153.4003 | 158.3333 | 158.3333 | 158.333 |
| mu assay (teach_gene) | 107.401 104.7743 108.5855 208.4326 212.4497 208.4417 | 106.6666 ×3, 210 ×3 | 106.6666 ×3, 210 ×3 | 106.667 / 210 |
| dispGeneEst | 0.02660597 | 0.02538756 | 0.02538756 | 0.025385 (CR) |
| dispFit (trend) | 0.1340589 | 0.133263 | 0.133263 | 0.08 (given) |
| dispPriorVar | 0.25 | 0.25 | 0.25 | 0.49 (= 0.7²) |
| dispMAP → dispersion | 0.09697603 | 0.0958387 | 0.0958387 → 0.05314734 (overwritten) | 0.053147 |
| dispOutlier / betaConv | FALSE / TRUE | FALSE / TRUE | FALSE / TRUE | – |
| log2FoldChange | 0.8741087 | 0.9772803 | 0.9772801 | 0.977280 |
| lfcSE | 0.3799469 | 0.377878 | 0.2890574 | 0.289057 |
| stat | 2.300608 | 2.586232 | 3.380921 | 3.380919 |
| pvalue | 0.02141381 | 0.009703152 | 0.0007224337 | 0.00072244 |
| padj (BH over 2001 genes) | 0.1177086 | 0.06112221 | 0.009247152 | (not computed) |

That the padj in the last row differs between (A), (B) and (C) is because it depends not only on teach_gene's p but also on the ranks of the p-values of the other 2000 genes. This is the same reason section 16.4 of the study text did not produce a padj for a single gene.

(C)'s stat, 3.3809208, differs from the hand calculation 3.3809197 in the 7th digit. It is a numerical difference coming from the ridge λ = 1e-6 and the convergence level of IRLS (the iteratively reweighted least squares calculation that finds the GLM coefficients). The ratio of lfcSE to the hand-computed SE is 0.9999999 in (C) and 0.9999998 in (B) (block below).

</details>

### Reconstructing DESeq2's prior by hand

Let me rebuild by hand how the prior of (B) differs from the study text's.

```r
df <- dispersionFunction(ddsB); cf <- attr(df, "coefficients"); pv <- attr(df, "dispPriorVar"); vl <- attr(df, "varLogDispEsts")
mcB <- mcols(ddsB)["teach_gene", ]
cat("cr(log 0.0253849) - cr(log 0.0253876) =", cr(log(0.0253849)) - cr(log(0.0253876)), "\n")   # cr() is the function from the 16.2 block
cat("dispersionFunction coefs:", round(cf, 5), "  dispPriorVar:", pv, "  varLogDispEsts:", vl, "\n")
cat("trend by hand:", cf["asymptDisp"] + cf["extraPois"]/mcB$baseMean,
    "  varLogDispEsts - trigamma((6-2)/2):", vl - trigamma(2), " -> pmax(., 0.25):", pmax(vl - trigamma(2), 0.25), "\n")
post_deseq <- function(theta) cr(theta) - 0.5*(theta - log(mcB$dispFit))^2/pv   # only the two prior values replaced with DESeq2's
cat("hand MAP with DESeq2 prior:", exp(optimize(post_deseq, log(c(1e-8, 10)), maximum=TRUE, tol=1e-12)$maximum), " vs dispMAP:", mcB$dispMAP, "\n")
```

```
cr(log 0.0253849) - cr(log 0.0253876) = 7.227747e-09 
dispersionFunction coefs: 0.09311 6.35753   dispPriorVar: 0.25   varLogDispEsts: 0.6885003 
trend by hand: 0.133263   varLogDispEsts - trigamma((6-2)/2): 0.04356624  -> pmax(., 0.25): 0.25 
hand MAP with DESeq2 prior: 0.09583826  vs dispMAP: 0.0958387 
```

First, the trend. DESeq2 fits the trend curve α = `asymptDisp` + `extraPois` / mean to the background genes, and plugging in gene A's mean of 158.333 gives 0.09311 + 6.35753 / 158.333 = 0.133263. It is higher than the 0.08 the study text set.

The prior variance differs too. The value estimated from the data is `varLogDispEsts − trigamma(2)` = 0.0436. `trigamma((6 − 2)/2)` is subtracted because it is the share due to the wobble of the gene-wise estimate itself. But DESeq2 sets a lower limit of 0.25 on the prior variance, so here, where the estimate is smaller than that, 0.25 (standard deviation 0.5) was used. It is a narrower prior than the study text's 0.7² = 0.49.

Putting just these two values into the study text's `post()` gives 0.09583826. This equals DESeq2's `dispMAP` 0.0958387 within a relative difference of 5 × 10⁻⁶, so the operation is the same and only the prior differs.

<details>
<summary>The ratio of lfcSE to the hand-computed SE, and the dispersions() assignment trap</summary>

```r
se_log2 <- function(a) se_from_alpha(a)$se_log2                                  # function from the 16.4 block
cat(sprintf("hand SE_log2 at dispMAP %.7f = %.6f ; at a_map = %.7f\n", mcB$dispMAP, se_log2(mcB$dispMAP), se_log2(a_map)))
cat(sprintf("lfcSE / hand SE: [B] %.7f  [C] %.7f\n", results(ddsB)["teach_gene","lfcSE"]/se_log2(mcB$dispMAP), rC$lfcSE/se_log2(a_map)))
```

```
hand SE_log2 at dispMAP 0.0958387 = 0.377878 ; at a_map = 0.2890574
lfcSE / hand SE: [B] 0.9999998  [C] 0.9999999
```

When making (C), assigning by name, as in `dispersions(ddsC)["teach_gene"] <- a_map`, fails. `dispersions()` returns an unnamed vector, so the character index appends a new named element, making the length 2002. So the (C) block above used a logical index (`rownames(ddsC) == "teach_gene"`). `which()` works too.

```r
tmp <- ddsB
cat("names(dispersions(tmp)) is NULL:", is.null(names(dispersions(tmp))), "\n")
err <- try(dispersions(tmp)["teach_gene"] <- a_map, silent = TRUE)
cat("name index ->", conditionMessage(attr(err, "condition")), "\n")
```

```
names(dispersions(tmp)) is NULL: TRUE 
name index -> 2002 elements in value to replace 2001 elements 
```

</details>

### When the background genes change

The trend of 0.133 above is not a property of gene A alone but a value set by the 2000 background genes. So I reran it with the background seed changed to 1–3, and with no background at all (code and output in the collapsible block).

With seeds 1–3, gene A's gene-wise value is the same, 0.025388, but the trend varies over 0.128–0.138, the MAP over 0.093–0.098 and p over 0.0086–0.0106. The prior variance was the lower limit 0.25 all three times.

With only gene A's single row there is no information to borrow at all. A message says the trend curve (parametric fit) could not be fitted and was replaced by a local fit, and the trend value becomes equal to the gene's own gene-wise value. The MAP is then ≈ gene-wise (0.02538672), with p = 4.1 × 10⁻⁶, almost the same as the Cox-Reid row of the table in section 3 (p 4.12 × 10⁻⁶). The warning of section 15.3 of the study text shows up in numbers.

<details>
<summary>Code and output of the runs with different backgrounds</summary>

```r
# the same gene A with different background genes (seed) (size factors fixed at 1)
for (s in 1:3) {
  set.seed(s); bg <- makeExampleDESeqDataSet(n = 2000, m = 6, betaSD = 1)
  d <- DESeqDataSetFromMatrix(rbind(counts(bg), teach_gene = as.integer(k)), colData(bg)[, "condition", drop=FALSE], ~ condition)
  sizeFactors(d) <- rep(1, 6); d <- DESeq(d, quiet = TRUE); m <- mcols(d)["teach_gene", ]; rr <- results(d)["teach_gene", ]
  cat(sprintf("seed %d: dispGeneEst=%.6f dispFit=%.4f priorVar=%.3f dispMAP=%.4f lfcSE=%.4f p=%.5f\n",
              s, m$dispGeneEst, m$dispFit, attr(dispersionFunction(d), "dispPriorVar"), m$dispMAP, rr$lfcSE, rr$pvalue))
}
# gene A's single row in DESeq() with no background
d1 <- DESeqDataSetFromMatrix(matrix(as.integer(k), 1, dimnames=list("teach_gene", colnames(dds2))), colData(dds2), ~ condition)
sizeFactors(d1) <- rep(1, 6); d1 <- suppressWarnings(DESeq(d1, quiet = TRUE))
cat("fitType used:", attr(dispersionFunction(d1), "fitType"), "\n")
print(as.data.frame(mcols(d1)[, c("dispGeneEst","dispFit","dispMAP")])); print(as.data.frame(results(d1))[, c("lfcSE","stat","pvalue")])
```

```
seed 1: dispGeneEst=0.025388 dispFit=0.1333 priorVar=0.250 dispMAP=0.0958 lfcSE=0.3779 p=0.00970
seed 2: dispGeneEst=0.025388 dispFit=0.1280 priorVar=0.250 dispMAP=0.0926 lfcSE=0.3719 p=0.00859
seed 3: dispGeneEst=0.025388 dispFit=0.1375 priorVar=0.250 dispMAP=0.0984 lfcSE=0.3826 p=0.01063
-- note: fitType='parametric', but the dispersion trend was not well captured by the
   function: y = a/x + b, and a local regression fit was automatically substituted.
   specify fitType='local' or 'mean' to avoid this message next time.
fitType used: local 
           dispGeneEst    dispFit    dispMAP
teach_gene  0.02538756 0.02538756 0.02538672
               lfcSE     stat       pvalue
teach_gene 0.2122124 4.605197 4.120762e-06
```

</details>

## 5. What does the three-group paired practice code actually output?

Chapter 15 of the study text does a three-group paired analysis with the practice file `DESeq2_workbook.R`. In this series' example the conditions are Ctrl, Starvation and Starvation+Glucose. The third condition is called Glucose for short, but it is not ordinary medium with extra glucose. So the Glucose − Starvation comparison is read as "the effect of giving glucose again in the starved state", that is, in the reversal (rescue) direction. All three conditions came from the same donors, so the design is `~ pair + condition`, where pair stands for the donor.

But the workbook file is not in the repository, so I simulated data in the input format the study text specifies myself (the format is in the 15.1 collapsible block below). Since the true effects are known, the output can be scored too.

### The simulated data

4 donors × 3 conditions gives 12 samples, with 2000 genes. Counts were drawn from an NB, and each gene got a donor effect (individual differences) shared by the three conditions. 1600 genes are nulls with no condition effect, and the other 400 are of one of the five types below. $m$ is an effect size (log2) drawn between 1 and 2.5 for each gene.

| Type | Starvation − Ctrl | Glucose − Ctrl | Glucose − Starvation | Number of genes |
|---|---:|---:|---:|---:|
| rescue (comes all the way back) | −m | 0 | +m | 60 (+ mirror 60) |
| partial (comes halfway back) | −m | −m/2 | +m/2 | 40 (+ 40) |
| overshoot (goes up past Ctrl) | −m | +m | +2m | 40 (+ 40) |
| noRescue (does not come back) | −m | −m | 0 | 30 (+ 30) |
| GlucoseOnly (changes only in Glucose) | 0 | +m | +m | 30 (+ 30) |
| null | 0 | 0 | 0 | 1600 |

For each type, half had their signs flipped (`_mirror`) to balance the numbers of genes going up and going down. The direction written in the table (going down under Starvation) is the direction to be fixed in advance and tested in section 6.

<details>
<summary>The simulation model and the code that made the data (15.1)</summary>

The true model is this.

$$K_{ij} \sim \mathrm{NB}(\mu_{ij}, \alpha_i), \qquad \log_2 \mu_{ij} = \log_2 s_j + \beta_{i0} + \gamma_{i,\mathrm{pair}(j)} + \beta_{iS}\,[j \in \text{Starvation}] + \beta_{iG}\,[j \in \text{Glucose}], \qquad \alpha_i = 4/2^{\beta_{i0}} + 0.1$$

$\beta_{i0} \sim N(4, 2^2)$ and the formula for $\alpha_i$ are the same as the defaults of `makeExampleDESeqDataSet()`. $\gamma$ is the pair (individual) effect $N(0, 0.5^2)$ drawn per gene and pair, shared by the three conditions, so there is no interaction. That is, `~ pair + condition` has the same structure as the true model. The true size factors $s_j$ were drawn between 0.7 and 1.4, and the $\beta_S$, $\beta_G$ of the table above are the Starvation − Ctrl and Glucose − Ctrl columns.

`_mirror` genes are genes of the same type with only the signs flipped. Exactly half of each type was flipped to balance the numbers of up and down changes (`Starvation down/up: 170 170` in the output), because the median-ratio size factor relies on the premise that "most genes don't change, or the ups and downs balance". On a first try, with all effects only on the Starvation-down side and no mirrors (the numbers of types and the effect sizes were also different), the ratio of estimated to true size factor was 0.90–0.92 for the four Starvation samples and 1.07–1.15 for the other eight. Starvation alone was estimated about 17% too small relative to the rest (the output of this first try is not included in the note).

The condition of `makeExampleDESeqDataSet()` has only two levels (A/B). Simply renaming the 12 samples in order as Ctrl, Starvation and Glucose splits Starvation across A and B, creating a pair × condition interaction, and `~ pair + condition` becomes a misspecified design. An earlier version of this note did exactly that. So I simulated the three conditions and the pair effects directly from an NB. In the misspecified data of the earlier version, only the genes with large effects had inflated dispersions (1.318 at |trueBeta| > 1; the earlier output is not included in the note), but there is no such thing in this data (the section "If pairs are attached wrongly" in the body).

The input format the study text specifies is this. In `counts.csv` the first column is a unique gene ID and the other column names are sample IDs, and `metadata.csv` has the columns `sample_id`, `pair` and `condition`. The study text says that if this format differs from the real data, the design and contrasts must be revised to fit the research question, and states that the workbook, if run without input files, only makes synthetic data and is not a file that gives results for real research data. The code below writes the two files in this format.

```r
suppressPackageStartupMessages(library(DESeq2)); options(digits = 6, width = 110)
set.seed(2024)
n <- 2000
meta <- data.frame(sample_id = paste0("P", rep(1:4, 3), "_", rep(c("Ctrl","Starvation","Glucose"), each = 4)),
                   pair = rep(1:4, 3), condition = rep(c("Ctrl","Starvation","Glucose"), each = 4))
# true response types: log2 effects of Starvation (bS) and Glucose=Starvation+Glucose (bG) vs Ctrl. Size 1–2.5 per gene.
# exactly half of each type has its signs flipped (mirror): ups and downs must balance so the median-ratio size factor is not biased
type <- sample(rep(c("null","rescue","partial","overshoot","noRescue","GlucoseOnly"), c(1600, 120, 80, 80, 60, 60)))
sgn  <- ave(rep(1, n), type, FUN = function(x) sample(rep(c(-1, 1), length.out = length(x))))
size <- ifelse(type == "null", 0, sgn * runif(n, 1, 2.5))
bS <- size * c(null=0, rescue=-1, partial=-1,   overshoot=-1, noRescue=-1, GlucoseOnly=0)[type]
bG <- size * c(null=0, rescue= 0, partial=-0.5, overshoot= 1, noRescue=-1, GlucoseOnly=1)[type]
b0 <- rnorm(n, 4, 2)                        # log2 baseline expression (same distribution as the makeExampleDESeqDataSet default)
bP <- matrix(rnorm(n * 4, 0, 0.5), n, 4)    # per-gene pair (individual) effects: shared by the three conditions, no pair x condition interaction
trueDisp <- 4 / 2^b0 + 0.1                  # dispMeanRel of makeExampleDESeqDataSet
sf <- round(runif(12, 0.7, 1.4), 2)         # true size factors
log2mu <- b0 + bP[, meta$pair] + outer(bS, meta$condition == "Starvation") + outer(bG, meta$condition == "Glucose")
cts <- matrix(rnbinom(n * 12, mu = sweep(2^log2mu, 2, sf, `*`), size = 1 / trueDisp), n, 12,
              dimnames = list(paste0("gene", 1:n), meta$sample_id))
type <- ifelse(type != "null" & sgn < 0, paste0(type, "_mirror"), type)   # the table's direction (Starvation-down side) keeps the original type name
truth <- data.frame(type, bS, bG, trueDisp, row.names = rownames(cts))
write.csv(cts, "counts.csv"); write.csv(meta, "metadata.csv", row.names = FALSE)
cat(readLines("counts.csv", 2), readLines("metadata.csv", 3), sep = "\n")
print(table(truth$type)); cat("Starvation down/up:", sum(bS < 0), sum(bS > 0), "  Glucose down/up:", sum(bG < 0), sum(bG > 0), "\n")
```

```
"","P1_Ctrl","P2_Ctrl","P3_Ctrl","P4_Ctrl","P1_Starvation","P2_Starvation","P3_Starvation","P4_Starvation","P1_Glucose","P2_Glucose","P3_Glucose","P4_Glucose"
"gene1",0,0,0,0,0,0,0,0,0,0,0,0
"sample_id","pair","condition"
"P1_Ctrl",1,"Ctrl"
"P2_Ctrl",2,"Ctrl"

       GlucoseOnly GlucoseOnly_mirror           noRescue    noRescue_mirror               null
                30                 30                 30                 30               1600
         overshoot   overshoot_mirror            partial     partial_mirror             rescue
                40                 40                 40                 40                 60
     rescue_mirror
                60
Starvation down/up: 170 170   Glucose down/up: 140 140
```

</details>

### Aligning the inputs and checking the design (15.2)

The code of section 15.2 of the study text aligns the order of the count columns with the metadata rows and checks with `stopifnot` that the design matrix is full rank. Full rank means that the columns of the design matrix do not overlap in information, so all coefficients can be estimated. Run as it is, all three checks passed, but so did `meta_bad`, a mechanical labelling that shifts the pair numbers by one position for each condition (collapsible block). The final cross-table shows that each label mixes three different donors, yet it still passes.

Passing the checks only means the column order is right and X is full rank. Section 15.2 of the study text warns the same: "That it runs does not by itself verify that the pairing is true. The user must confirm that the pair IDs in the metadata mean a real biological correspondence. IDs given the same numbers mechanically do not create a statistically valid pairing." How costly such labels are is seen in numbers a little later.

<details>
<summary>The code of 15.2 of the study text and a test of mechanical pair labels</summary>

```r
cts <- as.matrix(read.csv("counts.csv", row.names=1, check.names=FALSE))
meta <- read.csv("metadata.csv", stringsAsFactors=FALSE); rownames(meta) <- meta$sample_id
stopifnot(setequal(colnames(cts), rownames(meta)))
cts <- cts[, rownames(meta), drop=FALSE]; stopifnot(identical(colnames(cts), rownames(meta)))
meta$pair <- factor(meta$pair); meta$condition <- factor(meta$condition, levels=c("Ctrl","Starvation","Glucose"))
X <- model.matrix(~ pair + condition, data=meta); stopifnot(qr(X)$rank == ncol(X), nrow(X) > ncol(X))
print(table(meta$pair, meta$condition)); print(colnames(X)); cat("rank:", qr(X)$rank, " n:", nrow(X), "\n")
# a case of the 15.2 warning: a 'mechanical' labelling that shifts pair numbers by one position per condition also passes the same checks
meta_bad <- meta; meta_bad$pair <- factor((as.integer(meta$pair) + as.integer(meta$condition) - 2) %% 4 + 1)
Xb <- model.matrix(~ pair + condition, data=meta_bad); stopifnot(qr(Xb)$rank == ncol(Xb), nrow(Xb) > ncol(Xb))
print(table(true_pair = meta$pair, label = meta_bad$pair))
```

```

    Ctrl Starvation Glucose
  1    1          1       1
  2    1          1       1
  3    1          1       1
  4    1          1       1
[1] "(Intercept)"         "pair2"               "pair3"               "pair4"
[5] "conditionStarvation" "conditionGlucose"
rank: 6  n: 12
         label
true_pair 1 2 3 4
        1 1 1 1 0
        2 0 1 1 1
        3 1 0 1 1
        4 1 1 0 1
```

</details>

### The basic analysis and the three contrasts (15.3)

A contrast is a vector that sets which two conditions to compare. With the study text's code I made the three comparisons (Starvation − Ctrl, Glucose − Starvation, Glucose − Ctrl) and checked the genes with padj < 0.05 against the true effects (code and output in the collapsible block below).

The pre-filter `rowSums(cts >= 10) >= 3` is a criterion the study text describes as "for teaching, not a fixed standard for every study", and 1574 of the 2000 genes remained. The residual degrees of freedom (number of samples − number of coefficients) are 12 − 6 = 6, and `betaPrior=FALSE` is the same as the default.

Among the genes with padj < 0.05, the proportion whose true LFC is 0 is 11/179 = 6.1%, 7/174 = 4.0% and 6/97 = 6.2%. What BH (the Benjamini–Hochberg adjustment) holds at 5% is the FDR, the expected proportion of false discoveries in the list of discoveries, so a single realization can be a little above 5% ([06](06_multiple_testing.md)). Of the genes with a true effect, the numbers detected are 168/286 for Starvation − Ctrl, 167/282 for Glucose − Starvation and 91/228 for Glucose − Ctrl.

Glucose − Starvation is not in `resultsNames`, but it is computed as the difference of two coefficients (the `G-S check` line of the output). The number of NA padj differs between contrasts because the threshold of independent filtering (the step that removes genes with a low mean count from the BH targets) is set separately from each contrast's p distribution (section 3 of [06](06_multiple_testing.md)).

<details>
<summary>The code of 15.3 of the study text, checked against the true effects</summary>

```r
# teaching pre-filter: not a fixed standard for every study (study text 15.3 comment)
keep <- rowSums(cts >= 10) >= 3; cts <- cts[keep, , drop=FALSE]; truth <- truth[keep, ]
cat("genes kept:", sum(keep), "of", length(keep), "\n")
dds0 <- DESeqDataSetFromMatrix(countData=cts, colData=meta, design=~pair+condition)
dds <- DESeq(dds0, betaPrior=FALSE, quiet=TRUE); print(resultsNames(dds))
res_S_C <- results(dds, contrast=c("condition","Starvation","Ctrl"), alpha=0.05)
res_G_S <- results(dds, contrast=c("condition","Glucose","Starvation"), alpha=0.05)
res_G_C <- results(dds, contrast=c("condition","Glucose","Ctrl"), alpha=0.05)
true_lfc <- list(res_S_C = truth$bS, res_G_S = truth$bG - truth$bS, res_G_C = truth$bG)
for (nm in names(true_lfc)) { r <- get(nm); hit <- which(r$padj < 0.05)
  cat(sprintf("%s: padj<0.05 = %d (true LFC = 0: %d)  true LFC != 0: %d  NA padj = %d\n",
              nm, length(hit), sum(true_lfc[[nm]][hit] == 0), sum(true_lfc[[nm]] != 0), sum(is.na(r$padj)))) }
g <- rownames(dds)[1]
cat("G-S check:", res_G_S[g,"log2FoldChange"], "==", res_G_C[g,"log2FoldChange"] - res_S_C[g,"log2FoldChange"], "\n")
```

```
genes kept: 1574 of 2000
[1] "Intercept"                    "pair_2_vs_1"                  "pair_3_vs_1"
[4] "pair_4_vs_1"                  "condition_Starvation_vs_Ctrl" "condition_Glucose_vs_Ctrl"
res_S_C: padj<0.05 = 179 (true LFC = 0: 11)  true LFC != 0: 286  NA padj = 61
res_G_S: padj<0.05 = 174 (true LFC = 0: 7)  true LFC != 0: 282  NA padj = 0
res_G_C: padj<0.05 = 97 (true LFC = 0: 6)  true LFC != 0: 228  NA padj = 0
G-S check: -1.02534 == -1.02534
```

</details>

<details>
<summary>Do the pair coefficients catch signal too?</summary>

```r
for (nm in c("pair_2_vs_1","pair_3_vs_1","pair_4_vs_1")) { r <- results(dds, name=nm)
  cat(nm, ": p<0.05 fraction =", round(mean(r$pvalue < 0.05, na.rm=TRUE), 3), " padj<0.05 =", sum(r$padj < 0.05, na.rm=TRUE), "\n") }
```

```
pair_2_vs_1 : p<0.05 fraction = 0.194  padj<0.05 = 100 
pair_3_vs_1 : p<0.05 fraction = 0.199  padj<0.05 = 89 
pair_4_vs_1 : p<0.05 fraction = 0.217  padj<0.05 = 104 
```

The pair coefficients are real signal in this data. The proportion with p < 0.05 is 0.19–0.22, and padj < 0.05 holds for 89–104 genes per coefficient. The donor effect $\gamma$ put into the simulation shows up as it is.

</details>

### If pairs are attached wrongly (the warning of 15.2 in numbers)

I analyse the same data with three designs and compare the ratio of the gene-wise α to the true α.

```r
# design comparison: true pair / mechanical pair shifted by one (meta_bad of 15.2) / pair dropped
eff <- truth$bS != 0 | truth$bG != 0
for (des in c("true pair", "shifted pair", "no pair")) {
  d <- if (des == "true pair") dds else DESeq(DESeqDataSetFromMatrix(cts, if (des == "no pair") meta else meta_bad,
         if (des == "no pair") ~ condition else ~ pair + condition), quiet = TRUE)
  ratio <- mcols(d)$dispGeneEst / truth$trueDisp
  r <- results(d, contrast = c("condition","Starvation","Ctrl"), alpha = 0.05); hit <- which(r$padj < 0.05)
  cat(sprintf("%-13s median dispGeneEst/trueDisp: null %.3f  effect %.3f | Starvation-Ctrl padj<0.05: %d (true LFC = 0: %d)\n",
              des, median(ratio[!eff]), median(ratio[eff]), length(hit), sum(truth$bS[hit] == 0))) }
cat("reference, median of chi-square/df: df=6", qchisq(0.5, 6)/6, " df=9", qchisq(0.5, 9)/9, "\n")
```

```
true pair     median dispGeneEst/trueDisp: null 0.852  effect 0.942 | Starvation-Ctrl padj<0.05: 179 (true LFC = 0: 11)
shifted pair  median dispGeneEst/trueDisp: null 1.373  effect 1.372 | Starvation-Ctrl padj<0.05: 118 (true LFC = 0: 2)
no pair       median dispGeneEst/trueDisp: null 1.307  effect 1.288 | Starvation-Ctrl padj<0.05: 129 (true LFC = 0: 2)
reference, median of chi-square/df: df=6 0.891353  df=9 0.926981 
```

In the true pair design the median ratio is 0.852 for null genes and 0.942 for effect genes. A variance estimate with 6 residual degrees of freedom naturally has a median a little below 1 (the median of χ²₆ / 6 = 0.891), and it is right at that level. α is not inflated for the effect genes alone either, so the design absorbed the condition effects properly.

Dropping pair (`~ condition`, 9 residual degrees of freedom, reference value 0.927) puts the donor differences into the residuals, so the ratio grows to 1.307 / 1.288 and the Starvation − Ctrl detections drop from 179 to 129. The mechanical pairs shifted by one position are worse still: they use 3 degrees of freedom yet fail to explain the donor differences, so the ratio is 1.373 and the detections 118. The study text's warning that "mechanical IDs do not create a valid pairing" shows up here as a result worse than no pairing at all.

### Tracing the intermediate values: which way did the MAP move? (15.4)

With the study text's `trace` table, I put the gene-wise value, trend, MAP and final α side by side for each gene.

```r
trace <- data.frame(gene=rownames(dds), baseMean=mcols(dds)$baseMean, alpha_gene=mcols(dds)$dispGeneEst,
  alpha_trend=mcols(dds)$dispFit, alpha_MAP=mcols(dds)$dispMAP, alpha_final=dispersions(dds),
  dispersion_outlier=mcols(dds)$dispOutlier, beta_converged=mcols(dds)$betaConv)
print(head(trace))
out <- trace$dispersion_outlier
cat("outliers:", sum(out), " final==MAP for non-outlier:", all(trace$alpha_final[!out] == trace$alpha_MAP[!out]),
    " final==geneEst for outlier:", all(trace$alpha_final[out] == trace$alpha_gene[out]), "\n")
lo <- !out & trace$alpha_gene < trace$alpha_trend
cat("non-outlier, gene-wise < trend:", sum(lo), " of which MAP > gene-wise:", sum(trace$alpha_MAP[lo] > trace$alpha_gene[lo]), "\n")
```

```
    gene baseMean alpha_gene alpha_trend alpha_MAP alpha_final dispersion_outlier beta_converged
1  gene4 14.53436  0.2870371    0.452593  0.388855    0.388855              FALSE           TRUE
2  gene7 13.46715  0.4588111    0.480609  0.472677    0.472677              FALSE           TRUE
3  gene8  7.49448  3.3130931    0.784685  1.352995    3.313093               TRUE           TRUE
4 gene10 44.58308  0.0647908    0.214312  0.153293    0.153293              FALSE           TRUE
5 gene11 23.21980  0.8211376    0.320352  0.480467    0.480467              FALSE           TRUE
6 gene12 32.31258  0.1276545    0.258079  0.207501    0.207501              FALSE           TRUE
outliers: 5  final==MAP for non-outlier: TRUE  final==geneEst for outlier: TRUE 
non-outlier, gene-wise < trend: 958  of which MAP > gene-wise: 958 
```

`head(trace)` starts at gene4 because gene1–3 and others were removed by the pre-filter. gene4, gene10 and gene12 had gene-wise values below the trend, and their MAP got larger (0.287 → 0.389, 0.0648 → 0.153, 0.128 → 0.208). All 958 non-outliers with gene-wise < trend did this, and appendix Problem 9 is exactly this situation. gene11, conversely, had a gene-wise 0.821 above the trend 0.320, so it went down to 0.480.

gene8 is `dispOutlier = TRUE`. A gene flagged as an outlier because its gene-wise value is too far above the trend discards the MAP (1.353) and uses the gene-wise value (3.313) as its final α. There were 5 such genes (Problem 10). The decision rule is in [03](03_dispersion_estimation.md).

<details>
<summary>Size factors, the mu assay, column descriptions (the rest of the output of 15.4)</summary>

The study text's `sizeFactors(dds)` line was changed to compare with the true values. `ratio_rescaled` is the estimated/true ratio scaled to a geometric mean of 1.

```r
r_sf <- sizeFactors(dds) / sf; print(round(rbind(estimated = sizeFactors(dds), true = sf, ratio_rescaled = r_sf / exp(mean(log(r_sf)))), 3))
print(head(assays(dds)[["mu"]][, 1:6], 3))
print(mcols(mcols(dds))[c("baseMean","dispGeneEst","dispFit","dispMAP","dispersion","dispOutlier","betaConv"), ])
cat("betaConv all TRUE:", all(trace$beta_converged), "  assayNames(dds):", assayNames(dds), "\n")
```

```
               P1_Ctrl P2_Ctrl P3_Ctrl P4_Ctrl P1_Starvation P2_Starvation P3_Starvation P4_Starvation
estimated        1.297   1.199   0.756   1.324          0.69         1.137         1.082         0.797
true             1.360   1.280   0.780   1.400          0.72         1.160         1.100         0.810
ratio_rescaled   0.986   0.968   1.002   0.978          0.99         1.013         1.016         1.016
               P1_Glucose P2_Glucose P3_Glucose P4_Glucose
estimated           1.165      0.912      1.076      1.319
true                1.220      0.940      1.060      1.370
ratio_rescaled      0.987      1.002      1.049      0.995
      P1_Ctrl  P2_Ctrl  P3_Ctrl P4_Ctrl P1_Starvation P2_Starvation
gene4 48.3873  9.06572  9.26289 18.2875      27.13488       9.06367
gene7 18.9198 15.83541 16.53846 35.6956       4.02686       6.00876
gene8 49.5931  6.16720 10.06323 14.2651       6.37811       1.41405
DataFrame with 7 rows and 2 columns
                    type            description
             <character>            <character>
baseMean    intermediate mean of normalized c..
dispGeneEst intermediate gene-wise estimates ..
dispFit     intermediate fitted values of dis..
dispMAP     intermediate maximum a posteriori..
dispersion  intermediate final estimate of di..
dispOutlier intermediate dispersion flagged a..
betaConv         results   convergence of betas
betaConv all TRUE: TRUE   assayNames(dds): counts mu H cooks
```

Apart from a constant multiple, the size factors recover the true values within ±5% (`ratio_rescaled` 0.968–1.049).

</details>

### Running it split into step functions (15.5)

`DESeq()` calls the size factor, gene-wise dispersion, trend, MAP and Wald test functions in turn. Section 15.5 of the study text calls these functions one at a time, and while doing so I checked whether the mu assay changes at each step (code and output in the collapsible block below). The mu assay is the matrix holding, for each sample, the expected count μ_ij the model predicts.

mu is first created by `estimateDispersionsGeneEst()` when it fits with the initial α, and the `Fit` and `MAP` steps leave it alone. Then `nbinomWaldTest()` overwrites it when it refits with the final α. The median relative difference between the first and the last is 0.46%, and the absolute difference is at most 49.3. This is what section 15.4 of the study text means by "the final mu assay can be overwritten at a later step".

The biggest change is P1_Glucose of gene1559. As α changed from 0.381 to 0.241, mu went from 651.8 to 602.5. The size factors of the two steps are the same, so if mu = $s_j e^{Xb}$ changed, the coefficients $b$ (and hence β) moved with α. A design with a pair term does not solve as group means, so α affects the coefficients through the weights. In the chapter 16 example the means did not change only because the design was simple.

The step-by-step result `d4` matched `DESeq()` down to the dispersions and p-values. Section 15.5 of the study text writes that "listing the low-level calls does not give exactly the same result as the top-level function in every setting", which is right in general. `DESeq()` has an additional step that replaces count outliers found with Cook's distance (a measure of how much one count shakes the coefficient estimates) and refits.

But this step is considered only when one design cell has 7 or more samples. A design cell is a set of samples whose design matrix rows are completely identical. In `~ pair + condition` all 12 rows differ, so the cell size is 1, and no refit happens however many replicates there are per condition. `DESeq()` refits only when the same cell has 7 or more samples and there is a Cook's outlier (run examples: the m = 14 experiment of [00](00_overview.md), experiment 3 of [07](07_lfc_shrinkage_and_qc.md)).

<details>
<summary>The code of 15.5 of the study text, and the assays and columns each step creates</summary>

```r
d0 <- estimateSizeFactors(dds0); d1 <- estimateDispersionsGeneEst(d0); mu_initial <- assays(d1)[["mu"]]
d2 <- estimateDispersionsFit(d1); d3 <- estimateDispersionsMAP(d2)
d4 <- nbinomWaldTest(d3, betaPrior=FALSE); mu_final <- assays(d4)[["mu"]]
cat("mu_initial == mu_final ?", isTRUE(all.equal(mu_initial, mu_final)), " max|diff| =", max(abs(mu_initial - mu_final)),
    " median rel diff =", median(abs(mu_initial - mu_final)/mu_final), "\n")
ij <- which(abs(mu_initial - mu_final) == max(abs(mu_initial - mu_final)), arr.ind = TRUE); g <- rownames(dds)[ij[1, 1]]
cat("largest mu change:", g, colnames(dds)[ij[1, 2]], ": mu_initial =", mu_initial[ij], " mu_final =", mu_final[ij],
    " | alpha_gene =", mcols(d1)[g, "dispGeneEst"], " alpha_final =", dispersions(d4)[ij[1, 1]], "\n")
cat("d4 vs DESeq(): dispersions equal?", isTRUE(all.equal(dispersions(d4), dispersions(dds))),
    " pvalues equal?", isTRUE(all.equal(results(d4)$pvalue, results(dds)$pvalue)), "\n")
mm <- attr(dds, "modelMatrix")
cat("design cells (identical model-matrix rows): max size =", max(table(apply(mm, 1, paste, collapse = "_"))),
    "  any(nOrMoreInCell(mm, 7)) =", any(DESeq2:::nOrMoreInCell(mm, 7)), "\n")
```

```
mu_initial == mu_final ? FALSE  max|diff| = 49.3356  median rel diff = 0.00458844 
largest mu change: gene1559 P1_Glucose : mu_initial = 651.792  mu_final = 602.457  | alpha_gene = 0.38132  alpha_final = 0.240576 
d4 vs DESeq(): dispersions equal? TRUE  pvalues equal? TRUE 
design cells (identical model-matrix rows): max size = 1   any(nOrMoreInCell(mm, 7)) = FALSE 
```

```r
cat("assayNames d0:", assayNames(d0), "| d1:", assayNames(d1), "| d3:", assayNames(d3), "| d4:", assayNames(d4), "\n")
cat("d2 mu identical to d1?", identical(assays(d2)[["mu"]], mu_initial), " d3 mu identical to d1?", identical(assays(d3)[["mu"]], mu_initial), "\n")
cat("mcols(d1) cols:", colnames(mcols(d1)), "\n"); cat("mcols(d3) new cols:", setdiff(colnames(mcols(d3)), colnames(mcols(d1))), "\n")
cat("mcols(d4) new cols:", setdiff(colnames(mcols(d4)), colnames(mcols(d3))), "\n")
```

```
assayNames d0: counts | d1: counts mu | d3: counts mu | d4: counts mu H cooks 
d2 mu identical to d1? TRUE  d3 mu identical to d1? TRUE 
mcols(d1) cols: baseMean baseVar allZero dispGeneEst dispGeneIter 
mcols(d3) new cols: dispFit dispersion dispIter dispOutlier dispMAP 
mcols(d4) new cols: Intercept pair_2_vs_1 pair_3_vs_1 pair_4_vs_1 condition_Starvation_vs_Ctrl condition_Glucose_vs_Ctrl SE_Intercept SE_pair_2_vs_1 SE_pair_3_vs_1 SE_pair_4_vs_1 SE_condition_Starvation_vs_Ctrl SE_condition_Glucose_vs_Ctrl WaldStatistic_Intercept WaldStatistic_pair_2_vs_1 WaldStatistic_pair_3_vs_1 WaldStatistic_pair_4_vs_1 WaldStatistic_condition_Starvation_vs_Ctrl WaldStatistic_condition_Glucose_vs_Ctrl WaldPvalue_Intercept WaldPvalue_pair_2_vs_1 WaldPvalue_pair_3_vs_1 WaldPvalue_pair_4_vs_1 WaldPvalue_condition_Starvation_vs_Ctrl WaldPvalue_condition_Glucose_vs_Ctrl betaConv betaIter deviance maxCooks 
```

The mu of `d2` and `d3` is identical to that of `d1`. In `mcols`, `estimateDispersionsGeneEst()` adds `dispGeneEst` and others, `estimateDispersionsFit()` and `estimateDispersionsMAP()` add `dispFit`, `dispMAP`, `dispOutlier` and the final `dispersion`, and `nbinomWaldTest()` adds the coefficients, SEs, Wald statistics, p-values, `betaConv` and `maxCooks`.

</details>

## 6. How is the two-direction claim "it goes down and comes back" tested?

### To claim two conditions at once (15.6)

The claim "starvation lowers the expression and glucose brings it back" is really two claims combined. Starvation − Ctrl < 0 (it goes down under starvation) must hold, and at the same time Glucose − Starvation > 0 (it goes up when glucose is added).

Each claim has a one-sided p that looks in one direction only. `altHypothesis="less"` of `results()` tests the LFC < 0 side and `"greater"` the LFC > 0 side; with Z = LFC / SE, the p of less is Φ(Z) and the p of greater 1 − Φ(Z).

For example, as in appendix Problem 20, suppose the one-sided p of Starvation-down is 0.003 and the one-sided p of Glucose-up is 0.08. To say both conditions are true, both must pass, so the weaker one, 0.08, becomes the criterion. Using the larger of the two p-values like this is called the intersection–union test.

$$p_{\text{joint}} = \max(p_1, p_2)$$

$p_1$ is the one-sided p of Starvation-down, and $p_2$ the one-sided p of Glucose-up. After computing $p_{\text{joint}}$ for each gene, adjust with BH.

The study text pins the scope of this test down in two ways. One is that it tests directions fixed in advance (Starvation-down, Glucose-up); if the directions are chosen after looking at the data, this p loses its meaning. The other is that it "is not a test that automatically establishes complete normalization or mechanistic rescue". After the study text's code I appended a tally by true type.

```r
r_down <- results(dds, contrast=c("condition","Starvation","Ctrl"), altHypothesis="less", independentFiltering=FALSE)
r_up   <- results(dds, contrast=c("condition","Glucose","Starvation"), altHypothesis="greater", independentFiltering=FALSE)
p1 <- r_down$pvalue; p2 <- r_up$pvalue; p1[is.na(p1)] <- 1; p2[is.na(p2)] <- 1
p_joint <- pmax(p1, p2); padj_joint <- p.adjust(p_joint, method="BH")
cat("n p1<0.05:", sum(p1<0.05), " n p2<0.05:", sum(p2<0.05), " n p_joint<0.05:", sum(p_joint<0.05), " n padj_joint<0.05:", sum(padj_joint<0.05), "\n")
# passes by true response type: proper (BH after pmax) vs wrong (BH after pmin)
lv <- c(outer(c("rescue","partial","overshoot","noRescue","GlucoseOnly"), c("", "_mirror"), paste0), "null")
hit_max <- padj_joint < 0.05; hit_min <- p.adjust(pmin(p1, p2), "BH") < 0.05
print(cbind(all = table(factor(truth$type, lv)), joint_pmax = table(factor(truth$type[hit_max], lv)),
            wrong_pmin = table(factor(truth$type[hit_min], lv))))
```

```
n p1<0.05: 188  n p2<0.05: 176  n p_joint<0.05: 94  n padj_joint<0.05: 41
                    all joint_pmax wrong_pmin
rescue               52         24         42
partial              27          2         14
overshoot            36         15         32
noRescue             20          0         10
GlucoseOnly          27          0         15
rescue_mirror        50          0          0
partial_mirror       36          0          0
overshoot_mirror     37          0          0
noRescue_mirror      28          0          0
GlucoseOnly_mirror   17          0          0
null               1244          0          9
```

With the proper method (BH after pmax), 41 genes have padj_joint < 0.05: 24 rescue, 2 partial and 15 overshoot, all types for which "Starvation-down and Glucose-up" is true, and not a single noRescue, GlucoseOnly or null. The `_mirror` types also give 0; for example, rescue_mirror genes go up under starvation and come back with glucose, which is not the direction fixed in advance, so it is right that they are not caught.

The wrong method of choosing the smaller p (BH after pmin), on the other hand, passes 122. Of these, 34 (10 noRescue, 15 GlucoseOnly, 9 null) have only one of the two conditions true, or neither. The answer to appendix Problem 20 comes out as actual counts.

### Why the two-sided p must not be halved

I also ran what happens if the two-sided p is halved instead of using the one-sided p (code and output in the collapsible block below). 764 genes have a positive Starvation − Ctrl LFC, that is, actually went up, and 189 of them fall below 0.05 when the two-sided p is halved, looking like "significant Starvation-down". With the proper `altHypothesis="less"` there are 0. Exactly the warning of section 15.6 of the study text: "simply halving the two-sided p-value can wrongly make even observations of the opposite sign significant, so using the proper one-sided call is clearer."

<details>
<summary>Halving the two-sided p, what the one-sided p is, and the trap in the stat column</summary>

```r
# halving the two-sided p makes even genes of the opposite sign (LFC>0) 'Starvation-down significant'
r2 <- results(dds, contrast=c("condition","Starvation","Ctrl"), independentFiltering=FALSE)
ok <- !is.na(r2$pvalue); up <- ok & r2$log2FoldChange > 0
cat("LFC>0 genes:", sum(up), " p_two/2 < 0.05 among them:", sum(r2$pvalue[up]/2 < 0.05), " p_less < 0.05 among them:", sum(r_down$pvalue[up] < 0.05), "\n")
```

```
LFC>0 genes: 764  p_two/2 < 0.05 among them: 189  p_less < 0.05 among them: 0 
```

```r
r_gr <- results(dds, contrast=c("condition","Starvation","Ctrl"), altHypothesis="greater", independentFiltering=FALSE)
z <- (r2$log2FoldChange / r2$lfcSE)[ok]
cat("less p == pnorm(LFC/SE):", isTRUE(all.equal(r_down$pvalue[ok], pnorm(z))), "  less stat == pmin(LFC/SE, 0):", isTRUE(all.equal(r_down$stat[ok], pmin(z, 0))), "\n")
cat("less + greater == 1:", isTRUE(all.equal(r_down$pvalue[ok] + r_gr$pvalue[ok], rep(1, sum(ok)))),
    "  two-sided == 2*min(less, greater):", isTRUE(all.equal(r2$pvalue[ok], 2*pmin(r_down$pvalue[ok], r_gr$pvalue[ok]))), "\n")
```

```
less p == pnorm(LFC/SE): TRUE   less stat == pmin(LFC/SE, 0): TRUE 
less + greater == 1: TRUE   two-sided == 2*min(less, greater): TRUE 
```

The p of `altHypothesis="less"` is the lower tail `pnorm(LFC/SE)` (threshold 0). The one-sided p-values of the two directions add up to 1, and the two-sided p is twice the smaller one. These relations follow directly from the `results()` source ("What I checked in the DESeq2 source" below). But the `stat` column of the results table is clipped to `pmin(LFC/SE, 0)`, so for genes with a positive LFC, `pnorm(stat)` = 0.5. So rebuilding p from `stat` goes wrong. The p itself is the unclipped `pnorm(LFC/SE)`.

</details>

### If both directions are significant, is it normalization? (Problem 19)

Among the 41 that passed the joint test, this time I look at Glucose − Ctrl (two-sided), because Glucose − Ctrl being 0 would mean the expression came all the way back to the Ctrl level.

```r
# Problem 19: what does Glucose-Ctrl (two-sided) look like in the genes that passed the joint test?
r_gc <- results(dds, contrast=c("condition","Glucose","Ctrl"), independentFiltering=FALSE)
print(table(type = factor(truth$type[hit_max], lv[1:3]), Glucose_vs_Ctrl_p_lt_0.05 = r_gc$pvalue[hit_max] < 0.05))
```

```
           Glucose_vs_Ctrl_p_lt_0.05
type        FALSE TRUE
  rescue       23    1
  partial       2    0
  overshoot     3   12
```

Of the 41, 28 have a Glucose − Ctrl p ≥ 0.05 (23 rescue, 2 partial, 3 overshoot). Reading these 28 as "no difference, so normalized" would misread 5 of them, the 2 partial and 3 overshoot genes whose true Glucose − Ctrl is not 0, because p ≥ 0.05 is not evidence of no difference. Meanwhile, 12 of the 15 overshoot genes are significant for Glucose − Ctrl too, showing that they went up past the Ctrl level, and 1 of the 24 rescue genes had p < 0.05 even though its true Glucose − Ctrl = 0 (a type 1 error).

So, as in the study text's answer, claiming normalization needs a margin set in advance and an equivalence test (a test that shows the difference lies within the margin).

<details>
<summary>The three genes with the smallest joint p</summary>

```r
i <- order(p_joint)[1:3]
print(data.frame(gene=rownames(dds)[i], type=truth$type[i], lfc_S_C=round(r_down$log2FoldChange[i],3), p_down=signif(p1[i],3),
                 lfc_G_S=round(r_up$log2FoldChange[i],3), p_up=signif(p2[i],3), p_joint=signif(p_joint[i],3), padj_joint=signif(padj_joint[i],3)))
```

```
      gene      type lfc_S_C   p_down lfc_G_S     p_up  p_joint padj_joint
1  gene385 overshoot  -2.886 7.70e-12   4.096 3.78e-22 7.70e-12   1.21e-08
2 gene1309    rescue  -2.144 4.22e-11   2.154 3.36e-11 4.22e-11   3.32e-08
3  gene535    rescue  -2.528 9.50e-09   2.674 1.34e-09 9.50e-09   4.98e-06
```

The top 3 are 1 overshoot and 2 rescue. The joint p is the larger of the two one-sided p-values (gene385: max(7.70e-12, 3.78e-22) = 7.70e-12).

</details>

### What must be saved so it can be rebuilt? (15.7)

The code of section 15.7 of the study text saves four files: `dds` (RDS), the Starvation − Ctrl results table (CSV), the `trace` table (CSV) and `sessionInfo()`. I ran it as it is and opened the saved files to see what is left (code and output in the collapsible block). Checking the four files against the list section 15.7 of the study text says to preserve:

| Item to preserve in 15.7 of the study text | Are the four files enough? |
|---|---|
| Original counts and metadata | Only partly. `dds_wald.rds` has only the 1574 rows after the pre-filter, so the original `counts.csv` and `metadata.csv` must be kept separately |
| Pre-filter criterion | Missing. `rowSums(cts >= 10) >= 3` is only in the code |
| Design matrix | Present. The rds stores the design formula `~pair + condition` and the `modelMatrix` attribute |
| Direction of the contrast | Not in the CSV. Its columns are just `baseMean ... padj`, and the "condition Starvation vs Ctrl" description is only in the R object `mcols(res_S_C)`, so the file name is the only clue |
| Original p and padj | Only one contrast (Starvation vs Ctrl) is in a CSV |
| Effect estimator | Not stated in the CSV. This analysis is the MLE (`betaPrior=FALSE`, no shrinkage) |
| Package versions and session information | Present. `sessionInfo.txt` and `metadata(dds)$version` (1.50.2) |
| Annotation and upstream quantification versions | Missing |

It is just as section 15.7 of the study text says: "Saving dds as an RDS lets you check the intermediate estimates, but that alone does not account for every source of the raw data." In this simulation the true effect table `truth` was not saved to a file either, so it can only be obtained by rerunning the 15.1 block (seed 2024).

<details>
<summary>The code of 15.7 of the study text and a check of the saved files</summary>

```r
saveRDS(dds, "dds_wald.rds"); write.csv(as.data.frame(res_S_C), "Starvation_vs_Ctrl.csv")
write.csv(trace, "dispersion_trace.csv", row.names=FALSE); writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
print(file.info(c("dds_wald.rds","Starvation_vs_Ctrl.csv","dispersion_trace.csv","sessionInfo.txt"))[, "size", drop=FALSE])
cat(grep("DESeq2|^R version", readLines("sessionInfo.txt"), value=TRUE), sep="\n")
d <- readRDS("dds_wald.rds")
cat("rds: nrow =", nrow(d), " design =", deparse(design(d)), " modelMatrix stored:", !is.null(attr(d, "modelMatrix")),
    " metadata:", names(metadata(d)), as.character(metadata(d)$version), "\n")
cat("CSV columns:", names(read.csv("Starvation_vs_Ctrl.csv"))[-1], "\n")
cat("contrast (in R object only):", mcols(res_S_C)$description[2], "\n")
```

```
                         size
dds_wald.rds           977187
Starvation_vs_Ctrl.csv 184857
dispersion_trace.csv   172535
sessionInfo.txt          1820
R version 4.5.2 (2025-10-31)
 [1] DESeq2_1.50.2               SummarizedExperiment_1.40.0 Biobase_2.70.0
rds: nrow = 1574  design = ~pair + condition  modelMatrix stored: TRUE  metadata: version 1.50.2
CSV columns: baseMean log2FoldChange lfcSE stat pvalue padj
contrast (in R object only): log2 fold change (MLE): condition Starvation vs Ctrl
```

The file sizes (e.g. 1820 bytes for `sessionInfo.txt`) differ a little between environments. The table above is based on this output.

</details>

## Summary

- Gene A's α moved 0.0148 → 0.0254 → 0.0531 through $\ell$ → $\ell_{CR}$ → MAP, and computing the weights and SE with this α gave SE 0.289, Z 3.38 and p 0.00072. With two groups and size factors of 1, changing α leaves the LFC unchanged and changes only the SE and p, and assuming a Poisson makes p about 10¹⁹ times smaller.
- Actual DESeq2 does the same operations, so fixing α at the study text's value makes lfcSE, stat and p exactly equal to the hand calculation. What differs is the prior. The trend is estimated from the background genes (0.133) and the prior variance hits the lower limit 0.25, so the actual MAP is 0.0958, varying over 0.093–0.098 with the background.
- Dispersion shrinkage pulls toward the trend, so values below the trend go up and values above it go down. Genes flagged as outliers use their gene-wise value as it is.
- In a design with pairs, α moves the means and β too. Mechanically attached pair labels pass the checks but gave results worse even than dropping pair altogether.
- A two-direction claim must be adjusted with BH on the larger of two one-sided p-values. Choosing the smaller one or halving the two-sided p lets wrong genes in, and both directions being significant does not mean normalization.
- Even if `dds` is saved as an RDS, the original counts, the pre-filter rule and the direction of the contrast must be kept separately.

## Exercises

No appendix problems are assigned to this note specifically. Chapter 16 is itself a comprehensive example, so I chose six problems for which the calculations of this note serve as real examples. Concept-centred solutions are in [03](03_dispersion_estimation.md), [04](04_glm_condition_batch.md), [05](05_wald_vs_lrt.md) and [06](06_multiple_testing.md).

**Problem 9.** With gene-wise α=0.02 and trend α=0.1, the final α came out at 0.05. Since this is shrinkage yet the value grew, is it an error?

<details>
<summary>Solution</summary>

It is not an error. The target of dispersion shrinkage is not 0 but the trend, so gene-wise values below the trend are pulled up. 0.05 lies between 0.02 and 0.1. Gene A likewise had its gene-wise 0.0254 pulled toward the prior centre 0.08, becoming 0.0531 (section 2). In the simulation of section 5 too, all 958 non-outliers with gene-wise < trend got larger at the MAP.

</details>

**Problem 10.** Must `dispersions(dds)` agree with `mcols(dds)$dispMAP` for every gene?

<details>
<summary>Solution</summary>

No. Genes with `dispOutlier` TRUE keep their gene-wise value as the final α. gene8 of section 5 has a `dispMAP` of 1.353 but a final α equal to its gene-wise value, 3.313. 5 of the 1574 genes were like this, and for all the rest the final α = the MAP.

</details>

**Problem 11.** At the same expected count μ=100, α rose from 0.01 to 0.1. Compute the information weight for each and explain the direction of the effect on the SE.

<details>
<summary>Solution</summary>

At α = 0.01 it is 100 / (1 + 100 × 0.01) = 50, and at α = 0.1 it is 100 / (1 + 10) = 9.09. With the means and design unchanged, the information falls, so the covariance and SE grow. For gene A too, changing α from the MLE 0.0148 to the MAP 0.0531 lowered the Ctrl weight from 41.39 to 15.99 and raised the SE (log2) from 0.174 to 0.289 (the α comparison table in section 3).

</details>

**Problem 13.** log2FC=1.0, SE=0.25. Compute the default Wald Z and the nominal 95% CI. If the CI excludes 0, must BH padj<0.05 also hold?

<details>
<summary>Solution</summary>

Z = 1.0 / 0.25 = 4, and CI = 1.0 ± 1.96 × 0.25 = [0.51, 1.49]. But even if the CI excludes 0, padj < 0.05 is not guaranteed. The CI and Z correspond only to that one gene's raw p, while the BH padj depends on the ranks of the p-values of the other genes adjusted together. Section 4 (B) is a real example: gene A has p = 0.0097 and a nominal CI of 0.977 ± 1.96 × 0.378 = [0.237, 1.718], which excludes 0, but its BH padj over 2001 genes is 0.061.

</details>

**Problem 19.** Starvation−Ctrl<0 and Glucose−Starvation>0 are both significant, and Glucose−Ctrl has p=0.4. Is this enough to claim complete normalization?

<details>
<summary>Solution</summary>

It is not enough. Significance in both directions also arises from overshoot, and p = 0.4 is not evidence of no difference. In the simulation of section 6, of the 41 that passed the joint test, 28 had a Glucose − Ctrl p ≥ 0.05, and among them were 2 partial and 3 overshoot genes whose true Glucose − Ctrl is not 0. Claiming normalization needs a margin set in advance, an equivalence test to match it, and a multiple-testing adjustment fitting the joint claim.

</details>

**Problem 20.** The one-sided p-values of Starvation-down and Glucose-up are 0.003 and 0.08. What is the intersection–union joint p-value that both conditions hold at the same time? Why must the smaller one not be chosen?

<details>
<summary>Solution</summary>

The joint p is max(0.003, 0.08) = 0.08, because claiming "at the same time" requires passing even the weaker hypothesis. Choosing the smaller p would misread genes for which only one of the two holds as evidence for both conditions. In section 6, pmax passed 41 genes, all of types for which both conditions are true, while pmin passed 122, of which 34 had only one condition true or were null.

</details>

## Going deeper

<details>
<summary>What I checked in the DESeq2 source</summary>

The excerpts below were obtained from the installed DESeq2 1.50.2 by running `grep` on `deparse(body(...))` or printing whole functions with `deparse()`. They were run in order in one R session separate from the body. For parts whose full text was already seen in notes 01–07, I give only a summary and a link, and excerpt only the lines used directly in this note's calculations.

#### The betaPrior default of DESeq() and the refit condition

`betaPrior` is FALSE when the argument is not given. So `DESeq(dds0, betaPrior=FALSE)` in section 15.3 of the study text is the same as the default behaviour. Outlier replacement and refitting are considered only when at least one design cell (a set of samples whose model matrix rows are completely identical) has `minReplicatesForReplace = 7` or more. The criterion is the cell size, not the number of replicates per condition. With `fitType="glmGamPoi"`, `minReplicatesForReplace` changes to `Inf` and refitting is turned off.

```r
suppressPackageStartupMessages(library(DESeq2))
b <- deparse(body(DESeq))
cat(grep("minReplicatesForReplace <- Inf|missing\\(betaPrior\\)|betaPrior <- FALSE$|sufficientReps <-", b, value = TRUE), sep = "\n")
cat(deparse(DESeq2:::nOrMoreInCell), sep = "\n")
```

```
        minReplicatesForReplace <- Inf
    if (missing(betaPrior)) {
        betaPrior <- FALSE
    sufficientReps <- any(nOrMoreInCell(attr(object, "modelMatrix"), 
function (modelMatrix, n) 
{
    row_hash <- apply(modelMatrix, 1, paste0, collapse = "_")
    hash_table <- table(row_hash)
    numEqual <- as.vector(unname(hash_table[row_hash]))
    numEqual >= n
}
```

`nOrMoreInCell` groups each row of the model matrix as a string and counts how many times the same row appears. In a design with a pair term, every row has a different pair, so the cell size is 1 (see 15.5 in section 5). When a refit actually happens, and the whole flow, are in the collapsible blocks "Checked in the DESeq2 source (3): when do the refit and outlier replacement happen?" of [00](00_overview.md) and "The replacement conditions of DESeq() and the definition of maxCooks (source)" of [07](07_lfc_shrinkage_and_qc.md).

#### Gene-wise estimation: niter = 1, Cox-Reid on by default, and for group designs mu is "group mean of the normalized counts × s_j"

If `linearMu` is NULL, it becomes TRUE when "the number of unique rows (groups) of the model matrix == the number of columns" (FALSE with weights). Then, instead of a GLM, `linearModelMuNormalized` fits a linear model to the normalized counts and multiplies the size factors back in. For a group design this is the group mean of the normalized counts × $s_j$. In the chapter 16 example $s_j = 1$, so mu is just the arithmetic group mean.

```r
args(estimateDispersionsGeneEst)
g <- deparse(body(estimateDispersionsGeneEst))
cat(grep("linearMu <-|if \\(!linearMu\\)|linearModelMuNormalized", g, value = TRUE), sep = "\n")
cat(deparse(DESeq2:::linearModelMuNormalized), sep = "\n")
```

```
function (object, minDisp = 1e-08, kappa_0 = 1, dispTol = 1e-06, 
    maxit = 100, useCR = TRUE, weightThreshold = 0.01, quiet = FALSE, 
    modelMatrix = NULL, niter = 1, linearMu = NULL, minmu = if (type == 
        "glmGamPoi") 1e-06 else 0.5, alphaInit = NULL, type = c("DESeq2", 
        "glmGamPoi")) 
NULL
        linearMu <- nlevels(modelMatrixGroups) == ncol(modelMatrix)
            linearMu <- FALSE
        if (!linearMu) {
            fitMu <- linearModelMuNormalized(objectNZ[fitidx, 
function (object, x) 
{
    cts <- counts(object)
    norm.cts <- counts(object, normalized = TRUE)
    muhat <- linearModelMu(norm.cts, x)
    nf <- getSizeOrNormFactors(object)
    muhat * nf
}
```

#### The MAP and the outlier rule

The excerpt and the three branches are in the collapsible block "Reading the DESeq2 source (4): the prior width, the MAP and the outlier rule" of [03](03_dispersion_estimation.md). This note uses three facts: the prior mean is `log(dispFit)`, the prior variance is `pmax(varLogDispEsts - trigamma((m - p)/2), 0.25)` when the residual degrees of freedom $n-p>3$ (written `m - p` in the source), and genes flagged `dispOutlier` keep their gene-wise value as the final dispersion.

#### Wald: fitted in natural log and converted to log2, SE based on (XᵀWX + λ)⁻¹, λ = 1e-6

With `betaPrior=FALSE`, `betaPriorVar = 1e6`, that is, ridge λ = 1e-6 (on the log2 scale). The coefficients and SEs are multiplied by `log2(exp(1))` to convert to the log2 scale, and p is `2 * pnorm(abs(WaldStatistic), lower.tail = FALSE)` (`pt` only when `useT=TRUE`). The whole flow is in the collapsible block "Checked in the DESeq2 source: the Wald test (nbinomWaldTest)" of [05](05_wald_vs_lrt.md).

```r
f <- deparse(body(DESeq2:::fitNbinomGLMs))
cat(grep("lambda <- rep|betaMatrix <- log2|betaSE <- log2", f, value = TRUE), sep = "\n")
w <- deparse(body(nbinomWaldTest))
cat(grep("betaPriorVar <- rep\\(1e\\+06|WaldStatistic <-|WaldPvalue <-", w, value = TRUE), sep = "\n")
```

```
        lambda <- rep(1e-06, ncol(modelMatrix))
    betaMatrix <- log2(exp(1)) * betaRes$beta_mat
    betaSE <- log2(exp(1)) * sqrt(pmax(betaRes$beta_var_mat, 
        betaPriorVar <- rep(1e+06, ncol(fit$modelMatrix))
    WaldStatistic <- betaMatrix/betaSE
        WaldPvalue <- 2 * pt(abs(WaldStatistic), df = df, lower.tail = FALSE)
        WaldPvalue <- 2 * pnorm(abs(WaldStatistic), lower.tail = FALSE)
```

#### The one-sided p of results()

`altHypothesis="less"` computes p as the lower tail `pnorm(LFC/SE)` (threshold 0), but puts the `stat` column in clipped to `pmin(LFC/SE, 0)`. So reproducing p from the `stat` of the results table goes wrong for genes with a positive LFC (seen in numbers in the collapsible block of section 6). The same branch is also in the collapsible block "The p formula of lfcThreshold and the six altHypothesis options" of [05](05_wald_vs_lrt.md); here I show again only the lines used in section 6.

```r
r <- deparse(body(results))
cat(grep('pfunc <- function|altHypothesis == "(greater|less)"\\)|newStat <- pm|newPvalue <- pfunc\\(\\(', r, value = TRUE), sep = "\n")
```

```
            pfunc <- function(q) pt(q, df = df, lower.tail = FALSE)
            pfunc <- function(q) pnorm(q, lower.tail = FALSE)
            newStat <- pmin(newStatAbove, newStatBelow)
        else if (altHypothesis == "greater") {
            newStat <- pmax((LFC - T)/SE, 0)
            newPvalue <- pfunc((LFC - T)/SE)
        else if (altHypothesis == "less") {
            newStat <- pmin((LFC + T)/SE, 0)
            newPvalue <- pfunc((-T - LFC)/SE)
```

#### dispersions() and dispersions<-

`dispersions()` returns an unnamed numeric, and `dispersions<-` replaces `mcols(object)$dispersion` wholesale. So `dispersions(dds)["gene1"] <- x` appends a new named element, making the length n+1. The error occurs not in `validObject` but earlier, in the assignment to the mcols column (`[[<-`). Assign by row position (a logical index, `which()`).

```r
m <- deparse(getMethod("dispersions<-", c("DESeqDataSet", "numeric"))@.Data)
cat(grep("dispersion <- value|validObject", m, value = TRUE), sep = "\n")
set.seed(1); dd <- DESeq(makeExampleDESeqDataSet(n = 50, m = 6), quiet = TRUE)
is.null(names(dispersions(dd)))
e <- try(dispersions(dd)["gene1"] <- 0.1, silent = TRUE)
cat("error in", as.character(conditionCall(attr(e, "condition"))[[1]]), ":", conditionMessage(attr(e, "condition")), "\n")
```

```
        mcols(object)$dispersion <- value
        validObject(object)
[1] TRUE
error in [[<- : 51 elements in value to replace 50 elements 
```

</details>

<details>
<summary>Common misconceptions</summary>

- "The study text's 0.053147 is the dispersion DESeq2 would give this gene." No. The study text's prior centre 0.08 and SD 0.7 are given values. Actual DESeq2 estimates the trend from the other genes of the same object (0.133 with the seed 1 background), and the prior variance, whose data estimate of 0.0436 hits the lower limit, becomes 0.25 (SD 0.5). So the MAP is 0.0958, and even this varies over 0.093–0.098 with the background genes. But fixing α at the converged α_MAP of 16.3 makes lfcSE, stat and p equal to the hand calculation of 16.4. The latter part of the study text (16.4) is the same operation as DESeq2's Wald calculation.
- "When α grows, the log2FC changes too." With two groups and size factors of 1 it does not change (0.977280 in both (B) and (C)). In the chapter 15 design with a pair term, mu moves with α (median relative difference between mu_initial and mu_final 0.46%, maximum absolute difference 49.3). The size factors are the same, so β moved too.
- "The MAP is always smaller than the gene-wise value." Gene-wise values below the trend are pulled up. All 958 such non-outliers, including gene4, gene10 and gene12 of 15.4, do this, and the chapter 16 example itself goes 0.0254 → 0.0531.
- "p can be reproduced from the `stat` of `results(..., altHypothesis="less")`." The `stat` column is clipped to `pmin(LFC/SE, 0)`, so for genes with a positive LFC `pnorm(stat)` = 0.5. The p itself is the unclipped `pnorm(LFC/SE)`.
- "Halving the two-sided p gives the one-sided p." Wrong unless the sign is checked. In 15.6, 189 of the 764 genes with LFC > 0 looked significant for down with "p/2 < 0.05", but there were 0 with the proper `altHypothesis="less"`.
- "One gene can be fixed with `dispersions(dds)["gene"] <- x`." There are no names, so a new element is appended. So the column assignment in `mcols(object)$dispersion <- value` (`[[<-`) gives a length error (`51 elements in value to replace 50 elements`). Use a logical index or `which()`.
- "With 7 or more replicates per condition, `DESeq()` replaces outliers and refits (and with fewer than 7 it is always the same as the step-by-step listing)." The criterion is not the number of replicates per condition but the size of the design cell (a set of samples with the same model matrix row). In `~pair+condition` every cell has size 1, so there is no refit regardless of the number of replicates. That is why in this data the dispersions and p-values of the step-by-step listing and of `DESeq()` were `all.equal` TRUE. With `fitType="glmGamPoi"`, `minReplicatesForReplace` becomes `Inf`, so there is no refit either.
- "The `stopifnot` passed, so the pairing is right." It checked only the alignment and the rank. In 15.2, a mechanical labelling that shifted the pair numbers by one position per condition also passed. That design had dispersions inflated to 1.37 times the true values (median) and 118 Starvation − Ctrl hits, fewer even than the design with pair dropped altogether (129) (15.3).

</details>

<details>
<summary>Where the results differ from the study text</summary>

| Study text claim | Result | Evidence |
|---|---|---|
| 16.4 Z ≈ 3.380919 | Matches (differs by 1 in the last digit) | 3.3809197 at the converged α (3.380942 with the default tol). The fraction in the study text's formula, 0.977280/0.289057 = 3.380925, is divided by the SE rounded for display |
| 16.4 nominal 95% CI ≈ [0.410727, 1.543833] | Matches (the upper end differs by 1 in the last digit) | z = 1.96: [0.410727, 1.543832], qnorm(0.975): [0.410738, 1.543822]. The study text seems to have used 1.96. The study text's upper end equals 2·0.977280 − 0.410727 = 1.543833 |
| 15.3 `DESeq(dds0, betaPrior=FALSE)` | Matches (the same as the default) | `body(DESeq)`: `if (missing(betaPrior)) betaPrior <- FALSE` |
| 15.6 `altHypothesis="less"` tests a negative effect at threshold 0 | Matches (+ details) | p = `pnorm(LFC/SE)`. The `stat` column is clipped to `pmin(LFC/SE, 0)` |
| 15.1 "running LFC shrinkage also requires apeglm or ashr" | Partly matches | `requireNamespace`: apeglm FALSE, ashr FALSE. In this state `lfcShrink(dd, coef=2, type="apeglm")` fails with "type='apeglm' requires installing the Bioconductor package 'apeglm'", and `type="ashr"` with "... requires installing the CRAN package 'ashr'". But `type="normal"` runs without either (the run output is in the collapsible block "Calling the study text's example code as it is" of [note 07](07_lfc_shrinkage_and_qc.md)) |
| 15.1 the attached `DESeq2_workbook.R` and `DESeq2_single_gene_math.R` | Not confirmed | The files are not in the repository, so the study text's code in the body and the simulated data of section 5 were run instead |
| The figure of figure 3 (p.35) itself | Matches (compared by eye) | The figure redrawn with the same objectives (section 2, the figure of [03](03_dispersion_estimation.md)) and the grid values were viewed side by side with the study text's figure. −2.5 / −3.7 at α = 0.001, −7.38 / −4.23 / −7.88 at the right end α = 0.6, and the prior curve dropping below −12 at α ≈ 0.0033 match the study text's figure. The curve values are in the last two lines of the grid table block |

The rest of the study text's explanations, not in the table (the three α values and the LFC, SE and p of 16.4, the directions of Cox-Reid and shrinkage, the chapter 15 practice code and the warnings of each section, and so on), matched the results of actually running them.

</details>

<details>
<summary>What to look at after this note, and reproduction outputs</summary>

- This note ends the series. The list of notes is in [README.md](README.md), and the places where each note differs from the study text are collected in the collapsible block at the end of each note.
- In real research data, the trend and prior estimation of [03](03_dispersion_estimation.md) replaces the "given prior" of chapter 16. And the shrinkage and QC checks of [07](07_lfc_shrinkage_and_qc.md) (appendix C.2, C.3) come before the prespecified directional test of 15.6.
- The reproduction outputs (`counts.csv`, `metadata.csv`, `dds_wald.rds`, `dispersion_trace.csv`, `sessionInfo.txt`) are not in the repository. To add them, put them in `data/` and record the 15.1 block (seed 2024) and the pre-filter rule with them. The true effect table `truth` of the 15.1 block is not saved to a file, so it has to be regenerated as well.

</details>

---

← Previous: [07. LFC shrinkage and QC](07_lfc_shrinkage_and_qc.md) · Next: [List of notes](README.md) →

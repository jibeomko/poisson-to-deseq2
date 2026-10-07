# 03. How is the dispersion α set from three replicates?

[Note 02](02_negative_binomial.md) showed that counts spread more than the Poisson predicts, and I wondered how the size of that spread, $\alpha$, can be set from three replicates. In this note I compute $\alpha$ three times from the six counts of gene A (0.0148 → 0.0254 → 0.0531) and follow how DESeq2 borrows information from the other genes to set the final $\alpha$.

> Study text chapters 5–7 · Code: [03_dispersion_estimation.R](../../03_dispersion_estimation.R)

## 1. Can't α be computed straight from the sample variance?

Throughout this note the example is gene A from chapter 16 of the study text. Each of the two conditions has three replicates.

| Condition | count | Group mean |
|---|---|---|
| Ctrl | 100, 130, 90 | 106.667 |
| Starvation | 200, 250, 180 | 210 |

The count $K_{ij}$ is the number of reads assigned to gene $i$ in sample $j$. This note mostly deals with one gene, so I often drop the $i$ and write $K_j$. The size factor, the multiplier that corrects for sequencing depth differing between samples, is 1 for all six samples here.

In [note 02](02_negative_binomial.md) the variance of a count was written $\mu+\alpha\mu^2$. $\mu$ is the expected count, the mean count the model expects in that sample. The dispersion $\alpha$ is the size of the spread piled on top of the Poisson variance $\mu$, and $\alpha=0$ is the same as the Poisson. Be sure to remember that $\alpha$ is not the variance itself. So how is $\alpha$ set from six counts?

The first method that comes to mind is to solve the variance formula backwards. Taking the sample variance $s^2$ to be close to $\bar K+\alpha\bar K^2$ gives this.

$$\hat\alpha_{moment}\approx\frac{s^2-\bar K}{\bar K^2}$$

$\bar K$ is the mean of the counts and $s^2$ the sample variance, and the hat over a symbol ($\hat{\ }$) marks an estimate computed from data. The numerator is what is left of the observed spread after removing the Poisson share, and dividing it by the mean² turns it into $\alpha$. This kind of method is called moment estimation.

```r
K <- c(100, 130, 90, 200, 250, 180)
mom_g <- function(k) (var(k) - mean(k)) / mean(k)^2      # moment estimate
round(c(mean = mean(K), var = var(K), all_six = mom_g(K)), 6)
#>        mean         var     all_six 
#>  158.333333 3896.666667    0.149119
round(c(Ctrl = mom_g(K[1:3]), Starvation = mom_g(K[4:6])), 5)
#>       Ctrl Starvation 
#>    0.02871    0.02472
```

Putting all six values together gives 0.149, while splitting by group gives 0.029 and 0.025, a difference of more than five times. The variance of the six values mixed together contains the difference between the Ctrl and Starvation means (106.7 vs 210), so it counts even the treatment effect as spread between replicates.

Splitting by group avoids this problem. But in real experiments the size factors differ between samples and batch or donor effects are present too, so samples often cannot be grouped neatly into sets that share the same mean. A criterion that still works when the mean structure is complicated is needed, and that is the likelihood of the next section. This does not mean DESeq2 throws moment estimation away. It uses it to pick the starting point (initial value) when the computer searches for the best $\alpha$ ("Reading the DESeq2 source (1)" under "Going deeper").

## 2. Which α explains the six counts best?

The likelihood, with the observations held fixed, says how plausible a candidate parameter makes those observations. Here the observations are the six counts and the candidate parameter is $\alpha$.

All you need to grasp is that it runs in the opposite direction to probability. Probability asks "if $\alpha$ were 0.05, how likely would counts like these be?", while likelihood asks "these counts have already come out; which fits better, $\alpha=0.005$ or $\alpha=0.05$?" It is the same probability formula; the only difference is which side is treated as the variable.

Let me compute it. The means are fixed at the group means (106.667, 210). Then for each candidate $\alpha$, the probability of each of the six counts can be found from the negative binomial distribution. The negative binomial (NB) is a count distribution that allows overdispersion, that is, spread between replicates of the same condition larger than the Poisson predicts. In R, `dnbinom()` gives NB probabilities, and the sum of the logs of the six probabilities is the log-likelihood.

```r
mu <- rep(c(mean(K[1:3]), mean(K[4:6])), each = 3)    # fixed at the group means
ll_nb <- function(a) sum(dnbinom(K, mu = mu, size = 1/a, log = TRUE))
cand <- c(0.001, 0.005, 0.015, 0.05, 0.2)
setNames(round(sapply(cand, ll_nb), 3), cand)
#>   0.001   0.005   0.015    0.05     0.2 
#> -29.542 -27.739 -27.028 -28.079 -31.191
a_mle <- optimize(function(th) -ll_nb(exp(th)), c(log(1e-8), log(10)), tol = 1e-10)$minimum
round(exp(a_mle), 6)                                  # the α with the largest log-likelihood
#> [1] 0.014786
```

The value is largest (least negative) at 0.015, in the middle. If $\alpha$ is too small the distribution narrows like a Poisson and values like 130 or 250 become unlikely; if it is too large the distribution spreads too widely and the probability near the observed values drops. Finding the peak exactly with `optimize()` gives 0.014786. The parameter value with the largest likelihood like this is called the MLE (maximum likelihood estimate). The code uses `exp(th)` because it searches by moving $\log\alpha$ instead of $\alpha$; the reason is explained at the end of section 3.

As a formula:

$$\ell(\alpha)=\log L(\alpha)=\sum_{j=1}^{6}\log f_{NB}\big(K_j;\ \mu_j,\ \alpha\big)$$

$f_{NB}(K_j;\mu_j,\alpha)$ is the probability of the count $K_j$ under an NB with mean $\mu_j$ and dispersion $\alpha$, which is `dnbinom(K, mu = mu, size = 1/a)` in the code. $L(\alpha)$ is the product of the six probabilities; they can be multiplied because the samples are assumed independent (how this assumption is handled in designs such as pairs is in "A little more on the likelihood" under "Going deeper"). Multiplying many small probabilities drops to 0 on a computer (underflow), so in practice the sum of the logs, $\ell(\alpha)$, is used. Plugging in $\alpha=0.015$, the six terms sum to −27.028, the third value in the output above.

Two things to add. One is that the likelihood is not "the probability that this $\alpha$ is true". The MLE only picks the location of the peak of the curve; it does not assign probabilities to values of $\alpha$.

The other is that the calculation above fixed the means first and moved only $\alpha$. In theory one can also refit the means for each candidate $\alpha$ before comparing; this is called the profile likelihood. DESeq2's default does not do this. It computes the means once with an initial $\alpha$, fixes them ($\hat\mu^0$), and optimizes only $\alpha$ on top of them.

Computing the means is the job of the GLM (generalized linear model). The GLM is a model that writes the mean of the counts on the log scale as a sum of condition, batch and so on, and it is covered in detail in [note 04](04_glm_condition_batch.md). In a two-group design with all size factors equal, like gene A, the group means do not depend on $\alpha$, so the two methods give the same answer. But in designs where the size factors differ between samples or terms such as batch are mixed in, they can differ.

So when saying "I computed the MLE", the calculation becomes clear only if you also state whether the target is a mean coefficient or $\alpha$, and whether the other parameters were fixed or refitted. Chapter 5 of the study text sums it up the same way.

## 3. What if the mean was estimated from the same data? The Cox-Reid adjustment

Remember why the sample variance is divided by $n-1$ rather than $n$? The sample mean is fitted right to the middle of those very observations, so the distance between the observations and the sample mean is on average shorter than the distance to the true mean. Dividing by $n$ therefore sets the variance too small.

The same thing happens with gene A. The group means fixed in section 2, 106.667 and 210, were computed from these very counts. With two means fitted to the data, the remaining spread tends to look smaller than it really is, and the MLE of 0.014786 contains this effect too. The Cox-Reid adjustment is a correction term that reduces this underestimation of the dispersion caused by estimating the mean from the same data.

$$\ell_{CR}(\alpha)=\sum_{j}\log f_{NB}\big(K_j;\ \hat\mu_j,\ \alpha\big)\;-\;\frac12\log\det\big(X^\top W(\alpha)\,X\big)$$

The first term is the $\ell(\alpha)$ of section 2 unchanged, and the second term is the correction. Let me unpack the symbols one by one.

- $\hat\mu_j$ is the mean estimated from the data. For gene A it is the group mean.
- $X$ is the design matrix, the table that records in numbers which condition each sample belongs to. For gene A it has 6 rows and 2 columns; the first column is all 1 (the reference), and the second is 1 for the Starvation samples.
- $W(\alpha)$ is the matrix with $w_j=\hat\mu_j/(1+\alpha\hat\mu_j)$ on its diagonal. Read $w_j$ as the amount of information one sample $j$ contributes to estimating the mean. The larger $\alpha$, the more erratic the counts, so the information from one sample shrinks.
- $\det(X^\top W X)$ is the determinant of the 2×2 matrix that gathers that information according to the design, a single number summarizing how precisely the mean can be estimated.

Let me put in numbers. At $\alpha=0.025385$, $w_j$ for a Ctrl sample is $106.667/(1+0.025385\times106.667)\approx28.8$, and for a Starvation sample $210/(1+0.025385\times210)\approx33.2$. At $\alpha=0$ (Poisson) $w_j$ would have been the means themselves, 106.7 and 210, so overdispersion cuts the information from one sample a lot. In this design the determinant simplifies to $9\,w_{Ctrl}\,w_{Starvation}$, so the correction term is $-\tfrac12\log(9\times28.77\times33.17)\approx-4.53$.

As $\alpha$ grows, $w_j$ falls, the determinant falls, and the correction term $-\tfrac12\log\det$ grows. So adding this term moves the maximum toward larger $\alpha$.

```r
X  <- cbind(1, rep(0:1, each = 3))                    # reference column, Starvation indicator column
cr <- function(a) { W <- diag(mu/(1 + a*mu)); -0.5*log(det(t(X) %*% W %*% X)) }
ll_cr <- function(a) ll_nb(a) + cr(a)
setNames(round(sapply(cand, cr), 3), cand)            # correction term: grows with α
#>  0.001  0.005  0.015   0.05    0.2 
#> -5.961 -5.534 -4.918 -3.963 -2.673
a_cr <- optimize(function(th) -ll_cr(exp(th)), c(log(1e-8), log(10)), tol = 1e-10)$minimum
round(c(NB_MLE = exp(a_mle), CoxReid = exp(a_cr), ratio = exp(a_cr)/exp(a_mle)), 6)
#>   NB_MLE  CoxReid    ratio 
#> 0.014786 0.025385 1.716818
```

The last line shows that 0.014786 before the correction became 0.025385 after it, 1.72 times larger. These are the same numbers as in section 16.2 of the study text.

But just as the sample variance is divided by $n-1$ because one mean was estimated, couldn't we multiply the variance by $n/(n-p)=6/(6-2)=1.5$ since two means were estimated? Here $n$ is the number of samples and $p$ the number of estimated mean coefficients. I thought so at first too, but in this example the sizes just happen to come out similar; it is not the same operation. The Cox-Reid term is a function whose shape changes with $\alpha$, the design and the estimated means, so how many times it enlarges $\alpha$ differs from gene to gene. The numerical comparison is in "The difference between Cox-Reid and the $n/(n-p)$ correction" under "Going deeper".

Let me check that DESeq2 gives the same value. `estimateDispersionsGeneEst()` is the function in which DESeq2 first estimates the dispersion of each gene.

```r
cd   <- data.frame(condition = factor(rep(c("Ctrl", "Starvation"), each = 3), levels = c("Ctrl", "Starvation")))
dds1 <- DESeqDataSetFromMatrix(matrix(as.integer(K), nrow = 1), cd, ~condition)
sizeFactors(dds1) <- rep(1, 6)
g1 <- estimateDispersionsGeneEst(dds1, quiet = TRUE)
c(dispGeneEst = mcols(g1)$dispGeneEst, dispGeneIter = mcols(g1)$dispGeneIter)
#>  dispGeneEst dispGeneIter 
#>   0.02538756   5.00000000
assays(g1)[["mu"]]                                    # the fixed means
#>          [,1]     [,2]     [,3] [,4] [,5] [,6]
#> [1,] 106.6667 106.6667 106.6667  210  210  210
mcols(estimateDispersionsGeneEst(dds1, useCR = FALSE, quiet = TRUE))$dispGeneEst
#> [1] 0.01478521
```

`dispGeneEst` is 0.0253876, the same as the Cox-Reid maximum above (the difference is at the level of the optimization tolerance), and turning the correction off with `useCR = FALSE` gives 0.0147852, the same as the uncorrected MLE. The stored mean `mu` is also exactly the group means. So `dispGeneEst`, often called the "dispersion MLE" in DESeq2, is strictly speaking not a plain MLE but the maximum of the Cox-Reid adjusted likelihood. Because it is computed from the data of one gene alone, it is called the gene-wise estimate.

Why search over $\log\alpha$ rather than $\alpha$, then? DESeq2 finds the maximum by moving $\theta=\log\alpha$ instead of $\alpha$, for two reasons. First, $\alpha$ must be positive, but $\theta$ can be any real number. Also, a change of $\alpha$ from 0.001 to 0.002 and one from 1 to 2 are both a doubling, and on the log scale both are the same distance ($\log 2$). Meaningful changes in dispersion are relative changes like these. The actual implementation does set a search range for $\theta$ and lower and upper limits for $\alpha$ ("Edge cases of the implementation" under "Going deeper").

### Does α change if the design changes?

What is left as spread depends on what the mean structure is taken to be. Let me check with the two genes of section 6.5 of the study text. The study text calls them A and B; to avoid confusion with gene A of this note, I call them A6.5 and B6.5.

| Gene | Ctrl | Starvation | Feature |
|---|---|---|---|
| A6.5 | 100, 105, 95 | 200, 205, 195 | Large difference in means, small spread within conditions |
| B6.5 | 50, 140, 80 | 120, 300, 170 | A difference in means and a large spread within conditions |

I compare `dispGeneEst` under a design with the condition (`~condition`) and one without it (`~1`, all samples share one mean). The last column applies the moment estimate of section 1 to all six values.

```r
cts <- rbind(A6.5 = c(100,105,95, 200,205,195), B6.5 = c(50,140,80, 120,300,170))
storage.mode(cts) <- "integer"
fit <- function(design) { dds <- DESeqDataSetFromMatrix(cts, cd, design); sizeFactors(dds) <- rep(1, 6)
  mcols(estimateDispersionsGeneEst(dds, quiet = TRUE))$dispGeneEst }
signif(cbind(alpha_condition = fit(~condition), alpha_intercept = fit(~1),
             moment_pooled = apply(cts, 1, function(k) max(0, mom_g(k)))), 4)
#>      alpha_condition alpha_intercept moment_pooled
#> A6.5         1.0e-08          0.1323        0.1276
#> B6.5         2.2e-01          0.3407        0.3681
```

A6.5 gives 1e-08 under `~condition`. The variance within each condition is 25 in both groups, smaller even than the variance the Poisson predicts (100 and 200, equal to the means). So the best $\alpha$ heads toward 0 and hits DESeq2's lower limit `minDisp = 1e-8`. Switching to `~1` turns the whole difference between 100 and 200 into remaining spread, and it jumps to 0.132. B6.5 really has a large spread within conditions, so it is 0.220 even under `~condition`, and under `~1` the difference in means is piled on top, giving 0.341.

This is why changing the design means re-estimating the dispersions from scratch. If you changed the design with `design(dds) <-`, you have to run `estimateDispersions()` again.

## 4. How much does an α from three replicates wobble?

Let me make 2000 genes with a true $\alpha$ of 0.1. All have a mean of 100, 3 vs 3, with no condition difference. If the gene-wise value is computed for each gene, will they gather around 0.1?

```r
set.seed(2)
sim <- matrix(rnbinom(2000*6, mu = 100, size = 1/0.1), ncol = 6); storage.mode(sim) <- "integer"
ds  <- DESeqDataSetFromMatrix(sim, data.frame(condition = factor(rep(c("a", "b"), each = 3))), ~condition)
sizeFactors(ds) <- rep(1, 6)
gw  <- mcols(estimateDispersionsGeneEst(ds, quiet = TRUE))$dispGeneEst
signif(c(quantile(gw, c(0.01, 0.05, 0.5, 0.95, 0.99)), at_minDisp = mean(gw <= 1e-8),
         lt_0.01 = mean(gw < 0.01), gt_0.5 = mean(gw > 0.5), gt_1 = mean(gw > 1)), 3)
#>         1%         5%        50%        95%        99% at_minDisp    lt_0.01 
#>   1.00e-08   1.07e-02   8.25e-02   2.45e-01   3.70e-01   2.90e-02   4.90e-02 
#>     gt_0.5       gt_1 
#>   5.00e-04   0.00e+00
```

Even though every gene's true value is 0.1, the central 90% (the 5–95% interval) spreads from 0.011 to 0.25, that is, from 1/10 of the true value to 2.5 times it. About 3% stuck to the lower limit 1e-8, and at this mean level there were no values above 1. When the three values happen to come out similar the estimate is small; when one value jumps it is large. Without more replicates, this wobble is hard to reduce from the data of one gene alone.

What the shrinkage discussed below fixes is the uncertainty of this estimate. It is not erasing the fact that genes really differ in variability.

### Genes with similar expression: the trend

But there are thousands of genes. Genes with a similar mean expression also tend to have dispersions in a similar range, and the curve that summarizes this tendency is the trend. That is, the trend is a curve showing where the dispersion roughly lies as a function of mean expression. DESeq2's default trend has this shape.

$$\alpha_{tr}(\bar q_i)=a_0+\frac{a_1}{\bar q_i},\qquad \bar q_i=\frac1n\sum_{j}\frac{K_{ij}}{s_j}$$

$\bar q_i$ is the mean normalized count of gene $i$. A normalized count is the count divided by the size factor $s_j$, and $n$ is the number of samples. $a_0$ is the value $\alpha_{tr}$ approaches at very high expression, and $a_1/\bar q_i$ is the part that grows as expression gets lower. $a_0$ and $a_1$ are set by fitting the curve to the gene-wise values of all genes.

I make example data with 2000 genes and 6 samples (3 vs 3) using DESeq2's built-in `makeExampleDESeqDataSet()` and fit the trend. This data is used again in section 6.

```r
set.seed(1)
dds <- makeExampleDESeqDataSet(n = 2000, m = 6, betaSD = 1)   # design ~condition, 3 vs 3
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, quiet = TRUE)
fn  <- dispersionFunction(dds)                                # the fitted trend function
round(attr(fn, "coefficients"), 5)
#> asymptDisp  extraPois 
#>    0.09389    6.19823
sapply(c(q10 = 10, q1000 = 1000), function(q) c(alpha_tr = unname(fn(q)), Var = unname(q + fn(q)*q^2)))
#>                 q10        q1000
#> alpha_tr  0.7137142 1.000892e-01
#> Var      81.3714206 1.010892e+05
```

The fit gives $a_0=0.0939$ and $a_1=6.198$. Plugging in, $\alpha_{tr}=0.0939+6.198/10\approx0.714$ at $\bar q=10$, and about 0.100 at $\bar q=1000$.

The last output contains two things to be careful about when reading the trend. First, $\alpha$ fell 7-fold, yet the count variance $\mu+\alpha\mu^2$ rose more than 1000-fold, from 81 to about 101,000. This is why $\alpha$ and the variance must not be confused. And "higher expression means lower $\alpha$" is not a law of biology. The trend only summarizes where genes with similar means roughly lie in this data. If the shape of the curve does not fit the data, use `fitType = "local"` or `"mean"`.

### Turning the trend into a prior

The trend is an expectation: "at this expression level, $\alpha$ is roughly about here". Writing this expectation as a probability distribution gives the prior. A prior is a distribution for where a parameter is likely to be before the data are seen, and DESeq2 puts a normal prior on $\log\alpha$.

$$\theta_i=\log\alpha_i\ \sim\ N\big(m_i,\ \sigma_d^2\big),\qquad m_i=\log\alpha_{tr}(\bar q_i)$$

The centre of the prior, $m_i$, is the log of the trend value at gene $i$'s expression level, so it differs from gene to gene. The width (variance) of the prior, $\sigma_d^2$, on the other hand, represents how widely the genes' true $\log\alpha$ scatter around the trend, and there is one for all genes. Because both the centre and the width are estimated from all genes, this approach is called empirical Bayes.

The width $\sigma_d^2$ is set like this. If you measure how widely the gene-wise $\log\alpha$ scatter around the trend, it mixes real differences between genes with the estimation noise seen in the simulation just above. So the share expected from estimation noise is subtracted from the observed scatter. The expected noise is computed with 4 degrees of freedom, the number of samples, 6, minus the number of coefficients, 2; `trigamma()` in the code is the function that gives the size of this noise for a given number of degrees of freedom (the reasoning is in "How is the prior width σ_d² computed?" under "Going deeper").

```r
pv <- attr(fn, "dispPriorVar")
round(c(observed = attr(fn, "varLogDispEsts"), expected_noise = trigamma((6 - 2)/2), dispPriorVar = pv), 3)
#>       observed expected_noise   dispPriorVar 
#>          0.679          0.645          0.250
```

The first two numbers are almost the same. Subtracting the expected noise of 0.645 from the observed scatter of 0.679 leaves only 0.034. DESeq2 sets a lower limit of 0.25 so that the width does not get too narrow, so $\sigma_d^2$ in this example became 0.25. The scatter was almost all noise because this example data was made with the true $\alpha$ placed exactly on the curve (which can be checked in the same block).

## 5. How are a gene's own value and the trend mixed? MAP

Now there are two pieces of information: the $\ell_{CR}$ given by the gene's own counts, and the prior supplied by the other genes. By Bayes' rule, the distribution after seeing the data (the posterior) is proportional to the product of the likelihood and the prior, and taking logs turns the product into a sum. The value where that sum is largest is the MAP (maximum a posteriori), the most plausible value when the likelihood and the prior are considered together.

$$\hat\theta_{MAP}=\arg\max_{\theta}\Big[\ \ell_{CR}(e^{\theta})\;-\;\frac{(\theta-m_i)^2}{2\sigma_d^2}\ \Big],\qquad \hat\alpha_{MAP}=e^{\hat\theta_{MAP}}$$

The first term inside the brackets, $\ell_{CR}(e^\theta)$, is the Cox-Reid adjusted likelihood of section 3 with $\alpha=e^\theta$. The second term, $(\theta-m_i)^2/(2\sigma_d^2)$, is a penalty from the prior: it grows quadratically as $\theta$ moves away from the trend centre $m_i$, and it is stronger the smaller $\sigma_d^2$ is. $\arg\max_\theta$ means picking the $\theta$ that makes the bracket largest.

The trend cannot be estimated from gene A alone. So section 16.3 of the study text gives a teaching prior: trend value 0.08 and prior SD 0.7, that is, $\theta\sim N(\log 0.08,\ 0.7^2)$, which is not estimated from data.

```r
post  <- function(th, m = log(0.08), s2 = 0.7^2) ll_cr(exp(th)) - (th - m)^2/(2*s2)
a_map <- optimize(function(th) -post(th), c(log(1e-8), log(10)), tol = 1e-10)$minimum
round(c(NB_MLE = exp(a_mle), CoxReid = exp(a_cr), MAP = exp(a_map)), 6)
#>   NB_MLE  CoxReid      MAP 
#> 0.014786 0.025385 0.053147
```

The MAP is 0.053147. It lies between the gene's own value, 0.025385, and the prior centre, 0.08, and is larger than its own value. Pulling estimates that carry little information toward the overall tendency is called shrinkage, and the target it pulls toward is not 0 but the trend. So gene A, which was below the trend, went up.

![Gene A: three objective functions over alpha](../../figures/03_gene_a_objectives.png)

Figure 1. The three objective functions for gene A plotted against $\alpha$ (figure 3 of the study text). Each curve is shifted so that its own maximum is 0, so don't compare heights between curves; just watch the peak move to the right, 0.0148 → 0.0254 → 0.0531.

<details>
<summary>Code for the figure</summary>

```r
ag <- 10^seq(-3, log10(0.6), length.out = 400)
curves <- cbind(sapply(ag, ll_nb) - ll_nb(exp(a_mle)),
                sapply(ag, ll_cr) - ll_cr(exp(a_cr)),
                sapply(log(ag), post) - post(a_map))
peaks <- exp(c(a_mle, a_cr, a_map)); cols <- c("#2a78d6", "#eb6834", "#1baf7a"); ltys <- c(1, 2, 4)
out <- "../figures/03_gene_a_objectives.png"
png(out, width = 7, height = 4.5, units = "in", res = 150, bg = "white")
par(mar = c(4.2, 4.6, 1.2, 1), las = 1, col.axis = "#52514e", col.lab = "#0b0b0b", fg = "#8a8984")
matplot(ag, curves, type = "l", log = "x", lty = ltys, lwd = 2.2, col = cols, ylim = c(-12, 1.2),
        xlab = "dispersion alpha (log scale)", ylab = "objective minus its own maximum", xaxt = "n")
axis(1, at = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3), labels = c("0.001", "0.003", "0.01", "0.03", "0.1", "0.3"))
abline(v = 0.08, col = "grey55", lty = 3)
text(0.08, 1.0, "prior centre 0.08", pos = 4, cex = 0.75, col = "#52514e")
points(peaks, rep(0, 3), pch = 19, col = cols)
text(peaks, 0, sprintf("%.4f", peaks), pos = 3, cex = 0.72, col = "#0b0b0b")
legend("bottomright", lty = ltys, lwd = 2.2, col = cols, bty = "n", cex = 0.8, text.col = "#0b0b0b",
       legend = c("NB likelihood (MLE)", "Cox-Reid adjusted", "Cox-Reid adjusted + teaching prior (MAP)"))
invisible(dev.off())
file.exists(out)
#> [1] TRUE
```

</details>

How far it is pulled differs from gene to gene. If the $\ell_{CR}$ curve on the $\theta$ scale is roughly a parabola (a polynomial of degree 2) near its peak, the MAP is approximated by a weighted average of the two values. A parabolic log-likelihood means the likelihood has the shape of a normal distribution, and the prior is normal too, so combining the two gives a weighted average.

$$\hat\theta_{MAP}\approx\frac{\hat\theta_{gw}/v_i+m_i/\sigma_d^2}{1/v_i+1/\sigma_d^2}$$

$\hat\theta_{gw}$ is the log of the gene-wise value, $\log 0.025385$ for gene A. $v_i$ is the uncertainty of the gene-wise value, the reciprocal of how sharp the $\ell_{CR}$ curve on $\theta$ is at its peak. If the curve is sharp (lots of information), $v_i$ is small. The weights are $1/v_i$ and $1/\sigma_d^2$, so a gene with lots of information keeps its own value, and a gene with little information moves a long way toward the trend. A wide prior allows more differences between genes, and a narrow one increases the influence of the trend. The prior matters more with fewer samples, but this means shrinkage is not applied in the same proportion to every gene.

```r
th0 <- a_cr; h <- 1e-4
v <- -1 / ((ll_cr(exp(th0 + h)) - 2*ll_cr(exp(th0)) + ll_cr(exp(th0 - h))) / h^2)   # reciprocal of the curvature
th_approx <- (th0/v + log(0.08)/0.49) / (1/v + 1/0.49)
round(c(v = v, approx = exp(th_approx), exact_MAP = exp(a_map)), 6)
#>         v    approx exact_MAP 
#>  0.795046  0.051642  0.053147
round(c(w_gw_approx = (1/v)/(1/v + 1/0.49), w_gw_exact = (a_map - log(0.08))/(a_cr - log(0.08))), 4)
#> w_gw_approx  w_gw_exact 
#>      0.3813      0.3563
```

With $v_i=0.795$ and $\sigma_d^2=0.49$, the weight on the gene-wise side is $1.258/(1.258+2.041)\approx0.381$. Mixing with this weight gives 0.0516, close to the exact MAP of 0.0531. They are not identical because $\ell_{CR}$ is not exactly of degree 2 in $\theta$. The gene-wise weight back-calculated from the actual MAP location is 0.356. So this formula is only an approximation to aid understanding; DESeq2 maximizes the MAP objective above directly.

One more thing to watch: the MAP depends on the scale on which the maximum is sought. Turning the same prior into a density on the $\alpha$ scale and finding its maximum gives 0.0384 ("The scale of the MAP" under "Going deeper"). DESeq2 optimizes on the $\theta=\log\alpha$ scale, so the answer in DESeq2's manner is 0.0531.

It is also worth pointing out that Cox-Reid and shrinkage are different steps. Cox-Reid corrects for the effect of estimating the mean within the data of one gene (0.0148 → 0.0254), while shrinkage borrows information between genes (0.0254 → 0.0531).

## 6. Checking in the DESeq2 results table

`estimateDispersions()` splits the process so far into three functions and calls them in turn.

1. `estimateDispersionsGeneEst()` computes the mean $\hat\mu^0$ of each gene once, fixes it, and computes the gene-wise value from the Cox-Reid adjusted likelihood → `dispGeneEst`
2. `estimateDispersionsFit()` fits the trend to the gene-wise values → `dispFit`
3. `estimateDispersionsMAP()` sets the prior width, computes the MAP, and then applies an exception to genes far above the trend (outliers), explained later in this section → `dispMAP`, `dispOutlier`, `dispersion`

Let me pick a few genes from the example data made in section 4 and look at the numbers. `baseMean` is the mean normalized count ($\bar q_i$).

```r
df   <- as.data.frame(mcols(dds)); ok <- !df$allZero     # per-gene results table
pick <- c("gene733", "gene798", "gene12", "gene15", "gene1", "gene2")
round(df[pick, c("baseMean", "dispGeneEst", "dispFit", "dispMAP", "dispOutlier", "dispersion")], 4)
#>         baseMean dispGeneEst dispFit dispMAP dispOutlier dispersion
#> gene733  14.3625      3.9334  0.5254  1.0829           1     3.9334
#> gene798   6.2682      7.3963  1.0827  1.9315           1     7.3963
#> gene12   20.8986      1.2598  0.3905  0.5537           0     0.5537
#> gene15  113.6308      0.3119  0.1484  0.1985           0     0.1985
#> gene1     1.8177      1.2302  3.5039  2.9337           0     2.9337
#> gene2    17.2139      0.4213  0.4540  0.4464           0     0.4464
```

gene12 and gene15 have gene-wise values above the trend, so their MAP came down. Conversely, gene1 and gene2 have gene-wise values below the trend, so their MAP went up, the same case as gene A. gene733 and gene798 have gene-wise values far above the trend. These two have `dispOutlier` equal to 1, and their final `dispersion` is the gene-wise value, not the MAP.

Plotting the whole example data, including these six genes, in one figure gives this.

![Dispersion estimates of the example data](../../figures/03_dispersion_shrinkage.png)

Figure 2. Grey points are gene-wise values, the orange line is the trend, and blue points are the final values. The triangles at the bottom gather gene-wise values below 0.01 (including the lower limit 1e-8) at a height of 0.01. The blue points cluster closer to the trend than the grey points, and only the outliers, marked with circles, stay at their gene-wise values.

<details>
<summary>Code for the figure</summary>

```r
out <- "../figures/03_dispersion_shrinkage.png"
png(out, width = 7, height = 4.5, units = "in", res = 150, bg = "white")
par(mar = c(4.2, 4.6, 1.2, 1), las = 1, col.axis = "#52514e", fg = "#8a8984")
plotDispEsts(dds, ymin = 1e-2,   # gene-wise values below 0.01 as triangles at the bottom
             genecol = "grey70", fitcol = "#eb6834", finalcol = "#2a78d6",
             legend = FALSE, xlab = "mean of normalized counts (baseMean)", ylab = "dispersion", xaxt = "n")
axis(1, at = 10^(-1:3), labels = c("0.1", "1", "10", "100", "1000"))
g <- c("gene1", "gene12")                                # examples that moved up / down
arrows(df[g, "baseMean"], df[g, "dispGeneEst"], df[g, "baseMean"], df[g, "dispersion"],
       length = 0.06, lwd = 1.4, col = "#0b0b0b")
text(df[c(g, "gene733"), "baseMean"], df[c(g, "gene733"), "dispGeneEst"], c(g, "gene733"),
     pos = 4, cex = 0.75, col = "#0b0b0b")
legend("topright", bty = "n", cex = 0.8, text.col = "#0b0b0b", pch = c(20, 16, 16, 1), pt.cex = c(1, 0.8, 0.8, 1.3),
       col = c("grey70", "#eb6834", "#2a78d6", "#2a78d6"),
       legend = c("gene-wise (dispGeneEst)", "trend (dispFit)", "final (dispersion)", "outlier: final = gene-wise"))
invisible(dev.off())
file.exists(out)
#> [1] TRUE
```

</details>

There is a reason for not pulling down genes far above the trend. The SE (standard error), which says how much an estimate would move if the same experiment were repeated, comes out smaller the smaller $\alpha$ is, other things being equal. If that gene really is highly variable, then forcing $\alpha$ down makes the SE small and the p-value too small, risking significance that is not there. So DESeq2 makes an exception only for genes that deviate far upward and uses their gene-wise value as it is. The criterion is this.

$$\log(\texttt{dispGeneEst})>\log(\texttt{dispFit})+2\sqrt{\texttt{varLogDispEsts}}$$

In words, a gene is an outlier if the log of its gene-wise value is more than twice the SD of the observed scatter above the log of the trend. `varLogDispEsts` is the observed scatter of 0.679 seen in section 4 (about 0.82 as an SD), not the prior width $\sigma_d^2$.

```r
vld <- attr(fn, "varLogDispEsts")
myOut <- with(df[ok, ], log(dispGeneEst) > log(dispFit) + 2*sqrt(vld))
c(rule_reproduced = identical(myOut, df$dispOutlier[ok]), n_outlier = sum(myOut),
  if_prior_SD = sum(with(df[ok, ], log(dispGeneEst) > log(dispFit) + 2*sqrt(pv))))
#> rule_reproduced       n_outlier     if_prior_SD 
#>               1               8              79
with(df[ok & !df$dispOutlier, ], c(non_outlier = length(dispMAP), MAP_up = sum(dispMAP > dispGeneEst),
  below_trend = sum(dispGeneEst < dispFit), below_trend_down = sum(dispGeneEst < dispFit & dispMAP < dispGeneEst)))
#>      non_outlier           MAP_up      below_trend below_trend_down 
#>             1969             1326             1336                0
identical(dispersions(dds), df$dispersion)
#> [1] TRUE
c(allZero = sum(!ok), allZero_all_NA = all(is.na(df[!ok, c("dispGeneEst", "dispFit", "dispMAP", "dispersion", "dispOutlier")])))
#>        allZero allZero_all_NA 
#>             23              1
```

Reading the output in order: first, the rule is reproduced exactly, and there are 8 outliers. Had the prior SD ($\sqrt{0.25}$) been put into the same formula, 79 would have been caught.

Of the 1969 non-outliers, 1326 (67%) had a MAP larger than their gene-wise value. Shrinkage is not an operation that makes values smaller. And of the 1336 genes below the trend, not a single one went further down at the MAP step (the 10 that stayed put are genes at the upper limit 10). Genes below the trend move toward a larger or unchanged $\alpha$, that is, the conservative side where the p-value gets less small, so the exception only needs to be on the upper side.

`dispersions(dds)` pulls out the `mcols(dds)$dispersion` column, which starts from `dispMAP` and replaces only the outliers with `dispGeneEst`. The 23 genes whose counts are 0 in every sample are left out of the calculation, so all five columns are NA for them.

The stored `dispMAP` does not always lie between `dispGeneEst` and `dispFit`. 3 of the 1969 lay outside, because the stored `dispGeneEst` was not the exact maximum of $\ell_{CR}$ ("The 3 genes whose MAP fell outside the gene-wise value and the trend" under "Going deeper").

Dispersion outliers are a different concept from the Cook's distance outliers that come up later. A dispersion outlier is a judgement about the dispersion of a gene as a whole, while Cook's distance is a judgement about how much one sample's count influences the coefficient estimates ([note 07](07_lfc_shrinkage_and_qc.md)).

## Summary

- Moment estimation, which subtracts the Poisson share from the sample variance, works only when all samples have the same mean. For gene A it gave 0.149 with the six values mixed and 0.029 and 0.025 per group.
- So candidate values of $\alpha$ are compared with the likelihood. The NB MLE with the means fixed at the group means is 0.014786, and the Cox-Reid term $-\tfrac12\log\det(X^\top WX)$, which corrects for the share lost by estimating the means from the same data, makes it 0.025385. This is DESeq2's `dispGeneEst`, and it changes when the design changes.
- With few replicates, gene-wise values wobble a lot (true value 0.1, 5–95% interval 0.011–0.25). So DESeq2 puts a prior on $\log\alpha$ centred on the trend over expression levels, and estimates its width from all genes too (empirical Bayes).
- The MAP mixes the gene's own value and the trend according to how much information there is, so it can also move up (gene A with the teaching prior: 0.0254 → 0.0531). Only outliers far above the trend use the gene-wise value instead of the MAP, and the final values are in `dispersions(dds)`.
- Next comes the step that fixes this final $\alpha$ and computes the condition coefficient and its SE ([04](04_glm_condition_batch.md)). The calculation that carries gene A's 0.053147 through to the SE and p-value is in [08](08_one_gene_end_to_end.md).

## Exercises

Problems 6–10 of appendix A of the study text.

**Problem 6.** In what situation is the explanation 'the dispersion MLE is just the sample variance minus the Poisson variance' intuitively useful, and what does it leave out of the general DESeq2 process?

<details>
<summary>Solution</summary>

$\hat\alpha\approx(s^2-\bar K)/\bar K^2$ is meaningful only in a single group where every sample has the same mean and a size factor of 1. Compared with the general DESeq2 process, four things are missing.

1. The mean structure: the size factors (offsets) and the differences in means due to condition, pair and batch. For gene A, the pooled moment with the group difference mixed in is 0.149, the NB MLE with the condition included and the means fixed is 0.0148, and the within-group moments are 0.029/0.025.
2. The likelihood optimization itself.
3. The Cox-Reid adjustment (0.0148 → 0.0254).
4. The trend, prior and MAP that follow (with the teaching prior → 0.053).

In the DESeq2 source, the moment estimate `momentsDispEstimate` enters the initial value through `pmin(roughDisp, momentsDisp)`, and for gene A the rough side (0.0232) is chosen.

</details>

**Problem 7.** Explain the difference in how the fitted mean is handled between a strict profile likelihood and DESeq2's default gene-wise implementation.

<details>
<summary>Solution</summary>

The profile likelihood refits the mean coefficients $\hat b(\alpha)$ for each candidate $\alpha$ to build the curve $\ell(\hat b(\alpha),\alpha)$. DESeq2 computes $\hat\mu^0$ once with an initial $\alpha$ (`pmin(rough, moments)`) (via `linearModelMuNormalized` or `fitNbinomGLMs`) and passes it as a fixed value to `fitDispWrapper(mu_hatSEXP = fitMu, ...)`. The outer mean→dispersion cycle defaults to `niter = 1`. Even with `niter` set to 2 or more, from the second iteration on only genes with $|\Delta\log\alpha|>0.05$ are refitted. For gene A the group means do not depend on $\alpha$, so the two approaches agree, but in general designs with offsets mixed in they can differ.

</details>

**Problem 8.** What problem does each of the Cox–Reid term and the MAP prior penalty address? Can both be lumped together under the word shrinkage?

<details>
<summary>Solution</summary>

The Cox-Reid term $-\tfrac12\log\det(X^\top W X)$ corrects, within one gene, the underestimation of the dispersion that arises from estimating $p$ mean coefficients from the same data. With the mean fixed, it only moves $\alpha$ to an equal or larger value. The prior penalty $-(\theta-m_i)^2/(2\sigma_d^2)$ combines an unstable gene-wise estimate with the trend learned from all genes; its direction is toward the trend, so it can go up or down. In the source too, `useCR` and `usePrior` are separate switches, with `usePriorSEXP = FALSE` in the gene-wise step and `TRUE` in the MAP step. Calling both shrinkage hides the purpose and direction of the two operations. For gene A the former is 0.015 → 0.025 and the latter 0.025 → 0.053.

</details>

**Problem 9.** With gene-wise α=0.02 and trend α=0.1, the final α came out at 0.05. Since this is shrinkage yet the value grew, is it an error?

<details>
<summary>Solution</summary>

It is not an error. The prior centre is not 0 but $\log\alpha_{tr}$, so genes below the trend move up. Gene A is exactly this situation: gene-wise 0.0254, prior centre 0.08, MAP 0.0531, and the weighted-average approximation, with weights from $v_i=0.795$ and $\sigma_d^2=0.49$ on the $\theta$ scale, gives 0.0516.

The numbers in the problem fit this picture too. On the $\theta$ scale, $\log 0.05$ lies between $\log 0.02$ and $\log 0.1$ at a gene-wise weight of $(\log0.05-\log0.1)/(\log0.02-\log0.1)=0.43$ (`ex9_weight_gw` in "The scale of the MAP" under "Going deeper"). So it can be read as a gene with $1/v_i : 1/\sigma_d^2 \approx 0.43 : 0.57$. In the example data too, 1326 of the 1969 non-outliers moved up.

</details>

**Problem 10.** Must `dispersions(dds)` agree with `mcols(dds)$dispMAP` for every gene?

<details>
<summary>Solution</summary>

No. `dispersions(dds)` is `mcols(dds)$dispersion`, a column that starts from `dispMAP` and overwrites only the genes flagged `dispOutlier` with `dispGeneEst` (source: `dispersionFinal[dispOutlier] <- dispGeneEst[dispOutlier]`). The test is `log(dispGeneEst) > log(dispFit) + 2 * sqrt(varLogDispEsts)`. In the example, 8 of the 1977 genes were outliers; for those 8, `dispersion == dispGeneEst`, and for the rest `dispersion == dispMAP` ("Reading the DESeq2 source (4)" under "Going deeper"). For all-zero genes (23 in the example), `dispGeneEst`, `dispFit`, `dispMAP`, `dispersion` and `dispOutlier` are all NA (`allZero_all_NA` in section 6).

</details>

## Going deeper

<details>
<summary>Reading the DESeq2 source (1): the three steps of estimateDispersions and the gene-wise step</summary>

`estimateDispersions()` calls three functions in order. From the body of `getMethod("estimateDispersions", "DESeqDataSet")` I excerpted only the three call lines (with no arguments omitted).

```r
# excerpt (not run): the three call lines in the body of getMethod("estimateDispersions", "DESeqDataSet")
object <- estimateDispersionsGeneEst(object, maxit = maxit, useCR = useCR, weightThreshold = weightThreshold,
    quiet = quiet, modelMatrix = modelMatrix, minmu = minmu, type = dispersionEstimator)
object <- estimateDispersionsFit(object, fitType = fitType, quiet = quiet)
object <- estimateDispersionsMAP(object, maxit = maxit, useCR = useCR, weightThreshold = weightThreshold,
    quiet = quiet, modelMatrix = modelMatrix, type = dispersionEstimator)
```

The default arguments of the gene-wise step are these.

```r
args(DESeq2:::estimateDispersionsGeneEst)
#> function (object, minDisp = 1e-08, kappa_0 = 1, dispTol = 1e-06, 
#>     maxit = 100, useCR = TRUE, weightThreshold = 0.01, quiet = FALSE, 
#>     modelMatrix = NULL, niter = 1, linearMu = NULL, minmu = if (type == 
#>         "glmGamPoi") 1e-06 else 0.5, alphaInit = NULL, type = c("DESeq2", 
#>         "glmGamPoi")) 
#> NULL
```

`niter = 1` is the default, and `maxit = 100` is a separate argument. Below is an excerpted and condensed version of the body (`...` marks omitted arguments, and the comments are mine).

```r
# excerpt (not run): body of estimateDispersionsGeneEst, excerpted and condensed
if (nrow(modelMatrix) == ncol(modelMatrix)) { stop("the number of samples and the number of model coefficients are equal, ...") }
objectNZ <- object[!mcols(object)$allZero, , drop = FALSE]      # exclude all-zero genes
# initial dispersion: the smaller of the rough residual estimate and the moment estimate, clipped to [minDisp, maxDisp]
roughDisp <- roughDispEstimate(y = counts(objectNZ, normalized = TRUE), x = modelMatrix)
momentsDisp <- momentsDispEstimate(objectNZ)
alpha_hat <- pmin(roughDisp, momentsDisp)
alpha_hat <- alpha_hat_new <- alpha_init <- pmin(pmax(minDisp, alpha_hat), maxDisp)
# initial mean path: linear mean for a group design whose number of unique rows equals the number of columns, NB-GLM otherwise. Always GLM with weights
if (is.null(linearMu)) {
    modelMatrixGroups <- modelMatrixGroups(modelMatrix)
    linearMu <- nlevels(modelMatrixGroups) == ncol(modelMatrix)
    if (useWeights) { linearMu <- FALSE }
}
for (iter in seq_len(niter)) {
    if (!linearMu) { fit <- fitNbinomGLMs(..., alpha_hat = alpha_hat[fitidx], ...); fitMu <- fit$mu }
    else           { fitMu <- linearModelMuNormalized(objectNZ[fitidx, , drop = FALSE], modelMatrix) }
    fitMu[fitMu < minmu] <- minmu
    dispRes <- fitDispWrapper(ySEXP = counts(objectNZ)[fitidx, , drop = FALSE], xSEXP = modelMatrix,
        mu_hatSEXP = fitMu, log_alphaSEXP = log(alpha_hat)[fitidx], ...,
        min_log_alphaSEXP = log(minDisp/10), kappa_0SEXP = kappa_0, tolSEXP = dispTol,
        maxitSEXP = maxit, usePriorSEXP = FALSE, ..., useCRSEXP = useCR)
    fitidx <- abs(log(alpha_hat_new) - log(alpha_hat)) > 0.05
    ...
}
if (niter == 1) { noIncrease <- last_lp < initial_lp + abs(initial_lp)/1e+06
                  dispGeneEst[which(noIncrease)] <- alpha_init[which(noIncrease)] }
dispGeneEstConv <- dispIter < maxit & !(dispIter == 1)
refitDisp <- !dispGeneEstConv & dispGeneEst > minDisp * 10
if (sum(refitDisp) > 0) { dispGrid <- fitDispGridWrapper(...); dispGeneEst[refitDisp] <- dispGrid }
```

Matching the steps described in section 6.3 of the study text to the source gives this.

| Step in section 6.3 of the study text | What it is in the source |
|---|---|
| Preparation (size factors, design matrix, excluding all-zero genes) | `estimateDispersions()` stops if there are no size factors. If the number of samples = the number of coefficients, `stop()`. It computes with the all-zero genes removed via `objectNZ <- object[!mcols(object)$allZero, ]`, so those genes have NA for `dispGeneEst`, `dispFit`, `dispMAP`, `dispersion` and `dispOutlier` (checked in `allZero_all_NA` in section 6) |
| Initial dispersion | `pmin(roughDisp, momentsDisp)`, clipped to `[minDisp, max(10, ncol)]`. If the optimization fails to raise the objective (`noIncrease`), this initial value becomes `dispGeneEst` as it is |
| Initial mean $\hat\mu^0$ | Group designs such as 2-group → `linearModelMuNormalized` (QR linear mean × size factor); if there are more unique rows than columns, as in `~batch + condition`, `fitNbinomGLMs`; with sample weights, always `fitNbinomGLMs` |
| $\hat\mu^0$ fixed, optimize $\theta=\log\alpha$ | `fitDispWrapper(mu_hatSEXP = fitMu, log_alphaSEXP = log(alpha_hat), usePriorSEXP = FALSE)` |
| `niter` vs `maxit` | `niter` = number of outer mean→dispersion cycles (default 1), `maxit` = upper limit of the inner C++ line-search iterations (default 100) |
| Convergence handling | If `dispIter == maxit` or `== 1` and `dispGeneEst > 10·minDisp`, re-optimize on a grid (`fitDispGridWrapper`: a 20-point grid on the R side, and the C++ `fitDispGrid` sweeps a fine grid around the optimum once more). Genes at the lower limit (A6.5 of section 6.5 of the study text) are not grid targets. The MAP step also uses a grid fallback, but its condition is only `iter == maxit` (no `iter == 1` or `minDisp` condition), and the grid call is fixed at `useCRSEXP = TRUE` regardless of the `useCR` argument |

I checked the `linearMu` decision by running real designs. Only colData is needed, no `dds`.

```r
cd6 <- data.frame(condition = factor(rep(c("Ctrl", "Starvation"), each = 3)), batch = factor(rep(c("b1", "b2", "b3"), 2)))
for (f in list(~condition, ~batch + condition)) {
  mm <- model.matrix(f, cd6); g <- nlevels(DESeq2:::modelMatrixGroups(mm))
  cat(deparse(f), ": n =", nrow(mm), " p =", ncol(mm), " groups =", g, " linearMu =", g == ncol(mm), "\n")
}
#> ~condition : n = 6  p = 2  groups = 2  linearMu = TRUE 
#> ~batch + condition : n = 6  p = 4  groups = 6  linearMu = FALSE
```

</details>

<details>
<summary>Reading the DESeq2 source (2): the C++ objective function fitDisp</summary>

The installed package contains only the compiled `.so`, so I read `src/DESeq2.cpp` in the DESeq2 1.50.2 source tarball of Bioconductor 3.22.

```cpp
# excerpt (not run): src/DESeq2.cpp in the DESeq2 1.50.2 source tarball (fitDisp, fitDispGrid)
ll_part = sum(lgamma(y + alpha_neg1) - Rf_lgammafn(alpha_neg1) - y * log(mu + alpha_neg1) - alpha_neg1 * log(1.0 + mu * alpha));  // alpha_neg1 = 1/alpha
arma::vec w_diag = pow(pow(mu, -1) + alpha, -1);      // w_j = 1/(1/mu_j + alpha)
b = x.t() * (x.each_col() % w_diag);                   // X^T W X
cr_term = -0.5 * log(det(b));                          // Cox-Reid term (when useCR)
prior_part = -0.5 * R_pow_di(log_alpha - log_alpha_prior_mean,2)/log_alpha_prior_sigmasq;  // when usePrior
double res =  ll_part + prior_part + cr_term;
// iteration (fitDisp): a_propose = a + kappa * dlp;
//   theta_hat_kappa = -1.0 * lp - kappa * epsilon * R_pow_di(dlp, 2);  if the Armijo condition fails, kappa = kappa / 2.0;
//   break if change < tol
// bounds: if (a_propose < -30.0) { kappa = (-30.0 - a)/dlp; }   if (a_propose > 10.0) { kappa = (10.0 - a)/dlp; }
//       if (a < min_log_alpha) { break; }                        // min_log_alpha = log(minDisp/10)
// grid fallback (fitDispGrid): after finding the best a_hat on a 20-point grid
//   disp_grid_fine = arma::linspace<arma::vec>(a_hat - delta, a_hat + delta, disp_grid_n);  // delta = grid spacing
```

The $\ell_{CR}$ and penalty forms of sections 6.2 and 7.4 of the study text and the backtracking line search of 6.4 correspond directly. But `ll_part` lacks the terms that do not depend on $\alpha$, $-\log\Gamma(K_j+1)$ and $K_j\log\mu_j$. So the value of the C++ objective differs from the study text's $\ell_{CR}$ by a constant, but the maximum is the same (a difference of −919.43 for gene A; see "Solving gene A again with DESeq2's optimizer" below). The tolerance of the `noIncrease` test, `abs(initial_lp)/1e6`, is also based on this C++ value. In R this function can be called directly as `DESeq2:::fitDispWrapper()`.

The $w_j=1/(1/\hat\mu_j+\alpha)$ of `w_diag` is the Fisher information weight of a log-link NB-GLM, and as $\alpha\to0$ it goes back to the Poisson $\mu_j$.

</details>

<details>
<summary>Reading the DESeq2 source (3): fitting the trend, and varLogDispEsts</summary>

```r
args(DESeq2:::estimateDispersionsFit)
#> function (object, fitType = c("parametric", "local", "mean", 
#>     "glmGamPoi"), minDisp = 1e-08, quiet = FALSE) 
#> NULL
```

```r
# excerpt (not run): body of parametricDispersionFit. Iteratively fits a Gamma GLM (identity link) as disps ~ 1/means
fit <- glm(disps[good] ~ I(1/means[good]), family = Gamma(link = "identity"), start = coefs)
names(coefs) <- c("asymptDisp", "extraPois");  ans <- function(q) coefs[1] + coefs[2]/q
# on failure: "-- note: fitType='parametric', but ... a local regression fit was automatically substituted."
```

The moment the trend function is set with `dispersionFunction(dds) <-`, the setter creates the `dispFit` column and the residual spread. `varLogDispEsts` (a MAD-based robust variance), used in the outlier test, comes from here. An excerpt of the setter body:

```r
# excerpt (not run): body of getMethod("dispersionFunction<-", c("DESeqDataSet", "function"))
aboveMinDisp <- dispGeneEst >= minDisp * 100                 # gene-wise values below 1e-6 (near the lower limit) are left out of the spread
dispResiduals <- log(dispGeneEst) - log(dispFit)
varLogDispEsts <- mad(dispResiduals[aboveMinDisp], na.rm = TRUE)^2
attr(value, "varLogDispEsts") <- varLogDispEsts
```

</details>

<details>
<summary>Reading the DESeq2 source (4): the prior width, the MAP and the outlier rule</summary>

```r
args(DESeq2:::estimateDispersionsMAP)
#> function (object, outlierSD = 2, dispPriorVar, minDisp = 1e-08, 
#>     kappa_0 = 1, dispTol = 1e-06, maxit = 100, useCR = TRUE, 
#>     weightThreshold = 0.01, modelMatrix = NULL, type = c("DESeq2", 
#>         "glmGamPoi"), quiet = FALSE) 
#> NULL
```

```r
# excerpt (not run): body of estimateDispersionsPriorVar, three branches (m in the source = number of samples, p = ncol(X))
if (((m - p) <= 3) & (m > p)) {         # 1 <= m-p <= 3: minimize the KL divergence to a simulated log-chi^2 + normal noise distribution
    ...; dispPriorVar <- pmax(argminKL, 0.25); return(dispPriorVar) }
if (m > p) {                            # m-p > 3: observed spread - trigamma noise, lower limit 0.25
    expVarLogDisp <- trigamma((m - p)/2)
    dispPriorVar <- pmax((varLogDispEsts - expVarLogDisp), 0.25)
} else {                                # m == p: no noise subtraction or lower limit. Not reached in the standard workflow (checked below)
    dispPriorVar <- varLogDispEsts; expVarLogDisp <- 0 }
# --- estimateDispersionsMAP excerpt ---
dispInit <- ifelse(dispGeneEst > 0.1 * dispFit, dispGeneEst, dispFit)
dispInit[is.na(dispInit)] <- mcols(objectNZ)$dispFit[is.na(dispInit)]
dispResMAP <- fitDispWrapper(..., mu_hatSEXP = mu, log_alphaSEXP = log(dispInit),
    log_alpha_prior_meanSEXP = log(mcols(objectNZ)$dispFit),
    log_alpha_prior_sigmasqSEXP = log_alpha_prior_sigmasq, ..., usePriorSEXP = TRUE, ...)
dispConv <- dispResMAP$iter < maxit;  refitDisp <- !dispConv      # convergence is only iter < maxit (no iter == 1 or 10*minDisp condition as in gene-wise)
if (sum(refitDisp) > 0) { dispGrid <- fitDispGridWrapper(..., usePrior = TRUE, ..., useCRSEXP = TRUE); dispMAP[refitDisp] <- dispGrid }  # regardless of the useCR argument
# outlier exception rule
dispOutlier <- log(dispGeneEst) > log(dispFit) + outlierSD * sqrt(varLogDispEsts)
dispOutlier[is.na(dispOutlier)] <- FALSE
dispersionFinal <- dispMAP;  dispersionFinal[dispOutlier] <- dispGeneEst[dispOutlier]
```

The `m == p` branch of the source is never reached in the standard workflow. Running a design in which every sample has its own level (number of samples = number of coefficients) makes both functions stop before the prior variance is computed. So this branch is reached only when `estimateDispersionsMAP` and the like are called directly on input that has not gone through the gene-wise step; it is effectively dead code.

```r
dp <- makeExampleDESeqDataSet(n = 100, m = 6); dp$g <- factor(1:6); design(dp) <- ~g; dp <- estimateSizeFactors(dp)
for (f in list(estimateDispersions, estimateDispersionsGeneEst)) cat(tryCatch({ f(dp, quiet = TRUE); "ran" }, error = conditionMessage), "\n")
#> 
#> 
#>   The design matrix has the same number of samples and coefficients to fit,
#>   so estimation of dispersion is not possible. Treating samples
#>   as replicates was deprecated in v1.20 and no longer supported since v1.22.
#> 
#>  
#> the number of samples and the number of model coefficients are equal,
#>   i.e., there are no replicates to estimate the dispersion.
#>   use an alternate design formula
```

A few more details. The `mu` of the MAP step is `assays(objectNZ)[["mu"]]`, the $\hat\mu^0$ stored in the gene-wise step, reused as it is. The SD in the outlier criterion is not the prior SD $\sigma_d$ but `sqrt(varLogDispEsts)` (the MAD-based SD of the residuals), and the multiplier is `outlierSD = 2`. And `dispersions(dds)` is just an accessor that returns `mcols(dds)$dispersion`.

```r
getMethod("dispersions", "DESeqDataSet")
#> Method Definition:
#> 
#> function (object, ...) 
#> {
#>     .local <- function (object) 
#>     mcols(object)$dispersion
#>     .local(object, ...)
#> }
#> <bytecode: 0x5bc059769628>
#> <environment: namespace:DESeq2>
#> 
#> Signatures:
#>         object        
#> target  "DESeqDataSet"
#> defined "DESeqDataSet"
```

I also looked at the column names and descriptions of the example data, and at which column the final values came from.

```r
cat(names(df), fill = TRUE)
#> trueIntercept trueBeta trueDisp baseMean baseVar allZero dispGeneEst 
#> dispGeneIter dispFit dispersion dispIter dispOutlier dispMAP
as.data.frame(mcols(mcols(dds)))[c("dispGeneEst", "dispFit", "dispMAP", "dispOutlier", "dispersion"), ]
#>                     type                       description
#> dispGeneEst intermediate gene-wise estimates of dispersion
#> dispFit     intermediate       fitted values of dispersion
#> dispMAP     intermediate     maximum a posteriori estimate
#> dispOutlier intermediate     dispersion flagged as outlier
#> dispersion  intermediate      final estimate of dispersion
with(df[ok, ], c(nonoutlier_eq_MAP = all(dispersion[!dispOutlier] == dispMAP[!dispOutlier]),
                 outlier_eq_geneEst = all(dispersion[dispOutlier] == dispGeneEst[dispOutlier])))
#>  nonoutlier_eq_MAP outlier_eq_geneEst 
#>               TRUE               TRUE
```

</details>

<details>
<summary>How is the prior width σ_d² computed?</summary>

First, the symbols: $n$ is the number of samples (`m` in the DESeq2 source), $p$ the number of columns of the design matrix, and $n-p$ the residual degrees of freedom.

$$\sigma_d^2=\max\!\Big(\big(1.4826\cdot\mathrm{MAD}\big)^2\big[\log\hat\alpha_{gw}-\log\alpha_{tr}\big]-\psi_1\!\big(\tfrac{n-p}{2}\big),\;0.25\Big)\qquad(n-p>3;\ \text{only genes with }\hat\alpha_{gw}\ge 100\cdot\texttt{minDisp})$$

- $\psi_1$ is the trigamma function. For the sample variance $s^2$ of normal residuals with $\nu$ degrees of freedom, $\mathrm{Var}(\log s^2)=\psi_1(\nu/2)$. DESeq2 borrows this as an approximation to the estimation noise of the log dispersion: it subtracts this estimation noise from the observed spread of the log-dispersion residuals and takes what is left as the real spread between genes.
- The spread is measured with R's `mad()^2`, which is robust, instead of the variance. `mad()` already multiplies by the constant 1.4826 that converts to a normal SD, so it is about 2.2 times the square of the raw $\mathrm{median}|x-\mathrm{median}(x)|$.
- If $1\le n-p\le3$, it takes a separate (KL divergence) path. I actually ran it further down in this block. The source also has an $n=p$ branch, but the standard workflow stops with an error before reaching it ("Reading the DESeq2 source (4)").
- In the example, $n-p=4$, so subtracting trigamma(2) = 0.645 from the observed residual spread of 0.679 left only 0.034, and the lower limit 0.25 applied. This is not because of the number of samples but because of how the example data were constructed. `makeExampleDESeqDataSet` places the true dispersions exactly on the curve `dispMeanRel = function(x) 4/x + 0.1` (the `all.equal` below is TRUE). So the true spread between genes is 0, and the observed spread of 0.679 is almost entirely estimation noise.

```r
all.equal(df$trueDisp, 4/2^df$trueIntercept + 0.1)
#> [1] TRUE
round(c(varLogDispEsts = vld, trigamma2 = trigamma(2), dispPriorVar = pv), 5)
#> varLogDispEsts      trigamma2   dispPriorVar 
#>        0.67906        0.64493        0.25000
```

With 3 or fewer degrees of freedom the KL path is taken, but the examples in this note all have $n-p=4$, so they never pass through this branch. So I made data with 4 samples and ran `~condition` ($n-p=2$) and `~1` ($n-p=3$) separately. `klByHand()` is the branch body of the source copied by hand, and it uses the same `set.seed(2)` as the source. I compared the case where the true log dispersions lie exactly on the trend curve (`sd=0`, the same construction as this note's example) with the case where they scatter around the curve with SD 1 (`sd=1`). This block was run in a fresh R session.

```r
suppressMessages(library(DESeq2))
# reproduce the 1 <= m-p <= 3 branch of the source by hand. Also returns argminKL before pmax(, 0.25)
klByHand <- function(res, df) {
  set.seed(2); brks <- -20:20/2
  oh <- hist(res[res > min(brks) & res < max(brks)], breaks = brks, plot = FALSE)$density
  grid <- seq(0, 8, length = 200)
  kl <- sapply(grid, function(x) {
    r <- log(rchisq(10000, df = df)) + rnorm(10000, 0, sqrt(x)) - log(df)
    rh <- hist(r[r > min(brks) & r < max(brks)], breaks = brks, plot = FALSE)$density
    small <- min(c(oh, rh)[c(oh, rh) > 0])
    sum(oh * (log(oh + small) - log(rh + small))) })
  fine <- seq(0, 8, length = 1000)
  fine[which.min(predict(loess(kl ~ grid, span = 0.2), fine))]
}
run <- function(label, sdLog, d) {
  set.seed(1)  # sdLog: SD of the true log dispersion around the trend (0 = on the curve)
  dds <- makeExampleDESeqDataSet(m = 4, n = 2000,
           dispMeanRel = function(x) (4/x + 0.1) * exp(rnorm(length(x), 0, sdLog)))
  design(dds) <- as.formula(d)
  dds <- suppressMessages(estimateDispersions(estimateSizeFactors(dds)))
  mc <- mcols(dds)[!mcols(dds)$allZero, ]; ok <- mc$dispGeneEst >= 100 * 1e-8
  df <- nrow(colData(dds)) - ncol(model.matrix(as.formula(d), colData(dds)))
  f <- dispersionFunction(dds); v <- attr(f, "varLogDispEsts")
  kl <- klByHand(log(mc$dispGeneEst[ok]) - log(mc$dispFit[ok]), df)
  cat(sprintf("%-6s %-10s n-p=%d var=%.3f trigamma=%.3f | trigamma formula %.3f | argminKL %.3f -> pmax %.3f | DESeq2 %.3f\n",
    label, d, df, v, trigamma(df/2), max(v - trigamma(df/2), 0.25), kl, max(kl, 0.25),
    attr(f, "dispPriorVar")))
}
for (d in c("~condition", "~1")) for (s in c(0, 1)) run(paste0("sd=", s), s, d)
```

```
sd=0   ~condition n-p=2 var=0.802 trigamma=1.645 | trigamma formula 0.250 | argminKL 0.000 -> pmax 0.250 | DESeq2 0.250
sd=1   ~condition n-p=2 var=1.468 trigamma=1.645 | trigamma formula 0.250 | argminKL 0.488 -> pmax 0.488 | DESeq2 0.488
sd=0   ~1         n-p=3 var=0.868 trigamma=0.935 | trigamma formula 0.250 | argminKL 0.000 -> pmax 0.250 | DESeq2 0.250
sd=1   ~1         n-p=3 var=1.748 trigamma=0.935 | trigamma formula 0.813 | argminKL 0.809 -> pmax 0.809 | DESeq2 0.809
```

In all four cases the hand-copied calculation (the `pmax` column) equals DESeq2's `dispPriorVar`. With `sd=0`, `argminKL` is 0 so the lower limit 0.25 is used, and the trigamma formula is also 0.25 then, so the two paths cannot be told apart. The two formulas diverge at `sd=1` with $n-p=2$. The observed variance of 1.468 is smaller than $\psi_1(1)=1.645$, so the trigamma formula gives 0.25, but the KL path gives 0.488, and DESeq2 uses 0.488. At $n-p=3$ the two values are similar, 0.813 and 0.809, and DESeq2 uses the KL value, 0.809.

</details>

<details>
<summary>Edge cases of the implementation: lower and upper limits, noIncrease, grid re-optimization</summary>

- The range searched for $\theta$ is not "the unrestricted real line". The C++ `fitDisp` clips proposals to $\log\alpha\in[-30,10]$ and stops if $\log\alpha$ falls below `log(minDisp/10)`. The R side also clips the result again to `[minDisp, max(10, number of samples)]`. This is why A6.5 hits 1e-8 under `~condition`.
- The initial value is `pmin(roughDisp, momentsDisp)` clipped to `[minDisp, maxDisp]`. For gene A, rough 0.0232 < moment 0.149, so the actual initial value is the rough one.
- For genes where the optimization fails to raise the objective above the initial value (`noIncrease`: `last_lp < initial_lp + abs(initial_lp)/1e6`), that initial value becomes `dispGeneEst` as it is. So "moment estimation is used only for the initial value" is slightly too strong.
- In the gene-wise step, if `dispIter == maxit` or `dispIter == 1` and `dispGeneEst > 10·minDisp`, it is recomputed with `fitDispGridWrapper` (a 20-point grid on the R side, and the C++ `fitDispGrid` sweeps a fine grid around the optimum once more). Genes at the lower limit (A6.5) are not grid targets.
- The grid fallback condition of the MAP step is only `iter == maxit` (no `iter == 1` or `minDisp` condition). The grid call is fixed at `useCRSEXP = TRUE` regardless of the `useCR` argument.
- `niter` is the number of outer mean→dispersion cycles (default 1), and `maxit` the upper limit of the inner C++ line-search iterations (default 100). Raising `maxit` does not update the mean more often.
- The initial mean is computed as a linear mean (`linearModelMuNormalized`) for a group design whose number of unique rows equals the number of columns, and with an NB-GLM (`fitNbinomGLMs`) otherwise. With sample weights it is always the GLM.
- The upper limit `maxDisp` is `max(10, number of samples)`. In the example, 10 genes below the trend sat at the upper limit 10 in both the gene-wise and MAP values and did not move.

</details>

<details>
<summary>A little more on the likelihood: the independence assumption, the NB probability formula, the profile likelihood</summary>

Building $L$ by multiplying per-sample probabilities rests on the model assumption that the samples are independent conditional on the design. Pair and nested designs use this product only after absorbing that dependence into design terms.

With $r=1/\alpha$, the NB probability formula is the following, the same as R's `dnbinom(size = r, mu = mu)`.

$$f(k;\mu,\alpha)=\frac{\Gamma(k+r)}{\Gamma(r)\,\Gamma(k+1)}\left(\frac{r}{r+\mu}\right)^{r}\left(\frac{\mu}{r+\mu}\right)^{k}$$

$$\ell_j=\log\Gamma(K_j+r)-\log\Gamma(r)-\log\Gamma(K_j+1)+r\{\log r-\log(r+\mu_j)\}+K_j\{\log\mu_j-\log(r+\mu_j)\}$$

The gamma function extends the factorial to a real $r$. Computing the lgamma expansion and the `dnbinom` sum side by side at the Cox-Reid value $\alpha$ = 0.025385 shows that they agree.

```r
r <- 1/exp(a_cr)
c(lgamma = sum(lgamma(K+r) - lgamma(r) - lgamma(K+1) + r*(log(r)-log(r+mu)) + K*(log(mu)-log(r+mu))), dnbinom = ll_nb(exp(a_cr)))
#>    lgamma   dnbinom 
#> -27.24031 -27.24031
```

The strict profile likelihood is $\ell_{profile}(\alpha)=\ell(\hat b(\alpha),\alpha)$, the curve built by refitting $\hat b(\alpha)=\arg\max_b \ell(b,\alpha)$ for each candidate $\alpha$. DESeq2's default implementation fixes the fitted mean $\hat\mu^0$ obtained with an initial $\alpha$ and optimizes only $\alpha$, with one outer cycle by default (`niter = 1`). So `dispGeneEst` must be distinguished both from the plain MLE and from the maximum of the strict profile likelihood (section 6.6 of the study text).

Section 5.5 of the study text writes the pooled mean with an accent over $K$. In this note I wrote it as $\bar K$, and the normalized mean written the same way in 7.2 as $\bar q$.

</details>

<details>
<summary>The difference between Cox-Reid and the n/(n−p) correction, and the direction of the correction</summary>

As with Bessel's correction in the normal model, in regression the residual degrees of freedom shrink by the number of parameters $p$ used for the mean structure. The Cox-Reid adjustment is a way of handling this problem for the NB, but it is not the same as a formula that multiplies the variance by $n/(n-p)$. The adjustment term $-\tfrac12\log|X^\top W(\alpha)X|$ is a function that depends on $\alpha$ through $W(\alpha)$, and on $X$ and $\hat\mu$, not a fixed multiplier.

When comparing in numbers, don't miss that $n/(n-p)$ is a factor that multiplies the **variance**. Since $\alpha=(\mathrm{Var}-\mu)/\mu^2$, multiplying the variance by 1.5 makes $\alpha$, with the Poisson part removed, grow by more than 1.5 times.

```r
mom <- function(s, k = K, m = mu) sum((s*(k - m)^2 - m)/m^2)/length(k)   # moment α: multiply residuals² (variance) by s, then subtract the Poisson share
round(c(mom_n = mom(1), mom_var_x_6_4 = mom(6/4), alpha_ratio = mom(6/4)/mom(1), CR_over_MLE = exp(a_cr)/exp(a_mle)), 5)
#>         mom_n mom_var_x_6_4   alpha_ratio   CR_over_MLE 
#>       0.01545       0.02671       1.72871       1.71682
ag2 <- 10^seq(-4, 1, length.out = 200)
c(CR_term_increasing = all(diff(sapply(ag2, cr)) > 0))                     # -0.5 log|X'W(a)X| is increasing in a
#> CR_term_increasing 
#>               TRUE
```

For gene A, the moment estimate with the residuals² multiplied by 6/4 enlarges $\alpha$ 1.73 times, and Cox-Reid 1.72 times, so the two corrections are almost the same size. The numbers alone cannot tell them apart, and the point of section 6.1 of the study text is about structure. For B6.5, with the same `~condition` design, Cox-Reid gives 1.49 times and the variance correction 1.53 times, ratios different from gene A's (block below).

The direction is fixed mathematically. With $\hat\mu$ fixed, $\partial w_j/\partial\alpha=-\hat\mu_j^2/(1+\alpha\hat\mu_j)^2<0$, so $X^\top W(\alpha)X$ shrinks and its determinant shrinks, which makes the Cox-Reid term $-\tfrac12\log|X^\top WX|$ an increasing function of $\alpha$ (`CR_term_increasing` above). The maximum of an objective with an increasing function added cannot be smaller than the original maximum, so with the mean fixed, Cox-Reid makes $\alpha$ equal or larger (equal when it hits the lower limit, as with A6.5).

I solved the two genes of section 6.5 of the study text again with the function I implemented, and also printed the mean path DESeq2 used and the mean values.

```r
for (dsg in list(~condition, ~1)) {
  dd <- DESeqDataSetFromMatrix(cts, cd, dsg); sizeFactors(dd) <- rep(1, 6)
  gg <- estimateDispersionsGeneEst(dd, quiet = TRUE); mm <- model.matrix(dsg, cd)
  cat(deparse(dsg), ": ncol(X) =", ncol(mm), " linearMu =", nlevels(DESeq2:::modelMatrixGroups(mm)) == ncol(mm),
      " mu(A6.5) =", assays(gg)[["mu"]]["A6.5", ], "\n")
}
#> ~condition : ncol(X) = 2  linearMu = TRUE  mu(A6.5) = 100 100 100 200 200 200 
#> ~1 : ncol(X) = 1  linearMu = TRUE  mu(A6.5) = 150 150 150 150 150 150
own2 <- function(K, X) { mu <- as.vector(X %*% solve(crossprod(X), crossprod(X, K)))
  f <- function(th) { a <- exp(th); W <- diag(mu/(1+a*mu)); -(sum(dnbinom(K, mu=mu, size=1/a, log=TRUE)) - 0.5*log(det(t(X)%*%W%*%X))) }
  exp(optimize(f, c(log(1e-8), log(10)), tol = 1e-10)$minimum) }
Xc <- cbind(1, rep(0:1, each = 3)); Xi <- matrix(1, 6, 1)
signif(rbind(A6.5 = c(condition = own2(cts["A6.5", ], Xc), intercept = own2(cts["A6.5", ], Xi)),
             B6.5 = c(own2(cts["B6.5", ], Xc), own2(cts["B6.5", ], Xi))), 7)
#>         condition intercept
#> A6.5 1.000011e-08 0.1321920
#> B6.5 2.199911e-01 0.3407456
mle2 <- function(K, X) { mu <- as.vector(X %*% solve(crossprod(X), crossprod(X, K)))
  exp(optimize(function(th) -sum(dnbinom(K, mu=mu, size=1/exp(th), log=TRUE)), c(log(1e-8), log(10)), tol = 1e-10)$minimum) }
kB <- cts["B6.5", ]; muB <- as.vector(Xc %*% solve(crossprod(Xc), crossprod(Xc, kB)))
round(c(B_CR_over_MLE = own2(kB, Xc)/mle2(kB, Xc), B_var_x_6_4 = mom(6/4, kB, muB)/mom(1, kB, muB)), 3)   # gene A: 1.717 / 1.729
#> B_CR_over_MLE   B_var_x_6_4 
#>         1.486         1.527
```

The study text's claim that "dropping the condition changes the dispersion a lot" is exactly right, with the addition that for A6.5 the "small α" is actually the lower limit `minDisp = 1e-8`.

</details>

<details>
<summary>Solving gene A again with DESeq2's optimizer</summary>

Feeding the same data directly into DESeq2's C++ optimizer (`fitDispWrapper`) also gives the study text's numbers.

```r
Y <- matrix(K, nrow = 1); MU <- matrix(mu, nrow = 1); Wt <- matrix(1, 1, 6)
run <- function(useCR, usePrior = FALSE, pm = 0, ps = 1, init = 0.1)
  DESeq2:::fitDispWrapper(ySEXP = Y, xSEXP = X, mu_hatSEXP = MU, log_alphaSEXP = log(init),
     log_alpha_prior_meanSEXP = pm, log_alpha_prior_sigmasqSEXP = ps, min_log_alphaSEXP = log(1e-9),
     kappa_0SEXP = 1, tolSEXP = 1e-6, maxitSEXP = 100, usePriorSEXP = usePrior,
     weightsSEXP = Wt, useWeightsSEXP = FALSE, weightThresholdSEXP = 0.01, useCRSEXP = useCR)
for (r in list(run(FALSE), run(TRUE), run(TRUE, TRUE, log(0.08), 0.49))) cat(sprintf("alpha = %.8f  iter = %d\n", exp(r$log_alpha), r$iter))
#> alpha = 0.01478473  iter = 8
#> alpha = 0.02538751  iter = 5
#> alpha = 0.05314358  iter = 8
rc <- run(TRUE); a1 <- exp(rc$log_alpha)   # C++ objective vs study text l_CR: differ only by a constant independent of alpha
c(cpp_last_lp = rc$last_lp, textbook_ll_cr = ll_cr(a1), diff = rc$last_lp - ll_cr(a1), const = sum(lgamma(K+1)) - sum(K*log(mu)))
#>    cpp_last_lp textbook_ll_cr           diff          const 
#>     -951.20376      -31.76941     -919.43436     -919.43436
c(rough = DESeq2:::roughDispEstimate(counts(g1, normalized = TRUE), X), moment = DESeq2:::momentsDispEstimate(DESeq2:::getBaseMeansAndVariances(g1)))
#>      rough     moment 
#> 0.02317952 0.14911911
```

The initial value is `pmin(rough, moment)` = rough 0.0232 (not moment 0.149).

Running the whole of `estimateDispersions()` on this 1-gene object gives a result unrelated to the prior (0.08, 0.7) of section 16.3 of the study text. The parametric fit fails and is replaced by local; with only one point, `dispFit == dispGeneEst` and `varLogDispEsts = 0`, and `dispPriorVar` becomes the lower limit 0.25.

```r
e1 <- suppressWarnings(estimateDispersions(dds1, quiet = TRUE))   # warning: "Estimated rdf < 1.0; not estimating variance"
#> -- note: fitType='parametric', but the dispersion trend was not well captured by the
#>    function: y = a/x + b, and a local regression fit was automatically substituted.
#>    specify fitType='local' or 'mean' to avoid this message next time.
as.data.frame(mcols(e1))[, c("dispGeneEst", "dispFit", "dispMAP", "dispersion")]
#>   dispGeneEst    dispFit    dispMAP dispersion
#> 1  0.02538756 0.02538756 0.02538672 0.02538672
c(fitType = attr(dispersionFunction(e1), "fitType"), varLogDispEsts = attr(dispersionFunction(e1), "varLogDispEsts"),
  dispPriorVar = attr(dispersionFunction(e1), "dispPriorVar"))
#>        fitType varLogDispEsts   dispPriorVar 
#>        "local"            "0"         "0.25"
```

| Objective | Optimal α (recomputed) | Study text | DESeq2 function | Interpretation |
|---|---|---|---|---|
| $\ell$ (NB log-likelihood) | 0.014786 | 0.014786 | `useCR = FALSE` → 0.014785 | plain dispersion MLE (mean fixed) |
| $\ell_{CR}$ | 0.025385 | 0.025385 | `useCR = TRUE` → 0.025388 | adjusted for the bias from estimating the mean, `dispGeneEst` |
| $\ell_{CR}$ − penalty, $N(\log 0.08, 0.7^2)$ | 0.053147 | 0.053147 | `usePrior = TRUE` → 0.053144 | teaching MAP, between gene-wise 0.025 and the prior centre 0.08 (the prior is given) |
| pooled moment $(s^2-\bar K)/\bar K^2$ | 0.149 | — | candidate initial value (loses to rough 0.0232 here) | inflated because the difference between the two group means is mixed into the variance. Within-group moments are 0.029 (Ctrl) and 0.025 (Starvation) |

Figure 3 of the study text (Figure 1 of this note) shifts each of the three objectives so that its maximum is 0, so the absolute heights between the curves are not meant to be compared. It is the same reason that the C++ objective value and the study text's $\ell_{CR}$ differ by −919.43 yet have the same maximum. The calculation that turns figure 3 into a grid table is in [08](08_one_gene_end_to_end.md).

</details>

<details>
<summary>The scale of the MAP (Jacobian)</summary>

Even for the same lognormal prior, turning it into a density on $\alpha$ adds the Jacobian term $-\log\alpha$, which changes the mode.

```r
post_alpha <- function(a, m = log(0.08), s2 = 0.49) ll_cr(a) - (log(a) - m)^2/(2*s2) - log(a)   # the same prior as a density on alpha
c(mode_theta = exp(a_map), mode_alpha = optimize(function(a) -post_alpha(a), c(1e-6, 5), tol = 1e-12)$minimum)
#> mode_theta mode_alpha 
#> 0.05314735 0.03843084
c(ex9_weight_gw = log(0.05/0.1)/log(0.02/0.1))   # Problem 9: gene-wise weight that corresponds to 0.05 on the theta scale
#> ex9_weight_gw 
#>     0.4306766
```

The mode on the $\alpha$ scale is 0.038, and on the $\theta$ scale 0.053. DESeq2 optimizes in $\theta$ space, so 0.053 is the value in DESeq2's manner. The MAP depends on the parameterization, so don't assume the mode is the same on every scale (section 7.4 of the study text).

</details>

<details>
<summary>Computing dispGeneEst and dispMAP directly to match them</summary>

Taking only `assays(dds)[["mu"]]`, `dispFit` and `dispPriorVar` from the example object and maximizing $\ell_{CR}$ and the MAP objective directly with `optimize()` gives the same values as DESeq2's columns.

```r
KK <- counts(dds); MUe <- assays(dds)[["mu"]]; mm <- model.matrix(design(dds), colData(dds))
obj <- function(i, a, prior = FALSE) {
  m <- MUe[i, ]; W <- diag(m/(1 + a*m))
  v <- sum(dnbinom(KK[i, ], mu = m, size = 1/a, log = TRUE)) - 0.5*log(det(t(mm) %*% W %*% mm))
  if (prior) v <- v - (log(a) - log(df$dispFit[i]))^2/(2*pv)
  v
}
own <- function(i, prior) exp(optimize(function(th) -obj(i, exp(th), prior), c(log(1e-8), log(10)), tol = 1e-10)$minimum)
cmp <- t(sapply(pick, function(g) { i <- match(g, rownames(df))
  c(geneEst = df$dispGeneEst[i], own_CR = own(i, FALSE), MAP = df$dispMAP[i], own_MAP = own(i, TRUE)) }))
round(cmp, 5)
#>         geneEst  own_CR     MAP own_MAP
#> gene733 3.93338 3.93333 1.08290 1.08294
#> gene798 7.39634 7.39793 1.93150 1.93160
#> gene12  1.25983 1.25984 0.55373 0.55377
#> gene15  0.31192 0.31183 0.19847 0.19848
#> gene1   1.23019 1.23203 2.93373 2.93379
#> gene2   0.42131 0.42136 0.44639 0.44639
signif(c(max_relerr_CR = max(abs(cmp[, "own_CR"]/cmp[, "geneEst"] - 1)), max_relerr_MAP = max(abs(cmp[, "own_MAP"]/cmp[, "MAP"] - 1))), 3)
#>  max_relerr_CR max_relerr_MAP 
#>       1.50e-03       7.46e-05
```

The difference is at most a relative error of 0.15% (the gene-wise value of gene1). The C++ convergence criterion is an objective change `< dispTol = 1e-6`, so on a flat curve $\alpha$ stops with correspondingly less precision. In other words, the study text's $\ell_{CR}$ and MAP formulas differ from the functions the installed DESeq2 actually maximizes only by constants that do not depend on $\alpha$, and the maxima are the same.

</details>

<details>
<summary>The 3 genes whose MAP fell outside the gene-wise value and the trend</summary>

```r
nz <- df[ok & !df$dispOutlier, ]; lo <- pmin(nz$dispGeneEst, nz$dispFit); hi <- pmax(nz$dispGeneEst, nz$dispFit)
with(nz, c(between = sum(dispMAP >= lo & dispMAP <= hi), outside = sum(dispMAP < lo | dispMAP > hi), below_trend = sum(dispGeneEst < dispFit),
           below_up = sum(dispGeneEst < dispFit & dispMAP > dispGeneEst), below_down = sum(dispGeneEst < dispFit & dispMAP < dispGeneEst),
           below_same = sum(dispGeneEst < dispFit & dispMAP == dispGeneEst), below_same_at_10 = sum(dispGeneEst < dispFit & dispMAP == dispGeneEst & dispMAP == 10)))
#>          between          outside      below_trend         below_up 
#>             1966                3             1336             1326 
#>       below_down       below_same below_same_at_10 
#>                0               10               10
out3 <- rownames(nz)[nz$dispMAP < lo | nz$dispMAP > hi]
signif(cbind(df[out3, c("dispGeneEst", "dispGeneIter", "dispFit", "dispMAP")], own_CR = sapply(match(out3, rownames(df)), own, prior = FALSE)), 6)
#>         dispGeneEst dispGeneIter  dispFit  dispMAP   own_CR
#> gene4    0.12105700            7 0.122465 0.122834 0.123550
#> gene274  0.11263100            8 0.112523 0.112522 0.112628
#> gene367  0.00000001            1 0.397029 0.443694 0.549531
DESeq2:::roughDispEstimate(counts(dds, normalized = TRUE)["gene367", , drop = FALSE], mm)   # candidate initial value for gene367
#> gene367 
#>       0
```

For 1966 of the 1969 non-outliers, `dispMAP` lies between `dispGeneEst` and `dispFit`. Of the 1336 below the trend, none went down, 1326 went up, and 10 sat at the upper limit `maxDisp = 10` in both the gene-wise and MAP values.

For the 3 that fell outside, the stored `dispGeneEst` is not the exact maximum of $\ell_{CR}$.

- gene367 has `roughDisp = 0`, so its initial value became the lower limit 1e-8. $\ell_{CR}$ is flat around there, so the C++ code stopped after 1 iteration (`dispGeneIter = 1`; the slope is almost 0, so the change in the objective at the first step was below `dispTol`), and 1e-8 also fails the grid refit condition `> 10·minDisp`, so it stayed. The actual maximum of $\ell_{CR}$ is 0.55, and the MAP of 0.444 lies between 0.55 and the trend 0.397.
- gene4 is a case where `dispGeneEst` stopped 2% below the actual maximum, 0.1236. The MAP of 0.1228 lies between 0.1236 and the trend 0.1225.
- gene274 is a case at the level of numerical error, where the MAP and the trend agree within 1e-5.

So the "between" picture of section 7.5 of the study text is right relative to the actual maximum of $\ell_{CR}$, but relative to the stored `dispGeneEst` column there are exceptions.

</details>

<details>
<summary>Common misconceptions</summary>

| Misconception | In fact |
|---|---|
| "The dispersion MLE is just the sample variance minus the Poisson variance" | That is a moment approximation that holds only with a single mean and size factor 1. For gene A (chapter 16), the pooled moment is 0.149, the NB MLE with the condition included and the mean fixed is 0.0148, and `dispGeneEst`, which adds Cox-Reid on top, is 0.0254. `dispGeneEst`, often called the "dispersion MLE" in DESeq2, is also not a plain MLE but this CR-adjusted estimate. DESeq2 uses the moment estimate as a candidate initial value, and the initial value remains only for genes where the optimization fails to improve on it (`noIncrease`) |
| "DESeq2 refits β for every candidate α (profile likelihood)" | The default is `niter = 1`. It computes $\hat\mu^0$ once with an initial α and fixes it. The MAP step reuses the same `assays[["mu"]]` |
| "Raising `maxit` updates the mean more often" | `maxit` is the upper limit of the C++ line-search iterations; the number of mean updates is `niter` |
| "Cox-Reid is a kind of shrinkage too" | Cox-Reid deals with the bias from estimating the mean within one gene, and with the mean fixed it always makes α equal or larger (the CR term is an increasing function of α). Shrinkage is a prior across genes. In chapter 16 the former is 0.015→0.025 and the latter 0.025→0.053 |
| "Shrinkage makes the dispersion smaller" | The target is the trend. In the example data, 67% of the non-outliers (1326/1969) had a MAP larger than their gene-wise value |
| "Higher means should have lower α / if the trend goes down, the count variance goes down too" | The trend is a summary of a tendency, not a law. The variance is $\mu+\alpha\mu^2$, so in the example fit, as $q$ goes 10→1000, $\alpha_{tr}$ falls 0.714→0.100 while the variance rises 81→about 101,000 |
| "`dispersions(dds)` is always `dispMAP`" | For outlier genes (8 in the example) it is `dispGeneEst`. You have to look at the `dispOutlier` column |
| "The outlier criterion is 2 times the prior SD" | It is `2 * sqrt(varLogDispEsts)` (the MAD-based SD of the residuals). In the example, using the prior SD would have caught 79 genes instead of 8 |
| "Dispersion outliers and Cook's outliers are the same thing" | The former is a judgement on a gene's dispersion, the latter on how much one sample's count influences the coefficient fit ([07](07_lfc_shrinkage_and_qc.md)) |
| "The MAP has the same mode on every scale" | 0.0531 in $\theta$ space vs 0.0384 in $\alpha$ space. DESeq2 uses $\theta$ |

</details>

<details>
<summary>Where the results differ from the study text</summary>

I kept only the claims where the study text's explanation differs a little from the actual behaviour of DESeq2, or where there are details beyond what the study text gives.

| Study text claim | Result | Evidence |
|---|---|---|
| 5.1 The likelihood is the product of the per-sample pmfs (assuming independence conditional on the design) | Matches (conceptual) | The sum of `dnbinom(..., log=TRUE)` = $\log\prod$; DESeq2 also adds per-sample terms (the C++ `ll_part` is a sum over samples. But $-\log\Gamma(K+1)+K\log\mu$, which does not depend on α, is left out, so the value differs by a constant: −919.43 in chapter 16) |
| 5.5 Moment estimation is used for the initial value | Matches (with an addition) | `alpha_hat <- pmin(roughDisp, momentsDisp)`; in chapter 16, pooled moment 0.149 vs NB MLE with the mean fixed 0.0148 (within-group moments 0.029/0.025, `dispGeneEst` with CR 0.0254). Addition: for `noIncrease` genes the initial value becomes `dispGeneEst` as it is |
| 6.1 Cox-Reid is not a formula that multiplies by $n/(n-p)$ | Matches (conceptual). Not distinguishable in numbers | Structure: the CR term is a function of α that depends on $W(\alpha)$ and $X$. Numbers: the moment with the variance (residuals²) multiplied by 6/4 enlarges α 1.729 times, CR 1.717 times (chapter 16), and for B6.5 it is 1.53 vs 1.49, similar sizes. Comparing the ratio of α with 1.5 is not appropriate (1.5 is a factor on the variance) |
| 6.4 Optimizing in log α (removes the positivity constraint, relative change), backtracking line search | Matches (with an addition; the C++ checked in the Bioconductor 1.50.2 tarball) | R side `log_alphaSEXP = log(alpha_hat)`; C++ `a_propose = a + kappa*dlp`, the Armijo condition `theta_hat_kappa`, `kappa = kappa/2.0`, `change < tol`. Addition: it is not "the unrestricted real line". The C++ code clips log α to [−30, 10] and stops when `a < min_log_alpha`, and R clips the result to `[minDisp, max(10, n)]` |
| 6.5 For A (A6.5 in this note), α is small under `~condition` and changes a lot under `~1` | Matches (with an addition) | A6.5: 1e-08 (= the `minDisp` lower limit) vs 0.132; B6.5: 0.220 vs 0.341. The study text does not mention that it hits the lower limit |
| 7.3 Prior variance ≈ observed residual variance − trigamma-based noise, robust spread, a lower limit, a separate path for small degrees of freedom | Matches (with an addition) | `mad(dispResiduals[aboveMinDisp])^2` (only `dispGeneEst >= 100*minDisp`), `trigamma((m-p)/2)`, `pmax(..., 0.25)`, the KL path if `1 <= m-p <= 3` (m in the source = number of samples). I ran the KL path on data with 4 samples ($n-p=2, 3$) and confirmed that the hand-copied calculation equals DESeq2's `dispPriorVar` in all four cases ("How is the prior width σ_d² computed?" block). Addition: the source has an `m == p` branch (no noise subtraction, no lower limit), but `estimateDispersions()` and `estimateDispersionsGeneEst()` `stop()` first, so it is not reached in the standard workflow. The lower limit 0.25 applied in the example because of the simulation's construction, in which the true dispersions lie exactly on the trend curve (true spread 0) |
| 7.5 A wide prior allows more differences between genes, a narrow one increases the influence of the trend / 7.3 the shrinkage proportion is not the same for every gene | Matches (a property of the formula) | Weights $1/v_i$ vs $1/\sigma_d^2$. In the 1-gene example of chapter 16 ($v_i=0.795$, $\sigma_d^2=0.49$), the gene-wise weight is 0.381 by the approximation and 0.356 from the actual MAP location (the weighted-average approximation block in section 5) |
| 7.6 Outlier genes use the gene-wise value instead of the MAP (forcing them down exaggerates significance) | Matches (the study text does not give the criterion) | The criterion is `log(dispGeneEst) > log(dispFit) + outlierSD(2) * sqrt(varLogDispEsts)`; not the prior SD; the exception is only on the upper side. `dispersion == dispGeneEst` for the 8 genes |
| 16.1 (group means 106.67/210, "not output from running the package") | Outside the scope of this note → checked in [08](08_one_gene_end_to_end.md) | This note uses them only as the $\hat\mu$ input (same as `assays(g1)[["mu"]]`) |
| 16.3 The prior (0.08, 0.7) is given, and there is no claim that the trend was estimated from 1 gene | Matches (with an addition) | Running `estimateDispersions()` on a 1-gene object: parametric fails → local substituted, `dispFit == dispGeneEst`, `varLogDispEsts = 0`, `dispPriorVar = 0.25` (lower limit), `dispMAP = 0.02539`. Values unrelated to the study text's prior |

The rest of the study text's explanations (5.2–5.3, 6.2–6.3, 6.6, the other claims of 7.1–7.6, and the numbers and figure 3 of 16.2–16.3) matched the actual behaviour of DESeq2.

The installed package has only the compiled `.so` for the C++ `fitDisp`/`fitDispGrid` bodies, so I downloaded `DESeq2_1.50.2.tar.gz` from Bioconductor 3.22 and compared them line by line in `src/DESeq2.cpp` (DESCRIPTION `Version: 1.50.2`). Every line I quoted is there as is, and the numbers obtained by implementing the same objective myself also match the output of the installed `fitDispWrapper`.

</details>

<details>
<summary>How the code was run</summary>

Blocks whose first line is `# excerpt (not run)` are source excerpts and were not run. The other R code blocks were run on R 4.5.2 and DESeq2 1.50.2, after calling `library(DESeq2)`, in one session in document order, through both the body and the collapsible blocks. So later blocks use objects from earlier blocks as they are (`K`, `mu`, `X`, `ll_nb`, `cr`, `ll_cr`, `a_mle`, `a_cr`, `a_map`, `cd`, `dds1`, `g1`, `cts`, `dds`, `fn`, `pv`, `vld`, `df`, `ok`, `pick`, `own`, `mm` and others). The `#>` lines are the output of that run through knitr (`collapse = TRUE`, `comment = "#>"`), copied without changes, and values that differ from run to run, such as `<bytecode: ...>` addresses, were left as they were. The figure code was run with the `notes/` folder as the working directory. The matching script in the roadmap is `03_dispersion_estimation.R`, and it was checked on 2026-09-26.

</details>

---

← Previous: [02. The negative binomial and size factors](02_negative_binomial.md) · Next: [04. The NB-GLM](04_glm_condition_batch.md) →

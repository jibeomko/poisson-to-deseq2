# 02. How much do counts spread, and how is the difference in sequencing depth corrected? (The negative binomial and size factors)

In [note 01](01_poisson_simulation.md) we saw that the counts of replicates of the same condition spread more than the Poisson predicts. How, then, is that spread written as a number, and how is the sequencing depth, which differs between samples, corrected? In this note I compute in turn the dispersion α, which sets how much the counts spread, the negative binomial distribution, and the size factor, which corrects for depth by multiplying the expected count.

> Study text 2.3–2.4, chapter 3 · Code: [02_negative_binomial.ipynb](../../02_negative_binomial.ipynb)

All code was run from top to bottom in a single session on R 4.5.2 and DESeq2 1.50.2, after first running `suppressPackageStartupMessages(library(DESeq2)); options(digits = 6, width = 100)`. So later blocks reuse objects created in earlier ones, and the output keeps DESeq2 messages such as "converting counts to integer mode" as they were.

## 1. What does the dispersion α measure?

Let me start with gene A, the example gene from chapter 16 of the study text. The counts of the three Ctrl replicates are 100, 130 and 90. Here the count $K_{ij}$ is the number of reads (or fragments) assigned to gene $i$ in sample $j$.

The mean of the three values is 106.7, so under the Poisson model of note 01 the variance should also be near 106.7. But the sample variance comes out at 433. With only three values this alone settles nothing, but spread beyond the Poisson prediction like this is common between biological replicates. In note 01 we called it overdispersion.

DESeq2 writes this excess spread as one number per gene.

$$\mathrm{Var}(K_{ij}) = \mu_{ij} + \alpha_i\,\mu_{ij}^2$$

$\mu_{ij}$ is the expected count, the mean count the model expects for gene $i$ in sample $j$. The first term, $\mu_{ij}$, is the share the Poisson already predicts: sampling noise from drawing reads at random. The second term, $\alpha_i\mu_{ij}^2$, is added on top; it comes from the expression level itself differing a little between replicates. The coefficient $\alpha_i$ that sets the size of this second term is the dispersion of gene $i$. If $\alpha=0$, the second term disappears and the model is the same as the Poisson.

The final dispersion of gene A comes out in [note 03](03_dispersion_estimation.md) as $\alpha=0.053147$, obtained with the prior that the study text sets for teaching. Plugging this value into the Ctrl mean gives a variance of $106.667 + 0.053147\times106.667^2 \approx 711$.

```r
ctrl <- c(100, 130, 90)                  # the three Ctrl replicates of gene A
m <- mean(ctrl); a <- 0.053147           # a: gene A's final dispersion from note 03
cat("mean:", m, "  sample variance of the three:", var(ctrl), "\n")
cat("if Poisson      var:", m, "  SD:", sqrt(m), "\n")
cat("NB (alpha = a)  var:", m + a*m^2, "  SD:", sqrt(m + a*m^2), "\n")
```

```
mean: 106.667   sample variance of the three: 433.333 
if Poisson      var: 106.667   SD: 10.328 
NB (alpha = a)  var: 711.361   SD: 26.6714 
```

The NB (negative binomial, covered in detail in section 2) standard deviation of 26.7 is about 2.6 times wider than the Poisson 10.3. One thing here is easy to confuse: α (0.053) is neither the variance (711) nor the standard deviation (26.7). At first I too thought of the dispersion as another name for the variance. The model variance of 711 also differs from the sample variance of the three Ctrl values, 433. This is because α is not set from the three Ctrl values alone; it is set by looking at the three Starvation values as well and going through the correction and the prior of note 03.

### The same α with a different mean?

I recomputed the table on p.6 of the study text. Rows 1–4 are the study text's table, rows 5–6 are the legend of the figure in section 2, and row 7 is the value from Problem 3.

```r
tab <- data.frame(mu = c(100,100,100,1000, 30,30, 50), alpha = c(0,0.01,0.10,0.10, 0.05,0.2, 0.2))
tab$var <- tab$mu + tab$alpha*tab$mu^2; tab$CV2 <- 1/tab$mu + tab$alpha
tab$sd_over_mu <- sqrt(tab$var)/tab$mu; print(tab)
```

```
    mu alpha    var       CV2 sd_over_mu
1  100  0.00    100 0.0100000   0.100000
2  100  0.01    200 0.0200000   0.141421
3  100  0.10   1100 0.1100000   0.331662
4 1000  0.10 101000 0.1010000   0.317805
5   30  0.05     75 0.0833333   0.288675
6   30  0.20    210 0.2333333   0.483046
7   50  0.20    550 0.2200000   0.469042
```

Compare rows 3 and 4. Both have α = 0.1, but the variances are 1,100 and 101,000, about 92 times apart. The size of the variance is set mainly by μ.

So the spread relative to the mean is looked at separately. The standard deviation divided by the mean is called the coefficient of variation (CV), and its square is this.

$$CV^2 = \frac{\mathrm{Var}(K)}{\mu^2} = \frac{1}{\mu} + \alpha$$

$1/\mu$ is the Poisson share and $\alpha$ the share of the extra spread. Plugging in row 3 gives $1/100 + 0.1 = 0.11$, and row 4 gives $1/1000 + 0.1 = 0.101$. Even though the absolute variances are 92 times apart, $CV^2$ is almost the same. As μ grows, $1/\mu$ approaches 0, so α becomes the floor below which the relative variation cannot go.

### Are low-expression genes a problem because their variance is large?

This time I fixed α at 0.1 and changed only μ. `poisson_share` is the fraction of $CV^2$ taken up by the Poisson share $1/\mu$.

```r
lo <- data.frame(mu = c(5, 10, 100, 1000, 10000), alpha = 0.1)
lo$var <- lo$mu + lo$alpha*lo$mu^2; lo$CV2 <- 1/lo$mu + lo$alpha
lo$poisson_share <- (1/lo$mu)/lo$CV2; print(lo)
```

```
     mu alpha       var    CV2 poisson_share
1     5   0.1 7.500e+00 0.3000   0.666666667
2    10   0.1 2.000e+01 0.2000   0.500000000
3   100   0.1 1.100e+03 0.1100   0.090909091
4  1000   0.1 1.010e+05 0.1010   0.009900990
5 10000   0.1 1.001e+07 0.1001   0.000999001
```

At μ = 5 the absolute variance is only 7.5. Instead, $CV^2$ is large at 0.30, and 2/3 of it is the Poisson share. Conversely, once μ exceeds 1,000, $CV^2$ sits almost on α.

So the statement "low-expression genes have a large absolute variance" is wrong. Low-expression genes have few reads, so their relative variation is large, and the estimates of α and of the effect are correspondingly unstable. The ways of handling this instability are the shrinkage of note 03 (pulling estimates that carry little information toward the overall tendency) and the independent filtering of [note 06](06_multiple_testing.md) (a procedure that removes genes with a low mean count from the multiple-testing adjustment).

## 2. What does the negative binomial look like?

Consider a gene with a mean of 30. Plotting the probability of each count value for α of 0, 0.05 and 0.2 gives the figure below, a redrawing of figure 1 on p.7 of the study text.

![Probability distributions with mean 30 and alpha 0, 0.05, 0.2](../../figures/02_nb_pmf_mean30.png)

All three curves have a mean of 30, but as α grows the peak gets lower and both tails, especially the tail toward large counts, get longer.

<details>
<summary>Code for the figure</summary>

```r
library(ggplot2)
k <- 0:100
lv <- c("Poisson (alpha = 0), variance 30", "NB, alpha = 0.05, variance 75", "NB, alpha = 0.2, variance 210")
d <- rbind(data.frame(k, p = dpois(k, 30),                    model = lv[1]),
           data.frame(k, p = dnbinom(k, size = 1/0.05, mu = 30), model = lv[2]),
           data.frame(k, p = dnbinom(k, size = 1/0.2,  mu = 30), model = lv[3]))
d$model <- factor(d$model, levels = lv)
p <- ggplot(d, aes(k, p, colour = model)) +
  geom_vline(xintercept = 30, colour = "#c3c2b7", linetype = "dashed", linewidth = 0.4) +
  annotate("text", x = 31, y = 0.0005, label = "mean = 30", hjust = 0, vjust = 0, size = 3.4, colour = "#52514e") +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("#2a78d6", "#eb6834", "#1baf7a"), name = NULL) +
  labs(x = "Observed count", y = "Probability",
       title = "Same mean (30), different dispersion alpha",
       subtitle = "Probability of each count; variance = mu + alpha * mu^2") +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
        axis.text = element_text(colour = "#52514e"), plot.subtitle = element_text(colour = "#52514e"),
        legend.position = "inside", legend.position.inside = c(0.98, 0.95), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = NA))
ggsave("../figures/02_nb_pmf_mean30.png", p, width = 7, height = 4.5, dpi = 150, bg = "white")
cat("saved:", file.exists("../figures/02_nb_pmf_mean30.png"), "\n")
```

```
saved: TRUE 
```

</details>

These curves are the negative binomial (NB), in a word a count distribution that allows overdispersion. DESeq2 takes each count to follow $K_{ij}\sim NB(\mu_{ij},\alpha_i)$, and the probability that a count is exactly $k$ is written like this.

$$P(K=k)=\frac{\Gamma(k+r)}{\Gamma(r)\,k!}\Big(\frac{r}{r+\mu}\Big)^{r}\Big(\frac{\mu}{r+\mu}\Big)^{k},\qquad r=\frac{1}{\alpha}$$

$k$ is a count value, one of 0, 1, 2, …; $\mu$ is the mean, and $r$ is the reciprocal of the dispersion, $1/\alpha$. $\Gamma$ is the gamma function; since $\Gamma(n)=(n-1)!$ for an integer $n$, think of it as the factorial extended to real numbers. The mean of this distribution is $\mu$, and the variance is $\mu+\mu^2/r=\mu+\alpha\mu^2$.

For the orange curve in the figure, $\mu=30$, $\alpha=0.05$, $r=20$, and the variance is $30+30^2/20=75$. As $\alpha\to0$, $r\to\infty$ and it goes back to the Poisson.

R's `dnbinom()` takes the same distribution as `size` and `mu`. Be sure to remember that `size` here is $r=1/\alpha$. Passing `size = α` gives a completely different distribution.

```r
mu <- 30; alpha <- 0.05; r <- 1/alpha
nb_manual <- function(k, mu, r) exp(lgamma(k+r) - lgamma(r) - lgamma(k+1) + r*log(r/(r+mu)) + k*log(mu/(r+mu)))
cat("max |formula - dnbinom|, k = 0..500:", max(abs(nb_manual(0:500, mu, r) - dnbinom(0:500, size = r, mu = mu))), "\n")
kk <- 0:2000; p <- dnbinom(kk, size = r, mu = mu)
cat("sum p =", sum(p), " E[K] =", sum(kk*p), " Var[K] =", sum(kk^2*p) - sum(kk*p)^2, " mu+alpha*mu^2 =", mu + alpha*mu^2, "\n")
cat("alpha -> 0: dnbinom(30, size=1e8, mu=30) =", dnbinom(30, size = 1e8, mu = 30), " dpois(30, 30) =", dpois(30, 30), "\n")
```

```
max |formula - dnbinom|, k = 0..500: 2.08167e-15 
sum p = 1  E[K] = 30  Var[K] = 75  mu+alpha*mu^2 = 75 
alpha -> 0: dnbinom(30, size=1e8, mu=30) = 0.0726345  dpois(30, 30) = 0.0726345 
```

The difference between the formula above and `dnbinom(size = 1/α)` is only rounding error on the order of $10^{-15}$. The variance computed directly from the probabilities also comes out exactly 75, and setting α close to 0 gives the same values as the Poisson.

The differences seen in the figure can be written as numbers too: the 2.5 % and 97.5 % quantiles, and the two tail probabilities.

```r
fig1 <- t(sapply(c(0, 0.05, 0.2), function(a) {
  sz <- if (a == 0) Inf else 1/a                      # with size=Inf, dnbinom gives the Poisson
  c(alpha = a, var = 30 + a*30^2, q2.5 = qnbinom(0.025, size = sz, mu = 30), q97.5 = qnbinom(0.975, size = sz, mu = 30),
    "P(K<=20)" = pnbinom(20, size = sz, mu = 30), "P(K>=45)" = 1 - pnbinom(44, size = sz, mu = 30)) }))
print(round(fig1, 4))
cat("size=Inf check: dnbinom(30, size=Inf, mu=30) =", dnbinom(30, size = Inf, mu = 30), " dpois(30,30) =", dpois(30, 30), "\n")
```

```
     alpha var q2.5 q97.5 P(K<=20) P(K>=45)
[1,]  0.00  30   20    41   0.0353   0.0063
[2,]  0.05  75   15    49   0.1298   0.0579
[3,]  0.20 210    8    64   0.2811   0.1524
size=Inf check: dnbinom(30, size=Inf, mu=30) = 0.0726345  dpois(30,30) = 0.0726345 
```

With the same mean of 30, the central 95 % range widens $[20,41]\to[15,49]\to[8,64]$, and the probability of 20 or less grows from 3.5 % to 28 % (because the distribution is discrete, the probability of falling in this range is not exactly 0.95). This puts into numbers what the study text says with figure 1: "even at the same mean, a larger dispersion widens the range of possible counts". The last line shows that treating α = 0 as `size = Inf` also gives the same as the Poisson.

But where does the NB come from? In note 01 we let the Poisson rate itself wobble as a gamma distribution from sample to sample, and saw that the count then has mean μ and variance $\mu+\alpha\mu^2$. Integrating that mixture gives exactly the NB formula above. The derivation is written up in "Deriving the NB formula from the Poisson–gamma mixture" under "Going deeper". But the mixture is only a device for obtaining the NB; DESeq2 does not estimate and store a hidden rate per sample (implementation check (3) under "Going deeper" in [note 01](01_poisson_simulation.md)).

## 3. Why does one gene have only one α?

Let me look at gene A's six values again. Ctrl is 100, 130, 90 (mean 106.7) and Starvation is 200, 250, 180 (mean 210), so the Starvation mean is about twice the Ctrl mean.

The expected count $\mu_{ij}$ can differ between samples. In this example all size factors are 1, so μ is 106.7 for the Ctrl samples and 210 for the Starvation samples. The dispersion $\alpha_A$, on the other hand, is a single value for gene A. The difference in means between the two groups is explained by the mean coefficient β (section 5), and α takes charge only of the spread left after that, that is, how much each sample scatters around its own mean. Section 2.4 of the study text explains this with three groups (Ctrl, Starvation, Glucose), but the point does not depend on the number of groups. Glucose here is the short name for the condition where starved cells are given glucose again (Starvation+Glucose).

What happens, then, if the difference in means is left out of the model? I estimated gene A's α with two designs and compared them. A design is a formula stating what the means of the samples depend on. `~ condition` means "give Ctrl and Starvation separate means", and `~ 1` means "give all six samples one mean".

The estimate is the maximum likelihood estimate (MLE). The likelihood, with the observations held fixed, says how plausible each candidate α makes those observations. The MLE is the α that makes that likelihood largest. The code below does the same thing turned around: it minimizes the log likelihood with a minus sign in front (`nll`). The detailed calculation is in note 03.

```r
K_A <- c(100, 130, 90, 200, 250, 180)            # gene A: 3 Ctrl, 3 Starvation
grp <- factor(rep(c("Ctrl", "Starvation"), each = 3))
mu_cond <- ave(K_A, grp)                          # ~condition: each sample's mean = its group mean
mu_one  <- rep(mean(K_A), 6)                      # ~1: the six samples share one mean
nll <- function(logA, mu) -sum(dnbinom(K_A, size = exp(-logA), mu = mu, log = TRUE))
mle <- function(mu) exp(optimize(nll, log(c(1e-8, 10)), mu = mu, tol = 1e-12)$minimum)
cat("group means:", unique(mu_cond), "  overall mean:", mu_one[1], "\n")
cat("alpha MLE   ~condition:", mle(mu_cond), "   ~1:", mle(mu_one), "\n")
```

```
group means: 106.667 210   overall mean: 158.333 
alpha MLE   ~condition: 0.014786    ~1: 0.12509 
```

The 0.014786 of `~ condition` is the same as the value that appears in note 03 as the "plain NB MLE". But with `~ 1`, where the means are tied into one, it becomes 0.125, about 8.5 times larger. The twofold difference between the groups has leaked into α.

The same thing happens in a simulation with 1,000 genes. I picked the 332 genes whose true log2 fold change (LFC, the log2 of the ratio of two condition means) has an absolute value above 2 and compared the medians of the dispersions estimated separately per gene (gene-wise estimates). It was 0.19 with `~ condition` and 1.32 with `~ 1`, while the median of the true values is 0.37. `~ 1` inflated far above the truth. That the 0.19 of `~ condition` is below the truth is the underestimation that comes from estimating each gene separately from 6 samples, which is covered in note 03. The code for this simulation is under "Going deeper".

In one line: **β explains the structure of the mean, and α the overdispersion left around that structure.** "The difference in means is large" and "the variation within a condition is large" are different pieces of information.

## 4. What if the sequencing depth differs between samples?

Suppose sample 2 was sequenced twice as deeply as sample 1. Then even at the same expression level, the counts of sample 2 come out roughly twice as large. This difference must not be read as a condition effect, so the depth difference has to be corrected first. Before that, let me go over what goes into DESeq2.

### Why must the input be raw counts?

What goes into DESeq2 is a matrix of non-negative integer counts. These are the original counts, not divided or transformed in any way, so they are called raw counts. TPM, FPKM, log values, VST values (values transformed so that the variance depends less on the mean, [note 07](07_lfc_shrinkage_and_qc.md)) and values already divided by library size do not go in. The NB itself is a distribution on the integers, and DESeq2 re-estimates the depth correction from the counts and multiplies it into the expected counts inside the model. Entering already normalized values would apply the correction twice or break the integer assumption. So DESeq2 blocks some of these at the input stage.

```r
cd <- data.frame(row.names = c("a","b"), cond = factor(c("x","y")))
tpm_like <- matrix(c(10.5, 20.25, 3.7, 8.1), 2, 2, dimnames = list(c("g1","g2"), c("a","b")))
msg <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))
print(msg(DESeqDataSetFromMatrix(tpm_like, cd, ~cond)))
print(msg(DESeqDataSetFromMatrix(matrix(c(1,-1,2,3), 2, 2, dimnames = dimnames(tpm_like)), cd, ~cond)))
x <- DESeqDataSetFromMatrix(matrix(c(1,2,3,4), 2, 2, dimnames = dimnames(tpm_like)), cd, ~cond)
cat("double but whole numbers:", class(x), "/ storage mode:", storage.mode(counts(x)), "\n")
```

```
[1] "some values in assay are not integers"
[1] "some values in assay are negative"
converting counts to integer mode
double but whole numbers: DESeqDataSet / storage mode: integer 
```

A real-valued matrix is rejected with a "not integers" error, and negative values with a "negative" error. Values like 1, 2, 3, 4 stored as double are accepted and converted to integer as long as they are whole numbers.

But a rounded TPM passes this check. That does not make it a correct input. TPM is already divided by length and depth, so size factor estimation and the variance structure of the NB lose their meaning. The estimated counts of Salmon or kallisto can be imported through `tximport` → `DESeqDataSetFromTximport()`, or `tximeta` → `DESeqDataSet()`. For details, see the tximport part under "Going deeper".

These are also worth checking once before an analysis.

- whether the counts are in reads or in fragments
- whether gene IDs are unique
- whether the columns of the count matrix and the rows of the metadata are in the same order
- whether condition was not stored as a number

### Size factor: a multiplier on the expected count

Even for a gene at the same expression level, we expect about twice the count in a sample read twice as deeply. DESeq2 writes this as a multiplication.

$$\mu_{ij} = s_j\,q_{ij}$$

$s_j$ is the size factor of sample $j$, the multiplier that corrects for the difference in sequencing depth. $q_{ij}$ is the expression level with the depth removed, which the study text calls the "normalized expected expression".

For example, even for a gene with $q=100$, the expected count in a sample with $s=2$ is 200. Conversely, the count divided by the size factor, $K_{ij}/s_j$, is called the normalized count.

### Computing the size factors: median-of-ratios

Using the ratio of total read counts as the size factor would be simple, but the total read count is pulled strongly by a few highly expressed genes. So DESeq2 computes a ratio per gene and then takes their median, hence the name median-of-ratios. The calculation has three steps.

1. For each gene, compute the geometric mean $g_i$ of its counts over all samples. The geometric mean is the $n$-th root of the product of $n$ values. This becomes the "virtual average sample", the reference.
2. Divide each count by the reference to get the ratio $K_{ij}/g_i$.
3. For each sample, take the median of these ratios; that value is the size factor $\hat s_j$.

$$g_i=\Big(\prod_{j=1}^{n}K_{ij}\Big)^{1/n},\qquad \hat s_j=\operatorname{median}_i\frac{K_{ij}}{g_i}$$

$n$ is the number of samples, and the median is taken along the gene direction ($i$). I computed it directly with the table on p.8 of the study text. To avoid confusion with gene A of the shared example, the study text's genes A, B and C are called g1, g2 and g3 here.

```r
K <- matrix(c(100,200, 50,100, 200,400), nrow = 3, byrow = TRUE,
            dimnames = list(c("g1","g2","g3"), c("sample1","sample2")))
g <- exp(rowMeans(log(K)))              # per-gene geometric mean = reference
ratio <- K / g                          # each count as a multiple of the reference
print(cbind(K, g = g, ratio))
s_hat <- apply(ratio, 2, median)        # median per sample
cat("size factor by hand:", s_hat, "  ratio s2/s1 =", s_hat[2]/s_hat[1], "\n")
cat("estimateSizeFactorsForMatrix:", estimateSizeFactorsForMatrix(K), "\n")
cat("normalized count K/s:\n"); print(t(t(K) / s_hat))
```

```
   sample1 sample2        g  sample1 sample2
g1     100     200 141.4214 0.707107 1.41421
g2      50     100  70.7107 0.707107 1.41421
g3     200     400 282.8427 0.707107 1.41421
size factor by hand: 0.707107 1.41421   ratio s2/s1 = 2 
estimateSizeFactorsForMatrix: 0.707107 1.41421 
normalized count K/s:
    sample1  sample2
g1 141.4214 141.4214
g2  70.7107  70.7107
g3 282.8427 282.8427
```

The geometric mean of g1 is $\sqrt{100\times200}=141.42$, and its ratios are $100/141.42=0.707$ and $200/141.42=1.414$. All three genes have the same ratios, so the medians are also 0.707 and 1.414. The hand calculation and `estimateSizeFactorsForMatrix()` give the same values, and the normalized counts become equal within each row. What matters here is not the absolute value 0.707 but the ratio of 2 between the two samples.

In this small example, though, the ratio of total read counts (350 : 700) is also 2, so the difference between the two does not show. So I added a gene, g4, that is 20 times higher in sample 2 only.

```r
K4 <- rbind(K, g4 = c(1000, 20000))     # add one gene that is 20x higher in sample2 only
sf4 <- estimateSizeFactorsForMatrix(K4); tc4 <- colSums(K4)/mean(colSums(K4))
cat("median-of-ratios:", round(sf4, 4), " -> ratio", round(sf4[2]/sf4[1], 4), "\n")
cat("by total read count (colSums/mean):", round(tc4, 4), " -> ratio", round(tc4[2]/tc4[1], 4), "\n")
cat("ratio matrix:\n"); print(round(K4/exp(rowMeans(log(K4))), 4))
```

```
median-of-ratios: 0.7071 1.4142  -> ratio 2 
by total read count (colSums/mean): 0.1224 1.8776  -> ratio 15.3333 
ratio matrix:
   sample1 sample2
g1  0.7071  1.4142
g2  0.7071  1.4142
g3  0.7071  1.4142
g4  0.2236  4.4721
```

Because of g4 alone, the ratio based on total read counts jumps to 15.3, but the median-of-ratios ratio stays at 2. However extreme g4's ratios are (0.22, 4.47), the median is set by g1–g3.

The principle is the same for data with thousands of genes. The figure below is a histogram of the ratios $K_{ij}/g_i$ per gene in one sample of simulated data.

![Per-gene ratios in one sample and their median](../../figures/02_median_of_ratios.png)

The ratios of the individual genes are widely scattered; even the central 80 % spans 0.91–3.73. The size factor is the single median of this distribution (1.95), the same as the value from `sizeFactors()`. Genes with even one 0 have a geometric mean of 0 and drop out of the calculation, so only 1,469 of the 2,000 genes were used.

<details>
<summary>Code for the figure</summary>

```r
library(ggplot2)
set.seed(42)
dsim <- makeExampleDESeqDataSet(n = 2000, m = 6, sizeFactors = c(0.5, 1, 1.5, 0.8, 1.2, 2))
cnt <- counts(dsim)
use <- rowSums(cnt == 0) == 0                          # genes with a 0 drop out of the default ratio method
g <- exp(rowMeans(log(cnt[use, ])))
rat <- cnt[use, "sample6"] / g
med <- median(rat)
cat("genes used:", sum(use), "/", nrow(cnt), "  median =", med,
    "  sizeFactors(sample6) =", sizeFactors(estimateSizeFactors(dsim))["sample6"], "\n")
cat("10% and 90% quantiles of the ratios:", quantile(rat, c(0.1, 0.9)), "\n")
p <- ggplot(data.frame(rat), aes(rat)) +
  geom_histogram(bins = 60, fill = "#2a78d6", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = med, colour = "#0b0b0b", linewidth = 0.8) +
  annotate("text", x = med / 1.06, y = 100, hjust = 1, size = 3.8, colour = "#0b0b0b",
           label = sprintf("size factor of sample 6 = median = %.2f", med)) +
  scale_x_continuous(trans = "log2", breaks = 2^(-3:3), labels = c("1/8", "1/4", "1/2", "1", "2", "4", "8")) +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Ratio to reference, K / g  (log2 axis)", y = "Number of genes",
       title = "Median-of-ratios for one sample",
       subtitle = sprintf("Simulated data, sample 6: one ratio per gene (%d genes without zeros)", sum(use))) +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
        axis.text = element_text(colour = "#52514e"), plot.subtitle = element_text(colour = "#52514e"))
ggsave("../figures/02_median_of_ratios.png", p, width = 7, height = 4.5, dpi = 150, bg = "white")
cat("saved:", file.exists("../figures/02_median_of_ratios.png"), "\n")
```

```
genes used: 1469 / 2000   median = 1.94983   sizeFactors(sample6) = 1.94983 
10% and 90% quantiles of the ratios: 0.906005 3.72959 
saved: TRUE 
```

</details>

### When the median moves

Median-of-ratios relies on the assumption that "most genes do not change between samples". So p.8 of the study text warns that this assumption breaks down when most genes change strongly in one direction. Trying it in a simulation, that is exactly what happened. When I scaled the counts of donors p1, p2 and p3 by 1×, 2× and 4× in 80 % of the genes, log2 1.70 of the log2 2 difference between p1 and p3 was absorbed into the size factors. The remaining 20 % of genes, which had not changed, then looked after normalization as if they dropped about 2.7-fold from p1 to p3. Even with the same manipulation on only 20 % of the genes, part of it (log2 0.40) was absorbed.

For such cases the study text says to use an independent reference, such as spike-ins or known control genes. In DESeq2, `estimateSizeFactors(dds, controlGenes = ...)` plays that role, taking the median only over the specified genes. Giving the unchanged genes of the example in section 6 as controls brings the size factors back near 1, and the effect estimates also recover toward their theoretical values. The code is under "Going deeper".

### What if a count is 0?

The default method (`type = "ratio"`) drops genes with even one 0 from the median calculation. So if every gene has at least one 0, no genes are left to use and it fails with an error. The alternative for this case is `type = "poscounts"`. Changing the input with a pseudocount (adding a small number such as 1 to every count) is not an alternative. As the study text's phrase "alternatives such as" suggests, there is also `type = "iterate"`, which I did not run in this note. How `poscounts` builds its reference is followed by hand under "Going deeper".

## 5. Where does the size factor enter the model?

For gene A of the shared example, all size factors were set to 1. If one Ctrl sample had a size factor of 2, its expected count should be 213.3 instead of 106.7, even at the same expression level. DESeq2 handles this multiplier by putting it inside the GLM.

A GLM (generalized linear model) is a model that writes the mean of the counts on the log scale as a sum of factors such as condition and batch. DESeq2's formula looks like this.

$$\log\mu_{ij}=\log s_{j}+x_j^{T}b_i$$

$x_j$ is row $j$ of the design matrix $X$. The design matrix is a table that records in numbers which condition, batch or pair each sample belongs to. $b_i$ is the coefficient vector of gene $i$, in natural-log units, and $x_j^Tb_i$ is the sum of the products of sample $j$'s row and this coefficient vector. The coefficients converted to log2 units are written $\beta_i=b_i/\log 2$, and the same formula can then also be written $\mu_{ij}=s_j\,2^{x_j^T\beta_i}$.

The remaining $\log s_j$ is the offset. An offset is a term that is added to the formula but not estimated; in other words, a term whose coefficient is fixed at 1.

Let me plug in gene A. The design is `~ condition`, with Ctrl as the reference.

```r
X <- cbind(Intercept = 1, Starvation = rep(0:1, each = 3))   # design matrix: ~ condition
b <- c(log(320/3), log(210 / (320/3)))                    # natural-log coefficients: log(Ctrl mean), log(Starvation mean / Ctrl mean)
s_A <- rep(1, 6)                                          # size factors of the gene A example
mu_A <- s_A * exp(X %*% b)
print(cbind(X, mu = drop(mu_A)))
cat("log2 coefficients beta = b / log(2):", b / log(2), "\n")
cat("expected count of sample 1 if only its size factor were 2:", 2 * exp(sum(X[1, ] * b)), "\n")
```

```
     Intercept Starvation      mu
[1,]         1          0 106.667
[2,]         1          0 106.667
[3,]         1          0 106.667
[4,]         1          1 210.000
[5,]         1          1 210.000
[6,]         1          1 210.000
log2 coefficients beta = b / log(2): 6.73697 0.97728
expected count of sample 1 if only its size factor were 2: 213.333
```

The Ctrl rows have $x=(1,0)$, so $\mu=\exp(b_0)=106.667$; the Starvation rows have $x=(1,1)$, so $\mu=\exp(b_0+b_1)=210$. The Starvation coefficient converted to log2, $\beta_1=0.977280$, is exactly the value used as gene A's LFC in the shared example. If the size factor becomes 2, the offset $\log 2$ is added and the expected count becomes 213.3.

What "a coefficient fixed at 1" means becomes clear by changing the size factors. I ran `DESeq()` on simulated data, then, with the dispersions left as they were, doubled all the size factors and recomputed the coefficients. The coefficients in the output are in log2 units.

```r
set.seed(7)
d1 <- makeExampleDESeqDataSet(n = 300, m = 6, betaSD = 1, sizeFactors = c(0.5, 1, 1.5, 0.8, 1.2, 2))
d1 <- DESeq(d1, quiet = TRUE)
d2 <- d1; sizeFactors(d2) <- 2 * sizeFactors(d1)   # dispersions(d2) stay those of d1
d2 <- nbinomWaldTest(d2, quiet = TRUE)
print(head(cbind(coef(d1), coef(d2)), 3))
cat("range of intercept differences:", range(coef(d2)[,1] - coef(d1)[,1], na.rm = TRUE), "\n")
cat("max |difference| of condition coefficient:", max(abs(coef(d2)[,2] - coef(d1)[,2]), na.rm = TRUE), "\n")
mu1 <- assays(d1)[["mu"]]; mu2 <- assays(d2)[["mu"]]
cat("max relative difference of mu matrix:", max(abs(mu2/mu1 - 1), na.rm = TRUE), "\n")
```

```
      Intercept condition_B_vs_A Intercept condition_B_vs_A
gene1   8.68204         0.587260   7.68204         0.587260
gene2  -1.22403         2.942984  -2.22402         2.942979
gene3   3.11048        -0.510782   2.11048        -0.510783
range of intercept differences: -1.00004 -0.999976 
max |difference| of condition coefficient: 3.67796e-05 
max relative difference of mu matrix: 2.4514e-05 
```

The intercept drops by 1 (difference −1.00004 to −0.999976) and the condition coefficient stays within $4\times10^{-5}$. In $\log_2\mu=\log_2 s+\beta_0+\beta_1x$, when $s$ doubles, the only way to keep the same μ is to lower $\beta_0$ by 1. The remaining differences on the order of $10^{-5}$ are numerical error of the iterative calculation. Factors such as the convergence tolerance and the starting values are mixed in, and I did not split out the share of each. If it were a freely estimated coefficient, there would be no reason for it to move by exactly 1 like this.

The size factors are computed from the data too, but in the GLM step they are treated as already fixed constants. For comparison, putting $\log_2 s_j$ into the design as a column makes its coefficient get estimated separately for each gene, scattering from −2.9 to 4.0.

One thing to be careful about: "changing the multiplier leaves the condition coefficient unchanged" holds when the dispersions are fixed. Rerunning `DESeq()` from scratch, the starting values of the dispersion estimation respond to the scale of the size factors, and the results change a little. The two experiments are collected in "Changing the scale of the size factors and rerunning DESeq()" under "Going deeper".

### Gene length and normalization factors

DESeq2 compares the same gene across samples, so if the gene length is fixed, it cancels out in the ratio. But if the isoform composition changes so that the effective length differs between samples, it does not cancel. In that case a genes×samples matrix $s_{ij}$ is used instead of one $s_j$ per sample; DESeq2 calls this a normalization factor.

When importing with `tximport` and `countsFromAbundance = "no"` (the default), the transcript lengths are stored as `avgTxLength`, and `estimateSizeFactors()` builds normalization factors from them. A matrix you made yourself can be supplied with `estimateSizeFactors(normMatrix = ...)`. When normalization factors exist, the size factors are not used.

## 6. What do normalized counts not correct for?

`counts(dds, normalized = TRUE)` is the count divided by the size factor, $K_{ij}/s_j$ ($K_{ij}/s_{ij}$ if there are normalization factors). It is a value with only the depth difference corrected.

If, in the three-group example, Ctrl, Starvation and Glucose were obtained from the same donors, the design becomes `~ pair + condition`, where a pair is the set of samples from the same donor. Does the pair difference then also disappear from the normalized counts when pair is in the design? At first glance it seems it would.

I checked by simulation. Three donors (p1, p2, p3) each gave one sample in condition A and one in B, and a donor effect (p1 ×1, p2 ×2, p3 ×4) was multiplied into genes 1–200 only, out of 1,000 genes. The true condition effect is 0 for every gene.

```r
set.seed(3)
dp <- makeExampleDESeqDataSet(n = 1000, m = 6, betaSD = 0)   # true condition effect = 0 (all genes)
dp$pair <- factor(rep(c("p1","p2","p3"), each = 2)); dp$condition <- factor(rep(c("A","B"), 3))   # one A and one B in each pair
cnt <- counts(dp); eff <- c(p1 = 1, p2 = 2, p3 = 4)[as.character(dp$pair)]
cnt[1:200, ] <- round(t(t(cnt[1:200, ]) * eff)); storage.mode(cnt) <- "integer"; counts(dp) <- cnt
design(dp) <- ~ pair + condition; dp <- DESeq(dp, quiet = TRUE)
cat("sizeFactors:", round(sizeFactors(dp), 3), "\n")
nc <- counts(dp, normalized = TRUE)
pm <- sapply(split(seq_len(6), dp$pair), function(j) rowMeans(nc[, j, drop = FALSE]))
cat("per-pair geometric mean of normalized counts, gene 1-200 (pair effect 1, 2, 4x):\n"); print(round(exp(colMeans(log(pm[1:200, ] + 0.5))), 2))
cat("gene 201-1000 (no pair effect):\n"); print(round(exp(colMeans(log(pm[201:1000, ] + 0.5))), 2))
cat("resultsNames:", resultsNames(dp), "\n")
cat("median GLM pair coefficient (log2), gene 1-200:", round(apply(coef(dp)[1:200, 2:3], 2, median, na.rm = TRUE), 3), " (theory 1, 2)\n")
cat("gene 201-1000:", round(apply(coef(dp)[201:1000, 2:3], 2, median, na.rm = TRUE), 3), " (theory 0, 0)\n")
lb <- limma::removeBatchEffect(log2(nc + 1), batch = dp$pair, design = model.matrix(~ condition, colData(dp)))
cat("pair means of gene 1-200 after limma::removeBatchEffect:\n"); print(round(2^sapply(split(seq_len(6), dp$pair), function(j) mean(lb[1:200, j])), 2))
```

```
sizeFactors: 0.877 0.91 1.028 1.051 1.158 1.196 
per-pair geometric mean of normalized counts, gene 1-200 (pair effect 1, 2, 4x):
   p1    p2    p3 
16.75 27.03 47.75 
gene 201-1000 (no pair effect):
   p1    p2    p3 
15.73 14.60 12.95 
resultsNames: Intercept pair_p2_vs_p1 pair_p3_vs_p1 condition_B_vs_A 
median GLM pair coefficient (log2), gene 1-200: 0.67 1.652  (theory 1, 2)
gene 201-1000: -0.149 -0.361  (theory 0, 0)
pair means of gene 1-200 after limma::removeBatchEffect:
  p1   p2   p3 
26.2 26.2 26.2 
```

The pair effect is still there in the normalized counts. The per-pair geometric mean of genes 1–200 climbs 17 → 27 → 48. What absorbs this effect are the GLM coefficients `pair_p2_vs_p1` and `pair_p3_vs_p1`. To remove the pair effect for a figure, a separate operation such as `limma::removeBatchEffect` is needed, after which all three pair means become 26.2. Because this is a balanced design with one A and one B in each pair, subtracting the estimated pair means makes them exactly equal by construction.

That the medians of the pair coefficients come out as 0.67 and 1.65 instead of the theoretical 1 and 2 is the "median moves" case from section 4. When 20 % of the genes moved in one direction, the size factors absorbed part of it (log2 0.40), and as a result the genes without a pair effect (201–1000) were biased in the opposite direction by −0.15 and −0.36.

In the end, the model fits the raw counts and multiplies the expected counts by the size factors; on the log scale that is adding an offset. It is not a process of making normalized counts and then feeding those values back in as if they were integer observations. Normalized counts are values for figures and summaries.

## Summary

- α is not the variance but the size of the extra spread beyond the Poisson. The variance is $\mu+\alpha\mu^2$, so plugging α = 0.053147 into gene A's Ctrl (mean 106.7) gives 711.
- $CV^2=1/\mu+\alpha$. Low-expression genes have a small absolute variance and a large relative variation. The real problem is that they carry little information, so the estimates are unstable.
- In R, compute it with `dnbinom(k, size = 1/α, mu = μ)`. At the same mean, a larger α widens the range of counts that can appear.
- There is one α per gene, while μ differs between samples. Leaving the difference in means out of the design inflates α (gene A: 0.0148 → 0.125).
- The size factor is computed by median-of-ratios, not as a ratio of total read counts, and enters $\log\mu_{ij}=\log s_j+x_j^Tb_i$ as an offset with its coefficient fixed at 1. When most genes move in one direction, this median moves too.
- Normalized counts are only $K_{ij}/s_j$, so pair and batch effects put in the design remain in them.

## Exercises

**Problem 1.** The histogram of all gene counts in one sample does not look like a normal distribution. Is this a direct reason why DESeq2 cannot be used? Explain the direction in which DESeq2 assumes a distribution.

<details>
<summary>Solution</summary>

No. The NB assumption is that $K_{ij}\sim NB(\mu_{ij},\alpha_i)$ "with gene $i$ fixed, along the sample direction $j$". A histogram drawn across the genes of one sample mixes values drawn one each from thousands of distributions with different $\mu_{ij}$ and $\alpha_i$, so there is no reason for it to follow any particular distribution. The answer in appendix B of the study text is the same.

In addition, the intercepts of `makeExampleDESeqDataSet` are drawn from a normal distribution (log2 scale) with `interceptMean = 4, interceptSD = 2` (`beta <- cbind(rnorm(n, interceptMean, interceptSD), …)`, `mu <- t(2^(x %*% t(beta)) * sizeFactors)`; the source block under "Going deeper"). So in that simulation the per-gene log2 means follow a normal distribution, but this is only the simulator's choice and has nothing to do with the NB assumption along the sample direction.

</details>

**Problem 3.** Compute the variance and CV² of the count when μ=50 and α=0.2.

<details>
<summary>Solution</summary>

```r
cat("Var =", 50 + 0.2*50^2, " CV^2 =", 1/50 + 0.2, " CV =", sqrt(1/50+0.2), "\n")
cat("if Poisson, Var =", 50, " CV^2 =", 1/50, "\n")
```

```
Var = 550  CV^2 = 0.22  CV = 0.469042 
if Poisson, Var = 50  CV^2 = 0.02 
```

$\mathrm{Var}=50+0.2\times2500=550$ and $CV^2=0.02+0.2=0.22$. The same as the study text's answer (550, 0.22). Of the variance of 550, 500 is the $\alpha\mu^2$ share, and of the $CV^2$ of 0.22, 0.2 is the α share. So α = 0.2 is neither the variance of 550 nor the standard deviation of 23.5.

</details>

**Problem 4.** Correct the error in the explanation that one more β coefficient for the size factor is added to `~pair+condition`.

<details>
<summary>Solution</summary>

The size factor enters $\log\mu_{ij}=\log s_j+x_j^Tb_i$ as an offset, with its coefficient fixed at 1. Only the pair and condition coefficients are estimated. This can be seen from two directions.

1. With the dispersions fixed, doubling the size factors moves the intercept by −1 within numerical tolerance, and the other coefficients stay within $4\times10^{-5}$ (section 5). This is a property that comes from the coefficient being fixed.
2. Putting the same $\log_2 s_j$ into the design as a column adds `logsf` to `resultsNames`, and its coefficient is estimated separately for each gene, from −2.9 to 4.0 ("Changing the scale of the size factors and rerunning DESeq()" under "Going deeper"). "Adding one more coefficient" looks exactly like this.

In the code too, `normalizationFactors` is multiplied outside `modelMatrix`, as in line 128 of `fitNbinomGLMs`, `mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))`.

</details>

## Going deeper

<details>
<summary>What I looked up in the DESeq2 source</summary>

I unpacked the function bodies with `deparse()`, numbered the lines, and printed only the lines I needed. The line numbers are for DESeq2 1.50.2.

First, `estimateSizeFactors`. The default `type` is `"ratio"` and the default `locfunc` is `median`, and the path splits four ways.

```r
m <- deparse(getMethod("estimateSizeFactors", "DESeqDataSet")@.Data)
cat(sprintf("[%d] %s", seq_along(m), m)[c(3:5, 12:13, 16, 22, 25, 27:31, 36:37, 40:41)], sep = "\n")
```

```
[3]     .local <- function (object, type = c("ratio", "poscounts", 
[4]         "iterate"), locfunc = stats::median, geoMeans, controlGenes, 
[5]         normMatrix, quiet = FALSE) 
[12]         if (type == "iterate") {
[13]             sizeFactors(object) <- estimateSizeFactorsIterate(object)
[16]             if (type == "poscounts") {
[22]                     exp(sum(log(x[x > 0]))/length(x))
[25]                 geoMeans <- apply(counts(object), 1, geoMeanNZ)
[27]             if ("avgTxLength" %in% assayNames(object)) {
[28]                 nm <- assays(object)[["avgTxLength"]]
[29]                 nm <- nm/exp(MatrixGenerics::rowMeans(log(nm)))
[30]                 normalizationFactors(object) <- estimateNormFactors(counts(object), 
[31]                   normMatrix = nm, locfunc = locfunc, geoMeans = geoMeans, 
[36]             else if (missing(normMatrix)) {
[37]                 sizeFactors(object) <- estimateSizeFactorsForMatrix(counts(object), 
[40]             else {
[41]                 normalizationFactors(object) <- estimateNormFactors(counts(object), 
```

`type = "iterate"` goes to a separate path, `estimateSizeFactorsIterate`, which repeats NB fits (lines 12–13). If there is an `avgTxLength` assay (tximport, `countsFromAbundance = "no"`), it fills `normalizationFactors` (lines 27–31), and supplying `normMatrix =` directly also fills `normalizationFactors` (lines 40–41). The remaining default path fills `sizeFactors` (lines 36–37). `poscounts` builds a geometric mean by dividing the sum of logs, with the 0s dropped, by the total number of samples `length(x)`, and passes it on as `geoMeans` (lines 22, 25).

Next is `estimateSizeFactorsForMatrix`, which actually computes the median-of-ratios.

```r
e <- deparse(estimateSizeFactorsForMatrix)
cat(sprintf("[%d] %s", seq_along(e), e)[c(5:6, 8, 19, 23, 25:26, 28:31, 38:41, 44:45)], sep = "\n")
```

```
[5]     if (missing(geoMeans)) {
[6]         incomingGeoMeans <- FALSE
[8]             loggeomeans <- MatrixGenerics::rowMeans(log(counts))
[19]         incomingGeoMeans <- TRUE
[23]         loggeomeans <- log(geoMeans)
[25]     if (all(is.infinite(loggeomeans))) {
[26]         stop("every gene contains at least one zero, cannot compute log geometric means")
[28]     sf <- if (missing(controlGenes)) {
[29]         apply(counts, 2, function(cnts) {
[30]             exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & 
[31]                 cnts > 0]))
[38]         loggeomeansSub <- loggeomeans[controlGenes]
[39]         apply(counts[controlGenes, , drop = FALSE], 2, function(cnts) {
[40]             exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & 
[41]                 cnts > 0]))
[44]     if (incomingGeoMeans) {
[45]         sf <- sf/exp(mean(log(sf)))
```

It takes the median of $\log K_{ij}-\log g_i$ on the log scale and exponentiates it (lines 29–31). `is.finite(loggeomeans) & cnts > 0` is the filter that drops genes containing a 0. If `controlGenes` is given, the median is taken only over those genes (lines 38–41); this is where the spike-ins / control genes mentioned in the study text come in. The filter that drops genes with a 0 is the same on this path too.

Lines 44–45 rescale the geometric mean of the size factors to 1 only when the `geoMeans` argument was supplied. The `poscounts` path of the dds is the typical case, and it also applies when the user passes `geoMeans =` directly. The default ratio path has no such rescaling. Indeed, in "Size factor details" below, the geometric mean on the default path was 1.0375, and passing the same reference as `geoMeans` made it 1.

Which do `counts(normalized = TRUE)` and `getSizeOrNormFactors` look at first, normalization factors or size factors?

```r
cm <- deparse(getMethod("counts", "DESeqDataSet")@.Data)
cat(sprintf("[%d] %s", seq_along(cm), cm)[23:31], sep = "\n")
```

```
[23]             if (!is.null(normalizationFactors(object))) {
[24]                 return(cnts/normalizationFactors(object))
[25]             }
[26]             else if (is.null(sizeFactors(object)) || any(is.na(sizeFactors(object)))) {
[27]                 stop("first calculate size factors, add normalizationFactors, or set normalized=FALSE")
[28]             }
[29]             else {
[30]                 return(t(t(cnts)/sizeFactors(object)))
[31]             }
```

```r
cat(deparse(DESeq2:::getSizeOrNormFactors), sep = "\n")
```

```
function (object) 
{
    if (!is.null(normalizationFactors(object))) {
        return(normalizationFactors(object))
    }
    else {
        return(matrix(rep(sizeFactors(object), each = nrow(object)), 
            ncol = ncol(object)))
    }
}
```

Both ignore the size factors when normalization factors exist. `DESeq()` checks in the same order.

```r
d <- deparse(DESeq); cat(sprintf("[%d] %s", 96:110, d[96:110]), sep = "\n")
```

```
[96]     if (!is.null(sizeFactors(object)) || !is.null(normalizationFactors(object))) {
[97]         if (!quiet) {
[98]             if (!is.null(normalizationFactors(object))) {
[99]                 message("using pre-existing normalization factors")
[100]             }
[101]             else {
[102]                 message("using pre-existing size factors")
[103]             }
[104]         }
[105]     }
[106]     else {
[107]         if (!quiet) 
[108]             message("estimating size factors")
[109]         object <- estimateSizeFactors(object, type = sfType, 
[110]             quiet = quiet)
```

If there are normalization factors it prints "using pre-existing normalization factors", if there are only size factors "using pre-existing size factors", and uses those values as they are. Only when neither exists does it call `estimateSizeFactors(object, type = sfType)`. So `DESeq()` does not overwrite size factors set in advance. Note that normalization factors must be positive finite values.

```r
nm <- deparse(getMethod("normalizationFactors<-", c("DESeqDataSet", "matrix"))@.Data)
cat(sprintf("[%d] %s", seq_along(nm), nm)[5:7], sep = "\n")
```

```
[5]         stopifnot(all(!is.na(value)))
[6]         stopifnot(all(is.finite(value)))
[7]         stopifnot(all(value > 0))
```

Now the place where the offset actually enters, `DESeq2:::fitNbinomGLMs`.

```r
f <- deparse(DESeq2:::fitNbinomGLMs)
cat(sprintf("[%d] %s", seq_along(f), f)[c(3, 31, 44, 47:48, 123:128, 143:144, 153:154)], sep = "\n")
```

```
[3]     useOptim = TRUE, useQR = TRUE, forceOptim = FALSE, warnNonposVar = TRUE, 
[31]     normalizationFactors <- getSizeOrNormFactors(object)
[44]     if (type == "glmGamPoi") {
[47]         gp_res <- glmGamPoi::glm_gp(counts(object), design = modelMatrix, 
[48]             size_factors = FALSE, offset = log(normalizationFactors), 
[123]     betaRes <- fitBetaWrapper(ySEXP = counts(object), xSEXP = modelMatrix, 
[124]         nfSEXP = normalizationFactors, alpha_hatSEXP = alpha_hat, 
[125]         beta_matSEXP = beta_mat, lambdaSEXP = lambdaNatLogScale, 
[126]         weightsSEXP = weights, useWeightsSEXP = useWeights, tolSEXP = betaTol, 
[127]         maxitSEXP = maxit, useQRSEXP = useQR, minmuSEXP = minmu)
[128]     mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))
[143]     rowsForOptim <- if (useOptim) {
[144]         which(!betaConv | !rowStable | !rowVarPositive)
[153]         resOptim <- fitNbinomGLMsOptim(object, modelMatrix, lambda, 
[154]             rowsForOptim, rowStable, normalizationFactors, alpha_hat, 
```

Line 128 is $\mu_{ij}=s_{ij}\exp(x_j^Tb_i)$ translated directly into code. `normalizationFactors` is just a constant that gets multiplied in, not a column of `modelMatrix`. This is how "coefficient fixed at 1" is implemented in the code.

`getSizeOrNormFactors` on line 31 builds a genes×samples matrix even when there are only size factors. This matrix goes into `nfSEXP` of `fitBetaWrapper`, the default C++ IRLS routine (lines 123–124). IRLS is the iteratively reweighted least squares calculation that finds the coefficients. The same matrix is passed unchanged to two other paths: with `type = "glmGamPoi"` it goes into `offset = log(normalizationFactors)` of `glmGamPoi::glm_gp` (lines 44, 47–48), and rows where the C++ IRLS did not converge or was unstable go into the arguments of `fitNbinomGLMsOptim` (lines 143–144, 153–154; `useOptim = TRUE` is the default on line 3). And `ySEXP = counts(object)` on line 123 is the raw counts. This line is the evidence for the study text's statement that normalized counts are not fed back in.

The input checks are done in `DESeqDataSet`.

```r
s <- deparse(DESeqDataSet); cat(sprintf("[%d] %s", 20:41, s[20:41]), sep = "\n")
```

```
[20]     if ("tximetaInfo" %in% names(metadata(se))) {
[21]         se <- processTximeta(se)
[22]     }
[23]     if (any(is.na(assay(se)))) 
[24]         stop("NA values are not allowed in the count matrix")
[25]     if (any(assay(se) < 0)) {
[26]         stop("some values in assay are negative")
[27]     }
[28]     if (!skipIntegerMode & !is.integer(assay(se))) {
[29]         if (!is.numeric(assay(se))) {
[30]             stop(paste("counts matrix should be numeric, currently it has mode:", 
[31]                 mode(assay(se))))
[32]         }
[33]         if (any(round(assay(se)) != assay(se))) {
[34]             stop("some values in assay are not integers")
[35]         }
[36]         message("converting counts to integer mode")
[37]         mode(assay(se)) <- "integer"
[38]     }
[39]     if (all(assay(se) == 0)) {
[40]         stop("all samples have 0 counts for all genes. check the counting script.")
[41]     }
```

If `metadata(se)` has `tximetaInfo`, it is handled by `processTximeta(se)` (lines 20–22), and then NA (lines 23–24) and negative values (lines 25–27) are rejected. Inside the block entered only when `skipIntegerMode = FALSE` (the default, line 28), non-integer values are rejected (lines 33–34), and doubles that are whole numbers are converted to integer with the message "converting counts to integer mode" (lines 36–37). If every count is 0, it is an error (lines 39–40).

tximport results come in through `DESeqDataSetFromTximport`.

```r
tx <- deparse(DESeqDataSetFromTximport); cat(sprintf("[%d] %s", 4:18, tx[4:18]), sep = "\n")
```

```
[4]     counts <- round(txi$counts)
[5]     mode(counts) <- "integer"
[6]     object <- DESeqDataSetFromMatrix(countData = counts, colData = colData, 
[7]         design = design, ...)
[8]     stopifnot(txi$countsFromAbundance %in% c("no", "scaledTPM", 
[9]         "lengthScaledTPM"))
[10]     if (txi$countsFromAbundance %in% c("scaledTPM", "lengthScaledTPM")) {
[11]         message("using just counts from tximport")
[12]     }
[13]     else {
[14]         message("using counts and average transcript lengths from tximport")
[15]         lengths <- txi$length
[16]         stopifnot(all(lengths > 0))
[17]         dimnames(lengths) <- dimnames(object)
[18]         assays(object)[["avgTxLength"]] <- lengths
```

Lines 4–5 apply `round()` and then convert to integer mode, and lines 10–18 decide whether to store `avgTxLength` depending on `countsFromAbundance`. The transcript lengths are kept as the `avgTxLength` assay only with `"no"`.

Finally, how `makeExampleDESeqDataSet`, used for the simulations in this note, generates counts.

```r
mk <- deparse(makeExampleDESeqDataSet); cat(sprintf("[%d] %s", seq_along(mk), mk)[c(1:3, 5:7, 16:18)], sep = "\n")
```

```
[1] function (n = 1000, m = 12, betaSD = 0, interceptMean = 4, interceptSD = 2, 
[2]     dispMeanRel = function(x) 4/x + 0.1, sizeFactors = rep(1, 
[3]         m)) 
[5]     beta <- cbind(rnorm(n, interceptMean, interceptSD), rnorm(n, 
[6]         0, betaSD))
[7]     dispersion <- dispMeanRel(2^(beta[, 1]))
[16]     mu <- t(2^(x %*% t(beta)) * sizeFactors)
[17]     countData <- matrix(rnbinom(m * n, mu = mu, size = 1/dispersion), 
[18]         ncol = m)
```

The defaults are `n = 1000, m = 12, betaSD = 0, interceptMean = 4, interceptSD = 2, dispMeanRel = function(x) 4/x + 0.1, sizeFactors = rep(1, m)`. The intercepts are drawn from a normal distribution on the log2 scale (lines 5–6), the expected counts are `2^(x %*% t(beta)) * sizeFactors` (line 16), and the counts are generated with `rnbinom(mu = mu, size = 1/dispersion)` (line 17).

</details>

<details>
<summary>Where the results differ from the study text</summary>

While running the study text's explanations as code, I collected only what came out differently from the study text or was not in it.

| Study text claim (p.) | Result | Evidence |
|---|---|---|
| Correspondence between $NB(\mu,\alpha)$ and R's `dnbinom(size=1/α, mu=μ)` (not stated in the study text) | Matches; a detail not in the study text | Maximum difference of $2\times10^{-15}$ in section 2 |
| "Rounding TPM to make it look like counts is not a correct reconstruction" (p.8) | Not confirmed (conceptual) | Rounded TPM passes the input check, so the code does not block it. Why it is wrong is a matter of model assumptions (explained in section 4), which code cannot decide |
| The warning that the invariant-reference assumption breaks down when most genes change strongly in one direction (p.8) | Matches, with an addition | When 80 % of the genes moved 4-fold, 1.70 of log2 2 was absorbed into the size factors ("When some or most genes move"). Even at the 20 % of section 6, weaker than the study text's condition, 0.40 was absorbed |
| (Not mentioned in the study text) values of `estimateSizeFactors(type="poscounts")` and `estimateSizeFactorsForMatrix(type="poscounts")` | Not a discrepancy; an addition | The ratios are the same, and only the dds path rescales the geometric mean to 1 ("Counts of 0: the exclusion rule and poscounts"). A fact not in the study text, recorded here only |
| The coefficient of $\log s_j$ is fixed at 1 and is not an estimated coefficient (p.9) | Matches (with the dispersions fixed) | With the dispersions fixed, the intercept moves by $-1.00004$ to $-0.999976$, that is, by $-1$ within numerical tolerance (section 5, line 128 of `fitNbinomGLMs`). Re-estimating the dispersions scatters it to $-1.012$ to $-0.971$ and changes the condition coefficients and padj a little, not because of the offset but because the starting value of the gene-wise dispersion (`roughDispEstimate`) depends on the scale ("Changing the scale of the size factors and rerunning DESeq()") |
| Notation: right after the formula $\log\mu_{ij}=\log s_j+x_j^Tb_i$ (natural-log coefficients $b$), the estimated coefficients are called "β coefficients" (p.9) | Notation mismatch (content unaffected) | Appendix C (p.39) distinguishes $b$ in natural-log units from $\beta$ in log2 units. This note keeps $b$ and $\beta$ apart (section 5) and renders the corresponding sentence of section 3.4 of the study text as "condition and pair coefficients" |
| Appendix B, Problem 1: the model along the sample direction | Matches (conceptual) | The conceptual answer is the same as the study text's. The remark about `interceptSD=2` of `makeExampleDESeqDataSet` is a supplementary explanation |

The rest of the study text's explanations matched the actual behaviour of DESeq2.

</details>

<details>
<summary>Deriving the NB formula from the Poisson–gamma mixture</summary>

As in section 2.2 of the study text, let $K\mid\Lambda\sim\mathrm{Poisson}(\Lambda)$, $E[\Lambda]=\mu$, $\mathrm{Var}(\Lambda)=\alpha\mu^2$. Taking $\Lambda$ as a gamma gives shape $r=1/\alpha$ and scale $\alpha\mu$, since we need shape·scale $=\mu$ and shape·scale² $=\alpha\mu^2$. This is the same $r$ as the gamma shape in note 01. The marginal pmf is then:

$$P(K=k)=\int_0^\infty\frac{e^{-\lambda}\lambda^{k}}{k!}\cdot\frac{\lambda^{r-1}e^{-\lambda/(\alpha\mu)}}{\Gamma(r)(\alpha\mu)^{r}}\,d\lambda=\frac{1}{k!\,\Gamma(r)(\alpha\mu)^{r}}\int_0^\infty\lambda^{k+r-1}e^{-\lambda\,(1+\alpha\mu)/(\alpha\mu)}\,d\lambda .$$

The integral at the end is the definition of the gamma function, $\int_0^\infty\lambda^{a-1}e^{-b\lambda}d\lambda=\Gamma(a)/b^{a}$, with $a=k+r$ and $b=(1+\alpha\mu)/(\alpha\mu)$, which gives $\Gamma(k+r)\,(\alpha\mu)^{k+r}/(1+\alpha\mu)^{k+r}$. Simplifying:

$$P(K=k)=\frac{\Gamma(k+r)}{\Gamma(r)\,k!}\Big(\frac{1}{1+\alpha\mu}\Big)^{r}\Big(\frac{\alpha\mu}{1+\alpha\mu}\Big)^{k}$$

Since $1/(1+\alpha\mu)=r/(r+\mu)$ and $\alpha\mu/(1+\alpha\mu)=\mu/(r+\mu)$, this is exactly the NB pmf of section 2.

I also checked it in numbers. Continuing with `nb_manual` from the block in section 2, I put the formula and `dnbinom` side by side at several $k$, and also compared with the `prob` parameterization and with random draws.

```r
mu <- 30; alpha <- 0.05; r <- 1/alpha
k <- c(0, 10, 20, 30, 40, 60)
cmp <- data.frame(k, manual = nb_manual(k, mu, r), dnbinom = dnbinom(k, size = r, mu = mu))
cmp$abs_diff <- abs(cmp$manual - cmp$dnbinom); print(cmp, digits = 10)
cat("same via prob parameterization: p = r/(r+mu) =", r/(r+mu), "; dnbinom(30,size=r,prob=r/(r+mu)) =", dnbinom(30, size = r, prob = r/(r+mu)), "\n")
set.seed(1); x <- rnbinom(1e6, size = r, mu = mu)
cat("rnbinom 1e6: mean =", mean(x), " var =", var(x), "\n")
```

```
   k          manual         dnbinom        abs_diff
1  0 1.099511628e-08 1.099511628e-08 0.000000000e+00
2 10 1.331660435e-03 1.331660435e-03 1.561251128e-17
3 20 2.770707453e-02 2.770707453e-02 2.671474153e-16
4 30 4.582342113e-02 4.582342113e-02 2.498001805e-16
5 40 2.053684624e-02 2.053684624e-02 2.532696275e-16
6 60 4.749447753e-04 4.749447753e-04 4.521123059e-17
same via prob parameterization: p = r/(r+mu) = 0.4 ; dnbinom(30,size=r,prob=r/(r+mu)) = 0.0458234 
rnbinom 1e6: mean = 29.981  var = 74.9131 
```

`size = 1/α, mu = μ` and `size = r, prob = r/(r+μ)` are the same distribution. The sample variance of $10^6$ random draws, 74.9, also agrees well with the theoretical 75.

This time I compute the integral above numerically and compare it with `dnbinom`.

```r
gp <- function(k) integrate(function(l) dpois(k, l) * dgamma(l, shape = 1/alpha, scale = alpha*mu), 0, Inf, rel.tol = 1e-10, abs.tol = 0)$value
print(cbind(k = k, integral = sapply(k, gp), dnbinom = dnbinom(k, size = 1/alpha, mu = mu)), digits = 10)
cat("mean of gamma = shape*scale =", (1/alpha)*(alpha*mu), " variance = shape*scale^2 =", (1/alpha)*(alpha*mu)^2, " = alpha*mu^2 =", alpha*mu^2, "\n")
```

```
      k        integral         dnbinom
[1,]  0 1.099511628e-08 1.099511628e-08
[2,] 10 1.331660435e-03 1.331660435e-03
[3,] 20 2.770707453e-02 2.770707453e-02
[4,] 30 4.582342113e-02 4.582342113e-02
[5,] 40 2.053684624e-02 2.053684624e-02
[6,] 60 4.749447753e-04 4.749447753e-04
mean of gamma = shape*scale = 30  variance = shape*scale^2 = 45  = alpha*mu^2 = 45 
```

They agree to 10 digits. With the default tolerance of `integrate()` (an absolute error of about $10^{-4}$), $k=0$ (on the order of $10^{-8}$) and the tail at $k=60$ did not match, so I gave `rel.tol = 1e-10, abs.tol = 0`.

</details>

<details>
<summary>"One α" seen across many genes</summary>

In a simulation with large condition effects (`betaSD = 2`), I computed the gene-wise dispersions on the same counts, changing only the design.

```r
set.seed(1)
d24 <- estimateSizeFactors(makeExampleDESeqDataSet(n = 1000, m = 6, betaSD = 2))   # simulation with large condition effects
big <- abs(mcols(d24)$trueBeta) > 2                                                # genes with |true log2FC| above 2
ge <- function(des) { design(d24) <- des; mcols(estimateDispersionsGeneEst(d24, quiet = TRUE))$dispGeneEst[big] }
cat("number of genes:", sum(big), "\n")
cat("median dispGeneEst  ~condition:", median(ge(~ condition), na.rm = TRUE), "  ~1:", median(ge(~ 1), na.rm = TRUE),
    "  true dispersion (trueDisp):", median(mcols(d24)$trueDisp[big]), "\n")
```

```
number of genes: 332 
median dispGeneEst  ~condition: 0.187612   ~1: 1.31605   true dispersion (trueDisp): 0.370905 
```

With condition as the mean structure (`~condition`), the median gene-wise dispersion is 0.19; putting all samples around one mean (`~1`) makes it 1.32, about 7 times larger. The difference in means between the conditions leaked into $\alpha$. That the 0.19 of `~condition` is below the true median of 0.37 is a property of gene-wise estimation from 6 samples. This underestimation, the details of how the design affects dispersion estimation, and the implementation of the gene-wise step, which optimizes $\alpha$ with the mean fixed, are all covered in note 03.

</details>

<details>
<summary>Size factor details: the dds method, the geometric mean, the priority of normalization factors</summary>

First, does computing the study text's table from section 4 with the dds method give the same values?

```r
cat("all.equal(by hand, estimateSizeFactorsForMatrix):", all.equal(unname(s_hat), unname(estimateSizeFactorsForMatrix(K))), "\n")
dds3 <- DESeqDataSetFromMatrix(K, data.frame(row.names = colnames(K), cond = factor(c("a","b"))), ~cond)
dds3 <- estimateSizeFactors(dds3)
cat("sizeFactors(dds3):", sizeFactors(dds3), "\n"); print(counts(dds3, normalized = TRUE))
cat("colSums ratio:", colSums(K)/mean(colSums(K)), "\n")
```

```
all.equal(by hand, estimateSizeFactorsForMatrix): TRUE 
converting counts to integer mode
sizeFactors(dds3): 0.707107 1.41421 
    sample1  sample2
g1 141.4214 141.4214
g2  70.7107  70.7107
g3 282.8427 282.8427
colSums ratio: 0.666667 1.33333 
```

The study text's 141.42 / 70.71 / 282.84 / 0.707 / 1.414 come out identically from the hand calculation, the matrix function and the dds method, and the normalized matrix is equal within each row. In this toy data the colSums ratio (0.667, 1.333) is also 2-fold, so the result looks the same as total-count normalization, but changing just one gene, like g4 in section 4, makes them diverge.

How different are `sizeFactors(dds)` and the colSums ratio in a simulation with 2,000 genes?

```r
set.seed(42)
dds <- makeExampleDESeqDataSet(n = 2000, m = 6, sizeFactors = c(0.5, 1, 1.5, 0.8, 1.2, 2))
dds <- estimateSizeFactors(dds)
cs <- colSums(counts(dds)); sf <- sizeFactors(dds); gm <- function(x) x/exp(mean(log(x)))
out <- data.frame(true_sf = c(0.5,1,1.5,0.8,1.2,2), sizeFactors = sf, sf_scaled = gm(sf), colSums = cs, cs_scaled = gm(cs), ratio = gm(sf)/gm(cs))
print(out, digits = 4)
cat("max |sf_scaled/cs_scaled - 1| =", max(abs(out$ratio - 1)), "\n")
cat("identical(sizeFactors, colSums/mean(colSums))?", isTRUE(all.equal(unname(sf), unname(cs/mean(cs)))), "\n")
cat("geometric mean of sizeFactors(dds):", exp(mean(log(sf))), "\n")
cat("identity check sizeFactors == estimateSizeFactorsForMatrix(counts):", all.equal(sf, estimateSizeFactorsForMatrix(counts(dds))), "\n")
sf_g <- estimateSizeFactorsForMatrix(counts(dds), geoMeans = exp(rowMeans(log(counts(dds)))))   # pass the same reference as geoMeans
cat("passing geoMeans directly: same ratios?", all.equal(sf_g/sf_g[1], sf/sf[1]), " geometric mean =", exp(mean(log(sf_g))), "\n")
```

```
        true_sf sizeFactors sf_scaled colSums cs_scaled  ratio
sample1     0.5      0.4919    0.4741   40007    0.4566 1.0384
sample2     1.0      0.9533    0.9188   80379    0.9174 1.0016
sample3     1.5      1.4417    1.3896  126849    1.4478 0.9598
sample4     0.8      0.7880    0.7595   65961    0.7528 1.0088
sample5     1.2      1.2008    1.1574  100190    1.1435 1.0121
sample6     2.0      1.9498    1.8793  167833    1.9155 0.9811
max |sf_scaled/cs_scaled - 1| = 0.0401773 
identical(sizeFactors, colSums/mean(colSums))? FALSE 
geometric mean of sizeFactors(dds): 1.03751 
identity check sizeFactors == estimateSizeFactorsForMatrix(counts): TRUE 
passing geoMeans directly: same ratios? TRUE  geometric mean = 1 
```

Even with both scaled to a geometric mean of 1, they differ by up to 4 %. This simulation has no DE genes (`betaSD = 0`) and no extreme genes, so the difference is on the small side. When a few genes take a large share of the total reads, as in real data, they diverge further, like the g4 example in section 4. `sizeFactors(dds)` is exactly `estimateSizeFactorsForMatrix(counts(dds))`, and its geometric mean is not 1 (1.0375). Passing the same reference directly as `geoMeans =` leaves the ratios unchanged and rescales only the geometric mean to 1. So what decides the rescaling is not the type but whether the geoMeans argument came in.

What happens to the size factors when normalization factors are supplied?

```r
dds_nf <- dds; nf <- matrix(1, nrow(dds), ncol(dds)); normalizationFactors(dds_nf) <- nf
cat("sizeFactors still there:", round(sizeFactors(dds_nf), 4), "\n")
cat("counts(normalized=TRUE)[1,] == counts[1,] (nf=1 is used):", all.equal(counts(dds_nf, normalized = TRUE)[1,], counts(dds_nf)[1,]), "\n")
cat("getSizeOrNormFactors == nf:", all.equal(unname(DESeq2:::getSizeOrNormFactors(dds_nf)), nf), "\n")
cat("without nf, getSizeOrNormFactors copies sf into every row:", all.equal(DESeq2:::getSizeOrNormFactors(dds)[1:2,], matrix(rep(sf, each = 2), 2)), "\n")
cat("DESeq() message: "); invisible(capture.output(tmp <- DESeq(dds_nf[1:50,], quiet = FALSE), type = "output"))
cat("setting normalizationFactors to 0: "); print(tryCatch({normalizationFactors(dds_nf) <- nf*0; "ok"}, error = function(e) conditionMessage(e)))
```

```
sizeFactors still there: 0.4919 0.9533 1.4417 0.788 1.2008 1.9498 
counts(normalized=TRUE)[1,] == counts[1,] (nf=1 is used): TRUE 
getSizeOrNormFactors == nf: TRUE 
without nf, getSizeOrNormFactors copies sf into every row: TRUE 
DESeq() message: using pre-existing normalization factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
setting normalizationFactors to 0: [1] "all(value > 0) is not TRUE"
```

Even though the `sizeFactors` slot is still there, the normalized counts, the GLM offset and the `DESeq()` message all follow the normalization factors. As the rejection of 0 shows, normalization factors must be positive finite values (the `stopifnot` in `normalizationFactors<-`).

</details>

<details>
<summary>Counts of 0: the exclusion rule and poscounts</summary>

I tried a matrix with one 0 each in g2 and g4.

```r
K0 <- matrix(c(100,200, 0,100, 200,400, 50,0), nrow = 4, byrow = TRUE,
             dimnames = list(c("g1","g2","g3","g4"), c("s1","s2")))     # one 0 each in g2 and g4
cat("rowMeans(log(K0)):", rowMeans(log(K0)), "\n")
cat("ratio:", estimateSizeFactorsForMatrix(K0), " (only g1, g3 used -> by hand:", exp(median((log(K0[,1]) - rowMeans(log(K0)))[c(1,3)])), ")\n")
Kall0 <- K0[c("g2","g4"), ]
cat("one 0 in every gene: "); print(tryCatch(estimateSizeFactorsForMatrix(Kall0), error = function(e) conditionMessage(e)))
sf_mat <- estimateSizeFactorsForMatrix(Kall0, type = "poscounts"); cat("poscounts (same matrix):", sf_mat, "\n")
lc <- log(Kall0); lc[!is.finite(lc)] <- 0
cat("poscounts by hand: loggeomeans =", rowMeans(lc), "; s1 =", exp(median((log(Kall0[,1]) - rowMeans(lc))[Kall0[,1] > 0])),
    " s2 =", exp(median((log(Kall0[,2]) - rowMeans(lc))[Kall0[,2] > 0])), "\n")
```

```
rowMeans(log(K0)): 4.95174 -Inf 5.64489 -Inf 
ratio: 0.707107 1.41421  (only g1, g3 used -> by hand: 0.707107 )
one 0 in every gene: [1] "every gene contains at least one zero, cannot compute log geometric means"
poscounts (same matrix): 7.07107 10 
poscounts by hand: loggeomeans = 2.30259 1.95601 ; s1 = 7.07107  s2 = 10 
```

g2 and g4 have a geometric mean of 0 (a $\log$ mean of $-\infty$), so they drop out of the ratio calculation, and 0.707/1.414 comes from g1 and g3 alone. Keeping only g2 and g4 gives an error. `poscounts` uses as its reference a "shrunken" geometric mean that sets $\log 0$ to 0 and divides by the total number of samples $n=2$ ($e^{2.30}=10$, $e^{1.96}\approx7.07$). A geometric mean of only the values that are not 0 would have been 100 and 50.

But computing the same matrix with the dds method gives different values.

```r
dd0 <- DESeqDataSetFromMatrix(Kall0, data.frame(row.names = c("s1","s2"), c = factor(c("a","b"))), ~c)
sf_dds <- sizeFactors(estimateSizeFactors(dd0, type = "poscounts", quiet = TRUE))
cat("matrix poscounts:", sf_mat, " geomean =", exp(mean(log(sf_mat))), "\n")
cat("dds    poscounts:", sf_dds, " geomean =", exp(mean(log(sf_dds))), "\n")
cat("ratio s2/s1 is the same:", sf_mat[2]/sf_mat[1], sf_dds[2]/sf_dds[1], "; dds value == matrix value/geomean:",
    all.equal(unname(sf_dds), unname(sf_mat/exp(mean(log(sf_mat))))), "\n")
```

```
converting counts to integer mode
matrix poscounts: 7.07107 10  geomean = 8.40896 
dds    poscounts: 0.840896 1.18921  geomean = 1 
ratio s2/s1 is the same: 1.41421 1.41421 ; dds value == matrix value/geomean: TRUE 
```

`estimateSizeFactors(dds, type = "poscounts")` builds and passes `geoMeans`, so `incomingGeoMeans = TRUE` and it is rescaled to a geometric mean of 1, while `estimateSizeFactorsForMatrix(type = "poscounts")` does not rescale. The ratios between samples are the same, so with the dispersions fixed only the intercept moves (section 5). But running `DESeq()` from scratch, the starting value of the gene-wise dispersion depends on the scale of the size factors ("Changing the scale of the size factors and rerunning DESeq()"), so the dispersions, p-values and padj can differ a little. It is also easy to get confused when comparing the values by eye.

</details>

<details>
<summary>Importing with tximport and tximeta</summary>

tximport is not installed in my environment, so I built a list with the same shape as the output of `tximport()` by hand and passed it to `DESeqDataSetFromTximport()`.

```r
cat("DESeqDataSetFromTximport lines 4-5:", deparse(DESeqDataSetFromTximport)[4:5], "\n")
txi <- list(counts = matrix(c(10, 20, 30, 40), 2, 2, dimnames = dimnames(tpm_like)), length = matrix(1000, 2, 2), countsFromAbundance = "no")
cat("countsFromAbundance = 'no'        -> assayNames:", assayNames(DESeqDataSetFromTximport(txi, cd, ~cond)), "\n")
txi$countsFromAbundance <- "scaledTPM"
cat("countsFromAbundance = 'scaledTPM' -> assayNames:", assayNames(DESeqDataSetFromTximport(txi, cd, ~cond)), "\n")
```

```
DESeqDataSetFromTximport lines 4-5:     counts <- round(txi$counts)     mode(counts) <- "integer" 
using counts and average transcript lengths from tximport
countsFromAbundance = 'no'        -> assayNames: counts avgTxLength 
using just counts from tximport
countsFromAbundance = 'scaledTPM' -> assayNames: counts 
```

`estimateSizeFactors` goes through the normalization factor path (lines 27–31 of `estimateSizeFactors`) only if the `avgTxLength` assay exists. Counts made with `scaledTPM`/`lengthScaledTPM` have no such assay, so they are handled with ordinary size factors.

That does not mean the length correction is missing. In the tximport source (`makeCountsFromAbundance`), `scaledTPM` is TPM (already divided by effective length) scaled up to the library size, and `lengthScaledTPM` is TPM multiplied by each gene's average length over samples and then put on the same scale. Both build the effect of lengths that differ between samples into the counts in advance, so no offset is needed. But I only read this part of the tximport source and did not run it.

On the tximeta path, `DESeqDataSet()` sees `tximetaInfo` in `metadata(se)` and handles it with `processTximeta` (lines 20–22 of `DESeqDataSet` in the source block).

</details>

<details>
<summary>Changing the scale of the size factors and rerunning DESeq()</summary>

First, what happens if it goes in as a covariate instead of an offset? In `d1` from section 5, I set the size factors to 1 and put $\log_2 s_j$ into the design as a column.

```r
d3 <- d1; d3$logsf <- log2(sizeFactors(d1)); sizeFactors(d3) <- rep(1, 6)
design(d3) <- ~ logsf + condition; d3 <- DESeq(d3, quiet = TRUE)
print(summary(coef(d3)[, "logsf"]))
cat("resultsNames(d3):", resultsNames(d3), "\n")
```

```
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 -2.931   0.513   1.021   1.013   1.583   3.967       1 
resultsNames(d3): Intercept logsf condition_B_vs_A 
```

The median is near 1, but it scatters from −2.9 to 4.0 across genes. This experiment shows from the opposite direction that an offset is not "one more estimated coefficient".

This time I doubled the size factors and then re-estimated the dispersions as well.

```r
d2r <- d1; sizeFactors(d2r) <- 2 * sizeFactors(d1); d2r <- DESeq(d2r, quiet = TRUE)   # re-estimate the dispersions
cat("DESeq(d2r) range of intercept differences:", range(coef(d2r)[,1] - coef(d1)[,1], na.rm = TRUE), " median:", median(coef(d2r)[,1] - coef(d1)[,1], na.rm = TRUE), "\n")
dc <- coef(d2r)[,2] - coef(d1)[,2]
cat("max |difference| of condition coefficient:", max(abs(dc), na.rm = TRUE), "; genes changed by more than 1e-3:", sum(abs(dc) > 1e-3, na.rm = TRUE), "/", sum(!is.na(dc)),
    "; number with padj<0.1 d1, d2r:", sum(results(d1)$padj < 0.1, na.rm = TRUE), sum(results(d2r)$padj < 0.1, na.rm = TRUE), "\n")
g1 <- mcols(d1)$dispGeneEst; g2 <- mcols(d2r)$dispGeneEst; chg <- which(abs(g2/g1 - 1) > 1e-3)
cat("genes whose gene-wise dispersion changed by 0.1 % or more:", length(chg), "/", sum(!is.na(g1)), "; of these collapsed to 1e-8 (minDisp):", sum(g2[chg] <= 1e-7), "\n")
cat("trend coefficients (a0, a1)  d1:", attr(dispersionFunction(d1), "coefficients"), "  d2r:", attr(dispersionFunction(d2r), "coefficients"), "\n")
```

```
DESeq(d2r) range of intercept differences: -1.01175 -0.970659  median: -0.999992 
max |difference| of condition coefficient: 0.0298464 ; genes changed by more than 1e-3: 193 / 299 ; number with padj<0.1 d1, d2r: 53 49 
genes whose gene-wise dispersion changed by 0.1 % or more: 61 / 299 ; of these collapsed to 1e-8 (minDisp): 48 
trend coefficients (a0, a1)  d1: 0.115638 5.47721   d2r: 0.108227 3.46123 
```

The intercept shift scatters from −1.012 to −0.971 (up to 3 %). The median is −1.0000, so most genes are unaffected, but the condition coefficient also changes by up to 0.03 (193 of 299 genes by more than $10^{-3}$), and the number of genes with padj < 0.1 drops from 53 to 49. So "the scale of the size factors does not matter" is true only when the dispersions are fixed.

What changed is the dispersion. The gene-wise estimates change for 61 of 299 genes, and 48 of them collapse to `minDisp` = $10^{-8}$. The trend coefficients change too, and $a_1$ does not simply halve. I traced at which step the scale leaks in with the next block.

```r
gs <- deparse(estimateDispersionsGeneEst); fs <- deparse(estimateDispersionsFit)
cat(sprintf("[GeneEst %d] %s", 40:43, gs[40:43]), sprintf("[Fit %d] %s", 8, fs[8]), sep = "\n")   # print part of the source
nz1 <- d1[!mcols(d1)$allZero, ]; nz2 <- d2r[!mcols(d2r)$allZero, ]; X <- model.matrix(design(d1), colData(d1))
m1 <- DESeq2:::momentsDispEstimate(nz1); m2 <- DESeq2:::momentsDispEstimate(nz2)
r1 <- DESeq2:::roughDispEstimate(counts(nz1, normalized = TRUE), X); r2 <- DESeq2:::roughDispEstimate(counts(nz2, normalized = TRUE), X)
cat("momentsDispEstimate max relative difference:", max(abs(m2/m1 - 1)), "; genes with roughDispEstimate = 0 in d1, d2r:", sum(r1 == 0), sum(r2 == 0), "\n")
col <- match(rownames(d1)[chg][g2[chg] <= 1e-7], rownames(nz1))
cat("collapsed genes with roughDispEstimate = 0: d2r", sum(r2[col] == 0), "/", length(col), ", d1", sum(r1[col] == 0), "; their max d1 dispGeneEst:", max(g1[chg][g2[chg] <= 1e-7]), "\n")
d2i <- estimateDispersionsGeneEst(d2r, alphaInit = pmin(r1, m1), quiet = TRUE)   # give d2r the starting values of d1
cat("estimating d2r from d1's starting values, max relative difference of dispGeneEst:", max(abs(mcols(d2i)$dispGeneEst/g1 - 1), na.rm = TRUE), "\n")
```

```
[GeneEst 40]         roughDisp <- roughDispEstimate(y = counts(objectNZ, normalized = TRUE), 
[GeneEst 41]             x = modelMatrix)
[GeneEst 42]         momentsDisp <- momentsDispEstimate(objectNZ)
[GeneEst 43]         alpha_hat <- pmin(roughDisp, momentsDisp)
[Fit 8]     useForFit <- mcols(objectNZ)$dispGeneEst > 100 * minDisp
momentsDispEstimate max relative difference: 0 ; genes with roughDispEstimate = 0 in d1, d2r: 38 87 
collapsed genes with roughDispEstimate = 0: d2r 48 / 48 , d1 0 ; their max d1 dispGeneEst: 10 
estimating d2r from d1's starting values, max relative difference of dispGeneEst: 0 
```

The starting value of the gene-wise optimization is `pmin(roughDisp, momentsDisp)` (GeneEst lines 40–43). Of the two, `momentsDispEstimate` is $(\text{baseVar}-\bar{x}\,\text{baseMean})/\text{baseMean}^2$ with $\bar x=$`mean(1/s_j)`. Under $s\to2s$, baseMean and $\bar x$ halve and baseVar becomes 1/4, so it stays exactly the same (maximum difference 0).

The cause is on the `roughDispEstimate` side. This function treats the normalized counts $y$ as if they were raw counts, computes $\frac{1}{m-p}\sum_j\frac{(y_j-\hat\mu_j)^2-\hat\mu_j}{\hat\mu_j^2}$ and clips negative values to 0 ($\hat\mu$ also has a `pmax(1, mu)` lower bound). When $y$ halves, the first term stays the same but the Poisson term $-1/\hat\mu_j$ doubles. So the number of genes clipped to 0 grows from 38 to 87.

All 48 collapsed genes are genes whose rough estimate is 0 in d2r, so their starting value was clipped to `minDisp` (0 such genes in d1). They include a gene whose d1 estimate was 10, so these are not originally small values that got rounded. Passing d1's starting values as `alphaInit` makes d2r's gene-wise estimates exactly equal to d1's (difference 0). That means the whole difference comes from the starting values. Because the collapsed genes are left out of the trend fit (Fit line 8, `useForFit`), the trend coefficients change along with them.

In short, the cause is not the offset (the GLM step, $10^{-5}$ in section 5) but the starting value of the dispersion estimation. The details of that step are covered in note 03.

</details>

<details>
<summary>When some or most genes move: controlGenes</summary>

First, the setup of the simulation in section 6 again. A pair effect (p1 ×1, p2 ×2, p3 ×4) was multiplied into genes 1–200 only, out of 1,000. An effect present in only some genes is not library size, so it must not disappear into the size factors. I made it with no condition effect (`betaSD = 0`), and since `makeExampleDESeqDataSet` makes condition A A A B B B, I relabelled it so that each pair contains one A and one B. The true condition effect is 0 for every gene, so relabelling leaves the truth unchanged. Therefore the true values of the pair and condition coefficients of genes 201–1000 are 0 for every gene.

In the result of section 6, the share absorbed by the size factors is the mean $s_j$ of p3, 1.18, against that of p1, 0.89, that is, $\log_2(1.177/0.894)\approx0.40$. Looking at the difference between the two groups of genes (0.67 + 0.15 ≈ 0.82, 1.65 + 0.36 ≈ 2.01), p3 is close to the theoretical 2. p2 stays at 0.85 even after the reference is corrected below, so the median of the 200 genes itself seems to have moved. The same mechanism already shows at 20 %, weaker than the condition in section 3.3 of the study text ("most genes").

Following the study text's prescription, I fix the reference with control genes.

```r
dpc <- estimateSizeFactors(dp, controlGenes = 201:1000)   # only genes without a pair effect as the reference
cat("controlGenes=201:1000 sizeFactors:", round(sizeFactors(dpc), 3), "\n")
dpc <- DESeq(dpc, quiet = TRUE)
cat("median pair coefficient, gene 1-200:", round(apply(coef(dpc)[1:200, 2:3], 2, median, na.rm = TRUE), 3),
    " / gene 201-1000:", round(apply(coef(dpc)[201:1000, 2:3], 2, median, na.rm = TRUE), 3), "\n")
```

```
controlGenes=201:1000 sizeFactors: 0.996 1.015 1.024 1.047 1.037 1.058 
median pair coefficient, gene 1-200: 0.846 1.99  / gene 201-1000: 0.026 -0.024 
```

Taking the median only over unaffected genes makes the size factors independent of the pair (all ≈ 1). The medians of the pair coefficients also come back toward their theoretical values, 0.85/1.99 and 0.03/−0.02. In real data, though, you do not know which genes are "unaffected", so spike-ins or prior knowledge are needed. This is exactly what the study text is saying.

What happens, then, if most (80 %) of the genes move in one direction?

```r
dp8 <- dp; cnt8 <- counts(dp); cnt8[201:800, ] <- round(t(t(cnt8[201:800, ]) * eff)); storage.mode(cnt8) <- "integer"
counts(dp8) <- cnt8; dp8 <- estimateSizeFactors(dp8)   # now a pair effect in 1-800 (80 %)
lr <- function(d) round(log2(mean(sizeFactors(d)[5:6])/mean(sizeFactors(d)[1:2])), 3)
cat("20 % moved: sizeFactors", round(sizeFactors(dp), 3), " log2(s_p3/s_p1) =", lr(dp), "\n")
cat("80 % moved: sizeFactors", round(sizeFactors(dp8), 3), " log2(s_p3/s_p1) =", lr(dp8), " (2 if the whole effect were absorbed into the size factors)\n")
nc8 <- counts(dp8, normalized = TRUE); pm8 <- sapply(split(seq_len(6), dp8$pair), function(j) rowMeans(nc8[801:1000, j, drop = FALSE]))
g8 <- exp(colMeans(log(pm8 + 0.5)))   # same geometric-mean measure as before
cat("unaffected gene 801-1000, per-pair geometric mean of normalized counts:", round(g8, 2), " p1/p3 =", round(g8[1]/g8[3], 2), "; s_p3/s_p1 =", round(2^lr(dp8), 2), "\n")
```

```
20 % moved: sizeFactors 0.877 0.91 1.028 1.051 1.158 1.196  log2(s_p3/s_p1) = 0.397 
80 % moved: sizeFactors 0.559 0.579 1.028 1.051 1.816 1.887  log2(s_p3/s_p1) = 1.702  (2 if the whole effect were absorbed into the size factors)
unaffected gene 801-1000, per-pair geometric mean of normalized counts: 25.61 15.47 9.39  p1/p3 = 2.73 ; s_p3/s_p1 = 3.25 
```

When 80 % of the genes moved 4-fold, the median moved inside that group, and 1.70 of the pair effect of log2 2 became size factor. The "invariant reference" assumption has collapsed. The unaffected 20 % of genes (801–1000) then look after normalization as if they dropped about 2.7-fold from p1 → p3 (25.6 → 9.4 by the same geometric-mean measure as in section 6). In terms of the size factor ratio it is $2^{1.70}\approx3.25$-fold. The warning on p.8 of the study text shows up exactly.

</details>

<details>
<summary>Common misconceptions</summary>

| Misconception | In fact |
|---|---|
| A large dispersion means a large variance | True at the same $\mu$ (section 2: at $\mu=30$ the variance goes 30 → 75 → 210). But not if $\mu$ differs. The variance is $\mu+\alpha\mu^2$, so the variance at $\alpha=0.1,\mu=100$ (1,100) is smaller than at $\alpha=0.01,\mu=1000$ (11,000). $\alpha$ is the floor of $CV^2$. |
| Low-expression genes are a problem because their variance is large | Their absolute variance is actually small. In $CV^2=1/\mu+\alpha$, $1/\mu$ is large, so the relative variation is large and that is why the estimates are unstable. |
| Writing `dnbinom(size=α)` is fine | `size` is $1/\alpha$. The correspondence in section 2 holds only for `size=1/alpha, mu=mu`. |
| The size factor is the ratio of total read counts | It is median-of-ratios. In "Size factor details" the two differed by up to 4 % even after matching their geometric means, and in the extreme-gene example of section 4 by about 7.7-fold (ratio 2 vs 15.3). |
| Size factors are normalized to a geometric mean of 1 | Not on the default ratio path (1.0375). They are rescaled only when `geoMeans` is passed in (the poscounts path of the dds, or passed directly). |
| If every gene has a 0, add a pseudocount | `type="poscounts"` uses a shrunken geometric mean that sets $\log 0$ to 0 and divides by $n$. The input is left untouched. |
| The median holds up even if most genes move | Even with only 20 % moving, part of it is absorbed (0.40), and with 80 % most of it is (1.70/2). In that case `controlGenes=` or spike-ins are needed. |
| Normalized counts are fed back into DESeq2 | The model fits the raw counts and the size factor enters as an offset. Normalized counts are for visualization and summaries, and they are not integers either (rejected in the input block of section 4). |
| With `~pair+condition`, `counts(normalized=TRUE)` removes the pair | It does not (section 6). The GLM coefficients absorb the pair, and to remove it for a figure use `removeBatchEffect`. |
| `DESeq()` re-estimates and overwrites size factors set in advance | `DESeq()` uses existing sizeFactors/normalizationFactors as they are ("using pre-existing …"). |
| Even with normalizationFactors supplied, sizeFactors are used as well | If normalizationFactors exist, sizeFactors are ignored (`getSizeOrNormFactors`, "Size factor details"). |
| With any tximport option, DESeq2 applies a length offset | The length offset (`avgTxLength` → normalization factor) applies only with `countsFromAbundance="no"`. `scaledTPM`/`lengthScaledTPM` build the lengths into the counts in advance, so DESeq2 uses only ordinary size factors (from reading the tximport source; not run). |

</details>

<details>
<summary>How this connects to the next notes</summary>

I said $\alpha_i$ is one value per gene, but I have not yet said how that value is obtained. [Note 03](03_dispersion_estimation.md) covers the gene-wise MLE, the Cox-Reid adjustment, the trend and MAP shrinkage.

Changing the scale of the size factors made some gene-wise dispersion estimates collapse to `minDisp` and changed the trend coefficients. I traced the cause as far as the `roughDispEstimate` side (which treats normalized counts as raw counts) of the starting value of the gene-wise optimization, `pmin(roughDispEstimate, momentsDispEstimate)`. It comes up again when note 03 looks at the starting values and the optimization path.

How $\beta$ and the SE come out with the offset fixed is shown in [note 04](04_glm_condition_batch.md). Two things seen here, that the matrix built by `getSizeOrNormFactors` goes into `nfSEXP` of `fitBetaWrapper`, and that `DESeq()` does not overwrite existing size factors, are assumed in the reproduction experiments of notes 03 and 04.

</details>

---

← Previous: [01. Poisson and overdispersion](01_poisson_simulation.md) · Next: [03. Dispersion estimation](03_dispersion_estimation.md) →

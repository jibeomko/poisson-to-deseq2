# 01. Why do counts differ so much between replicates of the same condition?

I counted the reads of one gene in three replicates grown under the same condition, and the counts came out fairly far apart: 90, 100 and 130. Can the chance involved in counting reads alone make them wobble this much? In this note I first build the model that contains only that chance (the Poisson), and then follow by simulation why the variance becomes $\mu + \alpha\mu^2$ when the expression level itself wobbles from sample to sample. It helps to look at the flow of the whole series in [00. The overall map](00_overview.md) first.

> Study text 2.1–2.2 · Code: [01_poisson_simulation.ipynb](../../01_poisson_simulation.ipynb)

## 1. How much do the counts of gene A wobble?

This series follows one gene from beginning to end, gene A from chapter 16 of the study text. These are the reads counted in a control grown in ordinary medium (Ctrl) and a group grown in medium without glucose (Starvation), with three replicates each.

| Group | replicate 1 | replicate 2 | replicate 3 | Mean |
|---|---:|---:|---:|---:|
| Ctrl | 100 | 130 | 90 | 106.667 |
| Starvation | 200 | 250 | 180 | 210 |

Each number in the table is a count. The count $K_{ij}$ is the number of reads (or fragments) assigned to gene $i$ in sample $j$; since this note looks at only one gene, I drop the $i$ and write $K_j$. The size factor $s_j$, the multiplier that corrects for differences in sequencing depth between samples, is 1 for all six samples. In effect the depths are taken to be equal. Size factors are covered in detail in [note 02](02_negative_binomial.md).

The three Ctrl samples were grown under exactly the same condition, yet their counts range from 90 to 130. Why do they differ like this? The causes can be split into two broad kinds.

The first is chance in the counting process. Sequencing the same RNA library again does not give exactly the same number of reads. Which of the tens of millions of reads come from this gene is decided at random.

The second is biological difference. The three Ctrl samples are different biological replicates. If, for example, the donors or culture dishes differ, the cell state, tissue composition and preparation all differ a little, so the actual expression level of this gene also differs from sample to sample. Section 2.2 of the study text adds that "heterogeneity not explained by the covariates put in the model also remains". That is, some difference is left even after explaining what the factors in the design, such as condition or batch, can explain.

I first build the model with only the first cause (the Poisson), then see what changes when the second cause is added. Before that, here is how spread out the two groups are, in numbers (SD is the standard deviation).

```r
geneA <- list(Ctrl = c(100, 130, 90), Starvation = c(200, 250, 180))
for (grp in names(geneA)) {
  k <- geneA[[grp]]
  cat(sprintf("%-6s mean %.3f | sample var %6.1f | var/mean %.2f | Poisson SD %.1f vs actual SD %.1f\n",
              grp, mean(k), var(k), var(k) / mean(k), sqrt(mean(k)), sqrt(var(k))))
}
```

```text
Ctrl   mean 106.667 | sample var  433.3 | var/mean 4.06 | Poisson SD 10.3 vs actual SD 20.8
Starvation mean 210.000 | sample var 1300.0 | var/mean 6.19 | Poisson SD 14.5 vs actual SD 36.1
```

The column to watch is "var/mean". As sections 2 and 3 will show, this value should be near 1 if the Poisson is right, but for gene A it is 4–6 in both groups. The last two numbers are the standard deviation expected under a Poisson ($\sqrt{\text{mean}}$) and the actual standard deviation.

## 2. How would counts wobble if reads were only being counted?

Say one sample was sequenced and gave 20,000,000 reads, of which a little over 100 come from gene A. The probability that a given read comes from gene A is tiny, but the number of reads is huge. When a trial with a small success probability is repeated a very large number of times and the successes are counted, that count is well approximated by a Poisson distribution. This property is called the Poisson limit of the binomial (an explanation not in the study text).

$$K_j \sim \mathrm{Poisson}(\mu_j), \qquad E(K_j) = \mu_j, \qquad Var(K_j) = \mu_j$$

$K_j$ is the count of sample $j$, and $\mu_j$ is the mean count the model expects in that sample, the expected count. Think of $E(K_j)$ and $Var(K_j)$ as the mean and variance of the count if the same library were sequenced and counted again and again without end.

So the Poisson makes one promise: the variance equals the mean. For gene A's Ctrl, $\mu = 106.667$, so the variance is also 106.667 and the standard deviation is $\sqrt{106.667} = 10.3$. But the actual standard deviation of the three values is 20.8, about twice as large.

Section 2.1 of the study text points out one thing here. The Poisson is not a distribution newly invented for sequencing; an existing probability distribution is borrowed as an approximate model of this observation process. So it cannot capture all the complications of how reads are really assigned, such as reads that align to several places (multi-mapping), GC bias and PCR duplicates. Conversely, it has been observed that in technical replicates, where the same library is sequenced again, the count variation is close to Poisson (Marioni et al. 2008, *Genome Res.* 18:1509). These examples and the observation of Marioni et al. are additions not in the study text.

How common, then, is a standard deviation about twice the Poisson one, as in Ctrl, under a Poisson? I drew three values from a Poisson with the same mean as gene A, repeating it 100,000 times.

```r
# If the Poisson is right, how often is the var/mean of three replicates as large as gene A's?
set.seed(1)
for (grp in c("Ctrl", "Starvation")) {
  mu  <- c(Ctrl = 106.667, Starvation = 210)[grp]
  obs <- c(Ctrl = 433.333 / 106.667, Starvation = 1300 / 210)[grp]
  K   <- matrix(rpois(3 * 1e5, mu), ncol = 3)
  r   <- apply(K, 1, var) / rowMeans(K)
  cat(sprintf("%-6s observed ratio %.2f | share at least as large under Poisson %.4f\n", grp, obs, mean(r >= obs)))
}
```

```text
Ctrl   observed ratio 4.06 | share at least as large under Poisson 0.0168
Starvation observed ratio 6.19 | share at least as large under Poisson 0.0020
```

The last numbers show that if the Poisson is right, a spread as large as Ctrl's happens only 1.7% of the time, and one as large as Starvation's only 0.2%. These numbers are of the same nature as a p-value computed assuming "the Poisson is right". Gene A wobbles more than the Poisson predicts. But this evidence comes from one gene and three replicates, so how much it can be trusted on its own is taken up again in section 5.

## 3. Why is "the variance grows with the mean" not evidence?

The Starvation group has a mean of 210, about twice that of Ctrl, and its sample variance has also grown, from 433 to 1,300. Seeing this, one is tempted to say "the variance is larger where the mean is larger, so the Poisson won't do". I thought so at first too. But section 2.2 of the study text blocks exactly this inference, citing the DESeq2 paper (Love et al. 2014, reference [2] of the study text) as the source of the point.

Think about it: under the Poisson the variance equals the mean, so when the mean goes from 100 to 200, the variance goes from 100 to 200 as well. A variance that grows with the mean is exactly what the Poisson predicts. What the Poisson cannot explain is the part that exceeds the variance it predicts, the part where $Var(K) - \mu > 0$. So the thing to look at is not the size of the variance but **the ratio of the variance to the mean**.

$$\frac{Var(K_j)}{\mu_j} = 1 \;\;(\text{Poisson}), \qquad \frac{Var(K_j)}{\mu_j} > 1 \;\;(\text{overdispersion})$$

When this ratio is above 1, that is, when the spread between replicates of the same condition is larger than the Poisson predicts, it is called overdispersion. For gene A it is 4.06 in Ctrl and 6.19 in Starvation, above 1 in both groups.

There is one thing to be careful about. The $\mu_j$ in the denominator is the expected count of sample $j$, so it is set separately for each sample, reflecting its condition (Ctrl or Starvation) and its size factor. If samples with different means are pooled to compute the ratio, the ratio comes out far above 1 even when the counts are purely Poisson.

```r
# Pool two groups with different means and compute var/mean? (both pure Poisson)
set.seed(1)
k100 <- rpois(1e5, 100); k200 <- rpois(1e5, 200)
ratio <- function(k) var(k) / mean(k)
cat(sprintf("mu=100 only: %.2f | mu=200 only: %.2f | pooled: %.2f\n",
            ratio(k100), ratio(k200), ratio(c(k100, k200))))
```

```text
mu=100 only: 1.00 | mu=200 only: 0.99 | pooled: 17.66
```

Looked at separately, the two groups have a ratio of 1, but pooled it becomes 17.66. The difference between the group means (100 and 200) has leaked into the variance. So overdispersion has to be judged on the spread that is left after the condition difference has been explained. In DESeq2 this is the job of the GLM (generalized linear model), a model that writes the mean of the counts on the log scale as a sum of factors such as condition and batch. It is covered in detail in [note 04](04_glm_condition_batch.md).

## 4. What is the variance when the expression level itself wobbles?

Now I put the second cause, biological difference, into the model. It is easiest to think of it in two stages. First, each sample $j$ has an actual expression level of this gene, $\Lambda_j$ (capital Greek lambda), which differs a little around the mean $\mu_j$ depending on the donor or the cell state. Then sequencing counts reads with that $\Lambda_j$ as the expected value, and this stage is exactly the Poisson of section 2. As formulas:

$$K_j \mid \Lambda_j \sim \mathrm{Poisson}(\Lambda_j), \qquad E(\Lambda_j) = \mu_j, \qquad Var(\Lambda_j) = \alpha\mu_j^2$$

The vertical bar is read "given". So $K_j \mid \Lambda_j$ is the count when the expression level is fixed at $\Lambda_j$. $\Lambda_j$ is the actual expression level of sample $j$, the rate of the Poisson, and it is a random variable that takes a different value in each sample. $\alpha$ is the dispersion, which says how much the rate wobbles relative to its mean. Precisely, it is $Var(\Lambda)/E(\Lambda)^2$, the square of the coefficient of variation (CV, the standard deviation divided by the mean), so it has no units.

We also need to decide which distribution $\Lambda_j$ wobbles with, and the study text chooses the gamma distribution. The gamma is a continuous distribution that only takes values above 0. Its shape is set by two numbers, shape and scale; its mean is shape × scale and its variance shape × scale². The study text calls this model, with a gamma on the rate, the Poisson–gamma mixture. Setting the shape to $1/\alpha$ and the scale to $\alpha\mu$ gives exactly the mean and variance above.

$$\Lambda \sim \mathrm{Gamma}\left(\text{shape}=\tfrac{1}{\alpha},\ \text{scale}=\alpha\mu\right): \qquad E(\Lambda) = \tfrac{1}{\alpha}\cdot\alpha\mu = \mu, \qquad Var(\Lambda) = \tfrac{1}{\alpha}\cdot(\alpha\mu)^2 = \alpha\mu^2$$

Now let me find the mean and variance of the count. The tool is the law of total variance, which in words says "total variance = the mean of the variance with the rate fixed + the variance of the mean with the rate fixed". With the rate fixed, the count is Poisson, so $E(K\mid\Lambda)=\Lambda$ and $Var(K\mid\Lambda)=\Lambda$. In the formulas below, the $E$ and $Var$ outside the braces mean the mean and variance taken over $\Lambda$, which differs from sample to sample.

$$E(K) = E\{E(K\mid\Lambda)\} = E(\Lambda) = \mu \qquad (\text{law of total expectation})$$

$$Var(K) = \underbrace{E\{Var(K\mid\Lambda)\}}_{=\,E(\Lambda)\,=\,\mu} + \underbrace{Var\{E(K\mid\Lambda)\}}_{=\,Var(\Lambda)\,=\,\alpha\mu^2} = \mu + \alpha\mu^2$$

The mean stays $\mu$, and the variance now has two terms, each matching one of the causes in section 1. The first term, $\mu$, is the chance in counting, the Poisson noise that remains even with the rate fixed. The second term, $\alpha\mu^2$, is the share that comes from the expression level differing between samples. If α is 0, the second term disappears and we are back to the Poisson.

It is worth remembering here that the dispersion α is not the variance itself. The variance is $\mu + \alpha\mu^2$, so the same α gives a different variance when the mean changes. This relationship is laid out in a table in [note 02](02_negative_binomial.md).

Let me plug in gene A's Ctrl. $\mu = 106.667$, and for α I use gene A's final dispersion computed in [note 03](03_dispersion_estimation.md), $\alpha = 0.053147$. Note 03 attaches a prior to obtain this value; a prior is a distribution, set before seeing the data, for where α is likely to be (it comes up again at the end of section 5). This prior is not estimated from data but set by the study text for teaching.

```r
mu <- 106.667; alpha <- 0.053147
pois <- mu; bio <- alpha * mu^2
cat(sprintf("Poisson part %.1f + rate wobble part %.1f = variance %.1f (SD %.1f)\n",
            pois, bio, pois + bio, sqrt(pois + bio)))
cat(sprintf("SD of the rate = %.1f, %.1f%% of the mean\n", sqrt(bio), 100 * sqrt(alpha)))
```

```text
Poisson part 106.7 + rate wobble part 604.7 = variance 711.4 (SD 26.7)
SD of the rate = 24.6, 23.1% of the mean
```

Computed directly, only 106.7 of the variance of 711.4 is the Poisson share, and the remaining 604.7 comes from the wobble in the expression level. The second line of the output shows an easy way to read α. $\sqrt{\alpha} = 0.231$, so α = 0.053 means the actual expression level wobbles by about 23% of the mean (in standard deviation terms) from replicate to replicate.

![Count distributions under Poisson and Poisson-gamma](../../figures/01_same_mean_spread.png)

Both have a mean of 106.7, but when the expression level wobbles (orange), the counts spread much more widely. The dashed lines are the three Ctrl values of gene A; 130 sits in the tail of the Poisson (blue).

<details>
<summary>Code for the figure</summary>

```python
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(11)
mu, alpha, N = 106.667, 0.053147, 200_000          # gene A's Ctrl mean, final α
pois = rng.poisson(mu, N)
gp = rng.poisson(rng.gamma(1 / alpha, alpha * mu, N))   # draw the rate first, then Poisson with that rate

BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
bins = np.arange(0, 261, 4)
for k, color, label in [(pois, BLUE, f"Poisson (SD {pois.std():.1f})"),
                        (gp, ORANGE, f"Poisson-gamma, α = 0.053 (SD {gp.std():.1f})")]:
    ax.hist(k, bins=bins, weights=np.full(N, 1 / N), histtype="step", lw=2, color=color, label=label)
for x in (100, 130, 90):
    ax.axvline(x, color=MUTED, lw=1, ls=":")
ax.text(132, ax.get_ylim()[1] * 0.93, "Gene A, Ctrl: 90, 100, 130", color=MUTED, fontsize=9)
ax.set_xlabel("count in one replicate", color=INK)
ax.set_ylabel("fraction of replicates (bin width 4)", color=INK)
ax.set_title("Same mean (106.7), different spread", color=INK, loc="left", fontsize=11)
ax.legend(frameon=False, fontsize=9, loc="upper right", bbox_to_anchor=(1, 0.88))
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.set_xlim(0, 260)
fig.tight_layout()
fig.savefig("../figures/01_same_mean_spread.png", dpi=150, facecolor="white")
print(f"Poisson mean {pois.mean():.1f}, SD {pois.std():.1f} | Poisson-gamma mean {gp.mean():.1f}, SD {gp.std():.1f}")
```

```text
Poisson mean 106.7, SD 10.3 | Poisson-gamma mean 106.8, SD 26.7
```

Run from the `notes/` folder.

</details>

In fact, in the derivation above the choice of the gamma was used only to set the two values $E(\Lambda)$ and $Var(\Lambda)$. So $Var(K)=\mu+\alpha\mu^2$ holds whatever the shape of the distribution, as long as the rate wobbles with mean $\mu$ and variance $\alpha\mu^2$. What, then, is gained by choosing the gamma? The probability of a count can be written as a single formula without an integral. The resulting distribution is the negative binomial (NB), a count distribution that allows overdispersion; its definition and formula are covered in [note 02](02_negative_binomial.md).

The study text calls this mixture "one way of constructing the negative binomial". It does not mean that DESeq2 estimates and stores a hidden $\Lambda_j$ for each sample; section 6 opens a DESeq2 object to check this directly. Whether the integrated mixture really equals the NB probability is compared numerically under "Going deeper" below.

## 5. It is clear with many replicates, but can 3 tell?

First, the case with a very large sample. I counted a gene with a mean of 100 one hundred thousand times and compared three cases: a pure Poisson, a Poisson–gamma with α = 0.05, and a Poisson–gamma with α = 0.20. In the code and output, Poisson–gamma is shortened to GP.

```python
import numpy as np
rng = np.random.default_rng(1)

def gamma_poisson(mu, alpha, size):
    lam = rng.gamma(shape=1/alpha, scale=alpha*mu, size=size)  # E=mu, Var=alpha*mu^2
    return rng.poisson(lam)

n = 100_000
for label, k in [("Poisson",   rng.poisson(100, n)),
                 ("GP a=0.05", gamma_poisson(100, 0.05, n)),
                 ("GP a=0.20", gamma_poisson(100, 0.20, n))]:
    m, v = k.mean(), k.var(ddof=1)
    print(f"{label:10s} mean={m:7.2f} var={v:8.2f} var/mean={v/m:6.3f} (var-mean)/mean^2={(v-m)/m**2:7.4f}")
```

```text
Poisson    mean=  99.97 var=   99.61 var/mean= 0.996 (var-mean)/mean^2=-0.0000
GP a=0.05  mean=  99.97 var=  595.97 var/mean= 5.961 (var-mean)/mean^2= 0.0496
GP a=0.20  mean= 100.01 var= 2088.14 var/mean=20.879 (var-mean)/mean^2= 0.1988
```

| Model ($\mu=100$) | Theoretical variance | Observed $S^2$ | Theoretical var/mean $= 1+\alpha\mu$ | Observed $S^2/\bar K$ | Observed $\hat\alpha$ |
|---|---:|---:|---:|---:|---:|
| Poisson | 100 | 99.6 | 1 | 0.996 | −0.0000 |
| GP, $\alpha=0.05$ | 600 | 596.0 | 6 | 5.96 | 0.0496 |
| GP, $\alpha=0.20$ | 2,100 | 2,088.1 | 21 | 20.88 | 0.1988 |

All three have a mean of 100, but the variance is 6 and 21 times the Poisson one. The variance did not grow because the mean grew. This is exactly overdispersion.

The last column, $\hat\alpha$, is the variance formula solved backwards. Solving $Var(K)=\mu+\alpha\mu^2$ for α gives $\alpha = (Var(K)-\mu)/\mu^2$, and plugging in the sample mean $\bar K$ and the sample variance $S^2$ gives the simplest estimator.

$$\hat\alpha = \frac{S^2 - \bar K}{\bar K^2}$$

Solving for a parameter by matching the sample mean and variance to the theoretical formulas like this is called moment estimation. The sample variance is written $S^2$ rather than a lowercase $s^2$ because the lowercase $s_j$ is already the size factor (the symbols of appendix C of the study text). In the next code, `s2` and `m` are $S^2$ and $\bar K$. With 100,000 samples, $\hat\alpha$ recovers the true α with a relative error under 1% (0.0496 vs 0.05, 0.1988 vs 0.20).

But a real RNA-seq experiment has 3–6 replicates per condition. I repeated the same calculation 20,000 times each with 3, 6 and 30 replicates.

```python
import numpy as np
rng = np.random.default_rng(2)
reps, mu, a = 20_000, 100, 0.05
q = [2.5, 25, 50, 75, 97.5]
for n in (3, 6, 30):
    P  = rng.poisson(mu, (reps, n))
    GP = rng.poisson(rng.gamma(1/a, a*mu, (reps, n)))
    for label, K in (("Poisson", P), ("GP a=0.05", GP)):
        m, s2 = K.mean(1), K.var(1, ddof=1)
        print(f"n={n:2d} {label:9s} s2 quantiles={np.percentile(s2, q).round(0)} "
              f"s2/m={np.percentile(s2/m, q).round(2)} P(s2<m)={np.mean(s2 < m):.3f}")
```

```text
n= 3 Poisson   s2 quantiles=[  2.  30.  70. 140. 376.] s2/m=[0.02 0.3  0.71 1.4  3.77] P(s2<m)=0.626
n= 3 GP a=0.05 s2 quantiles=[  13.  165.  401.  812. 2305.] s2/m=[ 0.14  1.69  4.07  8.09 21.66] P(s2<m)=0.156
n= 6 Poisson   s2 quantiles=[ 17.  54.  87. 133. 254.] s2/m=[0.17 0.54 0.87 1.32 2.54] P(s2<m)=0.583
n= 6 GP a=0.05 s2 quantiles=[  96.  311.  509.  786. 1605.] s2/m=[ 1.01  3.18  5.15  7.86 15.32] P(s2<m)=0.025
n=30 Poisson   s2 quantiles=[ 55.  81.  98. 117. 157.] s2/m=[0.55 0.81 0.98 1.17 1.57] P(s2<m)=0.533
n=30 GP a=0.05 s2 quantiles=[325. 480. 581. 700. 981.] s2/m=[3.31 4.84 5.82 6.96 9.61] P(s2<m)=0.000
```

Summarizing where the sample variance $S^2$ scatters, in quantiles (20,000 repetitions, $\mu=100$):

| $n$ | Model (theoretical variance) | 2.5% | 25% | 50% | 75% | 97.5% | 95% range of $S^2/\bar K$ | $P(S^2 < \bar K)$ |
|---:|---|---:|---:|---:|---:|---:|---|---:|
| 3 | Poisson (100) | 2 | 30 | 70 | 140 | 376 | 0.02 – 3.77 | 0.626 |
| 3 | GP $\alpha=0.05$ (600) | 13 | 165 | 401 | 812 | 2,305 | 0.14 – 21.7 | 0.156 |
| 6 | Poisson (100) | 17 | 54 | 87 | 133 | 254 | 0.17 – 2.54 | 0.583 |
| 6 | GP $\alpha=0.05$ (600) | 96 | 311 | 509 | 786 | 1,605 | 1.01 – 15.3 | 0.025 |
| 30 | Poisson (100) | 55 | 81 | 98 | 117 | 157 | 0.55 – 1.57 | 0.533 |
| 30 | GP $\alpha=0.05$ (600) | 325 | 480 | 581 | 700 | 981 | 3.31 – 9.61 | 0.000 |

With 3 replicates, the $S^2$ of a Poisson gene spreads from 2 to 376, even though the true variance is 100. Wobble of this size is normal when a variance is measured from 3 values, because the sample variance wobbles roughly like a $\chi^2$ distribution with $n-1$ degrees of freedom (here 2). This relationship is exact only for the normal distribution, but at $\mu=100$ it is a good approximation.

So even a gene with α = 0.05, whose variance is 6 times the Poisson one, shows "variance smaller than the mean" 15.6% of the time when only 3 replicates are seen. Conversely, a Poisson gene with no overdispersion at all shows a variance larger than the mean 37.4% of the time (= 1 − 0.626). Only with 6 replicates does the lower end of the GP range of $S^2/\bar K$ reach 1 ($P(S^2<\bar K)=0.025$). Among the $n$ in the table, the two 95% ranges stop overlapping only at 30, and running $n$ more finely shows that they separate between 10 and 12 (see "Going deeper").

One more thing stands out: even for the Poisson gene, $P(S^2 < \bar K)$ is above 0.5 (0.626, 0.583, 0.533). The distribution of $S^2$ has a long right tail, so its median is below the mean (100). Using the moment estimator above on such genes gives a negative $\hat\alpha$; how DESeq2 handles these negative values is shown in section 6.

Plotting thousands of genes at once instead of one makes the situation clear at a glance.

![Sample variance vs sample mean per gene: 3 replicates and 30 replicates](../../figures/01_variance_vs_mean.png)

Each point is one gene. Blue points are Poisson genes and orange points are Poisson–gamma genes with α = 0.05 ("NB genes" in the figure); the two lines are their theoretical variances. With 3 replicates (left), the points scatter widely around their own line and the two clouds overlap a lot, but with 30 (right) they separate along the two lines. At small means the ratio of the two lines, $1+\alpha\mu$, is close to 1, so the lines themselves are close together.

<details>
<summary>Code for the figure</summary>

```python
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(12)
G, alpha = 1500, 0.05
mu = np.exp(rng.uniform(np.log(10), np.log(10_000), G))   # a different mean for each gene
BLUE, ORANGE, INK = "#2a78d6", "#eb6834", "#0b0b0b"
grid = np.geomspace(10, 10_000, 100)

fig, axes = plt.subplots(1, 2, figsize=(7, 4.5), dpi=150, sharex=True, sharey=True)
for ax, n in zip(axes, (3, 30)):
    P = rng.poisson(mu[:, None], (G, n))
    NB = rng.poisson(rng.gamma(1 / alpha, alpha * mu[:, None], (G, n)))
    for K, color, label in [(P, BLUE, "Poisson genes"), (NB, ORANGE, "NB genes, α = 0.05")]:
        m, s2 = K.mean(1), K.var(1, ddof=1)
        keep = s2 > 0                      # log axis: drop sample variance 0
        ax.scatter(m[keep], s2[keep], s=6, alpha=0.35, color=color, lw=0, label=label)
        print(f"n={n:2d} {label:18s} dropped (s2 = 0): {np.sum(~keep)}")
    ax.plot(grid, grid, color=BLUE, lw=1.8, label="variance = mean")
    ax.plot(grid, grid + alpha * grid**2, color=ORANGE, lw=1.8, label="variance = mean + 0.05·mean²")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title(f"{n} replicates per gene", color=INK, fontsize=11, loc="left")
    ax.set_xlabel("sample mean", color=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0].set_ylabel("sample variance", color=INK)
axes[0].set_ylim(1e-1, 1e8)
axes[1].legend(frameon=False, fontsize=8, loc="upper left", markerscale=2.5)
fig.tight_layout()
fig.savefig("../figures/01_variance_vs_mean.png", dpi=150, facecolor="white")
```

```text
n= 3 Poisson genes      dropped (s2 = 0): 0
n= 3 NB genes, α = 0.05 dropped (s2 = 0): 1
n=30 Poisson genes      dropped (s2 = 0): 0
n=30 NB genes, α = 0.05 dropped (s2 = 0): 0
```

Run from the `notes/` folder. The gene means were drawn evenly on the log scale between 10 and 10,000, and genes with a sample variance of 0 were dropped because they cannot be drawn on a log axis (the 1 in the output).

</details>

In the end, overdispersion is real and clearly visible with many samples, but with the 3–6 replicates of a real experiment, the sample variance of one gene cannot tell reliably whether there is overdispersion or how large α is. So DESeq2 borrows information across genes. It estimates a prior (a distribution for where a parameter is likely to be before the data are seen) from all genes and uses it; this approach is called empirical Bayes and is covered in [note 03](03_dispersion_estimation.md).

## 6. How does DESeq2 use this model?

There are three places where this note meets DESeq2 directly.

First, the variance formula is the same. DESeq2 models the count of gene $i$ in sample $j$ as $K_{ij} \sim NB(\mu_{ij}, \alpha_i)$ and writes the variance as $\mu_{ij} + \alpha_i\mu_{ij}^2$. The expected count $\mu_{ij}$ differs between samples, but there is one dispersion $\alpha_i$ per gene, shared by all samples of that gene. R's negative binomial functions (`dnbinom`, `rnbinom`) take an argument called `size` instead of α; passing `size = 1/α` gives the same variance. DESeq2 also uses `rnbinom(mu = mu, size = 1/dispersion)` when it generates example data.

That does not mean DESeq2 estimates a $\Lambda_j$ for each sample. The mixture is only a device for deriving the variance formula and the distribution, so even if you open a DESeq2 object after the analysis, none of its per-sample tables contains a rate.

Finally, the moment estimate of section 5, $\hat\alpha = (S^2 - \bar K)/\bar K^2$, sits inside DESeq2 almost unchanged. But it is used only as the starting point (initial value) of the dispersion calculation, so even if it comes out negative, the final dispersion used in the test is never negative.

The last two points are easiest to see by running it. I generate 200 genes and 6 samples (3 in each of two conditions) with DESeq2's fake-data function and run `DESeq()`.

```r
suppressPackageStartupMessages(library(DESeq2))
set.seed(1)
dds <- DESeq(makeExampleDESeqDataSet(n = 200, m = 6), quiet = TRUE)  # 200 genes, 6 samples
assayNames(dds)                            # tables with a value per sample: no rate table
mom <- DESeq2:::momentsDispEstimate(dds)   # (sample var - mean) / mean^2, on normalized counts
cat("genes with a negative moment estimate:", sum(mom < 0), "\n")
print(data.frame(moment = mom, final = dispersions(dds), row.names = rownames(dds))[mom < 0, ], digits = 3)
```

```text
[1] "counts" "mu"     "H"      "cooks" 
genes with a negative moment estimate: 4 
         moment  final
gene45  -0.0269  0.632
gene67  -0.1768 10.000
gene121 -0.1830  0.840
gene141 -0.4543 10.000
```

The first line lists the names of the tables with a value per sample (assays). `counts` holds the original counts, and `mu` the expected counts fitted by the model. `H` is each sample's influence on the fit (the diagonal of the hat matrix), and `cooks` is an outlier measure computed from it (Cook's distance). There is no table holding the rate $\Lambda_j$. `mu` too is the expected count $\mu_{ij}$ computed by the GLM, not a realization of $\Lambda_j$.

Comparing `moment` and `final` in the table at the bottom of the output, 4 of the 200 genes have a small sample variance and a negative moment estimate. Yet the final dispersion actually used in the test (`final`) is 0.632–10. 10 is the upper limit (`maxDisp`) DESeq2 sets for this data. So DESeq2 does not treat even genes with "variance smaller than the mean" as Poisson; it pulls them up toward the trend. The trend is a curve showing where the dispersion roughly lies as a function of mean expression. Pulling estimates that carry little information toward the overall tendency like this is called shrinkage, and the process is followed step by step in [note 03](03_dispersion_estimation.md).

Source excerpts, the number of genes that hit the lower limit `1e-8`, and how the assays change when there are outliers are collected under "Going deeper" below.

## Summary

- If the only source of variation were the chance in counting reads, counts would follow a Poisson, with variance equal to the mean. For gene A's Ctrl that is a variance of 106.7 and a standard deviation of 10.3.
- "The variance grows with the mean" happens under the Poisson too. What to check is whether the ratio $Var(K_j)/\mu_j$ to the expected count, which reflects the condition and the size factor, exceeds 1.
- If the expression level (rate) wobbles between biological replicates, the variance becomes $\mu + \alpha\mu^2$. $\mu$ is the counting chance and $\alpha\mu^2$ the share of biological differences between samples. α is the CV² of the rate, not the variance itself. Gene A's α = 0.053147 means the expression level wobbles by about 23% between replicates, and the Ctrl variance is then 711.4.
- Wobbling the rate with a gamma makes the count distribution a negative binomial, but this is only a way of deriving it; DESeq2 does not estimate and store a rate per sample.
- With 3–6 replicates, the sample variance of one gene alone can hardly tell overdispersion apart reliably. So DESeq2 borrows a trend and a prior from all genes.

The next note, [02](02_negative_binomial.md), derives the formula of the negative binomial from this mixture and covers the table of variance formulas, $CV^2 = 1/\mu + \alpha$, size factors and $\mu_{ij}=s_j q_{ij}$. Figure 1 of the study text is also recomputed there (study text 2.3–2.4, chapter 3).

## Exercises

### Problem 2

> As the mean count rose from 100 to 200, the variance also rose. Can you conclude from this alone that an NB is needed instead of a Poisson?

<details>
<summary>Solution</summary>

In short, no. Under a Poisson too, the variance grows along with the mean. What has to be checked is whether the variance exceeds the Poisson prediction $\mu$ relative to the expected mean, taking the conditions and offsets into account. Answer 2 in appendix B of the study text reaches the same conclusion. An offset here is a value that is not estimated but added as it is to the mean formula on the log scale; the size factor enters in this place in the form $\log s_j$.

Let me work it out in numbers. I computed the theoretical variances of the Poisson and the GP (α = 0.05) when the mean goes from 100 to 200, the range in which the sample variance of a 6-replicate Poisson sample falls, and the probability that the order of the two groups' sample variances flips.

```python
import numpy as np
rng = np.random.default_rng(4)
for mu in (100, 200):
    for a in (0.0, 0.05):
        v = mu + a*mu**2
        print(f"mu={mu:3d} alpha={a:.2f}  Var={v:6.0f}  Var/mu={v/mu:5.1f}  (Var-mu)/mu^2={(v-mu)/mu**2:.3f}")
n, reps = 6, 20_000
s2 = {}
for mu in (100, 200):
    s2[mu] = rng.poisson(mu, (reps, n)).var(1, ddof=1)
    print(f"Poisson n={n} mu={mu}: median s2={np.median(s2[mu]):5.1f}  95% band=({np.percentile(s2[mu],2.5):.0f}, {np.percentile(s2[mu],97.5):.0f})")
print(f"P(s2 at mu=200 > s2 at mu=100) = {np.mean(s2[200] > s2[100]):.3f}")
from scipy import stats
print(f"normal approx: P(F(5,5) > 0.5) = {stats.f.sf(0.5, 5, 5):.3f}")
```

```text
mu=100 alpha=0.00  Var=   100  Var/mu=  1.0  (Var-mu)/mu^2=0.000
mu=100 alpha=0.05  Var=   600  Var/mu=  6.0  (Var-mu)/mu^2=0.050
mu=200 alpha=0.00  Var=   200  Var/mu=  1.0  (Var-mu)/mu^2=0.000
mu=200 alpha=0.05  Var=  2200  Var/mu= 11.0  (Var-mu)/mu^2=0.050
Poisson n=6 mu=100: median s2= 87.4  95% band=(17, 255)
Poisson n=6 mu=200: median s2=175.0  95% band=(34, 523)
P(s2 at mu=200 > s2 at mu=100) = 0.768
normal approx: P(F(5,5) > 0.5) = 0.767
```

| Mean $\mu$ | Poisson variance | Poisson var/mean | GP (α = 0.05) variance | GP var/mean | 95% range of $S^2$ for an $n=6$ Poisson sample |
|---:|---:|---:|---:|---:|---|
| 100 | 100 | 1 | 600 | 6 | 17 – 255 |
| 200 | 200 | 1 | 2,200 | 11 | 34 – 523 |

Under the Poisson, when the mean goes from 100 to 200 the variance also doubles exactly, from 100 to 200. The observation that "the variance rose" fits the Poisson perfectly.

What is needed to tell them apart is the variance/mean. Under the Poisson it is 1 whatever the mean, while for a GP with α = 0.05 the ratio itself grows with the mean, from 6 to 11 ($1+\alpha\mu$). On the other hand, $(Var-\mu)/\mu^2 = \alpha$ stays at 0.05. So within the NB model, α is the overdispersion parameter with the size of the mean removed; unlike the variance/mean, it does not change along with $\mu$. Across real genes, though, α does tend to vary with the mean (a trend). The `dispMeanRel = 4/x + 0.1` of DESeq2's example-data function is such a relationship, and it is covered in [note 03](03_dispersion_estimation.md).

To say "the variance rose" from real data, the wobble of the sample variance itself must also be taken into account. At $n=6$, the 95% range of $S^2$ for a real Poisson(100) is 17–255, and for Poisson(200) it is 34–523. Even when the true variance is exactly double, the sample variance on the $\mu=200$ side comes out larger only 76.8% of the time. Approximating the counts by a normal distribution, the ratio of two sample variances can be computed with an $F_{5,5}$ distribution, and $P(F>0.5)=0.767$ obtained that way is almost the same. The order flips about 23% of the time, so a single observation that "the sample variance of the mean-200 group was larger" cannot even establish that the variance increased.

So the right question is whether, once the expected mean $\mu_{ij}$ that reflects the conditions and offsets is in place, the variance left around it systematically exceeds $\mu_{ij}$. DESeq2 judges this not from the sample variance of one gene but from the mean model (GLM) and a dispersion trend borrowed from all genes (notes [03](03_dispersion_estimation.md) and [04](04_glm_condition_batch.md)).

</details>

## Going deeper

The R code was run with R 4.5.2 and DESeq2 1.50.2, and the Python simulations with Python 3.13, numpy 2.2.6 and scipy 1.17.1, all with fixed seeds.

<details>
<summary>Why does the Poisson variance equal the mean?</summary>

If $K \sim \mathrm{Poisson}(\mu)$, then $P(K=k) = e^{-\mu}\mu^k / k!$. For this distribution $E(K)=\mu$ and $E[K(K-1)] = \mu^2$, so expanding the variance shows it equals the mean.

$$Var(K) = E[K(K-1)] + E(K) - E(K)^2 = \mu^2 + \mu - \mu^2 = \mu$$

</details>

<details>
<summary>Does integrating the mixture really give an NB? (Python, R)</summary>

I compared the integral over $\lambda$ of the product of the gamma density and the Poisson probability with scipy's NB probability `nbinom(n=1/α, p=1/(1+αμ)).pmf`, with $\mu = 30$ and $\alpha = 0.2$. Figure 1 of the study text (mean fixed at 30; probabilities for Poisson / α = 0.05 / α = 0.2) belongs to section 2.3 of the study text, so it is recomputed in [note 02](02_negative_binomial.md); here I only check in numbers the statement in 2.2 that "the Poisson–gamma mixture is one way of constructing the NB".

```python
import numpy as np
from scipy import stats, integrate
mu, a = 30, 0.2
nb = stats.nbinom(1/a, 1/(1 + a*mu))   # size=1/a, prob=1/(1+a*mu) -> mean mu, var mu+a*mu^2
for k in (5, 20, 30, 60):
    f = lambda lam: stats.poisson.pmf(k, lam) * stats.gamma.pdf(lam, 1/a, scale=a*mu)
    print(f"k={k:2d}  mixture integral={integrate.quad(f, 0, np.inf)[0]:.6f}  nbinom.pmf={nb.pmf(k):.6f}")
```

```text
k= 5  mixture integral=0.003469  nbinom.pmf=0.003469
k=20  mixture integral=0.028970  nbinom.pmf=0.028970
k=30  mixture integral=0.027064  nbinom.pmf=0.027064
k=60  mixture integral=0.003637  nbinom.pmf=0.003637
```

All four points agree to the sixth decimal place. The formal derivation is in note 02.

I did the same comparison in R, comparing the variance of values counted with `rpois` after drawing the rate with `rgamma` (the Poisson–gamma mixture) against R's `rnbinom(mu, size = 1/α)`.

```r
set.seed(2026); N <- 2e6
for (a in c(0.05, 0.2)) {
  lam  <- rgamma(N, shape = 1/a, scale = a*30)   # E=30, Var=a*30^2
  k_gp <- rpois(N, lam)                          # Poisson-gamma mixture
  k_nb <- rnbinom(N, mu = 30, size = 1/a)        # R's NB, mu/size parameterization
  cat(sprintf("alpha=%.2f  Var(Lambda)=%.1f (theory %.0f) | gamma-pois var=%.1f | rnbinom var=%.1f | theory %.0f\n",
              a, var(lam), a*900, var(k_gp), var(k_nb), 30 + a*900))
}
```

```text
alpha=0.05  Var(Lambda)=44.9 (theory 45) | gamma-pois var=74.9 | rnbinom var=75.0 | theory 75
alpha=0.20  Var(Lambda)=179.8 (theory 180) | gamma-pois var=209.9 | rnbinom var=210.7 | theory 210
```

$Var(\Lambda) = \alpha\mu^2$ (45, 180) and the sum of the two terms $\mu + \alpha\mu^2$ (75, 210) each match the derived formulas, and R's `rnbinom(size = 1/α)` gives the same variance. Incidentally, the GP samples in this note were made in two stages, exactly as in the formulas of section 4: drawing the rate from `rng.gamma(shape=1/α, scale=αμ)` and counting with `rng.poisson(λ)`. The NB probability function was used only in the integral comparison above.

</details>

<details>
<summary>What if the replicates go up to 8–20?</summary>

In the table of section 5 the two 95% ranges separated only at $n=30$. I ran the values between 6 and 30 more finely and compared the 97.5% quantile of $S^2/\bar K$ on the Poisson side with the 2.5% quantile on the GP (α = 0.05) side.

```python
import numpy as np
rng = np.random.default_rng(3)
reps, mu, a = 20_000, 100, 0.05
for n in (8, 10, 12, 15, 20):
    P  = rng.poisson(mu, (reps, n))
    GP = rng.poisson(rng.gamma(1/a, a*mu, (reps, n)))
    pq = np.percentile(P.var(1, ddof=1) / P.mean(1), 97.5)
    gq = np.percentile(GP.var(1, ddof=1) / GP.mean(1), 2.5)
    print(f"n={n:2d}  Poisson 97.5%={pq:.2f}  GP 2.5%={gq:.2f}  overlap={pq > gq}")
```

```text
n= 8  Poisson 97.5%=2.30  GP 2.5%=1.45  overlap=True
n=10  Poisson 97.5%=2.09  GP 2.5%=1.74  overlap=True
n=12  Poisson 97.5%=1.97  GP 2.5%=2.09  overlap=False
n=15  Poisson 97.5%=1.88  GP 2.5%=2.38  overlap=False
n=20  Poisson 97.5%=1.73  GP 2.5%=2.79  overlap=False
```

Even for a gene like $\alpha=0.05$, $\mu=100$, whose variance is 6 times the Poisson one, the replicates have to exceed 10 before the two 95% ranges separate. The difference at $n=12$ (1.97 vs 2.09) is borderline.

</details>

<details>
<summary>A look at the DESeq2 source: implementation checks (1)–(4)</summary>

Within the scope of this note, DESeq2 is touched directly in four places: (1) which parameterization R's NB uses and how it relates to the gamma mixture, (2) how DESeq2's simulator generates counts, (3) the study text's claim that DESeq2 does not store per-sample gamma rates, and (4) where in the code the moment estimate that inverts $Var=\mu+\alpha\mu^2$ lives.

**(1) The mu/size parameterization of `stats::dnbinom` and the gamma mixture**

```r
txt <- gsub("_\b", "", capture.output(tools::Rd2txt(utils:::.getHelpFile(help("dnbinom")))))
cat(grep("alternative parametrization|variance is|gamma mixture", txt, ignore.case = TRUE, value = TRUE), sep = "\n")
```

```text
      mu: alternative parametrization via mean: see ‘Details’.
     An alternative parametrization (often used in ecology) is by the
     where ‘prob’ = ‘size/(size+mu)’.  The variance is ‘mu + mu^2/size’
     ‘rnbinom’ uses the derivation as a gamma mixture of Poisson
     ## Alternative parametrization
```

The variance given in the R manual is `mu + mu^2/size`, so setting `size = 1/α` gives exactly $\mu + \alpha\mu^2$. `rnbinom` itself also generates random numbers as a "gamma mixture of Poisson". The [9] that the study text cites is this very R manual, and the mixture construction and parameterization taken from [9] match the manual text. But section 2.2 of the study text attaches [9] right after the sentence "this does not mean that DESeq2 estimates and stores a hidden gamma rate for each real sample". [9], an R manual, says nothing about how DESeq2 stores things, so the evidence for that part is what I checked directly in (3) below.

**(2) DESeq2's simulator `makeExampleDESeqDataSet` uses `rnbinom(mu, size = 1/dispersion)`**

```r
suppressPackageStartupMessages(library(DESeq2))
print(args(makeExampleDESeqDataSet))
b <- deparse(body(makeExampleDESeqDataSet))
cat(grep("dispersion <-|rnbinom", b, value = TRUE), sep = "\n")
```

```text
function (n = 1000, m = 12, betaSD = 0, interceptMean = 4, interceptSD = 2, 
    dispMeanRel = function(x) 4/x + 0.1, sizeFactors = rep(1, 
        m)) 
NULL
    dispersion <- dispMeanRel(2^(beta[, 1]))
    countData <- matrix(rnbinom(m * n, mu = mu, size = 1/dispersion), 
```

DESeq2 too passes `size = 1/dispersion` when it creates example data. So DESeq2's α is R's `1/size`. The default `dispMeanRel = 4/x + 0.1` is a relationship in which α grows as the mean gets smaller; what this trend means is covered in note 03.

**(3) `DESeq()` does not store per-sample gamma rates**

```r
suppressPackageStartupMessages(library(DESeq2))
set.seed(1)
dds <- DESeq(makeExampleDESeqDataSet(n = 200, m = 6), quiet = TRUE)
cat("assayNames:", assayNames(dds), "\n")
cat("mcols:", names(mcols(dds)), "\n")
cat("dim(assays(dds)$mu):", dim(assays(dds)$mu), "\n")
```

```text
assayNames: counts mu H cooks 
mcols: trueIntercept trueBeta trueDisp baseMean baseVar allZero dispGeneEst dispGeneIter dispFit dispersion dispIter dispOutlier dispMAP Intercept condition_B_vs_A SE_Intercept SE_condition_B_vs_A WaldStatistic_Intercept WaldStatistic_condition_B_vs_A WaldPvalue_Intercept WaldPvalue_condition_B_vs_A betaConv betaIter deviance maxCooks 
dim(assays(dds)$mu): 200 6 
```

In this run (m = 6, no outlier replacement), the assays with a sample dimension (columns) are the four `counts`, `mu`, `H` and `cooks`, and the gene-level `mcols` has only dispersion-related columns with one value per gene. `assays(dds)$mu` is a genes×samples matrix, but it is the expected count $\mu_{ij} = s_j q_{ij}$ fitted by the GLM, not a realization of $\Lambda_j$ ($q_{ij}$ is the expected expression without the size factor; note 02). There is no storage named `Lambda`, `gamma`, `rate` or the like either. The statement in section 2.2 of the study text, "this does not mean that DESeq2 estimates and stores a hidden gamma rate for each real sample", is in fact right. Note that `trueIntercept`, `trueBeta` and `trueDisp` are the answers planted by the simulator and do not exist in real data.

The assays are not always these four, though. If any design cell has 7 or more replicates (`minReplicatesForReplace = 7`, `any(nOrMoreInCell(...))`), `DESeq()` runs an outlier replacement step, and only samples in such cells can be replaced. So the assays come in three states.

| Situation | `assayNames` |
|---|---|
| No cell with 7 or more (m = 6 above) | `counts mu H cooks` |
| A cell with 7 or more, but no count replaced | `counts mu H cooks originalCounts` (+ a `replace` column in `mcols`) |
| A Cook's outlier in a sample of a cell with 7 or more, which is replaced | `counts mu H cooks replaceCounts replaceCooks` |

The third case also arises in an unbalanced design like 7 vs 3 if a sample in the cell with 7 has an outlier. In every case, the added assays are the original counts, the replaced counts and the Cook's distances at that point, not per-sample rates.

```r
suppressPackageStartupMessages(library(DESeq2))
cat("minReplicatesForReplace:", formals(DESeq)$minReplicatesForReplace, "\n")
cat(grep("sufficientReps <-", deparse(body(DESeq)), value = TRUE), "\n")
set.seed(1)
dds14 <- makeExampleDESeqDataSet(n = 200, m = 14)                   # 7 each in condition A/B
d0 <- DESeq(dds14, quiet = TRUE)
cat("m=14, no outlier:", assayNames(d0), "| sum(mcols$replace):", sum(mcols(d0)$replace), "\n")
counts(dds14)[1, 1] <- 50000L                                       # plant one outlier
cat("m=14, with outlier:", assayNames(DESeq(dds14, quiet = TRUE)), "\n")
dds10 <- makeExampleDESeqDataSet(n = 200, m = 10)
dds10$condition <- factor(rep(c("A", "B"), c(7, 3))); counts(dds10)[1, 1] <- 50000L   # 7 vs 3, outlier on the A side
cat("m=10 (7 vs 3), with outlier:", assayNames(DESeq(dds10, quiet = TRUE)), "\n")
```

```text
minReplicatesForReplace: 7 
    sufficientReps <- any(nOrMoreInCell(attr(object, "modelMatrix"),  
m=14, no outlier: counts mu H cooks originalCounts | sum(mcols$replace): 0 
m=14, with outlier: counts mu H cooks replaceCounts replaceCooks 
m=10 (7 vs 3), with outlier: counts mu H cooks replaceCounts replaceCooks 
```

**(4) $\hat\alpha = (S^2 - \bar K)/\bar K^2$ sits in `momentsDispEstimate` as it is**

```r
suppressPackageStartupMessages(library(DESeq2))
print(DESeq2:::momentsDispEstimate)
print(formals(estimateDispersionsGeneEst)$minDisp)
```

```text
function (object) 
{
    xim <- if (!is.null(normalizationFactors(object))) {
        mean(1/MatrixGenerics::colMeans(normalizationFactors(object)))
    }
    else {
        mean(1/sizeFactors(object))
    }
    bv <- mcols(object)$baseVar
    bm <- mcols(object)$baseMean
    (bv - xim * bm)/bm^2
}
<bytecode: 0x643f5c8db5a0>
<environment: namespace:DESeq2>
[1] 1e-08
```

`(bv - xim * bm)/bm^2` is $(S^2 - \overline{1/s}\,\bar K)/\bar K^2$. Here `bv` ($S^2$) and `bm` ($\bar K$) are the row-wise variance and mean of the normalized counts $K_{ij}/s_j$ (counts divided by the size factor), and `xim` is $\overline{1/s} = \mathrm{mean}_j(1/s_j)$. In the dds of section 6 the size factors are 0.944–1.122 and `xim` is 0.950, so it is almost the same as the formula of section 5.

Why subtract $\overline{1/s}\,\bar K$ instead of $\bar K$? If $K_{ij}\sim NB(\mu_{ij},\alpha)$ and $\mu_{ij}=s_j q_{ij}$, then $Var(K_{ij}/s_j) = (\mu_{ij}+\alpha\mu_{ij}^2)/s_j^2 = q_{ij}/s_j + \alpha q_{ij}^2$. On the normalized scale the Poisson term becomes $\mu_{ij}/s_j^2 = q_{ij}/s_j$, which is smaller than $q_{ij}$ if $s_j>1$ and larger if $s_j<1$. If every sample has the same $q_{ij}$ (intercept-only), $E(S^2) = \overline{1/s}\,q + \alpha q^2$, which is why $\overline{1/s}\,\bar K$ is subtracted. If all size factors are 1, `xim = 1` and it is the same as the formula of section 5.

If there is a condition effect, then as section 3 showed, the difference between the group means leaks into $S^2$ and inflates this value. So DESeq2 compares the moment estimate with `roughDispEstimate`, which reflects the design, takes the smaller one (`pmin`), clips it to `[minDisp, maxDisp]` and uses it only as the initial value `alpha_init` of the gene-wise optimization. The default of `minDisp` is `1e-08`, and `maxDisp` is `max(10, ncol(object))`.

```r
suppressPackageStartupMessages(library(DESeq2))
print(DESeq2:::roughDispEstimate)
b <- deparse(body(estimateDispersionsGeneEst))
i <- grep("momentsDisp|alpha_hat <- pmin|maxDisp <-|alpha_init <-|noIncrease", b)
cat(paste(i, b[i]), sep = "\n")
```

```text
function (y, x) 
{
    mu <- linearModelMu(y, x)
    mu <- matrix(pmax(1, mu), ncol = ncol(mu))
    m <- nrow(x)
    p <- ncol(x)
    est <- MatrixGenerics::rowSums(((y - mu)^2 - mu)/mu^2)/(m - 
        p)
    pmax(est, 0)
}
<bytecode: 0x58cc0ce1d5d8>
<environment: namespace:DESeq2>
37         momentsDisp <- momentsDispEstimate(objectNZ)
38         alpha_hat <- pmin(roughDisp, momentsDisp)
49     maxDisp <- max(10, ncol(object))
50     alpha_hat <- alpha_hat_new <- alpha_init <- pmin(pmax(minDisp, 
121         noIncrease <- last_lp < initial_lp + abs(initial_lp)/1e+06
122         dispGeneEst[which(noIncrease)] <- alpha_init[which(noIncrease)]
```

`roughDispEstimate` ends with `pmax(est, 0)`, so it never gives a negative value. The two `noIncrease` lines (121–122) reset the gene-wise value to the initial value for genes whose objective did not rise at `niter = 1`, so an initial value that came from the moment side can remain as the final gene-wise value of some genes.

After the initial value comes gene-wise estimation (the Cox-Reid adjusted likelihood) → trend → MAP. The likelihood, with the observations held fixed, measures how plausible a candidate parameter makes those observations, and the Cox-Reid adjustment is a correction term that reduces the underestimation of the dispersion caused by estimating the mean from the same data. The MAP (maximum a posteriori) is the most plausible value when the likelihood and the prior are considered together. This path and its details (the `noIncrease` reset, the grid refit for genes that did not converge) are looked at with the source in implementation check (1) of [note 03](03_dispersion_estimation.md).

What this note wants to know is how far it carries into the DESeq2 results when the sample variance is smaller than the sample mean and the moment estimate is negative. I counted on the same dds as in (3).

```r
suppressPackageStartupMessages(library(DESeq2))
set.seed(1)
dds <- DESeq(makeExampleDESeqDataSet(n = 200, m = 6), quiet = TRUE)   # same dds as in (3) (no allZero genes)
rough <- DESeq2:::roughDispEstimate(y = counts(dds, normalized = TRUE), x = model.matrix(design(dds), colData(dds)))
mom   <- DESeq2:::momentsDispEstimate(dds)
init  <- pmin(pmax(1e-8, pmin(rough, mom)), max(10, ncol(dds)))      # alpha_init
gw    <- mcols(dds)$dispGeneEst
cat("moment < 0:", sum(mom < 0), "| alpha_init = 1e-8:", sum(init <= 1e-8),
    "(rough = 0:", sum(init <= 1e-8 & rough == 0), ", moment < 0:", sum(init <= 1e-8 & mom < 0), ")\n")
cat("dispGeneEst = 1e-8:", sum(gw <= 1e-8), "| dispGeneEst == moment:", sum(abs(gw - mom) < 1e-12),
    "| dispersions(dds) = 1e-8:", sum(dispersions(dds) <= 1e-8), "\n")
print(data.frame(moments = mom, rough, dispGeneEst = gw, dispFit = mcols(dds)$dispFit,
                 dispersion = dispersions(dds), row.names = rownames(dds))[mom < 0, ], digits = 4)
```

```text
moment < 0: 4 | alpha_init = 1e-8: 17 (rough = 0: 17 , moment < 0: 4 )
dispGeneEst = 1e-8: 18 | dispGeneEst == moment: 1 | dispersions(dds) = 1e-8: 0 
        moments rough dispGeneEst dispFit dispersion
gene45  -0.0269     0       1e-08  0.9373     0.6320
gene67  -0.1768     0       1e-08 33.9399    10.0000
gene121 -0.1830     0       1e-08  1.3373     0.8405
gene141 -0.4543     0       1e-08 15.8081    10.0000
```

Of the 200 genes, 4 have a negative moment estimate, and 17 have an initial value at the lower limit `minDisp = 1e-08`. In all 17 cases `roughDispEstimate` gave 0, and the 4 negative moments also have rough = 0. In this example the lower limit is hit mainly through rough = 0, and there is no case where a negative moment alone caused it.

The gene-wise estimate `dispGeneEst` is 1e-08 for 18 genes. At the gene-wise stage DESeq2 sets "the Poisson limit α → 0" not to exactly 0 but to $10^{-8}$. For one gene (`dispGeneEst == moment`) the initial value came from the moment side, and the `noIncrease` reset leaves its gene-wise value equal to the moment estimate itself.

Still, none of the dispersions used in the test (`dispersions(dds)`, the result of trend + MAP) is at the lower limit. The final values of the 4 negative-moment genes are 0.63–10, and 10 is the `maxDisp = max(10, ncol)` upper limit. So DESeq2 does not treat even genes with "variance smaller than the mean" as Poisson; it pulls them up toward the trend of all genes (note 03).

</details>

<details>
<summary>Common misconceptions</summary>

- It is easy to think "the variance grows with the mean, so it isn't Poisson", but under the Poisson too $Var=\mu$, so the variance grows with the mean. The criterion is not an increase in variance but $Var(K_{ij})/\mu_{ij} > 1$ relative to the expected count (section 3, Problem 2).
- Nor is the NB a distribution made for RNA-seq. Both the Poisson and the NB existed before; the Poisson is an approximation of sampling, and the NB is the marginal distribution of a Poisson whose rate wobbles as a gamma (section 2, the `?dnbinom` quotation in implementation check (1)).
- Nor does DESeq2 estimate a gamma rate $\Lambda_j$ for each sample. The mixture is a derivation device, and no per-sample rate is stored anywhere in `dds` (section 6, implementation check (3)).
- A sample variance smaller than the sample mean does not mean the gene has no overdispersion either. With 3 replicates, even a gene with true α = 0.05 shows $S^2<\bar K$ 15.6% of the time, so one sample variance cannot tell (section 5).
- α is not the variance. It is the $CV^2$ of the rate, so it has no units, while the variance, $\mu + \alpha\mu^2$, depends on the mean. A detailed table is in note 02 (section 4).
- That does not make the Poisson simply a wrong model. For technical replicates, or when $\alpha\mu \ll 1$ (e.g. $\mu=10$ and $\alpha=0.01$ give $Var/\mu=1.1$), it is a good approximation. The problem is that for biological replicates the $\alpha\mu^2$ term becomes too large to ignore (the split into two terms in section 4).

</details>

<details>
<summary>Where this note differs from the study text, or adds to it</summary>

Explanations in the study text that are not in the table below matched this note's derivations, simulations and the actual behaviour of DESeq2. These include the Poisson's $E=Var=\mu$, $Var(K_j)=\mu_j+\alpha\mu_j^2$ from the law of total variance, the fact that DESeq2 does not store per-sample gamma rates (implementation check (3)), and answer 2 of appendix B, that "the variance grows with the mean" alone does not call for an NB.

| # | Study text claim (2.1–2.2, appendix B) | Result | Evidence |
|---|---|---|---|
| 1 | The Poisson is "not a newly invented distribution but an approximate application of an existing one" | Conceptual statement (not something to verify) | Nothing to check numerically. The Poisson limit of the binomial in section 2 is a supplementary explanation not in the study text, not a verification of this statement |
| 2 | $K_j\mid\Lambda_j \sim \mathrm{Poisson}(\Lambda_j)$, $E(\Lambda_j)=\mu_j$, $Var(\Lambda_j)=\alpha\mu_j^2$ | Model assumption (not something to verify) | The study text states this as an assumption. All I checked is that the shape $1/\alpha$ and scale $\alpha\mu$ chosen in this note satisfy it (gamma parameterization in section 4; $Var(\Lambda)$ = 44.9 / 179.8 with `rgamma`, theory 45 / 180) |
| 3 | "The Poisson–gamma mixture is one way of constructing the NB" [9] | Matches, but the placement of [9] is inaccurate | `?dnbinom`: "‘rnbinom’ uses the derivation as a gamma mixture of Poisson"; the integrated mixture = `nbinom.pmf` (four points). The study text attaches [9] to the end of the following sentence, "this does not mean that DESeq2 ... stores", but [9], an R manual, supports only the mixture and parameterization part. The evidence for how DESeq2 stores things is the direct check in implementation check (3) |
| 4 | (Implicit) the gamma parameterization in section 2.2 of the study text assumes shape $1/\alpha$, scale $\alpha\mu$ | Matches, but the study text never states the shape/scale explicitly | This is the only gamma parameterization satisfying $E(\Lambda)=\mu$, $Var(\Lambda)=\alpha\mu^2$, so there is no logical gap. This note states it explicitly |
| 5 | The study text does not write out the relationship between DESeq2's α and R's `size` in 2.2 | Not stated (not a discrepancy) | Confirmed $\alpha = 1/\text{size}$ from `rnbinom(mu = mu, size = 1/dispersion)` in the body of `makeExampleDESeqDataSet`. Covered formally in note 02 |

There were no formulas or numbers to record as discrepancies, and the one difference from the study text is the citation placement in #3. Figure 1 of the study text (p.7) belongs to 2.3, so its numbers are covered in note 02.

</details>

---

← Previous: [00. The overall map](00_overview.md) · Next: [02. The negative binomial and size factors](02_negative_binomial.md) →

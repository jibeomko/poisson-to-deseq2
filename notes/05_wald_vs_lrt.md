# 05. 같은 fold change인데 왜 p 값이 다를까? Wald 검정과 LRT

[04 노트](04_glm_condition_batch.md)에서 유전자 A의 LFC(0.977)와 SE(0.289)를 구했는데, DESeq2는 여기에 p = 0.00072를 붙여. 그런데 LFC가 똑같이 1인 유전자끼리도 p 값은 0.00000057부터 0.32까지 갈려. 이 노트에서는 그 p가 어떻게 나오는지 Wald 검정을 손으로 따라가 보고, 세 그룹 중 "어딘가" 다른지를 한 번에 묻는 LRT와 비교해 볼게.

> 교재 9장, 10장, 12장 · 코드: [05_wald_vs_lrt.R](../05_wald_vs_lrt.R)

코드는 접이식 블록 안의 것까지 모두 `options(width = 110)`을 둔 R 세션 하나에서 위에서부터 차례로 실행했어(기본 폭 80에서는 값은 같고 출력의 줄바꿈만 달라져).

## 1. 유전자 A의 p 값은 어디서 나올까?

유전자 A의 count는 Ctrl이 100, 130, 90이고 Starvation이 200, 250, 180이야. count는 한 sample에서 그 유전자에 배정된 read 수이고, sample마다 다른 sequencing 깊이를 맞춰 주는 배율인 size factor는 모두 1로 둘게. [04 노트](04_glm_condition_batch.md)에서는 이 유전자를 GLM(일반화 선형모형)으로 적합했어. GLM은 count의 평균을 log scale에서 조건 효과의 합으로 쓰는 모형인데, 결과로 두 숫자가 나왔지.

하나는 LFC = 0.977이야. log2 fold change(LFC)는 두 조건 평균의 비를 log2로 나타낸 값인데, 그룹 평균이 106.667과 210이니 비는 1.97배이고 log2로는 0.977이 돼. LFC가 1이면 2배, −1이면 절반이라는 뜻이야. 다른 하나는 SE = 0.289야. SE(표준오차)는 같은 실험을 다시 하면 추정치가 얼마나 흔들릴지를 나타내.

그러면 궁금한 건 하나야. Starvation이 실제로는 아무 효과가 없는데도, replicate 세 개씩의 우연한 차이만으로 LFC 0.977이 나올 수 있을까?

### 추정치를 그 흔들림으로 나누기

진짜 효과가 0이라고 해 볼게. 그래도 실험을 반복할 때마다 LFC 추정치는 0 근처에서 SE 정도의 폭으로 흩어져. 그러니 관측한 0.977이 0에서 SE 몇 개만큼 떨어져 있는지 세어 보면 돼. 0.977 ÷ 0.289 = 3.38, 즉 SE 3.38개만큼 떨어져 있어. 정규분포에서 평균으로부터 표준편차 3.38개 이상 벗어날 확률은 양쪽을 합쳐 약 0.07%인데, 이게 바로 p 값이야.

이 계산을 Wald 검정이라고 불러. 식으로 쓰면 이래.

$$Z=\frac{\hat\beta}{SE(\hat\beta)},\qquad p=2\,\Phi(-\lvert Z\rvert)$$

- $\hat\beta$: 추정한 LFC(log2 단위 계수). 모자(^)는 데이터에서 추정한 값이라는 표시야.
- $SE(\hat\beta)$: 그 추정치의 표준오차.
- $Z$: 추정치가 0에서 SE 몇 개만큼 떨어져 있는지. `results()` 표의 `stat` 열이야.
- $\Phi$: 표준정규분포의 누적분포함수. $\Phi(-3.38)$은 표준정규분포에서 −3.38보다 왼쪽의 넓이야.
- 앞의 2: 증가와 감소를 모두 극단적인 결과로 치기 때문에 붙어(양측 검정).

유전자 A에 대입하면 $Z=0.977280/0.289057=3.381$, $p=2\,\Phi(-3.381)=0.00072$야.

이번엔 DESeq2로 같은 계산을 해 볼게. 값이 하나 더 필요한데, 바로 dispersion α야. α는 같은 조건의 replicate가 Poisson이 예상하는 것보다 얼마나 더 퍼지는지를 나타내. 여기서는 [03 노트](03_dispersion_estimation.md)에서 구한 최종 α = 0.053147을 고정해 넣어. `nbinomWaldTest()`는 `DESeq()`가 차례로 하는 여러 단계 중 마지막 검정 단계만 따로 부르는 함수야.

```r
suppressPackageStartupMessages(library(DESeq2))
k  <- c(100L, 130L, 90L, 200L, 250L, 180L)                    # 유전자 A: Ctrl 3개, Starvation 3개
cd <- data.frame(condition = factor(rep(c("Ctrl", "Starvation"), each = 3)),
                 row.names = paste0(rep(c("Ctrl", "Starvation"), each = 3), "_", 1:3))
ddsA <- DESeqDataSetFromMatrix(matrix(k, nrow = 1, dimnames = list("geneA", rownames(cd))), cd, ~ condition)
sizeFactors(ddsA) <- rep(1, 6)                                # size factor 모두 1
dispersions(ddsA) <- 0.053147                                 # 03 노트의 최종 alpha
ddsA <- nbinomWaldTest(ddsA, quiet = TRUE)
resA <- results(ddsA)
print(as.data.frame(resA)[, c("log2FoldChange", "lfcSE", "stat", "pvalue")], digits = 6)
z <- resA$log2FoldChange / resA$lfcSE
cat("손 계산: Z = LFC/SE =", round(z, 4), "  p = 2*pnorm(-|Z|) =", signif(2 * pnorm(-abs(z)), 5), "\n")
cat("95% CI = LFC ± 1.96*SE:", round(resA$log2FoldChange + c(-1, 1) * qnorm(0.975) * resA$lfcSE, 4), "\n")
```

```
      log2FoldChange    lfcSE    stat      pvalue
geneA        0.97728 0.289057 3.38093 0.000722408
손 계산: Z = LFC/SE = 3.3809   p = 2*pnorm(-|Z|) = 0.00072241
95% CI = LFC ± 1.96*SE: 0.4107 1.5438
```

직접 계산해 보니 `stat`은 LFC를 lfcSE로 나눈 값 그대로이고, `pvalue`도 손으로 계산한 2Φ(−|Z|)와 똑같이 나왔어.

### 정규분포를 가정하는 대상은 count일까?

나도 처음엔 Wald 검정이 정규분포를 쓰니까 count도 정규분포여야 한다고 생각했어. 하지만 그렇지 않아. count는 음이항분포(NB)로 모형화해. NB는 과산포, 즉 같은 조건 replicate 사이의 퍼짐이 Poisson이 예상하는 것보다 큰 현상을 허용하는 count 분포야([02](02_negative_binomial.md)).

정규분포를 가정하는 대상은 추정치 $\hat\beta$야. 같은 실험을 여러 번 반복했다고 상상하면 $\hat\beta$도 매번 조금씩 다르게 나올 텐데, 그 $\hat\beta$들이 이루는 분포가 정규분포에 가깝다는 거야. 근거는 MLE의 일반 성질이야. $\hat\beta$는 NB likelihood가 가장 커지는 값, 즉 MLE(최대가능도추정)거든. likelihood는 관측값을 고정해 두고, parameter 후보가 그 관측값을 얼마나 그럴듯하게 만드는지 나타내는 값이야. 그리고 MLE는 데이터가 많아질수록 분포가 정규분포에 가까워져.

### SE와 기준 분포에 깔린 가정

앞에서 α를 Poisson보다 더 퍼지는 정도라고 했는데, 식으로 쓰면 분산 = μ + αμ²이야. μ는 모형이 그 sample에서 예상하는 평균 count(기대 count)야. α = 0이면 분산이 μ인 Poisson과 같아지니까, α는 분산 자체가 아니라 Poisson보다 추가로 퍼지는 정도인 셈이지. 그런데 SE를 계산할 때 DESeq2는 α와 size factor를 이미 정해진 값처럼 넣어(plug-in). 사실 α도 데이터에서 추정한 값이라 불확실한데, 기본 SE에는 그 불확실성이 들어 있지 않아. 그래서 여기서 나오는 SE와 신뢰구간은 "명목상(nominal)" 값이야.

Z를 비교하는 기준 분포도 기본은 표준정규분포야. `nbinomWaldTest(useT = TRUE)`를 주면 정규분포보다 꼬리가 두꺼운 t 분포와 비교하는데, sample이 적으면 두 결과가 크게 달라. 3절부터 쓰는 예제 데이터(12 sample, 계수 6개)의 Glucose–Ctrl 비교에서는 padj<0.05 유전자가 163개에서 35개로 줄었어([더 깊이 보기](#더-깊이-보기)).

### 95% 신뢰구간

$\hat\beta\pm1.96\times SE$가 nominal 95% 신뢰구간(CI)이야. 1.96은 표준정규분포에서 양쪽 끝 2.5%씩을 잘라내는 값이고. 유전자 A라면 0.977 ± 1.96 × 0.289 = [0.411, 1.544]야. CI가 0을 포함하지 않는 것과 p < 0.05는 같은 사건이야. 둘 다 "|Z| > 1.96"이라는 한 조건을 다르게 표현한 것이거든. 1.96으로 반올림한 값을 쓸 때 생기는 아주 작은 예외는 더 깊이 보기에 적어 두었어.

## 2. LFC가 같아도 p 값이 다른 이유는?

교재 9.2의 예를 그대로 계산해 볼게. LFC가 모두 1(2배)이고 SE만 다른 유전자 세 개야.

```r
tab <- data.frame(log2FC = 1, SE = c(0.2, 0.5, 1.0))
tab$Z <- tab$log2FC / tab$SE
tab$p <- 2 * pnorm(-abs(tab$Z))
tab$CI_low  <- tab$log2FC - qnorm(0.975) * tab$SE
tab$CI_high <- tab$log2FC + qnorm(0.975) * tab$SE
print(tab, digits = 4)
```

```
  log2FC  SE Z         p   CI_low CI_high
1      1 0.2 5 5.733e-07  0.60801   1.392
2      1 0.5 2 4.550e-02  0.02002   1.980
3      1 1.0 1 3.173e-01 -0.95996   2.960
```

SE가 커지면서 Z는 5, 2, 1로 줄고, p는 0.00000057에서 0.32까지 커져. CI가 0을 포함하는 건 SE가 1.0인 유전자 하나뿐이야.

![LFC가 1로 같고 SE만 다른 세 유전자의 점추정치와 95% 신뢰구간](../figures/05_same_lfc_different_se.png)

*세 유전자의 점은 모두 LFC = 1에 있어. SE가 커질수록 막대(95% CI)가 넓어지고, SE가 1.0인 유전자는 막대가 0선을 넘어서 p = 0.32가 돼.*

<details>
<summary>그림을 만든 코드</summary>

```r
suppressPackageStartupMessages(library(ggplot2))
fig <- transform(tab, yy = 3:1, lab = c("p = 5.7e-07", "p = 0.046", "p = 0.32"))    # SE 0.2 를 맨 위에
g1 <- ggplot(fig, aes(y = yy)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#8a8987") +
  annotate("text", x = 0, y = 3.45, label = "no change (LFC = 0)", hjust = -0.05, size = 3.2, colour = "#52514e") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.15, linewidth = 0.8, colour = "#2a78d6", orientation = "y") +
  geom_point(aes(x = log2FC), size = 3, colour = "#2a78d6") +
  geom_text(aes(x = 3.15, label = lab), hjust = 0, size = 3.5, colour = "#0b0b0b") +
  scale_x_continuous(limits = c(-1.2, 3.9), breaks = -1:3) +
  scale_y_continuous(breaks = fig$yy, labels = sprintf("SE = %.1f", fig$SE), limits = c(0.6, 3.6)) +
  labs(title = "Same log2 fold change, different uncertainty",
       subtitle = "Point = estimate (LFC = 1); bar = nominal 95% CI = LFC ± 1.96 × SE",
       x = "log2 fold change", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), plot.subtitle = element_text(colour = "#52514e"))
ggsave("../figures/05_same_lfc_different_se.png", g1, width = 7, height = 4.5, dpi = 150)
cat("saved:", file.exists("../figures/05_same_lfc_different_se.png"), "\n")
```

```
saved: TRUE
```

</details>

"몇 배 바뀌었나"와 "얼마나 정확히 추정했나"는 서로 다른 정보야. SE는 유전자마다 다른데, count가 작을수록(Poisson 잡음이 상대적으로 커서), dispersion α가 클수록(replicate끼리 더 들쭉날쭉해서) 커져. replicate가 적거나, 모형에 넣은 요인(batch, pair 등)이 많아 추정할 계수가 많을 때도 커지고.

dispersion의 영향은 유전자 A에서 바로 볼 수 있어. [03 노트](03_dispersion_estimation.md)에서 α를 세 가지 방법으로 구했거든. 보정 없는 NB MLE는 0.014786이고, 평균을 같은 데이터로 추정해서 생기는 과소추정을 줄여 주는 Cox-Reid 보정을 넣으면 0.025385야. 마지막은 prior(사전분포)를 붙인 값이야. prior는 데이터를 보기 전에 α가 어디쯤 있을지에 대한 분포인데, 여기서는 교재가 학습용으로 log α ~ N(log 0.08, 0.7²)를 지정했어. likelihood와 prior를 함께 고려했을 때 가장 그럴듯한 값(MAP)이 0.053147이야. 이제 count는 그대로 두고 α만 바꿔 넣어 볼게.

```r
for (a in c(0.014786, 0.025385, 0.053147)) {          # 03 노트: MLE, Cox-Reid, MAP
  d <- ddsA; dispersions(d) <- a
  r <- results(nbinomWaldTest(d, quiet = TRUE))
  cat(sprintf("alpha = %.6f   LFC = %.4f   SE = %.4f   Z = %.3f   p = %.2g\n",
              a, r$log2FoldChange, r$lfcSE, r$stat, r$pvalue))
}
```

```
alpha = 0.014786   LFC = 0.9773   SE = 0.1741   Z = 5.612   p = 2e-08
alpha = 0.025385   LFC = 0.9773   SE = 0.2122   Z = 4.605   p = 4.1e-06
alpha = 0.053147   LFC = 0.9773   SE = 0.2891   Z = 3.381   p = 0.00072
```

LFC는 세 줄 모두 0.9773으로 똑같아. 그런데 α가 0.015에서 0.053으로 커지자 SE는 0.174에서 0.289로 커지고, p는 2e-08에서 0.00072가 됐어. dispersion을 어떻게 추정하느냐가 p 값을 직접 바꾸는 거지.

## 3. 결과표 여섯 열은 각각 무엇을 말하나?

여기부터는 유전자 1000개짜리 연습 데이터를 써. DESeq2의 `makeExampleDESeqDataSet()`으로 만든 모의 데이터에 교재의 세 조건 이름을 붙인 거야. Ctrl, Starvation, Glucose가 4 sample씩이고, Glucose는 굶긴 세포에 glucose를 다시 넣은 조건(Starvation+Glucose)을 줄여 부른 이름이야. 같은 donor에서 세 조건을 얻었다고 보고 pair 1–4를 붙였기 때문에, design(모형에 넣을 요인을 적은 식)은 `~ pair + condition`이야.

그리고 데이터를 만든 방식 때문에 **Starvation–Ctrl의 참 효과는 모든 유전자에서 0이야.** Glucose–Ctrl에는 유전자마다 참 효과(`trueBeta`)가 있고, pair는 참 효과가 없지만 design에 넣어 추정해. 그러니 Starvation–Ctrl에서 유의하게 나오는 유전자는 전부 거짓 양성이야.

`results()`에 `contrast = c("condition", "Starvation", "Ctrl")`처럼 어떤 두 조건을 비교할지(contrast) 적어 주면 그 비교의 결과표가 나와.

```r
set.seed(1)
dds0 <- makeExampleDESeqDataSet(n = 1000, m = 16, betaSD = 1)   # A: 1-8, B: 9-16
keep <- c(1:8, 9:12)                                              # Ctrl=A1-4, Starvation=A5-8, Glucose=B1-4
cts  <- counts(dds0)[, keep]
meta <- data.frame(condition = factor(rep(c("Ctrl","Starvation","Glucose"), each = 4),
                                      levels = c("Ctrl","Starvation","Glucose")),
                   pair = factor(rep(1:4, 3)))
colnames(cts) <- rownames(meta) <- paste0(meta$condition, "_", meta$pair)
dds3 <- DESeqDataSetFromMatrix(cts, meta, ~ pair + condition)
dds_wald <- DESeq(dds3, quiet = TRUE)
res_wS <- results(dds_wald, contrast = c("condition", "Starvation", "Ctrl"))
res_wG <- results(dds_wald, name = "condition_Glucose_vs_Ctrl")
print(as.data.frame(res_wS[c("gene841", "gene420", "gene4"), ]), digits = 4)
print(mcols(res_wS)$description)
```

```
        baseMean log2FoldChange  lfcSE    stat   pvalue   padj
gene841    918.0        -0.3490 0.3258 -1.0710 0.284176 0.9988
gene420    373.6        -0.8044 0.3066 -2.6239 0.008693 0.9988
gene4      166.7         0.1590 0.3546  0.4484 0.653838 0.9988
[1] "mean of normalized counts for all samples"
[2] "log2 fold change (MLE): condition Starvation vs Ctrl"
[3] "standard error: condition Starvation vs Ctrl"
[4] "Wald statistic: condition Starvation vs Ctrl"
[5] "Wald test p-value: condition Starvation vs Ctrl"
[6] "BH adjusted p-values"
```

아래쪽 `description` 여섯 줄이 각 열이 무엇인지 적어 줘. 표를 읽다가 헷갈리면 이것부터 보면 돼.

| 열 | 뜻 | 헷갈리기 쉬운 점 |
|---|---|---|
| `baseMean` | 모든 sample의 정규화 count 평균 | 비교에 쓰지 않은 Glucose sample까지 들어가. gene841은 Glucose count가 커서 918.0인데, Ctrl과 Starvation만으로 새로 분석하면 174.8이야(6절). |
| `log2FoldChange` | 요청한 비교의 LFC | 분자가 Starvation, 분모가 Ctrl이야. 기본값(`betaPrior=FALSE`)에서는 LFC를 0 쪽으로 당기는 shrinkage를 하지 않은 MLE야. |
| `lfcSE` | LFC의 SE | `lfcShrink()`를 거치면 값과 뜻이 바뀌어([07](07_lfc_shrinkage_and_qc.md)). |
| `stat` | Wald Z = LFC / lfcSE | LRT 결과표에서는 전혀 다른 통계량 D야(4절). |
| `pvalue` | 원래 p 값 | "귀무가설이 참일 확률"도, "이 유전자가 거짓 양성일 확률"도 아니야. 효과가 0이라면 지금만큼 또는 더 극단적인 Z가 나올 확률이야. |
| `padj` | BH로 보정한 p 값 | 이 유전자 하나가 틀렸을 확률(posterior error probability)은 아니야. |

padj는 유전자 수천 개를 한꺼번에 검정하기 때문에 필요한 보정값이야. 목표는 FDR, 즉 발견 목록에서 거짓 발견이 차지하는 비율의 기대값을 조절하는 거야. DESeq2의 기본 보정법은 BH(Benjamini–Hochberg)이고, 자세한 계산은 [06](06_multiple_testing.md)에서 다뤄.

### CI가 0을 벗어났는데 padj 별표가 없는 이유

논문 그림에서 흔히 보는 조합이 있어. 막대는 95% CI이고, 별표는 padj<0.05인 경우지. 그런데 CI가 0을 벗어났는데도 별표가 없는 유전자가 생겨. Glucose–Ctrl 비교에서 세어 볼게.

```r
r  <- results(dds_wald, name = "condition_Glucose_vs_Ctrl", alpha = 0.05)
lo <- r$log2FoldChange - qnorm(0.975) * r$lfcSE; hi <- r$log2FoldChange + qnorm(0.975) * r$lfcSE
ci_ex <- (lo > 0 | hi < 0)
cat("non-NA pvalue:", sum(!is.na(r$pvalue)), "| CI excludes 0:", sum(ci_ex, na.rm = TRUE),
    " raw p<0.05:", sum(r$pvalue < 0.05, na.rm = TRUE), " padj<0.05:", sum(r$padj < 0.05, na.rm = TRUE), "\n")
cat("CI excludes 0 but padj>=0.05:", sum(ci_ex & !is.na(r$padj) & r$padj >= 0.05),
    "| CI excludes 0 but padj NA (independent filtering):", sum(ci_ex & is.na(r$padj), na.rm = TRUE), "\n")
cat("CI-excludes-0 <=> raw p<0.05:", identical(which(ci_ex), which(r$pvalue < 0.05)), "\n")
cat("padj<0.05 but CI includes 0:", sum(!is.na(r$padj) & r$padj < 0.05 & !ci_ex),
    "| padj >= pvalue for all non-NA:", all(r$padj >= r$pvalue, na.rm = TRUE), "\n")
```

```
non-NA pvalue: 999 | CI excludes 0: 277  raw p<0.05: 277  padj<0.05: 163
CI excludes 0 but padj>=0.05: 103 | CI excludes 0 but padj NA (independent filtering): 11
CI-excludes-0 <=> raw p<0.05: TRUE
padj<0.05 but CI includes 0: 0 | padj >= pvalue for all non-NA: TRUE
```

CI가 0을 벗어난 277개와 raw p<0.05인 277개는 정확히 같은 유전자야. 그런데 padj<0.05는 그중 163개뿐이야.

CI는 유전자 하나만 보고 만든 구간이지만, padj는 여러 유전자를 한꺼번에 검정했다는 사실을 반영해. 여기서는 p 값이 있는 999개 중 independent filtering을 통과한 766개를 함께 보정했어. 그래서 CI는 0을 벗어났는데 padj<0.05는 아닌 유전자가 114개 생겨. 103개는 padj가 0.05 이상이고, 11개는 independent filtering(평균 count가 낮은 유전자를 보정 대상에서 빼는 절차)에 걸려 padj가 NA야.

반대 방향은 성립해. BH padj는 언제나 raw p보다 크거나 같으니까, padj<0.05인 163개는 모두 CI가 0을 벗어나. 그림에 nominal CI 막대와 padj 별표를 함께 그린다면, 둘의 기준이 다르다는 걸 legend에 적어 둬.

### "|LFC| > 1인 유전자"를 제대로 고르는 법

"padj<0.05이고 |LFC|>1인 유전자"라는 기준을 많이 써. 그런데 이건 "효과가 0인가"를 검정한 다음 추정치가 1을 넘는 것만 남긴 필터야. "진짜 효과가 2배보다 크다"를 검정한 건 아니지. 그렇게 주장하고 싶다면 귀무가설 자체를 |β| ≤ 1로 바꿔서 검정해야 해. `results()`의 `lfcThreshold = 1`이 바로 그 검정이야.

```r
r0 <- results(dds_wald, name = "condition_Glucose_vs_Ctrl", alpha = 0.05)
r1 <- results(dds_wald, name = "condition_Glucose_vs_Ctrl", alpha = 0.05,
              lfcThreshold = 1, altHypothesis = "greaterAbs")
cat("필터 padj<0.05 & |LFC|>1:", sum(r0$padj < 0.05 & abs(r0$log2FoldChange) > 1, na.rm = TRUE),
    "| lfcThreshold = 1 검정 padj<0.05:", sum(r1$padj < 0.05, na.rm = TRUE), "\n")
gg <- c("gene8", "gene20", "gene38")
print(data.frame(LFC = r0[gg, "log2FoldChange"], SE = r0[gg, "lfcSE"],
                 p_vs0 = r0[gg, "pvalue"], padj_vs0 = r0[gg, "padj"],
                 p_vs1 = r1[gg, "pvalue"], padj_vs1 = r1[gg, "padj"], row.names = gg), digits = 3)
```

```
필터 padj<0.05 & |LFC|>1: 156 | lfcThreshold = 1 검정 padj<0.05: 18
         LFC    SE    p_vs0 padj_vs0    p_vs1 padj_vs1
gene8  -1.84 0.601 2.18e-03 1.49e-02 0.080679   0.6182
gene20  1.83 0.389 2.44e-06 5.84e-05 0.016130   0.2403
gene38  2.99 0.610 9.53e-07 2.71e-05 0.000554   0.0327
```

필터로는 156개가 남지만, "|LFC|가 1보다 크다"를 직접 검정하면 18개만 남아.

차이는 gene8에서 잘 보여. 추정 LFC가 −1.84이고 padj도 0.015로 작아서 필터는 통과해. 하지만 SE가 0.60이라 참값이 1 안쪽일 가능성을 배제하지 못해서, |β| ≤ 1에 대한 p는 0.081, padj는 0.62야. SE가 큰 유전자는 LFC 추정치가 커도 이 검정을 통과하지 못하는 거지. 18개는 모두 필터 156개 안에 있었어. 필터만 통과한 나머지 138개 중 135개는 이 검정의 padj가 0.05 이상이었고, 3개는 이 검정의 independent filtering에서 빠져 padj가 NA였어.

기준값 T와 검정 방향(`altHypothesis`)은 결과를 보기 전에 정해 둬야 해. 같은 T = 1이라도 DESeq2가 제공하는 계산 방식에 따라 결과가 18개와 12개로 갈렸거든(더 깊이 보기). 그래서 교재도 어떤 버전과 옵션을 썼는지 기록해 두라고 해.

정리하면 LFC, SE, raw p, CI, padj는 서로 다른 질문에 답해. LFC는 얼마나 큰지, SE는 얼마나 정밀한지, raw p와 CI는 이 유전자 하나만 봤을 때 0과 양립하는지, padj는 목록 전체의 FDR 아래에서도 남는지를 말해 줘. 그리고 통계적으로 유의하다고 해서 생물학적으로 중요하다거나 기전이 증명됐다는 뜻은 아니야.

## 4. 세 그룹 중 "어딘가" 다른지는 어떻게 묻나?

세 그룹 design에는 condition 계수가 두 개 있어. Starvation–Ctrl($b_S$)과 Glucose–Ctrl($b_G$)이지. 여기서 $b$는 자연로그 단위 계수이고, log2 단위로 바꾼 $\beta=b/\log 2$가 결과표의 LFC야. 그런데 "세 조건의 평균이 모두 같은가"는 계수 하나로는 답할 수 없어. 두 계수가 동시에 0인지 물어야 하거든. 이럴 때 쓰는 게 LRT(가능도비 검정)야. 계수를 뺀 모형과 넣은 모형의 likelihood를 비교해서 여러 계수를 한꺼번에 검정하는 방법이야.

### 먼저 유전자 A로: 계수를 빼면 likelihood가 얼마나 떨어지나

유전자 A에 두 모형을 맞춰 볼게.

- full: Ctrl과 Starvation이 각자 평균(106.667과 210)을 가져.
- reduced: condition 계수를 빼서, 여섯 sample이 공통 평균 158.333 하나를 써.

각 모형에서 log-likelihood를 구해. 관측한 count 여섯 개가 그 모형 아래에서 나올 NB 확률의 log를 더한 값인데, 클수록 모형이 데이터를 잘 설명한다는 뜻이야. reduced는 full에서 계수 하나를 0으로 고정한 특수한 경우라서, full의 log-likelihood가 reduced보다 작을 수는 없어. 문제는 얼마나 더 크냐는 거지.

```r
ddsA_lrt <- nbinomLRT(ddsA, reduced = ~ 1, quiet = TRUE)        # 같은 alpha, 같은 size factor
resA_lrt <- results(ddsA_lrt)
a <- 0.053147
mu_full <- rep(c(mean(k[1:3]), mean(k[4:6])), each = 3)        # full: 그룹마다 평균
mu_red  <- rep(mean(k), 6)                                      # reduced: 공통 평균 하나
l_full <- sum(dnbinom(k, mu = mu_full, size = 1 / a, log = TRUE))
l_red  <- sum(dnbinom(k, mu = mu_red,  size = 1 / a, log = TRUE))
cat("l_full =", round(l_full, 3), "  l_reduced =", round(l_red, 3), "  D = 2(l_full - l_reduced) =", round(2 * (l_full - l_red), 4), "\n")
cat("DESeq2 LRT : stat =", round(resA_lrt$stat, 4), "  p =", signif(resA_lrt$pvalue, 4), "\n")
cat("Wald 비교  : Z^2  =", round(resA$stat^2, 4), "  p =", signif(resA$pvalue, 4), "\n")
```

```
l_full = -28.182   l_reduced = -33.829   D = 2(l_full - l_reduced) = 11.2933
DESeq2 LRT : stat = 11.2933   p = 0.0007779
Wald 비교  : Z^2  = 11.4307   p = 0.0007224
```

손으로 구한 D가 DESeq2의 `stat`과 같고, Wald의 $Z^2$와도 가까워. `nbinomLRT()`는 `nbinomWaldTest()` 대신 LRT를 하는 검정 단계 함수야.

condition 계수를 빼자 log-likelihood가 −28.182에서 −33.829로 떨어졌어. 이 차이를 두 배 한 것이 LRT 통계량이야.

$$D=2\,(\ell_{full}-\ell_{reduced})\ \overset{H_0}{\approx}\ \chi^2_{df}$$

- $\ell_{full}$, $\ell_{reduced}$: 두 모형의 최대 log-likelihood.
- $D$: 계수를 뺐을 때 잃은 설명력. 결과표의 `stat` 열이야.
- $H_0$: 귀무가설. 뺀 계수가 모두 0이라는 가설이야.
- $\chi^2_{df}$: 자유도가 df인 카이제곱 분포. 귀무가설이 맞으면 D가 대략 이 분포를 따라.
- $df$: 뺀 계수의 수, 즉 두 design matrix의 열 수 차이야. design matrix는 각 sample이 어떤 조건·pair에 속하는지 숫자로 적은 표야.

유전자 A에 대입하면 $D=2\{-28.182-(-33.829)\}=11.29$이고, 뺀 계수가 하나라서 df = 1이야. $\chi^2_1$에서 11.29 이상이 나올 확률이 p = 0.00078인데, Wald의 $Z^2=11.43$, p = 0.00072와 거의 같아. 두 그룹 비교에서는 Wald와 LRT가 같은 질문을 다른 방법으로 묻는 셈이야.

"여러 평균이 모두 같은가"를 묻는다는 점에서 질문의 모양은 ANOVA와 닮았어. 다만 분포가 정규가 아니라 NB이고, 통계량이 F가 아니라 D이며, 기준 분포가 카이제곱이야. 그러니 raw count에 ANOVA를 먼저 돌릴 필요는 없어.

### 세 그룹: `~ pair + condition` 대 `~ pair`

이제 연습 데이터의 세 그룹으로 넘어갈게. full은 `~ pair + condition`, reduced는 `~ pair`야. 귀무가설은 $b_S=b_G=0$, 즉 pair 차이를 감안하고 나면 세 조건의 평균이 같다는 거야. pair는 검정 대상이 아니니까 reduced에도 남겨 둬. full의 계수는 6개(절편, pair 3개, Starvation, Glucose)이고 reduced는 4개라서 df = 2야.

```r
dds_lrt <- DESeq(dds3, test = "LRT", reduced = ~ pair, quiet = TRUE)
res_l <- results(dds_lrt)
nF <- ncol(attr(dds_lrt, "modelMatrix")); nR <- ncol(attr(dds_lrt, "reducedModelMatrix"))
cat("df =", nF, "-", nR, "=", nF - nR, "\n")
cat("dispersion이 Wald object와 같은가:", identical(dispersions(dds_wald), dispersions(dds_lrt)), "\n")
g3 <- c("gene841", "gene420", "gene4")
print(data.frame(Z_Starvation = res_wS[g3, "stat"], Z_Glucose = res_wG[g3, "stat"],
                 LRT_D = res_l[g3, "stat"], LRT_p = res_l[g3, "pvalue"], row.names = g3), digits = 4)
print(counts(dds_wald)["gene841", ])
```

```
df = 6 - 4 = 2
dispersion이 Wald object와 같은가: TRUE
        Z_Starvation Z_Glucose  LRT_D     LRT_p
gene841      -1.0710   11.2459 184.30 9.540e-41
gene420      -2.6239    5.2985  64.77 8.628e-15
gene4         0.4484    0.7308   0.53 7.672e-01
      Ctrl_1       Ctrl_2       Ctrl_3       Ctrl_4 Starvation_1 Starvation_2 Starvation_3 Starvation_4
         283          197          148          177          203          131          167          144
   Glucose_1    Glucose_2    Glucose_3    Glucose_4
        1379         2511         3770         2385
```

gene841의 D = 184.3은 두 계수를 함께 뺀 효과야. count를 보면 그 대부분은 Glucose 쪽 증가에서 와.

dispersion은 full design으로 한 번 추정한 값을 두 모형이 그대로 나눠 써. LRT는 reduced 모형용 dispersion을 따로 추정하지 않고 GLM만 두 번 적합하거든(`nbinomLRT` 소스, 더 깊이 보기). 출력에서 `identical`이 TRUE인 것도 LRT object의 dispersion이 Wald object와 한 자리도 다르지 않다는 뜻이야.

### 두 그룹이면 Wald와 LRT는 늘 같은 답을 주나?

유전자 A에서는 두 p가 0.00072와 0.00078로 비슷했어. 유전자 1000개에서도 그럴까? 연습 데이터에서 Ctrl과 Glucose만으로 새 object를 만들고(df = 1), 두 검정의 p를 나란히 놓아 볼게.

```r
keepG <- meta$condition %in% c("Ctrl", "Glucose")
metaG <- droplevels(meta[keepG, , drop = FALSE])
ddsG  <- DESeq(DESeqDataSetFromMatrix(cts[, rownames(metaG)], metaG, ~ pair + condition), quiet = TRUE)
ddsG_lrt <- nbinomLRT(ddsG, reduced = ~ pair, quiet = TRUE)
wG <- results(ddsG, name = "condition_Glucose_vs_Ctrl"); lG <- results(ddsG_lrt)
ok <- !is.na(wG$pvalue) & !is.na(lG$pvalue)
cat("genes:", sum(ok), "| cor(-log10 p):", round(cor(-log10(wG$pvalue[ok]), -log10(lG$pvalue[ok])), 4), "\n")
cat("p<0.05  Wald:", sum(wG$pvalue[ok] < 0.05), " LRT:", sum(lG$pvalue[ok] < 0.05),
    " both:", sum(wG$pvalue[ok] < 0.05 & lG$pvalue[ok] < 0.05), "\n")
cat("share of genes with D < Z^2:", round(mean(lG$stat[ok] < wG$stat[ok]^2), 3), "\n")
```

```
genes: 995 | cor(-log10 p): 0.9977
p<0.05  Wald: 255  LRT: 229  both: 229
share of genes with D < Z^2: 0.93
```

두 p의 상관은 0.998이야. 다만 p<0.05인 유전자 수는 Wald 255개, LRT 229개로 조금 달라.

![두 그룹 비교에서 유전자별 Wald p와 LRT p를 비교한 산점도](../figures/05_wald_vs_lrt_pvalues.png)

*점 하나가 유전자 하나야. 대부분 대각선 근처에 있어서 두 검정이 거의 같은 답을 줘. 오른쪽 확대 그림의 주황색 점 26개는 Wald로만 p<0.05이고 LRT로는 0.05를 넘는 유전자야.*

<details>
<summary>그림을 만든 코드</summary>

```r
wo <- ok & wG$pvalue < 0.05 & lG$pvalue >= 0.05                  # Wald로만 p<0.05
pp <- data.frame(x = -log10(wG$pvalue[ok]), y = -log10(lG$pvalue[ok]),
                 grp = ifelse(wo[ok], "p < 0.05 by Wald only", "other genes"))
zoom <- subset(pp, x < 3 & y < 3)
pp$panel <- "All genes"; zoom$panel <- "Zoom: both p > 0.001"
both <- rbind(pp, zoom); both$panel <- factor(both$panel, levels = c("All genes", "Zoom: both p > 0.001"))
thr <- data.frame(panel = factor("Zoom: both p > 0.001", levels = levels(both$panel)), v = -log10(0.05))
g2 <- ggplot(both, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#8a8987") +
  geom_vline(data = thr, aes(xintercept = v), linetype = "dotted", colour = "#8a8987") +
  geom_hline(data = thr, aes(yintercept = v), linetype = "dotted", colour = "#8a8987") +
  geom_point(aes(colour = grp), size = 1.3, alpha = 0.7) +
  scale_colour_manual(values = c("other genes" = "#2a78d6", "p < 0.05 by Wald only" = "#eb6834"), name = NULL) +
  facet_wrap(~ panel, scales = "free") +
  labs(title = "Two-group comparison (df = 1): Wald vs LRT",
       subtitle = sprintf("Ctrl vs Glucose, %d genes. Dashed: equal p-values. Dotted: p = 0.05", sum(ok)),
       x = "-log10(p), Wald", y = "-log10(p), LRT") +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA), panel.grid.minor = element_blank(),
        legend.position = "top", legend.justification = "left", plot.subtitle = element_text(colour = "#52514e"))
ggsave("../figures/05_wald_vs_lrt_pvalues.png", g2, width = 7, height = 4.5, dpi = 150)
cat("saved:", file.exists("../figures/05_wald_vs_lrt_pvalues.png"), "| Wald-only genes:", sum(wo),
    "| their Wald p:", signif(range(wG$pvalue[wo]), 3), "| LRT p:", signif(range(lG$pvalue[wo]), 3), "\n")
for (b in list(c(0, 2), c(2, 4), c(4, 8), c(8, Inf))) {                # |Z| 구간별 p 차이
  s <- abs(wG$stat[ok]) >= b[1] & abs(wG$stat[ok]) < b[2]
  cat(sprintf("|Z| in [%g, %g): n = %3d  median log10(p_LRT / p_Wald) = %.3f\n",
              b[1], b[2], sum(s), median(log10(lG$pvalue[ok][s] / wG$pvalue[ok][s]))))
}
```

```
saved: TRUE | Wald-only genes: 26 | their Wald p: 0.0258 0.0495 | LRT p: 0.0502 0.132
|Z| in [0, 2): n = 749  median log10(p_LRT / p_Wald) = 0.016
|Z| in [2, 4): n = 195  median log10(p_LRT / p_Wald) = 0.128
|Z| in [4, 8): n =  47  median log10(p_LRT / p_Wald) = 0.400
|Z| in [8, Inf): n =   4  median log10(p_LRT / p_Wald) = 3.658
```

</details>

두 검정은 대부분 같은 결론을 내. LRT로 p<0.05인 229개는 모두 Wald로도 p<0.05였어. 이 데이터에서는 유전자의 93%에서 D가 $Z^2$보다 조금 작았으니, LRT 쪽이 조금 더 보수적이었던 셈이야. 다만 이건 이 모의 데이터에서 본 결과이지 일반 법칙으로 확인한 건 아니야. 두 p의 차이는 |Z|가 클수록 벌어져서, |Z|가 8 이상인 유전자 4개에서는 p가 수천 배씩 차이 났어(중앙값 기준, 그림 코드의 출력). 그래도 둘 다 아주 작은 값이라 결론은 같아. 결론이 갈린 26개는 Wald p 0.026–0.0495, LRT p 0.050–0.13으로 경계 근처에 있었어. 어느 쪽이 "맞는" p라기보다, 두 근사가 경계 근처에서 조금 다른 답을 준 거야.

### LRT 결과표의 함정

LRT 결과표에도 `log2FoldChange` 열이 있어. 그래서 옆의 `pvalue`도 그 LFC에 대한 p로 읽기 쉬운데, 나도 처음엔 그렇게 생각했어. 하지만 그렇게 읽으면 안 돼.

```r
res_lS <- results(dds_lrt, contrast = c("condition", "Starvation", "Ctrl"))
res_lG <- results(dds_lrt, name = "condition_Glucose_vs_Ctrl")
cat("contrast를 바꿔도 pvalue가 그대로인가:", identical(res_l$pvalue, res_lS$pvalue), identical(res_l$pvalue, res_lG$pvalue), "\n")
cat("LRT 표의 LFC 열:", mcols(res_l)$description[2], "\n")
cat("LRT 표의 p 열  :", mcols(res_l)$description[5], "\n")
res_lS_wald <- results(dds_lrt, contrast = c("condition", "Starvation", "Ctrl"), test = "Wald")
cat("test = 'Wald'로 뽑은 p == Wald object의 p:", isTRUE(all.equal(res_lS_wald$pvalue, res_wS$pvalue)), "\n")
sig <- function(r) !is.na(r$padj) & r$padj < 0.05
cat("padj<0.05  Wald Glucose-Ctrl:", sum(sig(res_wG)), " LRT:", sum(sig(res_l)),
    " Wald만:", sum(sig(res_wG) & !sig(res_l)), " LRT만:", sum(sig(res_l) & !sig(res_wG)), "\n")
```

```
contrast를 바꿔도 pvalue가 그대로인가: TRUE TRUE
LRT 표의 LFC 열: log2 fold change (MLE): condition Glucose vs Ctrl
LRT 표의 p 열  : LRT p-value: '~ pair + condition' vs '~ pair'
test = 'Wald'로 뽑은 p == Wald object의 p: TRUE
padj<0.05  Wald Glucose-Ctrl: 163  LRT: 146  Wald만: 37  LRT만: 20
```

`contrast`를 바꿔도 `pvalue`는 한 글자도 바뀌지 않아. p 열의 설명도 "LRT p-value: '~ pair + condition' vs '~ pair'"라고 되어 있고.

LRT 결과표의 `pvalue`는 "Starvation이든 Glucose든 어딘가 Ctrl과 다른가"에 대한 p 하나야. 이렇게 여러 비교를 한꺼번에 묻는 p를 omnibus p라고 불러. `log2FoldChange`와 `lfcSE`는 full 모형의 계수 하나를 보여 주려고 붙여 둔 것이고, 기본으로는 마지막 계수(Glucose vs Ctrl)가 붙어. 특정 두 조건을 비교한 p가 필요하면 `test = "Wald"`를 주거나 Wald object를 따로 둬. 표에는 어느 검정의 p인지도 적어 두고.

padj<0.05로 세면 Wald(Glucose–Ctrl) 163개와 LRT 146개는 어느 쪽도 다른 쪽을 포함하지 않아(Wald만 37개, LRT만 20개). 일반적으로 LRT는 증거를 두 계수(df 2)에 나눠 쓰기 때문에, 한 비교에 몰린 효과에는 Wald보다 힘이 약할 수 있어. 반대로 두 비교에 조금씩 퍼진 효과는 LRT가 더 잘 잡을 수 있고.

하지만 이 데이터에서 LRT로만 잡힌 20개는 그런 퍼진 참 효과일 수 없어. Starvation–Ctrl의 참 효과가 모든 유전자에서 0이니까. 이 20개는 Glucose–Ctrl Wald padj의 중앙값이 0.105로 경계 근처에 있었고, 그중 5개는 Starvation 쪽 잡음(Starvation–Ctrl raw p<0.05)이 LRT 통계량에 더해진 경우였어. 여기에 두 검정이 서로 다른 p 분포에 BH를 적용하고 filtering 기준도 다르다는 점이 겹쳤어(더 깊이 보기).

그러면 LRT를 먼저 하고, 유의한 유전자에서만 pairwise 비교를 해야 할까? 그럴 필요는 없어. 질문이 처음부터 "Starvation–Ctrl"이라면 그 Wald 검정을 바로 하면 돼. LRT가 유의했다고 해서 뒤따르는 pairwise 비교의 오류율이 자동으로 보호되는 것도 아니고([06](06_multiple_testing.md)). LRT는 "어딘가 차이가 있는가", Wald는 "이 특정 차이가 있는가"를 물을 때 맞는 도구야.

### Wald가 크게 빗나가는 경우: 한 그룹 count가 모두 0

Wald와 LRT가 극단적으로 갈리는 유전자도 있어. gene818은 Starvation count 네 개가 모두 0이야(Ctrl 1, 13, 0, 0 / Glucose 0, 0, 4, 0). 그러면 Starvation 평균의 추정치는 0, 즉 LFC는 −∞로 가야 해. 계산은 도중에 멈추는데, 멈춘 자리의 LFC −21.6과 SE로 만든 Wald p가 $4\times10^{-14}$가 됐어. Starvation–Ctrl에서 padj<0.05인 유일한 유전자가 바로 이거야. 참 효과가 0인 비교에서 나온 거짓 양성이지. 같은 유전자의 omnibus LRT(두 계수를 함께 검정)는 D = 0.28, p = 0.87이야. likelihood는 계수가 −∞로 가도 유한한 값에 다가가기 때문에 D는 계산이 멈춘 위치에 거의 흔들리지 않지만, Wald Z는 멈춘 위치의 LFC와 SE에 그대로 달려 있거든. 이렇게 count가 작은 유전자는 `stat`보다 count 자체를 먼저 봐.

## 5. "A에서는 유의, B에서는 아니다"는 차이의 증거일까?

genotype A와 B에서 각각 처리(T)와 대조(C)를 비교했다고 해 볼게. A에서는 p<0.05, B에서는 p>0.05가 나왔어. 이때 "처리 효과가 genotype에 따라 다르다"고 말해도 될까? 모든 참 효과가 0인 모의 데이터로 확인해 봤어. genotype, condition, 둘의 상호작용 모두 참 효과가 없으니, 여기서 나오는 "차이"는 전부 거짓이야.

```r
set.seed(2); d <- makeExampleDESeqDataSet(n = 500, m = 16)        # betaSD 기본값 0: 모든 참 효과 0
d$genotype  <- factor(rep(c("A", "B"), each = 8))
d$condition <- factor(rep(c("C", "T"), 8))
design(d) <- ~ genotype + condition + genotype:condition
dw  <- DESeq(d, quiet = TRUE)
rA  <- results(dw, name = "condition_T_vs_C")                                        # A 안의 T-C
rB  <- results(dw, contrast = list(c("condition_T_vs_C", "genotypeB.conditionT")))    # B 안의 T-C
int <- results(dw, name = "genotypeB.conditionT")                                     # 두 효과의 차이
one <- xor(rA$pvalue < 0.05, rB$pvalue < 0.05)
cat("한쪽에서만 p<0.05:", sum(one, na.rm = TRUE),
    "| 그중 interaction p<0.05:", sum(one & int$pvalue < 0.05, na.rm = TRUE),
    "| 전체 interaction p<0.05:", sum(int$pvalue < 0.05, na.rm = TRUE), "/", sum(!is.na(int$pvalue)), "\n")
cat("interaction LFC == (B의 T-C) - (A의 T-C):",
    isTRUE(all.equal(int$log2FoldChange, rB$log2FoldChange - rA$log2FoldChange)), "\n")
```

```
한쪽에서만 p<0.05: 60 | 그중 interaction p<0.05: 21 | 전체 interaction p<0.05: 27 / 500
interaction LFC == (B의 T-C) - (A의 T-C): TRUE
```

참 효과가 하나도 없는데도 "한쪽에서만 유의"한 유전자가 60개나 나와.

"한쪽만 유의하니 효과가 다르다"고 읽으면 60개의 거짓 차이를 보고하게 되는 거지. 한쪽 p가 0.04이고 다른 쪽이 0.06이기만 해도 "유의/비유의"로 갈리니까. 유의성이 갈렸다는 건 차이의 증거가 아니야. 물어야 할 것은 두 효과의 차이 그 자체야.

$$\delta_{int}=\delta_B-\delta_A=(\log_2 q_{B,T}-\log_2 q_{B,C})-(\log_2 q_{A,T}-\log_2 q_{A,C})$$

- $q_{B,T}$: genotype B, 처리 T의 평균 발현. 나머지 $q$도 같은 식이야.
- $\delta_A$, $\delta_B$: 각 genotype 안에서의 처리 효과(LFC).
- $\delta_{int}$: 두 효과의 차이. DESeq2에서는 `genotypeB.conditionT` 계수인데, 이게 interaction(상호작용) 항이야.

올바른 검정은 이 계수를 그 SE로 나눈 Wald 검정이야. 위 시뮬레이션에서 interaction p<0.05는 500개 중 27개(5.4%)로, 참 효과가 없을 때 기대하는 5% 근처였어. "한쪽만 유의" 60개 중에서는 21개였는데, 이것도 모두 거짓 양성이야.

두 추정치를 더하거나 뺀 값의 SE를 $\sqrt{SE_1^2+SE_2^2}$로 구하는 공식은 두 추정치가 서로 상관이 없을 때만 맞아. 상관되어 있으면 공분산(두 추정치가 함께 흔들리는 정도) 항이 더 필요하거든(더 깊이 보기). 다행히 interaction 계수의 `lfcSE`는 공분산까지 포함해 계산된 값이라 따로 합칠 필요가 없어. LRT로 검정하려면 reduced 모형에서 interaction 항만 빼고 genotype과 condition은 남기면 돼(df = 1). 다만 두 방법의 p가 정확히 같지는 않아.

## 6. 세 번째 그룹을 넣고 빼면 Ctrl–Starvation 결과가 왜 바뀌나?

같은 Ctrl–Starvation 비교라도 얻는 방법은 두 가지야. Ctrl, Starvation, Glucose 세 그룹을 함께 적합한 object에서 `contrast`로 뽑을 수도 있고, Ctrl과 Starvation만 raw count에서 새 object로 만들어 처음부터 분석할 수도 있어. 같은 8 sample의 같은 count인데도 결과는 달라. 무엇이 다른지, 그리고 이미 적합한 object를 잘라 쓰는 흔한 방법이 왜 문제인지 같이 볼게.

```r
sub <- dds_wald[, dds_wald$condition != "Glucose"]                  # 적합된 3 그룹 object를 자른다
cat("(가) 그대로 DESeq():", tryCatch({ DESeq(sub, quiet = TRUE); "ok" }, error = function(e) conditionMessage(e)), "\n")
sub$condition <- droplevels(sub$condition)
msgs <- character()
sub_fit <- withCallingHandlers(DESeq(sub), message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
cat("(나) droplevels 후 첫 메시지:", trimws(msgs[1]), "\n")
keep2 <- meta$condition %in% c("Ctrl", "Starvation")                  # (다) 교재 12.4: raw count로 새 object
meta2 <- droplevels(meta[keep2, , drop = FALSE]); cts2 <- cts[, rownames(meta2), drop = FALSE]
dds2  <- DESeq(DESeqDataSetFromMatrix(cts2, meta2, ~ pair + condition), quiet = TRUE)
r3 <- results(dds_wald, contrast = c("condition", "Starvation", "Ctrl"), alpha = 0.05)
r2 <- results(dds2,     contrast = c("condition", "Starvation", "Ctrl"), alpha = 0.05)
cmp <- function(f) c(three_group = format(signif(unname(f(dds_wald, r3)), 4)), two_group = format(signif(unname(f(dds2, r2)), 4)))
print(noquote(rbind(
  size_factor_Ctrl_1 = cmp(function(d, r) sizeFactors(d)[["Ctrl_1"]]),
  residual_df        = cmp(function(d, r) ncol(d) - ncol(attr(d, "dispModelMatrix"))),
  median_dispersion  = cmp(function(d, r) median(dispersions(d), na.rm = TRUE)),
  median_lfcSE       = cmp(function(d, r) median(r$lfcSE, na.rm = TRUE)),
  baseMean_gene841   = cmp(function(d, r) r["gene841", "baseMean"]),
  filterThreshold    = cmp(function(d, r) metadata(r)$filterThreshold))))
cat("cor(LFC 3 그룹, LFC 2 그룹):", round(cor(r3$log2FoldChange, r2$log2FoldChange, use = "complete"), 3), "\n")
```

```
(가) 그대로 DESeq(): full model matrix is less than full rank
(나) droplevels 후 첫 메시지: using pre-existing size factors
                   three_group two_group
size_factor_Ctrl_1 1.037       1.042
residual_df        6           3
median_dispersion  0.4143      0.5139
median_lfcSE       0.7119      0.7852
baseMean_gene841   918         174.8
filterThreshold    0.07903     0.1194
cor(LFC 3 그룹, LFC 2 그룹): 0.82
```

(가)는 오류로 멈추고, (나)는 3 그룹의 size factor를 그대로 물려받아. 아래 표를 보면 size factor부터 filtering 기준까지 모든 줄이 달라.

먼저 잘라 쓰기의 함정부터 볼게. (가) 적합된 3 그룹 object를 `dds[, condition != "Glucose"]`로 자르면 factor level `Glucose`가 이름만 남아. 그러면 design matrix에 전부 0인 열이 생기고 `DESeq()`는 "full model matrix is less than full rank"로 멈춰. (나) `droplevels()`로 빈 level을 지우면 돌아가기는 해. 하지만 "using pre-existing size factors"라며 3 그룹에서 계산한 size factor를 그대로 써. (다) 교재 12.4처럼 raw count와 metadata로 새 object를 만들어야 처음부터 다시 추정돼.

(다)와 3 그룹 결과를 나란히 놓으면 이런 차이가 있어.

- size factor부터 달라. size factor를 구할 때 기준으로 삼는 유전자별 기하평균(geometric mean)을 12 sample로 구하느냐 8 sample로 구하느냐가 다르니까([02](02_negative_binomial.md)).
- dispersion도 유전자별 추정치, 평균에 따른 경향(trend), 최종 α가 모두 바뀌어.
- 잔여 자유도가 6 대 3이라, 3 그룹 쪽이 dispersion을 추정할 정보가 더 많아.
- 그래서 SE와 p가 바뀌어. LFC는 대체로 비슷하지만(gene841 −0.349 대 −0.318), count가 거의 0인 gene818 같은 유전자는 예외야(더 깊이 보기). 이 비교는 참 효과가 0이라 LFC가 대부분 잡음인데도 두 분석의 LFC 상관이 0.82야.
- `baseMean`도 포함한 sample이 달라서 바뀌어. gene841은 918.0 대 174.8이야.
- independent filtering 기준(`filterThreshold`)과 padj가 NA인 유전자 수도 바뀌어([06](06_multiple_testing.md)).

이 데이터의 Glucose는 Ctrl, Starvation과 같은 방식으로 만든 모의 그룹이라 변동성이 비슷해. 이런 replicate가 더 들어오면 dispersion 추정에 정보를 더 보태 주는데(교재 12.2), 3 그룹 쪽 median SE가 0.712로 2 그룹의 0.785보다 작은 것도 그 방향과 맞아. DESeq2는 유전자마다 dispersion α를 하나만 추정해서 포함된 모든 그룹에 써. 그렇다고 그룹별 분산이 같다는 뜻은 아니야. 평균이 다르면 분산 μ + αμ²도 다르니까. 다만 한 그룹만 replicate 사이 변동이 유독 크면, 그 영향이 공통 α를 통해 Ctrl–Starvation 비교의 SE에도 들어가.

그럼 함께 적합할까, 두 그룹만 따로 분석할까? DESeq2 공식 문서(vignette FAQ)는 보통 모든 그룹을 함께 적합하고 `contrast`로 필요한 비교를 뽑으라고 권해. 같은 FAQ에는 PCA 같은 탐색 분석에서 한 그룹의 그룹 내 변동이 극단적으로 클 때는 두 그룹만 따로 분석하는 편이 더 민감할 수 있다는 설명도 있어. 중요한 건 이 판단을 결과를 보기 전에 탐색 분석으로 내린다는 점이야. 결과가 더 좋게 나오는 쪽을 나중에 고르는 것과는 다르지.

두 분석이 왜 다른지 진단할 때는 순서가 있어. 먼저 같은 유전자 집합인지, size factor를 새로 추정했는지 물려받았는지, outlier와 filtering 옵션이 같은지를 맞추거나 기록해 둬. 그다음 LFC → SE → raw p → padj 순서로 비교해. padj만 비교하면 모형이 바뀐 효과(SE, p)와 보정이 바뀐 효과(filterThreshold, 검정 수)가 섞여 보이거든. 결국 2 그룹과 3 그룹 결과의 차이는 비교 수가 늘어서만 생기는 게 아니야. 모형이 학습한 정규화와 dispersion, 그리고 padj를 함께 계산하는 유전자 집합을 나눠서 확인해야 해.

## 정리

- Wald 검정은 $Z=\hat\beta/SE$를 표준정규분포와 비교해. 정규분포를 가정하는 건 count가 아니라 추정치 $\hat\beta$이고, 유전자 A는 Z = 3.381, p = 0.00072였어.
- 같은 LFC라도 SE가 다르면 p가 크게 달라. SE는 count 크기, dispersion, replicate 수에 달려 있어서, 유전자 A도 α를 0.015에서 0.053으로 바꾸자 SE가 0.174에서 0.289로 커졌지.
- CI가 0을 벗어나는 건 raw p<0.05와 같은 말이지 padj<0.05와는 달라(277개 대 163개). "|LFC|>1 필터"와 `lfcThreshold = 1` 검정도 달라(156개 대 18개).
- LRT는 계수를 뺀 모형과 넣은 모형의 likelihood 차이 $D=2(\ell_{full}-\ell_{reduced})$를 카이제곱 분포와 비교해. 두 모형은 같은 dispersion을 쓰고, 두 그룹 비교에서는 Wald와 대부분 같은 답을 줘. 다만 LRT 결과표의 `pvalue`는 `contrast`를 바꿔도 그대로인 omnibus p야. 특정 비교의 p는 Wald로 구해.
- "A에서 유의, B에서 비유의"는 차이의 증거가 아니야. 두 효과의 차이(interaction 계수)를 직접 검정해.
- 세 번째 그룹을 넣고 빼면 size factor, dispersion, SE, baseMean, filtering이 바뀌어. 두 그룹만 분석하려면 raw count로 새 object를 만들어.

p 값을 BH로 보정하고 independent filtering을 거치는 과정은 [06](06_multiple_testing.md)으로, LFC shrinkage와 QC는 [07](07_lfc_shrinkage_and_qc.md)로 이어져. 한 유전자를 count부터 padj까지 손으로 따라가는 흐름은 [08](08_one_gene_end_to_end.md)에 모아 두었어.

## 연습문제

교재 부록 A의 문제를 그대로 옮겼어.

> **문제 13.** log2FC=1.0, SE=0.25 이다. 기본 Wald Z 와 nominal 95% CI 를 계산하라. CI 가 0 을 제외하면 BH padj<0.05 도 반드시 성립하는가?

<details>
<summary>풀이</summary>

```r
cat(sprintf("Z = %.2f, p = %.3g, 95%% CI = [%.2f, %.2f]\n",
            1/0.25, 2*pnorm(-4), 1 - qnorm(0.975)*0.25, 1 + qnorm(0.975)*0.25))
```

```
Z = 4.00, p = 6.33e-05, 95% CI = [0.51, 1.49]
```

$Z=1.0/0.25=4.00$, $p=2\Phi(-4)=6.33\times10^{-5}$이고, CI $=1\pm1.95996\times0.25=[0.51,\,1.49]$야. 교재 정답 [0.51, 1.49]와 같아.

CI가 0을 제외한다고 BH padj<0.05가 보장되지는 않아. CI가 0을 제외하는 건 raw p<0.05와 같은 사건일 뿐이야. padj는 전체 검정 수와 다른 유전자들의 p 분포에 따라 달라지기 때문에, 이 유전자의 raw p가 $6\times10^{-5}$로 작아도 0.05를 넘을 수 있어.

관계는 한 방향이야. BH padj는 언제나 raw p 이상이라서 padj<0.05이면 CI는 반드시 0을 제외하지만, 역은 성립하지 않아. 3절의 Glucose–Ctrl 비교에서 CI가 0을 제외한 277개 중 114개는 padj<0.05가 아니었어(103개는 padj ≥ 0.05, 11개는 independent filtering으로 NA). padj<0.05인데 CI가 0을 포함한 유전자는 하나도 없었고.

</details>

> **문제 14.** paired 세 그룹에서 full=`~pair+condition`, reduced=`~pair`인 LRT 의 기본 자유도 차이와 귀무가설을 말하라.

<details>
<summary>풀이</summary>

design matrix의 열 수는 full이 6(절편, pair 2·3·4, Starvation, Glucose), reduced가 4야. 그래서 df = 6 − 4 = 2야(4절 출력 `df = 6 - 4 = 2`). DESeq2 소스도 `df <- ncol(fullModelMatrix) - ncol(reducedModelMatrix)`로 계산해.

귀무가설은 pair 효과를 통제한 상태에서 $b_S=b_G=0$이야. reference coding 아래에서는 세 condition의 보정된 평균이 모두 같다는 뜻이지. 교재 정답과 같아.

그룹이 더 있거나 covariate가 추가되면 "세 그룹이니 2"가 아니라 두 design matrix의 rank 차이를 따져야 해. DESeq2는 열 수 차이를 계산하는데, full rank 검사를 통과한 nested design에서는 열 수 차이와 rank 차이가 같아(더 깊이 보기의 LRT 소스).

</details>

> **문제 15.** LRT object 의 결과표에서 Glucose−Starvation LFC 를 표시했다. 그 행의 p 값을 Glucose−Starvation pairwise p 값으로 읽어도 되는가?

<details>
<summary>풀이</summary>

안 돼. `results(dds_lrt, contrast = ...)`는 `log2FoldChange`와 `lfcSE`만 요청한 비교로 바꾸고, `stat`과 `pvalue`에는 `getStat(object, test, name = NULL)`로 omnibus 값을 그대로 붙여(더 깊이 보기의 `results()` 소스).

더 깊이 보기의 'Wald와 LRT 전체 실행 기록'에서 문제가 묻는 Glucose–Starvation 추출을 그대로 해 봤어. `pvalue`는 omnibus와 `identical` TRUE였고, 표시된 LFC만 Wald object의 Glucose–Starvation LFC와 같았어. Starvation–Ctrl, Glucose–Ctrl 추출도 마찬가지였어(4절).

gene841을 예로 들면 표시된 LFC가 무엇이든 p는 $9.5\times10^{-41}$이야. 이 p는 "Starvation이든 Glucose든 어느 쪽이 Ctrl과 다르다"는 증거야.

pairwise p가 필요하면 `results(dds_lrt, contrast = c("condition", "Glucose", "Starvation"), test = "Wald")`로 Wald p를 다시 만들거나 Wald object를 따로 두면 돼. 이 호출의 p는 Wald object의 Glucose–Starvation p와 같았고, 설명도 "Wald test p-value: condition Glucose vs Starvation"으로 바뀌었어. 표에는 어느 검정의 p인지 적어 둬. 교재 정답과 같아.

</details>

## 더 깊이 보기

아래 코드는 본문 코드에 이어서 같은 R 세션에서 실행했어.

<details>
<summary>DESeq2 소스로 확인한 것: Wald 검정 (<code>nbinomWaldTest</code>)</summary>

계수 하나가 아니라 계수 조합(contrast)을 검정할 때의 식이야. $c$는 contrast 벡터이고 단위는 log2야.

$$Z=\frac{\hat\delta}{SE(\hat\delta)},\qquad \hat\delta=c^\top\hat\beta,\qquad SE(\hat\delta)=\sqrt{c^\top\,\widehat{\mathrm{Cov}}(\hat\beta)\,c},\qquad p=2\{1-\Phi(\lvert Z\rvert)\}=2\,\Phi(-\lvert Z\rvert).$$

nominal 95% CI는 $\hat\beta\pm z_{0.975}\,SE(\hat\beta)$이고, $z_{0.975}=1.95996$이야.

SE는 IRLS(계수를 찾는 반복 가중 최소제곱 계산)가 수렴한 자리의 공분산 $\widehat{\mathrm{Cov}}(\hat b)\approx(X^\top W X)^{-1}$에서 와. $X^\top WX$는 Fisher 정보, 즉 데이터가 계수에 대해 주는 정보의 양이야. 그 역행렬이 공분산이니까 정보가 많을수록 SE가 작아. 단위는 자연로그 계수 $b$이고, $W$는 [04](04_glm_condition_batch.md)에 있어. log2 계수로는 $SE(\hat\beta)=SE(\hat b)/\ln 2$야(소스 `betaSE <- log2(exp(1)) * sqrt(sigma)`). `betaPrior=FALSE`여도 `fitNbinomGLMs`는 `lambda <- rep(1e-06, ...)`의 아주 작은 ridge를 넣어.

$W$ 안의 dispersion $\alpha_i$와 size factor $s_j$는 이미 추정된 값을 고정한 거야(`fitNbinomGLMs` 안의 `alpha_hat <- dispersions(object)`). 그래서 nominal SE에는 upstream 추정의 불확실성이 담기지 않아.

`print(DESeq2::nbinomWaldTest)`에서 해당 부분만 옮기면 이래.

```r
# 발췌 (실행하지 않음)
WaldStatistic <- betaMatrix/betaSE
...
if (useT) {
    if (!missing(df)) {
        ...                                      # 사용자가 준 df (길이 1 이면 유전자 수만큼 복제)
        df <- rep(df, nrow(objectNZ))
    } else {
        ...                                      # weights 가 있으면 num.samps <- rowSums(weights)
        num.samps <- rep(ncol(object), nrow(objectNZ))
        df <- num.samps - ncol(dispModelMatrix)
    }
    df <- ifelse(df > 0, df, NA)
    WaldPvalue <- 2 * pt(abs(WaldStatistic), df = df, lower.tail = FALSE)
} else {
    WaldPvalue <- 2 * pnorm(abs(WaldStatistic), lower.tail = FALSE)
}
```

`args(nbinomWaldTest)`를 보면 `betaPrior = FALSE`, `useT = FALSE`, `useQR = TRUE`, `minmu = 0.5`가 기본값이고, 기본값 없는 `df` 인자도 있어. `2*pnorm(abs(Z), lower.tail=FALSE)`는 $2\Phi(-\lvert Z\rvert)$와 같아.

`useT = TRUE`인데 `df`를 주지 않고 weights도 없으면 자유도는 $n-p$야. $n$은 sample 수, $p$는 `ncol(dispModelMatrix)`, 즉 계수 수야. 교재 4.5의 기호라서 p 값의 p와는 달라. DESeq2 소스는 sample 수를 `m`으로 써. `df`를 직접 주면 그 값을 쓰고.

</details>

<details>
<summary>t 분포로 바꾸면 얼마나 달라질까? (<code>useT = TRUE</code>)</summary>

```r
dds_t <- nbinomWaldTest(dds_wald, useT = TRUE, quiet = TRUE)
res_t <- results(dds_t, name = "condition_Glucose_vs_Ctrl")
cat("tDegreesFreedom (unique):", unique(mcols(dds_t)$tDegreesFreedom),
    " = m - ncol(dispModelMatrix) =", ncol(dds_t), "-", ncol(attr(dds_t, "dispModelMatrix")), "\n")
z <- res_wG["gene841", "stat"]
cat(sprintf("gene841: Z=%.3f  normal p=%.3g  t(df=6) p=%.3g  DESeq2 useT p=%.3g\n",
            z, 2 * pnorm(-abs(z)), 2 * pt(abs(z), df = 6, lower.tail = FALSE), res_t["gene841", "pvalue"]))
cat("padj<0.05: normal", sum(res_wG$padj < 0.05, na.rm = TRUE), " t:", sum(res_t$padj < 0.05, na.rm = TRUE), "\n")
```

```
tDegreesFreedom (unique): 6 NA  = m - ncol(dispModelMatrix) = 12 - 6
gene841: Z=11.246  normal p=2.43e-29  t(df=6) p=2.95e-05  DESeq2 useT p=2.95e-05
padj<0.05: normal 163  t: 35
```

`useT = TRUE`로 돌리면 `mcols(dds)$tDegreesFreedom` 열이 생기고(all-zero 행은 NA), `results()`는 이 열이 있는지를 보고 `useT`를 판단해(`useT <- "tDegreesFreedom" %in% names(mcols(object))`). 12 sample, 계수 6개에서 정규분포 대 t(6)의 차이는 유의 유전자 163 → 35개로 꽤 커. gene841의 p도 $2.4\times10^{-29}$에서 $3.0\times10^{-5}$가 되고. 교재 9.1은 "선택적인 t-reference 설정도 표준 기본값과 구분해야 한다"고만 적었는데, 왜 구분해야 하는지 이 차이를 보면 알 수 있어.

</details>

<details>
<summary>DESeq2 소스로 확인한 것: LRT (<code>nbinomLRT</code>)</summary>

자연로그 계수 $b$와 offset $\log s_j$로 쓴 기본 `type="DESeq2"` 경로의 식이야.

$$\ell(b)=\sum_j \log f_{NB}\!\left(K_{ij};\;\mu_{ij}=s_j\exp(x_j^\top b),\;\alpha_i\right),\qquad D=2\{\ell_{full}-\ell_{reduced}\}\;\overset{H_0}{\approx}\;\chi^2_{p_F-p_R}.$$

$K_{ij}$는 count, $\mu_{ij}$는 기대 count, $x_j$는 design matrix의 $j$번째 행이야. $p_F$, $p_R$은 교재 표기로 두 모형의 모수 수(열 수)이고, 이것도 p 값과는 다른 기호야.

full: $\log\mu_j=\log s_j+b_0+b_{pair(j)}+b_S I(S_j)+b_G I(G_j)$, reduced: $\log\mu_j=\log s_j+b_0+b_{pair(j)}$. $H_0: b_S=b_G=0$이고 $p_F-p_R=6-4=2$야. 두 모형 모두 같은 $\alpha_i$를 써.

```r
# 발췌 (실행하지 않음)
fullModelMatrix    <- stats::model.matrix.default(full,    data = as.data.frame(colData(object)))
reducedModelMatrix <- stats::model.matrix.default(reduced, data = as.data.frame(colData(object)))
df <- ncol(fullModelMatrix) - ncol(reducedModelMatrix)
...
fullModel    <- fitNbinomGLMs(objectNZ, modelFormula = full, ...)
reducedModel <- fitNbinomGLMs(objectNZ, modelFormula = reduced, ...)
LRTStatistic <- (2 * (fullModel$logLike - reducedModel$logLike))
LRTPvalue    <- pchisq(LRTStatistic, df = df, lower.tail = FALSE)
deviance     <- -2 * fullModel$logLike
```

- `nbinomLRT` 안에는 `estimateDispersions` 호출이 없어. 두 `fitNbinomGLMs` 호출 모두 내부에서 `alpha_hat <- dispersions(object)`로 object에 저장된 $\alpha_i$를 읽어. `logLike`는 `nbinomLogLike()`, 즉 `rowSums(dnbinom(counts, mu = mu, size = 1/disp, log = TRUE))`야.
- `checkLRT`는 reduced의 변수가 full에 모두 있는지만 확인해.
- `DESeq()` 쪽에서는 `test="LRT"`인데 `reduced`가 없으면 stop, `betaPrior=TRUE`면 "test='LRT' does not support use of LFC shrinkage"로 stop, `test="Wald"`에 `reduced`를 주면 "'reduced' ignored when test='Wald'"로 stop 해.
- 자유도를 교재는 "rank 차이"라고 쓰고, 소스는 열 수 차이로 계산해. `DESeq()`가 full design을 `checkFullRank`로 검사하고(matrix로 주면 reduced도 검사), 보통의 nested formula에서는 reduced 열이 full 열의 일부야. 그래서 이 검사를 통과한 design에서는 열 수 차이와 rank 차이가 같아.
- dispersion은 두 모형이 공유해. Wald용 object에 `nbinomLRT()`를 바로 걸어도 `DESeq(test="LRT")`와 통계량이 소수점까지 같았어(차이 0, 아래 전체 실행 기록). 교재 10.3은 Wald object와 LRT object를 따로 두어 해석 실수를 줄이고, low-level 함수로 dispersion을 공유할 때는 무엇을 재사용하는지 적으라고 권해. 여기서 `nbinomLRT(dds_wald, ...)`가 재사용한 건 `dds_wald`의 size factor와 최종 $\alpha_i$야.

위 코드는 기본 `type = "DESeq2"` 분기야. `args(nbinomLRT)`의 `type = c("DESeq2", "glmGamPoi")` 중 `glmGamPoi` 분기는 `disp_trend <- mcols(objectNZ)$dispFit`으로 trend 값을 쓰고, `glmGamPoi::test_de`의 `LRTStatistic <- qlr$f_statistic`, 즉 quasi-likelihood F 통계량을 넣어. `DESeq()`는 `nbinomLRT(..., type = dispersionEstimator)`로 부르고, `fitType == "glmGamPoi"`일 때 이 분기로 가. glmGamPoi는 이 환경에 설치되어 있지 않아서 실행으로 확인하지는 못했어. 이 노트의 LRT 설명은 모두 기본 `type="DESeq2"` 경로 기준이야.

교재 10.5의 $\delta_{int}$ 식은 로그의 밑을 적지 않았어. log2로 읽어야 DESeq2의 `genotypeB.conditionT` 계수와 같은 수가 되고, 자연로그로 읽으면 그 $\ln 2$배야. LRT로 검정하려면 reduced에서 interaction 항만 빼고 main effect는 남겨(df = 1).

</details>

<details>
<summary>DESeq2 소스로 확인한 것: <code>results()</code>의 검정 선택과 LRT 결과표</summary>

```r
# 발췌 (실행하지 않음)
if (missing(test)) {
    test <- attr(object, "test")
} else if (test == "Wald" & attr(object, "test") == "LRT") {
    object <- makeWaldTest(object)          # 저장된 beta/SE 로 Z, p 를 새로 만든다
} else if (test == "LRT" & attr(object, "test") == "Wald") {
    stop("the LRT requires the user run nbinomLRT or DESeq(dds,test='LRT')")
}
...
if (!(lfcThreshold == 0 & altHypothesis == "greaterAbs")) {
    if (test == "LRT") {
        stop("tests of log fold change above or below a theshold must be Wald tests.")
    }
    ...
    if (altHypothesis == "greaterAbs") {
        if (!useT) {
            pfunc_lfc <- function(lfc_T, lfc, se) {
                pnorm(-abs(lfc) + lfc_T, sd = se) + pnorm(-abs(lfc) - lfc_T, sd = se)
            }
            newPvalue <- mapply(pfunc_lfc, T, LFC, SE)
        }
        ...
        newStat <- LFC/SE
    }
    ...
    else if (altHypothesis == "greaterAbs2014") {
        newStat <- sign(LFC) * pmax((abs(LFC) - T)/SE, 0)
        newPvalue <- pmin(1, 2 * pfunc((abs(LFC) - T)/SE))
    }
    else if (altHypothesis == "lessAbs") {
        newStatAbove <- pmax((T - LFC)/SE, 0)
        ...
        newStatBelow <- pmax((LFC + T)/SE, 0)
        ...
        newStat <- pmin(newStatAbove, newStatBelow)
        newPvalue <- pmax(pvalueAbove, pvalueBelow)
    }
    else if (altHypothesis == "greater") {
        newStat <- pmax((LFC - T)/SE, 0)
        newPvalue <- pfunc((LFC - T)/SE)
    }
    else if (altHypothesis == "less") {
        newStat <- pmin((LFC + T)/SE, 0)
        newPvalue <- pfunc((-T - LFC)/SE)
    }
    res$stat <- newStat
```

LRT object에서 `results()`가 표를 채우는 `cleanContrast` 코드는 이래. `stat`/`pvalue`는 `name=NULL`로, 즉 contrast와 상관없이 `LRTStatistic`/`LRTPvalue` 열을 가져와.

```r
# 발췌 (실행하지 않음)
if (test == "LRT") {
    stat   <- getStat(object, test, name = NULL)
    pvalue <- getPvalue(object, test, name = NULL)
    res <- cbind(res[c("baseMean", "log2FoldChange", "lfcSE")], stat, pvalue)
```

- `baseMean`은 `getBaseMeansAndVariances`의 `rowMeans(counts(object, normalized=TRUE))`이고, description도 "mean of normalized counts for all samples"야.
- Cook's cutoff 기본값은 `qf(0.99, p, m - p)`야. 소스의 `p`, `m`은 `dispModelMatrix`의 열·행 수(본문의 계수 수와 sample 수)야.
- padj의 기본 보정법은 `formals(results)$pAdjustMethod`로 보면 "BH"야(아래 전체 실행 기록).

</details>

<details>
<summary><code>lfcThreshold</code>의 p 공식과 <code>altHypothesis</code> 여섯 가지</summary>

1.50.2의 `altHypothesis`에는 `"greaterAbs"`, `"greaterAbsUPSHOT"`, `"lessAbs"`, `"greater"`, `"less"`, `"greaterAbs2014"` 여섯 가지가 있어. 기본값 `"greaterAbs"`의 p 공식은 2014 논문 방식과 다르고, 논문 방식은 `greaterAbs2014`로 따로 남아 있어. $T$ = `lfcThreshold`일 때 DESeq2 1.50.2 소스에서 읽은 공식은 이래.

$$p_{\text{greaterAbs}}=\Phi\!\left(\frac{T-\lvert\hat\beta\rvert}{SE}\right)+\Phi\!\left(\frac{-T-\lvert\hat\beta\rvert}{SE}\right),\qquad p_{\text{greaterAbs2014}}=\min\!\left\{1,\;2\,\Phi\!\left(-\frac{\lvert\hat\beta\rvert-T}{SE}\right)\right\}.$$

앞의 것은 경계 $\beta=T$에서 $\hat\beta\sim N(T,SE^2)$라고 놓고, $\lvert\hat\beta\rvert$가 관측값 이상일 확률을 양쪽 꼬리 그대로 더한 거야. 뒤의 것은 한쪽 꼬리를 두 배 한 거고.

$-T-\lvert\hat\beta\rvert\le T-\lvert\hat\beta\rvert$이므로 둘째 항은 첫째 항보다 크지 않아. 따라서 $p_{\text{greaterAbs}}\le2\Phi((T-\lvert\hat\beta\rvert)/SE)$야. 또 두 항의 합은 $\lvert\hat\beta\rvert$가 커질수록 줄고 $\hat\beta=0$에서 1이라서 1을 넘지도 않아. 그래서 $p_{\text{greaterAbs}}\le p_{\text{greaterAbs2014}}$가 모든 경우에 성립해. 현재 기본값은 p 기준으로 2014 방식보다 보수적일 수 없다는 뜻이야. 아래 실행에서도 모든 유전자에서 그랬어.

`stat` 열도 `altHypothesis`에 따라 달라. `greaterAbs`(와 `greaterAbsUPSHOT`)만 $\hat\beta/SE$를 그대로 두고 `pvalue`만 바꿔. `greaterAbs2014`는 `stat`을 $\mathrm{sign}(\hat\beta)\max\{(\lvert\hat\beta\rvert-T)/SE,\,0\}$로 바꾸고, `greater`는 $\max\{(\hat\beta-T)/SE,0\}$, `less`는 $\min\{(\hat\beta+T)/SE,0\}$, `lessAbs`는 $\min[\max\{(T-\hat\beta)/SE,0\},\max\{(\hat\beta+T)/SE,0\}]$로 바꿔.

```r
eq  <- function(a, b) isTRUE(all.equal(a, b))
r0  <- results(dds_wald, name = "condition_Glucose_vs_Ctrl", alpha = 0.05)
r1  <- results(dds_wald, name = "condition_Glucose_vs_Ctrl", lfcThreshold = 1, altHypothesis = "greaterAbs",     alpha = 0.05)
r14 <- results(dds_wald, name = "condition_Glucose_vs_Ctrl", lfcThreshold = 1, altHypothesis = "greaterAbs2014", alpha = 0.05)
lfc <- r0$log2FoldChange; se <- r0$lfcSE
p_hand <- pnorm(-abs(lfc) + 1, sd = se) + pnorm(-abs(lfc) - 1, sd = se)
p_2014 <- pmin(1, 2 * pnorm((abs(lfc) - 1)/se, lower.tail = FALSE))
cat("naive filter padj<0.05 & |LFC|>1:", sum(r0$padj < 0.05 & abs(lfc) > 1, na.rm = TRUE),
    "| lfcThreshold=1 greaterAbs padj<0.05:", sum(r1$padj < 0.05, na.rm = TRUE), "| greaterAbs2014:", sum(r14$padj < 0.05, na.rm = TRUE), "\n")
cat("greaterAbs:     stat == r0$stat (LFC/SE):", eq(r1$stat, r0$stat), "\n")
cat("greaterAbs2014: stat == r0$stat:", eq(r14$stat, r0$stat),
    " | stat == sign(LFC)*pmax((|LFC|-1)/SE, 0):", eq(r14$stat, sign(lfc) * pmax((abs(lfc) - 1)/se, 0)), "\n")
cat("greaterAbs p == Phi((T-|b|)/SE)+Phi((-T-|b|)/SE):", eq(r1$pvalue, p_hand),
    " | greaterAbs2014 p == 2*Phi(-(|b|-T)/SE):", eq(r14$pvalue, p_2014), "\n")
cat("stat == LFC/SE for lessAbs / greater / less:",
    sapply(c("lessAbs", "greater", "less"), function(h)
      eq(results(dds_wald, name = "condition_Glucose_vs_Ctrl", lfcThreshold = 1, altHypothesis = h)$stat, r0$stat)), "\n")
cat("p_greaterAbs <= p_greaterAbs2014 for every gene:", all(r1$pvalue <= r14$pvalue, na.rm = TRUE), "\n")
gg <- c("gene8", "gene20", "gene38")
print(data.frame(LFC = r0[gg, "log2FoldChange"], SE = r0[gg, "lfcSE"], stat0 = r0[gg, "stat"], p0 = r0[gg, "pvalue"],
                 p_T1 = r1[gg, "pvalue"], padj_T1 = r1[gg, "padj"], stat_2014 = r14[gg, "stat"], p_2014 = r14[gg, "pvalue"],
                 row.names = gg), digits = 3)
cat("metadata(r1)$lfcThreshold:", metadata(r1)$lfcThreshold, "\n")
cat("LRT + lfcThreshold ->", tryCatch(results(dds_lrt, lfcThreshold = 1), error = function(e) conditionMessage(e)), "\n")
```

```
naive filter padj<0.05 & |LFC|>1: 156 | lfcThreshold=1 greaterAbs padj<0.05: 18 | greaterAbs2014: 12
greaterAbs:     stat == r0$stat (LFC/SE): TRUE
greaterAbs2014: stat == r0$stat: FALSE  | stat == sign(LFC)*pmax((|LFC|-1)/SE, 0): TRUE
greaterAbs p == Phi((T-|b|)/SE)+Phi((-T-|b|)/SE): TRUE  | greaterAbs2014 p == 2*Phi(-(|b|-T)/SE): TRUE
stat == LFC/SE for lessAbs / greater / less: FALSE FALSE FALSE
p_greaterAbs <= p_greaterAbs2014 for every gene: TRUE
         LFC    SE stat0       p0     p_T1 padj_T1 stat_2014  p_2014
gene8  -1.84 0.601 -3.06 2.18e-03 0.080679  0.6182     -1.40 0.16136
gene20  1.83 0.389  4.71 2.44e-06 0.016130  0.2403      2.14 0.03226
gene38  2.99 0.610  4.90 9.53e-07 0.000554  0.0327      3.26 0.00111
metadata(r1)$lfcThreshold: 1
LRT + lfcThreshold -> tests of log fold change above or below a theshold must be Wald tests.
```

`stat`은 `greaterAbs`에서 LFC/SE 그대로야(gene8: −3.06). `greaterAbs2014`에서는 $\mathrm{sign}(\hat\beta)\max\{(\lvert\hat\beta\rvert-1)/SE,0\}=-(1.842-1)/0.601=-1.40$으로 바뀌고, `lessAbs`/`greater`/`less`도 `stat`을 바꿔(출력의 FALSE 세 개). 같은 T = 1에서 `greaterAbs`는 18개, `greaterAbs2014`는 12개야. LRT object에는 `lfcThreshold`를 아예 쓸 수 없어.

</details>

<details>
<summary>CI 기준과 raw p 기준이 정확히 같아지는 조건</summary>

CI가 0을 제외하는 것과 raw p<0.05가 같다는 말은 같은 임계값을 쓸 때만 정확해. 3절 코드는 `qnorm(0.975)` $=1.959964$를 썼어. 1.96으로 반올림하면 $1.959964<\lvert Z\rvert<1.96$인 얇은 띠에서 두 기준이 어긋나.

</details>

<details>
<summary>Wald와 LRT 전체 실행 기록 (3 그룹 연습 데이터)</summary>

먼저 연습 데이터를 만든 방식이야. `makeExampleDESeqDataSet(n=1000, m=16, betaSD=1)`은 condition A(1–8), B(9–16) 두 그룹만 만들어. A의 앞 4개를 Ctrl, 뒤 4개를 Starvation, B의 앞 4개를 Glucose로 이름 붙이고, 각 그룹 안에서 pair 1–4를 매겼어. 그래서 Starvation–Ctrl은 모든 유전자에서 참 효과가 0이고, Glucose–Ctrl의 참 효과는 `mcols(dds0)$trueBeta`야.

```r
eq <- function(a, b) isTRUE(all.equal(a, b))
cat("resultsNames(dds_wald):", resultsNames(dds_wald), "\n")
cat("attr test:", attr(dds_wald, "test"), "/", attr(dds_lrt, "test"), " betaPrior:", attr(dds_wald, "betaPrior"), "\n")
cat("identical dispersions(wald, lrt):", identical(dispersions(dds_wald), dispersions(dds_lrt)), "\n")
fmm <- attr(dds_lrt, "modelMatrix"); rmm <- attr(dds_lrt, "reducedModelMatrix")   # nbinomLRT 가 object 에 저장한 두 행렬
cat("df = ncol(full) - ncol(reduced) =", ncol(fmm), "-", ncol(rmm), "=", ncol(fmm) - ncol(rmm), "\n")
cat("colnames(attr(dds_lrt, 'reducedModelMatrix')):", colnames(rmm), "\n")
cat("colnames(model.matrix(~ pair, meta))        :", colnames(model.matrix(~ pair, meta)), "\n")
cat("LRT-only mcols:", setdiff(names(mcols(dds_lrt)), names(mcols(dds_wald))), "\n")
lrt2 <- nbinomLRT(dds_wald, reduced = ~ pair, quiet = TRUE)   # Wald object 에 바로 LRT: 저장된 dispersion 을 재사용한다
cat("max |LRTStatistic diff| (nbinomLRT on wald obj vs DESeq LRT):",
    max(abs(mcols(lrt2)$LRTStatistic - mcols(dds_lrt)$LRTStatistic), na.rm = TRUE), "\n")
cat("Wald p == 2*pnorm(-|stat|):", eq(res_wS$pvalue, 2 * pnorm(-abs(res_wS$stat))), "\n")
cat("LRT  p == pchisq(stat, 2, lower=FALSE):", eq(res_l$pvalue, pchisq(res_l$stat, 2, lower.tail = FALSE)), "\n")
cat("res_l$stat == mcols(dds_lrt)$LRTStatistic:", eq(res_l$stat, mcols(dds_lrt)$LRTStatistic), "\n")
cat("LRT: pvalue identical for omnibus / Starvation-Ctrl / Glucose-Ctrl extraction:",
    identical(res_l$pvalue, res_lS$pvalue), identical(res_l$pvalue, res_lG$pvalue), "\n")
cat("LRT: log2FoldChange shown by default (last coef):", mcols(res_l)$description[2], "\n")
cat("LRT: Starvation-Ctrl extraction LFC == Wald object's Starvation-Ctrl LFC:", eq(res_lS$log2FoldChange, res_wS$log2FoldChange),
    " lfcSE equal:", eq(res_lS$lfcSE, res_wS$lfcSE), "\n")
cat("results(dds_lrt, test='Wald') p == Wald object p:", eq(res_lS_wald$pvalue, res_wS$pvalue), "\n")
cat("mcols(res_wS)$description:\n"); print(mcols(res_wS)$description)
cat("mcols(res_l)$description:\n");  print(mcols(res_l)$description)
res_lGS      <- results(dds_lrt,  contrast = c("condition", "Glucose", "Starvation"))                  # 문제 15 의 경우
res_lGS_wald <- results(dds_lrt,  contrast = c("condition", "Glucose", "Starvation"), test = "Wald")
res_wGS      <- results(dds_wald, contrast = c("condition", "Glucose", "Starvation"))
cat("LRT Glucose-Starvation extraction: pvalue identical to omnibus:", identical(res_l$pvalue, res_lGS$pvalue),
    "| LFC == Wald object's Glucose-Starvation LFC:", eq(res_lGS$log2FoldChange, res_wGS$log2FoldChange), "\n")
cat("results(dds_lrt, Glucose-Starvation, test='Wald') p == Wald object p:", eq(res_lGS_wald$pvalue, res_wGS$pvalue),
    "|", mcols(res_lGS_wald)$description[5], "\n")
cat("formals(results)$pAdjustMethod:", formals(results)$pAdjustMethod, "\n")
sh <- lfcShrink(dds_wald, coef = "condition_Glucose_vs_Ctrl", type = "normal", quiet = TRUE)
cat("lfcShrink(type='normal') description[2:3]:", mcols(sh)$description[2:3], sep = "\n  ")
cat("apeglm / ashr installed:", requireNamespace("apeglm", quietly = TRUE), requireNamespace("ashr", quietly = TRUE), "\n")
```

```
resultsNames(dds_wald): Intercept pair_2_vs_1 pair_3_vs_1 pair_4_vs_1 condition_Starvation_vs_Ctrl condition_Glucose_vs_Ctrl
attr test: Wald / LRT  betaPrior: FALSE
identical dispersions(wald, lrt): TRUE
df = ncol(full) - ncol(reduced) = 6 - 4 = 2
colnames(attr(dds_lrt, 'reducedModelMatrix')): Intercept pair_2_vs_1 pair_3_vs_1 pair_4_vs_1
colnames(model.matrix(~ pair, meta))        : (Intercept) pair2 pair3 pair4
LRT-only mcols: LRTStatistic LRTPvalue fullBetaConv reducedBetaConv
max |LRTStatistic diff| (nbinomLRT on wald obj vs DESeq LRT): 0
Wald p == 2*pnorm(-|stat|): TRUE
LRT  p == pchisq(stat, 2, lower=FALSE): TRUE
res_l$stat == mcols(dds_lrt)$LRTStatistic: TRUE
LRT: pvalue identical for omnibus / Starvation-Ctrl / Glucose-Ctrl extraction: TRUE TRUE
LRT: log2FoldChange shown by default (last coef): log2 fold change (MLE): condition Glucose vs Ctrl
LRT: Starvation-Ctrl extraction LFC == Wald object's Starvation-Ctrl LFC: TRUE  lfcSE equal: TRUE
results(dds_lrt, test='Wald') p == Wald object p: TRUE
mcols(res_wS)$description:
[1] "mean of normalized counts for all samples"
[2] "log2 fold change (MLE): condition Starvation vs Ctrl"
[3] "standard error: condition Starvation vs Ctrl"
[4] "Wald statistic: condition Starvation vs Ctrl"
[5] "Wald test p-value: condition Starvation vs Ctrl"
[6] "BH adjusted p-values"
mcols(res_l)$description:
[1] "mean of normalized counts for all samples"         "log2 fold change (MLE): condition Glucose vs Ctrl"
[3] "standard error: condition Glucose vs Ctrl"         "LRT statistic: '~ pair + condition' vs '~ pair'"
[5] "LRT p-value: '~ pair + condition' vs '~ pair'"     "BH adjusted p-values"
LRT Glucose-Starvation extraction: pvalue identical to omnibus: TRUE | LFC == Wald object's Glucose-Starvation LFC: TRUE
results(dds_lrt, Glucose-Starvation, test='Wald') p == Wald object p: TRUE | Wald test p-value: condition Glucose vs Starvation
formals(results)$pAdjustMethod: BH
lfcShrink(type='normal') description[2:3]:
  log2 fold change (MAP): condition Glucose vs Ctrl
  standard error: condition Glucose vs Ctrl
apeglm / ashr installed: FALSE FALSE
```

`lfcSE` 열의 뜻도 shrinkage 방법에 따라 달라져. `lfcShrink`의 apeglm/ashr 분기는 값과 description을 "posterior SD"로 바꾸고, `type="normal"`은 description이 "standard error"로 남지만 값은 MAP의 sandwich SE야(위 출력, 값의 정체는 [07](07_lfc_shrinkage_and_qc.md)에서 확인). apeglm/ashr는 이 환경에 설치되어 있지 않아서 소스로만 확인했어.

이번엔 같은 유전자 세 개에서 세 결과표의 `stat`/`pvalue`를 나란히 놓고, contrast별 padj<0.05 개수를 세 볼게.

```r
g3 <- c("gene841", "gene420", "gene4")
wald_tab <- data.frame(baseMean = res_wG[g3, "baseMean"],
  LFC_Glucose = res_wG[g3, "log2FoldChange"], stat_Glucose = res_wG[g3, "stat"], p_Glucose = res_wG[g3, "pvalue"],
  LFC_Starvation = res_wS[g3, "log2FoldChange"], stat_Starvation = res_wS[g3, "stat"], p_Starvation = res_wS[g3, "pvalue"],
  trueBeta_Glucose = mcols(dds0)[g3, "trueBeta"], row.names = g3)
lrt_tab <- data.frame(LFC_shown = res_l[g3, "log2FoldChange"], lfcSE_shown = res_l[g3, "lfcSE"],
  stat_D = res_l[g3, "stat"], p_omnibus = res_l[g3, "pvalue"], row.names = g3)
cat("Wald object (Glucose = Glucose-Ctrl, Starvation = Starvation-Ctrl):\n"); print(wald_tab, digits = 4)
cat("LRT object (LFC_shown = 표시용 마지막 계수 Glucose-Ctrl):\n"); print(lrt_tab, digits = 4)
sig <- function(r) !is.na(r$padj) & r$padj < 0.05
cat("padj<0.05 counts: Wald Glucose-Ctrl", sum(sig(res_wG)), "| Wald Starvation-Ctrl", sum(sig(res_wS)), "| LRT omnibus", sum(sig(res_l)), "\n")
cat("Wald Starvation-Ctrl padj<0.05 gene:", rownames(res_wS)[sig(res_wS)], " allZero:", mcols(dds_wald)$allZero[sig(res_wS)],
    "| counts:", counts(dds_wald)["gene818", ], "\n")
cat("gene818: Wald Starvation-Ctrl p =", signif(res_wS["gene818", "pvalue"], 3), "| LRT D =", signif(res_l["gene818", "stat"], 3),
    " LRT p =", signif(res_l["gene818", "pvalue"], 3), "\n")
cat("Wald Glucose padj<0.05 but LRT padj>=0.05:", sum(sig(res_wG) & !sig(res_l)), "| LRT<0.05 but Wald Glucose>=0.05:", sum(sig(res_l) & !sig(res_wG)), "\n")
lo <- sig(res_l) & !sig(res_wG)
cat("filterThreshold (baseMean): Wald Glucose-Ctrl", metadata(res_wG)$filterThreshold, "| LRT", metadata(res_l)$filterThreshold, "\n")
cat("LRT-only genes: Wald Glucose-Ctrl padj median", signif(median(res_wG$padj[lo]), 3),
    "| Starvation-Ctrl raw p<0.05:", sum(res_wS$pvalue[lo] < 0.05), "of", sum(lo), "\n")
```

```
Wald object (Glucose = Glucose-Ctrl, Starvation = Starvation-Ctrl):
        baseMean LFC_Glucose stat_Glucose p_Glucose LFC_Starvation stat_Starvation p_Starvation
gene841    918.0      3.6113      11.2459 2.427e-29        -0.3490         -1.0710     0.284176
gene420    373.6      1.5965       5.2985 1.167e-07        -0.8044         -2.6239     0.008693
gene4      166.7      0.2589       0.7308 4.649e-01         0.1590          0.4484     0.653838
        trueBeta_Glucose
gene841           2.9192
gene420           1.2825
gene4             0.2107
LRT object (LFC_shown = 표시용 마지막 계수 Glucose-Ctrl):
        LFC_shown lfcSE_shown stat_D p_omnibus
gene841    3.6113      0.3211 184.30 9.540e-41
gene420    1.5965      0.3013  64.77 8.628e-15
gene4      0.2589      0.3543   0.53 7.672e-01
padj<0.05 counts: Wald Glucose-Ctrl 163 | Wald Starvation-Ctrl 1 | LRT omnibus 146
Wald Starvation-Ctrl padj<0.05 gene: gene818  allZero: FALSE | counts: 1 13 0 0 0 0 0 0 0 0 4 0
gene818: Wald Starvation-Ctrl p = 4.09e-14 | LRT D = 0.28  LRT p = 0.87
Wald Glucose padj<0.05 but LRT padj>=0.05: 37 | LRT<0.05 but Wald Glucose>=0.05: 20
filterThreshold (baseMean): Wald Glucose-Ctrl 5.299015 | LRT 4.067603
LRT-only genes: Wald Glucose-Ctrl padj median 0.105 | Starvation-Ctrl raw p<0.05: 5 of 20
```

위 출력을 표로 정리하면 이래(W = Wald, G = Glucose–Ctrl, S = Starvation–Ctrl).

| gene | baseMean | W LFC G | W stat G | W p G | W LFC S | W stat S | W p S | LRT LFC(표시) | LRT stat | LRT p | trueBeta |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gene841 | 918.0 | 3.611 | 11.25 | 2.43e-29 | -0.349 | -1.07 | 0.284 | 3.611 | 184.30 | 9.54e-41 | 2.919 |
| gene420 | 373.6 | 1.597 | 5.30 | 1.17e-07 | -0.804 | -2.62 | 0.00869 | 1.597 | 64.77 | 8.63e-15 | 1.283 |
| gene4 | 166.7 | 0.259 | 0.73 | 0.465 | 0.159 | 0.45 | 0.654 | 0.259 | 0.53 | 0.767 | 0.211 |

- LRT 표의 `log2FoldChange`는 마지막 계수(Glucose vs Ctrl)이고, `stat`은 3.611/0.321이 아니라 두 계수를 함께 지운 $D$야. gene841의 omnibus p는 거의 전적으로 Glucose 효과가 만든 건데, 표만 봐서는 알 수 없어.
- gene420은 Starvation–Ctrl의 참 효과가 0인데도 raw p = 0.0087이야. 귀무가설 아래에서도 1000개 중 9개쯤은 이만큼 작은 p가 나올 거라고 기대할 수 있는 값이야.
- BH 후 Starvation–Ctrl에서 padj<0.05는 gene818 하나뿐인데, Starvation count가 모두 0인 유전자야. Ctrl과 Glucose에는 count가 있어서 DESeq2의 `allZero`는 FALSE야. 같은 유전자의 LRT omnibus는 $D=0.28$, $p=0.87$이야.

omnibus LRT(146)와 Glucose–Ctrl Wald(163)가 서로 포함 관계가 아닌 이유(Wald만 37개, LRT만 20개)는 본문 4절에 적었어. 출력으로 다시 보면 LRT만 잡은 20개는 Wald Glucose padj 중앙값이 0.105로 경계 근처에 있고, 그중 5개는 Starvation–Ctrl raw $p<0.05$라서 Starvation 쪽 잡음이 df 2 통계량에 더해졌어. 여기에 서로 다른 p 분포에 적용한 BH와 서로 다른 `filterThreshold`(baseMean 5.30 대 4.07)가 겹친 결과야.

</details>

<details>
<summary>gene841의 LRT 통계량을 두 NB log-likelihood로 직접 계산하기</summary>

DESeq2가 저장한 $\alpha_i$와 $s_j$를 그대로 쓰고, 자연로그 계수 $b$를 `optim`으로 최대화했어.

```r
g <- "gene841"
K <- counts(dds_lrt)[g, ]; sf <- sizeFactors(dds_lrt); a <- mcols(dds_lrt)[g, "dispersion"]
XF <- model.matrix(~ pair + condition, meta); XR <- model.matrix(~ pair, meta)
nll <- function(b, X) { mu <- sf * exp(X %*% b); -sum(dnbinom(K, mu = mu, size = 1/a, log = TRUE)) }
fit <- function(X) { st <- c(log(mean(K/sf)), rep(0, ncol(X)-1))
  optim(st, nll, X = X, method = "BFGS", control = list(maxit = 1000, reltol = 1e-14)) }
fF <- fit(XF); fR <- fit(XR)
lF <- -fF$value; lR <- -fR$value; D <- 2 * (lF - lR)
cat("gene:", g, " alpha =", a, "\n"); print(K)
cat("l_full =", lF, " l_reduced =", lR, " D = 2(l_full - l_reduced) =", D, "\n")
cat("DESeq2 LRTStatistic:", mcols(dds_lrt)[g, "LRTStatistic"], " p(chisq,df=2):", pchisq(D, 2, lower.tail = FALSE),
    " DESeq2 LRTPvalue:", mcols(dds_lrt)[g, "LRTPvalue"], "\n")
cat("full-model beta (log2) optim vs DESeq2 Wald object:\n")
print(rbind(optim_log2 = fF$par / log(2), DESeq2 = coef(dds_wald)[g, ]), digits = 5)
cat("deviance check: mcols(dds_lrt)$deviance =", mcols(dds_lrt)[g, "deviance"], " -2*l_full =", -2 * lF, "\n")
df_ <- cbind(meta, K = K, lsf = log(sf))
mF <- glm(K ~ pair + condition + offset(lsf), family = MASS::negative.binomial(theta = 1/a), data = df_)
mR <- glm(K ~ pair + offset(lsf),             family = MASS::negative.binomial(theta = 1/a), data = df_)
cat("MASS glm: D =", 2 * (as.numeric(logLik(mF)) - as.numeric(logLik(mR))), "\n")
cat(sprintf("Wald by hand: Z = %.4f / %.4f = %.4f, p = %.3g ; DESeq2 stat = %.4f\n",
            res_wG[g, "log2FoldChange"], res_wG[g, "lfcSE"], res_wG[g, "log2FoldChange"] / res_wG[g, "lfcSE"],
            res_wG[g, "pvalue"], res_wG[g, "stat"]))
```

```
gene: gene841  alpha = 0.09638025
      Ctrl_1       Ctrl_2       Ctrl_3       Ctrl_4 Starvation_1 Starvation_2 Starvation_3 Starvation_4
         283          197          148          177          203          131          167          144
   Glucose_1    Glucose_2    Glucose_3    Glucose_4
        1379         2511         3770         2385
l_full = -73.25379  l_reduced = -165.4043  D = 2(l_full - l_reduced) = 184.3011
DESeq2 LRTStatistic: 184.3011  p(chisq,df=2): 9.53965e-41  DESeq2 LRTPvalue: 9.53965e-41
full-model beta (log2) optim vs DESeq2 Wald object:
           Intercept pair_2_vs_1 pair_3_vs_1 pair_4_vs_1 condition_Starvation_vs_Ctrl
optim_log2    7.7307    -0.19151    0.016532    -0.27648                     -0.34896
DESeq2        7.7307    -0.19149    0.016564    -0.27646                     -0.34896
           condition_Glucose_vs_Ctrl
optim_log2                    3.6113
DESeq2                        3.6113
deviance check: mcols(dds_lrt)$deviance = 146.5076  -2*l_full = 146.5076
MASS glm: D = 184.3011
Wald by hand: Z = 3.6113 / 0.3211 = 11.2459, p = 2.43e-29 ; DESeq2 stat = 11.2459
```

$D$가 소수점 셋째 자리까지 일치하고, `MASS::negative.binomial(theta=1/α)`로 적합한 `glm()`도 같은 값을 줘. `mcols(dds_lrt)$deviance`가 $-2\ell_{full}$이라는 것도 볼 수 있어. full 모형의 $\hat b/\log 2$는 Wald object의 `coef()`와 같아. 결국 Wald와 LRT는 같은 full 적합 위에서 서로 다른 통계량을 계산할 뿐이야.

</details>

<details>
<summary>2 그룹 대 3 그룹 전체 기록</summary>

```r
sub <- dds_wald[, dds_wald$condition != "Glucose"]           # 적합된 object 를 자른다
cat("levels after subset:", levels(sub$condition), "| sizeFactors kept:", !is.null(sizeFactors(sub)), "\n")
cat("DESeq(sub) ->", tryCatch(DESeq(sub, quiet = TRUE), error = function(e) conditionMessage(e)), "\n")
sub$condition <- droplevels(sub$condition)
msgs <- character()
sub_fit <- withCallingHandlers(DESeq(sub), message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
cat("after droplevels, DESeq(sub) messages:", trimws(msgs), sep = "\n  ")
# 교재 12.4: raw count 와 metadata 로 새 object
keep <- meta$condition %in% c("Ctrl", "Starvation")
meta2 <- droplevels(meta[keep, , drop = FALSE]); cts2 <- cts[, rownames(meta2), drop = FALSE]
dds2 <- DESeq(DESeqDataSetFromMatrix(cts2, meta2, ~ pair + condition), quiet = TRUE)
r3 <- results(dds_wald, contrast = c("condition", "Starvation", "Ctrl"), alpha = 0.05)
r2 <- results(dds2,     contrast = c("condition", "Starvation", "Ctrl"), alpha = 0.05)
rs <- results(sub_fit,  contrast = c("condition", "Starvation", "Ctrl"), alpha = 0.05)
cat("size factors (Ctrl/Starvation samples):\n")
print(rbind(three_group = sizeFactors(dds_wald)[1:8], two_group_fresh = sizeFactors(dds2), subset_kept = sizeFactors(sub_fit)), digits = 5)
trend <- function(d) unlist(attr(dispersionFunction(d), "coefficients"))
cat("dispersion trend coef (asymptDisp, extraPois): 3g", trend(dds_wald), "| 2g", trend(dds2), "\n")
cat("prior var (log dispersion): 3g", attr(dispersionFunction(dds_wald), "dispPriorVar"), "| 2g", attr(dispersionFunction(dds2), "dispPriorVar"), "\n")
p3 <- ncol(attr(dds_wald, "dispModelMatrix")); m3 <- ncol(dds_wald); p2 <- ncol(attr(dds2, "dispModelMatrix")); m2 <- ncol(dds2)
cat("residual df m - p: 3g", m3, "-", p3, "=", m3 - p3, "| 2g", m2, "-", p2, "=", m2 - p2, "\n")
cat("median dispersion 3g / 2g:", median(dispersions(dds_wald), na.rm = TRUE), "/", median(dispersions(dds2), na.rm = TRUE), "\n")
cat("median lfcSE 3g / 2g:", median(r3$lfcSE, na.rm = TRUE), "/", median(r2$lfcSE, na.rm = TRUE),
    "| cor(LFC):", cor(r3$log2FoldChange, r2$log2FoldChange, use = "complete"), "\n")
cat("baseMean equal?", eq(r3$baseMean, r2$baseMean), " (3g uses all 12 samples)\n")
cat("filterThreshold 3g / 2g:", metadata(r3)$filterThreshold, "/", metadata(r2)$filterThreshold, "| padj NA:", sum(is.na(r3$padj)), "/", sum(is.na(r2$padj)), "\n")
cat("Cook's cutoff qf(.99,p,m-p): 3g", qf(.99, p3, m3 - p3), " 2g", qf(.99, p2, m2 - p2),
    "| maxCooks all NA (3g, 2g):", all(is.na(mcols(dds_wald)$maxCooks)), all(is.na(mcols(dds2)$maxCooks)), "\n")
cat("padj<0.05 Starvation-Ctrl: 3g", sum(r3$padj < 0.05, na.rm = TRUE), " 2g fresh", sum(r2$padj < 0.05, na.rm = TRUE), " subset-kept-sf", sum(rs$padj < 0.05, na.rm = TRUE), "\n")
cat("dispGeneEst identical 2g fresh vs subset-kept-sf?", identical(mcols(dds2)$dispGeneEst, mcols(sub_fit)$dispGeneEst),
    "| lfcSE max abs diff:", max(abs(rs$lfcSE - r2$lfcSE), na.rm = TRUE), "\n")
g4 <- c("gene841", "gene420", "gene4", "gene818")
print(data.frame(LFC_3g = r3[g4, "log2FoldChange"], SE_3g = r3[g4, "lfcSE"], p_3g = r3[g4, "pvalue"],
                 LFC_2g = r2[g4, "log2FoldChange"], SE_2g = r2[g4, "lfcSE"], p_2g = r2[g4, "pvalue"],
                 alpha_3g = mcols(dds_wald)[g4, "dispersion"], alpha_2g = mcols(dds2)[g4, "dispersion"],
                 baseMean_3g = r3[g4, "baseMean"], baseMean_2g = r2[g4, "baseMean"], row.names = g4), digits = 4)
cat("gene818 counts:", counts(dds_wald)["gene818", ], "| allZero:", mcols(dds_wald)["gene818", "allZero"], "\n")
cat("gene818 betaConv/betaIter: 3g", unlist(mcols(dds_wald)["gene818", c("betaConv", "betaIter")]),
    "| 2g", unlist(mcols(dds2)["gene818", c("betaConv", "betaIter")]), "\n")
keepG <- meta$condition %in% c("Ctrl", "Glucose"); metaG <- droplevels(meta[keepG, , drop = FALSE])
ddsG <- DESeq(DESeqDataSetFromMatrix(cts[, rownames(metaG)], metaG, ~ pair + condition), quiet = TRUE)
rG <- results(ddsG, contrast = c("condition", "Glucose", "Ctrl"), alpha = 0.05)
cat("Glucose-Ctrl: padj<0.05 3g", sum(res_wG$padj < 0.05, na.rm = TRUE), " 2g", sum(rG$padj < 0.05, na.rm = TRUE),
    "| median lfcSE 3g/2g", median(res_wG$lfcSE, na.rm = TRUE), "/", median(rG$lfcSE, na.rm = TRUE),
    "| cor LFC", cor(res_wG$log2FoldChange, rG$log2FoldChange, use = "complete"), "\n")
```

```
levels after subset: Ctrl Starvation Glucose | sizeFactors kept: TRUE
DESeq(sub) -> full model matrix is less than full rank
after droplevels, DESeq(sub) messages:
  using pre-existing size factors
  estimating dispersions
  found already estimated dispersions, replacing these
  gene-wise dispersion estimates
  mean-dispersion relationship
  final dispersion estimates
  fitting model and testing
size factors (Ctrl/Starvation samples):
                Ctrl_1 Ctrl_2  Ctrl_3 Ctrl_4 Starvation_1 Starvation_2 Starvation_3 Starvation_4
three_group     1.0367 1.0254 0.99114 1.0231       1.0088       1.0356       1.0534       1.0600
two_group_fresh 1.0425 1.0250 1.01683 1.0372       1.0415       1.0427       1.0341       1.0571
subset_kept     1.0367 1.0254 0.99114 1.0231       1.0088       1.0356       1.0534       1.0600
dispersion trend coef (asymptDisp, extraPois): 3g 0.08018223 5.795031 | 2g 0.07442075 7.33897
prior var (log dispersion): 3g 0.25 | 2g 0.25
residual df m - p: 3g 12 - 6 = 6 | 2g 8 - 5 = 3
median dispersion 3g / 2g: 0.4142668 / 0.5139128
median lfcSE 3g / 2g: 0.7118803 / 0.7851631 | cor(LFC): 0.8197001
baseMean equal? FALSE  (3g uses all 12 samples)
filterThreshold 3g / 2g: 0.07902854 / 0.1194288 | padj NA: 1 / 4
Cook's cutoff qf(.99,p,m-p): 3g 8.466125  2g 28.23708 | maxCooks all NA (3g, 2g): TRUE TRUE
padj<0.05 Starvation-Ctrl: 3g 1  2g fresh 0  subset-kept-sf 0
dispGeneEst identical 2g fresh vs subset-kept-sf? FALSE | lfcSE max abs diff: 0.03222073
          LFC_3g  SE_3g      p_3g  LFC_2g  SE_2g     p_2g alpha_3g alpha_2g baseMean_3g baseMean_2g
gene841  -0.3490 0.3258 2.842e-01 -0.3181 0.3140 0.310934  0.09638  0.08898      917.99     174.769
gene420  -0.8044 0.3066 8.693e-03 -0.8198 0.3130 0.008811  0.08429  0.08808      373.58     184.286
gene4     0.1590 0.3546 6.538e-01  0.1436 0.3332 0.666441  0.11467  0.10043      166.72     158.169
gene818 -21.6388 2.8630 4.088e-14 -2.0068 2.8101 0.475141  6.47240  6.16654        1.45       1.705
gene818 counts: 1 13 0 0 0 0 0 0 0 0 4 0 | allZero: FALSE
gene818 betaConv/betaIter: 3g 1 100 | 2g 1 22
Glucose-Ctrl: padj<0.05 3g 163  2g 152 | median lfcSE 3g/2g 0.7195894 / 0.7479812 | cor LFC 0.9790813
```

| gene | LFC 3g | SE 3g | p 3g | LFC 2g | SE 2g | p 2g | $\alpha$ 3g | $\alpha$ 2g | baseMean 3g | baseMean 2g |
|---|---|---|---|---|---|---|---|---|---|---|
| gene841 | -0.349 | 0.326 | 0.284 | -0.318 | 0.314 | 0.311 | 0.0964 | 0.0890 | 918.0 | 174.8 |
| gene420 | -0.804 | 0.307 | 0.00869 | -0.820 | 0.313 | 0.00881 | 0.0843 | 0.0881 | 373.6 | 184.3 |
| gene4 | 0.159 | 0.355 | 0.654 | 0.144 | 0.333 | 0.666 | 0.1150 | 0.1000 | 166.7 | 158.2 |
| gene818 | -21.639 | 2.863 | 4.09e-14 | -2.007 | 2.810 | 0.475 | 6.47 | 6.17 | 1.5 | 1.7 |

- 같은 8개 sample의 같은 raw count인데도 size factor부터 달라(기준 geometric mean이 12개 대 8개 sample). `sub_fit`은 3 그룹 size factor를 물려받았지만(`subset_kept` 행) dispersion은 새로 추정돼서 `dds2`와 gene-wise 값이 달라.
- 교재의 "LFC는 비슷하고 SE가 달라진다"는 gene841/420/4에서 그대로 보여. dispersion이 두 object에서 다르게 shrink되고(trend 계수 $(a_0,a_1)$이 다름) 잔여 자유도도 다르기 때문이야. prior variance는 두 경우 모두 하한 0.25에 걸려 같았어(`estimateDispersionsPriorVar` 소스의 `dispPriorVar <- pmax((varLogDispEsts - expVarLogDisp), 0.25)`, 값은 출력의 `attr(dispersionFunction(.), "dispPriorVar")`).
- 3 그룹 쪽 SE가 전반적으로 작아. Starvation–Ctrl median lfcSE가 0.712 대 0.785, Glucose–Ctrl도 0.720 대 0.748이야. 잔여 df가 6 대 3이고 median $\alpha$도 0.414 대 0.514야. Glucose는 같은 모의 분포에서 나온 비슷한 변동성의 그룹이라, 교재 12.2의 "비슷한 수준의 변동성을 가진 replicate가 더 들어오면 dispersion 추정에 정보를 더 준다"는 방향과 맞아. 반대로 제3 그룹의 이질성이 크면 Ctrl–Starvation 비교의 SE가 그 영향을 받아.
- gene818은 Starvation count가 모두 0인 유전자야(Ctrl 1, 13, 0, 0 / Starvation 0, 0, 0, 0 / Glucose 0, 0, 4, 0). Starvation 계수의 MLE는 $-\infty$이고, IRLS가 멈춘 자리에 따라 3 그룹에서는 LFC $-21.6$, $p=4\times10^{-14}$, 2 그룹에서는 $-2.0$, $p=0.48$이 나왔어(`betaIter` 100 대 22). 이런 유전자는 정규근사가 성립하지 않아서 Wald p에 의미가 없는데, 3 그룹 object에서 Starvation–Ctrl padj<0.05인 유일한 유전자가 바로 이거야. count가 이 정도로 작으면 `stat`이 아니라 count 자체를 봐야 해(같은 유전자의 LRT omnibus는 $p=0.87$).
- `baseMean`은 정의상 포함 sample 전체의 평균이라, 3 그룹 object에서는 Glucose의 큰 count가 섞여 918이 돼. independent filtering의 threshold도 그래서 달라져(0.079 vs 0.119).
- 이 design에서는 Cook's 필터가 어느 쪽에도 적용되지 않았어(`maxCooks` 1000개 모두 NA). `recordMaxCooks`가 `nOrMoreInCell(modelMatrix, n = 3)`으로 "같은 design 행을 3개 이상 가진 sample"만 대상으로 삼는데, pair×condition 조합마다 sample이 하나뿐이거든. `qf(0.99, p, m-p)` 값 자체는 8.47과 28.2로 다르지만 여기서는 쓰이지 않았어. Cook's distance의 정의와 replacement 조건은 [07](07_lfc_shrinkage_and_qc.md)에서 다뤄.
- Glucose–Ctrl도 padj<0.05 수가 3 그룹 163개 대 2 그룹 152개로 다르고, LFC 상관은 0.979야.
- vignette FAQ 원문은 설치된 `DESeq2/doc/DESeq2.Rmd` 2803행에 있어. "If I have multiple groups, should I run all together or split into pairs of groups?"라는 질문에 "Typically, we recommend users to run samples from all groups together"라며 함께 적합한 뒤 `contrast`로 비교를 뽑으라고 권해. 이어서 EDA(PCA 등)에서 한 그룹의 within-group 변동이 극단적으로 클 때는 두 그룹만 subset해서 `DESeq`를 돌리는 편이 더 민감할 수 있다고 설명해.

</details>

<details>
<summary>interaction 전체 기록</summary>

```r
set.seed(2); d <- makeExampleDESeqDataSet(n = 500, m = 16)        # betaSD 기본값 0: 모든 유전자의 참 효과 0
cat("betaSD default:", formals(makeExampleDESeqDataSet)$betaSD, "| all trueBeta == 0:", all(mcols(d)$trueBeta == 0), "\n")
d$genotype  <- factor(rep(c("A","B"), each = 8)); d$condition <- factor(rep(c("C","T"), 8))
design(d) <- ~ genotype + condition + genotype:condition
dw <- DESeq(d, quiet = TRUE)
dl <- DESeq(d, test = "LRT", reduced = ~ genotype + condition, quiet = TRUE)
rA  <- results(dw, name = "condition_T_vs_C")                                        # T-C in A
rB  <- results(dw, contrast = list(c("condition_T_vs_C", "genotypeB.conditionT")))    # T-C in B = dA + int
int <- results(dw, name = "genotypeB.conditionT")
rl  <- results(dl)
cat("resultsNames:", resultsNames(dw), "\n")
cat("interaction LRT description:", mcols(rl)$description[5], "\n")
cat("df =", ncol(attr(dl, "modelMatrix")) - ncol(attr(dl, "reducedModelMatrix")), "\n")
cat("interaction LFC == (T-C in B) - (T-C in A):", eq(int$log2FoldChange, rB$log2FoldChange - rA$log2FoldChange), "\n")
d2 <- d; design(d2) <- ~ genotype + genotype:condition                    # 같은 모형, 다른 parametrization
dw2 <- DESeq(d2, quiet = TRUE)
cat("alt resultsNames:", resultsNames(dw2), "\n")
dA2 <- results(dw2, name = "genotypeA.conditionT")$log2FoldChange; dB2 <- results(dw2, name = "genotypeB.conditionT")$log2FoldChange
cat("interaction LFC == dB - dA from ~genotype+genotype:condition: max|diff| =",
    max(abs(int$log2FoldChange - (dB2 - dA2)), na.rm = TRUE), "\n")
cat("SE(int) vs sqrt(SE_A^2 + SE_B^2): max|diff| =", max(abs(int$lfcSE - sqrt(rA$lfcSE^2 + rB$lfcSE^2)), na.rm = TRUE),
    " (A 와 B 는 서로 다른 sample -> Cov(dA, dB) = 0)\n")
cat("SE(B)   vs sqrt(SE_A^2 + SE_int^2): max|diff| =", max(abs(rB$lfcSE - sqrt(rA$lfcSE^2 + int$lfcSE^2)), na.rm = TRUE),
    " (dB = dA + int 는 상관됨)\n")
cat("Wald p(interaction) == LRT p (df=1)? max|diff| =", max(abs(int$pvalue - rl$pvalue), na.rm = TRUE), "\n")
one <- xor(rA$pvalue < 0.05, rB$pvalue < 0.05)
cat("genes where exactly one of A/B effects has p<0.05:", sum(one, na.rm = TRUE),
    "; among them interaction p<0.05:", sum(one & int$pvalue < 0.05, na.rm = TRUE), "\n")
cat("interaction p<0.05 among all genes:", sum(int$pvalue < 0.05, na.rm = TRUE), "of", sum(!is.na(int$pvalue)), "\n")
```

```
betaSD default: 0 | all trueBeta == 0: TRUE
resultsNames: Intercept genotype_B_vs_A condition_T_vs_C genotypeB.conditionT
interaction LRT description: LRT p-value: '~ genotype + condition + genotype:condition' vs '~ genotype + condition'
df = 1
interaction LFC == (T-C in B) - (T-C in A): TRUE
alt resultsNames: Intercept genotype_B_vs_A genotypeA.conditionT genotypeB.conditionT
interaction LFC == dB - dA from ~genotype+genotype:condition: max|diff| = 0.0002154644
SE(int) vs sqrt(SE_A^2 + SE_B^2): max|diff| = 0.0001432716  (A 와 B 는 서로 다른 sample -> Cov(dA, dB) = 0)
SE(B)   vs sqrt(SE_A^2 + SE_int^2): max|diff| = 3.184341  (dB = dA + int 는 상관됨)
Wald p(interaction) == LRT p (df=1)? max|diff| = 0.1378705
genes where exactly one of A/B effects has p<0.05: 60 ; among them interaction p<0.05: 21
interaction p<0.05 among all genes: 27 of 500
```

- 이 시뮬레이션은 `betaSD` 기본값이 0이라, 모든 유전자에서 genotype·condition·interaction의 참 효과가 0인 순수 귀무 자료야(출력 첫 줄).
- `interaction LFC == (T-C in B) - (T-C in A)`는 `rB`를 interaction 계수로 만들었으니 정의상 참이야. 그래서 같은 모형을 `~ genotype + genotype:condition`으로 다시 parametrize해 각 genotype 안의 T−C 효과를 직접 추정한 뒤 빼 봤더니, 최대 차이 $2.2\times10^{-4}$(수치 오차 수준)로 맞았어.
- A와 B sample은 겹치지 않아 $\mathrm{Cov}(\hat\delta_A,\hat\delta_B)=0$이므로, $\sqrt{SE_A^2+SE_B^2}$가 interaction SE와 거의 같아(최대 차이 $1.4\times10^{-4}$). 하지만 $\hat\delta_B=\hat\delta_A+\hat\delta_{int}$에서 $\hat\delta_A$와 $\hat\delta_{int}$는 genotype A의 같은 sample을 함께 써서 상관되어 있고, 이 둘로 $SE_B$를 같은 공식으로 구하면 최대 3.18 어긋나. 일반적으로 차이의 SE는 $\sqrt{c^\top\widehat{\mathrm{Cov}}\,c}$이고, interaction 계수로 parametrize하면 그게 곧 그 계수의 `lfcSE`야.
- df가 1이어도 Wald와 LRT의 p는 같지 않아(최대 차이 0.138). 두 근사가 다르기 때문이야.

</details>

<details>
<summary>흔히 하는 오해</summary>

- "DESeq2 Wald는 count가 정규분포라고 가정한다"? 아니야. 정규근사는 $\hat\beta$의 표본분포에 대한 것이고, count는 NB로 모형화해.
- "LFC가 크면 p도 작다"? 같은 LFC 1.0도 SE에 따라 $p=5.7\times10^{-7}$부터 0.32까지 가(2절).
- "`stat`은 항상 LFC/lfcSE다"? Wald에서만, 그것도 `altHypothesis`가 기본 `greaterAbs`(또는 `greaterAbsUPSHOT`)일 때만 그래. LRT object의 `stat`은 $D$이고(gene841: 184.3 대 11.25), `greaterAbs2014`·`greater`·`less`·`lessAbs`는 `stat`도 바꿔(gene8: −3.06 → −1.40).
- "LRT object에 contrast를 주면 그 비교의 p가 나온다"? `pvalue` 열은 contrast와 무관한 omnibus p야(`identical()`로 확인). pairwise p가 필요하면 `test="Wald"`를 주거나 Wald object를 따로 만들어.
- "CI가 0을 제외하면 유의하다"? raw $p<0.05$와 같을 뿐, padj<0.05를 보장하지는 않아(277 대 163). 거꾸로 padj<0.05이면 CI가 0을 제외하는 건 항상 성립해.
- "`padj<0.05 & abs(LFC)>1`이 곧 $\lvert\beta\rvert>1$ 검정이다"? 그 가설은 `lfcThreshold=1`로 검정해야 하고, 결과도 156 대 18로 달라. `greaterAbs`와 `greaterAbs2014`도 18 대 12로 달라.
- "LRT는 reduced 모형에 맞는 dispersion을 따로 추정한다"? 기본 `type="DESeq2"`의 `nbinomLRT`는 저장된 $\alpha_i$를 두 모형이 공유하고 GLM만 두 번 적합해(glmGamPoi 경로는 trend `dispFit`과 QL F).
- "3 그룹 object를 subset하면 2 그룹 분석과 같다"? level이 남아 rank 오류가 나고, `droplevels` 후에도 size factor는 재사용돼. 새 object를 만들어야 처음부터 다시 추정돼.
- "A에서 유의, B에서 비유의 → 효과가 다르다"? interaction 계수 $\delta_B-\delta_A$를 그 SE로 검정해야 해. 참 효과가 모두 0인 시뮬레이션에서 "한쪽만 유의"는 60개의 거짓 차이를 만들었고, interaction 검정으로는 그중 21개(전체 500개 중 27개, 명목 수준 근처)였어.

</details>

<details>
<summary>교재와 다르게 나온 부분</summary>

근거 열에 적은 이름은 이 노트의 본문 절 번호와 위 접이식 블록 제목이야. 표에 없는 나머지 교재 설명은 실제 DESeq2 동작과 맞았어.

| 교재 설명 | 확인 결과 | 근거 |
|---|---|---|
| 9.4 CI는 $\hat\beta\pm1.96SE$, 0 제외 = raw p<0.05, padj와 다름 | 일치. 다만 1.96은 $z_{0.975}=1.959964$의 반올림이라 정확한 동치는 $z_{0.975}$를 쓸 때만 성립 | 3절(`qnorm(0.975)` 사용): CI 제외 277 = raw p<0.05 277, padj<0.05 163. padj<0.05이면서 CI가 0을 포함하는 유전자 0개 (padj $\ge$ p) |
| 9.5 threshold test의 세부 구현은 버전에 따라 다를 수 있음 | 일치, 교재보다 구체적인 내용 있음 | 1.50.2의 `altHypothesis`는 6종이야. 기본 `greaterAbs`의 p는 $\Phi((T-\lvert\hat\beta\rvert)/SE)+\Phi((-T-\lvert\hat\beta\rvert)/SE)$이고, 논문 방식은 `greaterAbs2014`로 따로 있어(18 vs 12개). `stat`도 `greaterAbs` 외의 네 경우에는 LFC/SE가 아니야('lfcThreshold' 블록) |
| 10.1 NB-GLM에서는 LRT가 omnibus 역할, raw count ANOVA는 불필요 | 개념상 일치 | `nbinomLRT`가 full/reduced NB 적합의 logLike 차이를 써 |
| 10.2 $D=2(\ell_{full}-\ell_{reduced})\sim\chi^2_{p_F-p_R}$, 세 그룹이면 df 2 | 일치, 표현만 달라. 교재는 rank 차이, 소스는 `ncol` 차이인데 `checkFullRank`를 통과한 nested design에서는 같아 | 소스 `LRTStatistic <- 2*(fullModel$logLike - reducedModel$logLike)`, `pchisq(..., df = df)`, `df <- ncol(fullModelMatrix) - ncol(reducedModelMatrix)`; df $=6-4=2$; gene841 직접 계산 $D=184.301$ 일치. 기본 `type="DESeq2"` 기준이고, `glmGamPoi` 분기는 QL F(`qlr$f_statistic`)인데 미설치라 실행하지 못했어 |
| 10.4 omnibus를 먼저 할 필요 없음 | 미확인 (통계 원칙) | 소프트웨어 동작이 아니라 설계 원칙이야. Wald 163 / LRT 146이 서로 포함 관계가 아니라는 것은 확인했어 |
| 10.5 식 $\delta_{int}=(\log q_{B,T}-\log q_{B,C})-(\log q_{A,T}-\log q_{A,C})$ | 일치, 다만 로그의 밑이 없음 | 교재는 밑을 적지 않았어. log2로 읽어야 DESeq2 `genotypeB.conditionT`(log2 계수, 부록 C의 $\beta$)와 같고, 자연로그($b$)로 읽으면 $\ln 2$배야 |
| 12.2 3 그룹을 넣고 빼면 size factor, gene-wise dispersion, trend, prior, SE, baseMean, filtering이 달라짐 | 대체로 일치 | '2 그룹 대 3 그룹 전체 기록': size factor·trend 계수·median $\alpha$·median SE·baseMean·filterThreshold가 모두 달랐어. prior variance는 두 경우 모두 하한 0.25로 같았어. 다른 데이터에서는 달라질 수 있어 |
| 12.2 outlier 처리가 달라질 수 있음 | 조건부 일치 | `qf(0.99,p,m-p)`는 8.47 vs 28.2로 다르지만, 이 design에서는 `maxCooks`가 양쪽 모두 NA라 적용되지 않았어(기전은 [07](07_lfc_shrinkage_and_qc.md)) |
| 12.2 비슷한 변동성의 replicate가 더 들어오면 dispersion 추정에 정보를 더 줌 | 방향은 일치 | 잔여 df 6 vs 3, median $\alpha$ 0.414 vs 0.514, Starvation–Ctrl median lfcSE 0.712 vs 0.785 |
| 12.4 적합된 dds를 subset하면 기존 size factor가 남음 | 일치, 보충 있음 | `sizeFactors(sub)` 유지, `DESeq()` 메시지 "using pre-existing size factors". 다만 dispersion은 새로 추정돼(`dispGeneEst`가 fresh 2g와 다름) |

</details>

---

← 이전: [04. NB-GLM](04_glm_condition_batch.md) · 다음: [06. 다중검정](06_multiple_testing.md) →

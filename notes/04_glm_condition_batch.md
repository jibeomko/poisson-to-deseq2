# 04. 조건과 batch를 넣은 모형에서 LFC와 SE는 어떻게 나올까?

[03 노트](03_dispersion_estimation.md)에서 유전자 A의 dispersion α를 0.053147로 정했어. 그럼 Starvation이 이 유전자를 몇 배 바꿨는지, 그 추정은 얼마나 흔들리는지는 어떻게 구할까? α를 고정하고 count 평균의 log를 조건·donor·batch의 합으로 쓰는 모형(GLM)을 맞추면, 계수가 LFC가 되고 그 모형의 정보량에서 SE가 나와. 이 노트에서는 그 과정을 손으로 한 번, DESeq2로 한 번 따라가 봐.

> 교재 4장, 8장 · 코드: [04_glm_condition_batch.R](../04_glm_condition_batch.R)

## 1. Starvation은 유전자 A를 몇 배 바꿨을까?

교재 16장의 유전자 A를 다시 꺼내 볼게. 여기서 count는 한 sample에서 그 유전자에 배정된 read 수야.

| 조건 | sample 1 | sample 2 | sample 3 | 평균 |
|---|---|---|---|---|
| Ctrl | 100 | 130 | 90 | 106.667 |
| Starvation | 200 | 250 | 180 | 210 |

size factor는 sample마다 sequencing 깊이가 다른 것을 맞춰 주는 배율인데([02 노트](02_negative_binomial.md)), 여기서는 여섯 sample 모두 1이라 count를 그대로 비교해도 돼. 평균의 비는 210 / 106.667 = 1.969, 그러니까 Starvation에서 거의 2배가 됐어. 이 비에 log2를 취하면 0.977인데, 교재 16.4의 LFC 0.977280과 같은 값이야.

그런데 실제 실험은 이렇게 깔끔하지 않아. size factor가 sample마다 다르고, donor나 batch가 섞이고, 그룹이 셋 이상일 수도 있지. 그래서 DESeq2는 평균을 직접 나누는 대신, 어떤 경우에나 같은 방식으로 쓸 수 있는 모형을 맞춰. 그 모형이 GLM(일반화 선형모형)이야.

GLM은 count의 평균을 log scale에서 조건·batch 같은 항의 합으로 표현하는 모형이야. 평균을 그대로 쓰지 않고 log를 취해서 항들과 잇는다는 뜻에서 이걸 log link라고 불러. DESeq2의 GLM은 여기에 count가 음이항분포(NB)를 따른다는 가정을 더해. NB는 과산포, 곧 같은 조건 replicate 사이의 퍼짐이 Poisson이 예상하는 것보다 큰 현상을 허용하는 count 분포야([02 노트](02_negative_binomial.md)). 식으로 쓰면 이래.

$$K_{ij}\sim \mathrm{NB}(\mu_{ij},\alpha_i),\qquad \log\mu_{ij}=\log s_j+x_j^{T}b_i$$

- $K_{ij}$는 유전자 $i$, sample $j$의 count이고, $\mu_{ij}$는 모형이 그 sample에서 예상하는 평균 count(기대 count)야.
- $\alpha_i$는 유전자 $i$의 dispersion이야. 분산 $=\mu+\alpha\mu^2$에서 Poisson보다 더 퍼지는 정도로, [03 노트](03_dispersion_estimation.md)에서 구했어.
- $s_j$는 sample $j$의 size factor야. $\log s_j$는 계수 없이 그대로 더해지는데, 이런 항을 offset이라고 불러.
- $x_j$는 sample $j$가 어느 그룹에 속하는지 적은 숫자 한 줄, 곧 아래 design matrix의 한 행이야.
- $b_i$는 유전자 $i$의 계수(자연로그 단위)이고, $x_j^{T}b_i$는 $x_j$의 각 숫자에 대응하는 계수를 곱해서 더한 값이야. 위첨자 $T$(전치)는 세로 벡터를 가로로 눕힌다는 표시야.

log를 쓰면 평균이 항상 양수가 되고, "몇 배"라는 곱셈 차이가 log scale에서는 덧셈이 돼.

### sample 표를 숫자로 바꾼 design matrix

design matrix(설계행렬) X는 각 sample이 어떤 조건·batch·pair에 속하는지 숫자로 적은 표야. R에서는 `model.matrix()`가 design formula를 이 표로 바꿔 줘. 이 노트의 코드는 R 4.5.2와 DESeq2 1.50.2로 돌렸고, 뒤 코드가 앞에서 만든 객체를 이어 쓰니 따라 할 때는 위에서부터 차례로 실행하면 돼.

```r
k  <- c(100, 130, 90, 200, 250, 180)         # 유전자 A. size factor는 모두 1
cd <- data.frame(condition = factor(rep(c("Ctrl", "Starvation"), each = 3),
                                    levels = c("Ctrl", "Starvation")))
X <- model.matrix(~ condition, cd)
t(X[, ])                                     # 보기 쉽게 가로로 눕혔다: 열 = sample
b <- c(b0 = log(mean(k[1:3])), bS = log(mean(k[4:6]) / mean(k[1:3])))
b
b[["bS"]] / log(2)                           # beta_S = log2 fold change
```

```
                    1 2 3 4 5 6
(Intercept)         1 1 1 1 1 1
conditionStarvation 0 0 0 1 1 1
       b0        bS
4.6697087 0.6773988
[1] 0.9772799
```

보기 쉽게 표를 가로로 눕혀 찍었으니 한 열이 한 sample이야. `(Intercept)` 줄은 모든 sample에서 1이고, `conditionStarvation` 줄은 Starvation sample(4–6번)에서만 1이지. 그래서 Ctrl에서는 $\log\mu = b_0$, Starvation에서는 $\log\mu = b_0 + b_S$가 되고, 둘을 빼면 $b_S=\log(\mu_{\mathrm{Starvation}}/\mu_{\mathrm{Ctrl}})$이야. 즉 $b_S$는 "Starvation이 Ctrl의 몇 배인가"를 자연로그로 적은 값이야. 이렇게 비교의 기준이 되는 그룹(여기서는 Ctrl)을 reference라고 불러.

그런데 DESeq2가 보고하는 `log2FoldChange`는 자연로그가 아니라 log2 단위야. log2 fold change(LFC)는 두 조건 평균 비의 log2라서, LFC 1은 2배, −1은 절반, 0.5는 $2^{0.5}\approx1.414$배야. 두 단위가 섞이지 않도록 이 노트는 교재를 따라 자연로그 단위 계수를 b, log2 단위 계수를 β로 구분해 써($\beta=b/\log 2$). 유전자 A에 대입하면 $b_0=\log 106.667=4.6697$, $b_S=\log 1.96875=0.6774$, $\beta_S=0.6774/0.6931=0.9773$이야.

출력에서도 자연로그 계수 `bS` 0.6774를 log 2로 나누니 LFC 0.9773이 나왔어. 자연로그 계수를 그대로 LFC라고 읽으면 log 2 ≈ 0.693배만큼 어긋나니 조심해야 해. 그리고 여기서 그룹 평균의 비를 모형의 답으로 써도 되는지는 3절에서 확인해 볼게.

## 2. 세 그룹, donor, batch는 모형에 어떻게 들어갈까?

### 세 그룹은 계수 세 개로

이 시리즈의 세 조건은 Ctrl, Starvation, Starvation+Glucose야. Starvation+Glucose는 보통 배지에 glucose를 더 넣은 그룹이 아니라, 굶긴 세포에 glucose를 다시 넣은 조건이야. 이 노트에서는 줄여서 Glucose라고 부를게. Ctrl을 reference로 두면 모형은 이렇게 돼.

$$\log\mu_j=\log s_j+b_0+b_S\,I(S_j)+b_G\,I(G_j)$$

$I(S_j)$는 sample $j$가 Starvation이면 1, 아니면 0이고, $I(G_j)$는 Glucose에 대해 같은 뜻이야. 보고 싶은 비교는 계수를 조합해서 써. 어떤 두 조건(또는 계수 조합)을 비교할지 정하는 벡터를 contrast라고 불러.

| 비교 | 계수로 쓰면 | contrast $c$ (Intercept, Starvation, Glucose 순) |
|---|---|---|
| Starvation vs Ctrl | $\beta_S$ | $(0,1,0)$ |
| Glucose vs Ctrl | $\beta_G$ | $(0,0,1)$ |
| Glucose vs Starvation (rescue) | $\beta_G-\beta_S$ | $(0,-1,1)$ |

rescue 비교는 "glucose를 다시 넣으면 Starvation 효과가 되돌아오는가"를 물어. 이걸 보려고 Ctrl–Starvation과 Starvation–Glucose를 서로 다른 dataset처럼 따로 돌릴 필요는 없어. 세 그룹을 한 모형에 넣으면 dispersion과 size factor를 모든 sample에서 함께 추정하거든. rescue의 SE는 6절에서 다시 볼게.

### 같은 donor에서 세 조건을 얻었다면 `~ pair + condition`

donor A와 B에서 각각 세 조건을 얻었다고 해 볼게. 이때 기본 design은 `~ pair + condition`이야.

```r
mdp <- data.frame(pair = factor(rep(c("A", "B"), each = 3)),
                  condition = factor(rep(c("Ctrl", "Starvation", "Glucose"), 2),
                                     levels = c("Ctrl", "Starvation", "Glucose")))
rownames(mdp) <- paste(mdp$pair, mdp$condition, sep = "-")
Xp <- model.matrix(~ pair + condition, mdp)
Xp[, ]
cat("rank =", qr(Xp)$rank, " ncol =", ncol(Xp), " n =", nrow(Xp),
    " residual df =", nrow(Xp) - qr(Xp)$rank, "\n")
```

```
             (Intercept) pairB conditionStarvation conditionGlucose
A-Ctrl                 1     0                   0                0
A-Starvation           1     0                   1                0
A-Glucose              1     0                   0                1
B-Ctrl                 1     1                   0                0
B-Starvation           1     1                   1                0
B-Glucose              1     1                   0                1
rank = 4  ncol = 4  n = 6  residual df = 2
```

`pairB` 열은 "donor B는 원래 발현이 높다" 같은 donor 사이의 차이를 흡수하고, condition 열은 같은 donor 안에서 처리 후 달라진 부분을 설명해. 그래서 condition 계수는 donor 차이를 뺀 처리 효과가 돼.

이 모형은 donor마다 열을 하나 두고 그 계수를 따로 추정해. 이런 방식을 고정효과(fixed effect)라고 하는데, sample 사이의 임의의 상관이나 donor마다 다른 처리 효과의 크기까지 학습하는 random-effects 모형과는 달라. DESeq2의 design은 고정효과 열로만 펼쳐져(더 깊이 보기 A).

pair ID는 실제로 짝이 있을 때만 붙여야 해. 독립인 sample에 순서가 같다는 이유만으로 붙이면 데이터에 없는 구조를 모형에 강제하게 되거든. 기술 반복이나 반복 측정처럼 의존 구조가 복잡할 때는 어떻게 하는지 더 깊이 보기 H에 적어 뒀어.

출력 마지막 줄의 residual df(잔차 자유도)는 sample 수 n에서 독립적인 계수 수를 뺀 값이야. 독립적인 계수 수는 X의 rank, 곧 X의 열 가운데 서로 다른 정보를 담은 열의 수야. 계수를 다 추정하고 남은 정보로 replicate 사이의 퍼짐, 즉 dispersion을 추정하는데, 여기서는 6 − 4 = 2가 남아. 이 행렬은 교재 4.5의 두 검사 `qr(X)$rank == ncol(X)`와 `nrow(X) > ncol(X)`도 모두 통과해.

### batch가 조건과 완전히 겹친다면?

이번에는 Ctrl이 전부 batch 1에서, Starvation이 전부 batch 2에서 나왔다고 해 볼게.

```r
suppressPackageStartupMessages(library(DESeq2))
mdc <- data.frame(batch = factor(c(1, 1, 1, 2, 2, 2)),
                  condition = factor(rep(c("Ctrl", "Starvation"), each = 3)))
Xc <- model.matrix(~ batch + condition, mdc)
t(Xc[, ])                                    # 열 = sample
cat("rank =", qr(Xc)$rank, " ncol =", ncol(Xc), "\n")
set.seed(1); cts <- matrix(rpois(20 * 6, 50), nrow = 20)
dds_c <- try(DESeqDataSetFromMatrix(cts, mdc, ~ batch + condition))
```

```
                    1 2 3 4 5 6
(Intercept)         1 1 1 1 1 1
batch2              0 0 0 1 1 1
conditionStarvation 0 0 0 1 1 1
rank = 2  ncol = 3
Error in checkFullRank(modelMatrix) :
  the model matrix is not full rank, so the model cannot be fit as specified.
  One or more variables or interaction terms in the design formula are linear
  combinations of the others and must be removed.

  Please read the vignette section 'Model matrix not full rank':

  vignette('DESeq2')
```

`batch2` 줄과 `conditionStarvation` 줄이 똑같지. X의 열은 3개인데 rank는 2야. 차이가 batch 때문인지 Starvation 때문인지 가를 정보가 데이터에 없는 거야. `~ batch + condition`이라고 적어 준다고 이 사실이 바뀌지는 않아서, DESeq2는 dataset 객체를 만드는 단계에서 바로 멈춰. 옵션으로 풀 문제가 아니라 실험 설계에서 정보가 빠진 거야.

반대쪽 극단도 있어. sample 수와 계수 수가 같으면(n = p, saturated design) rank는 온전하지만 residual df가 0이야. replicate 사이의 퍼짐을 볼 정보가 남지 않으니, 이번에는 dispersion 추정 단계에서 멈춰(더 깊이 보기 A).

교재 4장의 정리를 빌리면, design은 무엇을 설명하고 무엇을 남길지 정하는 일이야. 조건별로 나눈 평균만 비교하는 게 아니라, 전체 sample에 같은 설계식의 계수를 적합하는 거지.

## 3. DESeq2는 계수를 어떻게 찾을까?

### dispersion을 구했는데 왜 계수를 다시 맞출까?

[03 노트](03_dispersion_estimation.md)에서 gene-wise dispersion을 구할 때 쓴 평균은 초기 α로 얻은 임시값이었어. 그 뒤 전체 유전자의 경향(trend와 prior)을 반영해 최종 α가 정해졌지. 이 값은 MAP, 곧 likelihood와 prior를 함께 고려했을 때 가장 그럴듯한 값이야. 유전자 A의 최종 α는 교재가 학습용으로 지정한 prior N(log 0.08, 0.7²)를 붙여 얻은 0.053147이야. α가 바뀌면 sample마다의 가중치가 바뀌고 계수와 SE도 따라 바뀌니까, 최종 α를 고정한 채 계수 b만 다시 찾는 거야.

이 단계에서 LFC에도 prior를 붙일지는 `betaPrior` 옵션이 정하는데, 기본값이 `FALSE`야. 그래서 기본 설정에서는 이 단계의 LFC에 prior shrinkage가 들어가지 않아. shrinkage는 정보가 적은 추정치를 전체 경향 쪽으로 당기는 것이고, LFC shrinkage는 [07 노트](07_lfc_shrinkage_and_qc.md)에서 따로 다뤄.

다만 계산을 안정시키는 아주 작은 장치 몇 개는 남아 있어. 이 노트에서 숫자로 보게 되는 건 ridge인데, 계수가 0에서 멀어질수록 조금씩 벌점을 주는 항이야. DESeq2의 ridge는 아주 약해서($\lambda=10^{-6}$) 유전자 A에서는 LFC를 소수점 일곱째 자리에서, SE를 여덟째 자리에서 바꾸는 정도야(더 깊이 보기 A, D).

### 어떤 b가 가장 그럴듯할까?

likelihood는 관측값을 고정해 두고, parameter 후보가 그 관측값을 얼마나 그럴듯하게 만드는지 나타내는 값이야([03 노트](03_dispersion_estimation.md)). α를 고정했을 때 likelihood를 가장 크게 만드는 b가 계수의 MLE야. MLE는 log likelihood를 b로 미분한 기울기가 0인 곳에 있고, 이 기울기를 score라고 불러.

$$U(b)=\sum_j x_j\,\frac{K_j-\mu_j}{1+\alpha\mu_j}$$

$K_j-\mu_j$는 관측 count와 모형이 예상한 count의 차이야. 분모 $1+\alpha\mu_j$는 그 차이를 얼마나 놀랍게 볼지 정하는데, α가 크면 같은 차이도 덜 놀라워. 앞의 $x_j$는 그 차이를 어느 계수에 반영할지 정하고.

관측이 예상보다 크면 score가 양수가 되어 평균을 올리는 쪽으로 b를 밀고, 작으면 내리는 쪽으로 밀어. 예를 들어 유전자 A에서 아래 반복의 시작값을 넣으면 Starvation 계수 방향의 score가 0.47이야. 아직 0이 아니니 b를 더 옮겨야 하지.

단순 평균을 구할 때 미분을 0으로 두는 것과 같은 원리야. 하지만 log link와 NB 분산 때문에 이 식은 일반적으로 한 번에 풀리지 않아서, 반복 계산을 해.

<details>
<summary>score 식은 어떻게 나올까</summary>

log likelihood는 $\ell(b)=\sum_j\log \mathrm{NB}(K_j;\mu_j,\alpha)$야. NB의 log 확률을 $\mu_j$로 미분하면 $\partial\ell/\partial\mu_j=(K_j-\mu_j)/(\mu_j(1+\alpha\mu_j))$이고, log link에서는 $\mu_j=e^{x_j^{T}b}$(offset 생략)이니 $\partial\mu_j/\partial b=\mu_j x_j$야. 둘을 곱하면 $\mu_j$가 약분되어 $U(b)=\sum_j x_j(K_j-\mu_j)/(1+\alpha\mu_j)$가 되고, MLE에서는 $U(\hat b)=0$이야.

</details>

### IRLS: 가중치를 고쳐 가며 되풀이하는 최소제곱

IRLS(반복 가중 최소제곱)는 계수를 찾는 반복 계산이야. 현재 추정치 $b^{(t)}$에서 한 번 반복하는 과정은 이래.

$$\eta_j=\log s_j+x_j^{T}b^{(t)},\qquad \mu_j=e^{\eta_j}$$

$$z_j=\eta_j+\frac{K_j-\mu_j}{\mu_j},\qquad w_j=\frac{\mu_j}{1+\alpha\mu_j}$$

$$b^{(t+1)}=\left(X^{T}WX\right)^{-1}X^{T}W\,(z-o)$$

- $\eta_j$는 현재 계수로 계산한 log 기대 count야.
- $z_j$는 working response라고 부르는데, $\log K_j$를 현재 $\mu_j$ 주변에서 직선으로 근사한 값이야.
- $w_j$는 sample $j$의 가중치, 곧 이 sample이 계수에 주는 정보의 양이야. $(d\mu/d\eta)^2/\mathrm{Var}(K)=\mu^2/(\mu+\alpha\mu^2)$를 정리한 거야.
- $W$는 $w_j$를 대각선에 놓은 행렬이고, $o$는 offset $\log s_j$를 모은 벡터야.
- 마지막 식은 $z-o$를 X에 가중 최소제곱으로 맞춘 해야. 가중 최소제곱은 오차 제곱에 sample마다 $w_j$를 곱해서, 정보가 많은 sample을 더 믿는 최소제곱이야.

말로 풀면 한 번의 반복은 "지금 계수로 기대 count 계산 → 가중치와 working response 계산 → 가중 최소제곱으로 계수 갱신"이야. b가 더 이상 바뀌지 않으면($b^{(t+1)}=b^{(t)}$) 멈추는데, 그 자리에서는 $X^{T}W(z-\eta)=U(b)=0$, 곧 score가 0이야. DESeq2는 언제 멈출지를 deviance의 상대 변화로 판정해(더 깊이 보기 B).

유전자 A에 α = 0.053147을 넣고 직접 반복해 볼게. 시작값은 DESeq2처럼 log(count + 0.1)의 최소제곱해로 잡았어.

```r
alpha <- 0.053147                            # 03 노트의 최종 α (MAP)
o <- log(rep(1, 6))                          # offset = log(size factor) = 0
b <- qr.solve(X, log(k + 0.1))               # 시작값: log(count + 0.1)의 최소제곱해
for (t in 0:3) {
  eta <- as.numeric(o + X %*% b); mu <- exp(eta)
  U <- crossprod(X, (k - mu) / (1 + alpha * mu))            # score
  cat(sprintf("iter %d  b_S = %.9f  beta_S = %.9f  score_S = %9.2e\n",
              t, b[2], b[2] / log(2), U[2]))
  z <- eta + (k - mu) / mu                   # working response
  w <- mu / (1 + alpha * mu)                 # weight
  b <- as.numeric(solve(crossprod(X, X * w), crossprod(X, w * (z - o))))
}
exp(b)                                       # Ctrl 평균, Starvation/Ctrl 비
```

```
iter 0  b_S = 0.679599032  beta_S = 0.980454153  score_S =  4.70e-01
iter 1  b_S = 0.677376390  beta_S = 0.977247559  score_S = -2.13e-03
iter 2  b_S = 0.677398822  beta_S = 0.977279922  score_S = -4.36e-08
iter 3  b_S = 0.677398824  beta_S = 0.977279923  score_S =  1.41e-14
[1] 106.66667   1.96875
```

직접 돌려 보니 score가 반복마다 0에 가까워지다가 네 번째 줄(iter 3)에서 $10^{-14}$ 수준이 되고, β_S는 0.977279923에서 멈췄어. 마지막 줄의 exp(b)는 106.667과 1.96875(= 210 / 106.667)야. 두 그룹이고 size factor가 모두 같으면 MLE 평균이 그룹 산술평균과 같거든(교재 16.1). 1절의 단순 계산과 같은 답이 나온 이유가 이거야.

### DESeq2에게 같은 일을 시키면?

같은 α를 `dispersions()`로 고정하고 DESeq2의 Wald 단계(`nbinomWaldTest()`)만 돌려 봐. 이 함수가 계수를 적합하고 SE와 Wald 검정 결과까지 만들어. Wald 검정은 (추정치 − 0) / SE를 표준정규분포와 비교하는 검정이야.

```r
cts1 <- matrix(as.integer(k), nrow = 1, dimnames = list("geneA", paste0("s", 1:6)))
dds <- DESeqDataSetFromMatrix(cts1, cd, ~ condition)
sizeFactors(dds) <- rep(1, 6)
dispersions(dds) <- alpha                    # dispersion 추정을 건너뛰고 α를 고정
dds <- nbinomWaldTest(dds, quiet = TRUE)
res <- results(dds, independentFiltering = FALSE, cooksCutoff = FALSE)
print(as.data.frame(res)[, c("log2FoldChange", "lfcSE", "stat", "pvalue")], digits = 9)
cat("betaIter =", mcols(dds)$betaIter, "\n")
```

```
      log2FoldChange       lfcSE       stat         pvalue
geneA    0.977280132 0.289056569 3.38093037 0.000722408466
betaIter = 2
```

LFC 0.977280과 lfcSE 0.289057이 교재 16.4의 숫자 그대로 나왔어. 직접 계산한 β_S 0.977279923과는 $2.1\times10^{-7}$ 차이가 나는데, 앞에서 말한 아주 약한 ridge 때문이야(더 깊이 보기 D). `betaIter`가 2라는 건 DESeq2가 두 번 반복하고 멈췄다는 뜻이고, `stat`과 `pvalue`는 [05 노트](05_wald_vs_lrt.md)에서 다뤄.

## 4. SE는 어디서 나오고, 왜 α가 크면 커질까?

SE(표준오차)는 같은 실험을 반복하면 추정치가 얼마나 흔들릴지를 나타내는 값이야. 유전자 A의 lfcSE 0.289가 어디서 나왔는지 따라가 볼게.

### 가중치를 모은 정보량

IRLS에서 쓴 가중치 $w_j$가 그대로 SE의 재료야. 각 sample이 주는 정보를 design에 맞게 모은 것이 Fisher information이고.

$$I(b)=X^{T}WX,\qquad \widehat{\mathrm{Cov}}(\hat b)\approx\left(X^{T}WX\right)^{-1}$$

Fisher information $I(b)$는 likelihood의 봉우리가 계수 방향으로 얼마나 뾰족한지를 재는 양이야. 뾰족할수록 추정치가 덜 흔들리지. $(\cdot)^{-1}$은 역행렬인데, 숫자 하나라면 $1/x$에 해당하니 정보가 많을수록 분산이 작다는 뜻이야. 그렇게 얻은 $\widehat{\mathrm{Cov}}(\hat b)$는 계수 추정치의 공분산 행렬로, 대각선은 각 계수의 분산이고 대각선 밖은 두 계수가 함께 흔들리는 정도(공분산)야.

SE는 대각선 값의 제곱근이고, log2 단위로 바꾸려면 log 2로 나눠. 두 그룹 design이라면 역행렬을 손으로도 풀 수 있어.

$$\mathrm{Var}(\hat b_S)=\frac{1}{\sum_{j\in \mathrm{Ctrl}}w_j}+\frac{1}{\sum_{j\in \mathrm{Starvation}}w_j}$$

첫째 항은 Ctrl 평균의 log가 흔들리는 정도, 둘째 항은 Starvation 평균의 log가 흔들리는 정도야. $b_S$는 두 log 평균의 차이라서 두 흔들림이 더해지고, 그룹의 가중치 합이 클수록 그 항이 작아져.

유전자 A에 대입해 볼게. Ctrl sample 하나의 가중치는 $106.667/(1+0.053147\times106.667)=15.99$이고, Starvation은 $210/(1+0.053147\times210)=17.27$이야. 그룹별로 합하면 47.98과 51.81이야.

$$\mathrm{Var}(\hat b_S)=\frac{1}{47.98}+\frac{1}{51.81}=0.040144,\qquad SE(\hat\beta_S)=\frac{\sqrt{0.040144}}{\log 2}=\frac{0.200359}{0.693147}=0.289057$$

```r
mu <- exp(as.numeric(o + X %*% b))           # 수렴한 기대 count
w  <- mu / (1 + alpha * mu)
round(w, 4)
XtWX <- crossprod(X, X * w)                  # Fisher information  X^T W X
V <- solve(XtWX)                             # 역행렬 = 계수의 근사 공분산 (자연로그)
V
sqrt(diag(V)) / log(2)                       # SE (log2 단위)
1 / sum(w[1:3]) + 1 / sum(w[4:6])            # Var(b_S)를 그룹별 weight 합으로
```

```
[1] 15.9944 15.9944 15.9944 17.2685 17.2685 17.2685
                    (Intercept) conditionStarvation
(Intercept)          0.02084067         -0.02084067
conditionStarvation -0.02084067          0.04014363
        (Intercept) conditionStarvation
          0.2082717           0.2890566
[1] 0.04014363
```

공분산 행렬의 오른쪽 아래 값 0.04014363이 손계산과 같고, 그 제곱근을 log 2로 나눈 0.2890566이 바로 교재의 SE 0.289057이야.

### α만 바꾸면 SE는 어떻게 될까?

평균(106.667, 210)과 design은 그대로 두고 α만 바꿔 볼게. [03 노트](03_dispersion_estimation.md)(교재 16.2–16.3)에서 구한 세 후보, 곧 일반 NB MLE 0.014786, Cox-Reid 조정 0.025385, MAP 0.053147을 넣고, 비교용으로 0, 0.1, 0.5도 넣었어. Cox-Reid 조정값은 평균을 같은 데이터로 추정해서 생기는 과소추정을 줄인 값이야.

```r
mu_fix <- rep(c(320 / 3, 210), each = 3); bS <- log(210 / (320 / 3))
tab <- t(sapply(c(0, 0.014786, 0.025385, 0.053147, 0.1, 0.5), function(a) {
  w <- mu_fix / (1 + a * mu_fix); V <- solve(crossprod(X, X * w))
  c(alpha = a, w_Ctrl = w[1], w_Starvation = w[4], LFC = bS / log(2),
    SE_log2 = sqrt(V[2, 2]) / log(2), Z = bS / sqrt(V[2, 2]))
}))
round(tab, 6)
```

```
        alpha     w_Ctrl w_Starvation     LFC  SE_log2        Z
[1,] 0.000000 106.666667   210.000000 0.97728 0.099036 9.867972
[2,] 0.014786  41.389015    51.156378 0.97728 0.174140 5.612032
[3,] 0.025385  28.768700    33.170901 0.97728 0.212207 4.605318
[4,] 0.053147  15.994370    17.268501 0.97728 0.289057 3.380929
[5,] 0.100000   9.142857     9.545455 0.97728 0.385443 2.535475
[6,] 0.500000   1.963190     1.981132 0.97728 0.838807 1.165083
```

`LFC` 열은 그대로인데 `SE_log2` 열만 커지는 게 보이지. α = 0(Poisson)이면 SE가 0.099인데, MLE → Cox-Reid → MAP으로 가면 0.174 → 0.212 → 0.289로 커지고 Z는 5.61 → 4.61 → 3.38로 작아져. 03 노트에서 dispersion을 어떻게 추정했는지가 검정 결과에 이렇게 들어오는 거야. 평균과 design을 고정하면 이 사슬은 한 방향이야.

α 증가 → $w_j$ 감소 → 정보량 감소 → 계수 분산 증가 → SE 증가

`LFC` 열이 상수인 건 코드에서 평균을 고정했으니 당연해. 이 예제에서 평균을 고정해도 되는 건, 두 그룹이고 size factor가 모두 같아서 어떤 α에서도 MLE 평균이 그룹 산술평균이기 때문이야. 일반적으로는 그렇지 않아. size factor를 sample마다 다르게 바꾸면 α가 0.0148에서 0.5로 갈 때 LFC도 0.965에서 0.946으로 움직이거든(더 깊이 보기 E). 실제 데이터에서는 α가 바뀌면 LFC와 SE가 함께 바뀔 수 있어.

### 이름이 비슷한 네 가지 "퍼짐" (교재 8장 정리)

| 양 | 식 | 무엇의 변동인가 |
|---|---|---|
| count의 분산 | $\mathrm{Var}(K_{ij})=\mu_{ij}+\alpha_i\mu_{ij}^2$ | sample 하나의 count가 평균 주위로 흩어지는 정도 |
| dispersion | $\alpha_i$ | Poisson을 넘는 biological 변동의 크기 (유전자별 값) |
| 계수의 추정 분산 | $\mathrm{Var}(\hat b_k)\approx[(X^{T}WX)^{-1}]_{kk}$ | 실험을 반복할 때 $\hat b_k$가 흩어지는 정도 (자연로그) |
| log2 계수의 SE | $SE(\hat\beta_k)=\sqrt{\mathrm{Var}(\hat b_k)}/\log 2$ | `results()`의 `lfcSE` |

최종 α가 likelihood의 곡률과 정보량을 바꾸니까, Wald 검정 결과도 따라서 바뀌어.

## 5. 더 깊게 sequencing하면 SE가 계속 줄어들까?

나도 처음엔 깊게 읽을수록 SE가 계속 줄어들 거라고 생각했어. α = 0.1인 유전자로 계산해 볼게. 기대 count가 100이면 sample 하나의 가중치는 100 / (1 + 10) = 9.09야. 같은 sample을 50배 깊게 읽어 기대 count가 5000이 되어도 가중치는 9.98로, 10을 넘지 못해.

$$w=\frac{\mu}{1+\alpha\mu}<\frac{1}{\alpha}\quad(\alpha>0),\qquad \lim_{\mu\to\infty}\frac{\mu}{1+\alpha\mu}=\frac{1}{\alpha}$$

μ가 작으면 분모의 αμ가 1보다 훨씬 작아서 w ≈ μ야. 이 구간에서는 깊게 읽을수록 정보가 거의 비례해서 늘어나. 그런데 μ가 커지면 αμ가 분모를 차지해 w가 1/α에 붙어 버려. biological 변동(α)은 깊게 읽는다고 줄지 않으니까.

![Information weight vs expected count for three dispersion values](../figures/04_information_weight.png)

교재 그림 2를 다시 그린 거야. 점선(Poisson, α = 0)은 μ를 따라 계속 올라가지만, α가 양수인 세 곡선은 각자의 1/α(가는 파선)에서 평평해져.

<details>
<summary>그림을 만든 코드</summary>

```python
# notes/ 폴더에서 실행한다. 그림은 ../figures/에 저장된다.
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out = "../figures/04_information_weight.png"
mu = np.logspace(0, np.log10(5000), 400)
alphas = [(0.01, "#2a78d6"), (0.1, "#eb6834"), (0.5, "#1baf7a")]

fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150, facecolor="white")
ax.plot(mu, mu, color="#9a9994", lw=1.2, ls=":", label="alpha = 0 (Poisson): w = mu")
for a, col in alphas:
    ax.plot(mu, mu / (1 + a * mu), color=col, lw=2, label=f"alpha = {a}")
    ax.axhline(1 / a, color=col, lw=0.8, ls="--", alpha=0.6)
    ax.text(5200, 1 / a, f"1/alpha = {1 / a:g}", va="center", ha="left", fontsize=9, color="#52514e")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 5000); ax.set_ylim(0.5, 300)
ax.set_xlabel("Expected count mu")
ax.set_ylabel("Information weight  w = mu / (1 + alpha * mu)")
ax.grid(True, which="major", color="#e6e5e0", lw=0.6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(loc="upper left", frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(out, dpi=150, facecolor="white")
print("saved", out, plt.imread(out).shape)
```

```
saved ../figures/04_information_weight.png (675, 1050, 4)
```

</details>

그럼 정보를 늘리려면 뭘 해야 할까? 한 그룹 log 평균의 분산은 그 그룹 가중치 합의 역수($1/\sum_j w_j$)야. α = 0.1, 기대 count 100인 sample 3개가 있으면 가중치 합은 27.3이야. 같은 sample을 두 배 깊게 읽으면 28.6으로 약 5% 늘지만, replicate를 두 배로 늘리면 54.5로 두 배가 돼. biological dispersion이 있으면 한 sample을 더 깊게 읽는 것보다 독립 replicate를 더하는 쪽이 $\sum_j w_j$를 키우는 거지. 다만 이건 식의 성질이지, 특정 실험의 최적 depth를 정하는 공식은 아니야.

<details>
<summary>그림의 값과 depth·replicate 비교를 숫자로 확인하기</summary>

```r
alpha_grid <- c(0.01, 0.1, 0.5); mu_grid <- c(1, 3, 10, 30, 100, 300, 1000, 5000)
tab <- outer(mu_grid, alpha_grid, function(m, a) m / (1 + a * m))
dimnames(tab) <- list(paste0("mu=", mu_grid), paste0("alpha=", alpha_grid))
round(rbind(tab, `limit 1/alpha` = 1 / alpha_grid), 3)
```

```
              alpha=0.01 alpha=0.1 alpha=0.5
mu=1               0.990     0.909     0.667
mu=3               2.913     2.308     1.200
mu=10              9.091     5.000     1.667
mu=30             23.077     7.500     1.875
mu=100            50.000     9.091     1.961
mu=300            75.000     9.677     1.987
mu=1000           90.909     9.901     1.996
mu=5000           98.039     9.980     1.999
limit 1/alpha    100.000    10.000     2.000
```

α = 0.1 열을 보면 μ = 100에서 이미 상한 10의 91%(9.091)이고, μ를 50배 늘려도 9.980이야.

```r
a <- 0.1                                     # 한 그룹 3개 sample의 weight 합 (클수록 SE가 작다)
c(now           = 3 * 100 / (1 + a * 100),   # 지금: 기대 count 100
  depth_x2      = 3 * 200 / (1 + a * 200),   # 같은 sample을 두 배 깊게
  replicates_x2 = 6 * 100 / (1 + a * 100))   # replicate를 두 배로
```

```
          now      depth_x2 replicates_x2
     27.27273      28.57143      54.54545
```

</details>

## 6. rescue 비교(Glucose − Starvation)의 SE는 어떻게 구할까?

### 두 SE를 그냥 합치면?

이번엔 세 그룹 예제가 필요해. DESeq2의 `makeExampleDESeqDataSet()`은 두 그룹만 만들어 주기 때문에, 유전자 2000개, sample 9개의 가상 데이터를 만든 뒤 condition을 Ctrl / Starvation / Glucose 세 그룹(각 3개)으로 덮어썼어. 그다음 `DESeq()`로 size factor, dispersion, 계수와 SE까지 한 번에 구했어. 아래의 gene4는 이 가상 데이터의 유전자 하나야.

```r
set.seed(1)
dds3 <- makeExampleDESeqDataSet(n = 2000, m = 9, betaSD = 1)
dds3$condition <- factor(rep(c("Ctrl", "Starvation", "Glucose"), each = 3),
                         levels = c("Ctrl", "Starvation", "Glucose"))
design(dds3) <- ~ condition
dds3 <- DESeq(dds3, quiet = TRUE)
rS  <- results(dds3, name = "condition_Starvation_vs_Ctrl")
rG  <- results(dds3, name = "condition_Glucose_vs_Ctrl")
rGS <- results(dds3, contrast = c("condition", "Glucose", "Starvation"))
g <- "gene4"; se_S <- rS[g, "lfcSE"]; se_G <- rG[g, "lfcSE"]
round(c(SE_S = se_S, SE_G = se_G, SE_GvsS_results = rGS[g, "lfcSE"],
        naive_sum = se_S + se_G, naive_sqrt = sqrt(se_S^2 + se_G^2)), 4)
```

```
           SE_S            SE_G SE_GvsS_results       naive_sum      naive_sqrt
         0.5419          0.5419          0.5415          1.0838          0.7664
```

gene4에서 Starvation vs Ctrl의 SE와 Glucose vs Ctrl의 SE는 둘 다 0.5419야. 두 SE로 차이의 SE를 짐작할 때 흔히 쓰는 방법은 더하기(`naive_sum`, 1.084)와 제곱합의 제곱근(`naive_sqrt`, 0.766)이지. 솔직히 나도 처음엔 제곱합의 제곱근이면 될 줄 알았어. 그런데 `results(contrast = c("condition", "Glucose", "Starvation"))`가 준 답은 0.5415야. 두 방법 모두 틀린 거야.

### 빠진 조각은 공분산

contrast로 쓴 비교의 추정치와 SE는 이렇게 구해.

$$\hat\delta=c^{T}\hat b,\qquad SE(\hat\delta)=\sqrt{c^{T}\,\widehat{\mathrm{Cov}}(\hat b)\,c}$$

$c=(0,-1,1)$을 넣고 풀면 이렇게 돼.

$$\mathrm{Var}(\hat b_G-\hat b_S)=\mathrm{Var}(\hat b_G)+\mathrm{Var}(\hat b_S)-2\,\mathrm{Cov}(\hat b_G,\hat b_S)$$

공분산(covariance)은 두 추정치가 같은 방향으로 함께 흔들리는 정도야. 제곱합의 제곱근은 공분산이 0일 때만 맞고, 더하기는 두 추정치의 상관이 −1일 때만 도달하는 상한이야.

여기서 공분산이 양수인 이유는 간단해. $\hat\beta_S$와 $\hat\beta_G$는 둘 다 "그 그룹 − Ctrl"이라서 같은 Ctrl 평균 추정치를 빼. 그래서 Ctrl 추정치가 우연히 높게 나오면 두 계수가 함께 낮아지고, 둘의 차이를 구하면 이 공통 흔들림이 지워져. `~ condition`처럼 요인이 하나면 이 공분산은 절편의 분산과 정확히 같아.

$$\mathrm{Cov}(\hat b_G,\hat b_S)=\mathrm{Var}(\hat b_0)=\frac{1}{\sum_{j\in\mathrm{Ctrl}}w_j}\quad\Rightarrow\quad \mathrm{Var}(\hat b_G-\hat b_S)=\frac{1}{\sum_{j\in\mathrm{Glucose}}w_j}+\frac{1}{\sum_{j\in\mathrm{Starvation}}w_j}$$

Ctrl 항이 빠지는 거야. gene4의 공분산 행렬을 직접 계산해서 확인해 볼게. log2 단위로 바꿀 때 공분산은 $(\log 2)^2$로 나눠.

```r
X3 <- attr(dds3, "modelMatrix")
mu <- assays(dds3)[["mu"]][g, ]
w  <- mu / (1 + dispersions(dds3)[rownames(dds3) == g] * mu)
V3 <- solve(crossprod(X3, X3 * w)) / log(2)^2        # log2 단위 공분산 행렬
cc <- c(0, -1, 1)                                    # Glucose - Starvation
c(Var_S = V3[2, 2], Var_G = V3[3, 3], Cov_GS = V3[2, 3], Var_b0 = V3[1, 1],
  SE_contrast = sqrt(drop(t(cc) %*% V3 %*% cc)))
```

```
      Var_S       Var_G      Cov_GS      Var_b0 SE_contrast
  0.2937042   0.2936522   0.1470652   0.1470652   0.5415034
```

`Cov_GS`와 `Var_b0`가 0.1470652로 똑같이 나왔어. 대입하면 $\sqrt{0.2937+0.2937-2\times0.1471}=0.5415$이고, 직접 계산한 0.5415034는 `results()`의 lfcSE 0.5415033과 ridge 몫만큼만 달라.

그러니까 `results(contrast = ...)`는 저장된 두 LFC를 빼기만 하는 게 아니야. LFC는 빼지만, SE는 $c^{T}\widehat{\mathrm{Cov}}(\hat b)\,c$로 다시 계산해(더 깊이 보기 C). `c("condition", "Glucose", "Starvation")`, `list(...)`, `c(0, -1, 1)` 세 가지 표기 모두 같은 LFC와 SE를 내고(더 깊이 보기 F).

덧붙이면, "공분산 = 절편의 분산 > 0"은 같은 요인의 두 level이 같은 reference를 공유할 때의 성질이야. `~ pair + condition`처럼 다른 요인이 있으면 이 등식이 성립하지 않고(공분산 0.01, 절편 분산 0.0133), 서로 다른 요인의 계수 사이 공분산은 음수일 수도 있어. 두 경우 모두 연습문제 12에서 숫자로 볼 수 있어.

그리고 rescue SE 식은 같은 가중치로 Starvation과 Glucose 두 그룹만 비교했을 때의 식과 같아. 그러니 $-2\,\mathrm{Cov}$ 항은 Ctrl 기준 계수에서 rescue SE를 바르게 읽는 방법이지, 한 모형이 주는 이점은 아니야. 한 모형의 이점은 2절에서 말한 dispersion·size factor 공동 추정에 있어(더 깊이 보기 H).

## 정리

- DESeq2의 GLM은 count 평균의 log를 offset(log size factor)과 design matrix 열들의 합으로 써. `log2FoldChange`는 자연로그 계수 b를 log 2로 나눈 β야.
- 세 그룹, donor, batch는 모두 design matrix의 열로 들어가. `~ pair + condition`은 donor를 고정효과로 넣고 log scale에서 donor 효과와 처리 효과를 더하기만 해. batch가 조건과 완전히 겹치거나 sample 수가 계수 수와 같으면 어떤 옵션으로도 풀 수 없어.
- 계수는 최종 α를 고정한 채 IRLS로 찾고, SE는 가중치 $w=\mu/(1+\alpha\mu)$로 만든 정보량 $X^{T}WX$의 역행렬에서 나와. 유전자 A는 두 번 반복으로 LFC 0.977280, SE 0.289057에 닿아.
- 평균과 design을 고정하면 α가 클수록 SE가 커져. 유전자 A에서는 MLE, Cox-Reid, MAP α를 차례로 넣었을 때 0.174, 0.212, 0.289였어. 게다가 관측 하나의 가중치는 1/α를 넘지 못해서, 기대 count가 이미 1/α 근처까지 크면 깊게 읽기보다 replicate를 늘리는 쪽이 정보를 훨씬 많이 늘려(α = 0.1, μ = 100에서 depth 2배는 약 5%, replicate 2배는 2배).
- 계수 차이의 SE에는 공분산이 필요하고, `results(contrast=)`가 그걸 계산해 줘. gene4의 rescue SE는 0.766이나 1.084가 아니라 0.5415야.

## 연습문제

교재 부록 A의 문제를 옮겼어. 문장은 그대로 두고 조사 띄어쓰기만 붙였어.

**문제 5.** 3-group 모형에서 Ctrl을 reference로 두었고 $\beta_S=-1.2,\ \beta_G=-0.3$이다. Glucose−Starvation의 log2FC와 Glucose−Ctrl의 log2FC는 각각 얼마인가?

<details>
<summary>풀이</summary>

```r
cat("Glucose-Starvation =", -0.3 - (-1.2), " Glucose-Ctrl =", -0.3, " fold =", 2^(-0.3 - (-1.2)), 2^(-0.3), "\n")
```

```
Glucose-Starvation = 0.9  Glucose-Ctrl = -0.3  fold = 1.866066 0.8122524
```

Glucose−Starvation은 $c=(0,-1,1)$이니 $\beta_G-\beta_S=-0.3-(-1.2)=0.9$, 약 1.87배야. Glucose−Ctrl은 $\beta_G=-0.3$ 그대로라 약 0.81배고. rescue 비교가 양수여도 Glucose 그룹은 여전히 Ctrl보다 조금 낮다는 뜻이야. 교재 정답(0.9, −0.3)과 같아.

</details>

**문제 11.** 같은 expected count μ=100에서 α가 0.01에서 0.1로 증가했다. 정보량 weight를 각각 계산하고 SE에 미치는 방향을 설명하라.

<details>
<summary>풀이</summary>

```r
cat("w:", 100 / (1 + 0.01 * 100), 100 / (1 + 0.1 * 100), " ratio =", (100 / 2) / (100 / 11), " SE ratio =", sqrt((100 / 2) / (100 / 11)), "\n")
```

```
w: 50 9.090909  ratio = 5.5  SE ratio = 2.345208
```

가중치는 $100/(1+1)=50$에서 $100/(1+10)\approx9.09$로 5.5배 줄어. 평균과 design을 고정하고 모든 sample의 μ가 같으면 $\mathrm{Var}(\hat b)\propto1/w$이니, SE는 $\sqrt{5.5}\approx2.35$배 커져. 교재 정답(50, 9.09, SE 증가)과 같고, 4절의 표가 같은 현상을 유전자 A에서 보여 줘.

</details>

**문제 12.** 두 coefficient의 SE가 각각 0.3과 0.4라고 해서 그 차이의 SE가 항상 0.5인가?

<details>
<summary>풀이</summary>

본문 2절의 `Xp`(paired design matrix)를 이어 써.

```r
cat("sum =", 0.3 + 0.4, " sqrt sum sq =", sqrt(0.3^2 + 0.4^2), " Cov=+0.05 ->", sqrt(0.09 + 0.16 - 0.1), " Cov=-0.05 ->", sqrt(0.09 + 0.16 + 0.1), "\n")
cat("SE(diff) at rho = -1, 0, +1:", sqrt(0.25 + 0.24), sqrt(0.25), sqrt(0.25 - 0.24), "\n")
Xb <- model.matrix(~ batch + cond, data.frame(batch = factor(c(1, 1, 1, 2, 2, 2)), cond = factor(c("C", "C", "S", "C", "S", "S"))))
V <- solve(crossprod(Xb) * 50); cat("~batch+cond, w = 50: Cov(batch2, condS) =", V[2, 3], " cor =", V[2, 3] / sqrt(V[2, 2] * V[3, 3]), "\n")
V <- solve(crossprod(Xp) * 50); cat("~pair+condition, w = 50: Cov(Starvation, Glucose) =", V[3, 4], " Var(Intercept) =", V[1, 1], "\n")
```

```
sum = 0.7  sqrt sum sq = 0.5  Cov=+0.05 -> 0.3872983  Cov=-0.05 -> 0.591608
SE(diff) at rho = -1, 0, +1: 0.7 0.5 0.1
~batch+cond, w = 50: Cov(batch2, condS) = -0.005  cor = -0.3333333
~pair+condition, w = 50: Cov(Starvation, Glucose) = 0.01  Var(Intercept) = 0.01333333
```

아니야. $\sqrt{0.3^2+0.4^2}=0.5$는 공분산이 0일 때만 맞아. 일반적으로는 $\sqrt{0.09+0.16-2\,\mathrm{Cov}}$이고, 상관 $\rho$가 −1에서 1 사이를 움직이면 0.7에서 0.1까지 될 수 있어.

- 요인이 하나인 `~ condition`에서 같은 reference(Ctrl)를 공유하는 두 level의 계수는 $\mathrm{Cov}=\mathrm{Var}(\hat b_0)>0$이야. rescue 비교가 이 경우라서 0.5보다 작은 쪽이 돼(6절 gene4: 0.5419, 0.5419 → 0.5415).
- 그렇다고 이게 reference와 비교하는 계수 표현(treatment coding) 전체의 성질은 아니야. 위 `~ batch + cond` 예에서는 batch 계수와 condition 계수의 공분산이 음수(상관 −1/3)라서, 차이의 SE가 $\sqrt{SE_1^2+SE_2^2}$보다 커.
- `~ pair + condition`에서는 두 condition 계수의 공분산(0.01)이 절편의 분산(0.0133)과 달라.

교재 정답("covariance가 0인 경우에만 0.5, 일반적으로는 $-2\,\mathrm{Cov}$ 항이 필요")과 같아.

</details>

## 더 깊이 보기

<details>
<summary>A. DESeq2 소스에서 본 design, 기본값, 수렴 처리</summary>

먼저 design이 model matrix가 되는 부분이야. `Rscript -e 'print(DESeq2:::fitNbinomGLMs)'`의 앞부분을 옮겼어.

```r
# 발췌 (실행하지 않음): print(DESeq2:::fitNbinomGLMs) 앞부분
    if (is.null(modelMatrix)) {
        modelAsFormula <- TRUE
        modelMatrix <- stats::model.matrix.default(modelFormula,
            data = as.data.frame(colData(object)))
    }
    ...
    if (renameCols) {
        convertNames <- renameModelMatrixColumns(colData(object), modelFormula)
```

- `DESeq2:::renameModelMatrixColumns`가 `conditionStarvation`을 `condition_Starvation_vs_Ctrl`로 바꿔(`paste0(v, "_", levels[-1], "_vs_", levels[1])`). 교재 4.2 표의 열 이름 "Starvation"은 이 이름을 줄여 쓴 거야.
- rank는 `DESeq2:::checkFullRank`가 검사해. `qr(modelMatrix)$rank < ncol(modelMatrix)`이면 constructor 단계에서 멈추고, n = p이면 `DESeq2:::checkForExperimentalReplicates`가 `nrow(modelMatrix) == ncol(modelMatrix)`을 보고 dispersion 단계에서 멈춰.
- design은 `stats::model.matrix.default`를 거쳐 고정효과 열로만 펼쳐지고, `DESeq()`의 인자 목록에도 random effect를 받는 인자가 없어.
- 본문 2절에서는 교재 4.4의 paired 행렬과 4.5의 confounded 예를 그대로 실행했어. paired 행렬은 교재 4.4 표와 열 단위로 같아.

아래 코드는 세 그룹 2반복의 rank, 교재 4.5의 두 검사, `DESeq()` 인자, saturated design의 오류, 그리고 이 노트가 인용하는 기본값을 확인해. 본문 2–3절의 `Xp`, `cts`, `dds`, `res`를 이어 써.

```r
X3g <- model.matrix(~ condition, data.frame(condition = factor(rep(c("Ctrl", "Starvation", "Glucose"), each = 2))))
cat("~condition, 3 groups x 2: rank =", qr(X3g)$rank, " residual df =", nrow(X3g) - qr(X3g)$rank, "\n")
stopifnot(qr(Xp)$rank == ncol(Xp), nrow(Xp) > ncol(Xp))    # 교재 4.5의 두 검사 (본문 2절의 Xp)
print(names(formals(DESeq)))
mds <- data.frame(condition = factor(c("A", "B", "C")))    # n = p = 3 인 saturated design
try(DESeq(DESeqDataSetFromMatrix(cts[, 1:3], mds, ~ condition), quiet = TRUE))
grep("missing(betaPrior)", deparse(body(DESeq)), fixed = TRUE, value = TRUE)
cat("betaPrior default (nbinomWaldTest):", formals(nbinomWaldTest)$betaPrior, "\n")
print(formals(DESeq)$minmu)
cat("minmu default (nbinomWaldTest, results):", formals(nbinomWaldTest)$minmu, formals(results)$minmu, "\n")
cat("log2(exp(1)) == 1/log(2):", isTRUE(all.equal(log2(exp(1)), 1 / log(2))), "\n")
cat("gene A: betaPriorVar =", attr(dds, "betaPriorVar"), "| LFC description:", mcols(res)$description[2], "\n")
```

```
~condition, 3 groups x 2: rank = 3  residual df = 3
 [1] "object"                  "test"
 [3] "fitType"                 "sfType"
 [5] "betaPrior"               "full"
 [7] "reduced"                 "quiet"
 [9] "minReplicatesForReplace" "modelMatrixType"
[11] "useT"                    "minmu"
[13] "parallel"                "BPPARAM"
Error in checkForExperimentalReplicates(object, modelMatrix) :

  The design matrix has the same number of samples and coefficients to fit,
  so estimation of dispersion is not possible. Treating samples
  as replicates was deprecated in v1.20 and no longer supported since v1.22.


[1] "    if (missing(betaPrior)) {"
betaPrior default (nbinomWaldTest): FALSE
if (fitType == "glmGamPoi") 1e-06 else 0.5
minmu default (nbinomWaldTest, results): 0.5 0.5
log2(exp(1)) == 1/log(2): TRUE
gene A: betaPriorVar = 1e+06 1e+06 | LFC description: log2 fold change (MLE): condition Starvation vs Ctrl
```

세 그룹 2반복(`~ condition`, n = 6)은 rank 3, residual df 3이고, paired 행렬은 두 `stopifnot` 검사를 조용히 통과했어.

다음은 최종 α, log2 변환, betaPrior, ridge, 수렴 처리가 들어 있는 `fitNbinomGLMs` 발췌야.

```r
# 발췌 (실행하지 않음): DESeq2:::fitNbinomGLMs
    if (missing(alpha_hat)) { alpha_hat <- dispersions(object) }
    if (missing(lambda)) { lambda <- rep(1e-06, ncol(modelMatrix)) }
    ...
    justIntercept <- if (modelAsFormula) { modelFormula == formula(~1) } else { ncol(modelMatrix) == 1 & all(modelMatrix == 1) }
    if (justIntercept & all(lambda <= 1e-06)) {
        betaIter <- rep(1, nrow(object))
        betaMatrix <- ... matrix(log2(MatrixGenerics::rowMeans(counts(object, normalized = TRUE))), ncol = 1)
        ...
        xtwx <- MatrixGenerics::rowSums(w); sigma <- xtwx^-1
        ...
        return(res)
    }
    ...
    lambdaNatLogScale <- lambda/log(2)^2
    betaRes <- fitBetaWrapper(ySEXP = counts(object), xSEXP = modelMatrix,
        nfSEXP = normalizationFactors, alpha_hatSEXP = alpha_hat,
        beta_matSEXP = beta_mat, lambdaSEXP = lambdaNatLogScale, ...,
        tolSEXP = betaTol, maxitSEXP = maxit, useQRSEXP = useQR, minmuSEXP = minmu)
    ...
    betaConv <- betaRes$iter < maxit
    betaMatrix <- log2(exp(1)) * betaRes$beta_mat
    betaSE <- log2(exp(1)) * sqrt(pmax(betaRes$beta_var_mat, 0))
    rowsForOptim <- if (useOptim) { which(!betaConv | !rowStable | !rowVarPositive) } else { which(!rowStable | !rowVarPositive) }
    if (forceOptim) { rowsForOptim <- seq_along(betaConv) }
    if (length(rowsForOptim) > 0) { resOptim <- fitNbinomGLMsOptim(...) }
```

- `alpha_hat`의 기본값은 `dispersions(object)`, 곧 [03 노트](03_dispersion_estimation.md)의 최종 α야.
- IRLS 시작값은 `y <- t(log(counts(object, normalized = TRUE) + 0.1)); beta_mat <- t(solve(R, t(Q) %*% y))`, 곧 log(normalized count + 0.1)의 최소제곱해야. 본문 3절과 D의 직접 계산도 같은 시작값을 써.
- log2 변환은 `log2(exp(1)) *`로 해(`log2(exp(1)) == 1/log(2)`는 위 출력에 있어). 계수와 SE가 함께 자연로그에서 log2로 바뀌지. 함수 안에 `type = "glmGamPoi"` 분기(`betaMatrix = gp_res$Beta/log(2), betaSE = NULL`)도 있긴 한데, `fitNbinomGLMs`를 부르는 함수 가운데 `type=`을 넘기는 건 gene-wise dispersion 단계의 `estimateDispersionsGeneEst`뿐이야. `nbinomWaldTest`, `nbinomLRT`, `fitGLMsWithPrior`, `estimateMLEForBetaPriorVar`, `rlogData`는 넘기지 않아. 그래서 보고되는 `log2FoldChange`·`lfcSE`는 위 `log2(exp(1)) *` 경로에서 나와.
- `formals(nbinomWaldTest)$betaPrior`는 `FALSE`이고, `DESeq()` 본문에는 `if (missing(betaPrior)) betaPrior <- FALSE`가 있어. 이때 `attr(dds, "betaPriorVar")`는 `rep(1e+06, ncol)`로 기록되고, 결과 설명은 `log2 fold change (MLE)`가 돼.
- 절편 외 계수가 있는 design의 IRLS에는 ridge $\lambda=10^{-6}$(log2 척도)이 들어가. 자연로그 척도로는 $10^{-6}/\log^2 2$야. `~1` design이고 `lambda <= 1e-06`이면 IRLS를 건너뛰고 log2(normalized count 평균)을 닫힌 형태로 쓰는데(`betaIter = 1`, SE는 $1/\sum_j w_j$에서), 이 경로에는 ridge도 `minmu`도 없어. ridge가 LFC를 $2.1\times10^{-7}$, SE를 $2.9\times10^{-8}$ 옮기는 것과 `~1`의 닫힌 형태는 D에서 확인해.
- 수렴 판정은 `betaRes$iter < maxit`이야(기본 `maxit = 100`, `betaTol = 1e-08`). 불안정한 행(`rowStable`)과 분산이 양수가 아닌 행(`rowVarPositive`)은 항상, 수렴하지 못한 행은 `useOptim = TRUE`(기본)일 때 `fitNbinomGLMsOptim`으로 넘어가(`forceOptim = TRUE`면 모든 행). 이 함수는 log2 척도 계수 $p$에 대해 `mu_row <- nf * 2^(x %*% p)`, `optim(..., method = "L-BFGS-B", lower = -30, upper = 30)`으로 penalized log likelihood(`dnorm(p, 0, sqrt(1/lambda))`)를 최대화하고, `o$convergence == 0`이면 `betaConv[row] <- TRUE`로 바꿔. 끝내 수렴하지 못한 행 수는 `nbinomWaldTest`가 `"rows did not converge in beta, labelled in mcols(object)$betaConv"` 메시지로 알려 줘.
- `minmu` 기본값은 `formals(nbinomWaldTest)$minmu`와 `formals(results)$minmu` 모두 0.5야. `formals(DESeq)$minmu`는 `if (fitType == "glmGamPoi") 1e-06 else 0.5`이고, `DESeq()`는 이 값을 `nbinomWaldTest(..., minmu = minmu)`로 넘겨.

</details>

<details>
<summary>B. IRLS 본체인 C++ fitBeta 발췌</summary>

설치된 패키지에는 컴파일된 `libs/`만 있고 `.cpp`는 없어(`DESeq2:::fitBeta`는 `.Call("_DESeq2_fitBeta", ...)` 래퍼야). 그래서 아래는 Bioconductor 3.22의 소스 tarball `DESeq2_1.50.2.tar.gz`(DESCRIPTION `Version: 1.50.2`; [03 노트](03_dispersion_estimation.md)의 C++ 발췌와 같은 파일)에서 `src/DESeq2.cpp` 안의 `fitBeta`를 가져왔어. 처음에는 GitHub master(`Version: 1.53.4`, devel)에서 가져왔는데, `diff`로 비교해 보니 두 `DESeq2.cpp`는 머리 주석의 저장소 URL 한 줄만 다르고 `fitBeta`를 포함한 나머지는 같았어. 여러 줄에 걸친 `for` 루프와 `if (...) { break; }`는 한 줄로 합쳐 적었어. 이 식이 내는 숫자(계수, SE, contrast SE, `minmu` 효과)는 D와 F에서 설치된 1.50.2로 재현돼.

```cpp
# 발췌 (실행하지 않음): DESeq2_1.50.2 src/DESeq2.cpp, fitBeta
	  w_vec = mu_hat/(1.0 + alpha_hat(i) * mu_hat);
	...
	z = arma::log(mu_hat / nfrow) + (yrow - mu_hat) / mu_hat;
	weighted_x_ridge = join_cols(x.each_col() % w_sqrt_vec, sqrt(ridge));
	qr_econ(q, r, weighted_x_ridge);
	gamma_hat = q.t() * big_z_sqrt_w;
	solve(beta_hat, r, gamma_hat);
	mu_hat = nfrow % exp(x * beta_hat);
	for (int j = 0; j < y_m; j++) { mu_hat(j) = fmax(mu_hat(j), minmu); }
	conv_test = fabs(dev - dev_old)/(fabs(dev) + 0.1);
	if ((t > 0) & (conv_test < tol)) { break; }
    ...
    sigma = (x.t() * (x.each_col() % w_vec) + ridge).i() * x.t() * (x.each_col() % w_vec) * (x.t() * (x.each_col() % w_vec) + ridge).i();
    contrast_num.row(i) = contrast.t() * beta_hat;
    contrast_denom.row(i) = sqrt(contrast.t() * sigma * contrast);
    beta_var_mat.row(i) = diagvec(sigma).t();
```

교재 8.3의 식과는 이렇게 대응해.

- `w_vec`은 $w_j$야.
- `z`는 ($\mu$가 `minmu`로 잘리지 않은 경우) $\eta-o+(K-\mu)/\mu$, 곧 교재의 $z-o$야. $\mu$가 잘리면 `log(mu_hat/nfrow)`가 $x_j^{T}b$와 달라져.
- 갱신은 $[\sqrt W X;\ \sqrt\Lambda]$의 QR 분해로 풀어(교재의 "역행렬 대신 QR").
- 수렴 기준은 deviance의 상대 변화이고, $\mu$는 `minmu`(기본 0.5) 아래로 내려가지 않아.
- `sigma`는 ridge를 포함한 sandwich $(X^{T}WX+\Lambda)^{-1}X^{T}WX(X^{T}WX+\Lambda)^{-1}$이야. $\Lambda\to0$이면 $(X^{T}WX)^{-1}$이 되고.

</details>

<details>
<summary>C. results(contrast=)는 SE를 어떻게 계산할까</summary>

`deparse(body(results))`와 `DESeq2:::cleanContrast`에서 contrast를 처리하는 줄만 모았어.

```r
# 발췌 (실행하지 않음): results 본문과 DESeq2:::cleanContrast
        contrast <- checkContrast(contrast, resNames)
            res <- cleanContrast(object, contrast, expanded = isExpanded, ...
        else if (is.character(contrast)) {
            contrastNumeric <- rep(0, length(resNames))
            contrastNumeric[resNames == contrastNumColumn] <- 1
            contrastNumeric[resNames == contrastDenomColumn] <- -1
            contrast <- contrastNumeric
        contrastResults <- getContrast(object, contrast, useT = useT, minmu)
```

여기서 불리는 `DESeq2:::getContrast`는 이래.

```r
# 발췌 (실행하지 않음): DESeq2:::getContrast
    beta_mat <- log(2) * as.matrix(mcols(objectNZ)[, coefColumns, drop = FALSE])
    lambda = 1/(log(2)^2 * attr(object, "betaPriorVar"))
    betaRes <- fitBeta(ySEXP = countsMatrix, xSEXP = modelMatrix, nfSEXP = normalizationFactors,
        alpha_hatSEXP = alpha_hat, contrastSEXP = contrast, beta_matSEXP = beta_mat,
        lambdaSEXP = lambda, ..., tolSEXP = 1e-08, maxitSEXP = 0, useQRSEXP = FALSE, minmuSEXP = minmu)
    contrastEstimate <- log2(exp(1)) * betaRes$contrast_num
    contrastSE <- log2(exp(1)) * betaRes$contrast_denom
```

`contrast = c("condition", "Glucose", "Starvation")`는 이런 순서로 처리돼.

1. 두 level이 모두 reference가 아니라서 숫자 벡터 $(0,-1,1)$로 바뀌어.
2. 이미 적합된 계수를 `log(2)`배 해서 자연로그로 되돌린 뒤 `maxit = 0`으로 `fitBeta`를 불러. 반복은 하지 않고 B의 `sigma`와 `contrast_denom = sqrt(cᵀ sigma c)`만 계산하는데, 이게 교재 8.6의 $\sqrt{c^{T}\widehat{\mathrm{Cov}}(\hat b)c}$야.
3. 결과를 다시 `log2(exp(1))`배 해서 log2 단위로 돌려줘.

reference가 들어간 contrast(`c("condition", "Starvation", "Ctrl")`)는 `getCoef`/`getCoefSE`로 저장된 열을 그대로 읽어. 분자가 reference이면(예: `c("condition", "Ctrl", "Starvation")`) `swapName`으로 저장된 `condition_Starvation_vs_Ctrl` 열을 읽고 LFC와 stat에 −1을 곱해. SE는 그대로고.

`coef`는 `DESeq2:::coef.DESeqDataSet`으로, `mcols(object)[resultsNames]`를 행렬로 돌려줘. `coef(dds, SE = TRUE)`는 `paste0("SE_", resNms)` 열을 돌려주고, 둘 다 log2 단위야.

</details>

<details>
<summary>D. 유전자 A의 IRLS 전 과정과, ridge·수렴 기준이 만드는 아주 작은 차이</summary>

본문 3절의 반복을 여덟 번까지 늘리고 deviance도 함께 적었어. 본문 1·3절의 `k`, `X`, `o`, `alpha`를 이어 써.

```r
options(digits = 9, width = 200)             # 본문 1·3절의 k, X, o, alpha를 이어 쓴다
b <- qr.solve(X, log(k + 0.1))
tab <- NULL
for (t in 0:7) {
  eta <- as.numeric(o + X %*% b); mu <- exp(eta)
  z <- eta + (k - mu) / mu
  w <- mu / (1 + alpha * mu)
  U <- as.numeric(crossprod(X, (k - mu) / (1 + alpha * mu)))
  dev <- -2 * sum(dnbinom(k, mu = mu, size = 1 / alpha, log = TRUE))
  tab <- rbind(tab, data.frame(iter = t, b0 = b[[1]], bS = b[[2]], exp_bS = exp(b[[2]]), beta_S_log2 = b[[2]] / log(2),
                               score_bS = U[2], deviance = dev))
  b <- as.numeric(solve(crossprod(X, X * w), crossprod(X, w * (z - o))))
}
print(tab, row.names = FALSE)
```

```
 iter         b0          bS     exp_bS beta_S_log2        score_bS   deviance
    0 4.65846441 0.679599032 1.97308643 0.980454153  4.70307252e-01 56.3747900
    1 4.66977216 0.677376390 1.96870583 0.977247559 -2.12509410e-03 56.3644595
    2 4.66970871 0.677398822 1.96875000 0.977279922 -4.35871990e-08 56.3644593
    3 4.66970871 0.677398824 1.96875000 0.977279923  1.40998324e-14 56.3644593
    4 4.66970871 0.677398824 1.96875000 0.977279923 -3.53050922e-14 56.3644593
    5 4.66970871 0.677398824 1.96875000 0.977279923 -3.53050922e-14 56.3644593
    6 4.66970871 0.677398824 1.96875000 0.977279923 -3.53050922e-14 56.3644593
    7 4.66970871 0.677398824 1.96875000 0.977279923  1.40998324e-14 56.3644593
```

세 번째 반복 뒤 score가 $10^{-14}$ 수준으로 떨어지고, 그 뒤로는 움직이지 않아. $\exp(b_0)=106.666667$, $\exp(b_S)=1.96875=210/106.666667$이고, $\beta_S=b_S/\log 2=0.977279923$으로 교재의 0.977280과 같아.

다음은 ridge를 넣은 SE와 DESeq2가 저장한 값, 그리고 둘의 차이야. 본문 3절의 `dds`, `res`를 이어 써.

```r
eta <- as.numeric(o + X %*% b); mu <- exp(eta); w <- mu / (1 + alpha * mu)
XtWX <- crossprod(X, X * w); Cov_nat <- solve(XtWX)
cat("weights w     =", w, "\n"); print(XtWX); print(Cov_nat)
cat("SE(b) natural =", sqrt(diag(Cov_nat)), "\nSE(beta) log2 =", sqrt(diag(Cov_nat)) / log(2), "\n")
L <- diag(rep(1e-6 / log(2)^2, 2)); S <- solve(XtWX + L) %*% XtWX %*% solve(XtWX + L)
cat("with ridge:  SE(beta) log2 =", sqrt(diag(S)) / log(2), "\n")
print(as.data.frame(mcols(dds)[, c("Intercept", "condition_Starvation_vs_Ctrl", "SE_Intercept", "SE_condition_Starvation_vs_Ctrl", "betaConv", "betaIter")]))
cat("mu assay:", assays(dds)[["mu"]], "\n")
print(as.data.frame(res))
cat("manual - DESeq2: LFC diff =", b[2] / log(2) - res$log2FoldChange, "  SE diff =", sqrt(Cov_nat[2, 2]) / log(2) - res$lfcSE, "\n")
cat("with-ridge manual SE - DESeq2 lfcSE =", sqrt(S[2, 2]) / log(2) - res$lfcSE, "\n")
```

```
weights w     = 15.99437 15.99437 15.99437 17.2685013 17.2685013 17.2685013
                    (Intercept) conditionStarvation
(Intercept)           99.788614           51.805504
conditionStarvation   51.805504           51.805504
                      (Intercept) conditionStarvation
(Intercept)          0.0208406667       -0.0208406667
conditionStarvation -0.0208406667        0.0401436349
SE(b) natural = 0.144362968 0.200358766
SE(beta) log2 = 0.208271739 0.289056597
with ridge:  SE(beta) log2 = 0.208271721 0.289056567
       Intercept condition_Starvation_vs_Ctrl SE_Intercept SE_condition_Starvation_vs_Ctrl betaConv betaIter
geneA 6.73696535                  0.977280132  0.208271723                     0.289056569     TRUE        2
mu assay: 106.666648 106.666648 106.666648 209.999994 209.999994 209.999994
        baseMean log2FoldChange       lfcSE       stat         pvalue           padj
geneA 158.333333    0.977280132 0.289056569 3.38093037 0.000722408466 0.000722408466
manual - DESeq2: LFC diff = -2.08896144e-07   SE diff = 2.85845623e-08
with-ridge manual SE - DESeq2 lfcSE = -2.07655015e-09
```

$SE(\beta_S)=0.289056597$로 교재의 0.289057이 재현돼. Ctrl 한 관측의 weight는 $106.67/(1+0.053147\times106.67)=15.99$, Starvation은 17.27이라, $\mathrm{Var}(\hat b_S)=1/51.81+1/47.98=0.040144$라는 닫힌 형태와도 맞아. DESeq2의 ridge를 넣으면 SE가 $3\times10^{-8}$ 줄어들어.

DESeq2의 LFC는 직접 계산한 MLE보다 $2.1\times10^{-7}$ 커. 이 차이가 2회 만에 멈춘 탓인지 ridge 탓인지 가르려고, `fitBeta`처럼 $(X^{T}WX+\Lambda)$로 푸는 penalized IRLS를 끝까지 돌리고 DESeq2의 `betaTol`도 조여 봤어. 앞 블록의 `b`(MLE), `L`, `XtWX`, `Cov_nat`, `res`를 이어 써.

```r
bR <- b                                                  # b: 위 IRLS 의 MLE
for (t in 1:30) {                                        # fitBeta 처럼 (X'WX + Lambda) 로 푸는 penalized IRLS
  eta <- as.numeric(o + X %*% bR); mu <- exp(eta); w <- mu / (1 + alpha * mu); z <- eta + (k - mu) / mu
  bR <- as.numeric(solve(crossprod(X, X * w) + L, crossprod(X, w * (z - o))))
}
mu <- exp(as.numeric(o + X %*% bR)); w <- mu / (1 + alpha * mu); A <- crossprod(X, X * w)
SR <- solve(A + L) %*% A %*% solve(A + L)
cat("ridge IRLS: LFC =", bR[2] / log(2), " lfcSE =", sqrt(SR[2, 2]) / log(2), " mu =", mu[c(1, 4)], "\n")
cat("ridge IRLS - DESeq2: LFC", bR[2] / log(2) - res$log2FoldChange, "  SE", sqrt(SR[2, 2]) / log(2) - res$lfcSE, "\n")
cat("ridge effect (ridge - MLE): LFC", (bR[2] - b[2]) / log(2), "  SE", (sqrt(SR[2, 2]) - sqrt(Cov_nat[2, 2])) / log(2),
    "  1st-order -(X'WX)^-1 L b:", -solve(XtWX, L %*% b)[2] / log(2), "\n")
dT <- DESeqDataSetFromMatrix(cts1, cd, ~ condition); sizeFactors(dT) <- rep(1, 6); dispersions(dT) <- alpha
dT <- nbinomWaldTest(dT, betaTol = 1e-14, maxit = 1000, quiet = TRUE)
cat("DESeq2 betaTol=1e-14: LFC =", mcols(dT)$condition_Starvation_vs_Ctrl, " betaIter =", mcols(dT)$betaIter,
    " (vs default betaTol:", mcols(dT)$condition_Starvation_vs_Ctrl - res$log2FoldChange, ")\n")
dI <- DESeqDataSetFromMatrix(cts1, cd, ~ 1); sizeFactors(dI) <- rep(1, 6); dispersions(dI) <- alpha
dI <- nbinomWaldTest(dI, quiet = TRUE)
cat("~1: Intercept =", mcols(dI)$Intercept, " log2(mean(k)) =", log2(mean(k)), " betaIter =", mcols(dI)$betaIter, "\n")
```

```
ridge IRLS: LFC = 0.977280134  lfcSE = 0.289056569  mu = 106.666648 209.999994
ridge IRLS - DESeq2: LFC 1.6785906e-09   SE 2.73474021e-11
ridge effect (ridge - MLE): LFC 2.10574735e-07   SE -2.85572149e-08   1st-order -(X'WX)^-1 L b: 2.10574776e-07
DESeq2 betaTol=1e-14: LFC = 0.977280134  betaIter = 3  (vs default betaTol: 1.67858871e-09 )
~1: Intercept = 7.3068212  log2(mean(k)) = 7.3068212  betaIter = 1
```

| 항목 | 직접 IRLS (MLE) | 직접 IRLS + ridge | DESeq2 1.50.2 | 교재 16.4 |
|---|---|---|---|---|
| $\beta_S$ (log2FoldChange) | 0.977279923 | 0.977280134 | 0.977280132 | 0.977280 |
| $SE(\beta_S)$ (lfcSE) | 0.289056597 | 0.289056569 | 0.289056569 | 0.289057 |
| Ctrl `mu` | 106.666667 | 106.666648 | 106.666648 | 106.666667 |
| 반복 수 | score $<10^{-7}$까지 2회 | 30회(수렴 확인용) | `betaIter` = 2 | — |

결국 차이는 ridge에서 나와.

- 벌점 $\tfrac12 b^{T}\Lambda b$가 큰 절편 $b_0\approx4.67$을 조금 내려(Ctrl `mu` 106.666648, 산술평균과 $2\times10^{-5}$ 차이). 그러면 Starvation 평균을 맞추려고 $b_S$가 거의 그만큼 올라가.
- 그래서 ridge는 LFC를 $+2.1\times10^{-7}$, SE를 $-2.9\times10^{-8}$ 옮기고, 1차 근사 $-(X^{T}WX)^{-1}\Lambda b$가 LFC 이동을 그대로 설명해.
- ridge를 넣은 IRLS는 DESeq2와 LFC $1.7\times10^{-9}$, SE $3\times10^{-11}$ 안에서 같아.
- `betaTol = 1e-14`로 조이면 DESeq2도 3회 만에 ridge 해 0.977280134에 닿아. 기본 `betaTol = 1e-08`에서 2회 만에 멈춘 영향은 $1.7\times10^{-9}$뿐이야.
- 앞 블록의 `with-ridge manual SE - DESeq2 lfcSE` $=-2\times10^{-9}$는 ridge sandwich를 ridge 해가 아닌 MLE의 $\mu$에서 계산해서 남은 잔차야.
- 절편만 있는 `~1` design은 IRLS 없이 $\log_2(\bar K)$를 그대로 쓰고(`betaIter` = 1), ridge도 들어가지 않아(A).

</details>

<details>
<summary>E. α가 바뀌면 LFC도 움직이는 경우</summary>

본문 4절의 표는 평균을 고정하고 α만 바꿨어. 여기서는 DESeq2로 α마다 다시 적합해서, size factor가 모두 1인 경우와 sample마다 다른 경우(0.8, 1, 1.25, 0.9, 1.1, 1)를 비교해. 본문 3절의 `cts1`, `cd`를 이어 써.

```r
for (a in c(0.014786, 0.025385, 0.053147, 0.5)) {        # DESeq2 로 alpha 마다 재적합
  dA <- DESeqDataSetFromMatrix(cts1, cd, ~ condition); dispersions(dA) <- a
  sizeFactors(dA) <- rep(1, 6);                    dA1 <- nbinomWaldTest(dA, quiet = TRUE)
  sizeFactors(dA) <- c(0.8, 1, 1.25, 0.9, 1.1, 1); dA2 <- nbinomWaldTest(dA, quiet = TRUE)
  cat("alpha =", a, " LFC (s = 1) =", mcols(dA1)$condition_Starvation_vs_Ctrl, "  LFC (unequal s) =", mcols(dA2)$condition_Starvation_vs_Ctrl, "\n")
}
```

```
alpha = 0.014786  LFC (s = 1) = 0.977280007   LFC (unequal s) = 0.965107491
alpha = 0.025385  LFC (s = 1) = 0.97728004   LFC (unequal s) = 0.958733077
alpha = 0.053147  LFC (s = 1) = 0.977280132   LFC (unequal s) = 0.952472028
alpha = 0.5  LFC (s = 1) = 0.977281615   LFC (unequal s) = 0.945818515
```

- `s = 1`처럼 두 그룹이고 size factor가 모두 같으면 그룹별 score 식이 $\sum_{j\in g}(K_j-\mu_g)/(1+\alpha\mu_g)=0$이 돼. 그래서 어떤 α에서도 $\hat\mu_g$가 산술평균이고, LFC는 ridge 영향($<2\times10^{-6}$) 안에서 같아.
- `unequal s`에서는 같은 그룹 안에서도 weight가 달라져서, α가 0.0148 → 0.5로 갈 때 LFC가 0.965 → 0.946으로 움직여.
- 계수가 그룹 평균으로 바로 풀리지 않는 여러 요인의 additive design도 같은 이유로 α에 따라 $\hat\beta$가 바뀌어. 이때 결과가 한 방향으로 단조롭게 변한다고 가정하면 안 돼(교재 8.4, 16.5).

</details>

<details>
<summary>F. 세 그룹 예제 전체에서 본 수렴, 공분산, minmu, all-zero contrast</summary>

먼저 본문 6절의 `dds3`에서 기본 설정과 수렴 상태를 봐.

```r
options(digits = 7)                          # 본문 6절의 dds3를 이어 쓴다
print(resultsNames(dds3))
cat("attr betaPrior =", attr(dds3, "betaPrior"), " modelMatrixType =", attr(dds3, "modelMatrixType"), "\n")
print(head(coef(dds3), 3)); print(head(coef(dds3, SE = TRUE), 3))
print(table(mcols(dds3)$betaConv, useNA = "ifany")); print(summary(mcols(dds3)$betaIter))
```

```
[1] "Intercept"                    "condition_Starvation_vs_Ctrl" "condition_Glucose_vs_Ctrl"
attr betaPrior = FALSE  modelMatrixType = standard
      Intercept condition_Starvation_vs_Ctrl condition_Glucose_vs_Ctrl
gene1  1.435842                    0.2414319                 0.6950322
gene2  4.736040                   -0.4925971                -2.5054179
gene3  1.239198                   -2.7912796                 3.4807676
      SE_Intercept SE_condition_Starvation_vs_Ctrl SE_condition_Glucose_vs_Ctrl
gene1    1.0202100                       1.4243343                    1.4036809
gene2    0.5629422                       0.8014664                    0.8635503
gene3    0.8535847                       1.5960562                    1.0883241

TRUE <NA>
1986   14
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
  2.000   3.000   4.000   4.127   5.000  15.000      14
```

`NA` 14행은 모든 count가 0인 유전자이고, 나머지 1986행은 2–15회 안에 모두 수렴했어.

이어서 rescue contrast를 세 가지 표기로 구하고, 계수의 공분산을 $(X^{T}WX)^{-1}$로 직접 계산해 모든 유전자에서 비교해. `mu`는 `assays(dds3)[["mu"]]`, α는 `dispersions(dds3)`, X는 `attr(dds3, "modelMatrix")`에서 가져왔어. `fitBeta`의 `fmax(mu_hat(j), minmu)`를 재현하려고 `pmax(mu, 0.5)`를 넣은 경우와 넣지 않은 경우를 모두 계산했고, 본문 6절의 `rS`, `rG`, `rGS`를 이어 써.

```r
rGS_list <- results(dds3, contrast = list("condition_Glucose_vs_Ctrl", "condition_Starvation_vs_Ctrl"))
rGS_num  <- results(dds3, contrast = c(0, -1, 1))
cat("same LFC/lfcSE for the 3 spellings:", all.equal(rGS$log2FoldChange, rGS_list$log2FoldChange), all.equal(rGS$lfcSE, rGS_num$lfcSE), "\n")
X <- attr(dds3, "modelMatrix"); MU <- assays(dds3)[["mu"]]; A <- dispersions(dds3); lam <- 1e-6 / log(2)^2
man <- function(clamp) t(sapply(seq_len(nrow(dds3)), function(g) {
  mu <- MU[g, ]; if (any(is.na(mu))) return(rep(NA, 4))
  if (clamp) mu <- pmax(mu, 0.5)                     # fitBeta: mu_hat(j) = fmax(mu_hat(j), minmu)
  w <- mu / (1 + A[g] * mu); XtWX <- crossprod(X, X * w)
  Sg <- solve(XtWX + diag(lam, 3)) %*% XtWX %*% solve(XtWX + diag(lam, 3)); C2 <- Sg / log(2)^2
  c(sqrt(C2[2, 2]), sqrt(C2[3, 3]), C2[2, 3], sqrt(C2[2, 2] + C2[3, 3] - 2 * C2[2, 3]))
}))
m0 <- man(FALSE); m1 <- man(TRUE); ok <- !mcols(dds3)$allZero
cat("genes with any mu < 0.5:", sum(apply(MU < 0.5, 1, any), na.rm = TRUE), "\n")
cat("[no clamp]  max|SE_diff - lfcSE(GvsS)| =", max(abs(m0[ok, 4] - rGS$lfcSE[ok])), "\n")
cat("[clamp 0.5] max|SE_S - lfcSE(S)| =", max(abs(m1[ok, 1] - rS$lfcSE[ok])), " max|SE_G - lfcSE(G)| =", max(abs(m1[ok, 2] - rG$lfcSE[ok])),
    " max|SE_diff - lfcSE(GvsS)| =", max(abs(m1[ok, 4] - rGS$lfcSE[ok])), "\n")
zeroGS <- rowSums(counts(dds3)[, 4:9]) == 0
cat("all-zero in Starvation+Glucose samples:", sum(zeroGS), " -> rGS LFC/stat/pvalue:", unique(rGS$log2FoldChange[zeroGS]), unique(rGS$stat[zeroGS]), unique(rGS$pvalue[zeroGS]), "\n")
cat("others: max|LFC(GvsS) - (LFC_G - LFC_S)| =", max(abs(rGS$log2FoldChange - (rG$log2FoldChange - rS$log2FoldChange))[ok & !zeroGS]), "\n")
cat("max|sqrt(SE_S^2+SE_G^2) - lfcSE(GvsS)| =", max(abs(sqrt(rS$lfcSE^2 + rG$lfcSE^2) - rGS$lfcSE)[ok]), "  min Cov_GS =", min(m1[ok, 3]), "\n")
g <- which(ok & rS$baseMean > 50)[1:4]
print(data.frame(gene = rownames(dds3)[g], baseMean = rS$baseMean[g], alpha = A[g],
                 LFC_S = rS$log2FoldChange[g], LFC_G = rG$log2FoldChange[g], LFC_GvsS = rGS$log2FoldChange[g],
                 SE_S = rS$lfcSE[g], SE_G = rG$lfcSE[g], Cov_GS = m1[g, 3],
                 SE_sum = rS$lfcSE[g] + rG$lfcSE[g], SE_sqrtsumsq = sqrt(rS$lfcSE[g]^2 + rG$lfcSE[g]^2),
                 SE_with_cov = m1[g, 4], lfcSE_results = rGS$lfcSE[g]), digits = 4, row.names = FALSE)
g1 <- g[1]; mu <- MU[g1, ]; w <- mu / (1 + A[g1] * mu)
cat("gene4: Cov_GS =", m1[g1, 3], " 1/sum(w_Ctrl)/log(2)^2 =", 1 / sum(w[1:3]) / log(2)^2, "\n")
cat("gene4: sqrt(1/sum(w_Starvation)+1/sum(w_Glucose))/log(2) =", sqrt(1 / sum(w[4:6]) + 1 / sum(w[7:9])) / log(2), " lfcSE(GvsS) =", rGS$lfcSE[g1], "\n")
```

```
same LFC/lfcSE for the 3 spellings: TRUE TRUE
genes with any mu < 0.5: 184
[no clamp]  max|SE_diff - lfcSE(GvsS)| = 0.6087883
[clamp 0.5] max|SE_S - lfcSE(S)| = 4.440892e-15  max|SE_G - lfcSE(G)| = 3.552714e-15  max|SE_diff - lfcSE(GvsS)| = 3.552714e-15
all-zero in Starvation+Glucose samples: 22  -> rGS LFC/stat/pvalue: NA 0 NA 0 NA 1
others: max|LFC(GvsS) - (LFC_G - LFC_S)| = 8.881784e-16
max|sqrt(SE_S^2+SE_G^2) - lfcSE(GvsS)| = 1.729878   min Cov_GS = 0.06837815
   gene baseMean  alpha    LFC_S   LFC_G  LFC_GvsS   SE_S   SE_G Cov_GS SE_sum SE_sqrtsumsq SE_with_cov lfcSE_results
  gene4   178.54 0.2061  0.07263 0.06725 -0.005378 0.5419 0.5419 0.1471  1.084       0.7664      0.5415        0.5415
 gene11    96.67 0.1698  0.87821 0.54089 -0.337325 0.5014 0.5026 0.1282  1.004       0.7099      0.4976        0.4976
 gene15   103.86 0.2177 -0.13691 0.51642  0.653321 0.5634 0.5608 0.1586  1.124       0.7950      0.5610        0.5610
 gene20    54.78 0.3239  0.29618 1.37367  1.077493 0.6970 0.6909 0.2453  1.388       0.9814      0.6874        0.6874
gene4: Cov_GS = 0.147065  1/sum(w_Ctrl)/log(2)^2 = 0.1470652
gene4: sqrt(1/sum(w_Starvation)+1/sum(w_Glucose))/log(2) = 0.5415034  lfcSE(GvsS) = 0.5415033
```

- `c("condition", "Glucose", "Starvation")`, `list(...)`, `c(0, -1, 1)` 세 표기의 LFC와 lfcSE가 같아. LFC는 정확히 $\hat\beta_G-\hat\beta_S$야(차이 $9\times10^{-16}$).
- 예외는 비교 대상인 Starvation과 Glucose sample 여섯 개가 모두 0인 22개 유전자야. `cleanContrast`의 `contrastAllZero` 처리 때문에 LFC 0, stat 0, p 1로 고정되고, 그중 전체가 0인 행은 NA야.
- 계수의 SE와 contrast의 SE는 모두 log 2로 나눈 $(X^{T}WX+\Lambda)^{-1}X^{T}WX(X^{T}WX+\Lambda)^{-1}$에서 나와. `pmax(mu, 0.5)`를 넣으면 1986개 유전자에서 $4\times10^{-15}$ 안에서 맞고, 넣지 않으면 $\mu<0.5$인 184개 유전자에서 최대 0.61 어긋나. `minmu`가 SE 계산에도 실제로 들어간다는 뜻이야.
- gene4에서 $SE_S=SE_G=0.5419$인데 Glucose−Starvation의 SE는 0.5415야. $\mathrm{Cov}(\hat\beta_G,\hat\beta_S)=0.1471=\mathrm{Var}(\hat\beta_0)=1/(\sum_{\mathrm{Ctrl}}w_j\log^2 2)$로 양수라서, rescue SE는 $\sqrt{1/\sum_S w_j+1/\sum_G w_j}/\log 2=0.5415034$로 Ctrl의 weight와 무관해지고, `lfcSE` 0.5415033과 ridge 몫만큼만 달라.
- 1986개 유전자 모두에서 $\mathrm{Cov}>0$이었어(최소 0.068). `sqrt(SE_S^2 + SE_G^2)`로 계산하면 `lfcSE`와 최대 1.73까지 어긋나.

</details>

<details>
<summary>G. 95% CI, padj, 교재 Z와의 반올림 차이</summary>

먼저 padj야. 본문 3절에서 만든 `res`의 `padj`(D의 출력)가 `pvalue`와 같은 0.000722408466인 건, 유전자가 하나뿐이라 BH 보정의 대상이 자기 자신뿐이기 때문이야. 의미 있는 padj가 아니지. 교재 16.4도 단일 유전자에서는 padj를 제시하지 않아. 실제 padj는 다른 유전자들의 p와 보정 대상 집합이 있어야 계산돼([06 노트](06_multiple_testing.md)).

다음은 95% CI야. 교재 16.4의 nominal 95% CI는 `results()`가 기본으로 출력하지 않아서 직접 계산했어.

```r
options(digits = 9)                          # 본문 3절의 res를 이어 쓴다
cat("95% CI, DESeq2 LFC ± qnorm(0.975)*lfcSE:", res$log2FoldChange + c(-1, 1) * qnorm(0.975) * res$lfcSE, "\n")
cat("rounded 0.977280 ± z*0.289057: z = qnorm(0.975) ->", 0.977280 + c(-1, 1) * qnorm(0.975) * 0.289057,
    "  z = 1.96 ->", 0.977280 + c(-1, 1) * 1.96 * 0.289057, "\n")
cat("z implied by textbook CI [0.410727, 1.543833]:", (1.543833 - 0.410727) / 2 / 0.289057, "\n")
```

```
95% CI, DESeq2 LFC ± qnorm(0.975)*lfcSE: 0.410739668 1.5438206
rounded 0.977280 ± z*0.289057: z = qnorm(0.975) -> 0.410738691 1.54382131   z = 1.96 -> 0.41072828 1.54383172
z implied by textbook CI [0.410727, 1.543833]: 1.96000443
```

- DESeq2 값으로는 [0.410740, 1.543821]이고, 교재의 반올림 값에 `qnorm(0.975)`를 써도 [0.410739, 1.543821]이야.
- 교재의 [0.410727, 1.543833]은 $z=1.96$을 쓴 값 [0.410728, 1.543832]에 가까워. 교재 CI에서 역산한 $z$가 1.960004라서, 교재는 $z=1.96$을 쓴 것으로 보여.
- 남는 마지막 자리 차이는 반올림한 입력(교재 표시값 0.977280, 0.289057, 또는 이 노트가 고정한 6자리 반올림 α 0.053147)에서 와.
- [08 노트](08_one_gene_end_to_end.md)의 16.4 부분처럼 수렴한 $\alpha_{MAP}=0.053147342$를 쓰면 $Z=3.3809197$, $p=0.00072244$, $z=1.96$ CI $[0.410727, 1.543832]$로 교재와 반올림 수준에서 맞아(CI 상한 마지막 자리 1 차이). CI를 어떻게 해석하는지는 [05 노트](05_wald_vs_lrt.md)에 있어.

마지막으로 Z야. 본문 3절의 `stat` 3.38093037은 교재 16.4의 $Z=3.380919$와 여섯째 유효숫자에서 달라. 이 노트가 α를 반올림값 0.053147로 고정했기 때문이야. 자세한 비교는 [08 노트](08_one_gene_end_to_end.md)에서, Wald $Z$와 $p$의 일반 구조는 [05 노트](05_wald_vs_lrt.md)에서 다뤄.

</details>

<details>
<summary>H. 자주 헷갈리는 것들</summary>

- 자연로그 식의 $b_S$가 곧 log2FC일까? 아니야. `log2FoldChange` $=b_S/\log 2$이고, DESeq2는 `log2(exp(1)) *`로 변환해(A).
- Glucose 군은 보통 배지에 glucose를 더 넣은 조건일까? 이 노트의 Glucose는 Starvation+Glucose, 곧 굶긴 뒤 glucose를 다시 넣은 조건이야. rescue는 $\beta_G-\beta_S$이고, $\beta_G$는 Ctrl 대비야.
- Ctrl–Starvation, Starvation–Glucose를 두 dataset으로 따로 돌리는 것과 세 그룹 한 모형은 같을까? rescue SE의 식 자체는 같아. `~ condition`에서 $-2\,\mathrm{Cov}$가 Ctrl 항을 지워 $\mathrm{Var}(\hat b_G-\hat b_S)=1/\sum_G w_j+1/\sum_S w_j$가 되는데, 이건 같은 $w$로 Starvation vs Glucose만 적합한 식과 같거든(본문 6절, F). $-2\,\mathrm{Cov}$는 Ctrl 기준 계수에서 rescue SE를 바르게 읽는 방법이지 한 모형의 이점이 아니야. 실제 차이는 $\alpha_i$와 dispersion trend를 더 많은 sample·residual df로 추정하는지, size factor를 모든 sample로 함께 구하는지에 있어.
- 차이의 SE는 $SE_G+SE_S$나 $\sqrt{SE_G^2+SE_S^2}$일까? 전자는 상한이라 $\rho=-1$일 때만 같고, 후자는 $\mathrm{Cov}=0$일 때만 맞아. 요인이 하나인 `~ condition`에서 같은 reference를 공유하는 두 level은 $\mathrm{Cov}=\mathrm{Var}(\hat b_0)>0$이라 실제 SE가 더 작아. 서로 다른 요인의 계수 사이 공분산은 음수일 수도 있고(문제 12).
- `~ pair + condition`은 donor를 random effect로 다룰까? 아니야. donor dummy 열을 넣은 고정효과 additive 모형이고, design은 `stats::model.matrix.default`로만 펼쳐져(A).
- 독립 sample이라도 순서가 같으면 pair ID를 붙여도 될까? 실제 짝이 없는데 pair ID를 붙이면 데이터에 없는 구조를 모형에 강제하게 돼(교재 4.4).
- 기술 반복, 시간 반복도 `~ pair + condition`에 그대로 넣으면 될까? 의존 구조가 복잡하면 이 단순한 design으로 충분한지 따로 판단해야 해. 같은 library를 여러 번 sequencing한 기술 반복은 `collapseReplicates()`로 count를 합쳐 한 sample로 만들어. 반복 측정 사이의 상관이나 donor별 random effect가 필요하면 DESeq2의 고정효과 design으로는 표현되지 않아.
- `~ batch + condition`이라고 쓰면 batch 효과가 보정될까? batch와 condition이 완전히 겹치면 rank deficiency로 `checkFullRank`에서 멈춰. 정보가 부족한 거지 옵션 문제가 아니야(본문 2절).
- 깊게 읽으면 SE는 계속 줄어들까? 관측 하나의 $w_j<1/\alpha$이고, $\mu\to\infty$에서 $1/\alpha$로 수렴해. α = 0.1이면 μ = 100에서 이미 상한의 91%야(본문 5절).
- `betaPrior=FALSE` 결과는 순수 MLE일까? 실무적으로 MLE라 부르지만, 절편 외 계수가 있는 design의 IRLS에는 $\lambda=10^{-6}$ ridge가 있고, `minmu`(기본 fitType에서 0.5, glmGamPoi면 $10^{-6}$)와 optim 경계 ±30도 있어. 유전자 A에서 ridge는 LFC를 $2.1\times10^{-7}$, SE를 $2.9\times10^{-8}$ 옮겨. `~1` design은 ridge 없는 닫힌 형태야(A, D).
- `results(contrast=)`는 저장된 두 LFC를 빼고 SE를 조합할까? LFC는 빼지만, SE는 `fitBeta(maxit = 0)`로 $c^{T}\Sigma c$를 다시 계산해(C).

</details>

<details>
<summary>I. 교재와 다르게 나온 부분</summary>

로드맵의 `04_glm_condition_batch.R`에 해당하는 확인이고, 2026-09-26에 했어. 괄호 속 글자는 이 노트의 "더 깊이 보기" 블록이야. 아래 표에는 교재와 다르거나 교재에 없던 내용이 나온 항목만 남겼어. 나머지 교재 설명, 예를 들어 design matrix와 rank 검사, `betaPrior` 기본값, 아주 약한 ridge, IRLS와 SE 식, contrast SE, 16.4의 LFC와 SE는 실제 DESeq2 동작과 맞았어.

| 교재 주장 | 확인 결과 | 근거 |
|---|---|---|
| 4.2 두 그룹 design matrix (Intercept, Starvation) | 일치(표기 차이) | `model.matrix(~condition)` 열 이름은 `conditionStarvation`, DESeq2는 `renameModelMatrixColumns`로 `condition_Starvation_vs_Ctrl` (본문 1절, A) |
| 8.3 IRLS 식 ($\eta,\mu,z,w$, 가중 최소제곱) | 일치(형태 차이) | C++ `z = log(mu_hat/nfrow) + (yrow-mu_hat)/mu_hat`은 ($\mu$가 `minmu`로 잘리지 않으면) 교재의 $z-o$; `w_vec = mu_hat/(1+alpha*mu_hat)` (B) |
| 16.4 $Z=3.380919$, $p\approx0.00072244$ | 일치 (이 노트의 α 입력은 반올림값) | 이 노트는 α를 반올림값 0.053147로 고정해서 $Z=3.380929$ (DESeq2 `stat` 3.38093037, `pvalue` 0.000722408)가 나와. 교재 값은 수렴한 $\alpha_{MAP}=0.053147342$에서 나온 거야. [08 노트](08_one_gene_end_to_end.md) 16.4 부분의 손계산은 $Z=3.3809197$, $p=0.00072244$이고, 같은 α로 고정한 DESeq2 `stat`은 3.380921이야. $0.977280/0.289057=3.380925$는 반올림한 SE로 나눈 값이야 |
| 16.1 "DESeq2를 실제로 실행한 출력이 아니다" | 확인 보강 | `dispersions(dds) <- 0.053147`로 고정하면 DESeq2도 같은 LFC/SE를 내 (본문 3절) |
| (교재 미기재) `minmu`(기본 fitType에서 0.5)가 W와 SE에 들어가 | 추가 확인 | `formals(DESeq)$minmu`는 `if (fitType == "glmGamPoi") 1e-06 else 0.5`. `pmax(mu,0.5)` 없이는 184/1986 유전자에서 SE 최대 0.61 차이, 넣으면 $4\times10^{-15}$ (F) |
| (교재 미기재) 대상 sample이 전부 0인 contrast | 추가 확인 | `cleanContrast`의 `contrastAllZero`가 LFC 0, stat 0, p 1로 고정 (F) |
| 16.4 nominal 95% CI $\approx[0.410727, 1.543833]$ | 일치 (교재는 $z=1.96$ 사용; 상한 마지막 자리 1 차이) | 반올림 α(0.053147)의 DESeq2 값에 `qnorm(0.975)`면 [0.410740, 1.543821], 교재 반올림 값에 $z=1.96$이면 [0.410728, 1.543832], 교재 CI에서 역산한 $z$ = 1.960004 (G). 수렴 $\alpha_{MAP}$에서는 $z=1.96$ CI [0.410727, 1.543832] ([08 노트](08_one_gene_end_to_end.md) 16.4 부분). `results()`는 CI를 기본으로 출력하지 않아 |

</details>

<details>
<summary>J. 이 노트의 숫자가 이어지는 곳</summary>

- [05 노트](05_wald_vs_lrt.md)에서는 여기서 얻은 $\hat\beta$와 SE로 $Z=\hat\beta/SE$와 p를 만들고 LRT와 비교해. 본문 3절의 `stat` 3.38093037이 출발점이야.
- [06 노트](06_multiple_testing.md)는 유전자 수천 개의 p를 BH로 보정해.
- [07 노트](07_lfc_shrinkage_and_qc.md)는 `betaPrior=FALSE`의 MLE 계수에 사후적으로 shrinkage를 거는 `lfcShrink`와, 그 결과의 `lfcSE`가 type에 따라 무엇인지 다뤄(apeglm/ashr는 posterior SD, normal은 MAP의 sandwich SE).
- [08 노트](08_one_gene_end_to_end.md)는 16장 전체를 한 번에 이어.

</details>

---

← 이전: [03. dispersion 추정](03_dispersion_estimation.md) · 다음: [05. Wald 검정과 LRT](05_wald_vs_lrt.md) →

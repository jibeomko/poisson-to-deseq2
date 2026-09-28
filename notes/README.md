# DESeq2 공부 노트

DESeq2를 돌리면 결과표에 log2 fold change, lfcSE, pvalue, padj가 나온다. 그런데 이 값들이 raw count에서 어떻게 계산되는지, 어떤 수리통계적 과정을 거치는지 궁금했다. 그래서 개인 학습 노트(「DESeq2 학습 노트 · COUNT에서 결론까지」)를 장 순서대로 따라가며 하나씩 직접 계산해 봤고, 이 폴더는 그 기록이다. 노트에 나온 설명과 숫자는 설치된 DESeq2로 다시 돌려서 확인했다. 
평균, 분산, p 값은 알지만 likelihood나 GLM은 이름만 들어 본 사람을 떠올리며 썼다.

## 읽는 순서

처음이라면 00_overview.md부터 순서대로 읽는 것을 추천한다. 노트마다 질문 하나로 시작해서 숫자로 먼저 보여 주고, 식은 그다음에 나온다. 소스 코드를 뜯어본 내용이나 세부 확인 과정은 노트 끝 "더 깊이 보기"에 접어 두었으니 처음 읽을 때는 건너뛰어도 괜찮다.

모든 노트가 같은 예제를 쓴다. 교재 16장에 나오는 유전자 A인데, 세포를 보통 배지(Ctrl)와 glucose를 뺀 배지(Starvation)에서 세 번씩 키워서 센 count이다. size factor는 모두 1로 둔다. 예시 이름은 교재와 다르게 바꿔 썼고, 교재 문장을 인용한 곳도 이 이름에 맞췄다.

| 조건 | count | 평균 |
|---|---|---:|
| Ctrl | 100, 130, 90 | 106.667 |
| Starvation | 200, 250, 180 | 210 |

이 여섯 숫자가 노트를 따라가면서 이렇게 바뀐다.

| 단계 | 유전자 A | 노트 |
|---|---|---|
| dispersion: 보통의 NB 최대가능도 → Cox-Reid 보정 | 0.014786 → 0.025385 | [03](03_dispersion_estimation.md) |
| 최종 dispersion: 교재가 정해 둔 학습용 prior로 구한 MAP | 0.053147 | [03](03_dispersion_estimation.md) |
| log2 fold change와 표준오차 | 0.977280, 0.289057 | [04](04_glm_condition_batch.md) |
| Wald 통계량과 p 값 | 3.3809, 0.00072244 | [05](05_wald_vs_lrt.md) |

[08](08_one_gene_end_to_end.md)에서는 이 계산을 처음부터 끝까지 한 번에 다시 해 보고, 같은 count를 실제 `DESeq()`에 넣어 비교한다. 굶긴 세포에 glucose를 다시 넣은 세 번째 조건(Starvation+Glucose, 노트에서는 줄여서 Glucose)을 더한 예제도 나온다.

## 노트 목록

| 노트 | 궁금했던 것 | 교재 | 코드 |
|---|---|---|---|
| [00. DESeq2는 count 여섯 개로 무엇을 물을까?](00_overview.md) | DESeq2는 무엇을 비교하고, count에서 padj까지 어떤 단계를 거칠까 | 학습 안내, 1장, 부록 C·D | |
| [01. 같은 조건의 replicate인데 count가 왜 이렇게 다를까?](01_poisson_simulation.md) | Poisson 모형은 왜 biological replicate 사이의 퍼짐을 설명하지 못할까 | 2.1–2.2 | [노트북 (Python)](../01_poisson_simulation.ipynb) |
| [02. count는 얼마나 퍼지고, sequencing 깊이 차이는 어떻게 맞출까?](02_negative_binomial.md) | 음이항분포의 dispersion α와 size factor는 각각 무슨 일을 할까 | 2.3–2.4, 3장 | [노트북 (R)](../02_negative_binomial.ipynb) |
| [03. replicate 세 개로 dispersion α를 어떻게 정할까?](03_dispersion_estimation.md) | likelihood, Cox-Reid 보정, trend와 prior로 α를 정하는 과정 | 5–7장 | [R 스크립트](../03_dispersion_estimation.R) |
| [04. 조건과 batch를 넣은 모형에서 LFC와 SE는 어떻게 나올까?](04_glm_condition_batch.md) | design matrix, IRLS, 표준오차, 비교(contrast)의 공분산 | 4장, 8장 | [R 스크립트](../04_glm_condition_batch.R) |
| [05. 같은 fold change인데 왜 p 값이 다를까?](05_wald_vs_lrt.md) | Wald 검정과 LRT는 각각 무엇을 검정할까 | 9장, 10장, 12장 | [R 스크립트](../05_wald_vs_lrt.R) |
| [06. 유전자 수천 개를 한꺼번에 검정하면 p 값을 어떻게 읽어야 할까?](06_multiple_testing.md) | BH 보정, independent filtering, 여러 비교와 rescue 주장의 보정 | 11장, 14장 | [R 스크립트](../06_multiple_testing.R) |
| [07. 적은 count에서 나온 큰 fold change를 믿어도 될까?](07_lfc_shrinkage_and_qc.md) | LFC shrinkage, 결과표의 NA, QC 그림 | 13장 | |
| [08. count 여섯 개는 어떻게 p 값 하나가 될까?](08_one_gene_end_to_end.md) | 유전자 A를 끝까지 손으로 계산하고 실제 DESeq2와 비교하기 | 15장, 16장 | |

코드 열의 노트북과 스크립트는 노트의 주요 숫자만 따로 재현한다. 설명은 노트에 있다. 실행하면 결과를 노트에 적힌 숫자와 비교하고, 하나라도 다르면 오류를 내면서 멈춘다.

## 확인 방법

R 4.5.2와 DESeq2 1.50.2(Bioconductor 3.22)에서 돌렸고, 01의 시뮬레이션 일부는 Python(numpy, scipy)으로 했다. 노트에 있는 코드 블록은 전부 실제로 실행했고, 바로 아래 출력은 그때 나온 그대로다. 코드는 노트마다 위에서부터 한 세션에서 이어서 돌렸고, 세션을 나눈 곳은 그 노트에 적어 두었다. 그림 코드는 `notes/` 폴더에서 실행하면 `../figures/`에 PNG가 저장된다. 함수를 출력하면 찍히는 `<bytecode: 0x...>` 주소는 실행할 때마다 달라지니 신경 쓰지 않아도 된다.

나중에 저장소의 `environment.yml`로 새 환경을 만들어서 07을 뺀 노트들과 코드 파일을 다시 돌려 봤다. 1e-15 정도의 부동소수점 차이, 08이 저장하는 파일 크기, 일부 그림 파일의 렌더링 차이 말고는 결과가 같았다.

DESeq2가 어떻게 동작한다고 쓴 문장은 설치된 패키지 소스를 `deparse()`나 `args()`로 꺼내 보면서 확인했다. `fitDisp`, `fitBeta`처럼 C++로 된 부분은 설치본에는 컴파일된 파일만 있어서 1.50.2 소스 tarball에서 읽었다. apeglm, ashr, glmGamPoi, stageR, tximport는 설치하지 않았기 때문에, 이 패키지들에 대한 내용은 소스와 도움말만 보고 썼고 노트에 "미확인(실행)"이라고 표시해 두었다.

## 교재와 다르게 나온 것들

교재의 설명과 숫자는 대부분 실행 결과와 맞았다. 16장의 숫자도 전부 재현된다. 다르게 나온 것 중에서 기억해 둘 만한 것들을 모아 봤다.

- 16장의 Z, p, 신뢰구간을 마지막 자리까지 맞추려면 `optimize()`의 허용오차를 줄여서 MAP dispersion을 끝까지 수렴시켜야 함. 교재 코드를 그대로 쓰면 0.053146542에서 멈춤([08](08_one_gene_end_to_end.md)).
- 16장의 prior(중심 0.08, SD 0.7)는 교재가 설명하려고 정해 둔 값임. 같은 count를 실제 `DESeq()`에 넣으면 trend와 prior 폭을 다른 유전자들에서 추정하기 때문에, 최종 dispersion이 0.053이 아니라 배경 유전자에 따라 0.093–0.098로 나옴([08](08_one_gene_end_to_end.md)).
- dispersion outlier를 가르는 기준은 교재에 나오지 않음. 실제로는 log gene-wise 값이 log trend보다 `2 × sqrt(varLogDispEsts)` 넘게 클 때 outlier로 보고, 이 폭은 prior SD가 아니라 잔차를 robust하게 잰 산포임([03](03_dispersion_estimation.md)).
- 교재 6.5에는 조건 안에서 거의 흔들리지 않는 유전자가 나옴. 교재는 이 유전자의 α가 `~condition`에서 "작을 수 있다"고만 하는데, 실제로 돌려 보면 작아지는 정도가 아니라 DESeq2의 하한값 1e-8에 걸림([03](03_dispersion_estimation.md)).
- "초기 dispersion으로 평균을 먼저 추정한다"는 설명은 일부 design에서만 맞음. `~condition`처럼 그룹 평균만으로 풀리는 design에서는 α와 상관없는 그룹 평균을 씀([00](00_overview.md)).
- `lfcShrink(type = "normal")`이 돌려주는 `lfcSE`는 posterior SD가 아님. posterior SD라는 설명은 apeglm과 ashr에만 맞음([07](07_lfc_shrinkage_and_qc.md)).
- 부록 B 17의 "fitting 문제"는 두 경우로 나뉨. β가 수렴하지 않아도 p는 NA가 되지 않고, weights 때문에 계수를 추정할 수 없는 행만 p까지 NA가 됨([07](07_lfc_shrinkage_and_qc.md)).

노트마다 끝의 "더 깊이 보기"에 교재와 다르게 나온 부분을 빠짐없이 정리해 두었다.

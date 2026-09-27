# DESeq2 학습 노트

DESeq2가 유전자 하나의 count를 log2 fold change, 표준오차, p 값, 보정 p 값으로 바꾸는 과정을 따라간 노트다. 개인 학습 교재(「DESeq2 학습 교재 · COUNT에서 결론까지」)를 장 순서대로 따라가면서, 교재의 설명과 숫자를 설치된 DESeq2로 다시 계산해 확인했다. 교재 PDF는 저장소에 넣지 않았다.

평균, 분산, p 값은 알지만 likelihood나 GLM은 이름만 들어 본 사람을 기준으로 썼다.

## 읽는 법

처음이면 00부터 순서대로 읽는다. 각 노트는 질문 하나로 시작해서 구체적인 숫자로 보여 준 다음 수식을 붙인다. 소스 코드 발췌, 구현의 경계 사례, 교재와 대조한 메모는 노트 끝 "더 깊이 보기"의 접힌 블록에 있다. 처음 읽을 때는 건너뛰어도 흐름을 따라가는 데 지장이 없다.

모든 노트가 같은 예제를 쓴다. 교재 16장의 유전자 A이고, 세포를 보통 배지(Ctrl)와 glucose를 뺀 배지(Starvation)에서 세 번씩 키워 센 count다. size factor는 모두 1로 둔다. 조건 이름은 교재와 다르다. 교재 예시의 구조(대조군, 처리군, 처리 효과를 되돌리는 rescue 조건)는 그대로 두고 이름만 이 예시로 바꿨다. 교재 인용문과 검증 메모에 나오는 조건 이름도 같이 바꿨다.

| 조건 | count | 평균 |
|---|---|---:|
| Ctrl | 100, 130, 90 | 106.667 |
| Starvation | 200, 250, 180 | 210 |

이 여섯 숫자는 노트를 지나면서 다음 값이 된다.

| 단계 | 유전자 A의 값 | 노트 |
|---|---|---|
| dispersion: 일반 NB 최대가능도 → Cox-Reid 보정 | 0.014786 → 0.025385 | [03](03_dispersion_estimation.md) |
| 최종 dispersion: 교재가 정한 학습용 prior로 구한 MAP | 0.053147 | [03](03_dispersion_estimation.md) |
| log2 fold change와 표준오차 | 0.977280, 0.289057 | [04](04_glm_condition_batch.md) |
| Wald 통계량과 p 값 | 3.3809, 0.00072244 | [05](05_wald_vs_lrt.md) |

[08](08_one_gene_end_to_end.md)은 이 계산을 처음부터 끝까지 한 번에 다시 하고, 같은 count를 실제 `DESeq()`에 넣어 비교한다. 굶긴 세포에 glucose를 다시 넣은 세 번째 조건(Starvation+Glucose)을 더한 세 조건 예제도 나온다. 노트에서는 이 조건을 Glucose라고 줄여 부른다.

## 노트 목록

| 노트 | 답하는 질문 | 교재 | 실행 파일 |
|---|---|---|---|
| [00. DESeq2는 count 여섯 개로 무엇을 묻는가?](00_overview.md) | DESeq2는 무엇을 비교하고, count에서 padj까지 어떤 단계를 거치는가 | 학습 안내, 1장, 부록 C·D | |
| [01. 같은 조건의 replicate인데 count가 왜 이렇게 다를까?](01_poisson_simulation.md) | Poisson 모형은 biological replicate 사이의 퍼짐을 왜 설명하지 못하는가 | 2.1–2.2 | [노트북 (Python)](../01_poisson_simulation.ipynb) |
| [02. count는 얼마나 퍼지고, sequencing 깊이 차이는 어떻게 맞추나?](02_negative_binomial.md) | 음이항분포의 dispersion α와 size factor는 각각 무엇을 맡는가 | 2.3–2.4, 3장 | [노트북 (R)](../02_negative_binomial.ipynb) |
| [03. replicate 세 개로 dispersion α를 어떻게 정할까?](03_dispersion_estimation.md) | likelihood, Cox-Reid 보정, trend와 prior로 α를 정하는 과정 | 5–7장 | [R 스크립트](../03_dispersion_estimation.R) |
| [04. 조건과 batch를 넣은 모형에서 LFC와 SE는 어떻게 나오는가?](04_glm_condition_batch.md) | design matrix, IRLS, 표준오차, 비교(contrast)의 공분산 | 4장, 8장 | [R 스크립트](../04_glm_condition_batch.R) |
| [05. 같은 fold change인데 왜 p 값이 다를까?](05_wald_vs_lrt.md) | Wald 검정과 LRT는 각각 무엇을 검정하는가 | 9장, 10장, 12장 | [R 스크립트](../05_wald_vs_lrt.R) |
| [06. 유전자 수천 개를 한꺼번에 검정하면 p 값을 어떻게 읽어야 하나?](06_multiple_testing.md) | BH 보정, independent filtering, 여러 비교와 rescue 주장의 보정 | 11장, 14장 | [R 스크립트](../06_multiple_testing.R) |
| [07. 적은 count에서 나온 큰 fold change를 믿어도 될까?](07_lfc_shrinkage_and_qc.md) | LFC shrinkage, 결과표의 NA, QC 그림 | 13장 | |
| [08. count 여섯 개는 어떻게 p 값 하나가 되는가?](08_one_gene_end_to_end.md) | 유전자 A를 끝까지 손으로 계산하고 실제 DESeq2와 비교한다 | 15장, 16장 | |

실행 파일은 해당 노트 본문의 핵심 숫자만 따로 재현한다. 설명은 노트에 있다. 실행하면 결과를 노트에 적힌 숫자와 비교하고, 다르면 오류를 내며 멈춘다.

## 어떻게 확인했는가

R 4.5.2와 DESeq2 1.50.2(Bioconductor 3.22)로 확인했다. 01의 시뮬레이션 일부는 Python(numpy, scipy)으로 했다. 노트의 코드 블록은 모두 실제로 실행했고, 바로 아래의 출력은 그 실행 결과를 그대로 옮긴 것이다. 코드는 노트마다 위에서 아래 순서로 한 세션에서 실행했고, 세션을 나눈 곳은 그 노트에 적었다. 그림 코드는 `notes/` 폴더에서 실행하면 `../figures/`에 PNG를 쓴다. 함수를 출력할 때 찍히는 `<bytecode: 0x...>` 주소는 실행할 때마다 달라진다.

나중에 저장소의 `environment.yml`로 새 환경을 만들어 07을 뺀 노트들과 실행 파일을 다시 돌렸다. 1e-15 수준의 부동소수점 차이, 08이 저장하는 파일 크기, 일부 그림 파일의 렌더링 차이만 있었다.

DESeq2의 동작에 관한 문장은 설치된 패키지의 소스를 `deparse()`나 `args()`로 꺼내 대조했다. 설치본에는 컴파일된 C++만 있어서 `fitDisp`, `fitBeta`의 본문은 1.50.2 소스 tarball에서 읽었다. apeglm, ashr, glmGamPoi, stageR, tximport는 설치하지 않았다. 이 패키지들에 관한 내용은 소스와 도움말로만 확인했고, 해당 메모에 "미확인(실행)"으로 적었다.

## 교재와 다르게 확인된 점

교재는 R 없이 작성되어 스크립트를 실행해 보지 못했다고 스스로 밝힌다(p.3). 대부분의 설명과 숫자는 실행 결과와 맞았다. 16장의 숫자도 모두 재현된다. 다르게 나온 것 가운데 중요한 것은 다음과 같다.

- 16장의 Z, p, 신뢰구간을 마지막 자리까지 맞추려면 `optimize()`의 허용오차를 줄여 MAP dispersion을 끝까지 수렴시켜야 한다. 교재 코드를 그대로 쓰면 0.053146542에서 멈춘다([08](08_one_gene_end_to_end.md)).
- 16장의 prior(중심 0.08, SD 0.7)는 교재가 학습용으로 정한 값이다. 같은 count를 실제 `DESeq()`에 넣으면 trend와 prior 폭을 다른 유전자에서 추정하므로 최종 dispersion은 0.053이 아니라 배경 유전자에 따라 0.093–0.098이 된다([08](08_one_gene_end_to_end.md)).
- dispersion outlier를 가르는 기준은 교재에 없다. 실제 규칙은 log gene-wise 값이 log trend보다 `2 × sqrt(varLogDispEsts)` 넘게 클 때이고, 이 폭은 prior SD가 아니라 잔차의 robust한 산포로 잰다([03](03_dispersion_estimation.md)).
- 교재 6.5는 유전자 A의 α가 `~condition`에서 "작을 수 있다"고 적는다. 실제로는 작아지는 정도가 아니라 DESeq2의 하한값 1e-8에 걸린다([03](03_dispersion_estimation.md)).
- "초기 dispersion으로 평균을 먼저 추정한다"(학습 안내)는 일부 design에만 맞다. `~condition`처럼 그룹 평균만으로 풀리는 design에서는 α와 무관한 그룹 평균을 쓴다([00](00_overview.md)).
- `lfcShrink(type = "normal")`이 돌려주는 `lfcSE`는 posterior SD가 아니다. posterior SD라는 설명은 apeglm과 ashr에만 맞다([07](07_lfc_shrinkage_and_qc.md)).
- 부록 B 17의 "fitting 문제"는 두 경우로 갈린다. β가 수렴하지 않아도 p는 NA가 되지 않는다. weights 때문에 계수를 추정할 수 없는 행만 p까지 NA가 된다([07](07_lfc_shrinkage_and_qc.md)).

전체 목록은 각 노트 끝 "더 깊이 보기"의 "교재와 다른 점 (검증 메모)" 블록에 있다.

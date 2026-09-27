# 03 · How is the dispersion alpha estimated from three replicates?
#
# Reproduces the numbers in notes/03_dispersion_estimation.md (Korean), which
# explains each step. Run from the repository root:
#
#     Rscript 03_dispersion_estimation.R
#
# The stopifnot() checks compare against the values printed in the note
# (R 4.5.2, DESeq2 1.50.2).

suppressPackageStartupMessages(suppressWarnings(library(DESeq2)))
section <- function(title) cat("\n==", title, "==\n")
near <- function(x, expected, digits) all(abs(unname(x) - expected) <= 0.5 * 10^-digits + 1e-9)  # x rounds to expected

# Gene A (study text, ch. 16): three control and three glucose-starved replicates, size factors all 1.
K  <- c(100, 130, 90, 200, 250, 180)
cd <- data.frame(condition = factor(rep(c("Ctrl", "Starvation"), each = 3), levels = c("Ctrl", "Starvation")))
maximize <- function(f) exp(optimize(f, log(c(1e-8, 10)), maximum = TRUE, tol = 1e-10)$maximum)  # over log(alpha)

section("1. Moment estimate: pooling the two groups inflates alpha")
moment <- function(k) (var(k) - mean(k)) / mean(k)^2
cat(sprintf("all six counts %.6f | Ctrl %.5f | Starvation %.5f\n", moment(K), moment(K[1:3]), moment(K[4:6])))
stopifnot(near(moment(K), 0.149119, 6), near(moment(K[1:3]), 0.02871, 5), near(moment(K[4:6]), 0.02472, 5))

section("2. Maximum likelihood with the group means held fixed")
mu    <- rep(c(mean(K[1:3]), mean(K[4:6])), each = 3)
ll_nb <- function(a) sum(dnbinom(K, mu = mu, size = 1 / a, log = TRUE))
cand  <- c(0.001, 0.005, 0.015, 0.05, 0.2)
print(setNames(round(sapply(cand, ll_nb), 3), cand))
a_mle <- maximize(function(th) ll_nb(exp(th)))
cat(sprintf("alpha, NB maximum likelihood: %.6f\n", a_mle))
stopifnot(near(a_mle, 0.014786, 6))

section("3. Cox-Reid adjustment for estimating the means from the same counts")
# l_CR(alpha) = l(alpha) - 1/2 log det(X' W X),  W = diag(mu / (1 + alpha mu))
X     <- model.matrix(~ condition, cd)
cr    <- function(a) -0.5 * log(det(crossprod(X, X * (mu / (1 + a * mu)))))
ll_cr <- function(a) ll_nb(a) + cr(a)
a_cr  <- maximize(function(th) ll_cr(exp(th)))
cat(sprintf("alpha, Cox-Reid adjusted: %.6f (%.2f times the MLE)\n", a_cr, a_cr / a_mle))
stopifnot(near(a_cr, 0.025385, 6))

# DESeq2's gene-wise estimate is this adjusted maximum; useCR = FALSE gives back the plain MLE.
dds1 <- DESeqDataSetFromMatrix(matrix(as.integer(K), nrow = 1), cd, ~ condition)
sizeFactors(dds1) <- rep(1, 6)
gw_cr <- mcols(estimateDispersionsGeneEst(dds1, quiet = TRUE))$dispGeneEst
gw_ml <- mcols(estimateDispersionsGeneEst(dds1, useCR = FALSE, quiet = TRUE))$dispGeneEst
cat(sprintf("DESeq2 dispGeneEst %.8f | with useCR = FALSE %.8f\n", gw_cr, gw_ml))
stopifnot(near(gw_cr, 0.02538756, 8), near(gw_ml, 0.01478521, 8))

# The design decides what spread is left over (the two genes of section 6.5 in the study text).
cts <- rbind(A6.5 = c(100, 105, 95, 200, 205, 195), B6.5 = c(50, 140, 80, 120, 300, 170))
storage.mode(cts) <- "integer"
gene_wise <- function(design) {
  d <- DESeqDataSetFromMatrix(cts, cd, design); sizeFactors(d) <- rep(1, 6)
  mcols(estimateDispersionsGeneEst(d, quiet = TRUE))$dispGeneEst
}
by_design <- cbind(with_condition = gene_wise(~ condition), intercept_only = gene_wise(~ 1))
rownames(by_design) <- rownames(cts)
print(signif(by_design, 4))
stopifnot(by_design["A6.5", "with_condition"] == 1e-8,           # hits DESeq2's lower bound minDisp
          near(by_design[, "intercept_only"], c(0.1323, 0.3407), 4))

section("4. Three replicates give noisy gene-wise estimates")
set.seed(2)
sim <- matrix(rnbinom(2000 * 6, mu = 100, size = 1 / 0.1), ncol = 6); storage.mode(sim) <- "integer"
ds  <- DESeqDataSetFromMatrix(sim, data.frame(condition = factor(rep(c("a", "b"), each = 3))), ~ condition)
sizeFactors(ds) <- rep(1, 6)
gw  <- mcols(estimateDispersionsGeneEst(ds, quiet = TRUE))$dispGeneEst
q   <- quantile(gw, c(0.05, 0.5, 0.95))
cat(sprintf("true alpha 0.1 for all 2000 genes; gene-wise 5%% / 50%% / 95%%: %s; share at the lower bound %.3f\n",
            paste(signif(q, 3), collapse = " / "), mean(gw <= 1e-8)))
stopifnot(near(unname(q), c(0.0107, 0.0825, 0.245), c(4, 4, 3)))

section("5. Trend and prior from all genes (empirical Bayes)")
set.seed(1)
dds <- makeExampleDESeqDataSet(n = 2000, m = 6, betaSD = 1)      # ~ condition, 3 vs 3
dds <- estimateDispersions(estimateSizeFactors(dds), quiet = TRUE)
fn  <- dispersionFunction(dds)
trend_coef <- attr(fn, "coefficients")                           # alpha_tr(q) = asymptDisp + extraPois / q
vld <- attr(fn, "varLogDispEsts"); prior_var <- attr(fn, "dispPriorVar")
cat(sprintf("trend: %.5f + %.5f / mean\n", trend_coef[["asymptDisp"]], trend_coef[["extraPois"]]))
cat(sprintf("spread of log gene-wise estimates around the trend %.3f, expected from estimation noise %.3f,",
            vld, trigamma((6 - 2) / 2)), "prior variance", prior_var, "(the lower bound)\n")
stopifnot(near(unname(trend_coef), c(0.09389, 6.19823), 5), near(vld, 0.679, 3), prior_var == 0.25)

section("6. MAP for gene A with the study text's teaching prior")
# The prior log(alpha) ~ N(log 0.08, 0.7^2) is fixed by the study text, not estimated from data.
post  <- function(th) ll_cr(exp(th)) - (th - log(0.08))^2 / (2 * 0.7^2)
a_map <- maximize(post)
cat(sprintf("MLE %.6f -> Cox-Reid %.6f -> MAP %.6f\n", a_mle, a_cr, a_map))
stopifnot(near(a_map, 0.053147, 6), a_cr < a_map, a_map < 0.08)   # pulled toward the trend, not toward 0

section("7. The dispersion columns DESeq2 stores")
df <- as.data.frame(mcols(dds)); ok <- !df$allZero
print(round(df[c("gene733", "gene798", "gene12", "gene15", "gene1", "gene2"),
               c("baseMean", "dispGeneEst", "dispFit", "dispMAP", "dispOutlier", "dispersion")], 4))
# Outlier rule: log gene-wise > log trend + 2 * sqrt(varLogDispEsts), not the prior SD
outlier <- with(df[ok, ], log(dispGeneEst) > log(dispFit) + 2 * sqrt(vld))
kept <- df[ok & !df$dispOutlier, ]
cat(sprintf("outliers %d (with the prior SD instead: %d) | non-outliers %d, MAP above gene-wise %d, below trend %d, moved down %d\n",
            sum(outlier), sum(with(df[ok, ], log(dispGeneEst) > log(dispFit) + 2 * sqrt(prior_var))),
            nrow(kept), sum(kept$dispMAP > kept$dispGeneEst), sum(kept$dispGeneEst < kept$dispFit),
            sum(kept$dispGeneEst < kept$dispFit & kept$dispMAP < kept$dispGeneEst)))
stopifnot(identical(outlier, df$dispOutlier[ok]), sum(outlier) == 8,
          sum(with(df[ok, ], log(dispGeneEst) > log(dispFit) + 2 * sqrt(prior_var))) == 79,
          nrow(kept) == 1969, sum(kept$dispMAP > kept$dispGeneEst) == 1326,
          identical(dispersions(dds), df$dispersion),
          all(with(df[ok, ], ifelse(dispOutlier, dispersion == dispGeneEst, dispersion == dispMAP))))

cat("\nAll checks passed.\n")

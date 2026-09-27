# 04 - The negative binomial GLM: design matrix, IRLS and the standard error
#
# The numbers from notes/04_glm_condition_batch.md, without the explanations
# (those are in the note, in Korean). Run it from the repo root:
#
#     Rscript 04_glm_condition_batch.R
#
# Each stopifnot() checks a result against the value printed in the note.
# I ran it with R 4.5.2 and DESeq2 1.50.2.

suppressPackageStartupMessages(suppressWarnings(library(DESeq2)))
section <- function(title) cat("\n==", title, "==\n")
near <- function(x, expected, digits) all(abs(unname(x) - expected) <= 0.5 * 10^-digits + 1e-9)  # x rounds to expected

k  <- c(100, 130, 90, 200, 250, 180)                  # gene A, size factors all 1
cd <- data.frame(condition = factor(rep(c("Ctrl", "Starvation"), each = 3), levels = c("Ctrl", "Starvation")))
alpha <- 0.053147                                     # gene A's final dispersion (note 03)

section("1. The design matrix")
X <- model.matrix(~ condition, cd)
print(t(X[, ]))                                       # one column per sample
# Donors A and B each give Ctrl, Starvation and Starvation+Glucose ("Glucose"): design ~ pair + condition
mdp <- data.frame(pair = factor(rep(c("A", "B"), each = 3)),
                  condition = factor(rep(c("Ctrl", "Starvation", "Glucose"), 2), levels = c("Ctrl", "Starvation", "Glucose")))
Xp <- model.matrix(~ pair + condition, mdp)
cat("~ pair + condition: rank", qr(Xp)$rank, "of", ncol(Xp), "columns, residual df", nrow(Xp) - qr(Xp)$rank, "\n")
stopifnot(qr(Xp)$rank == 4, nrow(Xp) - qr(Xp)$rank == 2)
# Batch fully confounded with condition: no option can separate the two
mdc <- data.frame(batch = factor(c(1, 1, 1, 2, 2, 2)), condition = factor(rep(c("Ctrl", "Starvation"), each = 3)))
set.seed(1)
err <- tryCatch(DESeqDataSetFromMatrix(matrix(rpois(20 * 6, 50), nrow = 20), mdc, ~ batch + condition),
                error = conditionMessage)
cat("~ batch + condition with batch = condition:", sub("\n.*", "", err), "\n")
stopifnot(grepl("not full rank", err))

section("2. IRLS by hand for gene A, dispersion fixed")
# eta = offset + X b, mu = exp(eta), z = eta + (K - mu)/mu, w = mu / (1 + alpha mu), b <- (X'WX)^-1 X'W (z - offset)
o <- log(rep(1, 6))                                   # offset = log(size factor)
b <- qr.solve(X, log(k + 0.1))                        # starting value, as in DESeq2
for (t in 0:3) {
  eta <- as.numeric(o + X %*% b); mu <- exp(eta)
  score <- crossprod(X, (k - mu) / (1 + alpha * mu))
  cat(sprintf("iter %d  beta_S = %.9f  score_S = %9.2e\n", t, b[2] / log(2), score[2]))
  z <- eta + (k - mu) / mu; w <- mu / (1 + alpha * mu)
  b <- as.numeric(solve(crossprod(X, X * w), crossprod(X, w * (z - o))))
}
cat("exp(b):", exp(b), "(Ctrl mean, Starvation/Ctrl ratio)\n")
stopifnot(near(b[2] / log(2), 0.977279923, 9), all.equal(exp(b), c(320 / 3, 1.96875)))

section("3. The same fit in DESeq2")
dds <- DESeqDataSetFromMatrix(matrix(as.integer(k), nrow = 1, dimnames = list("geneA", paste0("s", 1:6))), cd, ~ condition)
sizeFactors(dds) <- rep(1, 6)
dispersions(dds) <- alpha                             # skip dispersion estimation, use alpha
dds <- nbinomWaldTest(dds, quiet = TRUE)
res <- results(dds, independentFiltering = FALSE, cooksCutoff = FALSE)
print(as.data.frame(res)[, c("log2FoldChange", "lfcSE")], digits = 9)
stopifnot(near(res$log2FoldChange, 0.97728, 6), near(res$lfcSE, 0.289057, 6))

section("4. The standard error comes from X'WX")
mu <- exp(as.numeric(o + X %*% b)); w <- mu / (1 + alpha * mu)
V  <- solve(crossprod(X, X * w))                      # approximate covariance of b (natural log)
cat(sprintf("weights Ctrl %.4f, Starvation %.4f | Var(b_S) %.8f = 1/%.2f + 1/%.2f | SE(beta_S) %.7f\n",
            w[1], w[4], V[2, 2], sum(w[1:3]), sum(w[4:6]), sqrt(V[2, 2]) / log(2)))
stopifnot(near(V[2, 2], 0.04014363, 8), all.equal(V[2, 2], 1 / sum(w[1:3]) + 1 / sum(w[4:6])),
          near(sqrt(V[2, 2]) / log(2), 0.2890566, 7))
# Same means, different alpha: the estimate stays, the standard error grows
mu_fix <- rep(c(320 / 3, 210), each = 3); bS <- log(210 / (320 / 3))
sweep_alpha <- t(sapply(c(0, 0.014786, 0.025385, 0.053147, 0.1, 0.5), function(a) {
  w <- mu_fix / (1 + a * mu_fix); V <- solve(crossprod(X, X * w))
  c(alpha = a, LFC = bS / log(2), SE_log2 = sqrt(V[2, 2]) / log(2), Z = bS / sqrt(V[2, 2]))
}))
print(round(sweep_alpha, 6))
stopifnot(near(sweep_alpha[, "SE_log2"], c(0.099036, 0.17414, 0.212207, 0.289057, 0.385443, 0.838807), 6))

section("5. Deeper sequencing cannot push a sample's weight past 1/alpha")
mu_grid <- c(1, 10, 100, 1000, 5000)
weights <- outer(mu_grid, c(0.01, 0.1, 0.5), function(m, a) m / (1 + a * m))
dimnames(weights) <- list(paste0("mu=", mu_grid), paste0("alpha=", c(0.01, 0.1, 0.5)))
print(round(weights, 3))
a <- 0.1                                              # three samples at expected count 100
depth <- c(now = 3 * 100 / (1 + a * 100), depth_x2 = 3 * 200 / (1 + a * 200), replicates_x2 = 6 * 100 / (1 + a * 100))
print(round(depth, 5))
stopifnot(near(depth, c(27.27273, 28.57143, 54.54545), 5))

section("6. The rescue contrast Glucose - Starvation needs the covariance")
set.seed(1)
dds3 <- makeExampleDESeqDataSet(n = 2000, m = 9, betaSD = 1)
dds3$condition <- factor(rep(c("Ctrl", "Starvation", "Glucose"), each = 3), levels = c("Ctrl", "Starvation", "Glucose"))
design(dds3) <- ~ condition
dds3 <- DESeq(dds3, quiet = TRUE)
g <- "gene4"
se_S  <- results(dds3, name = "condition_Starvation_vs_Ctrl")[g, "lfcSE"]
se_G  <- results(dds3, name = "condition_Glucose_vs_Ctrl")[g, "lfcSE"]
se_GS <- results(dds3, contrast = c("condition", "Glucose", "Starvation"))[g, "lfcSE"]
cat(sprintf("%s: SE_S %.4f, SE_G %.4f | results(contrast) %.4f | naive sum %.4f, naive sqrt(sum of squares) %.4f\n",
            g, se_S, se_G, se_GS, se_S + se_G, sqrt(se_S^2 + se_G^2)))
X3 <- attr(dds3, "modelMatrix"); mu3 <- assays(dds3)[["mu"]][g, ]
w3 <- mu3 / (1 + dispersions(dds3)[rownames(dds3) == g] * mu3)
V3 <- solve(crossprod(X3, X3 * w3)) / log(2)^2       # log2-scale covariance
cc <- c(0, -1, 1)
cat(sprintf("Cov(b_G, b_S) %.7f = Var(b_0) %.7f | sqrt(c' V c) %.7f\n", V3[2, 3], V3[1, 1], sqrt(drop(t(cc) %*% V3 %*% cc))))
stopifnot(near(c(se_S, se_G, se_GS), c(0.5419, 0.5419, 0.5415), 4),
          all.equal(V3[2, 3], V3[1, 1]), near(sqrt(drop(t(cc) %*% V3 %*% cc)), 0.5415, 4))

cat("\nAll checks passed.\n")

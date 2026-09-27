# 06 · Which genes is padj adjusted over, and how do you test "down under starvation, back up once glucose returns"?
#
# Reproduces the numbers in notes/06_multiple_testing.md (Korean), which explains
# each step. Run from the repository root:
#
#     Rscript 06_multiple_testing.R
#
# The stopifnot() checks compare against the values printed in the note
# (R 4.5.2, DESeq2 1.50.2).

suppressPackageStartupMessages(suppressWarnings(library(DESeq2)))
options(width = 110)
section <- function(title) cat("\n==", title, "==\n")
near <- function(x, expected, digits) all(abs(unname(x) - expected) <= 0.5 * 10^-digits + 1e-9)  # x rounds to expected

section("1. Picking genes by p < 0.05: how many are false?")
# 1000 null genes (Z ~ N(0, 1)) and 100 changed genes (Z ~ N(3, 1)), repeated 200 times
set.seed(20260926); B <- 200; m0 <- 1000; m1 <- 100; alpha <- 0.05
fdp_bh <- fdp_bonf <- fdp_raw <- power_bh <- fwer_holm <- numeric(B)
for (b in seq_len(B)) {
  z <- c(rnorm(m0, 0), rnorm(m1, 3)); truth <- rep(c(FALSE, TRUE), c(m0, m1))
  p <- 2 * pnorm(-abs(z))
  R  <- p.adjust(p, "BH") < alpha;         fdp_bh[b] <- sum(R & !truth) / max(sum(R), 1); power_bh[b] <- mean(R[truth])
  Rb <- p.adjust(p, "bonferroni") < alpha; fdp_bonf[b] <- sum(Rb & !truth) / max(sum(Rb), 1)
  Rh <- p.adjust(p, "holm") < alpha;       fwer_holm[b] <- any(Rh & !truth)
  Rr <- p < alpha;                         fdp_raw[b] <- sum(Rr & !truth) / max(sum(Rr), 1)
}
cat(sprintf("raw p < 0.05: FDR %.3f | BH: FDR %.4f, power %.3f, P(any false) %.2f | Bonferroni FDR %.4f | Holm P(any false) %.3f\n",
            mean(fdp_raw), mean(fdp_bh), mean(power_bh), mean(fdp_bh > 0), mean(fdp_bonf), mean(fwer_holm)))
cat(sprintf("BH false discovery proportion in single data sets ranges from %.3f to %.3f\n", min(fdp_bh), max(fdp_bh)))
stopifnot(near(mean(fdp_raw), 0.372, 3), near(mean(fdp_bh), 0.0459, 4), near(mean(fwer_holm), 0.045, 3))

section("2. Benjamini-Hochberg by hand (study text 11.2)")
p <- c(0.001, 0.010, 0.030, 0.040, 0.200)
m <- length(p); ps <- sort(p)
bh_hand <- pmin(1, rev(cummin(rev(m * ps / seq_len(m)))))   # m p_(i) / i, then the running minimum from the end
print(data.frame(rank = 1:m, p = ps, m_p_over_i = m * ps / seq_len(m), BH_hand = bh_hand, p.adjust = p.adjust(ps, "BH")))
stopifnot(all.equal(bh_hand, p.adjust(ps, "BH")), all.equal(bh_hand, c(0.005, 0.025, 0.05, 0.05, 0.2)))
# The same p-value gets a different padj in a different family
set.seed(1)
fam <- c(A = p.adjust(p, "BH")[3], B = p.adjust(c(0.03, 0.5, 0.6, 0.7, 0.9), "BH")[1], C = p.adjust(c(0.03, runif(1000)), "BH")[1])
cat("padj of p = 0.03 in three families:", round(fam, 4), "\n")
stopifnot(near(fam, c(0.05, 0.15, 0.9418), 4))

section("3. DESeq2's padj family: one results() call, minus filtered genes")
# Ctrl / Starvation / Glucose (starved, then glucose added back), four samples each, 2000 genes with known effects (log2 scale)
set.seed(2026); n <- 2000; m <- 12
base <- makeExampleDESeqDataSet(n = n, m = m)
cond <- factor(rep(c("Ctrl", "Starvation", "Glucose"), each = 4), levels = c("Ctrl", "Starvation", "Glucose"))
X    <- model.matrix(~ cond)
cls  <- rep(c("null", "starvation_only", "rescued", "partial", "glucose_only"), c(1500, 150, 150, 100, 100))
eff  <- sample(c(-1, 1), n, TRUE) * runif(n, 1, 2.5)
bS   <- ifelse(cls %in% c("starvation_only", "rescued", "partial"), eff, 0)
bG   <- ifelse(cls == "starvation_only", bS, ifelse(cls == "partial", 0.5 * bS, ifelse(cls == "glucose_only", eff, 0)))
mu   <- t(2^(X %*% t(cbind(mcols(base)$trueIntercept, bS, bG))))
cnt  <- matrix(rnbinom(m * n, mu = mu, size = 1 / mcols(base)$trueDisp), ncol = m,
               dimnames = list(paste0("gene", 1:n), paste0("sample", 1:m)))
cnt[5, 1] <- 50000L; cnt[7, 6] <- 80000L                     # two planted count outliers
dds <- suppressMessages(DESeqDataSetFromMatrix(cnt, DataFrame(condition = cond), ~ condition))
dds <- DESeq(dds, quiet = TRUE)

r_if   <- results(dds, contrast = c("condition", "Starvation", "Ctrl"))                 # default: filtering on, alpha = 0.1
r_if05 <- results(dds, contrast = c("condition", "Starvation", "Ctrl"), alpha = 0.05)
r_no   <- results(dds, contrast = c("condition", "Starvation", "Ctrl"), independentFiltering = FALSE)
summarize <- function(r) c(NA_pvalue = sum(is.na(r$pvalue)), NA_padj = sum(is.na(r$padj)),
                           padj_0.05 = sum(r$padj < 0.05, na.rm = TRUE), cutoff = round(c(metadata(r)$filterThreshold, NA)[[1]], 2))
tab <- rbind(filter_alpha_0.1 = summarize(r_if), filter_alpha_0.05 = summarize(r_if05), no_filter = summarize(r_no))
print(tab)
filtered <- is.na(r_if$padj) & !is.na(r_if$pvalue)
cat("stat and pvalue identical across the three calls:", identical(r_if$pvalue, r_no$pvalue), identical(r_if05$stat, r_no$stat),
    "| padj NA only:", sum(filtered), "genes, all below the baseMean cutoff:", all(r_if$baseMean[filtered] < metadata(r_if)$filterThreshold), "\n")
stopifnot(tab[, "NA_pvalue"] == 9, tab[, "NA_padj"] == c(125, 472, 9), tab[, "padj_0.05"] == c(161, 163, 158),
          identical(r_if$pvalue, r_no$pvalue), sum(filtered) == 116)
# Filtering on an independent statistic keeps the FDR; pre-selecting by p does not (all genes null)
set.seed(5); B <- 1000; m <- 2000
any_disc <- matrix(NA, B, 3, dimnames = list(NULL, c("no_filter", "filter_by_baseMean", "filter_by_p")))
for (b in seq_len(B)) {
  mu_b <- rlnorm(m, 3, 1.5); p <- runif(m)
  any_disc[b, 1] <- any(p.adjust(p, "BH") < 0.05)
  keep <- mu_b > quantile(mu_b, 0.5); any_disc[b, 2] <- any(p.adjust(p[keep], "BH") < 0.05)
  keep <- p < 0.5;                     any_disc[b, 3] <- any(p.adjust(p[keep], "BH") < 0.05)
}
cat("all-null FDR at BH 0.05:", round(colMeans(any_disc), 3), "\n")
stopifnot(near(colMeans(any_disc), c(0.041, 0.046, 0.091), 3))

section("4. Three comparisons: which family?")
# Global null: correcting each comparison separately does not control the FDR of the combined list
set.seed(16); B <- 2000; mg <- 1000; se_g <- 0.3; se_d <- sqrt(2) * se_g
any_sep <- any_glob <- one <- numeric(B)
for (b in seq_len(B)) {
  Cm <- rnorm(mg, 0, se_g); Sm <- rnorm(mg, 0, se_g); Gm <- rnorm(mg, 0, se_g)
  Pm <- 2 * pnorm(-abs(cbind(Sm - Cm, Gm - Sm, Gm - Cm) / se_d))
  sep <- apply(Pm, 2, p.adjust, method = "BH") < 0.05
  any_sep[b] <- any(sep); one[b] <- any(sep[, 1]); any_glob[b] <- any(p.adjust(as.vector(Pm), "BH") < 0.05)
}
cat(sprintf("FDR: one comparison %.3f | three separately corrected lists combined %.3f | global BH %.3f\n",
            mean(one), mean(any_sep), mean(any_glob)))
stopifnot(near(c(mean(one), mean(any_sep), mean(any_glob)), c(0.045, 0.128, 0.046), 3))

# Global BH on the DESeq2 results: pool the raw p-values of the three contrasts (study text 11.4)
res_list <- list(
  S_C = results(dds, contrast = c("condition", "Starvation", "Ctrl"), independentFiltering = FALSE),
  G_S = results(dds, contrast = c("condition", "Glucose", "Starvation"), independentFiltering = FALSE),
  G_C = results(dds, contrast = c("condition", "Glucose", "Ctrl"),   independentFiltering = FALSE))
P <- do.call(cbind, lapply(res_list, function(x) x$pvalue)); P[is.na(P)] <- 1
Q_sep    <- sapply(res_list, function(x) x$padj)
Q_global <- matrix(p.adjust(as.vector(P), "BH"), nrow = nrow(P))
Q_BY     <- matrix(p.adjust(as.vector(P), "BY"), nrow = nrow(P))
Q_twice  <- p.adjust(as.vector(Q_sep), "BH")                  # wrong: BH applied to padj again
cells <- c(per_contrast = sum(Q_sep < 0.05, na.rm = TRUE), global_BH = sum(Q_global < 0.05),
           global_BY = sum(Q_BY < 0.05), BH_of_padj = sum(Q_twice < 0.05, na.rm = TRUE))
print(cells)
stopifnot(cells == c(360, 367, 211, 162))

section("5. 'Starvation lowers it and Glucose brings it back': test both parts together")
# Two one-sided p-values; the joint (intersection-union) p-value is their maximum.
set.seed(4); N <- 1e6; se_g <- 0.3; se_d <- sqrt(2) * se_g
Cm <- rnorm(N, 0, se_g); Sm <- rnorm(N, -1.2, se_g); Gm <- rnorm(N, -1.2, se_g)   # starvation effect real, no rescue
pS <- pnorm((Sm - Cm) / se_d); pR <- pnorm((Gm - Sm) / se_d, lower.tail = FALSE)
cat(sprintf("only one part true: P(max(pS, pR) <= 0.05) = %.4f, P(min(pS, pR) <= 0.05) = %.4f\n",
            mean(pmax(pS, pR) <= 0.05), mean(pmin(pS, pR) <= 0.05)))
stopifnot(near(mean(pmax(pS, pR) <= 0.05), 0.0498, 4), near(mean(pmin(pS, pR) <= 0.05), 0.8817, 4))

# Intersecting two significant lists vs the joint p-value, 200 simulated experiments
set.seed(1); B <- 200; se_g <- 0.30
n_cls <- c(A = 100, B = 400, C = 400, D = 1100); cls2 <- rep(names(n_cls), n_cls)   # A: both parts true
dS_true <- c(A = -1.2, B = -1.2, C = 0, D = 0)[cls2]; dR_true <- c(A = 1.2, B = 0, C = 1.2, D = 0)[cls2]
joint_true <- cls2 == "A"; fdp <- function(R) sum(R & !joint_true) / max(sum(R), 1)
out <- matrix(NA, B, 4, dimnames = list(NULL, c("intersect_one_sided", "intersect_two_sided_signs", "joint_max_p",
                                                "intersect_independent_Starvation")))
for (b in seq_len(B)) {
  Cm <- rnorm(length(cls2), 0, se_g); Sm <- rnorm(length(cls2), dS_true, se_g); Gm <- rnorm(length(cls2), dS_true + dR_true, se_g)
  dS <- Sm - Cm; dR <- Gm - Sm; se_d <- sqrt(2) * se_g
  pS <- pnorm(dS / se_d); pR <- pnorm(dR / se_d, lower.tail = FALSE)
  p2S <- 2 * pnorm(-abs(dS / se_d)); p2R <- 2 * pnorm(-abs(dR / se_d))
  Sm2 <- rnorm(length(cls2), dS_true, se_g)                   # a second, independent Starvation group
  pR_ind <- pnorm((Gm - Sm2) / se_d, lower.tail = FALSE)
  out[b, ] <- c(fdp(p.adjust(pS, "BH") < 0.05 & p.adjust(pR, "BH") < 0.05),
                fdp(p.adjust(p2S, "BH") < 0.05 & p.adjust(p2R, "BH") < 0.05 & dS < 0 & dR > 0),
                fdp(p.adjust(pmax(pS, pR), "BH") < 0.05),
                fdp(p.adjust(pS, "BH") < 0.05 & p.adjust(pR_ind, "BH") < 0.05))
}
print(round(colMeans(out), 3))   # sharing one Starvation group makes the intersection worse (0.132 vs 0.096)
stopifnot(near(colMeans(out), c(0.132, 0.077, 0.021, 0.096), 3))

# The two comparisons share the Starvation mean, so they are negatively correlated even for null genes
lfc <- sapply(res_list, function(x) x$log2FoldChange)
null_genes <- cls == "null" & !is.na(lfc[, 1])
r_SR <- cor(lfc[null_genes, "S_C"], lfc[null_genes, "G_S"])
cat(sprintf("null genes: cor(Starvation - Ctrl, Glucose - Starvation) = %.3f (theory -0.5)\n", r_SR))
stopifnot(near(r_SR, -0.565, 3))

section("6. 'Not significant' is not 'the same': equivalence test with lessAbs")
Tt <- log2(1.2)                                               # tolerance: within 1.2-fold of Ctrl
r_eq <- results(dds, contrast = c("condition", "Glucose", "Ctrl"), lfcThreshold = Tt, altHypothesis = "lessAbs", alpha = 0.05)
gc   <- res_list$G_C
tost <- function(b, s) max(pnorm((b - Tt) / s), pnorm((b + Tt) / s, lower.tail = FALSE))
cat(sprintf("Glucose vs Ctrl: two-sided p > 0.3 for %d genes, equivalent within 1.2-fold for %d genes\n",
            sum(gc$pvalue > 0.3, na.rm = TRUE), sum(r_eq$padj < 0.05, na.rm = TRUE)))
cat(sprintf("gene A's SE (0.289057) with an estimate of exactly 0: equivalence p = %.3f\n", tost(0, 0.289057)))
stopifnot(sum(gc$pvalue > 0.3, na.rm = TRUE) == 1292, sum(r_eq$padj < 0.05, na.rm = TRUE) == 0,
          near(tost(0, 0.289057), 0.181, 3))

cat("\nAll checks passed.\n")

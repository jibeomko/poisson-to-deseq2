# 05 · Why can the same fold change get different p-values, and what does the LRT test instead?
#
# Reproduces the numbers in notes/05_wald_vs_lrt.md (Korean), which explains
# each step. Run from the repository root:
#
#     Rscript 05_wald_vs_lrt.R
#
# The stopifnot() checks compare against the values printed in the note
# (R 4.5.2, DESeq2 1.50.2).

suppressPackageStartupMessages(suppressWarnings(library(DESeq2)))
options(width = 110)
section <- function(title) cat("\n==", title, "==\n")
near <- function(x, expected, digits) all(abs(unname(x) - expected) <= 0.5 * 10^-digits + 1e-9)  # x rounds to expected

section("1. Wald test for gene A: estimate divided by its standard error")
k  <- c(100L, 130L, 90L, 200L, 250L, 180L)
cd <- data.frame(condition = factor(rep(c("Ctrl", "Starvation"), each = 3)),
                 row.names = paste0(rep(c("Ctrl", "Starvation"), each = 3), "_", 1:3))
ddsA <- DESeqDataSetFromMatrix(matrix(k, nrow = 1, dimnames = list("geneA", rownames(cd))), cd, ~ condition)
sizeFactors(ddsA) <- rep(1, 6)
dispersions(ddsA) <- 0.053147                         # gene A's final dispersion (note 03)
ddsA <- nbinomWaldTest(ddsA, quiet = TRUE)
resA <- results(ddsA)
z <- resA$log2FoldChange / resA$lfcSE
ci <- resA$log2FoldChange + c(-1, 1) * qnorm(0.975) * resA$lfcSE
cat(sprintf("LFC %.5f  SE %.6f  Z %.5f  p %.6g | by hand: Z = LFC/SE = %.4f, p = 2*pnorm(-|Z|) = %.5g | 95%% CI [%.4f, %.4f]\n",
            resA$log2FoldChange, resA$lfcSE, resA$stat, resA$pvalue, z, 2 * pnorm(-abs(z)), ci[1], ci[2]))
stopifnot(near(resA$stat, 3.38093, 5), near(resA$pvalue, 0.000722408, 9),
          all.equal(resA$pvalue, 2 * pnorm(-abs(z))), near(ci, c(0.4107, 1.5438), 4))

section("2. Same fold change, different standard errors (study text 9.2)")
tab <- data.frame(log2FC = 1, SE = c(0.2, 0.5, 1.0))
tab$Z <- tab$log2FC / tab$SE
tab$p <- 2 * pnorm(-abs(tab$Z))
print(tab, digits = 4)
stopifnot(near(tab$p, c(5.733e-07, 0.0455, 0.3173), c(10, 5, 4)))
# The dispersion estimate changes the standard error and the p-value, not the fold change
for (a in c(0.014786, 0.025385, 0.053147)) {          # note 03: MLE, Cox-Reid, MAP
  d <- ddsA; dispersions(d) <- a
  r <- results(nbinomWaldTest(d, quiet = TRUE))
  cat(sprintf("alpha = %.6f   LFC = %.4f   SE = %.4f   Z = %.3f   p = %.2g\n", a, r$log2FoldChange, r$lfcSE, r$stat, r$pvalue))
}

section("3. Practice data: Ctrl / Starvation / Glucose, four donors, design ~ pair + condition")
# The true Starvation - Ctrl effect is 0 for every gene; Glucose - Ctrl has true effects.
set.seed(1)
dds0 <- makeExampleDESeqDataSet(n = 1000, m = 16, betaSD = 1)
keep <- c(1:8, 9:12)                                  # Ctrl = A1-4, Starvation = A5-8, Glucose = B1-4
cts  <- counts(dds0)[, keep]
meta <- data.frame(condition = factor(rep(c("Ctrl", "Starvation", "Glucose"), each = 4), levels = c("Ctrl", "Starvation", "Glucose")),
                   pair = factor(rep(1:4, 3)))
colnames(cts) <- rownames(meta) <- paste0(meta$condition, "_", meta$pair)
dds3 <- DESeqDataSetFromMatrix(cts, meta, ~ pair + condition)
dds_wald <- DESeq(dds3, quiet = TRUE)
res_wS <- results(dds_wald, contrast = c("condition", "Starvation", "Ctrl"))
res_wG <- results(dds_wald, name = "condition_Glucose_vs_Ctrl")
print(as.data.frame(res_wS[c("gene841", "gene420", "gene4"), ]), digits = 4)

# A nominal 95% CI that excludes 0 is the same event as raw p < 0.05; padj also counts the other genes
r  <- results(dds_wald, name = "condition_Glucose_vs_Ctrl", alpha = 0.05)
lo <- r$log2FoldChange - qnorm(0.975) * r$lfcSE; hi <- r$log2FoldChange + qnorm(0.975) * r$lfcSE
ci_excludes_0 <- lo > 0 | hi < 0
n_ci <- sum(ci_excludes_0, na.rm = TRUE); n_p <- sum(r$pvalue < 0.05, na.rm = TRUE); n_padj <- sum(r$padj < 0.05, na.rm = TRUE)
cat(sprintf("Glucose vs Ctrl: CI excludes 0 for %d genes, raw p < 0.05 for %d, padj < 0.05 for %d\n", n_ci, n_p, n_padj))
stopifnot(identical(which(ci_excludes_0), which(r$pvalue < 0.05)), c(n_ci, n_padj) == c(277, 163))

# "padj < 0.05 and |LFC| > 1" is a filter; testing |beta| > 1 needs lfcThreshold
r1 <- results(dds_wald, name = "condition_Glucose_vs_Ctrl", alpha = 0.05, lfcThreshold = 1, altHypothesis = "greaterAbs")
n_filter <- sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE); n_test <- sum(r1$padj < 0.05, na.rm = TRUE)
cat(sprintf("filter padj < 0.05 & |LFC| > 1: %d genes | test of |beta| > 1 (lfcThreshold = 1): %d genes\n", n_filter, n_test))
stopifnot(n_filter == 156, n_test == 18)

section("4. Likelihood ratio test")
# Gene A: drop the condition coefficient and see how much the log-likelihood falls
a <- 0.053147
l_full <- sum(dnbinom(k, mu = rep(c(mean(k[1:3]), mean(k[4:6])), each = 3), size = 1 / a, log = TRUE))
l_red  <- sum(dnbinom(k, mu = rep(mean(k), 6), size = 1 / a, log = TRUE))
resA_lrt <- results(nbinomLRT(ddsA, reduced = ~ 1, quiet = TRUE))
cat(sprintf("gene A: D = 2 (l_full - l_reduced) = %.4f | DESeq2 LRT stat %.4f, p %.4g | Wald Z^2 %.4f, p %.4g\n",
            2 * (l_full - l_red), resA_lrt$stat, resA_lrt$pvalue, resA$stat^2, resA$pvalue))
stopifnot(near(2 * (l_full - l_red), 11.2933, 4), near(resA_lrt$stat, 11.2933, 4))

# Three groups: ~ pair + condition against ~ pair, two coefficients at once (df = 2)
dds_lrt <- DESeq(dds3, test = "LRT", reduced = ~ pair, quiet = TRUE)
res_l <- results(dds_lrt)
df_lrt <- ncol(attr(dds_lrt, "modelMatrix")) - ncol(attr(dds_lrt, "reducedModelMatrix"))
cat("df =", df_lrt, "| dispersions shared with the Wald fit:", identical(dispersions(dds_wald), dispersions(dds_lrt)), "\n")
print(data.frame(Z_Starvation = res_wS[c("gene841", "gene420"), "stat"], Z_Glucose = res_wG[c("gene841", "gene420"), "stat"],
                 LRT_D = res_l[c("gene841", "gene420"), "stat"], row.names = c("gene841", "gene420")), digits = 4)
stopifnot(df_lrt == 2, identical(dispersions(dds_wald), dispersions(dds_lrt)), near(res_l["gene841", "stat"], 184.3, 1))

# The LRT p-value is an omnibus p-value: changing the displayed contrast does not change it
res_lS <- results(dds_lrt, contrast = c("condition", "Starvation", "Ctrl"))
cat("LRT p unchanged by contrast:", identical(res_l$pvalue, res_lS$pvalue), "|", mcols(res_l)$description[5], "\n")
sig <- function(x) !is.na(x$padj) & x$padj < 0.05
cat(sprintf("padj < 0.05: Wald Glucose-Ctrl %d, LRT %d, Wald only %d, LRT only %d\n",
            sum(sig(res_wG)), sum(sig(res_l)), sum(sig(res_wG) & !sig(res_l)), sum(sig(res_l) & !sig(res_wG))))
stopifnot(identical(res_l$pvalue, res_lS$pvalue), c(sum(sig(res_wG)), sum(sig(res_l))) == c(163, 146))

# Two groups (df = 1): Wald and LRT mostly agree
keepG <- meta$condition %in% c("Ctrl", "Glucose")
metaG <- droplevels(meta[keepG, , drop = FALSE])
ddsG  <- DESeq(DESeqDataSetFromMatrix(cts[, rownames(metaG)], metaG, ~ pair + condition), quiet = TRUE)
wG <- results(ddsG, name = "condition_Glucose_vs_Ctrl"); lG <- results(nbinomLRT(ddsG, reduced = ~ pair, quiet = TRUE))
ok <- !is.na(wG$pvalue) & !is.na(lG$pvalue)
cat(sprintf("Ctrl vs Glucose only: %d genes, cor(-log10 p) %.4f, p < 0.05 by Wald %d, by LRT %d\n", sum(ok),
            cor(-log10(wG$pvalue[ok]), -log10(lG$pvalue[ok])), sum(wG$pvalue[ok] < 0.05), sum(lG$pvalue[ok] < 0.05)))
stopifnot(sum(ok) == 995, sum(wG$pvalue[ok] < 0.05) == 255, sum(lG$pvalue[ok] < 0.05) == 229)

section("5. 'Significant in A, not in B' is not a difference")
set.seed(2)
d <- makeExampleDESeqDataSet(n = 500, m = 16)         # every true effect is 0
d$genotype  <- factor(rep(c("A", "B"), each = 8))
d$condition <- factor(rep(c("C", "T"), 8))
design(d) <- ~ genotype + condition + genotype:condition
dw  <- DESeq(d, quiet = TRUE)
rA  <- results(dw, name = "condition_T_vs_C")                                       # T - C within A
rB  <- results(dw, contrast = list(c("condition_T_vs_C", "genotypeB.conditionT")))   # T - C within B
int <- results(dw, name = "genotypeB.conditionT")                                    # the difference itself
one_side <- xor(rA$pvalue < 0.05, rB$pvalue < 0.05)
cat(sprintf("p < 0.05 in only one genotype: %d genes | interaction p < 0.05: %d of %d\n",
            sum(one_side, na.rm = TRUE), sum(int$pvalue < 0.05, na.rm = TRUE), sum(!is.na(int$pvalue))))
stopifnot(sum(one_side, na.rm = TRUE) == 60, sum(int$pvalue < 0.05, na.rm = TRUE) == 27,
          isTRUE(all.equal(int$log2FoldChange, rB$log2FoldChange - rA$log2FoldChange)))

section("6. Refitting Ctrl vs Starvation without the Glucose samples changes more than the number of groups")
keep2 <- meta$condition %in% c("Ctrl", "Starvation")
meta2 <- droplevels(meta[keep2, , drop = FALSE])
dds2  <- DESeq(DESeqDataSetFromMatrix(cts[, rownames(meta2), drop = FALSE], meta2, ~ pair + condition), quiet = TRUE)
r3 <- results(dds_wald, contrast = c("condition", "Starvation", "Ctrl"), alpha = 0.05)
r2 <- results(dds2, contrast = c("condition", "Starvation", "Ctrl"), alpha = 0.05)
cmp <- rbind(residual_df       = c(ncol(dds_wald) - ncol(attr(dds_wald, "dispModelMatrix")), ncol(dds2) - ncol(attr(dds2, "dispModelMatrix"))),
             median_dispersion = c(median(dispersions(dds_wald), na.rm = TRUE), median(dispersions(dds2), na.rm = TRUE)),
             median_lfcSE      = c(median(r3$lfcSE, na.rm = TRUE), median(r2$lfcSE, na.rm = TRUE)),
             baseMean_gene841  = c(r3["gene841", "baseMean"], r2["gene841", "baseMean"]))
colnames(cmp) <- c("three_groups", "two_groups")
print(noquote(format(signif(cmp, 4), drop0trailing = TRUE)))
stopifnot(cmp["residual_df", ] == c(6, 3), near(cmp["median_dispersion", ], c(0.4143, 0.5139), 4),
          near(cmp["baseMean_gene841", ], c(918.0, 174.8), 1))

cat("\nAll checks passed.\n")

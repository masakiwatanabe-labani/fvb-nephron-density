#!/usr/bin/env Rscript
## DESeq2 per stage on the hybrid matrix (C57BL/6J from --b6-arm, FVB/N from --fvb-arm).
## Same model and filters as run_deseq2.R (design ~ strain; rowSums >= 10; results(alpha = 0.1)).
## Optional: exclude genes on given contigs from the count matrix BEFORE any filtering.
suppressMessages(library(DESeq2))
## arguments: --key=value (no extra packages needed)
parse_kv <- function(defaults) {
  a <- commandArgs(trailingOnly = TRUE)
  for (x in a) { kv <- regmatches(x, regexpr("=", x), invert = TRUE)[[1]]; defaults[[sub("^--", "", kv[1])]] <- kv[2] }
  miss <- names(defaults)[vapply(defaults, is.null, logical(1))]
  if (length(miss)) stop("missing arguments: ", paste(miss, collapse = ", "))
  defaults }
opt <- parse_kv(list(matdir = NULL, `b6-arm` = "uncorrected", `fvb-arm` = "condB", `gene-chr` = NULL,
                     `exclude-chroms` = "", outdir = NULL, `count-suffix` = "_counts.tsv"))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
excl <- strsplit(opt$`exclude-chroms`, ",")[[1]]; excl <- excl[nzchar(excl)]
gchr <- read.table(opt$`gene-chr`, header = TRUE, sep = "\t", quote = "", comment.char = "")[, c("gene_id", "chr")]
rd <- function(arm) read.table(file.path(opt$matdir, paste0(arm, opt$`count-suffix`)), header = TRUE, sep = "\t",
                               row.names = 1, check.names = FALSE, quote = "", comment.char = "")
b6arm <- rd(opt$`b6-arm`); fvbarm <- rd(opt$`fvb-arm`)
groups <- list(E13.5 = list(B6 = c("B6135-1","B6135-2","B6135-3"), FVB = c("FVB135-1","FVB135-2","FVB135-3")),
               P1    = list(B6 = c("B6P1-1","B6P1-2","B6P1-3"),    FVB = c("FVBP1-1","FVBP1-2","FVBP1-3")))
log <- c()
for (tp in names(groups)) {
  g <- groups[[tp]]
  counts <- cbind(b6arm[, g$B6], fvbarm[rownames(b6arm), g$FVB])
  n0 <- nrow(counts); chr <- gchr$chr[match(rownames(counts), gchr$gene_id)]
  drop <- chr %in% excl
  excluded <- data.frame(gene_id = rownames(counts)[drop], chr = chr[drop])
  write.table(excluded, file.path(opt$outdir, paste0("excluded_genes_", tp, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  counts <- counts[!drop, ]
  meta <- data.frame(strain = factor(rep(c("B6","FVB"), each = 3), levels = c("B6","FVB")), row.names = colnames(counts))
  counts <- round(as.matrix(counts)); mode(counts) <- "integer"
  write.table(counts, file.path(opt$outdir, paste0("hybrid_counts_", tp, ".tsv")), sep = "\t", quote = FALSE)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ strain)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("strain", "FVB", "B6"), alpha = 0.1)
  res <- res[order(res$padj), ]
  out <- as.data.frame(res); out$gene_id <- rownames(out)
  write.table(out, file.path(opt$outdir, paste0("DESeq2_", tp, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(counts(dds, normalized = TRUE), file.path(opt$outdir, paste0("normcounts_", tp, ".tsv")), sep = "\t", quote = FALSE)
  write.table(data.frame(sample = colnames(dds), sizeFactor = sizeFactors(dds)), file.path(opt$outdir, paste0("sizeFactors_", tp, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  md <- metadata(res)
  line <- sprintf("%s\tgenes_in_matrix=%d\texcluded=%d (chrM=%d, chrY=%d)\tafter_rowSums>=10=%d\tnonNA_padj=%d\tpadj<0.1=%d\tpadj<0.1&|lfc|>0.5=%d\tfilterThreshold=%s",
                  tp, n0, sum(drop), sum(excluded$chr == "chrM"), sum(excluded$chr == "chrY"), nrow(dds), sum(!is.na(res$padj)),
                  sum(res$padj < 0.1, na.rm = TRUE), sum(res$padj < 0.1 & abs(res$log2FoldChange) > 0.5, na.rm = TRUE),
                  paste(signif(md$filterThreshold, 6), collapse = ","))
  cat(line, "\n"); log <- c(log, line)
}
writeLines(log, file.path(opt$outdir, "summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))

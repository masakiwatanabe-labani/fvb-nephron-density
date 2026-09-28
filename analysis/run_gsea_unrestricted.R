#!/usr/bin/env Rscript
# Unrestricted gene-set enrichment on the developing-kidney DESeq2 results.
#
# The submitted analysis tested only gene sets whose NAME matched a kidney regular expression
# (106 sets). That answers "is the nephrogenic programme altered?" but cannot answer "is anything
# else altered?" - in particular the postnatal-growth programmes that the phenotype (kidney mass per
# glomerulus) points at. This script tests the collections without a name filter.
#
# Genes flagged by the hybrid-design reverse control (D11: B6 read count shifts >2-fold when only the
# genome changes) are removed first, because the extreme tail of the DEG list is dominated by them.
#
# Usage: run_gsea_unrestricted.R --deseq-dir DIR --id2name FILE --d11-dir DIR --outdir DIR [--minsize 15]
#                                [--maxsize 500] [--samplesize 1001] [--seed 42]
suppressPackageStartupMessages({library(fgsea); library(msigdbr); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, d = NULL) { i <- which(args == paste0("--", k)); if (length(i)) args[i + 1] else d }
deseq_dir <- getarg("deseq-dir"); d11_dir <- getarg("d11-dir"); outdir <- getarg("outdir")
id2name <- getarg("id2name")
minsize <- as.integer(getarg("minsize", "15")); maxsize <- as.integer(getarg("maxsize", "500"))
samplesize <- as.integer(getarg("samplesize", "1001")); seed <- as.integer(getarg("seed", "42"))
stopifnot(!is.null(deseq_dir), !is.null(outdir), !is.null(id2name))
map <- fread(id2name, header = FALSE, col.names = c("gene_id", "gene_name"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- gene sets: no name filter ----
coll <- list(c("C5", "GO:BP"), c("C2", "CP:REACTOME"), c("C2", "CP:KEGG_MEDICUS"), c("H", NA))
sets <- list()
for (cc in coll) {
  m <- if (is.na(cc[2])) msigdbr(species = "Mus musculus", collection = cc[1])
       else msigdbr(species = "Mus musculus", collection = cc[1], subcollection = cc[2])
  if (!nrow(m)) next
  s <- split(m$gene_symbol, m$gs_name)
  sets <- c(sets, s)
  cat(sprintf("%s %s: %d sets\n", cc[1], ifelse(is.na(cc[2]), "", cc[2]), length(s)))
}
sets <- sets[!duplicated(names(sets))]
cat(sprintf("total gene sets loaded: %d\n", length(sets)))

# ---- D11 artefact genes (union over stages) ----
d11 <- character(0)
if (!is.null(d11_dir)) {
  for (f in list.files(d11_dir, pattern = "^screen_condB_(E13\\.5|P1)\\.tsv$", full.names = TRUE)) {
    t <- fread(f); d11 <- union(d11, t$Geneid)
  }
  cat(sprintf("D11-flagged genes excluded: %d\n", length(d11)))
}

res_all <- list()
for (tp in c("E13.5", "P1")) {
  d <- fread(file.path(deseq_dir, paste0("DESeq2_", tp, ".tsv")))
  d <- merge(d, map, by = "gene_id", all.x = TRUE)
  n0 <- nrow(d)
  d <- d[!(gene_id %in% d11)]
  d <- d[!is.na(stat) & !is.na(gene_name) & gene_name != ""]
  # one gene per symbol: keep the first in identifier order, as in the submitted analysis
  setorder(d, gene_name, gene_id)
  d <- d[!duplicated(gene_name)]
  r <- setNames(d$stat, d$gene_name)
  cat(sprintf("\n%s: %d genes in, %d after D11 removal and symbol collapse\n", tp, n0, length(r)))
  set.seed(seed)
  fg <- fgsea(pathways = sets, stats = r, minSize = minsize, maxSize = maxsize,
              eps = 0, sampleSize = samplesize)
  fg <- fg[order(padj)]
  fg[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
  fwrite(fg, file.path(outdir, paste0("gsea_unrestricted_", tp, ".tsv")), sep = "\t")
  cat(sprintf("%s: %d sets tested, min padj = %.4g, sets with padj<0.05 = %d, <0.1 = %d\n",
              tp, nrow(fg), min(fg$padj, na.rm = TRUE), sum(fg$padj < 0.05, na.rm = TRUE),
              sum(fg$padj < 0.1, na.rm = TRUE)))
  print(head(fg[, .(pathway, NES, padj, size)], 15))
  res_all[[tp]] <- fg
}
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("\nUNRESTRICTED_GSEA_DONE\n")

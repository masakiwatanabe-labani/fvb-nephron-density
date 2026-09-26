#!/usr/bin/env Rscript
## fgsea as in run_fgsea.R, with sampleSize as an argument (default 1001) and the corrected set-name regex as default; (ranks = DESeq2 stat; custom 4 sets + msigdbr C5 GO:BP + C2 CP:REACTOME filtered by
## a name regex; set.seed(42); fgsea(minSize = 3, maxSize = 500, eps = 0)). Paths and regex are arguments.
suppressMessages({ library(fgsea); library(msigdbr); library(dplyr) })
parse_kv <- function(defaults) {
  for (x in commandArgs(trailingOnly = TRUE)) { kv <- regmatches(x, regexpr("=", x), invert = TRUE)[[1]]; defaults[[sub("^--", "", kv[1])]] <- kv[2] }
  miss <- names(defaults)[vapply(defaults, is.null, logical(1))]; if (length(miss)) stop("missing: ", paste(miss, collapse = ", ")); defaults }
opt <- parse_kv(list(`deseq-dir` = NULL, id2name = NULL, outdir = NULL,
                     pattern = "(?i)kidney|nephro|(?<!ad)renal|metanephr|ureter|glomerul|podocyte", seed = "42", samplesize = "1001"))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
custom_sets <- list(
  nephron_progenitor = c("Six2","Cited1","Eya1","Meox1","Sall1","Wt1","Osr1"),
  ureteric_bud       = c("Ret","Wnt11","Calb1","Gfra1","Gata3","Wnt9b"),
  differentiation    = c("Lhx1","Pax8","Jag1","Hnf1b","Podxl","Nphs1","Nphs2"),
  stroma             = c("Foxd1","Meis1","Pdgfrb","Tcf21"))
mm_go <- suppressWarnings(msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP"))
mm_re <- suppressWarnings(msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME"))
go <- mm_go %>% filter(grepl(opt$pattern, gs_name, perl = TRUE)); re <- mm_re %>% filter(grepl(opt$pattern, gs_name, perl = TRUE))
all_sets <- c(custom_sets, split(go$gene_symbol, go$gs_name), split(re$gene_symbol, re$gs_name))
writeLines(c(sprintf("pattern\t%s", opt$pattern), sprintf("n_GO\t%d", length(unique(go$gs_name))), sprintf("n_Reactome\t%d", length(unique(re$gs_name))),
             sprintf("n_total\t%d", length(all_sets)), paste0("set\t", names(all_sets))), file.path(opt$outdir, "gene_sets.txt"))
id2name <- read.table(opt$id2name, header = FALSE, sep = "\t", col.names = c("gene_id", "gene_name"), quote = "", comment.char = "")
for (tp in c("E13.5", "P1")) {
  deres <- read.table(file.path(opt$`deseq-dir`, paste0("DESeq2_", tp, ".tsv")), header = TRUE, sep = "\t")
  deres <- merge(deres, id2name, by = "gene_id"); deres <- deres[!is.na(deres$stat), ]; deres <- deres[!duplicated(deres$gene_name), ]
  ranks <- sort(setNames(deres$stat, deres$gene_name), decreasing = TRUE)
  set.seed(as.integer(opt$seed))
  res <- fgsea(pathways = all_sets, stats = ranks, minSize = 3, maxSize = 500, eps = 0, sampleSize = as.integer(opt$samplesize))
  res <- res[order(res$pval), ]; res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ";")
  write.table(res, file.path(opt$outdir, paste0("fgsea_", tp, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%s: sets tested=%d, min padj=%.4f, n padj<0.25=%d, n padj<0.1=%d\n", tp, nrow(res), min(res$padj), sum(res$padj < 0.25), sum(res$padj < 0.1)))
}
writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))

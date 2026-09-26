## motifbreakR v3
## v2 からの変更点: 入力VCFを bcftools norm -m -any で biallelic に分解し、
## GT="alt" で FVB が実際に保有する ALT アレルのみを残したものを使う。
## (v1/v2 では多アレル部位で motifbreakR が最後の ALT を評価しており、
##  FVB が実際に持つアレルとは異なる予測が全体の約半数で生じていた)

suppressMessages({
  library(motifbreakR)
  library(BSgenome.Mmusculus.UCSC.mm39)
  library(MotifDb)
})

## arguments: --vcf=<decomposed, GT="alt"-filtered VCF> --outdir=<dir>
.kv <- list(vcf = NA, outdir = NA)
for (x in commandArgs(TRUE)) { kv <- regmatches(x, regexpr("=", x), invert = TRUE)[[1]]; .kv[[sub("^--", "", kv[1])]] <- kv[2] }
if (anyNA(unlist(.kv))) stop("usage: --vcf=... --outdir=...")
vcf_path <- .kv$vcf
outdir <- .kv$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

priority_tfs <- c("SIX2","WT1","PAX2","PAX8","HNF1B","OSR1","SALL1","TCF21","FOXD1","EYA1")

mm <- query(MotifDb, "Mmusculus")
gn <- toupper(values(mm)$geneSymbol)
motifs <- mm[gn %in% priority_tfs]
cat(sprintf("Matched %d priority-TF motifs (MotifDb): %s\n", length(motifs),
            paste(unique(values(motifs)$geneSymbol), collapse=", ")))
missing_tfs <- setdiff(priority_tfs, unique(gn))
if (length(missing_tfs) > 0) {
  cat(sprintf("No motif available in MotifDb for: %s (not curated/no direct DNA-binding motif)\n",
              paste(missing_tfs, collapse=", ")))
}

if (length(motifs) == 0) {
  stop("No priority TF motifs matched in MotifDb")
}

snps.mb <- snps.from.file(vcf_path, search.genome = BSgenome.Mmusculus.UCSC.mm39,
                           format = "vcf")
cat(sprintf("Loaded %d variants for motifbreakR\n", length(snps.mb)))

## sanity check: SNP_id 内の ALT と実際に評価される ALT が一致しているか
id_alt <- sub(".*/", "", as.character(snps.mb$SNP_id))
obs_alt <- as.character(snps.mb$ALT)
n_mismatch <- sum(id_alt != obs_alt)
cat(sprintf("ALT consistency check: %d / %d mismatched (expect 0)\n", n_mismatch, length(snps.mb)))
## hard stop: every record's identifier allele must equal the allele that will be scored
stopifnot("SNP_id allele differs from the scored ALT allele; decompose multi-allelic records first" = n_mismatch == 0)

results <- motifbreakR(snpList = snps.mb, pwmList = motifs,
                        threshold = 1e-3, method = "ic",
                        BPPARAM = BiocParallel::SerialParam())

cat(sprintf("motifbreakR found %d disruption events\n", length(results)))
saveRDS(results, file.path(outdir, "motifbreakr_results_raw.rds"))
cat("Saved raw GRanges to motifbreakr_results_raw.rds\n")

if (length(results) > 0) {
  df <- as.data.frame(results, row.names = NULL)
  df[] <- lapply(df, function(col) {
    if (is.list(col)) {
      vapply(col, function(x) paste(as.character(x), collapse = ";"), character(1))
    } else {
      col
    }
  })
  write.table(df, file.path(outdir, "motifbreakr_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Wrote motifbreakr_results.tsv\n")
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
} else {
  cat("No motif disruptions found at this threshold\n")
  write.table(data.frame(), file.path(outdir, "motifbreakr_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

## Per-gene log2 fold changes of the two pre-specified MSigDB Hallmark interferon sets,
## at E13.5 and P1. Replaces the 18-gene panel, which was assembled from the P1 leading edge.
suppressPackageStartupMessages({library(msigdbr); library(data.table)})
a <- commandArgs(TRUE); g <- function(k,d=NULL){i<-which(a==paste0("--",k)); if(length(i)) a[i+1] else d}
deseq <- g("deseq-dir"); id2 <- g("id2name"); out <- g("out")
m <- msigdbr(species="Mus musculus", collection="H")
sets <- split(m$gene_symbol, m$gs_name)
keep <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE")
map <- fread(id2, header=FALSE, col.names=c("gene_id","sym")); map <- map[!duplicated(gene_id)]
res <- list()
for (tp in c("E13.5","P1")) {
  d <- fread(file.path(deseq, paste0("DESeq2_", tp, ".tsv")))
  d <- merge(d, map, by="gene_id", all.x=TRUE)
  setorder(d, sym, gene_id); d <- d[!duplicated(sym) & !is.na(sym) & sym != ""]
  for (s in keep) {
    x <- d[sym %in% sets[[s]], .(set=s, timepoint=tp, gene=sym, log2FoldChange, baseMean, padj)]
    res[[length(res)+1]] <- x
  }
}
r <- rbindlist(res); fwrite(r, out, sep="\t")
cat(sprintf("%d rows -> %s\n", nrow(r), out))
print(r[, .(n_tested=sum(!is.na(log2FoldChange)), median_log2FC=round(median(log2FoldChange, na.rm=TRUE),4),
            median_baseMean=round(median(baseMean, na.rm=TRUE),0),
            frac_padj_lt_0.1=round(sum(padj<0.1,na.rm=TRUE)/sum(!is.na(padj)),3)), by=.(set,timepoint)])

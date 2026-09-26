## デコンボリューション v2 — v1失敗の修正版
##
## v1の失敗原因:
##   (1) NNLS: 特徴選択なしで全20,347遺伝子を使用した。腎の細胞型は転写産物の大半を共有する
##       ためシグネチャ行列の列が強く共線となり、解が1細胞型に退化した(全検体 Cap_Mesenchyme
##       100%)。CIBERSORT等が必ずマーカー遺伝子に限定するのはこのため。
##   (2) Bisque: 参照が単一ライブラリのため、擬似個体へのランダム分割で個体間分散を捏造して
##       いた。Bisqueは真の個体間分散から bulk 変換を推定するため、この操作では変換が退化し
##       全検体で同一の推定値を返した(差は小数第4位)。→ 本版では Bisque を使用しない。
##
## v2の設計:
##   - 細胞型特異的マーカー遺伝子による特徴選択
##   - **擬似バルクによる検証を必須化**: 既知の混合比を復元できるか確認してから実バルクに適用。
##     復元できない場合は実バルクの結果を採用しない。
##   - シグネチャ構築用と擬似バルク生成用に細胞を二分割し、循環参照を回避

.libPaths(c("/usr/local/jupyter/FVB_B6_glo/Rlib_deconv", .libPaths()))
suppressMessages({ library(nnls) })
set.seed(42)

## paths as --key=value arguments (defaults reproduce the original layout)
.kv <- list(base = "/usr/local/jupyter/FVB_B6_glo", sc = NA, out = NA, bulkdir = NA, id2name = NA)
for (x in commandArgs(TRUE)) { kv <- regmatches(x, regexpr("=", x), invert = TRUE)[[1]]; .kv[[sub("^--", "", kv[1])]] <- kv[2] }
BASE <- .kv$base
SC   <- if (is.na(.kv$sc)) file.path(BASE, "refs/scrna_kidney") else .kv$sc
OUT  <- .kv$out
BULKDIR <- .kv$bulkdir
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CELL_TYPES <- c("Cap_Mesenchyme","Distal_Tubule","Endothelial","Loop_of_Henle",
                "Nephron_Progenitor","Podocytes","Proximal_Tubule","Stromal","Ureteric_Bud")

## ---------- 参照の読み込み ----------
cat("参照 scRNA-seq (GSE94333 P1) 読み込み...\n")
dge <- read.table(gzfile(file.path(SC,"cold_dge.txt.gz")), header=TRUE, row.names=1,
                  sep="\t", check.names=FALSE)
labs <- list()
for (ct in CELL_TYPES) {
  bc <- trimws(readLines(gzfile(file.path(SC, sprintf("GSE94333_P1_Cold_%s.txt.gz", ct)))))
  labs[[ct]] <- bc[nzchar(bc)]
}
bc  <- unlist(labs, use.names=FALSE); ct <- rep(names(labs), lengths(labs))
keep <- bc %in% colnames(dge); bc <- bc[keep]; ct <- ct[keep]
sc <- as.matrix(dge[, bc, drop=FALSE])
cat(sprintf("  %d genes x %d labeled cells\n", nrow(sc), ncol(sc)))

## 細胞を2分割（A: シグネチャ用 / B: 擬似バルク用）
half <- unlist(lapply(CELL_TYPES, function(k){
  i <- which(ct==k); s <- sample(i); ifelse(seq_along(s) %% 2 == 0, "A","B") -> h; setNames(h, s)
}))
split_lab <- rep("A", length(ct))
for (k in CELL_TYPES) { i <- which(ct==k); s <- sample(i); split_lab[s[seq(2, length(s), by=2)]] <- "B" }

cpm <- function(m) sweep(m, 2, pmax(colSums(m),1), "/") * 1e6

## ---------- マーカー遺伝子選択 ----------
## 各細胞型について「その型での平均発現 / 全型平均」の比が高く、かつ十分発現する遺伝子を選ぶ
pick_markers <- function(mat, ctv, n_per_type = 60) {
  prof <- sapply(CELL_TYPES, function(k) rowMeans(cpm(mat)[, ctv==k, drop=FALSE]))
  markers <- unique(unlist(lapply(CELL_TYPES, function(k){
    spec <- (prof[,k] + 1) / (rowMeans(prof) + 1)
    ok <- prof[,k] > 20                     # 十分に発現している
    g <- names(sort(spec[ok], decreasing=TRUE))[seq_len(min(n_per_type, sum(ok)))]
    g[!is.na(g)]
  })))
  list(markers = markers, profile = prof)
}

idxA <- which(split_lab=="A"); idxB <- which(split_lab=="B")
mk <- pick_markers(sc[, idxA, drop=FALSE], ct[idxA])
cat(sprintf("  マーカー遺伝子: %d\n", length(mk$markers)))
SIG <- mk$profile[mk$markers, CELL_TYPES, drop=FALSE]

solve_props <- function(bulk_cpm, sig) {
  g <- intersect(rownames(sig), rownames(bulk_cpm))
  S <- sig[g,,drop=FALSE]; Y <- bulk_cpm[g,,drop=FALSE]
  ## 遺伝子ごとにスケールを揃える（高発現遺伝子の支配を防ぐ）
  sc_f <- apply(S, 1, max); sc_f[sc_f<=0] <- 1
  S <- S / sc_f; Y <- Y / sc_f
  t(apply(Y, 2, function(y){
    co <- nnls(S, y)$x
    if (sum(co) > 0) co/sum(co) else rep(NA_real_, ncol(S))
  })) -> P
  colnames(P) <- colnames(sig); P
}

## ---------- 擬似バルク検証（必須ゲート） ----------
cat("\n=== 擬似バルク検証 (既知混合比を復元できるか) ===\n")
make_pseudo <- function(props, cells_idx, n_cells = 2000) {
  pick <- unlist(lapply(CELL_TYPES, function(k){
    pool <- cells_idx[ct[cells_idx]==k]
    n <- round(props[k]*n_cells)
    if (n<=0 || length(pool)==0) return(integer(0))
    sample(pool, n, replace=TRUE)
  }))
  rowSums(sc[, pick, drop=FALSE])
}
truth <- rbind(
  rep(1/9, 9),
  c(.05,.05,.05,.05,.45,.05,.10,.10,.10),
  c(.15,.05,.05,.10,.25,.05,.15,.10,.10),
  c(.10,.10,.05,.10,.15,.05,.20,.15,.10)
)
colnames(truth) <- CELL_TYPES
pb <- sapply(seq_len(nrow(truth)), function(i) make_pseudo(truth[i,], idxB))
colnames(pb) <- paste0("pseudo", seq_len(ncol(pb)))
est <- solve_props(cpm(pb), SIG)
cat("真の混合比:\n"); print(round(truth,3))
cat("推定値:\n");   print(round(est,3))
err <- abs(est - truth)
r_all <- cor(as.vector(est), as.vector(truth))
cat(sprintf("\n全体相関 r = %.3f / 平均絶対誤差 = %.3f / 最大誤差 = %.3f\n",
            r_all, mean(err), max(err)))
npc_err <- abs(est[,"Nephron_Progenitor"] - truth[,"Nephron_Progenitor"])
cat(sprintf("Nephron_Progenitor の絶対誤差: %s\n", paste(round(npc_err,3), collapse=", ")))

GATE <- (r_all > 0.8) && (mean(err) < 0.05)
cat(sprintf("\n検証ゲート (r>0.8 かつ 平均絶対誤差<0.05): %s\n", ifelse(GATE,"合格","不合格")))
write.csv(cbind(truth=truth, est=est), file.path(OUT,"pseudobulk_validation.csv"))

if (!GATE) {
  cat("\n*** 検証不合格のため、実バルクへの適用は行わない。 ***\n")
  cat("*** この参照・手法の組み合わせでは細胞組成を定量できないと結論する。 ***\n")
  quit(save="no", status=0)
}

## ---------- 実バルクへの適用 ----------
id2name <- read.table(.kv$id2name, sep="\t", header=FALSE,
                      col.names=c("gene_id","gene_name"), quote="", comment.char="")
for (tp in c("P1","E13.5")) {
  cat(sprintf("\n===== %s =====\n", tp))
  if (tp=="E13.5") cat("  注意: 参照はP1由来。E13.5には成熟尿細管がほぼ存在せず参照が不一致。\n")
  bulk <- read.table(file.path(BULKDIR, sprintf("hybrid_counts_%s.tsv", tp)),
                     header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  meta <- data.frame(strain = ifelse(grepl("^B6", colnames(bulk)), "B6", "FVB"), row.names = colnames(bulk))
  sym <- id2name$gene_name[match(rownames(bulk), id2name$gene_id)]
  ok <- !is.na(sym) & !duplicated(sym)
  bulk <- bulk[ok,,drop=FALSE]; rownames(bulk) <- sym[ok]
  P <- solve_props(cpm(as.matrix(bulk)), SIG)
  print(round(P,4))
  strain <- meta[rownames(P),"strain"]
  stats <- t(sapply(colnames(P), function(k){
    b <- P[strain=="B6",k]; f <- P[strain=="FVB",k]
    c(B6=mean(b), FVB=mean(f), diff=mean(f)-mean(b),
      p=tryCatch(t.test(f,b)$p.value, error=function(e) NA_real_))
  }))
  cat("\n群間比較:\n"); print(round(stats,4))
  write.csv(data.frame(sample=rownames(P), strain=strain, P, check.names=FALSE),
            file.path(OUT, sprintf("proportions_v2_%s.csv", tp)), row.names=FALSE)
  write.csv(stats, file.path(OUT, sprintf("stats_v2_%s.csv", tp)))
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("\n完了\n")

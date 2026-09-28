# Code for the FVB/N nephron-density manuscript

Analysis code for

> **Recovered glomerular counts relative to kidney mass distinguish FVB/N from three other inbred
> mouse strains: phenotypic characterisation and bias-corrected developmental transcriptome
> analysis**
> (submitted manuscript, *International Journal of Molecular Sciences*)

This repository holds the scripts that produced the figures, tables and numbers reported in the
submitted manuscript. It does not hold raw or intermediate data.

---

## What is here

| Directory | Contents |
|---|---|
| `pipeline/` | Trimming, library-type inference, alignment, personalised-genome construction, quantification, differential expression |
| `g2gtools_patch/` | The patched g2gtools module and the wrapper used for reverse BAM conversion, with a minimal reproducible example |
| `analysis/` | Mapping-bias QC, variant extraction and annotation, candidate prioritisation, allele-dosage concordance, gene-set enrichment, deconvolution, the liver replication, table builders |
| `figures/` | One script per figure. Figure 1 marks significance with Tukey HSD across all four strains, matching the P values quoted in the text; Figure 7B shows the two predefined MSigDB Hallmark interferon sets in full rather than a post-hoc selection of genes |
| `source_data/` | The aggregated tables actually plotted, one per figure or table panel |
| `environment/` | Conda environment exports, R `sessionInfo()` per analysis, tool versions, reference-data provenance |
| `MANIFEST.tsv` | Figure or table → script → output, with the path of each script in the original analysis directory and its md5 |

## Reading the MANIFEST

`MANIFEST.tsv` has one row per script:

| column | meaning |
|---|---|
| `category` | which directory it lives in here |
| `target` | the figure, table or quantity it produces |
| `source_path_in_analysis_dir` | where the file sits in the original analysis directory |
| `repo_path` | where it sits here |
| `outputs` | the files it writes |
| `md5` | md5 of the file as it stands here |
| `md5_matches_source` | whether that md5 equals the md5 of the file in the original analysis directory |

Every file was copied byte-for-byte and the md5 column was verified against the original after
copying. Four files were subsequently edited here for wording only — `analysis/isg_predefined.R`,
`analysis/make_figure7B_data.R` and `figures/plot_figure7.py`, where a comment said
"pre-specified" instead of "predefined", and `analysis/run_gsea_unrestricted.R`, where a comment
said "renal mass" instead of "kidney mass". Their `md5_matches_source` is `False`; no executable
line was changed, so they still reproduce the same output as the files that were run.

## Important: paths are not rewritten

**The scripts take their inputs from the directory layout of the original analysis and those paths
were deliberately left untouched**, so that what is published is exactly what was run. Most scripts
take their paths as arguments; some of the earlier pipeline scripts contain an absolute path to the
analysis root. Running them elsewhere means supplying the same inputs, not editing the scripts.

## Input data

| Input | Where to get it |
|---|---|
| RNA-seq reads and processed count matrices | Gene Expression Omnibus, accession *(to be filled in at submission)* |
| Mouse reference genome and annotation | GENCODE release M39 (GRCm39). File names, URLs, download dates and md5 checksums are in `environment/reference_data.tsv` |
| Strain variants | Mouse Genomes Project release REL2021 |
| Gene sets | MSigDB release 2026.1, retrieved at run time with the R package **msigdbr 26.1.0** (`msigdbr(species = "Mus musculus", collection = ...)`). **No gene-set file is redistributed here**: the scripts fetch the sets themselves, so that the KEGG-derived collections, whose redistribution MSigDB restricts, are never copied into this repository. Human sets are mapped to mouse orthologues by msigdbr's built-in mapping (babelgene 22.9) |
| Kidney ChIP-seq and H3K27ac peaks | GEO accessions listed in `environment/reference_data.tsv` |
| Single-cell reference for deconvolution | GEO GSE94333 |
| Liver replication dataset | GEO GSE167328 |

`environment/reference_data.tsv` holds the provenance of every reference file used here: its
exact name, source URL, download date and md5. For the full analysis settings and the complete
reference-data table, see Supplementary Table S1 of the paper; no draft of that table is kept in
this repository.

## Order of execution

1. `pipeline/run_fastp*.sh` — trimming
2. `pipeline/run_salmon_*.sh` — library type inferred from the data (`--libType A`)
3. `pipeline/run_phase1_variants_v2.sh`, `run_phase2_norm.sh` — FVB/N variants
4. `pipeline/run_phase2_cond{A,B}.sh` — personalised genomes with g2gtools
5. `pipeline/run_star_*_index.sh` — STAR indices
6. `pipeline/run_phase2_s2pairs.sh` — alignment, reverse coordinate conversion and quantification of all three arms, the final settings (`featureCounts -p --countReadPairs -s 2`)
7. `pipeline/build_matrix_param.py` — count matrices
8. `pipeline/run_deseq2_param.R` — differential expression on the hybrid matrix
9. `analysis/*` — everything downstream
10. `figures/*` — figures

Steps 1–8 need the sequencing data; steps 9–10 can be re-run from the outputs of step 8 together
with the tables in `source_data/`.

## Environment

```bash
conda env create -n fvb_pipeline -f environment/conda_fvb_pipeline.yml
conda env create -n genomics     -f environment/conda_genomics.yml
```

`environment/tool_versions.tsv` lists the version of every external tool, and
`environment/sessionInfo/` the R session for each analysis that used R.

## The g2gtools wrapper

g2gtools 0.2.7 cannot convert alignments back from a personalised genome to reference coordinates:
`bsam.convert_bam_file` never populates the `@SQ` header or the contig-name lookup, so every read
fails. A second problem appears once that is fixed: contigs declared in the VCI header but carrying
no variant records are absent from the coordinate mapping tree, and every alignment on them is
written to the `.unmapped` file. For FVB/NJ these are chrM, chrY and 39 unplaced scaffolds.

`g2gtools_patch/` contains:

| File | What it is |
|---|---|
| `bsam_fixed.py` | the patched module, rebuilding `@SQ` from the VCI `##CONTIG` lines |
| `g2g_bam_convert_v2.py` | the wrapper, which additionally inserts an identity interval for every contig declared in the VCI header but holding no records |
| `upstream_0.2.7.diff` | unified diff against the pristine upstream module |
| `UPSTREAM_LICENSE` | the upstream MIT licence, which governs these files |
| `example/` | a minimal reproducible example |

Run the example with

```bash
bash g2gtools_patch/example/make_and_run_mre_v2.sh
```

It builds a small synthetic genome with one contig carrying variants and one without, converts a
BAM both ways, and shows that the unpatched path loses every read on the variant-free contig while
the wrapper keeps them. `example/expected_output.log` is the output we obtain; the run ends with
`ALL_OK` and an empty `.unmapped` file.

The conversion is valid only with `--reverse`, because the rebuilt `@SQ` lengths are those of the
reference genome.

## Licence

The code written for this study is under the MIT Licence (`LICENSE`). The files under
`g2gtools_patch/` derive from g2gtools 0.2.7, which is itself MIT-licensed; that licence is
reproduced in `g2gtools_patch/UPSTREAM_LICENSE` and governs those files.

## Citation

See `CITATION.cff`. The journal reference will be added when the paper is published.

## Repository

https://github.com/masakiwatanabe-labani/fvb-nephron-density

Tag `v1.0` marks the state of the code used for the submitted manuscript.

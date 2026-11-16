#!/usr/bin/env Rscript
# GSE145926 Figure 7c,d displayed cells -> Supplementary Data 6.
# The original, separate Control and ALI computational/selection methods are
# retained. The only post-selection label changes are the approved Control
# B cell -> B cell-like and Neutrophil -> Unresolved. No Figure 7e cells added.
# Usage: Rscript 01_GEO_to_SuppData6.R GEO_H5_DIR NEW_OUTPUT_DIR

repair_xlsx <- function(src,dst) {
if(!requireNamespace('zip',quietly=TRUE) || !requireNamespace('openxlsx',quietly=TRUE))stop('Packages zip and openxlsx are required.')
src<-normalizePath(src,winslash='/',mustWork=TRUE)
dst<-file.path(normalizePath(dirname(dst),winslash='/',mustWork=TRUE),basename(dst))
if(file.exists(dst))stop('Refusing to overwrite an existing workbook: ',dst)
entries<-utils::unzip(src,list=TRUE)$Name
if(any(grepl('(^/|^[A-Za-z]:|(^|/)[.][.](/|$))',entries)))stop('Unsafe path in workbook archive.')
tmp<-tempfile('fig7e_ooxml_');dir.create(tmp);on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
utils::unzip(src,exdir=tmp)
read_text<-function(path){f<-file(path,'rb');on.exit(close(f));rawToChar(readBin(f,'raw',n=file.info(path)$size))}
write_text<-function(path,text){f<-file(path,'wb');on.exit(close(f));writeBin(charToRaw(text),f)}
sheet_names<-grep('^xl/worksheets/sheet[0-9]+[.]xml$',entries,value=TRUE)
if(length(sheet_names)!=2L)stop('Unexpected number of worksheets: ',length(sheet_names))
changes<-list()
for(sheet in sheet_names) {
 sp<-file.path(tmp,sheet);xml<-read_text(sp)
 refs<-regmatches(xml,gregexpr('<c r="[A-Z]+[0-9]+"',xml,perl=TRUE))[[1]]
 if(!length(refs))stop('Worksheet has no cell references: ',sheet)
 refs<-sub('^<c r="','',refs);refs<-sub('"$','',refs)
 letters<-sub('[0-9]+$','',refs);nums<-as.integer(sub('^[A-Z]+','',refs)); u<-unique(letters)
 value<-vapply(strsplit(u,''),function(ch)sum((match(ch,LETTERS))*26^rev(seq_along(ch)-1)),numeric(1))
 maxcol<-u[which.max(value)];new_dim<-paste0('A1:',maxcol,max(nums));old_dim<-regmatches(xml,regexpr('<dimension ref="[A-Z0-9:]+"[[:space:]]*/>',xml,perl=TRUE))
 if(length(old_dim)!=1L || !nzchar(old_dim))stop('Worksheet dimension not found: ',sheet)
 fixed_dim<-paste0('<dimension ref="',new_dim,'"/>')
 xml2<-sub(old_dim,fixed_dim,xml,fixed=TRUE);stopifnot(identical(sub(fixed_dim,old_dim,xml2,fixed=TRUE),xml))
 write_text(sp,xml2)
 rel<-paste0('xl/worksheets/_rels/',basename(sheet),'.rels');removed<-character()
 if(rel %in% entries) {
  rp<-file.path(tmp,rel);rels<-read_text(rp);links<-regmatches(rels,gregexpr('<Relationship [^>]+/>',rels,perl=TRUE))[[1]]
  for(link in links) {
   if(!grepl('/(drawing|vmlDrawing)"',link,perl=TRUE))next
   id<-sub('.* Id="([^"]+)".*','\\1',link);target<-sub('.* Target="([^"]+)".*','\\1',link)
   if(identical(id,link)||identical(target,link))stop('Could not parse drawing relationship: ',rel)
   if(!file.exists(file.path(dirname(sp),target))) {
    if(grepl(paste0('r:id="',id,'"'),xml,fixed=TRUE))stop('Missing drawing is used by worksheet: ',sheet,' ',id)
    rels<-sub(link,'',rels,fixed=TRUE);removed<-c(removed,id)
   }
  }
  write_text(rp,rels)
 }
 changes[[length(changes)+1L]]<-data.frame(sheet=sheet,dimension=new_dim,removed_unused_drawing_relationships=paste(removed,collapse=';'))
}
out_entries<-list.files(tmp,all.files=TRUE,recursive=TRUE,include.dirs=FALSE)
if(!setequal(entries[!grepl('/$',entries)],out_entries))stop('Archive entry list changed during repair.')
zip::zipr(dst,files=out_entries,root=tmp,include_directories=FALSE,mode='mirror')
if(!file.exists(dst))stop('Repair output was not created.')
# Compare every visible worksheet value in both files, including headers.
sheets<-openxlsx::getSheetNames(src);stopifnot(identical(sheets,openxlsx::getSheetNames(dst)))
for(s in sheets){a<-openxlsx::read.xlsx(src,sheet=s,colNames=FALSE,skipEmptyRows=FALSE,skipEmptyCols=FALSE);b<-openxlsx::read.xlsx(dst,sheet=s,colNames=FALSE,skipEmptyRows=FALSE,skipEmptyCols=FALSE);if(!identical(a,b))stop('Worksheet values changed during repair: ',s)}
print(do.call(rbind,changes),row.names=FALSE);cat('Workbook structure repaired; all worksheet values identical.\n')

}

run_control_supp <- function(geo_dir,work) {
  CTRL_ROOT<-file.path(work,'control');dir.create(CTRL_ROOT,recursive=TRUE)
  FINAL_RESULTS<-CTRL_ROOT
RULE_LINES <- c(
  "gene,marker_group,detection_threshold,interpretation,evidence_status",
  "MARCO,AM_core,raw_UMI_count > 0,\"AM evidence, lineage positive >=2; cluster context required\",DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "FABP4,AM_core,raw_UMI_count > 0,\"AM evidence, lineage positive >=2; cluster context required\",DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "C1QA,AM_core,raw_UMI_count > 0,\"AM evidence, lineage positive >=2; cluster context required\",DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "PPARG,AM_support,raw_UMI_count > 0,Context only; cannot alone assign AM,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CD3D,T_core,raw_UMI_count > 0,\"T lineage; 2 of 3 supports, macrophage conflict excluded\",DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CD3E,T_core,raw_UMI_count > 0,\"T lineage; 2 of 3 supports, macrophage conflict excluded\",DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "TRAC,T_core,raw_UMI_count > 0,\"T lineage; 2 of 3 supports, macrophage conflict excluded\",DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "TRBC1,T_support,raw_UMI_count > 0,\"Supportive, not required for primary identity\",DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "IL7R,T_support,raw_UMI_count > 0,\"Supportive, not required for primary identity\",DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "LTB,T_support,raw_UMI_count > 0,\"Supportive, not required for primary identity\",DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "MS4A1,B_primary,raw_UMI_count > 0,>=2 detected and low AM/T conflict; B-like candidate only,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CD79A,B_primary,raw_UMI_count > 0,>=2 detected and low AM/T conflict; B-like candidate only,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CD79B,B_primary,raw_UMI_count > 0,>=2 detected and low AM/T conflict; B-like candidate only,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CD19,B_primary,raw_UMI_count > 0,>=2 detected and low AM/T conflict; B-like candidate only,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "PAX5,B_primary,raw_UMI_count > 0,>=2 detected and low AM/T conflict; B-like candidate only,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "BANK1,B_support_specific,raw_UMI_count > 0,Supportive only; not decisive without B primary markers,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "VPREB3,B_support_specific,raw_UMI_count > 0,Supportive only; not decisive without B primary markers,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CD37,B_support_nonspecific,raw_UMI_count > 0,NON-SPECIFIC: never use by itself to declare B cell,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CD74,B_support_nonspecific,raw_UMI_count > 0,NON-SPECIFIC: never use by itself to declare B cell,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "HLA-DRA,B_support_nonspecific,raw_UMI_count > 0,NON-SPECIFIC: never use by itself to declare B cell,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "FCGR3B,NEUT_specific,raw_UMI_count > 0,Specific lineage support; at least one with other evidence,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CEACAM8,NEUT_specific,raw_UMI_count > 0,Specific lineage support; at least one with other evidence,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CXCR2,NEUT_secondary,raw_UMI_count > 0,Support only; S100A8/S100A9/CSF3R alone not sufficient,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CSF3R,NEUT_secondary,raw_UMI_count > 0,Support only; S100A8/S100A9/CSF3R alone not sufficient,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "S100A8,NEUT_secondary,raw_UMI_count > 0,Support only; S100A8/S100A9/CSF3R alone not sufficient,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "S100A9,NEUT_secondary,raw_UMI_count > 0,Support only; S100A8/S100A9/CSF3R alone not sufficient,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "FPR1,NEUT_secondary,raw_UMI_count > 0,Support only; S100A8/S100A9/CSF3R alone not sufficient,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CTSG,NEUT_secondary,raw_UMI_count > 0,Support only; S100A8/S100A9/CSF3R alone not sufficient,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "MPO,NEUT_secondary,raw_UMI_count > 0,Support only; S100A8/S100A9/CSF3R alone not sufficient,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "ELANE,NEUT_secondary,raw_UMI_count > 0,Support only; S100A8/S100A9/CSF3R alone not sufficient,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CCL3,target_expression,raw_UMI_count > 0,raw UMI >0 = detection; log-normalized values for plotted intensity,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED",
  "CCL4,target_expression,raw_UMI_count > 0,raw UMI >0 = detection; log-normalized values for plotted intensity,DISCOVERY_RULE_NOT_AUTHOR_VALIDATED"
)
writeLines(RULE_LINES,file.path(FINAL_RESULTS,"ANNOTATION_RULES_FROZEN_FOR_REANALYSIS.csv"),useBytes=TRUE)
rm(RULE_LINES)

# Input choices (optional, not required):
# Sys.setenv(CONTROL_DOWNLOAD_MISSING="false")  # stop instead of GEO download
# Final release disables QC cache: full raw H5 read on each execution.
# Default ALWAYS runs H5->QC on fresh run, and redoes normalization/PCA/UMAP.
cat("\n========== STEP 1/2: ORIGINAL H5 -> QC -> PCA -> UMAP ==========\n")

# ---------- AUTHORITATIVE RAW-H5 ALI_v1.3 STEP ----------
# Rebuild the 3-sample unsupervised embedding with ALI v1.3 computation
# (NOT the same cell-annotation strategy as ALI).
# Original GEO GSE145926 C51 C52 C100 only; original H5 checksums enforced.
# DO NOT ALTER this file to make the scatter plot look like the manuscript.
# This step computes UMAP once and caches the exact Seurat object for step 02.
# ============================================================================
# CONTROL GSE145926 C51 + C52 + C100 | ALI v1.3 METHOD-MATCHED REPRODUCTION
# Reference: GSE145926_Figure7cd_GEO_ONLY_REPRODUCE.R (ALI v1.3)
# Single self-contained script: original GEO H5 -> matched QC/PCA/UMAP ->
# cluster-based Figure subset -> Figure/Source Data/validation reports.
#
# SAME ALI COMPUTATIONAL PROCEDURE, changing cohort (9 ALI -> 3 Control),
# input checksums, and independently reviewed CONTROL cluster->type map.
# NOTE: ALI's cluster numbers (0/3/4/6) are NOT transferable to Control.
# Reviewed Control map supports candidate macrophage-like/T only.
# Clusters failing the predeclared marker gate are excluded and explicitly reported.
# B and neutrophil clusters are NOT supported by Control marker evidence.
# This script will NOT falsely label any cluster just to draw four colors.
# Historical manuscript UMAP coordinates are not an input or an output claim.
#
# Strict execution needs the ALI v1.3 package versions (especially Seurat 5.4.0).
# ============================================================================
options(stringsAsFactors = FALSE, timeout = max(1800L, getOption("timeout")))

ROOT <- CTRL_ROOT
IN_DIR <- geo_dir
OUT <- FINAL_RESULTS
DOWNLOAD_MISSING <- tolower(Sys.getenv("CONTROL_DOWNLOAD_MISSING", "true")) %in% c("1", "yes", "true")
REUSE_QC <- FALSE # Final release: ALWAYS re-read original H5 and redo QC
REQUIRE_ALI_RUNTIME <- tolower(Sys.getenv("CONTROL_REQUIRE_ALI_RUNTIME", "true")) %in% c("1", "yes", "true")
stopifnot(normalizePath(IN_DIR, mustWork = FALSE) != normalizePath(OUT, mustWork = FALSE))
dir.create(IN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
STATUS <- file.path(OUT, "99_STATUS.txt")
writeLines("STARTED: GSE145926 C51 C52 C100; not a completed result", STATUS)
fail <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
wc <- function(x, name) utils::write.csv(x, file.path(OUT, name), row.names = FALSE, na = "")
logmsg <- function(...) { cat(..., "\n", sep = "") }

# ONLY the three Control samples in GSE145926. MD5 hashes were observed in
# the earlier three-H5 GEO audit, not inferred from 28,056-cell Source Data.
manifest <- data.frame(
  sample = c("C51", "C52", "C100"),
  gsm = c("GSM4475048", "GSM4475049", "GSM4475050"),
  h5 = c("GSM4475048_C51_filtered_feature_bc_matrix.h5",
         "GSM4475049_C52_filtered_feature_bc_matrix.h5",
         "GSM4475050_C100_filtered_feature_bc_matrix.h5"),
  expected_md5 = c("7ba33981a06a212d0b4cd8625b65eca8",
                   "d6532007fcb7691073158617ad79cbb0",
                   "8f2dfb95d992608b0ec1400619b687f5"),
  expected_qc = c(11112L, 10361L, 6583L),
  stringsAsFactors = FALSE)
geo_url <- function(gsm, h5) paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/",
  sub("[0-9]{3}$", "nnn", gsm), "/", gsm, "/suppl/", h5)
find_input <- function(h5,gsm) {
  p <- file.path(geo_dir,h5)
  if (!file.exists(p)) stop("Missing original GEO H5: ",h5)
  p
}
ali_versions <- c(R="4.5.2",Seurat="5.4.0",SeuratObject="5.2.0",
                  uwot="0.2.4",RSpectra="0.16.2",RcppAnnoy="0.0.22")
required <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "hdf5r",
              "uwot", "RSpectra", "RcppAnnoy")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
actual_versions <- vapply(names(ali_versions), function(p) {
  if (p == "R") return(as.character(getRversion()))
  as.character(utils::packageVersion(p))
}, character(1))
version_audit <- data.frame(component = names(ali_versions),
  ALI_v13 = unname(ali_versions), Control_runtime = unname(actual_versions),
  matches = unname(actual_versions) == unname(ali_versions))
wc(version_audit, "00_ALI_RUNTIME_COMPARISON.csv")
writeLines(capture.output(sessionInfo()), file.path(OUT, "00_SESSION_INFO.txt"))
if (!all(version_audit$matches)) {
  warn <- paste(version_audit$component[!version_audit$matches], collapse = ", ")
  if (REQUIRE_ALI_RUNTIME) {
    writeLines(paste("STOP: ALI v1.3 environment differs:", warn), STATUS)
    stop("EXACT ALI ENVIRONMENT REQUIRED; mismatched components: ", warn,
         "\nDo not replace your ALI environment. Use a separate R library/project with",
         " the ALI versions listed in 00_ALI_RUNTIME_COMPARISON.csv.",
         "\nFor exploratory-only calculations set CONTROL_REQUIRE_ALI_RUNTIME=false;",
         " that run must not be called an exact reproduction.", call. = FALSE)
  } else warning("EXPLORATORY runtime differs from ALI: ", warn)
}
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})
options(Seurat.object.assay.version = "v5")
set.seed(42)

# Prevent local precomputed results from becoming an undisclosed input.
manifest$filepath <- mapply(find_input, manifest$h5, manifest$gsm, USE.NAMES = FALSE)
manifest$observed_md5 <- unname(tools::md5sum(manifest$filepath))
manifest$md5_pass <- tolower(manifest$observed_md5) == tolower(manifest$expected_md5)
wc(manifest, "00_GEO_INPUT_MANIFEST.csv")
fail(all(manifest$md5_pass), "GEO H5 MD5 does not match the frozen audit: STOP")
logmsg("[PASS] Three original GEO control H5 inputs and MD5 checks")

# ALI code: CreateSeuratObject(min.features=200,min.cells=0), nCount>=500,
# no mitochondrial exclusion; merge raw counts, not normalized values.
checkpoint <- file.path(OUT, "CACHE_QC_COUNTS_ONLY.rds")
checkpoint_metadata <- file.path(OUT, "CACHE_QC_COUNTS_ONLY_META.rds")
signature <- list(h5 = manifest$h5, md5 = manifest$observed_md5,
  min_features = 200L, min_counts = 500L, mitochondrial_filter = FALSE,
  seurat = as.character(packageVersion("Seurat")),
  seurat_object = as.character(packageVersion("SeuratObject")),
  procedure = "ALI_v1.3_CreateSeuratObject_merge")
cache_ok <- REUSE_QC && file.exists(checkpoint) && file.exists(checkpoint_metadata) &&
  identical(tryCatch(readRDS(checkpoint_metadata), error = function(e) NULL), signature)
if (cache_ok) {
  obj <- readRDS(checkpoint)
  qc <- read.csv(file.path(OUT, "02_PER_SAMPLE_QC.csv"), check.names = FALSE)
  fail(inherits(obj, "Seurat") && nrow(qc) == 3L, "QC cache is invalid")
  logmsg("[CACHE] Verified original GEO counts-only checkpoint")
} else {
  obj <- NULL
  qclist <- vector("list", 3L)
  for (i in seq_len(nrow(manifest))) {
    sid <- manifest$sample[[i]]
    logmsg("[GEO] ", sid)
    mat <- Seurat::Read10X_h5(manifest$filepath[[i]])
    if (is.list(mat)) {
      fail("Gene Expression" %in% names(mat), paste("Missing Gene Expression assay:", sid))
      mat <- mat[["Gene Expression"]]
    }
    raw_bc <- colnames(mat)
    pre <- ncol(mat)
    so <- Seurat::CreateSeuratObject(counts = mat, project = sid,
                                     min.features = 200L, min.cells = 0L)
    rm(mat)
    so$sample_id <- sid
    so$GSM <- manifest$gsm[[i]]
    so$barcode_raw <- colnames(so)
    so <- subset(so, subset = nFeature_RNA >= 200 & nCount_RNA >= 500)
    qclist[[i]] <- data.frame(sample = sid, gsm = manifest$gsm[[i]],
      input_matrix_columns = pre, qc_cells = ncol(so), removed = pre - ncol(so))
    fail(ncol(so) == manifest$expected_qc[[i]],
         paste("ALI QC cell-count gate differs for", sid, ":", ncol(so)))
    if (is.null(obj)) obj <- so else obj <- merge(obj, y = so, merge.data = FALSE)
    rm(so, raw_bc); invisible(gc())
  }
  qc <- do.call(rbind, qclist)
  wc(qc, "02_PER_SAMPLE_QC.csv")
  if (REUSE_QC) {
    tmp <- paste0(checkpoint, ".partial")
    saveRDS(obj, tmp, compress = FALSE)
    if (file.exists(checkpoint)) unlink(checkpoint)
    fail(file.rename(tmp, checkpoint), "Cannot finalize QC checkpoint")
    saveRDS(signature, checkpoint_metadata)
  }
}
fail(ncol(obj) == 28056L && !anyDuplicated(colnames(obj)), "QC cell membership incorrect")
obs <- as.integer(table(factor(obj$sample_id, levels = manifest$sample)))
fail(identical(obs, as.integer(manifest$expected_qc)), "Per-sample QC membership incorrect")
logmsg("[QC] 28,056 cells; ", nrow(obj), " genes")

# ALI code: Join COUNTS only; split per sample BEFORE NormalizeData and VST.
DefaultAssay(obj) <- "RNA"
count_layers <- grep("^counts(\\.|$)", Layers(obj[["RNA"]]), value = TRUE)
if (length(count_layers) > 1L)
  obj <- JoinLayers(obj, assay = "RNA", layers = "counts")
obj[["RNA"]] <- split(obj[["RNA"]],
  f = factor(obj$sample_id, levels = manifest$sample), layers = "counts")
fail(length(grep("^counts(\\.|$)", Layers(obj[["RNA"]]), value = TRUE)) == 3L,
     "RNA count-layer count != three")

# ALI v1.3 computational calls and arguments: intentionally identical.
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = 3000, loess.span = 0.3,
                            verbose = FALSE)
hvgs <- VariableFeatures(obj)
fail(length(hvgs) == 3000L, "Variable-feature gate: expected 3000")
writeLines(hvgs, file.path(OUT, "03_HVG3000.txt"))
obj <- ScaleData(obj, features = hvgs, scale.max = 10, verbose = FALSE)
obj <- RunPCA(obj, features = hvgs, npcs = 40, seed.use = 42,
              verbose = FALSE)
pca <- Embeddings(obj, reduction = "pca")
fail(nrow(pca) == 28056L && ncol(pca) >= 40L, "PCA40 gate failed")
saveRDS(pca[, 1:40, drop = FALSE], file.path(OUT, "04_PCA40.rds"))
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30,
                     k.param = 20, nn.method = "annoy", n.trees = 50,
                     annoy.metric = "euclidean", verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, algorithm = 1,
                    random.seed = 0, verbose = FALSE)
suppressPackageStartupMessages(library(RSpectra))
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30,
               n.neighbors = 30L, n.components = 2L,
               metric = "cosine", min.dist = 0.3, seed.use = 42,
               uwot.sgd = FALSE, angular.rp.forest = FALSE,
               verbose = FALSE)

# Frozen post-UMAP result: subsequent cell typing/figures do not rerun clustering.
cache <- file.path(OUT, "CACHE_POST_UMAP_ALI_V13_SEURAT.rds")
cache_meta <- file.path(OUT, "CACHE_POST_UMAP_ALI_V13_META.csv")
stopifnot(nrow(Embeddings(obj, "umap")) == 28056L)
tmp_cache <- paste0(cache,".partial")
if (file.exists(tmp_cache)) unlink(tmp_cache)
saveRDS(obj, tmp_cache, compress="gzip")
stopifnot(file.exists(tmp_cache), file.size(tmp_cache)>1000000)
if (file.exists(cache)) unlink(cache)
if (!file.rename(tmp_cache,cache)) stop("Cannot finalize post-UMAP cache",call.=FALSE)
write.csv(data.frame(sample=manifest$sample, gsm=manifest$gsm,
                     h5=manifest$h5, md5=manifest$observed_md5,
                     runtime_R=as.character(getRversion()),
                     runtime_Seurat=as.character(packageVersion("Seurat")),
                     runtime_SeuratObject=as.character(packageVersion("SeuratObject")),
                     runtime_uwot=as.character(packageVersion("uwot")),
                     method="ALI_v1.3_3controls_joint_UMAP",
                     stringsAsFactors=FALSE), cache_meta, row.names=FALSE)
write.csv(data.frame(cell_barcode=colnames(obj), sample_id=as.character(obj$sample_id),
                     gsm=as.character(obj$GSM), barcode_raw=as.character(obj$barcode_raw),
                     cluster_id=as.character(obj$seurat_clusters),
                     UMAP_1=Embeddings(obj,"umap")[colnames(obj),1],
                     UMAP_2=Embeddings(obj,"umap")[colnames(obj),2]),
          file.path(OUT,"01_ALL_28056_UNSUPERVISED_EMBEDDING.csv"), row.names=FALSE)
writeLines("ALI_v1.3_COMPUTATION_COMPLETE; post-UMAP Seurat cache saved; annotations not yet checked",
           STATUS)
cat("[PASS] Post-UMAP Seurat object cached. Next run step 02 annotation.\n")
rm(obj); invisible(gc())

# ---------- RAW-DERIVED MARKER ANNOTATION AND PRIMARY SOURCE DATA ----------
# ============================================================================
# HEALTHY BALF Control | GSE145926 C51/C52/C100 | PROVISIONAL annotation
# Uses saved 28,056-cell ALI-v1.3-matched Seurat object (NO re-UMAP here).
# Marker sets are PROPOSED and FROZEN for this run in ANNOTATION_RULES_FROZEN_FOR_REANALYSIS.csv.
# This is an independent NEW annotation of healthy cells, not a recreation of
# historical C51-derived manuscript coordinates or its 7,285-cell selection.
# ---------------------------------------------------------------------------
# Safety: no cell identities are fabricated to match an expected figure;
#          zeros stay zeros, missing genes are not treated as expression=0;
#          ambiguous/unassigned cells stay in the full Source Data.
# ============================================================================
options(stringsAsFactors = FALSE)
root <- CTRL_ROOT
output <- FINAL_RESULTS
cache <- file.path(output, "CACHE_POST_UMAP_ALI_V13_SEURAT.rds")
cache_meta <- file.path(output, "CACHE_POST_UMAP_ALI_V13_META.csv")
marker_file <- file.path(FINAL_RESULTS, "ANNOTATION_RULES_FROZEN_FOR_REANALYSIS.csv")
status_file <- file.path(output, "99_STATUS.txt")
stop0 <- function(...) stop(paste0(...), call. = FALSE)
if (!file.exists(cache) || !file.exists(cache_meta)) {
  stop0("MISSING ALI_v1.3_UMAP_CACHE. Run 01_GEO_TO_ALI_V13_UMAP.R once first.\n", cache)
}
if (!file.exists(marker_file)) stop0("Missing FROZEN marker rules file: ", marker_file)
require_p <- c("Seurat","SeuratObject","ggplot2","Matrix","patchwork")
for (p in require_p) if (!requireNamespace(p,quietly=TRUE)) stop0("Missing required R package: ",p)
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(ggplot2); library(patchwork)})
expected_pkg <- c(Seurat="5.4.0", SeuratObject="5.2.0", uwot="0.2.4")
for (p in names(expected_pkg)) {
  cur <- as.character(utils::packageVersion(p))
  if (cur != expected_pkg[[p]]) stop0("Package mismatch for ",p, ": ",cur,
      "; expected ",expected_pkg[[p]]," (the audited ALI v1.3 environment).")
}
write_csv <- function(x, name) utils::write.csv(x, file.path(output,name), row.names=FALSE, na="")
writeLines("ANNOTATION_STARTED; NOT FOR PUBLICATION YET", status_file)

# Verify input GEO checksum identities and recorded software versions.
input_meta <- read.csv(cache_meta, check.names=FALSE, stringsAsFactors=FALSE)
if (!identical(as.character(input_meta$sample),c("C51","C52","C100")))
  stop0("Bad cache sample IDs: only C51,C52,C100 are allowed.")
if (!all(input_meta$runtime_R == "4.5.2") ||
    !all(input_meta$runtime_Seurat == "5.4.0") ||
    !all(input_meta$runtime_SeuratObject == "5.2.0") ||
    !all(input_meta$runtime_uwot == "0.2.4")) stop0("Cache was not computed under ALI v1.3 runtime")
manifest <- read.csv(file.path(output,"00_GEO_INPUT_MANIFEST.csv"), stringsAsFactors=FALSE)
if (nrow(manifest)!=3L || !all(manifest$md5_pass) ||
    !identical(as.character(manifest$sample),c("C51","C52","C100")))
  stop0("Original GEO input manifest missing, wrong, or failed MD5")
if (!identical(as.character(manifest$observed_md5),as.character(input_meta$md5)))
  stop0("GEO H5 checksum metadata and saved Seurat cache metadata disagree")
if (!all(file.exists(manifest$filepath))) stop0("Original GEO input H5 files no longer exist")
if (!all(unname(tools::md5sum(manifest$filepath)) == manifest$observed_md5))
  stop0("Raw GEO H5 changed since the saved UMAP was computed")

obj <- readRDS(cache)
if (!inherits(obj,"Seurat") || ncol(obj)!=28056L ||
    !identical(as.integer(table(factor(obj$sample_id,levels=c("C51","C52","C100")))),
               c(11112L,10361L,6583L)) ||
    !"umap" %in% Seurat::Reductions(obj)) stop0("Seurat cache failed integrity gate")
ids <- colnames(obj)
if (anyDuplicated(ids)) stop0("Duplicated Seurat barcodes")
if (anyNA(obj$barcode_raw) || anyNA(obj$sample_id)) stop0("Missing original barcode/sample mapping")
original_barcode <- as.character(obj$barcode_raw)
if (anyDuplicated(paste(obj$sample_id,original_barcode,sep="__")))
  stop0("Non-unique sample-scoped original barcodes")

rules <- read.csv(marker_file, stringsAsFactors=FALSE, check.names=FALSE, fileEncoding="UTF-8")
required_cols <- c("gene","marker_group","detection_threshold","interpretation","evidence_status")
if (!all(required_cols %in% colnames(rules)) || anyDuplicated(rules$gene))
  stop0("Rule file malformed or contains duplicated gene definitions")
expected_groups <- c("AM_core","AM_support","T_core","T_support","B_primary",
 "B_support_specific","B_support_nonspecific","NEUT_specific","NEUT_secondary","target_expression")
if (!all(expected_groups %in% unique(rules$marker_group))) stop0("Missing required marker groups")
get_genes <- function(group) rules$gene[rules$marker_group == group]
all_genes <- unique(as.character(rules$gene))
missing_genes <- setdiff(all_genes,rownames(obj))
write_csv(data.frame(gene=all_genes, present_in_matrix=all_genes %in% rownames(obj),
                     group=rules$marker_group[match(all_genes,rules$gene)]),
          "02_MARKER_AVAILABILITY.csv")
required_bio <- c("MARCO","FABP4","C1QA","CD3D","CD3E","TRAC",
                  "MS4A1","CD79A","CD79B","FCGR3B","CEACAM8",
                  "S100A8","S100A9","CCL3","CCL4")
if (any(!required_bio %in% rownames(obj))) stop0("Missing required genes in raw GEO matrix: ",
                                      paste(setdiff(required_bio,rownames(obj)),collapse=", "))
# Snapshot the rule file that actually generated the Source Data.
file.copy(marker_file,file.path(output,"ANNOTATION_RULES_APPLIED.csv"),overwrite=TRUE)
writeLines(unname(tools::md5sum(marker_file)),file.path(output,"02_APPLIED_RULE_MD5.txt"))

# Per-sample Seurat v5 layers are explicitly traversed. Never silently drop
# cells using FetchData on a single layer or join mismatched rows by position.
extract_layer <- function(prefix, genes) {
  lyr <- grep(paste0("^",prefix,"(\\.|$)"),SeuratObject::Layers(obj[["RNA"]]),value=TRUE)
  if (!length(lyr)) stop0("Missing RNA ",prefix," layers")
  ans <- matrix(NA_real_, nrow=length(genes),ncol=length(ids),dimnames=list(genes,ids))
  filled <- rep(FALSE,length(ids))
  avail <- intersect(genes,rownames(obj))
  for (one in lyr) {
    m <- SeuratObject::LayerData(obj[["RNA"]], layer=one)
    ix <- match(colnames(m),ids)
    if (anyNA(ix) || any(filled[ix])) stop0("Overlapping or unknown cells in layer: ",one)
    if (length(avail)) {
      sub <- as.matrix(m[avail,,drop=FALSE])
      ans[avail,ix] <- sub
    }
    filled[ix] <- TRUE
    rm(m)
  }
  if (!all(filled)) stop0("Some cells missing ",prefix," layer values: ",sum(!filled))
  if (anyNA(ans[avail,,drop=FALSE])) stop0("Missing present-gene data from layer ",prefix)
  ans
}
raw <- extract_layer("counts",all_genes)
expr <- extract_layer("data",c("CCL3","CCL4"))
if (!all(is.finite(expr))) stop0("Non-finite CCL3/CCL4 normalized expression")

score_detected <- function(group) {
  gs <- intersect(get_genes(group),rownames(raw))
  gs <- intersect(gs,rownames(obj))
  if (!length(gs)) return(rep(0L,length(ids)))
  as.integer(colSums(raw[gs,,drop=FALSE] > 0, na.rm=TRUE))
}
pos <- function(gene) if (gene %in% rownames(obj)) as.logical(raw[gene,] > 0) else rep(FALSE,length(ids))
am_n <- score_detected("AM_core")
t_n <- score_detected("T_core")
b_n <- score_detected("B_primary")
b_spec <- score_detected("B_support_specific")
b_non <- score_detected("B_support_nonspecific")
n_spec <- score_detected("NEUT_specific")
n_second <- score_detected("NEUT_secondary")

clusters <- as.character(obj$seurat_clusters)
cl_levels <- sort(unique(clusters))
if (all(grepl("^[0-9]+$",cl_levels))) cl_levels <- cl_levels[order(as.integer(cl_levels))]
# Identify AM and T cluster CONTEXT without borrowing ALI's arbitrary cluster IDs.
cl_frac <- do.call(rbind,lapply(cl_levels, function(cl) {
  sel <- clusters == cl
  gs <- c("MARCO","FABP4","C1QA","CD3D","CD3E","TRAC",
          "MS4A1","CD79A","CD79B","FCGR3B","CEACAM8","S100A8","S100A9")
  data.frame(cluster_id=cl,n_cells=sum(sel),gene=gs,
             pct_detected=100*vapply(gs,function(g)mean(pos(g)[sel]),numeric(1)),
             stringsAsFactors=FALSE)
}))
write_csv(cl_frac,"03_CLUSTER_MARKER_FRACTIONS.csv")
fraction <- function(cl,g) {
  z <- cl_frac$pct_detected[cl_frac$cluster_id==cl & cl_frac$gene==g]
  if (length(z)!=1L) stop0("Internal cluster fraction error ",cl," ",g)
  z
}
cluster_type <- setNames(rep("Unassigned",length(cl_levels)),cl_levels)
context <- do.call(rbind,lapply(cl_levels,function(cl) {
  am_support <- fraction(cl,"MARCO")>=60 && fraction(cl,"FABP4")>=60 &&
                fraction(cl,"C1QA")>=80
  t_support <- fraction(cl,"CD3D")>=25 && fraction(cl,"TRAC")>=25 &&
               fraction(cl,"C1QA")<60
  label <- if(am_support && !t_support) "AM-like" else
           if(t_support && !am_support) "T cell" else "Unassigned"
  data.frame(cluster_id=cl,n_cells=sum(clusters==cl),
             MARCO_pct=fraction(cl,"MARCO"),FABP4_pct=fraction(cl,"FABP4"),
             C1QA_pct=fraction(cl,"C1QA"),CD3D_pct=fraction(cl,"CD3D"),
             TRAC_pct=fraction(cl,"TRAC"),
             AM_cluster_gate=am_support,T_cluster_gate=t_support,
             cluster_context_label=label)
}))
cluster_type <- stats::setNames(context$cluster_context_label, context$cluster_id)
write_csv(context,"03_CLUSTER_CONTEXT_DECISIONS.csv")

# Strict lineage gates: B and neutrophil must have positive, complementary
# evidence and be incompatible with AM/T lineage. B-support CD74/HLA-DRA/CD37
# is deliberately NON-SPECIFIC and never enough by itself to claim B identity.
is_b <- (b_n>=2L) & (am_n<=1L) & (t_n<=1L) & (n_spec==0L)
is_neut <- (n_spec>=1L) & (n_second>=1L) & (am_n<=1L) &
           (t_n<=1L) & (b_n<=1L)
# Prevent double classification in case of a multi-lineage/doublet signal.
conflict <- is_b & is_neut
is_b[conflict] <- FALSE
is_neut[conflict] <- FALSE
# We also record weaker rare-lineage signals as NOT confidently annotated.
b_review <- (b_n>=1L & (b_spec>=1L | b_non>=2L) & am_n<=1L & t_n<=1L & !is_b)
neut_review <- (n_second>=3L & am_n<=1L & t_n<=1L & !is_neut)
ann <- rep("Unassigned",length(ids))
cluster_ctx <- unname(cluster_type[clusters])
ann[cluster_ctx=="AM-like" & am_n>=1L & t_n<=1L & !is_b & !is_neut] <- "Alveolar macrophage-like"
ann[cluster_ctx=="T cell" & t_n>=1L & am_n<=1L & !is_b & !is_neut] <- "T cell"
ann[is_b] <- "B cell"
ann[is_neut] <- "Neutrophil"
ann[conflict] <- "Unassigned"
confidence <- ifelse(ann=="Unassigned","UNRESOLVED", "RULE_SUPPORTED_PROVISIONAL")
confidence[ann=="B cell" | ann=="Neutrophil"] <- "RARE_LINEAGE_CANDIDATE_REQUIRES_EXPERT_REVIEW"
confidence[conflict] <- "MULTILINEAGE_CONFLICT"
confidence[ann=="Unassigned" & (b_review | neut_review)] <- "RARE_LINEAGE_REVIEW_ONLY"

um <- Seurat::Embeddings(obj,"umap")
all <- data.frame(cell_barcode=ids,
     sample_id=as.character(obj$sample_id),gsm=as.character(obj$GSM),
     barcode_raw=original_barcode,cluster_id=clusters,
     UMAP_1=as.numeric(um[ids,1]), UMAP_2=as.numeric(um[ids,2]),
     cluster_context_label=cluster_ctx,
     annotation_final=ann,annotation_status=confidence,
     display_in_figure=ann!="Unassigned",
     AM_primary_count=am_n,T_primary_count=t_n,
     B_primary_count=b_n,B_specific_support_count=b_spec,
     B_nonspecific_support_count=b_non,
     Neut_specific_count=n_spec,Neut_support_count=n_second,
     B_review_only=b_review,Neut_review_only=neut_review,
     multilineage_conflict=conflict,
     CCL3_raw=as.numeric(raw["CCL3",]),CCL4_raw=as.numeric(raw["CCL4",]),
     CCL3_expr=as.numeric(expr["CCL3",]),CCL4_expr=as.numeric(expr["CCL4",]),
     stringsAsFactors=FALSE)
for (g in all_genes) {
  if (g %in% c("CCL3","CCL4")) next
  all[[paste0("marker_",g,"_raw")]] <- as.numeric(raw[g,])
}
all$CCL3_positive <- all$CCL3_raw>0
all$CCL4_positive <- all$CCL4_raw>0
all$CCL3_CCL4_both_positive <- all$CCL3_positive & all$CCL4_positive
if (nrow(all)!=28056L || anyDuplicated(all$cell_barcode) ||
    any(!is.finite(all$UMAP_1)) || any(!is.finite(all$UMAP_2)))
  stop0("Invalid all-cell output")

fig <- all[all$display_in_figure,,drop=FALSE]
if (!nrow(fig)) stop0("No marker-supported cells, so no meaningful figure can be generated")
# ALI v1.3 plotting normalization: q05 and q95 on Figure-selected cells.
# When 95% quantile is zero (very sparse data), flag it instead of inventing
# a positive intensity range. Such scores stay zero; expression plot is raw
# log-normalized for the respective gene, not a falsely nonzero normalized map.
qscore <- function(z) {
  if (!length(z)) return(list(x=numeric(0),q05=NA_real_,q95=NA_real_,status="EMPTY"))
  q <- as.numeric(quantile(z,c(.05,.95),na.rm=TRUE,names=FALSE))
  if (q[2]<=q[1]) return(list(x=rep(0,length(z)),q05=q[1],q95=q[2],status="ZERO_DYNAMIC_RANGE"))
  list(x=pmin(pmax((z-q[1])/(q[2]-q[1]),0),1),q05=q[1],q95=q[2],status="OK")
}
s3 <- qscore(fig$CCL3_expr);s4 <- qscore(fig$CCL4_expr)
fig$CCL3_scaled_0to1 <- s3$x
fig$CCL4_scaled_0to1 <- s4$x
fig$CCL3_CCL4_product <- s3$x*s4$x
for (name in c("CCL3_scaled_0to1","CCL4_scaled_0to1","CCL3_CCL4_product"))
  all[[name]] <- NA_real_
ix <- match(fig$cell_barcode,all$cell_barcode)
for (name in c("CCL3_scaled_0to1","CCL4_scaled_0to1","CCL3_CCL4_product")) all[[name]][ix] <- fig[[name]]
write_csv(all,"Supplementary_Source_Data_ALL_28056.csv")
write_csv(fig,"Figure7cd_DISPLAYED_CELL_SOURCE.csv")
write_csv(data.frame(gene=c("CCL3","CCL4"),q05=c(s3$q05,s4$q05),
                     q95=c(s3$q95,s4$q95),status=c(s3$status,s4$status),
                     scope="Figure-selected cells only"),"04_CCL3_CCL4_QUANTILE_AUDIT.csv")
write_csv(as.data.frame(table(sample=all$sample_id,
                              annotation=all$annotation_final)),
          "05_ANNOTATION_BY_SAMPLE.csv")
write_csv(as.data.frame(table(cluster=all$cluster_id,
                              annotation=all$annotation_final)),
          "05_ANNOTATION_BY_CLUSTER.csv")

# Ensure source data are self-consistent after writing: no off-by-row or
# barcode/UMAP/expression mismatches between the all-cell and displayed files.
all_disk <- read.csv(file.path(output,"Supplementary_Source_Data_ALL_28056.csv"),
                     check.names=FALSE,stringsAsFactors=FALSE)
fig_disk <- read.csv(file.path(output,"Figure7cd_DISPLAYED_CELL_SOURCE.csv"),
                     check.names=FALSE,stringsAsFactors=FALSE)
match_ix <- match(fig_disk$cell_barcode,all_disk$cell_barcode)
if (anyNA(match_ix) || nrow(fig_disk)!=sum(all_disk$display_in_figure) ||
    anyDuplicated(fig_disk$cell_barcode) ||
    !identical(as.character(fig_disk$cell_barcode),as.character(all_disk$cell_barcode[match_ix])) ||
    max(abs(fig_disk$UMAP_1-all_disk$UMAP_1[match_ix]),na.rm=TRUE)>1e-10 ||
    max(abs(fig_disk$CCL3_raw-all_disk$CCL3_raw[match_ix]),na.rm=TRUE)>1e-10 ||
    max(abs(fig_disk$CCL4_raw-all_disk$CCL4_raw[match_ix]),na.rm=TRUE)>1e-10)
  stop0("Figure-vs-full-cell source data barcode/UMAP/expression audit failed")


  fig$annotation_final[fig$annotation_final=='B cell'] <- 'B cell-like'
  fig$annotation_final[fig$annotation_final=='Neutrophil'] <- 'Unresolved'
  stopifnot(nrow(fig)==20327L,sum(fig$annotation_final=='B cell-like')==14L,sum(fig$annotation_final=='Unresolved')==12L)
  control_plot <- data.frame(barcode=fig$cell_barcode,
    sample=paste0(fig$gsm,'_',fig$sample_id,'_filtered_feature_bc_matrix'),
    celltype=fig$annotation_final, UMAP_1=fig$UMAP_1, UMAP_2=fig$UMAP_2,
    CCL3_expr=fig$CCL3_expr, CCL4_expr=fig$CCL4_expr,
    CCL3_scaled_0to1=fig$CCL3_scaled_0to1, CCL4_scaled_0to1=fig$CCL4_scaled_0to1,
    CCL3_CCL4_coexpression_score=fig$CCL3_CCL4_product, stringsAsFactors=FALSE)
  saveRDS(control_plot,file.path(work,'control_plot_data.rds'))
  data.frame(barcode=control_plot$barcode,sample=control_plot$sample,CellType=control_plot$celltype,
    CCL3_expr=control_plot$CCL3_expr,CCL4_expr=control_plot$CCL4_expr,
    CCL3_CCL4_coexpr=control_plot$CCL3_CCL4_coexpression_score,stringsAsFactors=FALSE)
}

run_ali_supp <- function(geo_dir,work) {
  in_dir<-geo_dir;out_dir<-file.path(work,'ali');dir.create(out_dir,recursive=TRUE);run_mode<-'full';allow_download<-FALSE
  assert<-function(cond,msg)if(!isTRUE(cond))stop(msg,call.=FALSE)
manifest <- data.frame(
  gsm = c("GSM4339769", "GSM4339770", "GSM4339771", "GSM4339772", "GSM4339773",
          "GSM4339774", "GSM4475051", "GSM4475052", "GSM4475053"),
  donor = c("C141", "C142", "C143", "C144", "C145", "C146", "C148", "C149", "C152"),
  expected_cells = c(4233L, 4562L, 20289L, 915L, 17396L, 3401L, 2397L, 2485L, 6962L),
  expected_md5 = c("d21a4a65983ac7407dede9918a146bdd",
                   "f3ed7fda2ee390063fdbbe2ccc3cd927",
                   "ee95b9984bec80d4f8c67fe42a5bf408",
                   "5d278893ddad4c96a6cb8b6b6cd261bb",
                   "953c90404ecbd1f69b60a79c2a9ca8d1",
                   "591810707452f2378d7410449bb7cdce",
                   "6255440e85a6b8da8c647dfc6f3641d9",
                   "4a6c6304fe51ece1ccec9dccdce2079c",
                   "e1791ff67a806db99193dd8b4511f166"),
  stringsAsFactors = FALSE
)
# IMPORTANT: The manifest's MD5 values are verification references only.
# If a later source audit detects an incorrect manifest entry, stop and fix it;
# do not change a reference to force a mismatch to pass.
manifest$sample_id <- paste0(manifest$gsm, "_", manifest$donor, "_filtered_feature_bc_matrix")
manifest$filename <- paste0(manifest$sample_id, ".h5")
manifest$filepath <- file.path(in_dir, manifest$filename)

if(any(!file.exists(manifest$filepath)))stop("Missing original GEO H5: ",paste(manifest$filename[!file.exists(manifest$filepath)],collapse=", "))
manifest$observed_md5 <- unname(tools::md5sum(manifest$filepath))
manifest$md5_pass <- tolower(manifest$observed_md5) == tolower(manifest$expected_md5)
write.csv(manifest, file.path(out_dir, "00_GEO_H5_MANIFEST_AND_MD5.csv"), row.names = FALSE)
if (!all(manifest$md5_pass)) {
  bad <- paste(manifest$filename[!manifest$md5_pass], collapse = "\n  ")
  stop("PRE-FLIGHT FAIL: GEO H5 MD5 mismatch:\n  ", bad,
       "\nExpected hashes are frozen from previously audited original H5 files.\n",
       "Check the manifest against the original source before continuing.", call. = FALSE)
}
assert(sum(manifest$expected_cells) == 62640L, "Internal manifest error: cell sum")
message("[PASS] Nine GEO H5 files and MD5 checks. Expected 62640 cells.")
if (run_mode == "preflight") {
  message("[STOP BY DESIGN] Preflight complete. No Seurat analysis was run.")
  message("Only after reviewing this check, run with ALI_RUN_MODE=full.")
} else {

  ## 2 | PACKAGE PREFLIGHT ----------------------------------------------------
  # Dependencies are installed by researchers BEFORE running this script.
  # No auto-installation; do not alter the researcher's R environment silently.
  required_packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2",
                         "hdf5r", "uwot", "RSpectra", "RcppAnnoy")
  missing_packages <- required_packages[!vapply(
    required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
         "\nInstall the missing packages using install.packages(), ",
         "then rerun this script. See README.md.\n",
         "In a renv project, install packages into the project library ",
         "or use a standard RStudio session.", call. = FALSE)
  }
  for (p in required_packages)
    assert(requireNamespace(p, quietly = TRUE), paste("Missing R package:", p))
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(ggplot2)
  })
  vers <- vapply(c("R", "Seurat", "SeuratObject", "uwot", "RSpectra", "RcppAnnoy"),
    function(pkg) {
      if (pkg == "R") {
        as.character(getRversion())
      } else if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(packageVersion(pkg))
      } else {
        "NOT_INSTALLED"
      }
    }, character(1))
  write.csv(data.frame(package = names(vers), version = unname(vers)),
            file.path(out_dir, "01_RUNTIME_VERSIONS.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()),
             file.path(out_dir, "01_R_SESSION_INFO.txt"))
  message("[INFO] Versions and sessionInfo saved. Historical exact UMAP is NOT validated.")
  expected_runtime <- c(R = "4.5.2", Seurat = "5.4.0",
                        SeuratObject = "5.2.0", uwot = "0.2.4",
                        RSpectra = "0.16.2", RcppAnnoy = "0.0.22")
  different_versions <- names(vers)[vers != expected_runtime[names(vers)]]
  if (length(different_versions)) {
    warning("Runtime versions differ from audited environment: ",
            paste(different_versions, collapse = ", "),
            ". Reproduction gates will still be checked; consult 01_RUNTIME_VERSIONS.csv.")
  }
  assert(startsWith(as.character(utils::packageVersion("Seurat")), "5."),
         "Requires Seurat 5.x; do not silently change assay behavior")

  ## 3 | READ EACH H5, QC BEFORE ANY PCA ------------------------------------
  # Durable counts-only checkpoint: preserves the slow H5/QC/merge stage.
  # Its companion metadata binds it to these nine verified H5 MD5s and versions.
  reuse_checkpoint <- FALSE
  so_all <- NULL
  qc <- vector("list", nrow(manifest))
  {
  for (i in seq_len(nrow(manifest))) {
    sid <- manifest$sample_id[i]
    message("[H5] Reading ", i, "/9: ", sid)
    m <- Seurat::Read10X_h5(manifest$filepath[i])
    if (is.list(m)) {
      assert("Gene Expression" %in% names(m), paste("No Gene Expression assay:", sid))
      m <- m[["Gene Expression"]]
    }
    so <- Seurat::CreateSeuratObject(counts = m, project = sid,
                                     min.features = 200L, min.cells = 0L)
    rm(m)
    so$sample_id <- sid
    so$donor <- manifest$donor[i]
    # Historically identical CELL MEMBERSHIP: nFeature>=200, nCount>=500,
    # NO mitochondrial-percentage exclusion (2620 cells at mt>=25 retained).
    so <- subset(so, subset = nFeature_RNA >= 200 & nCount_RNA >= 500)
    qc[[i]] <- data.frame(sample_id = sid, cells = ncol(so), genes = nrow(so))
    if (ncol(so) != manifest$expected_cells[i]) {
      write.csv(do.call(rbind, qc[seq_len(i)]),
                file.path(out_dir, "02_PER_SAMPLE_QC_FAILED.csv"), row.names = FALSE)
      stop("QC GATE FAIL for ", sid, ": observed ", ncol(so),
           "; expected ", manifest$expected_cells[i], call. = FALSE)
    }
    if (is.null(so_all)) so_all <- so else so_all <- merge(so_all, y = so,
                                                              merge.data = FALSE)
    rm(so); invisible(gc())
  }
    qc <- do.call(rbind, qc)
    write.csv(qc, file.path(out_dir, "02_PER_SAMPLE_QC_PASS.csv"), row.names = FALSE)
  }
  assert(ncol(so_all) == 62640L, "QC GATE: merged cell count is not 62640")
  assert(nrow(so_all) == 33539L, "QC GATE: merged gene count is not 33539")
  assert(!anyDuplicated(colnames(so_all)), "QC GATE: duplicated cell barcodes")
  # Historical identifiers can have repeated _1 suffixes from sequential merges.
  # Fail early on source-proven barcode examples before expensive PCA.
  expected_examples <- c("AAACCTGAGATGTCGG-1_1_1_1_1",
                         "AAACCTGCATGGTCAT-1_1_1_1_1")
  assert(all(expected_examples %in% colnames(so_all)),
         "BARCODE GATE FAIL: merged barcode naming differs from the Figure 7 source")
  message("[PASS] 62640 cells / 33539 genes / exact per-sample counts / barcode examples")

  ## 4 | MAINTAIN NINE RNA LAYERS FOR VST -----------------------------------
  DefaultAssay(so_all) <- "RNA"
  get_count_layers <- function(obj) grep("^counts(\\.|$)", Layers(obj[["RNA"]]), value = TRUE)
  if (length(get_count_layers(so_all)) > 1L) {
    # Join COUNTS ONLY, immediately split COUNTS by original sample.
    # Do NOT JoinLayers on the normalized data: doing so previously gave
    # 2329/3000 historical HVGs instead of 3000/3000.
    so_all <- JoinLayers(so_all, assay = "RNA", layers = "counts")
  }
  sample_factor <- factor(so_all$sample_id, levels = manifest$sample_id)
  assert(!anyNA(sample_factor), "LAYER GATE FAIL: unknown sample_id")
  so_all[["RNA"]] <- split(so_all[["RNA"]], f = sample_factor, layers = "counts")
  assert(length(get_count_layers(so_all)) == 9L,
         paste0("LAYER GATE FAIL: expected nine count layers, found ",
                length(get_count_layers(so_all))))
  message("[PASS] Nine unjoined sample-specific RNA count layers")

  ## 5 | PROVEN ANALYSIS PARAMETERS -----------------------------------------
  so_all <- NormalizeData(so_all, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = FALSE)
  so_all <- FindVariableFeatures(so_all, selection.method = "vst",
                                 nfeatures = 3000, loess.span = 0.3,
                                 verbose = FALSE)
  hvgs <- VariableFeatures(so_all)
  assert(length(hvgs) == 3000L, "HVG GATE FAIL: not exactly 3000 variable genes")
  assert(all(c("SCGB1A1", "SCGB3A1", "SERPINB3", "MT1G") == head(hvgs, 4)),
         "HVG GATE FAIL: first four original HVGs differ; stop before PCA")
  writeLines(hvgs, file.path(out_dir, "03_GEO_HVG3000.txt"))
  so_all <- ScaleData(so_all, features = hvgs, scale.max = 10, verbose = FALSE)
  so_all <- RunPCA(so_all, features = hvgs, npcs = 40, seed.use = 42,
                   verbose = FALSE)
  emb_pca <- Embeddings(so_all, reduction = "pca")
  assert(nrow(emb_pca) == 62640L && ncol(emb_pca) >= 40L,
         "PCA GATE FAIL: dimension mismatch")
  saveRDS(emb_pca[, 1:40, drop = FALSE],
          file.path(out_dir, "04_GEO_PCA40.rds"))
  message("[PASS] PCA40 created; note PCA numerical identity requires external audit")

  so_all <- FindNeighbors(so_all, reduction = "pca", dims = 1:30,
                          k.param = 20, nn.method = "annoy", n.trees = 50,
                          annoy.metric = "euclidean", verbose = FALSE)
  so_all <- FindClusters(so_all, resolution = 0.4, algorithm = 1,
                         random.seed = 0, verbose = FALSE)
  cluster_counts <- table(as.character(so_all$seurat_clusters))
  expected_clusters <- c(`0`=11048, `1`=9261, `2`=8349, `3`=6767,
    `4`=3778, `5`=3317, `6`=2325, `7`=2277, `8`=2261, `9`=2043,
    `10`=1797, `11`=1757, `12`=1254, `13`=1237, `14`=1223,
    `15`=992, `16`=851, `17`=713, `18`=458, `19`=368,
    `20`=168, `21`=143, `22`=107, `23`=105, `24`=41)
  assert(identical(as.integer(cluster_counts[names(expected_clusters)]),
                   as.integer(expected_clusters)),
         "CLUSTER GATE FAIL: exact 25-cluster sizes differ")
  write.csv(data.frame(cluster = names(cluster_counts), count = as.integer(cluster_counts)),
            file.path(out_dir, "05_GEO_CLUSTER_COUNTS_PASS.csv"), row.names = FALSE)
  write.csv(data.frame(cell_barcode = colnames(so_all),
                       cluster_id = as.character(so_all$seurat_clusters)),
            file.path(out_dir, "05_GEO_ALL_62640_CELL_CLUSTER_MAP.csv"), row.names = FALSE)
  message("[PASS] Exact 25 cluster sizes; all 62640 cell-level cluster IDs exported")

  ## 6 | GENES: EACH DATA LAYER -> CELL-LEVEL VALUES -------------------------
  # Reading data layers directly avoids silent partial FetchData behavior when
  # Seurat v5's RNA assay has nine data.* layers.
  expr <- matrix(NA_real_, nrow = ncol(so_all), ncol = 2L,
                 dimnames = list(colnames(so_all), c("CCL3", "CCL4")))
  data_layers <- grep("^data(\\.|$)", Layers(so_all[["RNA"]]), value = TRUE)
  assert(length(data_layers) == 9L, "EXPRESSION GATE FAIL: not nine data layers")
  for (lay in data_layers) {
    m <- LayerData(so_all[["RNA"]], layer = lay)
    assert(all(c("CCL3", "CCL4") %in% rownames(m)),
           paste("EXPRESSION GATE FAIL: missing CCL3/CCL4 in", lay))
    idx <- match(colnames(m), rownames(expr))
    assert(!anyNA(idx), paste("EXPRESSION GATE FAIL: unrecognized cell in", lay))
    assert(all(is.na(expr[idx, 1])), "EXPRESSION GATE FAIL: duplicate assignment")
    expr[idx, "CCL3"] <- as.numeric(m["CCL3", , drop = TRUE])
    expr[idx, "CCL4"] <- as.numeric(m["CCL4", , drop = TRUE])
    rm(m)
  }
  assert(all(is.finite(expr)), "EXPRESSION GATE FAIL: incomplete expression")

  cell_map <- c(`0`="Alveolar macrophage", `3`="T cell", `4`="Neutrophil",
                `6`="B cell")
  cl <- as.character(so_all$seurat_clusters)
  keep <- cl %in% names(cell_map)
  assert(sum(keep) == 23918L, "FIGURE GATE FAIL: expected 23918 selected cells")
  fig_cells <- colnames(so_all)[keep]
  fig_expr <- expr[fig_cells, , drop = FALSE]
  fig <- data.frame(cell_barcode = fig_cells,
                    celltype = unname(cell_map[cl[keep]]),
                    CCL3_expr = fig_expr[, 1],
                    CCL4_expr = fig_expr[, 2],
                    stringsAsFactors = FALSE)
  # Original documented method: 5th and 95th percentiles on THE SELECTED 23918.
  scale_q <- function(x) {
    q <- as.numeric(quantile(x, c(0.05, 0.95), names = FALSE, type = 7))
    assert(is.finite(q[1]) && is.finite(q[2]) && q[2] > q[1],
           "EXPRESSION GATE FAIL: invalid q5/q95 thresholds")
    list(value = pmin(pmax((x - q[1]) / (q[2] - q[1]), 0), 1), q = q)
  }
  s3 <- scale_q(fig$CCL3_expr)
  s4 <- scale_q(fig$CCL4_expr)
  fig$CCL3_scaled_0to1 <- s3$value
  fig$CCL4_scaled_0to1 <- s4$value
  fig$CCL3_CCL4_coexpression_score <- s3$value * s4$value
  fig$CCL3_CCL4_sum <- fig$CCL3_expr + fig$CCL4_expr
  thresholds <- data.frame(gene = c("CCL3", "CCL4"),
                           q05 = c(s3$q[1], s4$q[1]),
                           q95 = c(s3$q[2], s4$q[2]),
                           manuscript_q95 = c(4.684894666805451,
                                              5.157831318585424))
  thresholds$delta <- abs(thresholds$q95 - thresholds$manuscript_q95)
  write.csv(thresholds, file.path(out_dir, "06_CCL3_CCL4_QUANTILE_AUDIT.csv"),
            row.names = FALSE)
  assert(all(thresholds$delta < 1e-8),
         "EXPRESSION GATE FAIL: q95 differs from historical Figure 7 source")
  message("[PASS] 23918 selected cells; CCL3/CCL4 quantiles and co-expression")


  # Original GEO-only v1.3 UMAP, after expression and selection gates.
  suppressPackageStartupMessages(library(RSpectra))
  so_all <- RunUMAP(so_all, reduction='pca', dims=1:30,
                    n.neighbors=30L, n.components=2L, metric='cosine',
                    min.dist=0.3, seed.use=42, uwot.sgd=FALSE,
                    angular.rp.forest=FALSE, verbose=FALSE)
  um <- Embeddings(so_all, reduction='umap')[fig_cells,1:2,drop=FALSE]
  idx<-match(fig$cell_barcode,colnames(so_all))
  stopifnot(!anyNA(idx),length(idx)==23918L,identical(rownames(um),fig$cell_barcode))
  sample_name<-as.character(so_all$sample_id[idx])
  ali_plot <- data.frame(barcode=paste0('ALI::',fig$cell_barcode),sample=sample_name,
    celltype=fig$celltype, UMAP_1=um[,1], UMAP_2=um[,2],
    CCL3_expr=fig$CCL3_expr, CCL4_expr=fig$CCL4_expr,
    CCL3_scaled_0to1=fig$CCL3_scaled_0to1, CCL4_scaled_0to1=fig$CCL4_scaled_0to1,
    CCL3_CCL4_coexpression_score=fig$CCL3_CCL4_coexpression_score,stringsAsFactors=FALSE)
  saveRDS(ali_plot,file.path(work,'ali_plot_data.rds'))
  result<-data.frame(barcode=paste0('ALI::',fig$cell_barcode),sample=sample_name,CellType=fig$celltype,CCL3_expr=fig$CCL3_expr,CCL4_expr=fig$CCL4_expr,CCL3_CCL4_coexpr=fig$CCL3_CCL4_coexpression_score,stringsAsFactors=FALSE)
  }
  result
}

run_supp_data6 <- function(geo_dir,output_dir) {
  options(stringsAsFactors=FALSE,timeout=max(1800L,getOption('timeout')))
  req<-c('Seurat','SeuratObject','Matrix','hdf5r','ggplot2','patchwork','uwot','RSpectra','RcppAnnoy','openxlsx','zip')
  absent<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent))stop('Missing packages: ',paste(absent,collapse=', '))
  want<-c(R='4.5.2',Seurat='5.4.0',SeuratObject='5.2.0',uwot='0.2.4',RSpectra='0.16.2',RcppAnnoy='0.0.22')
  got<-vapply(names(want),function(x)if(x=='R')as.character(getRversion()) else as.character(utils::packageVersion(x)),character(1))
  if(any(got!=want))stop('Version mismatch: ',paste(sprintf('%s %s != %s',names(got),got,want),collapse=', '))
  geo_dir<-normalizePath(geo_dir,winslash='/',mustWork=TRUE)
  if(!dir.exists(output_dir))dir.create(output_dir,recursive=TRUE)
  output_dir<-normalizePath(output_dir,winslash='/',mustWork=TRUE)
  if(tolower(output_dir)==tolower(geo_dir)||startsWith(tolower(output_dir),paste0(tolower(sub('/+$','',geo_dir)),'/')))stop('Output directory must be outside raw input directory.')
  destination<-file.path(output_dir,'Supplementary_Data_6.xlsx');work<-file.path(output_dir,'internal_qc_fig7cd')
  if(file.exists(destination)||dir.exists(work))stop('Output already exists; choose a new directory.')
  dir.create(work)
  # Preserve the audited workflow's ALI-before-Control package/load order.
  ali<-run_ali_supp(geo_dir,work);invisible(gc())
  suppressPackageStartupMessages(library(RSpectra))
  control<-run_control_supp(geo_dir,work);invisible(gc())
  stopifnot(nrow(control)==20327L,nrow(ali)==23918L,!anyDuplicated(control$barcode),!anyDuplicated(ali$barcode),identical(names(control),names(ali)))
  wb<-openxlsx::createWorkbook();openxlsx::addWorksheet(wb,'per_cell_H.C');openxlsx::addWorksheet(wb,'per_cell_ALI');openxlsx::writeData(wb,'per_cell_H.C',control);openxlsx::writeData(wb,'per_cell_ALI',ali)
  tmp<-tempfile(fileext='.xlsx');on.exit(unlink(tmp),add=TRUE);openxlsx::saveWorkbook(wb,tmp,overwrite=FALSE);repair_xlsx(tmp,destination)
  write.csv(data.frame(package=names(got),version=unname(got)),file.path(work,'package_versions.csv'),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(work,'session_info.txt'))
  message('Supplementary Data 6: ',destination)
  invisible(destination)
}
if(sys.nframe()==0L){args<-commandArgs(trailingOnly=TRUE);if(length(args)!=2L)stop('Usage: Rscript 01_GEO_to_SuppData6.R GEO_H5_DIR OUTPUT_DIR');run_supp_data6(args[1],args[2])}

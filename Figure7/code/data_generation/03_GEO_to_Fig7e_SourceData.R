#!/usr/bin/env Rscript
# Figure 7e: twelve original GSE145926 10x H5 files -> per-cell data,
# donor medians and Welch/BH statistics. No Figure 7c,d selection is reused.
# Usage: Rscript 03_GEO_to_Fig7e_SourceData.R GEO_H5_DIR OUTPUT_DIR
# OUTPUT_DIR must not already contain the result workbook; H5 is read only.

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
if(length(sheet_names)!=3L)stop('Unexpected number of worksheets: ',length(sheet_names))
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

run_fig7e <- function(geo_dir,output_dir) {
  options(Seurat.object.assay.version='v3',stringsAsFactors=FALSE)
  required <- c('Seurat','SeuratObject','Matrix','hdf5r','dplyr','tidyr','openxlsx','zip')
  absent <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent)) stop('Missing packages: ',paste(absent,collapse=', '))
  required_versions <- c(R='4.5.2',Seurat='5.4.0',SeuratObject='5.2.0',Matrix='1.7.4',hdf5r='1.3.16',dplyr='1.1.4',tidyr='1.3.1',openxlsx='4.2.8.1',zip='2.3.3')
  actual <- vapply(names(required_versions),function(p) if(p=='R')as.character(getRversion()) else as.character(utils::packageVersion(p)),character(1))
  if(any(actual!=required_versions)) stop('Package versions differ from the verified environment: ',paste(sprintf('%s expected=%s actual=%s',names(actual),required_versions,actual),collapse='; '))
  suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(Matrix);library(dplyr);library(tidyr)})
  geo_dir <- normalizePath(geo_dir,winslash='/',mustWork=TRUE)
  if(!dir.exists(output_dir))dir.create(output_dir,recursive=TRUE)
  output_dir <- normalizePath(output_dir,winslash='/',mustWork=TRUE)
  if(tolower(output_dir)==tolower(geo_dir) || startsWith(tolower(output_dir),paste0(tolower(sub('/+$','',geo_dir)),'/'))) stop('Output directory must be outside raw data directory.')
  output_file <- file.path(output_dir,'Figure7e_Source_Data.xlsx')
  audit_dir <- file.path(output_dir,'internal_qc')
  if(file.exists(output_file)||dir.exists(audit_dir))stop('Output already exists; choose a new output directory.')
  dir.create(audit_dir)
QC_MIN_FEATURES <- 200L      # CreateSeuratObject min.features
QC_MIN_CELLS <- 3L          # CreateSeuratObject min.cells per gene
QC_MAX_MT_FRAC <- 0.15      # percent.mt < 15% (strict inequality)
AM_Q <- 0.70                # AM_delta top 30% cutoff
AM_MIN_KEEP <- 30L
AM_POS_EXPR_THRESHOLD <- 0.25
EPITHELIAL_EXPR_THRESHOLD <- 0.05
AM_MEAN_DIFF_THRESHOLD <- 0.20
AM_FALLBACK_Q <- c(0.65, 0.60, 0.55, 0.50)
MODULE_SCORE_SEED <- 1L     # Seurat AddModuleScore default in source-era code
PLOT_JITTER_SEED <- 7L
NORMALIZE_SCALE_FACTOR <- 10000

AM_POS <- c(
  "MARCO", "MRC1", "PPARG", "APOE", "APOC1",
  "CD68", "CSF1R", "MERTK", "LPL", "SIGLEC1"
)
MONO_POS <- c("FCN1", "LST1", "S100A8", "S100A9", "VCAN", "LYZ", "CCR2")
NEG_EPI <- c("EPCAM", "KRT19", "KRT8", "KRT18")

manifest <- data.frame(
  gsm = c("GSM4475048","GSM4475049","GSM4475050",
          "GSM4339769","GSM4339770","GSM4339771","GSM4339772",
          "GSM4339773","GSM4339774","GSM4475051","GSM4475052",
          "GSM4475053"),
  donor = c("C51","C52","C100",
            "C141","C142","C143","C144","C145","C146",
            "C148","C149","C152"),
  group = c(rep("Control",3L), rep("ALI",9L)),
  stringsAsFactors = FALSE
)
manifest$sample <- paste(manifest$gsm, manifest$donor, sep="_")
manifest$filename <- paste0(
  manifest$sample,"_filtered_feature_bc_matrix.h5")
needed_files <- manifest$filename
stopifnot(nrow(manifest)==12L,
          sum(manifest$group=="Control")==3L,
          sum(manifest$group=="ALI")==9L,
          !anyDuplicated(manifest$gsm),
          length(needed_files)==12L,
          !anyDuplicated(needed_files),
          !"GSM3660650" %in% manifest$gsm)


geo_paths <- file.path(geo_dir,manifest$filename)
names(geo_paths)<-manifest$filename
if(any(!file.exists(geo_paths)))stop('Missing GEO H5: ',paste(manifest$filename[!file.exists(geo_paths)],collapse=', '))
input_manifest<-data.frame(gsm=manifest$gsm,donor=manifest$donor,group=manifest$group,filename=manifest$filename,md5=unname(tools::md5sum(geo_paths)))
## ===== 2. DATA ACCESS HELPERS =====

get_layer <- function(so, layer = c("data", "counts")) {
  layer <- match.arg(layer)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    SeuratObject::LayerData(so, assay = "RNA", layer = layer)
  } else {
    Seurat::GetAssayData(so, assay = "RNA", slot = layer)
  }
}

get_gene <- function(so, gene, layer = "data", required = FALSE) {
  rn <- rownames(so)
  hit <- rn[match(toupper(gene), toupper(rn))]
  if (is.na(hit)) {
    if (required) stop("Required gene absent: ", gene, call. = FALSE)
    return(rep(0, ncol(so))) # same fallback as original for selection markers
  }
  as.numeric(get_layer(so, layer)[hit, , drop = TRUE])
}

present_markers <- function(rn, genes) {
  hits <- rn[match(toupper(genes), toupper(rn))]
  unique(hits[!is.na(hits)])
}

read_geo_matrix <- function(f) {
  m <- Seurat::Read10X_h5(f)
  if (is.list(m)) {
    if ("Gene Expression" %in% names(m)) {
      m <- m[["Gene Expression"]]
    } else {
      stop("H5 contains multiple matrices and no 'Gene Expression': ", f)
    }
  }
  if (!inherits(m, "Matrix") && !is.matrix(m)) {
    stop("Invalid count matrix: ", f)
  }
  m
}

## ===== 3. PER-SAMPLE QC, AM IDENTIFICATION, EXPRESSION =====

process_one <- function(filepath, sample, donor, gsm, group) {
  message("\n[", group, "] ", donor, " / ", gsm)

  counts <- read_geo_matrix(filepath)
  if (!inherits(counts,"Matrix") && !is.matrix(counts))
    stop("Invalid sparse count matrix: ",sample)
  so <- Seurat::CreateSeuratObject(
    counts = counts,
    project = group,
    min.cells = QC_MIN_CELLS,
    min.features = QC_MIN_FEATURES
  )
  rm(counts)
  gc(verbose = FALSE)
  if (ncol(so) == 0L) stop("No cells after initial QC: ", sample)

  DefaultAssay(so) <- "RNA"
  so$sample <- sample
  so$donor <- donor
  so$gsm <- gsm
  so$group <- group

  # Legacy percent.mt calculation: fraction of reads in ^MT- genes.
  mt_genes <- rownames(so)[grepl("^MT-", rownames(so))]
  totals <- Matrix::colSums(get_layer(so, "counts"))
  mt_fraction <- if (length(mt_genes)) {
    Matrix::colSums(get_layer(so, "counts")[mt_genes, , drop = FALSE]) /
      totals
  } else {
    rep(0, ncol(so))
  }
  so[["percent.mt"]] <- mt_fraction * 100
  qc_before_mt <- ncol(so)
  so <- subset(so, subset = percent.mt < (QC_MAX_MT_FRAC * 100))
  if (ncol(so) == 0L) stop("No cells remain after MT filtering: ", sample)
  qc_after_mt <- ncol(so)

  so <- Seurat::NormalizeData(
    so, normalization.method = "LogNormalize",
    scale.factor = NORMALIZE_SCALE_FACTOR, verbose = FALSE
  )

  rn <- rownames(so)
  am_mks <- present_markers(rn, AM_POS)
  mono_mks <- present_markers(rn, MONO_POS)
  if (!length(am_mks) || !length(mono_mks)) {
    stop("Insufficient AM/monocyte markers in sample: ", sample)
  }

  # Legacy AddModuleScore parameters/default seed made explicit.
  so <- Seurat::AddModuleScore(
    so, features = list(am_mks), name = "AM_Score",
    nbin = 24, ctrl = 100, seed = MODULE_SCORE_SEED, verbose = FALSE
  )
  so <- Seurat::AddModuleScore(
    so, features = list(mono_mks), name = "Mono_Score",
    nbin = 24, ctrl = 100, seed = MODULE_SCORE_SEED, verbose = FALSE
  )
  so$AM_delta <- so$AM_Score1 - so$Mono_Score1

  threshold <- stats::quantile(so$AM_delta, probs = AM_Q, na.rm = TRUE)
  idxA <- SeuratObject::WhichCells(so, expression = AM_delta >= threshold)

  assay_data <- get_layer(so, "data")
  am_use <- intersect(rownames(assay_data), am_mks)
  mono_use <- intersect(rownames(assay_data), mono_mks)
  idxB <- if (length(am_use) && length(mono_use)) {
    am_mean <- Matrix::colMeans(assay_data[am_use, , drop = FALSE])
    mono_mean <- Matrix::colMeans(assay_data[mono_use, , drop = FALSE])
    colnames(so)[(am_mean - mono_mean) >= AM_MEAN_DIFF_THRESHOLD]
  } else {
    character(0L)
  }

  pos_any <- (get_gene(so, "MARCO") >= AM_POS_EXPR_THRESHOLD) |
    (get_gene(so, "MRC1") >= AM_POS_EXPR_THRESHOLD) |
    (get_gene(so, "CSF1R") >= AM_POS_EXPR_THRESHOLD)
  epi_positive <- sapply(
    NEG_EPI,
    function(g) get_gene(so, g) >= EPITHELIAL_EXPR_THRESHOLD
  )
  no_epi <- rowSums(epi_positive) == 0
  idxC <- colnames(so)[pos_any & no_epi]

  nonempty <- Filter(length, list(idxA, idxB, idxC))
  idx <- if (length(nonempty)) Reduce(intersect, nonempty) else character()
  selection <- "intersection"
  if (length(idx) < AM_MIN_KEEP) {
    idx <- unique(c(idxA, idxB, idxC))
    selection <- "union_fallback"
  }
  if (length(idx) < AM_MIN_KEEP) {
    for (q in AM_FALLBACK_Q) {
      thrA <- stats::quantile(so$AM_delta, probs = q, na.rm = TRUE)
      idxA <- SeuratObject::WhichCells(so, expression = AM_delta >= thrA)
      idx <- unique(c(idxA, idxB, idxC))
      selection <- paste0("union_fallback_q", q)
      if (length(idx) >= AM_MIN_KEEP) break
    }
  }
  message("  AM retained: ", length(idx), " [", selection, "]")
  if (!length(idx)) stop("No AM cells selected: ", sample)

  am <- subset(so, cells = idx)
  DefaultAssay(am) <- "RNA"

  # Legacy pseudobulk CPM from counts of selected AMs.
  ccl3_counts <- sum(get_gene(am, "CCL3", "counts", required=TRUE))
  ccl4_counts <- sum(get_gene(am, "CCL4", "counts", required=TRUE))
  all_counts <- sum(Matrix::colSums(get_layer(am, "counts")))
  cpm3 <- if (all_counts > 0) ccl3_counts / all_counts * 1e6 else NA_real_
  cpm4 <- if (all_counts > 0) ccl4_counts / all_counts * 1e6 else NA_real_

  # Legacy cell-level normalized expression and arithmetic SUM.
  # Do not confuse the sum with a co-expression *product*.
  met <- data.frame(
    cell = paste(sample, colnames(am), sep = "::"),
    sample = sample,
    donor = donor,
    gsm = gsm,
    group = group,
    CCL3_expr = get_gene(am, "CCL3", "data", required=TRUE),
    CCL4_expr = get_gene(am, "CCL4", "data", required=TRUE),
    stringsAsFactors = FALSE
  )
  met$CCL_sum <- met$CCL3_expr + met$CCL4_expr

  summ <- data.frame(
    sample = sample,
    donor = donor,
    gsm = gsm,
    group = group,
    qc_after_initial = qc_before_mt,
    qc_after_mt = qc_after_mt,
    am_n = nrow(met),
    am_selection = selection,
    CPM_CCL3 = cpm3,
    CPM_CCL4 = cpm4,
    mean_CCL3 = mean(met$CCL3_expr),
    mean_CCL4 = mean(met$CCL4_expr),
    mean_CCLsum = mean(met$CCL_sum),
    median_CCL3 = stats::median(met$CCL3_expr),
    median_CCL4 = stats::median(met$CCL4_expr),
    median_CCLsum = stats::median(met$CCL_sum),
    stringsAsFactors = FALSE
  )

  rm(so, am, assay_data)
  gc(verbose = FALSE)
  list(cells = met, sample = summ)
}

processed <- lapply(seq_len(nrow(manifest)), function(i) {
  process_one(
    filepath = geo_paths[[manifest$filename[i]]],
    sample = tools::file_path_sans_ext(manifest$filename[i]),
    donor = manifest$donor[i],
    gsm = manifest$gsm[i],
    group = manifest$group[i]
  )
})
am_cells <- dplyr::bind_rows(lapply(processed, `[[`, "cells"))
am_samples <- dplyr::bind_rows(lapply(processed, `[[`, "sample"))
am_cells$group <- factor(am_cells$group, levels = c("Control", "ALI"))
am_samples$group <- factor(am_samples$group, levels = c("Control", "ALI"))

stopifnot(
  nrow(am_cells) > 0L,
  !anyDuplicated(am_cells$cell),
  n_distinct(am_samples$sample) == 12L,
  all(abs(am_cells$CCL_sum - am_cells$CCL3_expr -
            am_cells$CCL4_expr) < 1e-8)
)

p_to_star <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}
median_columns <- c(CCL3="median_CCL3", CCL4="median_CCL4",
                    "CCL3+CCL4"="median_CCLsum")
main_welch <- dplyr::bind_rows(lapply(names(median_columns), function(nm) {
  col <- median_columns[[nm]]
  x <- as.numeric(am_samples[[col]][am_samples$group=="Control"])
  y <- as.numeric(am_samples[[col]][am_samples$group=="ALI"])
  stopifnot(length(x)==3L, length(y)==9L,
            all(is.finite(x)), all(is.finite(y)))
  tt <- stats::t.test(x,y,var.equal=FALSE, alternative="two.sided")
  data.frame(feature=nm, n_control=length(x), n_ALI=length(y),
             control_mean_of_sample_medians=mean(x),
             ALI_mean_of_sample_medians=mean(y),
             mean_difference_ALI_minus_Control=mean(y)-mean(x),
             welch_t=unname(tt$statistic),
             welch_df=unname(tt$parameter),
             p_value=tt$p.value, stringsAsFactors=FALSE)
}))
main_welch$p_BH <- stats::p.adjust(main_welch$p_value,method="BH")
main_welch$raw_star <- vapply(main_welch$p_value,p_to_star,character(1L))
main_welch$BH_star <- vapply(main_welch$p_BH,p_to_star,character(1L))

# Release workbook contains figure data and donor-level statistics only.
release_sample <- am_samples[,c('sample','donor','gsm','group','am_n','median_CCL3','median_CCL4','median_CCLsum'),drop=FALSE]
wb<-openxlsx::createWorkbook()
for(s in c('Fig7e_Cell','Fig7e_Sample_Medians','Fig7e_Welch_BH'))openxlsx::addWorksheet(wb,s)
openxlsx::writeData(wb,'Fig7e_Cell',am_cells)
openxlsx::writeData(wb,'Fig7e_Sample_Medians',release_sample)
openxlsx::writeData(wb,'Fig7e_Welch_BH',main_welch)
tmp<-tempfile(fileext='.xlsx');on.exit(unlink(tmp),add=TRUE)
openxlsx::saveWorkbook(wb,tmp,overwrite=FALSE)
repair_xlsx(tmp,output_file)
write.csv(input_manifest,file.path(audit_dir,'input_manifest_md5.csv'),row.names=FALSE)
write.csv(am_samples,file.path(audit_dir,'per_sample_QC.csv'),row.names=FALSE)
write.csv(data.frame(package=names(actual),version=unname(actual)),file.path(audit_dir,'package_versions.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(audit_dir,'session_info.txt'))
message('Figure 7e output: ',output_file)
invisible(output_file)
}

if(sys.nframe()==0L){
 args<-commandArgs(trailingOnly=TRUE)
 if(length(args)!=2L)stop('Usage: Rscript 03_GEO_to_Fig7e_SourceData.R GEO_H5_DIR OUTPUT_DIR')
 run_fig7e(args[1],args[2])
}

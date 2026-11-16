#!/usr/bin/env Rscript
# GEO-derived Figure 7c,d Source Data, reusing the same-run analysis in 01.
# Usage: Rscript 02_GEO_to_Fig7cd_SourceData.R GEO_H5_DIR OUTPUT_DIR
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Usage: Rscript 02_GEO_to_Fig7cd_SourceData.R GEO_H5_DIR OUTPUT_DIR')
scriptarg <- grep('^--file=',commandArgs(),value=TRUE)
if(length(scriptarg)!=1L)stop('Run this script with Rscript.')
scriptdir <- dirname(normalizePath(sub('^--file=','',scriptarg),winslash='/',mustWork=TRUE))
source(file.path(scriptdir,'01_GEO_to_SuppData6.R'),local=TRUE)
geo <- normalizePath(args[1],winslash='/',mustWork=TRUE)
out <- args[2]
if(!dir.exists(out))dir.create(out,recursive=TRUE)
out <- normalizePath(out,winslash='/',mustWork=TRUE)
dest <- file.path(out,'Figure7cd_Source_Data.xlsx')
if(file.exists(dest))stop('Refusing to overwrite: ',dest)
work <- file.path(out,'internal_qc_fig7cd')
if(!dir.exists(work)) run_supp_data6(geo,out)
# Confirm that reused intermediate tables belong to these original H5 inputs.
a<-read.csv(file.path(work,'ali','00_GEO_H5_MANIFEST_AND_MD5.csv'),stringsAsFactors=FALSE)
c<-read.csv(file.path(work,'control','00_GEO_INPUT_MANIFEST.csv'),stringsAsFactors=FALSE)
paths<-file.path(geo,c(a$filename,c$h5)); hashes<-c(a$observed_md5,c$observed_md5)
if(length(paths)!=12L || any(!file.exists(paths)) || any(tolower(unname(tools::md5sum(paths)))!=tolower(hashes)))stop('GEO files differ from the inputs used in 01; use a fresh output directory.')
# No workbook or historical coordinates are read to create these sheets.
files <- file.path(work,c('control_plot_data.rds','ali_plot_data.rds'))
if(!all(file.exists(files)))stop('Incomplete prior 01 run; use a fresh output directory.')
control <- readRDS(files[1]); ali <- readRDS(files[2])
expected <- c('barcode','sample','celltype','UMAP_1','UMAP_2','CCL3_expr','CCL4_expr','CCL3_scaled_0to1','CCL4_scaled_0to1','CCL3_CCL4_coexpression_score')
stopifnot(identical(names(control),expected),identical(names(ali),expected),nrow(control)==20327L,nrow(ali)==23918L,!anyDuplicated(control$barcode),!anyDuplicated(ali$barcode),all(vapply(control[4:10],function(x)all(is.finite(x)),logical(1))),all(vapply(ali[4:10],function(x)all(is.finite(x)),logical(1))))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb,'Control');openxlsx::addWorksheet(wb,'ALI')
openxlsx::writeData(wb,'Control',control);openxlsx::writeData(wb,'ALI',ali)
tmp <- tempfile(fileext='.xlsx');on.exit(unlink(tmp),add=TRUE)
openxlsx::saveWorkbook(wb,tmp,overwrite=FALSE);repair_xlsx(tmp,dest)
message('Figure 7c,d Source Data: ',dest)

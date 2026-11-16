#!/usr/bin/env Rscript
# Run the exact release scripts in fresh R sessions, from raw H5.
# Usage: Rscript 00_RUN_ALL.R GEO_H5_DIR NEW_OUTPUT_DIR
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Usage: Rscript 00_RUN_ALL.R GEO_H5_DIR NEW_OUTPUT_DIR')
scriptarg<-grep('^--file=',commandArgs(),value=TRUE)
if(length(scriptarg)!=1L)stop('Run this script with Rscript.')
scriptdir<-dirname(normalizePath(sub('^--file=','',scriptarg),winslash='/',mustWork=TRUE))
geo<-normalizePath(args[1],winslash='/',mustWork=TRUE)
out<-args[2]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Use a new empty output directory.')
if(!dir.exists(out))dir.create(out,recursive=TRUE)
out<-normalizePath(out,winslash='/',mustWork=TRUE)
if(tolower(out)==tolower(geo)||startsWith(tolower(out),paste0(tolower(sub('/+$','',geo)),'/')))stop('Output must be outside GEO input.')
rscript<-file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
for(s in c('01_GEO_to_SuppData6.R','02_GEO_to_Fig7cd_SourceData.R','03_GEO_to_Fig7e_SourceData.R')) {
 message('\n========== ',s,' ==========')
 status<-system2(rscript,args=c(shQuote(file.path(scriptdir,s)),shQuote(geo),shQuote(out)))
 if(!identical(status,0L))stop(s,' failed with exit status ',status,'; inspect output and rerun in a new directory.')
}
message('Completed all three workbooks in ',out)

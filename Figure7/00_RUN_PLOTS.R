#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=FALSE)
f <- grep('^--file=',args,value=TRUE)
if(length(f)!=1L) stop('Run with Rscript 00_RUN_PLOTS.R')
root <- dirname(normalizePath(sub('^--file=','',f),winslash='/',mustWork=TRUE))
cli <- commandArgs(trailingOnly=TRUE)
flags <- cli[startsWith(cli,'--')]
paths <- cli[!startsWith(cli,'--')]
allowed <- c('--include-control-source','--control-source-only')
if(length(paths)>1L || length(flags)>1L || any(!flags %in% allowed))
  stop('Usage: Rscript 00_RUN_PLOTS.R [NEW_OUTPUT_DIR] [--include-control-source | --control-source-only]')
out <- if(length(paths)) paths[1] else file.path(root,'figures')
if(dir.exists(out) && length(list.files(out, all.files=TRUE, no..=TRUE))) stop('Use a new empty output directory: ',out)
if(!dir.exists(out))dir.create(out,recursive=TRUE)
out <- normalizePath(out,winslash='/',mustWork=TRUE)
Sys.setenv(FIG7_OUTPUT_DIR=out, FIG7_CONTROL_COORD_MODE='manuscript',
  FIG7_CONTROL_CSV=file.path(root,'data/manuscript_plot_inputs/Figure7cd_Control_Manuscript.csv'),
  FIG7_ALI_CSV=file.path(root,'data/manuscript_plot_inputs/Figure7cd_ALI_Manuscript.csv'),
  FIG7E_INPUT_XLSX=file.path(root,'data/reference_workbooks/Figure7e_Source_Data.xlsx'),
  FIG7_INPUT_XLSX=file.path(root,'data/reference_workbooks/Figure_7f_7g_plotting_data_GSE58682.xlsx'))
rs <- file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
run_plot <- function(script) {
  path<-file.path(root,'code','plotting',script)
  message('Drawing: ', script, ' [Control mode: ', Sys.getenv('FIG7_CONTROL_COORD_MODE'), ']')
  rc<-system2(rs,shQuote(path))
  if(!identical(rc,0L))stop('Plot failed: ',script,' (exit ',rc,')')
}
if(!('--control-source-only' %in% flags))
  for(script in c('04_Plot_Fig7cd_Control.R','05_Plot_Fig7cd_ALI.R',
                  '06_Plot_Fig7e.R','07_Plot_Fig7fg.R')) run_plot(script)
if(any(flags %in% allowed)) {
  Sys.setenv(FIG7_CONTROL_COORD_MODE='source')
  run_plot('04_Plot_Fig7cd_Control.R')
}
message('Completed requested Figure 7 plots in: ',out)

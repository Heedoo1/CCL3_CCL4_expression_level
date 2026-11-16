#!/usr/bin/env Rscript
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
if(length(sheet_names)!=6L)stop('Unexpected number of worksheets: ',length(sheet_names))
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


# ============================================================
# FIGURE 7e - SOURCE DATA VIOLIN GRAPH + DONOR-LEVEL WELCH/BH
# Header detection corrected
# BH-adjusted significance annotations drawn on each facet
# ============================================================

# 1. Packages
pkgs <- c("readxl", "ggplot2", "openxlsx")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1),
                        quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
library(readxl)
library(ggplot2)
library(openxlsx)

# 2. Paths
input <- Sys.getenv('FIG7E_INPUT_XLSX', unset = file.path(getwd(),
  'data', 'reference_workbooks', 'Figure7e_Source_Data.xlsx'))
if (!file.exists(input)) stop('Missing Figure 7e data: ', input)
out <- file.path(Sys.getenv('FIG7_OUTPUT_DIR', unset = file.path(getwd(),'figures')), 'Fig7e')
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cat("INPUT:", normalizePath(input), "\n")
cat("SHEETS:", paste(excel_sheets(input), collapse = ", "), "\n")

# 3. Helpers
norm <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("[^a-z0-9]", "", tolower(trimws(x)))
}

find_col <- function(header, candidates) {
  ix <- which(header %in% candidates)
  if (length(ix)) ix[1] else NA_integer_
}

get_col <- function(dat, ix) {
  if (is.na(ix)) return(rep(NA_character_, nrow(dat)))
  as.character(dat[[ix]])
}

to_num <- function(x) suppressWarnings(as.numeric(x))

group_name <- function(x) {
  z <- norm(x)
  ifelse(
    z %in% c("control", "healthy", "healthycontrol", "hc", "normal"),
    "Control",
    ifelse(
      z %in% c("ali", "acutelunginjury"),
      "ALI",
      NA_character_
    )
  )
}

group_keys <- c(
  "v7egroup", "group", "condition", "samplegroup",
  "cellgroup", "status", "diseasegroup"
)

c3_keys <- c(
  "v7eccl3expr", "ccl3", "ccl3expr",
  "ccl3expression", "ccl3normalizedexpression",
  "normalizedccl3", "exprccl3"
)

c4_keys <- c(
  "v7eccl4expr", "ccl4", "ccl4expr",
  "ccl4expression", "ccl4normalizedexpression",
  "normalizedccl4", "exprccl4"
)

sum_keys <- c(
  "v7eccl3sum", "ccl3ccl4", "ccl34",
  "ccl3plusccl4", "cclsum",
  "ccl3ccl4expr", "ccl3ccl4expression"
)

sample_keys <- c(
  "v7esample", "sample", "sampleid",
  "gsm", "gsmid", "samplename"
)

# 4. Automatically identify source columns
cells <- NULL
support <- NULL
metadata <- NULL
detected_sheet <- NA_character_
detected_row <- NA_integer_

for (sheet in excel_sheets(input)) {
  
  grid <- suppressWarnings(
    read_excel(
      input,
      sheet = sheet,
      col_names = FALSE,
      col_types = "text",
      .name_repair = "minimal"
    )
  )
  
  if (nrow(grid) < 2) next
  
  max_row <- min(100L, nrow(grid) - 1L)
  
  for (r in seq_len(max_row)) {
    
    header <- norm(unlist(grid[r, , drop = FALSE],
                          use.names = FALSE))
    
    ig <- find_col(header, group_keys)
    i3 <- find_col(header, c3_keys)
    i4 <- find_col(header, c4_keys)
    is <- find_col(header, sum_keys)
    iid <- find_col(header, sample_keys)
    
    if (is.na(ig) || is.na(i3) || is.na(i4)) next
    
    dat <- as.data.frame(
      grid[seq.int(r + 1L, nrow(grid)), , drop = FALSE]
    )
    
    z <- data.frame(
      Group = group_name(get_col(dat, ig)),
      Sample = get_col(dat, iid),
      CCL3 = to_num(get_col(dat, i3)),
      CCL4 = to_num(get_col(dat, i4)),
      stringsAsFactors = FALSE
    )
    
    if (!is.na(is)) {
      z$CCL3_4 <- to_num(get_col(dat, is))
    } else {
      z$CCL3_4 <- z$CCL3 + z$CCL4
      message("CCL3+4 column absent: calculated CCL3 + CCL4.")
    }
    
    z <- z[
      !is.na(z$Group) &
        is.finite(z$CCL3) &
        is.finite(z$CCL4) &
        is.finite(z$CCL3_4),
      ,
      drop = FALSE
    ]
    
    if (nrow(z) < 10 ||
        length(unique(z$Group)) != 2) next
    
    # Exclude the extra historical Control sample
    if (any(grepl("GSM3660650", z$Sample))) {
      z <- z[!grepl("GSM3660650", z$Sample), ,
             drop = FALSE]
    }
    
    cells <- z
    detected_sheet <- sheet
    detected_row <- r
    
    # Optional supporting sample-level means
    sx <- find_col(header, c("vsamplesample"))
    sg <- find_col(header, c("vsamplegroup"))
    s3 <- find_col(header, c("vsamplemeanccl3"))
    s4 <- find_col(header, c("vsamplemeanccl4"))
    ss <- find_col(header, c("vsamplemeancclsum"))
    
    if (all(!is.na(c(sx, sg, s3, s4, ss)))) {
      support <- data.frame(
        Sample = get_col(dat, sx),
        Group = group_name(get_col(dat, sg)),
        CCL3 = to_num(get_col(dat, s3)),
        CCL4 = to_num(get_col(dat, s4)),
        CCL3_4 = to_num(get_col(dat, ss)),
        stringsAsFactors = FALSE
      )
      
      support <- support[
        !is.na(support$Sample) &
          nzchar(support$Sample) &
          !is.na(support$Group) &
          !grepl("GSM3660650", support$Sample),
        ,
        drop = FALSE
      ]
      
      support <- support[!duplicated(support$Sample), ,
                         drop = FALSE]
    }
    
    # Optional original manuscript metadata
    mf <- find_col(header, c("meta7efeature"))
    ms <- find_col(header, c("meta7edisplaysignif"))
    mp <- find_col(header, c("meta7ecelllevelp"))
    
    if (!is.na(mf) && !is.na(ms)) {
      metadata <- data.frame(
        Gene = get_col(dat, mf),
        Stars = get_col(dat, ms),
        Source_p = get_col(dat, mp),
        stringsAsFactors = FALSE
      )
      metadata <- metadata[
        metadata$Gene %in% c("CCL3", "CCL4", "CCL3+4"),
        ,
        drop = FALSE
      ]
    }
    
    break
  }
  
  if (!is.null(cells)) break
}

# 5. Diagnostics if structure is different
if (is.null(cells)) {
  
  cat("\nSOURCE DATA HEADER DIAGNOSTICS\n")
  
  for (sheet in excel_sheets(input)) {
    tmp <- suppressWarnings(
      read_excel(input, sheet = sheet,
                 col_names = FALSE, col_types = "text",
                 n_max = 15)
    )
    
    cat("\nSHEET:", sheet, "\n")
    
    for (r in seq_len(nrow(tmp))) {
      values <- as.character(
        unlist(tmp[r, , drop = FALSE], use.names = FALSE)
      )
      values <- values[!is.na(values) & nzchar(values)]
      
      hits <- values[grepl(
        "CCL3|CCL4|GROUP|CONTROL|ALI|EXPRESSION",
        values, ignore.case = TRUE
      )]
      
      if (length(hits)) {
        cat("ROW", r, ":",
            paste(head(hits, 15), collapse = " | "), "\n")
      }
    }
  }
  
  stop(
    "Automatic header detection failed. ",
    "See SOURCE DATA HEADER DIAGNOSTICS above."
  )
}

cat("\nDETECTED SHEET:", detected_sheet, "\n")
cat("HEADER ROW:", detected_row, "\n")
cat("\nCELL COUNTS:\n")
print(table(cells$Group))
cat("TOTAL:", nrow(cells), "\n")

if (sum(cells$Group == "Control") != 6333 ||
    sum(cells$Group == "ALI") != 6604) {
  warning(
    "Counts differ from the selected 12-sample reference ",
    "(Control 6333 / ALI 6604). Check source data."
  )
}

# 6. Long-format figure data
make_long <- function(dat) {
  do.call(rbind, lapply(
    c("CCL3", "CCL4", "CCL3_4"),
    function(g) {
      data.frame(
        Group = dat$Group,
        Sample = dat$Sample,
        Gene = if (g == "CCL3_4") "CCL3+4" else g,
        Expression = dat[[g]],
        stringsAsFactors = FALSE
      )
    }
  ))
}

plot_data <- make_long(cells)
plot_data$Group <- factor(
  plot_data$Group, levels = c("Control", "ALI")
)
plot_data$Gene <- factor(
  plot_data$Gene,
  levels = c("CCL3", "CCL4", "CCL3+4")
)

# 7. Figure - cell-level violin plots and BH-adjusted significance marks
fill_col <- c(Control = "#F3B0B0", ALI = "#A5D1F3")
point_col <- c(Control = "#ED6868", ALI = "#2E83C7")

y_top <- max(13.4, max(plot_data$Expression) + 0.3)

p <- ggplot(plot_data,
            aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(
    width = 0.82, scale = "width", trim = TRUE,
    color = "#696969", linewidth = 0.32, alpha = 0.95
  ) +
  geom_boxplot(
    width = 0.085, outlier.shape = NA,
    fill = "#D9D9D9", color = "#666666",
    alpha = 0.65, linewidth = 0.28
  ) +
  geom_point(
    aes(color = Group),
    position = position_jitter(
      width = 0.065, height = 0, seed = 7
    ),
    size = 0.32, alpha = 0.62,
    show.legend = FALSE
  ) +
  facet_wrap(~Gene, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = fill_col) +
  scale_color_manual(values = point_col) +
  scale_y_continuous(
    breaks = c(0, 4, 8, 12),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(0, y_top)) +
  scale_x_discrete(expand = expansion(add = 0.36)) +
  labs(
    x = NULL,
    y = "Normalized expression in AM cells",
    tag = "e"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.background = element_rect(
      fill = "#FBFBFB", color = NA
    ),
    panel.border = element_rect(
      fill = NA, color = "#777777", linewidth = 0.55
    ),
    panel.grid.major = element_line(
      color = "#EBEBEB", linewidth = 0.42
    ),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = "#EEEEEE", color = "#777777", linewidth = 0.55
    ),
    strip.text = element_text(
      face = "italic", size = 11, color = "black",
      margin = margin(2, 0, 2, 0)
    ),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title.y = element_text(color = "black", size = 12),
    panel.spacing.x = grid::unit(0.045, "in"),
    plot.tag = element_text(face = "bold", size = 13),
    plot.tag.position = c(0.015, 0.995),
    plot.margin = margin(8, 7, 5, 7)
  )

# Figure significance comes from the donor-median Welch/BH table, not cell-level tests.
bh <- as.data.frame(readxl::read_excel(input, sheet = 'Fig7e_Welch_BH'))
stopifnot(nrow(bh)==3L, all(c('feature','p_BH','BH_star') %in% names(bh)))
stopifnot(identical(as.character(bh$feature), c('CCL3','CCL4','CCL3+CCL4')))
stopifnot(all(is.finite(bh$p_BH)))
expect_star <- function(p) ifelse(p < .001,'***',ifelse(p < .01,'**',ifelse(p < .05,'*','ns')))
stopifnot(identical(as.character(bh$BH_star), as.character(expect_star(bh$p_BH))))
height_by_gene <- tapply(plot_data$Expression, plot_data$Gene, max)
labels <- data.frame(Gene=factor(c('CCL3','CCL4','CCL3+4'),levels=levels(plot_data$Gene)),
  star=as.character(bh$BH_star),
  line_y=pmin(y_top-.9, as.numeric(height_by_gene[c('CCL3','CCL4','CCL3+4')])+.48))
labels$text_y <- labels$line_y+.27
p <- p +
  ggplot2::geom_segment(data=labels, ggplot2::aes(x=1, xend=2, y=line_y, yend=line_y),
                        inherit.aes=FALSE, linewidth=.35, colour='black') +
  ggplot2::geom_text(data=labels, ggplot2::aes(x=1.5, y=text_y, label=star),
                     inherit.aes=FALSE, size=3.5, colour='black')
if (interactive()) print(p)

ggsave(
  file.path(out, "Figure7e_FINAL.png"),
  plot = p, width = 5.6, height = 5.0,
  dpi = 600, bg = "white"
)

ggsave(
  file.path(out, "Figure7e_FINAL.pdf"),
  plot = p, width = 5.6, height = 5.0,
  bg = "white"
)

# 8. Manuscript significance labels are taken directly from Welch/BH results.
metadata <- data.frame(Gene=bh$feature, Stars=as.character(bh$BH_star),
                       Source_p=as.numeric(bh$p_BH),
                       Provenance='Fig7e_Welch_BH',stringsAsFactors=FALSE)

# 9. Sample-level statistics
method <- "Not calculated: no sample-level values"

if (any(!is.na(cells$Sample) & nzchar(cells$Sample))) {
  # Per-sample medians from individual cells
  sample_values <- aggregate(
    cbind(CCL3, CCL4, CCL3_4) ~ Sample + Group,
    data = cells,
    FUN = median
  )
  method <- "Welch test on per-sample MEDIAN expression"
  
} else if (!is.null(support) && nrow(support)) {
  # Separate sample-level summary, if supplied
  sample_values <- support
  method <- "Welch test on per-sample MEAN expression"
  
} else {
  sample_values <- data.frame(
    Note = "No sample-level identifiers or means in this workbook"
  )
}

sample_stats <- data.frame(
  Gene = c("CCL3", "CCL4", "CCL3+4"),
  Control_n = NA_integer_,
  ALI_n = NA_integer_,
  Welch_p = NA_real_,
  BH_adjusted_p = NA_real_,
  Method = method,
  stringsAsFactors = FALSE
)

if (all(c("Sample", "Group", "CCL3", "CCL4", "CCL3_4") %in%
        names(sample_values))) {
  
  gene_cols <- c("CCL3", "CCL4", "CCL3_4")
  
  for (j in seq_along(gene_cols)) {
    g <- gene_cols[j]
    
    ctl <- sample_values[[g]][
      sample_values$Group == "Control"
    ]
    ali <- sample_values[[g]][
      sample_values$Group == "ALI"
    ]
    
    ctl <- ctl[is.finite(ctl)]
    ali <- ali[is.finite(ali)]
    
    sample_stats$Control_n[j] <- length(ctl)
    sample_stats$ALI_n[j] <- length(ali)
    
    if (length(ctl) >= 2 && length(ali) >= 2) {
      sample_stats$Welch_p[j] <- tryCatch(
        t.test(ctl, ali, var.equal = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }
  }
  
  sample_stats$BH_adjusted_p <- p.adjust(
    sample_stats$Welch_p, method = "BH"
  )
}

# 10. Descriptive statistics
cell_summary <- aggregate(
  Expression ~ Gene + Group,
  data = plot_data,
  FUN = function(x) {
    c(
      N = length(x),
      Mean = mean(x),
      Median = median(x),
      SD = sd(x),
      Maximum = max(x)
    )
  }
)

cell_summary <- data.frame(
  Gene = cell_summary$Gene,
  Group = cell_summary$Group,
  cell_summary$Expression
)

# 11. Export Excel
notes <- data.frame(
  Item = c(
    "Figure", "Figure significance", "Manuscript_Labels",
    "Sample_Welch_BH", "Statistical method",
    "Cell_summary", "Source"
  ),
  Explanation = c(
    "Cell-level violin plot of CCL3, CCL4 and CCL3+4",
    "BH-adjusted Welch significance marks drawn",
    "Original manuscript annotation metadata",
    "Welch test and BH correction, when sample values exist",
    method,
    "Descriptive cell-level statistics",
    basename(input)
  )
)

excel_file <- file.path(out, "Figure7e_STATISTICS.xlsx")

# Write a temporary workbook and repair openxlsx worksheet dimensions and
# unused dangling drawing relationships without changing any worksheet values.
if (file.exists(excel_file)) stop("Refusing to overwrite: ", excel_file)
tmp_excel <- tempfile(fileext = ".xlsx")
openxlsx::write.xlsx(
  list(
    Manuscript_Labels = metadata,
    Sample_Welch_BH = sample_stats,
    Sample_Values = sample_values,
    Cell_summary = cell_summary,
    Figure_source = plot_data,
    README = notes
  ),
  file = tmp_excel,
  overwrite = FALSE
)
repair_xlsx(tmp_excel, excel_file)
unlink(tmp_excel)

cat("\n====================================\n")
cat("FIGURE 7e COMPLETED\n")
cat("====================================\n")
cat("PNG:", file.path(out, "Figure7e_FINAL.png"), "\n")
cat("PDF:", file.path(out, "Figure7e_FINAL.pdf"), "\n")
cat("XLSX:", excel_file, "\n")
cat("\nSTATISTICS METHOD:", method, "\n")
print(sample_stats)

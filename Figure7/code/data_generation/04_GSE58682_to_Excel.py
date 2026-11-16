#!/usr/bin/env python3
"""Reproduce Figure 7f/7g data from GEO's deposited GSE58682 series matrix.

This reproduces GEO-deposited processed transcript-cluster expression, NOT raw CEL RMA.
Dependencies: numpy, scipy. Run:
python 04_GSE58682_to_Excel.py path/to/GSE58682_series_matrix.txt.gz output_directory
"""
import csv
import gzip
import math
import re
import sys
from pathlib import Path
from importlib.metadata import version
import numpy as np
from scipy import stats
import zipfile
from xml.sax.saxutils import escape, quoteattr

GENES = {'CCL3': '8014369', 'CCL4': '8006602', 'TNF': '8177983',
         'IL1B': '8054722', 'CCL2': '8006433', 'CD68': '8004510',
         'ITGAM': '7995096', 'SPP1': '8096301'}
GEO_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58682'
COLUMNS_SUB7 = ['sample','group','dose_ppb','CCL3_log2','CCL4_log2','zCCL3','zCCL4',
                'CoExpr_Zmean','AMInflam_Score','AMInflam_nGenes']

def read_series_matrix(path):
    metadata = {}
    expr_rows = {}
    in_table = False
    sample_gsms = None
    with gzip.open(path, 'rt') as f:
        reader = csv.reader(f, delimiter='\t', quotechar='"')
        for row in reader:
            if not row: continue
            if row[0] == '!series_matrix_table_begin':
                in_table = True
                continue
            if row[0] == '!series_matrix_table_end':
                in_table = False
                break
            if not in_table:
                if row[0].startswith('!Sample_'):
                    metadata.setdefault(row[0], []).append(row[1:])
                continue
            if row[0] == 'ID_REF':
                sample_gsms = row[1:]
                continue
            if row[0] in set(GENES.values()):
                if row[0] in expr_rows: raise ValueError(f'Duplicate probe ID: {row[0]}')
                expr_rows[row[0]] = {gsm: float(value) for gsm, value in zip(sample_gsms, row[1:])}
    assert sample_gsms is not None and len(sample_gsms) == 56
    assert set(expr_rows) == set(GENES.values())
    gsm_metadata = metadata['!Sample_geo_accession'][0]
    assert gsm_metadata == sample_gsms, 'Matrix columns != GEO metadata GSM order'
    def metarow(key, prefix):
        matched = [r for r in metadata[key] if all(s.startswith(prefix) for s in r)]
        if len(matched) != 1: raise ValueError(f'No unique {key}/{prefix}')
        return matched[0]
    subjects = metarow('!Sample_characteristics_ch1', 'subject id:')
    doses = metarow('!Sample_characteristics_ch1', 'ozone exposure (ppb):')
    entries = []
    for gsm, subj, dose in zip(sample_gsms, subjects, doses):
        subj_number = int(re.search(r'(\d+)$', subj).group(1))
        dose_number = int(re.search(r'(\d+)$', dose).group(1))
        entries.append((gsm, subj_number, dose_number,
                        f'{gsm}_{subj_number:02d}_{dose_number:03d}'))
    entries.sort(key=lambda s: (s[2], int(s[0][3:])))
    assert len({e[3] for e in entries}) == 56
    assert [sum(e[2] == q for e in entries) for q in (0,100,200)] == [18,19,19]
    ex = np.array([[expr_rows[probe][e[0]] for probe in GENES.values()]
                   for e in entries], dtype=float)
    assert ex.shape == (56, 8) and np.isfinite(ex).all()
    return entries, ex

def zscore(x):
    return (x - np.mean(x, axis=0)) / np.std(x, axis=0, ddof=1)

# XLSX export using the Python standard library only (ECMA-376 OOXML).
def _excel_col(n):
    s = ''
    while n > 0:
        n, rem = divmod(n - 1, 26)
        s = chr(65 + rem) + s
    return s

def _cell(ref, value, is_header=False):
    if isinstance(value, np.generic):
        value = value.item()
    if value is None:
        return ''
    sf = ' s="1"' if is_header else ''
    if isinstance(value, bool):
        return f'<c r="{ref}"{sf} t="b"><v>{int(value)}</v></c>'
    if isinstance(value, (float, int)):
        if not math.isfinite(float(value)):
            raise ValueError(f'Non-finite number at {ref}')
        return f'<c r="{ref}"{sf}><v>{value!r}</v></c>'
    return f'<c r="{ref}"{sf} t="inlineStr"><is><t>{escape(str(value))}</t></is></c>'

def save_xlsx(path, sheets):
    ns = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    typ = 'http://schemas.openxmlformats.org/package/2006/content-types'
    relns = 'http://schemas.openxmlformats.org/package/2006/relationships'
    odrel = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
    wb_sheets, wb_rels, overrides = [], [], []
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for i, (name, headers, rows) in enumerate(sheets, 1):
            if len(name) > 31 or any(len(row) != len(headers) for row in rows):
                raise ValueError('Invalid worksheet: ' + name)
            wb_sheets.append(f'<sheet name={quoteattr(name)} sheetId="{i}" r:id="rId{i}"/>')
            wb_rels.append(f'<Relationship Id="rId{i}" Type="{odrel}/worksheet" Target="worksheets/sheet{i}.xml"/>')
            overrides.append(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
            rr = []
            for row_i, line in enumerate([headers] + list(rows), 1):
                cells = ''.join(_cell(f'{_excel_col(col_i)}{row_i}', v, row_i == 1)
                                for col_i, v in enumerate(line, 1))
                rr.append(f'<row r="{row_i}">{cells}</row>')
            content = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                       f'<worksheet xmlns="{ns}"><dimension ref="A1:{_excel_col(len(headers))}{len(rows)+1}"/>'
                       '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
                       'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
                       '<sheetFormatPr defaultRowHeight="15"/><sheetData>'
                       + ''.join(rr) + '</sheetData></worksheet>')
            z.writestr(f'xl/worksheets/sheet{i}.xml', content)
        wb_rels.append(f'<Relationship Id="rId{len(sheets)+1}" Type="{odrel}/styles" Target="styles.xml"/>')
        z.writestr('xl/workbook.xml', f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            f'<workbook xmlns="{ns}" xmlns:r="{odrel}"><sheets>{"".join(wb_sheets)}</sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels', f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            f'<Relationships xmlns="{relns}">{"".join(wb_rels)}</Relationships>')
        z.writestr('xl/styles.xml', f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            f'<styleSheet xmlns="{ns}"><fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
            '<font><b/><color rgb="FFFFFFFF"/><sz val="11"/><name val="Calibri"/></font></fonts>'
            '<fills count="3"><fill><patternFill patternType="none"/></fill>'
            '<fill><patternFill patternType="gray125"/></fill>'
            '<fill><patternFill patternType="solid"><fgColor rgb="FF1C3859"/><bgColor indexed="64"/>'
            '</patternFill></fill></fills><borders count="1"><border><left/><right/><top/><bottom/><diagonal/>'
            '</border></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/>'
            '</cellStyleXfs><cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
            '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
            '</cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/>'
            '</cellStyles></styleSheet>')
        z.writestr('[Content_Types].xml', f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            f'<Types xmlns="{typ}"><Default Extension="rels" '
            'ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/styles.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            + ''.join(overrides) + '</Types>')
        z.writestr('_rels/.rels', f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            f'<Relationships xmlns="{relns}"><Relationship Id="rId1" '
            f'Type="{odrel}/officeDocument" Target="xl/workbook.xml"/></Relationships>')

def main():
    source = Path(sys.argv[1]) if len(sys.argv)>1 else Path('GSE58682_series_matrix.txt.gz')
    out = Path(sys.argv[2]) if len(sys.argv)>2 else Path('generated')
    out.mkdir(exist_ok=True, parents=True)
    outputs = ('Supplementary_Data_7_GSE58682_reproduced.xlsx',
               'Source_Data_Figure_7fg_extract.xlsx',
               'Figure_7f_7g_plotting_data_GSE58682.xlsx', 'execution_report.txt')
    existing = [str(out / name) for name in outputs if (out / name).exists()]
    if existing:
        raise FileExistsError('Use a new output folder; refusing to overwrite: ' + ', '.join(existing))
    entries, expr = read_series_matrix(source)
    z = zscore(expr)
    coexpr = z[:,:2].mean(axis=1)
    activation = z.mean(axis=1)
    zact = zscore(activation)
    zco = zscore(coexpr)
    r9 = lambda x: np.round(x,9)
    g = np.array([e[2] for e in entries])
    ids = [e[3] for e in entries]
    clean = g==0; ozone = g>0
    assert (clean.sum(), ozone.sum()) == (18,38)
    assert abs(expr[0,0]-8.077609)<1e-10 and abs(expr[0,1]-6.040576)<1e-10
    assert abs(activation[0]-0.210292418)<1e-8
    r_sub = []
    r_z = []
    r_plot = []
    r_sf = []; r_sg=[]
    for i in range(56):
        r_sub.append([ids[i],'CleanAir' if clean[i] else 'Ozone', None if clean[i] else int(g[i]),
                      float(expr[i,0]),float(expr[i,1]),float(r9(z[i,0])),float(r9(z[i,1])),
                      float(r9(coexpr[i])),float(r9(activation[i])),8])
        r_z.append([*map(float,r9(z[i])), ids[i]])
        r_plot.append([ids[i], 'Clean Air' if clean[i] else 'Ozone',int(g[i]),
                       float(expr[i,0]),float(expr[i,1]),float(r9(coexpr[i])),float(r9(activation[i])),
                       float(zact[i]),float(r9(z[i,0])),float(r9(z[i,1])),float(zco[i])])
        for j, gene in enumerate(['CCL3','CCL4','CCL3+CCL4']):
            y = [expr[i,0], expr[i,1], coexpr[i]][j]
            yn = [z[i,0], z[i,1], zco[i]][j]
            r_sf.append(['7f',gene,ids[i],'CleanAir' if clean[i] else 'Ozone',int(g[i]),
                         float(r9(y)) if j==2 else float(y),'mean z-score' if j==2 else 'log2 normalized expression'])
            r_sg.append(['7g',gene,ids[i],float(r9(activation[i])),float(r9(y)) if j==2 else float(y),
                         float(zact[i]),float(r9(yn)) if j<2 else float(yn)])
    feature_order={'CCL3':0,'CCL4':1,'CCL3+CCL4':2}
    r_sf.sort(key=lambda row: feature_order[row[1]])
    r_sg.sort(key=lambda row: feature_order[row[1]])
    save_xlsx(out/'Supplementary_Data_7_GSE58682_reproduced.xlsx', [
        ('Raw data',COLUMNS_SUB7,r_sub),
        ('AM_genes_found',['gene'],[[gene] for gene in GENES]),
        ('AMgene_zscores',['z.'+gene for gene in GENES]+['sample'],r_z)])
    clean_idx=np.where(clean)[0]; ozone_idx=np.where(ozone)[0]
    f_names=['CCL3_CleanAir_log2','CCL3_Ozone_log2','CCL4_CleanAir_log2','CCL4_Ozone_log2',
             'CCL3+CCL4_CleanAir_meanZ','CCL3+CCL4_Ozone_meanZ']
    f_rows=[]
    for i in range(38):
        def take(x, idx): return float(x[idx[i]]) if i < len(idx) else None
        f_rows.append([take(expr[:,0],clean_idx),take(expr[:,0],ozone_idx),
                       take(expr[:,1],clean_idx),take(expr[:,1],ozone_idx),
                       take(r9(coexpr),clean_idx),take(r9(coexpr),ozone_idx)])
    g_headers=['X_AMactivation_z_CCL3','Y_CCL3_z','X_AMactivation_z_CCL4','Y_CCL4_z',
               'X_AMactivation_z_CCL3+CCL4','Y_CCL3+CCL4_z']
    g_rows=[[float(zact[i]),float(r9(z[i,0])),float(zact[i]),float(r9(z[i,1])),
             float(zact[i]),float(zco[i])] for i in range(56)]
    stats_rows=[]
    for name, y in [('CCL3',expr[:,0]),('CCL4',expr[:,1]),('CCL3+CCL4 mean z',coexpr)]:
        t = stats.ttest_ind(y[ozone],y[clean],equal_var=False)
        va=np.var(y[ozone],ddof=1)/38; vb=np.var(y[clean],ddof=1)/18
        df=(va+vb)**2/(va**2/37+vb**2/17)
        md=float(np.mean(y[ozone])-np.mean(y[clean]))
        ci=stats.t.interval(0.95,df,loc=md,scale=math.sqrt(va+vb))
        stats_rows.append(['7f',name,"Welch's two-sided t-test",f't({df:.3f}) = {t.statistic:.4f}',
                           float(t.pvalue),'*' if t.pvalue<.05 else 'ns','Clean Air 18; Ozone 38',
                           f'Difference Ozone-Clean Air = {md:.6f}; 95% CI [{ci[0]:.6f}, {ci[1]:.6f}]'])
    # Two-response one-way Pillai MANOVA, explanatory only; not plotted on mean-z boxplot.
    Y=expr[:,:2]; ym=Y.mean(axis=0)
    E=np.zeros((2,2)); H=np.zeros((2,2))
    for mask in (clean,ozone):
        V=Y[mask]; mu=V.mean(axis=0)
        E += (V-mu).T@(V-mu)
        d=(mu-ym).reshape(2,1)
        H += mask.sum()*(d@d.T)
    pillai=float(np.trace(np.linalg.solve(E+H,H)))
    F=float((53/2)*pillai/(1-pillai)); p_man=float(stats.f.sf(F,2,53))
    stats_rows.append(['7f','CCL3 & CCL4 joint','Two-response MANOVA, Pillai',
                       f'V = {pillai:.6f}; F(2,53) = {F:.4f}',p_man,'not the mean-z boxplot test',
                       'Clean Air 18; Ozone 38', 'Multivariate test; separate from univariate boxplot'])
    corr=[]
    for name, y in [('CCL3',z[:,0]),('CCL4',z[:,1]),('CCL3+CCL4',zco)]:
        p=stats.pearsonr(zact,y)
        r=float(p.statistic); pv=float(p.pvalue)
        ci=p.confidence_interval(0.95)
        corr.append((name,r,pv))
        stats_rows.append(['7g','AM activation vs '+name,'Pearson r, two-sided',f'r = {r:.6f}',pv,
                           f'R = {r:.2f}; P = {pv:.3g}','56 profiles',
                           f'95% CI of r [{ci.low:.4f}, {ci.high:.4f}]'])
    notes=[
       ['METHOD','Input','GEO processed Series Matrix','NOT raw CEL RMA',None,None,None,GEO_URL],
       ['METHOD','CCL3+CCL4','mean(zCCL3,zCCL4)','56-sample z scores',None,None,None,'Boxplot P refers to mean-z t test; MANOVA not a boxplot P'],
       ['METHOD','AM activation score','mean 8 gene z scores','56-sample z scores',None,None,None,'CCL3,CCL4,TNF,IL1B,CCL2,CD68,ITGAM,SPP1'],
       ['LIMITATION','Part-whole relation','CCL3/4 included in score',None,None,None,None,'Pearson correlation is not independent of gene inclusion']
    ]
    save_xlsx(out/'Figure_7f_7g_plotting_data_GSE58682.xlsx', [
      ('Fig7f_Boxplot',f_names,f_rows),('Fig7g_Scatter',g_headers,g_rows),
      ('Sample_Data',['Sample_ID','Exposure_group','Ozone_ppb','CCL3_log2','CCL4_log2',
                      'CCL3+CCL4_meanZ','AM_activation_score','AM_activation_score_z','CCL3_z',
                      'CCL4_z','CCL3+CCL4_z'],r_plot),
      ('Stats_Notes',['Panel','Comparison','Test','Statistic','P-value','Figure label','N','Interpretation'],
          stats_rows+notes)])
    save_xlsx(out/'Source_Data_Figure_7fg_extract.xlsx', [
       ('Fig7f',['figure','feature','sample','group','dose_ppb','expression','unit'],r_sf),
       ('Fig7g',['figure','feature','sample','activation_score','gene_expression','activation_score_z','expression_z'],r_sg),
       ('Statistics',['Panel','Comparison','Test','Statistic','P-value','Figure label','N','Interpretation'],stats_rows)])
    # Core explicit reference checks: source values, tests, every profile gene score.
    assert np.max(abs(r9(z).mean(axis=1)-r9(activation))) < 1e-8
    assert abs(stats_rows[0][4]-0.0335066754)<1e-6
    assert abs(stats_rows[1][4]-0.2233744639)<1e-6
    assert abs(stats_rows[2][4]-0.4833900961)<1e-6
    assert abs(p_man-0.008428645)<1e-5
    assert np.max(abs(np.array([x[7] for x in r_sub])-r9(coexpr))) < 1e-10
    assert all(abs(r-expected)<1e-6 for (n,r,p),expected in zip(corr,[.2752292803,.5744039878,.4970710287]))
    (out/'execution_report.txt').write_text(
        'EXECUTED GSE58682 data generation successfully\nInput: '+str(source)+'\nGEO profiles: 56 (CleanAir=18, Ozone=38)\n'
        +'CCL3='+str(expr[0,0])+', CCL4='+str(expr[0,1])+' (GSM1416848)\n'
        +'Figure 7f Welch p-values: '+', '.join(f'{r[1]}={r[4]:.9g}' for r in stats_rows[:3])+'\n'
        +f'Joint Pillai MANOVA p={p_man:.9g}\n'
        +'Figure 7g Pearson: '+', '.join(f'{n}: r={r:.9g} p={p:.9g}' for n,r,p in corr)+'\n'
        +f'numpy={np.__version__}; scipy={version("scipy")}\n',
        encoding='utf-8')
    print((out/'execution_report.txt').read_text())

if __name__=='__main__': main()

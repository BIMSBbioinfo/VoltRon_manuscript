"""Compile per-step time / memory results from the method results.txt files
into a single .xlsx workbook (long table + per-run summary).

Values are read manually from:
  wsireg/{LungTMA,BreastCancer,Xeniumpups}/results*.txt
  STalign/{LungTMA,BreastCancer,Xeniumpups}/results.txt
  {PASTE,PASTE2,SANTO}/DLPFC/results.txt

No third-party xlsx writer is available in this environment, so the workbook
is assembled directly as the (zip + XML) OOXML container.
"""

import zipfile
from xml.sax.saxutils import escape

# Method, Dataset, Run, Step order, Step, Time (s), Current mem (MiB), Peak mem (MiB), Source
ROWS = [
    # --- wsireg: image-based registration (elastix) -------------------------
    ("wsireg", "LungTMA",      "cell coordinates",     1, "Image registration", 105.3640, 192.97,  519.40, "wsireg/LungTMA/results.txt"),
    ("wsireg", "LungTMA",      "cell coordinates",     2, "Transform images",     9.1819, 173.25,  178.00, "wsireg/LungTMA/results.txt"),
    ("wsireg", "LungTMA",      "cell coordinates",     3, "Transform points",    21.8412,   0.01,    7.42, "wsireg/LungTMA/results.txt"),
    ("wsireg", "LungTMA",      "molecule coordinates", 1, "Image registration", 226.0471, 192.41,  519.41, "wsireg/LungTMA/results_mol.txt"),
    ("wsireg", "LungTMA",      "molecule coordinates", 2, "Transform images",     9.8675, 173.25,  178.01, "wsireg/LungTMA/results_mol.txt"),
    ("wsireg", "LungTMA",      "molecule coordinates", 3, "Transform points",    65.9195,   0.03,  413.28, "wsireg/LungTMA/results_mol.txt"),
    ("wsireg", "BreastCancer", "default",              1, "Image registration", 138.1781, 195.73,  339.63, "wsireg/BreastCancer/results.txt"),
    ("wsireg", "BreastCancer", "default",              2, "Transform images",    14.1244,  89.59,  338.56, "wsireg/BreastCancer/results.txt"),
    ("wsireg", "BreastCancer", "default",              3, "Transform points",    38.9999,   0.03,   77.21, "wsireg/BreastCancer/results.txt"),
    ("wsireg", "Xeniumpups",   "default",              1, "Image registration", 168.2173, 199.80, 1462.46, "wsireg/Xeniumpups/results.txt"),
    ("wsireg", "Xeniumpups",   "default",              2, "Transform images",    91.6777, 497.35, 1452.56, "wsireg/Xeniumpups/results.txt"),
    ("wsireg", "Xeniumpups",   "default",              3, "Transform points",   263.7011,   0.03,  595.75, "wsireg/Xeniumpups/results.txt"),

    # --- STalign: landmark/LDDMM point registration -------------------------
    ("STalign", "LungTMA",      "default", 1, "Point registration", 146.7296, 35.81,  54.80, "STalign/LungTMA/results.txt"),
    ("STalign", "LungTMA",      "default", 2, "Transform points",     0.0033,  0.00,   0.19, "STalign/LungTMA/results.txt"),
    ("STalign", "BreastCancer", "default", 1, "Point registration", 168.9887, 57.83,  64.07, "STalign/BreastCancer/results.txt"),
    ("STalign", "BreastCancer", "default", 2, "Transform points",     0.0208,  0.00,   2.56, "STalign/BreastCancer/results.txt"),
    ("STalign", "Xeniumpups",   "default", 1, "Point registration", 336.5115, 73.34, 206.58, "STalign/Xeniumpups/results.txt"),
    ("STalign", "Xeniumpups",   "default", 2, "Transform points",     0.1556,  0.00,  19.82, "STalign/Xeniumpups/results.txt"),

    # --- spot-level slice alignment ----------------------------------------
    ("PASTE",  "DLPFC", "default", 1, "Slice alignment",  27.2988, 102.16, 1933.54, "PASTE/DLPFC/results.txt"),
    ("PASTE2", "DLPFC", "default", 1, "Slice alignment", 787.6818, 101.94, 1719.27, "PASTE2/DLPFC/results.txt"),
    ("SANTO",  "DLPFC", "default", 1, "Slice alignment", 798.4196,   2.30, 1661.58, "SANTO/DLPFC/results.txt"),
]

# Dataset / method context used as extra grouping columns.
MODALITY = {
    "LungTMA":      "Xenium + H&E / IF (LungTMA assay 4)",
    "BreastCancer": "Xenium + H&E (10x breast cancer)",
    "Xeniumpups":   "Xenium + H&E (mouse pups)",
    "DLPFC":        "Visium (DLPFC adjacent sections)",
}
APPROACH = {
    "wsireg":  "Image-based registration",
    "STalign": "Image/point-based registration",
    "PASTE":   "Spot-level slice alignment",
    "PASTE2":  "Spot-level slice alignment",
    "SANTO":   "Spot-level slice alignment",
}

LONG_HEADER = [
    "Method", "Approach", "Dataset", "Modality", "Run", "Step order", "Step",
    "Time (s)", "Time (min)", "Current memory (MiB)", "Peak memory (MiB)",
    "Source file",
]

SUMMARY_HEADER = [
    "Method", "Approach", "Dataset", "Modality", "Run", "Number of steps",
    "Total time (s)", "Total time (min)", "Max peak memory (MiB)", "Source file",
]


def long_rows():
    out = []
    for m, ds, run, order, step, t, cur, peak, src in ROWS:
        out.append([m, APPROACH[m], ds, MODALITY[ds], run, order, step,
                    t, round(t / 60.0, 4), cur, peak, src])
    return out


def summary_rows():
    groups = {}
    order = []
    for m, ds, run, _o, _s, t, _c, peak, src in ROWS:
        key = (m, ds, run)
        if key not in groups:
            groups[key] = {"n": 0, "t": 0.0, "peak": 0.0, "src": set()}
            order.append(key)
        g = groups[key]
        g["n"] += 1
        g["t"] += t
        g["peak"] = max(g["peak"], peak)
        g["src"].add(src)
    out = []
    for m, ds, run in order:
        g = groups[(m, ds, run)]
        out.append([m, APPROACH[m], ds, MODALITY[ds], run, g["n"],
                    round(g["t"], 4), round(g["t"] / 60.0, 4), g["peak"],
                    "; ".join(sorted(g["src"]))])
    return out


# --------------------------------------------------------------------------
# minimal xlsx writer
# --------------------------------------------------------------------------
def col_name(i):
    name = ""
    while i >= 0:
        name = chr(ord("A") + i % 26) + name
        i = i // 26 - 1
    return name


def sheet_xml(header, rows):
    parts = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
        '<sheetViews><sheetView workbookViewId="0">'
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
        '</sheetView></sheetViews>',
        "<sheetData>",
    ]
    all_rows = [header] + rows
    for r, row in enumerate(all_rows, start=1):
        cells = []
        for c, val in enumerate(row):
            ref = "%s%d" % (col_name(c), r)
            style = ' s="1"' if r == 1 else ""
            if isinstance(val, (int, float)) and not isinstance(val, bool):
                cells.append('<c r="%s"%s><v>%s</v></c>' % (ref, style, val))
            else:
                cells.append('<c r="%s"%s t="inlineStr"><is><t>%s</t></is></c>'
                             % (ref, style, escape(str(val))))
        parts.append("<row r=\"%d\">%s</row>" % (r, "".join(cells)))
    parts.append("</sheetData></worksheet>")
    return "".join(parts)


def write_xlsx(path, sheets):
    content_types = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
    ]
    for i in range(len(sheets)):
        content_types.append(
            '<Override PartName="/xl/worksheets/sheet%d.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' % (i + 1))
    content_types.append("</Types>")

    wb = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
          '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>']
    for i, (name, _h, _r) in enumerate(sheets):
        wb.append('<sheet name="%s" sheetId="%d" r:id="rId%d"/>' % (escape(name), i + 1, i + 1))
    wb.append("</sheets></workbook>")

    wb_rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
               '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
    for i in range(len(sheets)):
        wb_rels.append('<Relationship Id="rId%d" '
                       'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                       'Target="worksheets/sheet%d.xml"/>' % (i + 1, i + 1))
    wb_rels.append('<Relationship Id="rId%d" '
                   'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
                   'Target="styles.xml"/>' % (len(sheets) + 1))
    wb_rels.append("</Relationships>")

    root_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                 '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                 '<Relationship Id="rId1" '
                 'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
                 'Target="xl/workbook.xml"/></Relationships>')

    styles = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
              '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
              '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
              '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
              '<fills count="2"><fill><patternFill patternType="none"/></fill>'
              '<fill><patternFill patternType="gray125"/></fill></fills>'
              '<borders count="1"><border/></borders>'
              '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
              '<cellXfs count="2">'
              '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
              '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/>'
              '</cellXfs></styleSheet>')

    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", "".join(content_types))
        z.writestr("_rels/.rels", root_rels)
        z.writestr("xl/workbook.xml", "".join(wb))
        z.writestr("xl/_rels/workbook.xml.rels", "".join(wb_rels))
        z.writestr("xl/styles.xml", styles)
        for i, (_n, header, rows) in enumerate(sheets):
            z.writestr("xl/worksheets/sheet%d.xml" % (i + 1), sheet_xml(header, rows))


if __name__ == "__main__":
    out = "registration_speed_benchmark_other_methods.xlsx"
    write_xlsx(out, [
        ("per_step", LONG_HEADER, long_rows()),
        ("per_run_summary", SUMMARY_HEADER, summary_rows()),
    ])
    print("wrote", out, "-", len(long_rows()), "step rows,", len(summary_rows()), "summary rows")

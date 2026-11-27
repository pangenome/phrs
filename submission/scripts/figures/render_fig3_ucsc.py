#!/usr/bin/env python3
"""Fig. 3 UCSC panels re-rendered as in May 2026 (slides/chm13-phr-ucsc-browser, render_ucsc_browser_suite.py), with a wider
label column so that the side label "SEDEF Segmental Dups" (20 characters) is no longer cut ("EDEF Segmental Dups").

usage (repo root): python3 submission/scripts/figures/render_fig3_ucsc.py OUT_DIR [--label-width N] [--text-font NAME]
A new UCSC session (the May session 3966979908_... was reused since and no longer has the custom track or the same tracks):
  step 1: load data/chm13.phrs.bed (checked equal to the public GitHub copy) as a custom track named chm13.phrs.bed, as in May;
  step 2: per panel, the May configuration read from the May HTML (hgTracks JSON): ruler, chm13.phrs.bed, GC Percent,
  RepeatMasker, SEDEF Segmental Dups dense, CAT/Liftoff Genes pack, in that image order (imgOrd 0-5), same windows,
  hgt.labelWidth N (default 22), pix adjusted so that the data image keeps the May width of 1645 px; center labels off, Helvetica
  normal, text size 12 (May: label width 20 and text size 12 give its row heights 35/16/176 and data width 1645 exactly).
Panel = side-label crop + data crop of the UCSC sprites (CSS offsets), appended, 1 px #cfcfcf border (as the May script).
Checks: track rows and order as in May, data width 1645, then a pixel comparison of the data part with the May panel.
out: OUT_DIR/<Fig3 panel name>.png, OUT_DIR/compare.tsv; log on stderr.
"""
import html, json, re, subprocess, sys, tempfile
from pathlib import Path
from PIL import Image, ImageChops
OUT = Path(sys.argv[1]); LW = int(sys.argv[sys.argv.index("--label-width") + 1]) if "--label-width" in sys.argv else 22
FONT = sys.argv[sys.argv.index("--text-font") + 1] if "--text-font" in sys.argv else "Helvetica"
OUT.mkdir(parents=True, exist_ok=False)
BASE = "https://genome.ucsc.edu/cgi-bin/hgTracks"; DB = "hub_3671779_hs1"; PIX = 1800 + 7 * (LW - 20)
RAW = "https://raw.githubusercontent.com/pangenome/phrs/main/data/chm13.phrs.bed"
MAY = Path("slides/chm13-phr-ucsc-browser/_assets/ucsc/panels")
PANELS = [("Fig3a_C1_chr4q", "chr4", 193304946, 193574945, "chr4q_chr4_193304946_193574945"),
          ("Fig3a_C1_chr10q", "chr10", 134488135, 134758134, "chr10q_chr10_134488135_134758134"),
          ("Fig3b_C2_chr10p", "chr10", 1, 105002, "chr10p_chr10_1_105002"),
          ("Fig3b_C2_chr18p", "chr18", 1, 397502, "chr18p_chr18_1_397502"),
          ("Fig3c_C11_chr5q", "chr5", 181752940, 182045439, "chr5q_chr5_181752940_182045439"),
          ("Fig3c_C11_chr6q", "chr6", 171909129, 172126628, "chr6q_chr6_171909129_172126628")]
HUB = ["hub_3671779_gc5Base", "hub_3671779_t2tRepeatMasker", "hub_3671779_sedefSegDups", "hub_3671779_catLiftOffGenesV1"]
def log(m): print(m, file=sys.stderr, flush=True)
tmp = Path(tempfile.mkdtemp(prefix="fig3ucsc_", dir="/tmp")); jar = tmp / "cookies.txt"
def post(params, dest):
    cmd = ["curl", "-fsSL", "-A", "Mozilla/5.0", "--compressed", "--retry", "3", "--connect-timeout", "30", "-c", str(jar), "-b", str(jar), BASE]
    for k, v in params.items(): cmd += ["--data-urlencode", f"{k}={v}"]
    subprocess.run(cmd + ["-o", str(dest)], check=True); return dest.read_text(errors="replace")
def get(url, dest): subprocess.run(["curl", "-fsSL", "-A", "Mozilla/5.0", "--retry", "3", "-o", str(dest), url], check=True)
bed = subprocess.run(["curl", "-fsSL", RAW], check=True, capture_output=True, text=True).stdout
if bed != Path("data/chm13.phrs.bed").read_text(): sys.exit("ERROR the GitHub copy differs from data/chm13.phrs.bed")
t = post({"db": DB, "position": "chr4:193304946-193574945", "hgt.customText": "track name=chm13.phrs.bed\n" + bed}, tmp / "step1.html")  # named as in May
sid = re.search(r"hgsid=(\d+_[A-Za-z0-9]+)", t).group(1)
ct = sorted(set(re.findall(r"tr_(ct_chm13phrsbed_\d+)", t)))
if len(ct) != 1: sys.exit(f"ERROR custom track rows {ct}")
ct = ct[0]; log(f"session created, custom track {ct}, label width {LW}, pix {PIX}")
order = ["ruler", ct] + HUB
rows_out = ["panel\tmay_size\tnew_size\tdata_pixels\tdata_identical\tdata_identical_fraction"]
for name, chrom, a, b, may in PANELS:
    p = {"hgsid": sid, "db": DB, "position": f"{chrom}:{a}-{b}", "pix": PIX, "hideTracks": 1, "hgt.labelWidth": LW,
         "ruler": "dense", ct: "dense", "hub_3671779_gc5Base": "dense", "hub_3671779_t2tRepeatMasker": "dense",
         "hub_3671779_sedefSegDups": "dense", "hub_3671779_catLiftOffGenesV1": "pack", "centerLabels": 0, "textFont": FONT, "textStyle": "Normal", "textSize": 12}
    p.update({f"{tr}_imgOrd": i for i, tr in enumerate(order)})
    for _ in range(3):  # pix is the whole image width: adjust it until the data image is 1645 px wide, as in May
        t = post(p, tmp / f"{name}.html")
        dw = int(re.search(r"width:(\d+)px; height:\d+px;' class='sliceDiv dataImg", t).group(1))
        if dw == 1645: break
        p["pix"] = int(p["pix"]) + 1645 - dw; log(f"{name}: data width {dw}, pix -> {p['pix']}")
    got = re.findall(r"<TR id='tr_([A-Za-z0-9_]+)'", t)
    if got != order: sys.exit(f"ERROR {name}: rows {got}, expected {order}")
    sl = {}
    for kind in ("sideLab", "dataImg"):
        m = re.findall(rf"<div style='width:(\d+)px; height:(\d+)px;' class='sliceDiv {kind}[^']*'>.*?<IMG\b[^>]*\bsrc='([^']+)'[^>]*\bstyle='([^']+)'", t, re.S)
        w = {int(x[0]) for x in m}; src = {x[2] for x in m}
        xs = {-int(re.search(r"left:\s*(-?\d+)px", html.unescape(x[3])).group(1)) for x in m}
        hs = [-int(re.search(r"top:\s*(-?\d+)px", html.unescape(x[3])).group(1)) + int(x[1]) for x in m]
        if len(w) != 1 or len(src) != 1 or len(xs) != 1: sys.exit(f"ERROR {name} {kind}: widths {w}, sources {src}, x {xs}")
        f = tmp / f"{name}_{kind}.png"; get("https://genome.ucsc.edu/cgi-bin/" + html.unescape(src.pop()), f)
        sl[kind] = (Image.open(f).convert("RGB"), xs.pop(), w.pop(), max(hs))
    if sl["dataImg"][2] != 1645: sys.exit(f"ERROR {name}: data width {sl['dataImg'][2]}, May 1645")
    h = max(sl["sideLab"][3], sl["dataImg"][3])
    side = sl["sideLab"][0].crop((sl["sideLab"][1], 0, sl["sideLab"][1] + sl["sideLab"][2], h))
    data = sl["dataImg"][0].crop((sl["dataImg"][1], 0, sl["dataImg"][1] + sl["dataImg"][2], h))
    pan = Image.new("RGB", (side.width + data.width + 2, h + 2), (0xcf, 0xcf, 0xcf))
    pan.paste(side, (1, 1)); pan.paste(data, (1 + side.width, 1)); pan.save(OUT / f"{name}.png")
    old = Image.open(MAY / f"{may}.png").convert("RGB"); ow = old.width - 2 - 1645
    od = old.crop((1 + ow, 1, 1 + ow + 1645, old.height - 1)); nd = data
    n = min(od.height, nd.height); diff = ImageChops.difference(od.crop((0, 0, 1645, n)), nd.crop((0, 0, 1645, n)))
    same = sum(1 for px in diff.get_flattened_data() if px == (0, 0, 0)); tot = 1645 * n
    rows_out.append(f"{name}\t{old.width}x{old.height}\t{pan.width}x{pan.height}\t{tot}\t{same}\t{same / tot:.6f}")
    log(f"{name}: May {old.width}x{old.height}, new {pan.width}x{pan.height}; data part identical pixels {same}/{tot} ({same / tot:.4%}), heights May {od.height} new {nd.height}")
(OUT / "compare.tsv").write_text("\n".join(rows_out) + "\n")
log(f"wrote {OUT}")

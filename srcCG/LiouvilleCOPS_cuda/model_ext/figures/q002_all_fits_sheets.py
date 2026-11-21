# viewing contact sheets (scratch) for the full-set clusters of q002_all_fits_summary_docx.py: one JPEG per cluster,
# ==== version-controlled copy: header added 2026-11-20, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/q002_all_fits_20261120/sheets_fullsets.py (on Voyager: /data001/highresolu_vt/cops_output/tests/q002_all_fits_20261120/sheets_fullsets.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: b4c6b7a2bf6257f6659835e547a2acd039d92ce81184619b8084ae2ed7da0ddc
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: viewing contact sheets of the full-set clusters for q002_all_fits_summary_docx.py (observations from the images) (2026-11-20)
# RUN: srun ... docx_venv/bin/python -I <this> <out dir> [sets spec json] (HOW-TO-RERUN.md of figures_cluster/12M/q002_all_fits_summary/)
# ==== end of header ====
# the ten runs in the fixed order (uses the generator's discovery / catalogue)
#   python -I sheets_fullsets.py <out dir> [sets spec json]
import importlib.util, os, sys
from concurrent.futures import ProcessPoolExecutor
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None
spec = importlib.util.spec_from_file_location("gen", "/data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_all_fits_summary_docx.py")
gen = importlib.util.module_from_spec(spec); spec.loader.exec_module(gen)
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
class A: pass
a = A(); a.sets_spec = sys.argv[2] if len(sys.argv) > 2 else None
groups = [g for g in gen.setup(a) if g["section"] in gen.FS_SECTIONS]
def one(g, W=2500, cols=5):
    tw = W // cols; tiles = []
    for r in gen.RK:
        p = g["paths"].get(r)
        if p:
            im = Image.open(p).convert("RGB"); tiles.append((r, im.resize((tw, round(im.size[1] * tw / im.size[0])), Image.LANCZOS)))
        else:
            tiles.append((r + " (none)", Image.new("RGB", (tw, 60), "white")))
    th = max(t[1].size[1] for t in tiles) + 40; rows = (len(tiles) + cols - 1) // cols
    sheet = Image.new("RGB", (W, rows * th), "white"); dr = ImageDraw.Draw(sheet)
    font = ImageFont.truetype("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans-Bold.ttf", 28)
    for i, (k, im) in enumerate(tiles):
        x = (i % cols) * tw; y = (i // cols) * th; dr.text((x + 10, y + 4), k, fill="red", font=font); sheet.paste(im, (x, y + 40))
    name = g["cid"].replace("fs:", "").replace("/", "__").replace("<STEM>", "STEM").replace(".png", ".jpg")
    sheet.save(os.path.join(OUT, name), "JPEG", quality=80); return name
if __name__ == "__main__":
    with ProcessPoolExecutor(8) as ex:
        for n in ex.map(one, groups): print(n)

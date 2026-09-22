"""Render papers/paperF_hubble_meter.md to HTML + PDF for the send
package (the P1/P2 colleague_package template, same CSS). The render
strips ONLY internal build apparatus, mechanically and disclosed here:
the bold draft-status banner, [ASSEMBLY ...] tags, the bracketed
in-text verification note, and the reference-list (v)/[...] tags.
Nothing else is altered; the committed markdown stays the working copy.
Outputs: papers/paperF_render.html (figs/ relative paths resolve there)
+ colleague_package/hajek2026d_hubble_meter.pdf via Edge headless.
"""
import os, re, subprocess, sys
import markdown

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

src = open('papers/paperF_hubble_meter.md', encoding='utf-8').read()

# 1. strip the draft-status banner (the first ** ... ** block after the title)
src = re.sub(r"\n\*\*Draft 0\.3[^*]*\*\*\n", "\n", src, count=1, flags=re.S)
# 2. internal tags
src = re.sub(r"\s*\[ASSEMBLY[^\]]*\]", "", src)
src = re.sub(r"\s*\[All quotes verified[^\]]*\]", "", src)
# 3. reference-list verification tags: "(v)" with optional bracket notes
src = re.sub(r"\s*\(v\)(\s*\[[^\]]*\])?", "", src)
src = re.sub(r"\s*\[verify[^\]]*\]", "", src)
src = re.sub(r"\s*\[ID corrected[^\]]*\]", "", src)
# 4. the reference-list ASSEMBLY preamble line
src = re.sub(r"\n\[ASSEMBLY: verify every entry[^\]]*\]\n", "\n", src)
# 5. version line under the author
src = src.replace(
    "Author: Filip Hájek (independent researcher).",
    "Author: Filip Hájek (independent researcher).\n\n"
    "*Preprint draft, 2026-09-23. Full analysis chain: "
    "github.com/TheCake/thermal-horizon-rar; archived record: "
    "doi.org/10.5281/zenodo.22050990 (concept).*")

assert '[ASSEMBLY' not in src and '(v)' not in src, "internal tags remain"

body = markdown.markdown(src, extensions=['tables'])

CSS = """
@page { size: A4; margin: 20mm 18mm 20mm 18mm; }
html { -webkit-print-color-adjust: exact; }
body { font-family: "Charter", "Georgia", "Cambria", serif;
  font-size: 10.5pt; line-height: 1.5; color: #111; max-width: 46em;
  margin: 0 auto; padding: 1em; }
h1 { font-size: 17pt; line-height: 1.25; margin: 0 0 0.6em 0; font-weight: 600; }
h2 { font-size: 13pt; margin: 1.6em 0 0.5em 0; font-weight: 600;
     border-bottom: 1px solid #ccc; padding-bottom: 0.15em; page-break-after: avoid; }
h3 { font-size: 11.5pt; margin: 1.2em 0 0.4em 0; font-weight: 600; page-break-after: avoid; }
p { margin: 0 0 0.7em 0; text-align: justify; hyphens: auto; }
em { color: #333; }
hr { border: none; border-top: 1px solid #ccc; margin: 1.4em 0; }
ul, ol { margin: 0 0 0.8em 0; padding-left: 1.6em; }
li { margin-bottom: 0.25em; }
code { font-family: "Consolas", monospace; font-size: 0.88em; background: #f4f4f4;
       padding: 0.08em 0.25em; border-radius: 2px; }
table { border-collapse: collapse; width: 100%; margin: 0.6em 0 1em 0;
        font-size: 9pt; page-break-inside: avoid; }
th, td { border: 1px solid #bbb; padding: 0.32em 0.5em; text-align: left;
         vertical-align: top; }
th { background: #eee; font-weight: 600; }
img { max-width: 100%; margin: 0.6em 0 0.2em 0; page-break-inside: avoid; }
"""

html = ("<!doctype html>\n<html lang='en'><head><meta charset='utf-8'>"
        "<title>Hajek 2026d - A Hubble meter from galaxy rotation curves"
        "</title><style>" + CSS + "</style></head><body>\n"
        + body + "\n</body></html>\n")
out_html = os.path.abspath('papers/paperF_render.html')
with open(out_html, 'w', encoding='utf-8') as f:
    f.write(html)
print("wrote", out_html, len(html), "bytes")

edges = [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
         r"C:\Program Files\Microsoft\Edge\Application\msedge.exe"]
edge = next((e for e in edges if os.path.exists(e)), None)
assert edge, "Edge not found"
out_pdf = os.path.abspath('colleague_package/hajek2026d_hubble_meter.pdf')
if os.path.exists(out_pdf):
    os.remove(out_pdf)
subprocess.run([edge, "--headless", "--disable-gpu",
                f"--print-to-pdf={out_pdf}", "--no-pdf-header-footer",
                "file:///" + out_html.replace('\\', '/')],
               check=True, timeout=120)
sz = os.path.getsize(out_pdf)
print("wrote", out_pdf, sz, "bytes")
assert sz > 200000, "PDF suspiciously small"

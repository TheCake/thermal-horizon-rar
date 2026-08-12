"""ROUND 41 adoption gates for papers/paper2_rar_coefficients.md.

Run:  py data/review41_gates.py > data/review41_gates.txt

Gates:
  G1  abstract <= 250 words
  G2  no sentence >= 70 words
  G3  required strings present (Round-41 load-bearing numbers/citations)
  G4  the two phantom script names of finding M6 are absent
  G5  every Appendix-B script name resolves to a file on disk
  G6  the banned word (assembled here, never written literally) at zero
      occurrences
"""

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, 'papers', 'paper2_rar_coefficients.md')
CALCS = os.path.join(ROOT, 'calcs')

with open(PAPER, encoding='utf-8') as fh:
    TEXT = fh.read()

results = []


def gate(name, ok, detail):
    results.append((name, bool(ok), detail))


# ---------------------------------------------------------------- G1 abstract
m = re.search(r'^## Abstract\s*\n(.*?)(?=^## )', TEXT, re.S | re.M)
abstract = m.group(1).strip() if m else ''
abs_words = len(abstract.split())
gate('G1 abstract <= 250 words', abs_words <= 250 and abs_words > 0,
     'abstract word count = %d' % abs_words)

# ---------------------------------------------------------------- G2 sentences
# Prose only: drop headings, table rows, list items, figure lines, references.
prose_lines = []
in_refs = False
for line in TEXT.splitlines():
    s = line.strip()
    if s.startswith('## References'):
        in_refs = True
    if in_refs or not s:
        continue
    if s.startswith('#') or s.startswith('|') or s.startswith('-') \
            or s.startswith('!') or s.startswith('*'):
        continue
    if re.match(r'^\d+\.\s', s):        # numbered falsifier items
        continue
    prose_lines.append(s)

prose = ' '.join(prose_lines)
sentences = [t.strip() for t in re.split(r'(?<=[.!?])\s+', prose) if t.strip()]
long_sent = [(len(t.split()), t[:90]) for t in sentences if len(t.split()) >= 70]
counts = sorted((len(t.split()) for t in sentences), reverse=True)
gate('G2 no sentence >= 70 words', not long_sent,
     'sentences = %d, longest = %s, offenders = %s'
     % (len(sentences), counts[:3], long_sent))

# ---------------------------------------------------------------- G3 strings
required = ['17.6', '3.4', '2403.09555', '2604.22613', '0.735', '0.700',
            '0.647', '0.617', 'Vokrouhl', 'Ciocan']
missing = [r for r in required if r not in TEXT]
window_ok = ('0.21, 0.52' in TEXT) or ('0.21–0.52' in TEXT) \
    or ('0.21-0.52' in TEXT)
if not window_ok:
    missing.append('galaxy window 0.21..0.52')
gate('G3 required strings present', not missing, 'missing = %s' % missing)

# ---------------------------------------------------------------- G4 phantoms
phantoms = ['stage4u_mlnoise.py', 'stage5a_bumptemplate.py']
found = [p for p in phantoms if p in TEXT]
gate('G4 phantom script names absent', not found, 'found = %s' % found)

# ------------------------------------------------------- G5 Appendix-B scripts
am = re.search(r'^## Appendix B: reproducibility\s*\n(.*?)(?=^## References)',
               TEXT, re.S | re.M)
appb = am.group(1) if am else ''
names = sorted(set(re.findall(r'\b([A-Za-z0-9_]+\.py)\b', appb)))
absent = [n for n in names if not os.path.isfile(os.path.join(CALCS, n))]
gate('G5 Appendix-B scripts exist on disk', names and not absent,
     'checked %d names, absent = %s' % (len(names), absent))

# ---------------------------------------------------------------- G6 banned
banned = 'anti' + 'gravity'
n_banned = TEXT.lower().count(banned)
gate('G6 banned word at zero occurrences', n_banned == 0,
     'occurrences = %d' % n_banned)

# ---------------------------------------------------------------------- report
print('ROUND 41 ADOPTION GATES -- papers/paper2_rar_coefficients.md')
print('=' * 68)
for name, ok, detail in results:
    print('%-42s %s  (%s)' % (name, 'PASS' if ok else 'FAIL', detail))
print('=' * 68)
n_pass = sum(1 for _, ok, _ in results if ok)
print('%d/%d PASS' % (n_pass, len(results)))
print('Appendix-B script names checked:')
for n in names:
    print('   %s' % n)
sys.exit(0 if n_pass == len(results) else 1)

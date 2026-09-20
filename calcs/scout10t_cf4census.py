# SCOUT for stage 10T (desk census, instrument-side, 2026-09-20 late):
# how many of SPARC's 97 flow galaxies now hold anchor-grade distances
# in CosmicFlows-4 (Tully+23, J/ApJ/944/94; fetched data/cf4/ with
# SHA256 manifest), and how many LITTLE THINGS dwarfs (Oh+15, already
# on disk from the 9R arc) genuinely extend the anchored leg?
# PGC-keyed via the SIMBAD resolver (canonical PGC alias = "LEDA n" --
# the trap that zeroed the first run). SCOUT-GRADE: the 21 unresolved
# F5xx/D5xx LSB names need alias work at stage time; findings feed the
# 10T pre-registration, they are NOT sky reads (no fit anywhere here).
# Result summary in data/scout10t_cf4census.txt; name->PGC cache in
# data/scout10t_pgc_cache.json.
import json, os, re, time, urllib.parse, urllib.request

os.chdir(r'c:\Users\hfili\Projects\PhysicsResearch')
CACHE = r'data\scout10t_pgc_cache.json'
SIM = 'https://simbad.cds.unistra.fr/simbad/sim-id?'

def pgc_of(name):
    # SIMBAD canonical PGC alias is "LEDA NNNNN" (PGC == LEDA number)
    url = SIM + urllib.parse.urlencode(
        {'output.format': 'ASCII', 'Ident': name})
    try:
        req = urllib.request.Request(
            url, headers={'User-Agent': 'physres-xmatch/1.0'})
        txt = urllib.request.urlopen(req, timeout=30).read().decode(
            'utf-8', 'replace')
        m = re.search(r'(?:LEDA|PGC)\s+(\d+)', txt)
        return int(m.group(1)) if m else None
    except Exception:
        return None

def sparc_query_names(n):
    outs = []
    m = re.match(r'^([A-Za-z]+?)(\d.*)$', n)
    if m:
        pre, num = m.groups()
        if pre in ('UGC', 'NGC', 'IC', 'DDO', 'UGCA', 'PGC'):
            outs.append(f'{pre} {num.lstrip("0")}')
        elif pre == 'ESO':
            outs.append('ESO ' + num.replace('-G', '-'))
            outs.append(n)
        elif pre == 'F':
            outs.append('LSBC F' + num)
    outs.append(n)
    outs.append(re.sub(r'(\D)(\d)', r'\1 \2', n, count=1))
    return outs

# SPARC names + f_D + D
sparc = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        sparc[t[0]] = (int(t[4]), float(t[2]), float(t[3]))
    except ValueError:
        continue
print(f'SPARC rows: {len(sparc)}')

cache = {}
if os.path.exists(CACHE):
    cache = json.load(open(CACHE))
for i, n in enumerate(sparc):
    if n in cache: continue
    p = None
    for qn in sparc_query_names(n):
        p = pgc_of(qn)
        if p: break
        time.sleep(0.1)
    cache[n] = p
    time.sleep(0.1)
    if (i+1) % 40 == 0:
        print(f'  resolved {i+1}/{len(sparc)}')
        json.dump(cache, open(CACHE, 'w'))
json.dump(cache, open(CACHE, 'w'))
unres = [n for n, p in cache.items() if not p]
print(f'PGC resolved: {len(sparc)-len(unres)}/{len(sparc)}; '
      f'unresolved: {unres}')

# CF4 table2 parse (fixed-width per ReadMe; verified on first rows)
cf4 = {}
with open('data/cf4/table2.dat') as f:
    for l in f:
        try:
            pgc = int(l[0:7])
        except ValueError:
            continue
        def fld(a, b):
            s = l[a:b].strip()
            return float(s) if s else None
        cf4[pgc] = dict(
            DM=fld(28, 34), eDM=fld(35, 40),
            snia=fld(41, 47), tf=fld(53, 59), fp=fld(65, 71),
            sbf=fld(77, 83), snII=fld(90, 96), trgb=fld(102, 107),
            ceph=fld(113, 119), mas=fld(126, 131))
print(f'CF4 galaxies: {len(cf4)}')
spot = cf4.get(2557)   # NGC 253-adjacent sanity? just print one known row
print('spot PGC 2557:', spot)

ANCH = ('trgb', 'ceph', 'mas', 'sbf', 'snia', 'snII')
res = {m: [] for m in ANCH}
flow_hit, flow_miss, flow_nopgc = [], [], []
for n, (fd, D, eD) in sparc.items():
    if fd != 1: continue
    p = cache.get(n)
    if not p:
        flow_nopgc.append(n); continue
    row = cf4.get(p)
    if not row:
        flow_miss.append(n); continue
    hits = [m for m in ANCH if row[m] is not None]
    if hits:
        flow_hit.append((n, p, D, hits,
                         {m: row[m] for m in hits}))
        for m in hits: res[m].append(n)
    else:
        flow_miss.append(n)
print()
print(f'FLOW GALAXIES (97): anchor-grade CF4 hit = {len(flow_hit)}; '
      f'no anchor method = {len(flow_miss)}; unresolved PGC = '
      f'{len(flow_nopgc)} {flow_nopgc}')
for m in ANCH:
    print(f'  {m:5s}: {len(res[m]):3d}  {sorted(res[m])[:12]}'
          f'{" ..." if len(res[m]) > 12 else ""}')
print()
print('per-galaxy anchor hits (name, PGC, SPARC flow D, methods, DMs):')
for n, p, D, hits, dms in sorted(flow_hit):
    dstr = ', '.join(f'{m}={10**((v-25)/5):.1f}Mpc' for m, v in dms.items())
    print(f'  {n:11s} PGC {p:7d}  D_flow={D:6.1f}  {dstr}')

# also: UMa + anchored upgrades (better e_D or new methods)
uma_up, anc_up = [], []
for n, (fd, D, eD) in sparc.items():
    p = cache.get(n)
    row = cf4.get(p) if p else None
    if not row: continue
    hits = [m for m in ('trgb', 'ceph', 'mas') if row[m] is not None]
    if fd == 4 and hits: uma_up.append((n, hits))
    if fd in (2, 3, 5) and 'mas' in hits: anc_up.append((n, hits))
print()
print(f'UMa galaxies with direct TRGB/Ceph/maser in CF4: {len(uma_up)} '
      f'{uma_up}')
print(f'anchored galaxies gaining a MASER distance: {anc_up}')

# ---- LITTLE THINGS overlap vs SPARC (PGC-keyed, fixed resolver) ----
lt = []
for l in open('data/littlethings/table1.dat'):
    nm = l.split('|')[0].strip()
    if nm: lt.append(nm)
lt_query = {n: n.replace('_', ' ') for n in lt}
lt_query['CVnIdwA'] = 'CVn I dwA'
lt_query['F564-V3'] = 'LSBC F564-V03'
sp_pgcs = {p for p in cache.values() if p}
print()
print('LITTLE THINGS overlap (PGC-keyed):')
ov, new, un = [], [], []
for n in lt:
    key = 'LT:' + n
    if key not in cache:
        cache[key] = pgc_of(lt_query[n])
        time.sleep(0.1)
    p = cache[key]
    if not p: un.append(n)
    elif p in sp_pgcs: ov.append((n, p))
    else: new.append((n, p))
json.dump(cache, open(CACHE, 'w'))
print(f'  overlap ({len(ov)}): {ov}')
print(f'  genuinely new ({len(new)}): {[n for n, _ in new]}')
print(f'  unresolved ({len(un)}): {un}')
# CF4 anchor methods for the genuinely-new LT dwarfs
print('  CF4 anchor methods for the new LT dwarfs:')
for n, p in new:
    row = cf4.get(p)
    if row:
        hits = {m: f'{10**((row[m]-25)/5):.1f}Mpc'
                for m in ANCH if row[m] is not None}
        print(f'    {n:10s} PGC {p:7d}: {hits if hits else "no anchor method in CF4"}')
    else:
        print(f'    {n:10s} PGC {p:7d}: not in CF4 table2')

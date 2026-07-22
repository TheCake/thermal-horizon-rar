"""
GAIA WIDE-BINARY APSIDAL-ORIENTATION TEST (first pass).
Data: El-Badry & Rix 2018 catalog (Gaia DR2, d<200 pc, 55k pairs), VizieR TAP.
Statistic: axial (180-degree-degenerate) orientation of the on-sky separation
vector relative to the PROJECTED galactic field direction (toward the galactic
center), weighted by the projected field magnitude. Rayleigh test on 2*psi.
Newtonian-regime bins (s < ~2 kAU) serve as the selection-systematics control:
any instrumental anisotropy (scanning law, crowding) should appear in ALL bins;
new physics only in the wide ones.
Secondary statistic: angle between relative proper motion and separation.
"""
import csv, math

D2R = math.pi/180
# galactic center, ICRS
A_GC, D_GC = 266.4051*D2R, -28.936175*D2R
GC = (math.cos(D_GC)*math.cos(A_GC), math.cos(D_GC)*math.sin(A_GC), math.sin(D_GC))

def tangent_frame(a, d):
    east  = (-math.sin(a), math.cos(a), 0.0)
    north = (-math.sin(d)*math.cos(a), -math.sin(d)*math.sin(a), math.cos(d))
    return east, north

def dot(u, v): return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

bins = [(100, 1000), (1000, 3000), (3000, 8000), (8000, 50000)]
stats = {b: {'C': 0.0, 'S': 0.0, 'W': 0.0, 'W2': 0.0, 'N': 0,
             'gC': 0.0, 'gS': 0.0} for b in bins}

n_total = n_used = 0
with open('data/widebinaries_elbadry_rix2018.csv', newline='') as f:
    for row in csv.DictReader(f):
        n_total += 1
        try:
            a1, d1 = float(row['RA_ICRS'])*D2R, float(row['DE_ICRS'])*D2R
            a2, d2 = float(row['RA2deg'])*D2R, float(row['DE2deg'])*D2R
            s_au   = float(row['sAU'])
            pm1 = (float(row['pmRA']), float(row['pmDE']))
            pm2 = (float(row['pmRA2']), float(row['pmDE2']))
        except (ValueError, KeyError):
            continue
        b = next((bb for bb in bins if bb[0] <= s_au < bb[1]), None)
        if b is None: continue
        n_used += 1
        # separation position angle (east-of-north), small-angle tangent plane
        dx = (a2-a1)*math.cos(0.5*(d1+d2))     # east component
        dy = (d2-d1)                            # north component
        pa_sep = math.atan2(dx, dy)
        # projected galactic-field direction & weight at star 1
        east, north = tangent_frame(a1, d1)
        gE, gN = dot(GC, east), dot(GC, north)
        w = math.hypot(gE, gN)                  # in-plane field fraction (0..1)
        pa_ref = math.atan2(gE, gN)
        th = 2.0*(pa_sep - pa_ref)              # axial statistic
        st = stats[b]
        st['C'] += w*math.cos(th); st['S'] += w*math.sin(th)
        st['W'] += w; st['W2'] += w*w; st['N'] += 1
        # secondary: relative-PM vs separation angle (0-90 deg)
        dpm = (pm2[0]-pm1[0], pm2[1]-pm1[1])
        npm = math.hypot(*dpm); nsep = math.hypot(dx, dy)
        if npm > 0 and nsep > 0:
            cg = abs((dpm[0]*dx + dpm[1]*dy)/(npm*nsep))
            st['gC'] += min(cg, 1.0); st['gS'] += 1

print(f"pairs read: {n_total}, used in bins: {n_used}\n")
print(f"{'separation':>16} {'N':>6} {'R (axial)':>10} {'p(Rayleigh)':>12} "
      f"{'mean axis':>10} {'<v-r angle>':>12}")
for b in bins:
    st = stats[b]
    if st['W'] <= 0: continue
    R = math.hypot(st['C'], st['S'])/st['W']
    Neff = st['W']**2/st['W2'] if st['W2'] > 0 else 0
    p = math.exp(-Neff*R*R) if Neff > 0 else 1.0
    mean_axis = 0.5*math.degrees(math.atan2(st['S'], st['C']))
    gbar = math.degrees(math.acos(min(st['gC']/max(st['gS'],1),1.0)))
    print(f"{b[0]:>7,}-{b[1]:<8,} {st['N']:>6} {R:>10.4f} {p:>12.4g} "
          f"{mean_axis:>9.1f}deg {gbar:>11.1f}deg")
print("""
Reading guide:
  R ~ 0, p ~ 1        -> orientations uniform (no signal)
  p small in ALL bins -> selection systematic (scanning law etc.), not physics
  p small ONLY in wide bins, mean axis ~0 or ~90 deg to the field
                      -> the anomaly's signature (libration axes lock to field)
  <v-r angle>: mean angle between relative velocity and separation; tracks the
  eccentricity distribution -- MOND-regime orbits predicted rounder (lower e)
  by some completions -> higher mean angle in wide bins.""")

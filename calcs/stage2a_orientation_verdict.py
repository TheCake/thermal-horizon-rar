"""
STAGE 2A: does our Gaia orientation null (R < 2.2% @ 2sigma, 8-50 kAU) bite?

Chain: measured force law (Stage 1: p=0.443, a0=1.03e-10) -> precession law
A(a,e) + B(a,e) sin(2 psi) per orbit -> stationary distribution of the
apsis-vs-field angle psi under field rotation (230 Myr) -> Monte Carlo of the
OBSERVABLE sky statistic exactly as measured -> predicted R vs data bound.
All numpy, CPU, seconds.
"""
import numpy as np

A0   = 1.03e-10 * 6657.0          # Stage-1 a0 -> AU/yr^2  (1 m/s^2 = 6657 AU/yr^2 /... )
# careful: 1.2e-10 m/s^2 = 7.99e-7 AU/yr^2  =>  1 m/s^2 = 6.658e3 AU/yr^2
A0   = 1.03e-10 * 6.658e3         # = 6.86e-7 AU/yr^2
P_SCREEN = 0.443
GM   = 4*np.pi**2                  # per Msun
GEXT = 1.8 * 1.2e-10 * 6.658e3    # galactic field in AU/yr^2 (1.8 a0_canonical)
T_FIELD = 230e6                    # yr, galactic rotation period
M_TOT = 1.5                        # Msun typical

def nu(y):
    y = np.clip(y, 1e-12, None)
    return (1 - np.exp(-y**P_SCREEN))**(-1/(2*P_SCREEN))

def accel(pos, mond=True):
    # pos: (N,2) in orbit plane; field along +x (in-plane component)
    r = np.linalg.norm(pos, axis=1, keepdims=True)
    gs = -GM*M_TOT*pos/r**3
    if not mond: return gs
    tot = gs + np.array([GEXT, 0.0])
    y = np.linalg.norm(tot, axis=1, keepdims=True)/A0
    return nu(y)*gs

def drift_batch(a_kau, e, om0_deg, n_orbits=12, spo=4500, M=M_TOT):
    """LSQ apsidal drift [deg/orbit] for a batch of om0 at one (a,e)."""
    a = a_kau*1e3
    P = np.sqrt(a**3/M)
    dt = P/spo
    om0 = np.radians(np.asarray(om0_deg, float))
    nB = len(om0)
    rp, vp = a*(1-e), np.sqrt(GM*M*(1+e)/(a*(1-e)))
    pos = np.stack([rp*np.cos(om0), rp*np.sin(om0)], 1)
    vel = np.stack([-vp*np.sin(om0), vp*np.cos(om0)], 1)
    rprev = np.linalg.norm(pos, axis=1); dec = np.zeros(nB, bool)
    aps = [[] for _ in range(nB)]
    for _ in range(int(n_orbits*spo)):
        k1v, k1p = accel(pos), vel
        k2v = accel(pos+0.5*dt*k1p); k2p = vel+0.5*dt*k1v
        k3v = accel(pos+0.5*dt*k2p); k3p = vel+0.5*dt*k2v
        k4v = accel(pos+dt*k3p);     k4p = vel+dt*k3v
        pos = pos + dt/6*(k1p+2*k2p+2*k3p+k4p)
        vel = vel + dt/6*(k1v+2*k2v+2*k3v+k4v)
        r = np.linalg.norm(pos, axis=1)
        peri = dec & (r > rprev)
        if peri.any():
            ang = np.degrees(np.arctan2(pos[:,1], pos[:,0]))
            for i in np.where(peri)[0]:
                if aps[i]:                       # unwrap vs previous
                    d = ang[i]-aps[i][-1]
                    while d > 180: d -= 360
                    while d < -180: d += 360
                    aps[i].append(aps[i][-1]+d)
                else:
                    aps[i].append(ang[i])
        dec = r < rprev; rprev = r
    out = np.full(nB, np.nan)
    for i in range(nB):
        y = np.array(aps[i])
        if len(y) >= 4:
            x = np.arange(len(y))
            out[i] = np.polyfit(x, y, 1)[0]      # robust slope, deg/orbit
    return out

def precession_law(a_kau, e):
    om = [0., 22.5, 45., 67.5, 90.]
    d = drift_batch(a_kau, e, om)
    # fit d(om) = A + B sin(2 om)
    s2 = np.sin(2*np.radians(om))
    ok = ~np.isnan(d)
    X = np.stack([np.ones(ok.sum()), s2[ok]], 1)
    coef, *_ = np.linalg.lstsq(X, d[ok], rcond=None)
    return coef[0], coef[1], np.sqrt((a_kau*1e3)**3/M_TOT)

# validation against the earlier space_model result (~ -3.8 deg/orbit)
_sav = (P_SCREEN, A0, M_TOT)
P_SCREEN, A0, M_TOT = 0.5, 7.99e-7, 2.0
val = drift_batch(10, 0.5, [0.])[0]
P_SCREEN, A0, M_TOT = _sav
print(f"validation (10 kAU, e=0.5, M=2, RAR p=0.5): drift = {val:.2f} deg/orbit "
      f"(space_model got -3.76)")

print("precession law on the (a,e) grid  [deg/orbit]")
print(f"{'a[kAU]':>7} {'e':>5} {'A':>8} {'B':>8} {'Omega_f':>8}")
grid = []
for a_kau in (8, 12, 18, 27, 40):
    for e in (0.3, 0.5, 0.7, 0.85):
        A, B, P = precession_law(a_kau, e)
        Om = 360.0*P/T_FIELD
        grid.append((a_kau, e, A, B, Om))
        print(f"{a_kau:>7} {e:>5} {A:>8.2f} {B:>8.2f} {Om:>8.2f}")

# ---- Monte Carlo of the observable ----
rng = np.random.default_rng(7)
N = 200_000
# separations p(a) ~ a^-1.6 over 8-50 kAU; e thermal p(e)=2e
u = rng.random(N)
lo, hi = 8.0, 50.0; g = -0.6
a_s = (lo**g + u*(hi**g - lo**g))**(1/g)          # inverse CDF of a^-1.6
e_s = 0.95*rng.random(N)**(1/2.25)   # superthermal p(e)~e^1.25 (Hwang+22, s>1kAU)
# interpolate A,B,Omega from grid (nearest in e, linear in log a)
ga = np.array(sorted(set(x[0] for x in grid)), float)
ge = np.array(sorted(set(x[1] for x in grid)))
Agrid = np.zeros((len(ga), len(ge))); Bgrid = np.zeros_like(Agrid); Ogrid = np.zeros_like(Agrid)
for (ak, ek, A, B, Om) in grid:
    Agrid[list(ga).index(ak), list(ge).index(ek)] = A
    Bgrid[list(ga).index(ak), list(ge).index(ek)] = B
    Ogrid[list(ga).index(ak), list(ge).index(ek)] = Om
ia = np.clip(np.searchsorted(ga, a_s)-1, 0, len(ga)-2)
fa = np.clip((np.log(a_s)-np.log(ga[ia]))/(np.log(ga[ia+1])-np.log(ga[ia])), 0, 1)
ie = np.clip(np.searchsorted(ge, e_s)-1, 0, len(ge)-2)
fe = np.clip((e_s-ge[ie])/(ge[ie+1]-ge[ie]), 0, 1)
def interp(G):
    return (G[ia,ie]*(1-fa)*(1-fe) + G[ia+1,ie]*fa*(1-fe)
          + G[ia,ie+1]*(1-fa)*fe + G[ia+1,ie+1]*fa*fe)
A_s, B_s, O_s = interp(Agrid), interp(Bgrid), interp(Ogrid)
# random orbit plane: in-plane field fraction scales the torque term
nrm = rng.normal(size=(N,3)); nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
f_ip = np.sqrt(np.clip(1-nrm[:,0]**2, 0, 1))       # field along x-hat
B_eff = B_s*f_ip
# evolve psi (apsis angle relative to rotating field) to stationarity
psi = rng.random(N)*np.pi
v = np.radians(A_s - O_s)                            # rad/orbit
b = np.radians(B_eff)
for _ in range(3000):
    psi = psi + (v + b*np.sin(2*psi))
psi %= np.pi
# orbital phase, time-uniform: solve Kepler
Mean = rng.random(N)*2*np.pi
E = Mean.copy()
for _ in range(8):
    E -= (E - e_s*np.sin(E) - Mean)/(1 - e_s*np.cos(E))
theta = 2*np.arctan2(np.sqrt(1+e_s)*np.sin(E/2), np.sqrt(1-e_s)*np.cos(E/2))
r_s = a_s*(1 - e_s*np.cos(E))                        # kAU
# 3D geometry: orbit basis (e_f = in-plane field dir, e2 = n x e_f)
xhat = np.zeros((N,3)); xhat[:,0] = 1
ef = xhat - nrm*nrm[:,[0]]
ef /= np.maximum(np.linalg.norm(ef, axis=1, keepdims=True), 1e-12)
e2 = np.cross(nrm, ef)
Phi = psi + theta
sep = r_s[:,None]*(np.cos(Phi)[:,None]*ef + np.sin(Phi)[:,None]*e2)
# line of sight, projections, the measured statistic
los = rng.normal(size=(N,3)); los /= np.linalg.norm(los, axis=1, keepdims=True)
def proj(vv):
    return vv - los*np.sum(vv*los, axis=1, keepdims=True)
sp, fp = proj(sep), proj(xhat)
w = np.linalg.norm(fp, axis=1)
ns = np.linalg.norm(sp, axis=1)
# SIGNED relative angle in the sky plane (fold only by 180 via doubling):
b1 = fp/np.maximum(w[:,None], 1e-12)
b2 = np.cross(los, b1)
theta_rel = np.arctan2(np.sum(sp*b2, axis=1), np.sum(sp*b1, axis=1))
# projected separation cut: keep the 8-50 kAU window like the data bin
keep = (ns > 8.0) & (ns < 50.0)
th2 = 2*theta_rel[keep]; wk = w[keep]
C = np.sum(wk*np.cos(th2)); S = np.sum(wk*np.sin(th2))
R_pred = np.hypot(C, S)/np.sum(wk)
frac_libr = np.mean(np.abs(v) < np.abs(b))
print(f"\nMC binaries kept after projection cut: {keep.sum()}")
print(f"fraction in LIBRATION (|A-Omega| < |B_eff|): {frac_libr:.2%}")
print(f"PREDICTED observable anisotropy R = {R_pred:.4f}")
print(f"our Gaia measurement (8-50 kAU):   R = 0.0084, N=5917 "
      f"-> 2-sigma bound R < {np.sqrt(3/5917):.4f}")
if R_pred > np.sqrt(3/5917):
    print("VERDICT: the null BITES -- this EFE+rotation scenario is excluded")
else:
    print(f"VERDICT: prediction is BELOW current sensitivity "
          f"(need N ~ {int(3/R_pred**2):,} pairs -> EDR3-scale catalog)")

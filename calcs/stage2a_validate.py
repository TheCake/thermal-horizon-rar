"""Isolated validation of the drift integrator (consistent parameters)."""
import numpy as np

GM = 4*np.pi**2

def make_accel(p_screen, a0, gext, M):
    def nu(y):
        y = np.clip(y, 1e-12, None)
        return (1 - np.exp(-y**p_screen))**(-1/(2*p_screen))
    def accel(pos):
        r = np.linalg.norm(pos, axis=1, keepdims=True)
        gs = -GM*M*pos/r**3
        tot = gs + np.array([gext, 0.0])
        return nu(np.linalg.norm(tot, axis=1, keepdims=True)/a0)*gs
    return accel

def drift(a_kau, e, om0_deg, accel, M, n_orbits=12, spo=9000):
    a = a_kau*1e3
    P = np.sqrt(a**3/M); dt = P/spo
    om0 = np.radians(np.asarray(om0_deg, float)); nB = len(om0)
    rp, vp = a*(1-e), np.sqrt(GM*M*(1+e)/(a*(1-e)))
    pos = np.stack([rp*np.cos(om0), rp*np.sin(om0)], 1)
    vel = np.stack([-vp*np.sin(om0), vp*np.cos(om0)], 1)
    vr_prev = np.sum(pos*vel, axis=1)
    aps = [[] for _ in range(nB)]; rmin = np.full(nB, 1e18); rmax = np.zeros(nB)
    ecc_track = []
    for step in range(int(n_orbits*spo)):
        k1v, k1p = accel(pos), vel
        k2v = accel(pos+0.5*dt*k1p); k2p = vel+0.5*dt*k1v
        k3v = accel(pos+0.5*dt*k2p); k3p = vel+0.5*dt*k2v
        k4v = accel(pos+dt*k3p);     k4p = vel+dt*k3v
        pos = pos + dt/6*(k1p+2*k2p+2*k3p+k4p)
        vel = vel + dt/6*(k1v+2*k2v+2*k3v+k4v)
        r = np.linalg.norm(pos, axis=1)
        rmin = np.minimum(rmin, r); rmax = np.maximum(rmax, r)
        vr = np.sum(pos*vel, axis=1)
        peri = (vr_prev < 0) & (vr >= 0)          # radial velocity sign change
        if peri.any():
            ang = np.degrees(np.arctan2(pos[:,1], pos[:,0]))
            for i in np.where(peri)[0]:
                if aps[i]:
                    d = ang[i]-aps[i][-1]
                    while d > 180: d -= 360
                    while d < -180: d += 360
                    aps[i].append(aps[i][-1]+d)
                else:
                    aps[i].append(ang[i])
        vr_prev = vr
    out = []
    for i in range(nB):
        y = np.array(aps[i])
        out.append(np.polyfit(np.arange(len(y)), y, 1)[0] if len(y) >= 4 else np.nan)
    return np.array(out), rmin, rmax

# Case: 10 kAU, e=0.5, M=2, p=0.5 (RAR), canonical a0, gext=1.8 a0
A0, GEXT = 7.99e-7, 1.8*7.99e-7
acc = make_accel(0.5, A0, GEXT, 2.0)
d, rmin, rmax = drift(10, 0.5, [0., 45., 90.], acc, 2.0)
print("modified drift [deg/orbit] at om0=0/45/90:", np.round(d, 2))
print("r range [kAU]:", np.round(rmin/1e3, 2), "-", np.round(rmax/1e3, 2))
print("(initial rp-ra would be 5.0 - 15.0; big deviation = orbit reshaped)")

# Newtonian control (same integrator, nu=1)
accN = lambda pos: -GM*2.0*pos/np.linalg.norm(pos,axis=1,keepdims=True)**3
dN, rminN, rmaxN = drift(10, 0.5, [0.], accN, 2.0)
print("Newtonian control drift:", np.round(dN, 4), " (numerical noise floor)")

# convergence check: same modified case at double resolution
d2, _, _ = drift(10, 0.5, [0., 45., 90.], acc, 2.0, spo=18000)
print("modified drift at 2x resolution:", np.round(d2, 2), " (should match)")

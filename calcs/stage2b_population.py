"""
STAGE 2B: GPU population synthesis of wide binaries under the measured law.
One CUDA thread = one binary, integrated start-to-finish in registers (fp64).
Physics: g = nu(|g_int + g_ext(t)|/a0) * g_int, nu = (1-exp(-y^p))^(-1/2p),
p = 0.443, a0 = 1.03e-10 m/s^2 (Stage 1); in-plane external field of magnitude
gext*f_ip rotating with the galactic period. Snapshot after burn-in at random
phase -> orientation statistic + v-tilde, vs a Newtonian twin population.
"""
import time
import numpy as np
import cupy as cp

KERNEL = r'''
extern "C" __global__
void evolve(const double* a_, const double* e_, const double* psi0_,
            const double* fip_, const double* Mtot_, const double* uphase_,
            const double* vp_, const double* tab_,
            double a0, double pscr, double gext, double om_rot,
            double lny0, double dlny, int ntab,
            int n_orbits, int spo, int newton, int N,
            double* out)   // per binary: x,y,vx,vy,phi_field,rmin,rmax
{
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    if (i >= N) return;
    const double GM = 4.0*9.869604401089358;      // 4 pi^2
    double a = a_[i], e = e_[i], psi = psi0_[i], fip = fip_[i], M = Mtot_[i];
    double P = sqrt(a*a*a/M);
    double dt = P/spo;
    double rp = a*(1.0-e);
    double vp = (vp_[i] > 0.0) ? vp_[i]
              : sqrt(GM*M*(1.0+e)/(a*(1.0-e)));
    double x = rp*cos(psi),  y = rp*sin(psi);
    double vx = -vp*sin(psi), vy = vp*cos(psi);
    double t = 0.0, rmin = 1e30, rmax = 0.0;
    double ge = gext*fip;
    long long nsteps = (long long)n_orbits*spo + (long long)(uphase_[i]*spo);
    for (long long s = 0; s < nsteps; s++) {
        // RK4 on (x,y,vx,vy)
        double kx[4], ky[4], kvx[4], kvy[4];
        double cx = x, cy = y, cvx = vx, cvy = vy;
        double phi = om_rot*t;
        double gex = ge*cos(phi), gey = ge*sin(phi);
        for (int k = 0; k < 4; k++) {
            double r2 = cx*cx + cy*cy;
            double r  = sqrt(r2);
            double gs = GM*M/r2;
            double gsx = -gs*cx/r, gsy = -gs*cy/r;
            double ax, ay;
            if (newton == 1) { ax = gsx; ay = gsy; }
            else if (newton == 2) {
                // Chae-Milgrom 2022 AQUAL fitting formula (radial, conservative)
                double yy = gs/a0;
                double eN = gext/a0;
                double be = 1.1*eN;
                double yb = sqrt(yy*yy + be*be);
                double sq = sqrt(0.25 + 1.0/yb);
                double nus = 0.5 + sq;
                double nuhat = (1.0/yb)/(2.0*nus*sq);
                double bst = nus*(1.0 + tanh(pow(be/yy, 1.2))*nuhat/3.0);
                ax = bst*gsx; ay = bst*gsy;
            }
            else if (newton == 3) {
                // C&M EFE structure, but with OUR measured screening family
                double yy = gs/a0;
                double eN = gext/a0;
                double be = 1.1*eN;
                double yb = sqrt(yy*yy + be*be);
                double tt = pow(yb, pscr);
                double et = exp(-tt);
                double nus = pow(1.0 - et, -0.5/pscr);
                double nuhat = tt*et/(2.0*(1.0 - et));
                double bst = nus*(1.0 + tanh(pow(be/yy, 1.2))*nuhat/3.0);
                ax = bst*gsx; ay = bst*gsy;
            }
            else if (newton == 5) {
                // radial boost from a numerically solved EFE table (ln y grid)
                double yy = gs/a0;
                double u = (log(yy) - lny0)/dlny;
                double bst;
                if (u <= 0.0) bst = tab_[0];
                else if (u >= ntab-1) bst = 1.0;
                else {
                    int i0 = (int)u;
                    double fr = u - i0;
                    bst = tab_[i0]*(1.0-fr) + tab_[i0+1]*fr;
                }
                ax = bst*gsx; ay = bst*gsy;
            }
            else if (newton == 4) {
                // Bose-Einstein reading: g = gN * (1 + n_BE(x)),
                // x = sqrt(|g_tot|/a0), external field inside the mode energy
                double yy = sqrt(gs*gs + ge*ge)/a0;
                double xx = sqrt(yy);
                double nocc = (xx > 40.0) ? 0.0 : 1.0/(exp(xx) - 1.0);
                double bst = 1.0 + nocc;
                ax = bst*gsx; ay = bst*gsy;
            }
            else if (newton == 0) {
                // CONSERVATIVE radial recipe: nu on quadrature combination,
                // force stays central -> energy-conserving by construction
                double yy = sqrt(gs*gs + ge*ge)/a0;
                if (yy < 1e-12) yy = 1e-12;
                double nu = pow(1.0 - exp(-pow(yy, pscr)), -0.5/pscr);
                ax = nu*gsx; ay = nu*gsy;
            }
            else {
                // legacy anisotropic recipe (newton == -1): NON-conservative,
                // valid only for few-orbit runs
                double tx = gsx + gex, ty = gsy + gey;
                double yy = sqrt(tx*tx + ty*ty)/a0;
                if (yy < 1e-12) yy = 1e-12;
                double nu = pow(1.0 - exp(-pow(yy, pscr)), -0.5/pscr);
                ax = nu*gsx; ay = nu*gsy;
            }
            kx[k] = cvx; ky[k] = cvy; kvx[k] = ax; kvy[k] = ay;
            double h = (k == 0 || k == 1) ? 0.5*dt : dt;
            if (k < 3) {
                cx = x + h*kx[k];  cy = y + h*ky[k];
                cvx = vx + h*kvx[k]; cvy = vy + h*kvy[k];
            }
        }
        x  += dt/6.0*(kx[0] + 2.0*kx[1] + 2.0*kx[2] + kx[3]);
        y  += dt/6.0*(ky[0] + 2.0*ky[1] + 2.0*ky[2] + ky[3]);
        vx += dt/6.0*(kvx[0] + 2.0*kvx[1] + 2.0*kvx[2] + kvx[3]);
        vy += dt/6.0*(kvy[0] + 2.0*kvy[1] + 2.0*kvy[2] + kvy[3]);
        t  += dt;
        double r = sqrt(x*x + y*y);
        if (r < rmin) rmin = r;
        if (r > rmax) rmax = r;
    }
    out[7*i+0] = x;  out[7*i+1] = y;
    out[7*i+2] = vx; out[7*i+3] = vy;
    out[7*i+4] = om_rot*t;
    out[7*i+5] = rmin; out[7*i+6] = rmax;
}
'''
mod = cp.RawModule(code=KERNEL, options=('--std=c++11',))
evolve = mod.get_function('evolve')

A0_CAN = 7.99e-7                    # canonical a0 in AU/yr^2
A0     = 1.03e-10*6.658e3           # Stage-1 a0
P_SCR  = 0.443
GEXT   = 1.9*A0_CAN                 # Chae-style external field incl. vertical
OM_ROT = 2*np.pi/230e6              # rad/yr

def run(a, e, psi0, fip, M, uph, n_orbits, spo, newton, a0=A0, pscr=P_SCR,
        gext=GEXT, om_rot=OM_ROT, vp=None, tab=None, lny0=0.0, dlny=1.0):
    N = len(a)
    if vp is None: vp = np.zeros(N)
    if tab is None: tab = np.ones(2)
    ntab = len(tab)
    dev = [cp.asarray(np.ascontiguousarray(v, dtype=np.float64))
           for v in (a, e, psi0, fip, M, uph, vp, tab)]
    out = cp.zeros(7*N, dtype=cp.float64)
    tpb = 128
    evolve(((N+tpb-1)//tpb,), (tpb,),
           (*dev, np.float64(a0), np.float64(pscr), np.float64(gext),
            np.float64(om_rot), np.float64(lny0), np.float64(dlny),
            np.int32(ntab), np.int32(n_orbits), np.int32(spo),
            np.int32(newton), np.int32(N), out))
    cp.cuda.Stream.null.synchronize()
    return cp.asnumpy(out).reshape(N, 7)

# ---------- validation vs CPU integrator ----------
print("validation vs CPU (10 kAU, e=0.5, M=2, p=0.5, a0 canonical, static field)")
t0 = time.time()
o = run(np.array([1e4]*3), np.array([0.5]*3), np.radians([0., 45., 90.]),
        np.array([1.0]*3), np.array([2.0]*3), np.zeros(3),
        12, 9000, 0, a0=A0_CAN, pscr=0.5, gext=1.8*A0_CAN, om_rot=0.0)
print(f"  GPU r-range [kAU]: "
      f"{np.round(o[:,5]/1e3,2)} - {np.round(o[:,6]/1e3,2)}")
print(f"  CPU reference:     [4.23 3.97 3.49] - [11.73 14.39 29.19]")
print(f"  ({time.time()-t0:.1f} s)")

# ---------- production populations ----------
rng = np.random.default_rng(11)
N = 1_000_000
u = rng.random(N); lo, hi, g = 2.0, 50.0, -0.6
a_s   = ((lo**g + u*(hi**g - lo**g))**(1/g))*1e3          # AU
e_s   = 0.95*rng.random(N)**(1/2.25)
psi0  = rng.random(N)*2*np.pi
nrm   = rng.normal(size=(N,3)); nrm /= np.linalg.norm(nrm,axis=1,keepdims=True)
f_ip  = np.sqrt(np.clip(1-nrm[:,0]**2, 0, 1))
M_s   = 0.6 + 1.8*rng.random(N)
uph   = rng.random(N)

print("\nproduction: 1e6 binaries, CONSERVATIVE modified law (10 orbits)...")
t0 = time.time()
om = run(a_s, e_s, psi0, f_ip, M_s, uph, 10, 3000, 0)
t_mod = time.time()-t0
print(f"  modified population done in {t_mod:.0f} s")
t0 = time.time()
on = run(a_s, e_s, psi0, f_ip, M_s, uph, 10, 3000, 1)
print(f"  Newtonian twin done in {(time.time()-t0):.0f} s")
np.save('data/pop_modified.npy', om)
np.save('data/pop_newton.npy', on)
np.save('data/pop_params.npy',
        np.stack([a_s, e_s, psi0, f_ip, M_s, nrm[:,0], nrm[:,1], nrm[:,2]],1))
print("saved to data/pop_*.npy")

# ---------- quick in-plane summaries (full sky projection in analysis pass) ----------
for tag, o in (("modified", om), ("newton", on)):
    xs, ys, vxs, vys, phi = o[:,0], o[:,1], o[:,2], o[:,3], o[:,4]
    r = np.hypot(xs, ys)
    # in-plane separation angle relative to current field direction (axial)
    th = (np.arctan2(ys, xs) - phi)
    R2 = np.abs(np.mean(np.exp(2j*th)))
    # in-plane v-tilde: |v| relative to circular at r
    vt = np.hypot(vxs, vys)/(2*np.pi*np.sqrt(M_s/r))
    wide = a_s > 8e3
    print(f"{tag:>9}: in-plane R2 = {R2:.4f}; "
          f"median v-tilde (a>8kAU) = {np.median(vt[wide]):.3f}; "
          f"e-pump check: median rmax/a = {np.median(o[:,6]/a_s):.2f}")

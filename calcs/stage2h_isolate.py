"""Isolate: 10k identical binaries (a=10 kAU, e=0.5), table mode vs Newton."""
import numpy as np
src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2; A0_CAN = 7.99e-7

tabf = np.load('data/efe_boost_simple.npy')
y_t, b_t = tabf[0][::-1], tabf[1][::-1]
lny = np.log(y_t); lny_u = np.linspace(lny[0], lny[-1], 512)
tab = np.interp(lny_u, lny, b_t)
lny0, dlny = lny_u[0], lny_u[1]-lny_u[0]
def boost_tab(y):
    return np.interp(np.log(np.clip(y,1e-12,None)), lny_u, tab, right=1.0)
def vp_c(a, e, M):
    rp, ra = a*(1-e), a*(1+e)
    xg, wg = np.polynomial.legendre.leggauss(48)
    lo_, hi_ = np.log(rp), np.log(ra)
    lr = 0.5*(hi_-lo_)[:,None]*xg[None,:] + 0.5*(hi_+lo_)[:,None]
    r = np.exp(lr)
    gN = GM*M[:,None]/r**2
    gg = boost_tab(gN/A0_CAN)*gN
    dPhi = np.sum(wg[None,:]*gg*r, axis=1)*0.5*(hi_-lo_)
    return np.sqrt(np.maximum(2*dPhi/(1-(rp/ra)**2), 0))

N = 10_000
rng = np.random.default_rng(1)
a = np.full(N, 1e4); e = np.full(N, 0.5); M = np.full(N, 1.5)
psi0 = rng.random(N)*2*np.pi
fip = np.ones(N); uph = rng.random(N)

om = run(a,e,psi0,fip,M,uph,10,3000,5,a0=A0_CAN,vp=vp_c(a,e,M),
         tab=tab,lny0=lny0,dlny=dlny)
on = run(a,e,psi0,fip,M,uph,10,3000,1)
for tag, o in (("table-mode", om), ("newton", on)):
    r = np.hypot(o[:,0], o[:,1]); v = np.hypot(o[:,2], o[:,3])
    vk = 2*np.pi*np.sqrt(1.5/r)
    print(f"{tag:>10}: median r={np.median(r)/1e3:.2f} kAU, "
          f"median v/vK(r)={np.median(v/vk):.3f}, "
          f"rmax fidelity={np.median(o[:,6]/(a*(1+e))):.3f}")

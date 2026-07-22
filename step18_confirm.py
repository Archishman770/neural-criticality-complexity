import numpy as np, pandas as pd, logging
import mrestimator as mre
logging.getLogger("mrestimator").setLevel(logging.CRITICAL)
import warnings; warnings.filterwarnings("ignore")


TAU_STEPS_GRID = (40, 100, 200)
BIN_STEPS = (1, 2, 4)
KMAX_FIXED = 60
MAX_TAU_FRAC = 0.8
KMAX_STEPS_TIME = 960
N_TRIALS, T, BURN = 12, 60000, 6000
SEED = 20260718
CEIL = np.exp(-1/(MAX_TAU_FRAC*KMAX_FIXED))

def ar1(tau_steps, rng):
    a = np.exp(-1.0/tau_steps)
    x = np.zeros((N_TRIALS, T)); noise = rng.normal(0,1,(N_TRIALS,T))
    for t in range(1,T): x[:,t] = a*x[:,t-1] + noise[:,t]
    x = x[:,BURN:]; return x - x.min() + 1.0

def rebin(x,b):
    if b<=1: return x
    n,Tt=x.shape; T2=(Tt//b)*b
    return x[:,:T2].reshape(n,T2//b,b).sum(2)

def fit(x, kmax):
    nb = x.shape[1]
    if nb <= kmax+2: return np.nan, np.nan, "short"
    try:
        rk = mre.coefficients(x, steps=(1,kmax), dt=1, method="trialseparated")
        ft = mre.fit(rk, fitfunc="exponential"); m,tau = float(ft.mre), float(ft.tau)
    except Exception as e: return np.nan, np.nan, f"fail:{type(e).__name__}"
    if not (np.isfinite(m) and np.isfinite(tau)): return m,tau,"nonfinite"
    if tau < 2: return m,tau,"tau_small"
    if tau > MAX_TAU_FRAC*kmax: return m,tau,"CENSORED"
    return m,tau,"ok"

rng = np.random.default_rng(SEED)
print(f"induced ceiling exp(-1/48) = {CEIL:.5f}\n")
print(f"{'tau_s':>5} {'bin':>3} {'tau_bin_true':>12} {'m_true':>7} | "
      f"{'FIXED_m':>8} {'f_tau':>6} {'f_stat':>9} | {'TIME_m':>8} {'t_tau':>6} {'gap':>7}")
print("-"*95)
rows=[]
for ts in TAU_STEPS_GRID:
    x = ar1(ts, rng)
    for b in BIN_STEPS:
        xb = rebin(x,b)
        tbt = ts/b; mt = np.exp(-1/tbt)
        mf,tf,sf = fit(xb, KMAX_FIXED)
        kt = max(5, round(KMAX_STEPS_TIME/b))
        mtm,ttm,stm = fit(xb, kt)
        feff = mf if sf=="ok" else (CEIL if sf=="CENSORED" else mf)
        gap = (mtm - feff) if np.isfinite(mtm) and np.isfinite(feff) else np.nan
        rows.append(dict(tau_steps=ts,bin=b,tau_bin_true=tbt,m_true=mt,
                         fixed_m=mf,fixed_tau=tf,fixed_stat=sf,time_m=mtm,time_tau=ttm,gap=gap))
        print(f"{ts:>5} {b:>3} {tbt:>12.1f} {mt:>7.4f} | "
              f"{mf:>8.4f} {tf:>6.1f} {sf:>9} | {mtm:>8.4f} {ttm:>6.1f} {gap:>+7.4f}")

df = pd.DataFrame(rows)
print("\n--- PREDICTION CHECK ---")
finite = df[df.gap.notna()]
maxgap = finite.gap.abs().max()
print(f"max |FIXED vs TIME gap| across all cells: {maxgap:.4f}  (threshold for 'filter bites' was 0.02)")
print(f"any cell CENSORED by fixed rule: {(df.fixed_stat=='CENSORED').any()} "
      f"-> only where tau_bin_true > 48: {sorted(df[df.fixed_stat=='CENSORED'].tau_bin_true.round(1).tolist())}")
print(f"\nVERDICT: {'single-mechanism demo FAILS as predicted -> two-mechanism 3.3a' if maxgap < 0.02 else 'gap reached 0.02 -> reconsider'}")

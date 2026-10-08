#!/usr/bin/env python3
"""WAVE 20 analysis-only: k-DPP (determinantal point process) MAP selection + ordering,
head-to-head vs the shipped greedy-submodular portfolio.  NO GPU, NO writes outside
results/kdpp_* and the scratch dir.  Reuses the EXACT quality priors from wave19_portfolio.py.

Kernel  L = diag(q) @ S @ diag(q).  q = p_success (switch*bind*express).  S = similarity in [0,1]:
 selection (36 pool): S = mean(struct-TM, epitope-Jaccard, seq k-mer cosine)   [task's 3 axes]
 ordering (20)      : Tier-B/knockout pairs = same 3-axis mean; locked rows fall back to
                      mean(seq cosine, lineage/fold similarity) because the locked ORIGINAL
                      complexes are not under runs/rfd_obj1_w18 (only partial-diffusion cousins are).
k-DPP MAP = greedy max log det(L_S): start empty, add argmax marginal log-det gain.
"""
import csv, math, os, subprocess, sys
import numpy as np
from scipy.spatial import cKDTree

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"
sys.path.insert(0, os.path.join(REPO, "scripts"))
import wave19_portfolio as wp            # exact priors: P_EXPRESS, p_bind, p_switch, load_metrics, load_dev, portfolio_order
import developability_scan as devscan    # scan(name, seq) -> {... 'expr_risk': CLEAN/MODERATE/SEVERE}

W18 = os.path.join(REPO, "runs", "rfd_obj1_w18")
FOLDSEEK = os.path.join(REPO, "scratchpad", "foldseek", "bin", "foldseek")
SCRATCH = "/home/prachh/work/tmp/claude-1004/-home-prachh-v/80165829-a188-4ed2-ae27-49adbdc8ef4f/scratchpad/kdpp"
MONO = os.path.join(SCRATCH, "monomers")
POOL_CSV = os.path.join(REPO, "results", "md_selective_w19.csv")
METRICS_FINAL = os.path.join(REPO, "results", "metrics_final_w19.csv")
METRICS_20 = os.path.join(REPO, "results", "metrics_20.csv")
RESUBMIT = os.path.join(REPO, "results", "adaptyv_challenge1_resubmit_W19_FINAL.csv")
DEV_W19 = os.path.join(REPO, "results", "developability_adaptyv_challenge1_resubmit_W19.csv")
OUT_CSV = os.path.join(REPO, "results", "kdpp_selection_w20.csv")
OUT_MD = os.path.join(REPO, "results", "kdpp_report.md")
KNOCKOUT = "egfr_phsw_32"
LOCKED = ["egfr_phsw_02", "egfr_phsw_04", "egfr_phsw_06", "egfr_phsw_08",
          "egfr_phsw_09", "egfr_phsw_12", "egfr_phsw_16", "egfr_phsw_19"]

def log(*a): print(*a, flush=True)

# ---------- load tables ----------
def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))

def load_all():
    pool = read_csv(POOL_CSV)                 # 36 candidates: name(pdnov), hip,hie,docc,tier,fold,lin,pae,tm
    mf = {r["name"]: r for r in read_csv(METRICS_FINAL)}          # keyed by pdnov name
    mf_by_new = {r["new_name"]: r for r in read_csv(METRICS_FINAL) if r.get("new_name")}
    m20 = {r["name"]: r for r in read_csv(METRICS_20)}            # locked + knockout + (stale) tierB
    resub = read_csv(RESUBMIT)                 # 20 rows in portfolio order: name(egfr), sequence
    return pool, mf, mf_by_new, m20, resub

# ---------- quality q (exact wave19 priors) ----------
def q_from(d_occ, occ_HIE, bc, af2_pae, expr_risk):
    return wp.P_EXPRESS.get(expr_risk, wp.P_EXPRESS[""]) * wp.p_bind(bc, af2_pae) * wp.p_switch(d_occ, occ_HIE)

def ff(x):
    return wp.ff(x)

def bind_consensus_of(row, chai_key, boltz_key, bc_key=None):
    bc = ff(row.get(bc_key)) if bc_key else None
    if bc is None:
        chai = ff(row.get(chai_key)); boltz = ff(row.get(boltz_key))
        vals = [v for v in (chai, boltz) if v is not None]
        bc = min(vals) if vals else None
    return bc

# ---------- features ----------
def parse_pdb(path):
    """return (chainA_xyz Nx3, chainB_resatoms dict resnum->list xyz)."""
    A = []; B = {}
    with open(path) as fh:
        for ln in fh:
            if not ln.startswith("ATOM"): continue
            ch = ln[21]
            try:
                x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            except ValueError:
                continue
            if ch == "A":
                A.append((x, y, z))
            elif ch == "B":
                rn = int(ln[22:26])
                B.setdefault(rn, []).append((x, y, z))
    return np.array(A, float), B

def write_monomer(path, out):
    with open(path) as fh, open(out, "w") as o:
        for ln in fh:
            if ln.startswith("ATOM") and ln[21] == "A":
                o.write(ln)
        o.write("END\n")

def epitope_contacts(path, cutoff=5.0):
    A, B = parse_pdb(path)
    if len(A) == 0 or not B:
        return set()
    tree = cKDTree(A)
    hit = set()
    for rn, atoms in B.items():
        d, _ = tree.query(np.array(atoms, float), k=1)
        if np.min(d) <= cutoff:
            hit.add(rn)
    return hit

def jaccard(a, b):
    if not a and not b: return 0.0
    u = len(a | b)
    return len(a & b) / u if u else 0.0

def kmer_vec(seq, k=3):
    d = {}
    for i in range(len(seq) - k + 1):
        km = seq[i:i+k]
        d[km] = d.get(km, 0) + 1
    return d

def cosine(d1, d2):
    if not d1 or not d2: return 0.0
    dot = sum(v * d2.get(km, 0) for km, v in d1.items())
    n1 = math.sqrt(sum(v*v for v in d1.values())); n2 = math.sqrt(sum(v*v for v in d2.values()))
    return dot / (n1 * n2) if n1 and n2 else 0.0

def foldseek_tm(mono_dir, names):
    """all-vs-all TMalign; return symmetric TM matrix (query-normalized, averaged both directions)."""
    out = os.path.join(SCRATCH, "fs_allvall.m8"); tmp = os.path.join(SCRATCH, "fstmp")
    cmd = [FOLDSEEK, "easy-search", mono_dir, mono_dir, out, tmp,
           "--alignment-type", "1", "--tmscore-threshold", "0.0", "-e", "1000",
           "--max-seqs", "2000", "--threads", "16",
           "--format-output", "query,target,qtmscore,qcov"]
    r = subprocess.run(cmd, cwd=REPO, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    if r.returncode != 0:
        log("FOLDSEEK STDERR:", r.stderr.decode()[-800:])
        raise RuntimeError("foldseek failed")
    idx = {n: i for i, n in enumerate(names)}
    def resolve(tok):
        t = tok
        if t.endswith(".pdb"): t = t[:-4]
        if t in idx: return t
        if t.endswith("_A") and t[:-2] in idx: return t[:-2]
        for n in names:                         # prefix fallback
            if t == n or t.startswith(n): return n
        return None
    n = len(names); raw = np.full((n, n), np.nan)
    hits = 0
    for ln in open(out):
        p = ln.rstrip("\n").split("\t")
        if len(p) < 3: continue
        qi = resolve(p[0]); ti = resolve(p[1])
        if qi is None or ti is None: continue
        tm = float(p[2]); raw[idx[qi], idx[ti]] = tm; hits += 1
    log(f"  foldseek hit rows parsed: {hits}")
    sim = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            a = raw[i, j]; b = raw[j, i]
            vals = [v for v in (a, b) if not np.isnan(v)]
            sim[i, j] = float(np.mean(vals)) if vals else 0.0
    np.fill_diagonal(sim, 1.0)
    sim = 0.5 * (sim + sim.T)
    return sim

# ---------- PSD projection ----------
def psd_correlation(S, eps=1e-8):
    S = 0.5 * (S + S.T)
    w, V = np.linalg.eigh(S)
    neg = float(-w[w < 0].sum())
    w = np.clip(w, eps, None)
    Sp = (V * w) @ V.T
    d = np.sqrt(np.clip(np.diag(Sp), eps, None))
    Sc = Sp / np.outer(d, d)
    Sc = 0.5 * (Sc + Sc.T)
    np.fill_diagonal(Sc, 1.0)
    return Sc, neg

# ---------- log det helpers ----------
def slogdet_sub(L, idxs):
    if len(idxs) == 0: return 0.0
    M = L[np.ix_(idxs, idxs)]
    M = M + 1e-9 * np.eye(len(idxs))
    s, ld = np.linalg.slogdet(M)
    if s <= 0:
        M = M + 1e-6 * np.eye(len(idxs))
        s, ld = np.linalg.slogdet(M)
    return float(ld)

def greedy_kdpp(L, k, forced_last=None, candidates=None):
    """greedy MAP: pick argmax marginal log-det gain. candidates = indices eligible; forced_last appended."""
    n = L.shape[0]
    cand = list(range(n)) if candidates is None else list(candidates)
    if forced_last is not None and forced_last in cand:
        cand.remove(forced_last)
    chosen = []
    while len(chosen) < k and cand:
        base = slogdet_sub(L, chosen)
        best, bestgain = None, -1e18
        for c in cand:
            g = slogdet_sub(L, chosen + [c]) - base
            if g > bestgain + 1e-12:
                bestgain, best = g, c
        chosen.append(best); cand.remove(best)
    if forced_last is not None:
        chosen.append(forced_last)
    return chosen

# =====================================================================================
def main():
    os.makedirs(MONO, exist_ok=True)
    pool, mf, mf_by_new, m20, resub = load_all()
    pool_names = [r["name"] for r in pool]                 # 36 pdnov names
    assert len(pool_names) == 36, len(pool_names)

    # ---- map pdnov<->egfr for selected tierB; collect sequences & lineage/fold ----
    pdnov2new = {r["name"]: r["new_name"] for r in read_csv(METRICS_FINAL) if r.get("new_name")}
    seq_pool = {r["name"]: r["sequence"] for r in read_csv(METRICS_FINAL)}
    lin_pool = {r["name"]: r["lineage"] for r in read_csv(METRICS_FINAL)}
    fold_pool = {r["name"]: r["fold"] for r in read_csv(METRICS_FINAL)}

    # ---- portfolio's 11 selected tierB (selected==yes) ----
    port11_pdnov = [r["name"] for r in read_csv(METRICS_FINAL) if r.get("selected") == "yes"]
    assert len(port11_pdnov) == 11, (len(port11_pdnov), port11_pdnov)
    log("portfolio-11 (pdnov):", [pdnov2new[p] for p in port11_pdnov])

    # ---- developability expr_risk for all 36 pool via scan() (fair); plus W19 for the 20 ----
    dev_pool = {p: devscan.scan(p, seq_pool[p])["expr_risk"] for p in pool_names}
    dev20 = wp.load_dev(DEV_W19)              # keyed by egfr name
    # self-consistency: scan vs W19 for the 11 selected
    mism = []
    for p in port11_pdnov:
        new = pdnov2new[p]
        if new in dev20 and dev20[new] != dev_pool[p]:
            mism.append((new, dev20[new], dev_pool[p]))
    log("expr_risk scan-vs-W19 mismatches (selected-11):", mism or "none")

    # ---- quality q for the 36 pool (exact priors) ----
    qp = {}
    for r in pool:
        p = r["name"]
        d_occ = ff(r["docc"]); occ_HIE = ff(r["hie"])   # md_selective hie/docc = pool MD
        mfr = mf[p]
        bc = bind_consensus_of(mfr, "chai_iptm_human", "boltz2_iptm", "bind_consensus")
        af2 = ff(mfr.get("af2_pae"))
        qp[p] = q_from(d_occ, occ_HIE, bc, af2, dev_pool[p])
    log("pool q range:", round(min(qp.values()),4), "-", round(max(qp.values()),4))

    # ---- features for the 36 pool ----
    log("extracting monomers + epitope contacts (36 pool)...")
    contacts = {}
    kvecs = {}
    for p in pool_names:
        pdb = os.path.join(W18, p + "__best.pdb")
        write_monomer(pdb, os.path.join(MONO, p + ".pdb"))
        contacts[p] = epitope_contacts(pdb)
        kvecs[p] = kmer_vec(seq_pool[p])
    csz = {p: len(contacts[p]) for p in pool_names}
    log("  epitope contact-set sizes: min %d max %d mean %.1f" % (min(csz.values()), max(csz.values()), np.mean(list(csz.values()))))

    log("foldseek all-vs-all TM (36 monomers)...")
    TM = foldseek_tm(MONO, pool_names)        # 36x36

    n = 36; idx = {p: i for i, p in enumerate(pool_names)}
    EPI = np.zeros((n, n)); SEQ = np.zeros((n, n))
    for i, pi in enumerate(pool_names):
        for j, pj in enumerate(pool_names):
            EPI[i, j] = jaccard(contacts[pi], contacts[pj]) if i != j else 1.0
            SEQ[i, j] = cosine(kvecs[pi], kvecs[pj]) if i != j else 1.0
    S_pool = (TM + EPI + SEQ) / 3.0
    np.fill_diagonal(S_pool, 1.0)
    log("  mean off-diag: TM %.3f  EPI %.3f  SEQ %.3f  S %.3f" %
        (offmean(TM), offmean(EPI), offmean(SEQ), offmean(S_pool)))
    Sc_pool, neg = psd_correlation(S_pool)
    log("  PSD projection clipped negative eigen-mass: %.4g" % neg)

    qv = np.array([qp[p] for p in pool_names])
    L_pool = (qv[:, None] * Sc_pool) * qv[None, :]

    # ---- STEP A: k-DPP MAP select 11 from 36 ----
    log("\n=== STEP A: k-DPP MAP selection (11 of 36) ===")
    sel_idx = greedy_kdpp(L_pool, 11)
    kdpp11 = [pool_names[i] for i in sel_idx]
    log("k-DPP-11 (pdnov):", [pdnov2new.get(p, p) for p in kdpp11])

    port_idx = [idx[p] for p in port11_pdnov]
    kdpp_set, port_set = set(kdpp11), set(port11_pdnov)
    vol_kdpp = slogdet_sub(L_pool, [idx[p] for p in kdpp11])
    vol_port = slogdet_sub(L_pool, port_idx)
    meanq_kdpp = float(np.mean([qp[p] for p in kdpp11]))
    meanq_port = float(np.mean([qp[p] for p in port11_pdnov]))
    log("volume logdet(L_S): k-DPP %.4f   portfolio %.4f   delta %.4f" % (vol_kdpp, vol_port, vol_kdpp - vol_port))
    log("mean q:             k-DPP %.4f   portfolio %.4f" % (meanq_kdpp, meanq_port))
    added = kdpp_set - port_set; dropped = port_set - kdpp_set
    log("k-DPP adds:", [pdnov2new.get(p, p) for p in added])
    log("k-DPP drops:", [pdnov2new.get(p, p) for p in dropped])

    # ---- STEP A2: LOCKED-CONDITIONED k-DPP (the correct best-of objective) ----
    # The 8 locked are ALREADY in the submission. Tier-B should add volume ON TOP of that floor.
    # Pool-only Step A cannot see that (e.g. lineage-09 is already covered by locked egfr_phsw_09).
    log("\n=== STEP A2: locked-conditioned k-DPP (floor = 8 locked) ===")
    mm = wp.load_metrics()
    seq_locked = {r["name"]: r["sequence"] for r in resub if r["name"] in LOCKED}
    q_locked = {}
    for l in LOCKED:
        m = mm.get(l, {})
        q_locked[l] = q_from(m.get("d_occ"), m.get("occ_HIE"), m.get("bind_consensus"), m.get("af2_pae"), dev20.get(l, ""))
    lin44 = {p: lin_pool[p] for p in pool_names}; lin44.update({l: l for l in LOCKED})
    fold44 = {p: fold_pool[p] for p in pool_names}; fold44.update({l: l for l in LOCKED})
    kv44 = dict(kvecs); kv44.update({l: kmer_vec(seq_locked[l]) for l in LOCKED})
    names44 = pool_names + LOCKED
    i44 = {nm: i for i, nm in enumerate(names44)}
    def linfold44(a, b):
        if lin44[a] == lin44[b]: return 1.0
        if fold44[a] and fold44[a] == fold44[b]: return 0.6
        return 0.0
    N = 44; S44 = np.eye(N)
    for a in names44:
        for b in names44:
            ia, ib = i44[a], i44[b]
            if ia >= ib: continue
            if a in idx and b in idx:                      # both pool -> 3-axis
                s = (TM[idx[a], idx[b]] + EPI[idx[a], idx[b]] + SEQ[idx[a], idx[b]]) / 3.0
            else:                                          # locked-involved -> fallback
                s = (cosine(kv44[a], kv44[b]) + linfold44(a, b)) / 2.0
            S44[ia, ib] = S44[ib, ia] = s
    Sc44, neg44 = psd_correlation(S44)
    q44 = np.array([(qp[nm] if nm in idx else q_locked[nm]) for nm in names44])
    L44 = (q44[:, None] * Sc44) * q44[None, :]
    locked_idx = [i44[l] for l in LOCKED]
    pool_idx44 = [i44[p] for p in pool_names]
    # conditioned greedy: floor = locked, then add 11 from pool
    chosen = list(locked_idx)
    cond_sel = []
    while len(cond_sel) < 11:
        base = slogdet_sub(L44, chosen)
        best, bg = None, -1e18
        for c in pool_idx44:
            if c in chosen: continue
            g = slogdet_sub(L44, chosen + [c]) - base
            if g > bg + 1e-12: bg, best = g, c
        chosen.append(best); cond_sel.append(best)
    cond11 = [names44[i] for i in cond_sel]
    cond_set = set(cond11)
    log("conditioned-11 (pdnov):", [pdnov2new.get(p, p) for p in cond11])
    log("conditioned adds vs portfolio:", [pdnov2new.get(p, p) for p in (cond_set - port_set)])
    log("conditioned drops vs portfolio:", [pdnov2new.get(p, p) for p in (port_set - cond_set)])
    # marginal volume added ON TOP of the locked floor
    floor_vol = slogdet_sub(L44, locked_idx)
    mv_cond = slogdet_sub(L44, locked_idx + [i44[p] for p in cond11]) - floor_vol
    mv_port = slogdet_sub(L44, locked_idx + [i44[p] for p in port11_pdnov]) - floor_vol
    mv_poolA = slogdet_sub(L44, locked_idx + [i44[p] for p in kdpp11]) - floor_vol
    log("marginal volume on locked floor: conditioned %.4f  portfolio %.4f  poolA %.4f" % (mv_cond, mv_port, mv_poolA))
    meanq_cond = float(np.mean([(qp[p] if p in idx else q_locked[p]) for p in cond11]))
    log("conditioned mean q: %.4f" % meanq_cond)

    # ---- STEP B: k-DPP-MAP ORDER the VALIDATED portfolio-20 (zero new-design risk) ----
    # Recommended action is order-only, so we order the CURRENT validated 20 (8 locked + portfolio-11
    # + knockout); the selection swap (A2) is reported separately as an optional, user-gated upgrade.
    log("\n=== STEP B: k-DPP-MAP ordering of the validated 20 ===")
    seq20 = {r["name"]: r["sequence"] for r in resub}      # egfr -> seq (authoritative live)
    order20_names = LOCKED + [pdnov2new[p] for p in port11_pdnov] + [KNOCKOUT]
    assert len(order20_names) == 20 and len(set(order20_names)) == 20, order20_names

    recs = {}
    m20_lin = {r["name"]: r.get("lineage", "") for r in read_csv(METRICS_20)}
    metrics_merged = wp.load_metrics()
    for name in order20_names:
        if name == KNOCKOUT:
            role = "knockout"; seq = seq20[name]; lin = m20_lin.get(name, name); fold = lin
            m = metrics_merged.get(name, {})
            q = q_from(m.get("d_occ"), m.get("occ_HIE"), m.get("bind_consensus"), m.get("af2_pae"), dev20.get(name, ""))
            cpath = os.path.join(W18, "pdnov_egfr_phsw_06__des_65_2__best.pdb")  # E66Q 1-mut of egfr_phsw_23 (=this pdnov); proxy
            recs[name] = dict(role=role, seq=seq, lin=lin, fold=fold, q=q, cpath=cpath, has_cx=True, pdnov="(E66Q knockout of egfr_phsw_23)")
        elif name in LOCKED:
            role = "locked"; seq = seq20[name]; lin = m20_lin.get(name, name); fold = lin
            m = metrics_merged.get(name, {})
            q = q_from(m.get("d_occ"), m.get("occ_HIE"), m.get("bind_consensus"), m.get("af2_pae"), dev20.get(name, ""))
            recs[name] = dict(role=role, seq=seq, lin=lin, fold=fold, q=q, cpath=None, has_cx=False, pdnov="(locked original; no w18 complex)")
        else:  # tierB (portfolio-selected)
            role = "tierB"
            pdnov = [p for p in port11_pdnov if pdnov2new[p] == name][0]
            recs[name] = dict(role=role, seq=seq_pool[pdnov], lin=lin_pool[pdnov], fold=fold_pool[pdnov], q=qp[pdnov],
                              cpath=os.path.join(W18, pdnov + "__best.pdb"), has_cx=True, pdnov=pdnov)
    for name in order20_names:
        r = recs[name]
        r["kvec"] = kmer_vec(r["seq"])
        r["contacts"] = epitope_contacts(r["cpath"]) if r["has_cx"] else None
    cx_names = [nm for nm in order20_names if recs[nm]["has_cx"]]
    cxdir = os.path.join(SCRATCH, "mono20"); os.makedirs(cxdir, exist_ok=True)
    for f in os.listdir(cxdir):
        os.remove(os.path.join(cxdir, f))
    for nm in cx_names:
        write_monomer(recs[nm]["cpath"], os.path.join(cxdir, nm + ".pdb"))
    TM20 = foldseek_tm(cxdir, cx_names)
    tmidx = {nm: i for i, nm in enumerate(cx_names)}

    def linfold_sim(a, b):
        if recs[a]["lin"] == recs[b]["lin"]: return 1.0
        if recs[a]["fold"] and recs[a]["fold"] == recs[b]["fold"]: return 0.6
        return 0.0

    m = 20; S20 = np.eye(m)
    for i, a in enumerate(order20_names):
        for j, b in enumerate(order20_names):
            if i >= j: continue
            sq = cosine(recs[a]["kvec"], recs[b]["kvec"])
            if recs[a]["has_cx"] and recs[b]["has_cx"]:
                s = (TM20[tmidx[a], tmidx[b]] + jaccard(recs[a]["contacts"], recs[b]["contacts"]) + sq) / 3.0
            else:
                s = (sq + linfold_sim(a, b)) / 2.0     # locked fallback
            S20[i, j] = S20[j, i] = s
    Sc20, neg20 = psd_correlation(S20)
    log("  ordering PSD clipped negative mass: %.4g" % neg20)
    q20 = np.array([recs[nm]["q"] for nm in order20_names])
    L20 = (q20[:, None] * Sc20) * q20[None, :]
    o_idx = {nm: i for i, nm in enumerate(order20_names)}

    ko_i = order20_names.index(KNOCKOUT)
    ord_idx = greedy_kdpp(L20, k=19, forced_last=ko_i)
    ordered = [order20_names[i] for i in ord_idx]
    port_order = [r["name"] for r in resub]  # live portfolio order of the 20 (egfr names)
    log("k-DPP order :", ordered)
    log("portfolio   :", port_order)

    # ---- best-of top-T coverage diagnostics (volume + max q + mean q) ----
    def topT_stats(order, T):
        names = order[:T]; idxs = [o_idx[nm] for nm in names]
        return slogdet_sub(L20, idxs), max(recs[nm]["q"] for nm in names), float(np.mean([recs[nm]["q"] for nm in names]))
    log("\n  best-of top-T coverage (L20 logdet / max q / mean q):")
    topT = {}
    for T in (3, 5, 8):
        vk, mk, ak = topT_stats(ordered, T); vp, mp, ap = topT_stats(port_order, T)
        log("   T=%d  k-DPP: vol %.3f maxq %.3f meanq %.3f | portfolio: vol %.3f maxq %.3f meanq %.3f"
            % (T, vk, mk, ak, vp, mp, ap))
        topT[T] = dict(kdpp_vol=vk, kdpp_maxq=mk, kdpp_meanq=ak, port_vol=vp, port_maxq=mp, port_meanq=ap)

    # ---- STEP D: self-check ----
    ok_unique = len(set(ordered)) == 20
    ok_ko_last = ordered[-1] == KNOCKOUT
    ok_locked = all(l in ordered for l in LOCKED)
    ok_count_tierB = sum(1 for nm in ordered if recs[nm]["role"] == "tierB") == 11
    log("\n=== STEP D self-check ===")
    log("  20 unique: %s | knockout last: %s | all 8 locked present: %s | exactly 11 tierB: %s"
        % (ok_unique, ok_ko_last, ok_locked, ok_count_tierB))
    assert ok_unique and ok_ko_last and ok_locked and ok_count_tierB, "SELF-CHECK FAILED"

    # ---- selection flags from the PRINCIPLED (locked-conditioned) k-DPP ----
    cond_dropped = port_set - cond_set          # pdnov dropped by conditioned k-DPP
    cond_added = cond_set - port_set            # pdnov added
    poolA_added = kdpp_set - port_set; poolA_dropped = port_set - kdpp_set

    # ---- write results/kdpp_selection_w20.csv (validated 20 in k-DPP order) ----
    rank_of = {nm: i + 1 for i, nm in enumerate(ordered)}
    rows_out = []
    for nm in ordered:
        r = recs[nm]; is_tb = r["role"] == "tierB"; pdnov = r["pdnov"]
        in_port = "yes" if is_tb else "-"
        in_kd = ("no" if pdnov in cond_dropped else "yes") if is_tb else "-"
        rows_out.append(dict(rank=rank_of[nm], name=nm, pdnov_name=pdnov, role=r["role"],
                             q=round(r["q"], 4), in_portfolio_11=in_port, in_kdpp_11=in_kd))
    for p in sorted(cond_added):                # conditioned-add (optional swap-in), no rank/slot yet
        rows_out.append(dict(rank="", name="(swap-in candidate)", pdnov_name=p, role="tierB(kdpp-add)",
                             q=round(qp[p], 4), in_portfolio_11="no", in_kdpp_11="yes"))
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rank", "name", "pdnov_name", "role", "q", "in_portfolio_11", "in_kdpp_11"])
        w.writeheader(); w.writerows(rows_out)
    log("\nwrote", OUT_CSV)

    res = dict(
        port11=[pdnov2new[p] for p in port11_pdnov], port11_pdnov=list(port11_pdnov),
        poolA_11=[pdnov2new.get(p, p) for p in kdpp11], poolA_added=sorted(poolA_added), poolA_dropped=sorted(poolA_dropped),
        cond11=[pdnov2new.get(p, p) for p in cond11], cond11_pdnov=list(cond11),
        cond_added=sorted(cond_added), cond_dropped=sorted(cond_dropped),
        cond_dropped_new=[pdnov2new[p] for p in cond_dropped],
        vol_kdpp_poolA=vol_kdpp, vol_port=vol_port, meanq_poolA=meanq_kdpp, meanq_port=meanq_port,
        mv_cond=mv_cond, mv_port=mv_port, mv_poolA=mv_poolA, meanq_cond=meanq_cond,
        ordered=ordered, port_order=port_order, neg=neg, neg20=neg20, neg44=neg44,
        offTM=offmean(TM), offEPI=offmean(EPI), offSEQ=offmean(SEQ), offS=offmean(S_pool),
        qp={pdnov2new.get(p, p): round(qp[p], 4) for p in pool_names},
        q20={nm: round(recs[nm]["q"], 4) for nm in order20_names},
        role20={nm: recs[nm]["role"] for nm in order20_names},
        pdnov_of={nm: recs[nm]["pdnov"] for nm in order20_names},
        lin20={nm: recs[nm]["lin"] for nm in order20_names},
        topT=topT, mism=mism)
    import json
    with open(os.path.join(SCRATCH, "res.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    log("stashed res.json")
    return res

def offmean(M):
    n = M.shape[0]
    return float((M.sum() - np.trace(M)) / (n * (n - 1)))

if __name__ == "__main__":
    main()

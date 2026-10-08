#!/usr/bin/env python3
"""Structure-based surface-hydrophobic-patch (spatial-aggregation-propensity) scan.

WHY: developability_scan.py flags aggregation with whole-sequence GRAVY + hydrophobic-run
proxies, which miss what actually depresses expression/solubility: SOLVENT-EXPOSED,
SPATIALLY-CLUSTERED hydrophobic surface (Chennamsetty 2009 PNAS SAP; Lauer 2012 J Pharm Sci
developability index; the Therapeutic Antibody Profiler hydrophobic-patch score). A buried
hydrophobic is harmless; a contiguous exposed patch nucleates aggregation.

WHAT: scores each binder in ISOLATION (target stripped) -- aggregation is a property of the
free binder in solution, so the interface face is legitimately exposed. Per binder:
  max_patch_hSASA  largest contiguous exposed-hydrophobic patch, A^2 (CB-CB < PATCH_D graph)
  patch_frac_surf  max_patch_hSASA / total molecular SASA, %  <-- PRIMARY, size-INVARIANT flag metric
  exphyd_frac      fraction of residues that are exposed hydrophobics
  patch_per_res    max_patch_hSASA / nres (kept for continuity; size-BIASED toward small binders -> not the flag)

VERDICT is RELATIVE and flagged on the size-invariant patch_frac_surf (% of the molecular surface
taken up by the single largest hydrophobic patch): a design is flagged only if patch_frac_surf
exceeds the 95th percentile of a control distribution of OTHER de-novo pool-passer binders (same
design class). The natural well-folded EGFR Domain III is scored too but used only as a labeled
CONTEXT FLOOR (a natural domain minimises exposed hydrophobic surface) -- it is NOT the flag
threshold. An absolute patch number alone is meaningless for a small helix bundle whose core reads
as "exposed" once the partner is removed, hence the class-relative percentile.

Matches structures to submitted sequences BY SEQUENCE (exact, else Hamming<=2 to tolerate the
WAVE-16d Cys-clean point mutations vs a pre-clean fold). Reads .pdb and .cif. PyRosetta+numpy.
Read-only; writes results/sap_scan.csv. Never modifies a submission.

Usage: envs/pyrosetta/bin/python scripts/sap_scan.py
"""
import os
import csv
from pathlib import Path
from collections import OrderedDict
import numpy as np

REPO = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SUB  = REPO/"results/submission_track3_v2.csv"
HYD  = set("AVLIMFWY")
EXP_THRESH = 10.0     # A^2 per-res hydrophobic SASA => "exposed enough to join a patch"
PATCH_D    = 8.0      # CB-CB (A) linking two exposed hydrophobics into one patch
# dirs holding submit==folded structures. af2ig folds FIRST (consistent with the af2ig reference
# pool -> removes the boltz-vs-af2ig folding-method confound for the 4 rfd3 designs); boltz co-folds
# are a last-resort fallback only if a design has no exact af2ig fold.
FOLD_DIRS = ["runs/sap15_fix_af2ig",  # WAVE-16f: #15 surface-fix (v3) re-fold — submit==folded for #15
             "runs/cys_clean_af2ig01","runs/cys_clean_af2ig04","runs/cys_clean_af2ig_batch",
             "runs/cys_clean_af2ig_iface","runs/cys_clean_af2ig_iface2","runs/md_pd3_in",
             "runs/eng_foldin","runs/rfd3/af2ig_c700","runs/rfd3/af2ig_c800",
             "runs/rfd3/af2ig_switch2_0cys","runs/rfd3/af2ig_switch3_0cys",
             "runs/boltz_switch2full","runs/boltz_campaign",
             "runs/boltz_c800_boltz_g0","runs/boltz_c800_boltz_g1","runs/boltz_c800_boltz_g2"]
REF_DIR = "runs/eng_foldin_pd3"   # reference pool: other de-novo pd3 binder chains (class control)
AA3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
          'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
          'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def chain_seqs_pdb(pdb):
    ch = OrderedDict()
    for ln in open(pdb):
        if ln.startswith("ENDMDL"): break                       # first model only (multi-model guard)
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA" and ln[16] in " A":  # skip altloc B/C...
            ch.setdefault(ln[21], []).append(AA3to1.get(ln[17:20].strip(), 'X'))
    return {c: "".join(s) for c, s in ch.items()}

def chain_seqs_cif(cif):
    out = OrderedDict()
    lines = open(cif).read().splitlines(); i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j = i+1; tags = []
            while j < len(lines) and lines[j].strip().startswith("_atom_site."):
                tags.append(lines[j].strip()); j += 1
            if tags:
                idx = {t: k for k, t in enumerate(tags)}
                ca = idx.get("_atom_site.label_atom_id"); co = idx.get("_atom_site.label_comp_id")
                asym = idx.get("_atom_site.auth_asym_id", idx.get("_atom_site.label_asym_id"))
                alt = idx.get("_atom_site.label_alt_id"); modn = idx.get("_atom_site.pdbx_PDB_model_num")
                first_model = None; k = j
                while k < len(lines) and lines[k].strip() and not lines[k].strip().startswith(("_","loop_","#")):
                    f = lines[k].split()
                    if len(f) >= len(tags) and f[ca].strip('"') == "CA":
                        if modn is not None and first_model is None: first_model = f[modn]
                        ok_alt = alt is None or f[alt] in (".", "A")        # skip altloc B/C...
                        ok_mod = modn is None or f[modn] == first_model     # first model only
                        if ok_alt and ok_mod:
                            out.setdefault(f[asym], []).append(AA3to1.get(f[co], 'X'))
                    k += 1
                i = k; continue
        i += 1
    return {c: "".join(s) for c, s in out.items()}

def chain_seqs(p):
    return chain_seqs_cif(p) if p.suffix == ".cif" else chain_seqs_pdb(p)

def ham(a, b):
    return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 999

def build_index(dirs):
    """seq -> (path, chain); keep exact-preferred. Returns list of (path,chain,seq) too for fuzzy."""
    exact = {}; allrecs = []
    for d in dirs:
        for p in sorted((REPO/d).rglob("*.pdb")) + sorted((REPO/d).rglob("*.cif")):
            for c, seq in chain_seqs(p).items():
                if len(seq) < 10: continue
                allrecs.append((p, c, seq))
                exact.setdefault(seq, (p, c))
    return exact, allrecs

def locate(sub_seq, exact, allrecs):
    if sub_seq in exact:
        return (*exact[sub_seq], 0)
    # fuzzy: FIRST record within Ham<=2 in FOLD_DIRS order -> honors af2ig-first (allrecs is dir-ordered),
    # so a Cys-clean near-twin resolves to an af2ig fold before any boltz co-fold (keeps the confound control).
    for p, c, seq in allrecs:
        if ham(seq, sub_seq) <= 2:
            return (p, c, ham(seq, sub_seq))
    return None  # (path,chain,ham) or None

def score_binder(pose):
    from pyrosetta.rosetta.core.scoring.sasa import SasaCalc
    n = pose.total_residue()
    if n == 0: raise ValueError("empty pose")
    sc = SasaCalc(); sc.calculate(pose)
    hsasa = list(sc.get_residue_hsasa()); tot_sasa = float(sc.get_total_sasa())
    exposed = []
    for i in range(1, n+1):
        aa = pose.residue(i).name1(); h = hsasa[i-1]
        if aa in HYD and h > EXP_THRESH:
            r = pose.residue(i)
            xyz = r.xyz("CB") if r.has("CB") else r.xyz("CA")
            exposed.append((i, aa, h, np.array([xyz.x, xyz.y, xyz.z])))
    m = len(exposed); parent = list(range(m))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a in range(m):
        for b in range(a+1, m):
            if np.linalg.norm(exposed[a][3]-exposed[b][3]) < PATCH_D: parent[find(a)] = find(b)
    comp = {}
    for a in range(m): comp.setdefault(find(a), []).append(a)
    patches = sorted(([sum(exposed[k][2] for k in mem),
                       [f"{exposed[k][1]}{exposed[k][0]}" for k in mem]] for mem in comp.values()),
                     reverse=True)
    maxp = patches[0] if patches else [0.0, []]
    toth = float(sum(hsasa))
    # PRIMARY size-invariant metric: % of total molecular surface in the single largest hydrophobic patch
    # (fixes patch_per_res, which divides by nres and biases flags toward small binders). patch_per_res kept
    # for continuity only.
    patch_frac = round(100.0 * maxp[0] / tot_sasa, 2) if tot_sasa > 0 else 0.0
    return dict(nres=n, tot_sasa=round(tot_sasa, 0), hSASA_per_res=round(toth/n, 2), exphyd_frac=round(m/n, 3),
                max_patch_hSASA=round(maxp[0], 1), patch_frac_surf=patch_frac, patch_per_res=round(maxp[0]/n, 2),
                n_patch_res=len(maxp[1]), patch_members=";".join(maxp[1]))

def load_matching_subpose(path, target_seq):
    """Load structure, return (sub-pose, ham) for the chain whose sequence best matches target_seq,
    or None if PyRosetta built no usable chain. Caller must gate on ham (a huge ham = wrong chain)."""
    from pyrosetta import pose_from_file
    pose = pose_from_file(str(path))
    best = None
    for sub in pose.split_by_chain():
        if sub.total_residue() == 0: continue
        h = ham(sub.sequence(), target_seq)
        if best is None or h < best[1]: best = (sub, h)
    return best  # (subpose, ham) or None

def main():
    from pyrosetta import init
    init("-mute all -ignore_unrecognized_res -ignore_zero_occupancy false")
    rows = list(csv.DictReader(open(SUB)))
    exact, allrecs = build_index(FOLD_DIRS)
    print(f"[index] {len(allrecs)} chains across {len(FOLD_DIRS)} fold dirs\n")

    results = []
    for r in rows:
        name, seq = r["name"], r["sequence"].strip().upper()
        hit = locate(seq, exact, allrecs)
        if not hit:
            print(f"[MISS] {name}  (no fold within Ham<=2 of the submitted seq)")
            results.append(dict(name=name, status="NO_STRUCT")); continue
        path, chain, h = hit
        try:
            m = load_matching_subpose(path, seq)
            if m is None or m[1] > 2:                       # guard: wrong/empty chain, don't score garbage
                raise ValueError(f"no chain within Ham<=2 (best {None if m is None else m[1]})")
            sub, sh = m
            rec = dict(name=name, status="OK", seq_ham=sh, src=str(path.relative_to(REPO)))
            rec.update(score_binder(sub))
        except Exception as e:                              # one bad structure must not abort all 20
            print(f"[SCORE-FAIL] {name}: {e}")
            rec = dict(name=name, status="SCORE_FAIL")
        results.append(rec)

    # reference pool: other de-novo pd3 binder chains (class control) + natural DIII (labeled context FLOOR)
    print("[ref] scoring class-control pool (pd3 pool binders + natural EGFR DIII)...")
    refs = []; diii = []; pool_fail = 0
    pool_pdbs = sorted((REPO/REF_DIR).glob("*.pdb"))
    sub_seqs = [r["sequence"].strip().upper() for r in rows]
    step = max(1, len(pool_pdbs)//50)   # ~50 pool structures
    for p in pool_pdbs[::step]:
        cs = chain_seqs(p)
        if not cs: continue
        items = sorted(cs.items(), key=lambda kv: len(kv[1]))   # binder = shortest chain; DIII = longest
        bseq = items[0][1]
        # exclude any pool binder within Ham<=2 of a submitted seq — the Cys-clean near-twins would otherwise
        # pollute the control (inflate p95 -> under-flag, the UNSAFE direction for a developability gate).
        if any(ham(bseq, s) <= 2 for s in sub_seqs): continue
        try:
            sb = load_matching_subpose(p, bseq)
            if sb is None or sb[1] > 2: raise ValueError("pool chain mismatch")
            refs.append(score_binder(sb[0])["patch_frac_surf"])
        except Exception:
            pool_fail += 1                                   # count, don't silently swallow (#7)
        if not diii:
            try:
                sd = load_matching_subpose(p, items[-1][1])
                if sd is not None and sd[1] <= 2: diii.append(score_binder(sd[0]))
            except Exception: pass
    if pool_fail: print(f"[ref] WARNING: {pool_fail} pool structures failed scoring (excluded)")
    if len(refs) < 20:                                       # empty/shrunken-pool guard (#1): fail loud
        raise SystemExit(f"[FATAL] reference pool too small (n={len(refs)}; need >=20). "
                         f"runs/ is gitignored/local-only — run this in the populated workspace.")
    refs = np.array(refs)
    print(f"[ref] n={len(refs)} pool binders; patch_frac_surf(%) mean {refs.mean():.2f} sd {refs.std():.2f} "
          f"p50 {np.percentile(refs,50):.2f} p90 {np.percentile(refs,90):.2f} p95 {np.percentile(refs,95):.2f}")
    if diii: print(f"[ref] natural EGFR DIII patch_frac_surf = {diii[0]['patch_frac_surf']}% "
                   f"(well-folded-domain FLOOR for context, NOT the flag threshold)")

    ok = [r for r in results if r.get("status") == "OK"]
    for r in ok:
        r["pool_pctl"] = round(float((refs < r["patch_frac_surf"]).mean()*100), 0)
    ok.sort(key=lambda x: -x["patch_frac_surf"])
    print(f"\n{'name':<14}{'n':>4}{'totSASA':>8}{'maxPatch':>9}{'patch%surf':>11}{'pool%ile':>9}  patch")
    for r in ok:
        print(f"{r['name']:<14}{r['nres']:>4}{r['tot_sasa']:>8.0f}{r['max_patch_hSASA']:>9}"
              f"{r['patch_frac_surf']:>11}{r['pool_pctl']:>8.0f}%  {r['patch_members'][:44]}")
    p95 = float(np.percentile(refs, 95))
    flagged = [r["name"] for r in ok if r["patch_frac_surf"] > p95]
    print(f"\nflag basis = patch_frac_surf (size-invariant: % of surface in largest hydrophobic patch); pool p95 = {p95:.2f}%")
    print(f"OUTLIERS above pool p95: {flagged or 'NONE -- all 20 within de-novo-binder class norms'}")
    cols = ['name','status','nres','tot_sasa','hSASA_per_res','exphyd_frac','max_patch_hSASA',
            'patch_frac_surf','patch_per_res','pool_pctl','n_patch_res','patch_members','seq_ham','src']
    out = REPO/"results/sap_scan.csv"
    with open(out,'w',newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore'); w.writeheader()
        for r in results: w.writerow(r)
    print(f"\nwrote {out}")

if __name__ == "__main__":
    main()

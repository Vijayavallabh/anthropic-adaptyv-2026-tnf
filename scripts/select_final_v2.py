#!/usr/bin/env python3
"""Final Track-3 selection v2 — merge ALL candidate sources into the best <=20, ranked by the
competition priority (pH-selective > cross-reactive > affinity), His433-driven + dual-confirmed
(PROPKA linkage AND MD salt-bridge) preferred.

Sources (each candidate kept as its best VIABLE version = binds AND carries the strongest switch):
  1. baseline  : results/submission_manifest.csv (already selected 20; their seqs/metrics)
  2. rfd_denovo: de novo RFd binders  (rfd_af2_passers_w1fa.csv [orig pae] + rfd_obj1_switch_w1fa.csv)
  3. pd        : partial-diffusion binders (rfd_af2_passers_pd.csv + rfd_obj1_switch_pd.csv)
For RFd/PD, each design has an ENGINEERED version (obj1 binder_seq, ddg_pH, re-af2ig pae from
rfd_af2ig_eng_<tag>) and an ORIGINAL version (af2pred chain-A seq, orig_ddg_pH, orig pae). We keep the
higher-ddg_pH version that still binds (pae < --max-pae). MD d_occupancy (runs/md_*) and cross-react
are attached when available.

Usage: select_final_v2.py [--max-pae 12] [--n 20] [--out ...] [--manifest ...]
"""
from __future__ import annotations
import argparse, csv, glob, json, os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
D3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H",
      "ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
      "TYR":"Y","VAL":"V","MSE":"M"}


def ff(d, k, dv=99.0):
    try: return float(d[k])
    except (KeyError, ValueError, TypeError): return dv


def family_of(design):
    """Backbone family for the diversity cap: A_s/A_l/C_s/C_l/H_s/conserved, or denovo. pd3/amplified
    children inherit their parent's family via the embedded egfr_dom3_<fam> tag; de-novo-lineage
    designs (w1_ext/w2_ext/..., and their pd3_w* children) all map to 'denovo' so the whole lineage is
    capped together -- otherwise pd3_w* children label as 'other' and escape the de-novo cap."""
    import re
    m = re.search(r"egfr_dom3_([A-Za-z]+(?:_[sl])?)", design)
    if m: return m.group(1)
    if "7d12" in design: return "nanobody"   # 7D12-seeded nanobody-derived backbones (distinct fold family)
    if re.search(r"w\d_(ext|tight)", design) or design.startswith("w"): return "denovo"
    return "other"


def chainA_seq(pdb):
    seq = []
    try:
        for l in open(pdb):
            if l[:4] == "ATOM" and l[21] == "A" and l[12:16].strip() == "CA":
                seq.append(D3.get(l[17:20].strip(), "X"))
    except OSError:
        return ""
    return "".join(seq)


def load_csv_map(path, key="design"):
    if not os.path.exists(path): return {}
    return {r[key]: r for r in csv.DictReader(open(path))}


def load_sc_map(pats):
    out = {}
    for pat in pats:
        for f in glob.glob(pat):
            hdr = None
            for ln in open(f):
                t = ln.split()
                if ln.startswith("SCORE:") and "description" in ln: hdr = t[1:]
                elif ln.startswith("SCORE:") and hdr:
                    d = dict(zip(hdr, t[1:]))
                    out[d["description"].replace("_af2pred","")] = d
    return out


def load_md_map():
    """design(stem) -> d_occupancy (>0.1 = MD switch). Covers baseline + candidate MD.
    md_eng_denovo is globbed LAST so the engineered-structure MD (run on the ACTUAL submitted
    engineered de-novo sequence) overrides the original-structure MD for those designs -- we must
    not claim a switch measured on a different sequence than the one we submit."""
    m = {}
    for j in (glob.glob(str(REPO/"runs"/"md_submission"/"*.json"))
              + glob.glob(str(REPO/"runs"/"md_candidates"/"*.json"))
              + glob.glob(str(REPO/"runs"/"md_final2"/"*.json"))
              + glob.glob(str(REPO/"runs"/"md_eng_denovo"/"*.json"))
              + glob.glob(str(REPO/"runs"/"md_pd3"/"*.json"))
              + glob.glob(str(REPO/"runs"/"md_7d12"/"*.json"))):
        # NOTE: runs/md_net (network 2-mutation seqs) is deliberately NOT globbed here -- its design
        # keys collide with the pd2/de-novo PARENTS, so it would overwrite the single-switch MD. The
        # network source loads md_net separately and attaches it only to its "__net" candidates.
        try:
            d = json.load(open(j))
            stem = d["name"].replace("__best","").replace("_af2pred","")
            if d.get("d_occupancy") is not None:
                m[stem] = d["d_occupancy"]
        except Exception: pass
    return m


def load_xreact_map():
    """design -> (cross_reactive_bool, mouse_pae, human_pae). From cross-react af2ig:
    mouse = design's divergent positions mutated human->mouse + repacked; humanrepack = repack control.
    Cross-reactive = the mouse ortholog still binds (mouse pae_int < 12) AND is not much worse than
    the repacked human control (mouse - human < 4)."""
    def scmap(run):
        out = {}
        for f in glob.glob(str(REPO/"runs"/run/"out.sc")) + glob.glob(str(REPO/"runs"/run/"gpu*"/"out.sc")):
            hdr = None
            for ln in open(f):
                t = ln.split()
                if ln.startswith("SCORE:") and "description" in ln: hdr = t[1:]
                elif ln.startswith("SCORE:") and hdr:
                    d = dict(zip(hdr, t[1:])); out[d["description"]] = d
        return out
    # _eng_ merged LAST: cross-react re-checked on the engineered de-novo sequence overrides the
    # original-sequence cross-react for those designs (consistency with the submitted sequence).
    mouse = {**scmap("af2ig_xreact_cand_mouse"), **scmap("af2ig_xreact_xr2_mouse"), **scmap("af2ig_xreact_eng_mouse"),
             **scmap("af2ig_xreact_pd3_mouse"), **scmap("af2ig_xreact_netwin_mouse"), **scmap("af2ig_xreact_7d12_mouse")}
    human = {**scmap("af2ig_xreact_cand_human"), **scmap("af2ig_xreact_xr2_human"), **scmap("af2ig_xreact_eng_human"),
             **scmap("af2ig_xreact_pd3_human"), **scmap("af2ig_xreact_netwin_human"), **scmap("af2ig_xreact_7d12_human")}
    res = {}
    for desc, d in mouse.items():
        if not desc.endswith("_mouse_af2pred"): continue
        base = desc[:-len("_mouse_af2pred")]
        mp = ff(d, "pae_interaction")
        hd = human.get(base + "_humanrepack_af2pred")
        if hd is None:
            # No human-repack control -> we CANNOT assert "mouse binds about as well as human".
            # Conservatively NOT cross-reactive (never default the missing control to 99, which would
            # make (mp-99)<4 trivially true and fabricate a cross-react positive on mouse-binding alone).
            res[base] = (False, round(mp, 1), None)
            continue
        hp = ff(hd, "pae_interaction")
        res[base] = ((mp < 12.0) and ((mp - hp) < 4.0), round(mp, 1), round(hp, 1))
    return res


def load_boltz_map():
    """design -> iptm (Boltz-2 single-seq cofold, independent binding predictor).
    Path shape: runs/boltz_shard*/out/by_seed/<DESIGN>/s0/confidence_<DESIGN>_model_0.json
    so the design name is the grandparent directory of the json."""
    m = {}
    for j in (glob.glob(str(REPO/"runs"/"boltz_shard*"/"out"/"**"/"confidence_*.json"), recursive=True)
              + glob.glob(str(REPO/"runs"/"boltz_fin2"/"out"/"**"/"confidence_*.json"), recursive=True)
              + glob.glob(str(REPO/"runs"/"boltz_eng"/"out"/"**"/"confidence_*.json"), recursive=True)
              + glob.glob(str(REPO/"runs"/"boltz_pd3"/"out"/"**"/"confidence_*.json"), recursive=True)
              + glob.glob(str(REPO/"runs"/"boltz_netwin"/"out"/"**"/"confidence_*.json"), recursive=True)
              + glob.glob(str(REPO/"runs"/"boltz_7d12"/"out"/"**"/"confidence_*.json"), recursive=True)):
        try:
            name = os.path.basename(os.path.dirname(os.path.dirname(j)))
            m[name] = round(float(json.load(open(j))["iptm"]), 3)
        except Exception: pass
    return m


def load_iface_hyd():
    """design -> interface hydrophobic fraction (precomputed cache). Lower = more polar interface =
    better wet-lab solubility/expression and fewer AF2 false-positives (the GuideFlip / direct-AF2
    hydrophobic-bias concern). Used as a gentle tiebreak within a switch-strength band, not a gate."""
    p = REPO/"results"/"iface_hyd.csv"
    if not p.exists(): return {}
    out = {}
    for r in csv.DictReader(open(p)):
        try: out[r["design"]] = float(r["iface_hyd"])
        except (ValueError, KeyError): pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-pae", type=float, default=12.0)
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--max-per-seed", type=int, default=3, help="cap variants per PD/RFd seed (diversity)")
    ap.add_argument("--max-per-family", type=int, default=99, help="cap picks per baseline family (A_s/C_s/C_l/A_l/H_s/conserved/denovo) for backbone diversity")
    ap.add_argument("--guarantee-families", action="store_true", help="ensure every validated backbone family (with a dual-confirmed, cross-reactive, Boltz>=0.5 member) is represented, swapping the weakest over-represented pick within --family-tol ddg_pH")
    ap.add_argument("--family-tol", type=float, default=0.2, help="max ddg_pH (kcal/mol) to sacrifice when swapping in a missing family's representative (PROPKA/MD noise band)")
    ap.add_argument("--out", default=str(REPO/"results"/"submission_track3_v2.csv"))
    ap.add_argument("--manifest", default=str(REPO/"results"/"submission_manifest_v2.csv"))
    ap.add_argument("--exclude", default="", help="comma-sep EXACT design names (or 'name*' for an explicit prefix) to drop before selection. NOTE: cross_reactive / cys are only SOFT sort tiebreaks here, not hard filters, so --exclude is the manual way to drop a design that should fail a gate (non-cross-reactive, 3+ Cys); for a real filter use --require-xreact / --max-cys.")
    ap.add_argument("--require-xreact", action="store_true", help="HARD filter: drop designs not predicted cross-reactive (human+mouse) before selection")
    ap.add_argument("--max-cys", type=int, default=None, help="HARD filter: drop designs with more than this many cysteines (e.g. --max-cys 2 forbids 3+ Cys)")
    a = ap.parse_args()

    md = load_md_map()
    xreact = load_xreact_map()
    boltz = load_boltz_map()
    iface = load_iface_hyd()
    cand = []

    # 1. baseline (already-selected 20): carry their metrics
    for r in csv.DictReader(open(REPO/"results"/"submission_manifest.csv")):
        cand.append({
            "source":"baseline","id":r["name"],"design":r["design"],"version":r.get("version",""),
            "sequence":r["sequence"],"length":int(r.get("length") or len(r["sequence"])),
            "ddg_pH":ff(r,"ddg_pH",0),"his433_dpKa":ff(r,"his433_term",0),
            "af2_pae":ff(r,"af2_pae"),"iddg":ff(r,"iddg",0),
            "cross_reactive":str(r.get("cross_reactive"))=="True","boltz2":ff(r,"boltz2_iptm",0),
            "md_docc": md.get(r["design"].replace("__best",""), md.get(r["design"],"")),
            "iface_hyd": iface.get(r["design"], ""),
        })

    # 2+3. RFd de novo and PD
    for tag, src, orig_dir in [("w1fa","rfd_denovo","rfd_passers_w1fa"), ("pd","pd","rfd_passers_pd"),
                               ("pd2","pd_diverse","rfd_passers_pd2"), ("w2fa","rfd_denovo","rfd_passers_w2fa"),
                               ("pd3","pd3_amp","rfd_passers_pd3"), ("7d12","7d12_amp","rfd_passers_7d12")]:
        obj1 = load_csv_map(str(REPO/"results"/f"rfd_obj1_switch_{tag}.csv"))
        origaf2 = load_csv_map(str(REPO/"results"/f"rfd_af2_passers_{tag}.csv"))
        engaf2 = load_sc_map([str(REPO/"runs"/f"rfd_af2ig_eng_{tag}"/"gpu*"/"out.sc")])
        for design, o in obj1.items():
            if o.get("error"): continue
            versions = []
            # engineered: needs re-af2ig binding confirmation
            e_pae = ff(engaf2.get(design,{}), "pae_interaction")
            eng_seq = o.get("binder_seq","")
            if eng_seq and e_pae < a.max_pae:
                versions.append(("engineered", eng_seq, e_pae, ff(o,"ddg_pH",-9), ff(o,"ddg_pH_his433",0), ff(o,"his433_dpKa",0)))
            # original: af2-validated already
            o_pae = ff(origaf2.get(design,{}), "pae_interaction")
            o_seq = chainA_seq(REPO/"runs"/orig_dir/f"{design}.pdb")
            if o_seq and o_pae < a.max_pae:
                versions.append(("original", o_seq, o_pae, ff(o,"orig_ddg_pH",-9), 0.0, ff(o,"his433_dpKa",0)))
            versions = [v for v in versions if v[3] > 0]   # must have positive switch
            if not versions: continue
            versions.sort(key=lambda v:-v[3])
            kind, seq, pae, ddgph, histerm, dpka = versions[0]
            xr = xreact.get(design)
            cand.append({
                "source":src,"id":design,"design":design,"version":kind,
                "sequence":seq,"length":len(seq),"ddg_pH":round(ddgph,3),"his433_dpKa":round(dpka,2),
                "af2_pae":round(pae,1),"iddg":ff(o,"interface_ddg",0),
                "cross_reactive":(xr[0] if xr else None),
                "cross_mouse_pae":(xr[1] if xr else ""),"cross_human_pae":(xr[2] if xr else ""),
                "boltz2":boltz.get(design,0),
                "md_docc": md.get(design, md.get(design+"__best","")),
                "iface_hyd": iface.get(design, ""),
                "his433_term_eng":round(histerm,2),
            })

    # 4. NETWORK source: bidentate His433 ionizable-network switches (two binder carboxylates
    #    co-coordinating protonated His433). Kept ONLY when the network sequence still BINDS
    #    (rfd_af2ig_eng_net pae<max_pae) AND the network bridge is MD-confirmed (md_net d_occ>0.1) --
    #    installing two charges usually breaks binding (16/17 did), so the gate is strict. Named
    #    "__net" to compete as a distinct sequence against the backbone's single-switch version.
    net_obj = load_csv_map(str(REPO/"results"/"rfd_obj1_network.csv"))
    net_pae = load_sc_map([str(REPO/"runs"/"rfd_af2ig_eng_net"/"gpu*"/"out.sc")])
    net_md = {}   # loaded separately: md_net keys collide with pd2/de-novo parents (see load_md_map)
    for j in glob.glob(str(REPO/"runs"/"md_net"/"*.json")):
        try:
            dd = json.load(open(j))
            if dd.get("d_occupancy") is not None:
                net_md[dd["name"]] = dd["d_occupancy"]
        except Exception:
            pass
    for design, o in net_obj.items():
        if o.get("error") or str(o.get("best_is_network")) != "True":
            continue
        pae = ff(net_pae.get(design, {}), "pae_interaction")
        occ = net_md.get(design)
        try:
            is_switch = occ is not None and float(occ) > 0.1
        except (TypeError, ValueError):
            is_switch = False
        if pae >= a.max_pae or not is_switch or ff(o, "ddg_pH", -9) <= 0:
            continue
        xr = xreact.get(design)
        cand.append({
            "source": "network", "id": design + "__net", "design": design + "__net", "version": "network",
            "sequence": o.get("binder_seq", ""), "length": int(o.get("binder_len") or len(o.get("binder_seq", ""))),
            "ddg_pH": round(ff(o, "ddg_pH", 0), 3), "his433_dpKa": round(ff(o, "his433_dpKa", 0), 2),
            "af2_pae": round(pae, 1), "iddg": ff(o, "interface_ddg", 0),
            "cross_reactive": (xr[0] if xr else None),
            "cross_mouse_pae": (xr[1] if xr else ""), "cross_human_pae": (xr[2] if xr else ""),
            "boltz2": boltz.get(design, 0), "md_docc": occ,
            "iface_hyd": iface.get(design, ""), "his433_term_eng": round(ff(o, "ddg_pH_his433", 0), 2),
        })

    # rank: dual-confirmed His433 switch > ddg_pH > cross-react > affinity(pae,iddg)
    def md_switch(c):
        try: return float(c["md_docc"]) > 0.1
        except (ValueError, TypeError): return False
    def his_driven(c):
        return c["his433_dpKa"] > 0.5 or c.get("his433_term_eng",0) > 0.3
    for c in cand:
        c["his433_driven"] = his_driven(c); c["md_switch"] = md_switch(c)
        # Triple-confirm doctrine (2026-10-02): a design whose pKAI His433 shift is REVERSE (<0) is NOT a
        # confirmed switch even when PROPKA+MD agree -- do not let the sort/guarantee bless a pKAI-reverse
        # (blood-ON) design (cf. live #20, pKAI -0.76). Missing/blank pKAI => 0.0 => no veto. See AGENTS.md.
        pkai_not_reverse = ff(c, "pkai_his433_dpKa", 0.0) >= 0.0
        c["dual_confirmed"] = c["his433_driven"] and c["md_switch"] and pkai_not_reverse
        c["cys"] = c["sequence"].count("C")
        # developability penalty (Adaptyv finding: composition predicts expression): >=3 Cys risks an
        # unpaired free thiol (aggregation/expression); gross hydrophobicity hurts solubility. D/E
        # acidity is by-design (His433 switch) so it is NOT penalised.
        c["dev_penalty"] = (1 if c["cys"] >= 3 else 0) + (1 if ff(c,"iface_hyd",0.45) > 0.52 else 0)
    def sort_key(c):
        return (
            0 if c["dual_confirmed"] else (1 if c["his433_driven"] else 2),
            -round(c["ddg_pH"], 1),            # switch strength, bucketed to 0.1 (within-band -> solubility/affinity break ties)
            0 if c["cross_reactive"] else 1,
            -round(ff(c,"boltz2",0.0), 2),     # independent Boltz-2 binding confirmation (higher iptm first)
            c["dev_penalty"],                  # developability: fewer unpaired Cys / not grossly hydrophobic (expression)
            round(ff(c,"iface_hyd",0.45), 2),  # prefer polar interface (wet-lab solubility; GuideFlip AF2 hydrophobic-bias)
            c["af2_pae"], c["iddg"])
    if a.exclude:
        ex = [e.strip() for e in a.exclude.split(",") if e.strip()]   # strip: tolerate "a, b" spacing
        before = len(cand)
        # Exact design name, or an EXPLICIT trailing-'*' prefix (e.g. "c800*"). Exact match deliberately
        # does NOT drop a distinct "<name>__net" network sibling (different name, still-valid sequence) --
        # use "<name>*" to also drop siblings. (code-review 2026-10-02: the old startswith(e+"__") silently
        # dropped the __net sibling, contradicting intent; a bare "c800" matched nothing.)
        def _excluded(name):
            for e in ex:
                if e.endswith("*"):
                    if name.startswith(e[:-1]): return True
                elif name == e:
                    return True
            return False
        cand = [c for c in cand if not _excluded(c["design"])]
        print(f"[exclude] dropped {before - len(cand)} candidate(s) matching {ex}")
    # HARD gates (code-review 2026-10-02): enforce in code, not just as soft sort keys, so a strong-switch
    # design that fails cross-reactivity or carries too many Cys cannot rank into the final 20 by omission.
    if a.require_xreact:
        before = len(cand); cand = [c for c in cand if c.get("cross_reactive")]
        print(f"[require-xreact] dropped {before - len(cand)} non-cross-reactive candidate(s)")
    if a.max_cys is not None:
        before = len(cand); cand = [c for c in cand if c["cys"] <= a.max_cys]
        print(f"[max-cys {a.max_cys}] dropped {before - len(cand)} candidate(s) with > {a.max_cys} Cys")
    cand.sort(key=sort_key)

    seen=set(); final=[]; seedcount={}; famcount={}
    for c in cand:
        if c["sequence"] in seen or not (10 <= len(c["sequence"]) <= 250): continue
        seed = c["design"].split("__des")[0]   # PD/RFd parent seed; baseline design is its own
        fam = family_of(c["design"])
        if seedcount.get(seed, 0) >= a.max_per_seed: continue
        if famcount.get(fam, 0) >= a.max_per_family: continue
        seen.add(c["sequence"]); seedcount[seed] = seedcount.get(seed, 0) + 1
        famcount[fam] = famcount.get(fam, 0) + 1; final.append(c)
        if len(final) >= a.n: break

    # Diversity guarantee: every validated backbone family that HAS a fully-confirmed member
    # (dual-confirmed His433 + cross-reactive + Boltz iptm>=0.5 + binds) should appear in the final
    # set. When the quality sort drops a whole family, swap the weakest over-represented pick (worst
    # dev_penalty, then lowest Boltz, then lowest ddg_pH) for that family's best member -- but only if
    # the switch sacrificed is within --family-tol (PROPKA/MD noise). Single auditable swap per family.
    if a.guarantee_families:
        def qualifies(c):
            return (c["dual_confirmed"] and bool(c.get("cross_reactive"))
                    and ff(c, "boltz2", 0) >= 0.5 and c["af2_pae"] < a.max_pae
                    and 10 <= len(c["sequence"]) <= 250)
        best_by_fam = {}
        for c in cand:  # cand is pre-sorted best-first, so first qualifier per family is its best
            if qualifies(c):
                best_by_fam.setdefault(family_of(c["design"]), c)
        for fam_missing, incoming in best_by_fam.items():
            if any(family_of(c["design"]) == fam_missing for c in final): continue
            if incoming["sequence"] in seen: continue
            cnt = {}
            for c in final: cnt[family_of(c["design"])] = cnt.get(family_of(c["design"]), 0) + 1
            evictable = [c for c in final if cnt[family_of(c["design"])] > 1
                         and c["ddg_pH"] <= incoming["ddg_pH"] + a.family_tol]
            if not evictable: continue
            evictable.sort(key=lambda c: (-c["dev_penalty"], ff(c, "boltz2", 0), round(c["ddg_pH"], 2)))
            victim = evictable[0]
            print(f"[guarantee-families] +{fam_missing}: swap in {incoming['design']} "
                  f"(ddg {incoming['ddg_pH']:.2f}, boltz {ff(incoming,'boltz2',0):.2f}) "
                  f"<- evict {victim['design']} (ddg {victim['ddg_pH']:.2f}, "
                  f"boltz {ff(victim,'boltz2',0):.2f}, dev {victim['dev_penalty']})")
            final[final.index(victim)] = incoming
            seen.discard(victim["sequence"]); seen.add(incoming["sequence"])
    final.sort(key=sort_key)   # re-rank after any swap so ranks reflect the sort

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out,"w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","sequence","molecule_class"])
        # molecule_class: Adaptyv portal renamed the single-chain class 'protein' -> 'single_chain' (2026-10)
        for i,c in enumerate(final,1): w.writerow([f"egfr_phsw_{i:02d}", c["sequence"], "single_chain"])
    cols=["rank","name","source","design","version","length","ddg_pH","his433_dpKa","his433_driven",
          "md_docc","md_switch","dual_confirmed","af2_pae","iddg","cross_reactive","cross_mouse_pae",
          "cross_human_pae","boltz2","iface_hyd","cys","dev_penalty","sequence"]
    # Fail-loud (code-review 2026-10-02): downstream validation enriches the manifest with extra columns
    # (pkai_*, triple_confirmed, chai_*, af3_*, prodigy_dG_kcal, ...). This selector only emits the base
    # columns above, so re-running it over the enriched live manifest would SILENTLY DROP those columns.
    # Refuse instead; write to a different --manifest or back up/remove the enriched file deliberately.
    if os.path.exists(a.manifest):
        with open(a.manifest) as _fh:
            _existing = next(csv.reader(_fh), [])
        _extra = [c for c in _existing if c not in cols]
        if _extra:
            raise SystemExit(
                f"select_final_v2: refusing to overwrite {a.manifest} -- it has {len(_extra)} enriched "
                f"column(s) this selector does not produce ({_extra[:6]}{'...' if len(_extra)>6 else ''}). "
                f"Re-running would drop them. Use --manifest <other path>, or back up/remove the file first.")
    with open(a.manifest,"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols); w.writeheader()
        for i,c in enumerate(final,1):
            w.writerow({"rank":i,"name":f"egfr_phsw_{i:02d}",**{k:c.get(k,"") for k in cols if k not in("rank","name")}})
    bycount={}
    for c in final: bycount[c["source"]]=bycount.get(c["source"],0)+1
    print(f"candidates={len(cand)} selected={len(final)} -> {a.out}")
    print("by source:", bycount)
    print("dual-confirmed His433:", sum(1 for c in final if c["dual_confirmed"]),
          "| His433-driven:", sum(1 for c in final if c["his433_driven"]),
          "| cross-react:", sum(1 for c in final if c["cross_reactive"]))
    print(f"{'#':>2} {'source':9s} {'ver':4s} {'ddgpH':>6s} {'dpKa':>5s} {'pae':>5s} {'mdocc':>5s} L  design")
    for i,c in enumerate(final,1):
        print(f"{i:2d} {c['source'][:9]:9s} {c.get('version','')[:4]:4s} {c['ddg_pH']:6.2f} {c['his433_dpKa']:5.2f} "
              f"{c['af2_pae']:5.1f} {str(c['md_docc'])[:5]:>5s} {c['length']:3d} {c['design'][:34]}")


if __name__ == "__main__":
    main()

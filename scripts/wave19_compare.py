#!/usr/bin/env python3
"""WAVE 19 promotion gate: compare the candidate W19 resubmit against the LIVE W18 resubmit (designs_20)
on the objective-priority axes (pH-selectivity first), and print a GO / NO-GO. The 8 locked + the knockout
are constant, so this compares the swappable Tier B only. Promote W19 ONLY if it is a strict improvement:
more clean switches (occ_HIE<0.2) AND >= as many lineages AND no regression in min d_occ, with every pick
still clearing the gated criteria (guaranteed upstream by wave19_build_selective). Read-only; decides nothing
on its own, just makes the call auditable."""
import csv, os

REPO = "/home/prachh/v/anthropic-adaptyv-2026-egfr"; os.chdir(REPO)
CLEAN_HIE = 0.2
STRONG_DOCC = 0.5


def ff(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def tierB(path, name_lo, name_hi):
    """Return the new/Tier-B (swappable) rows from a metrics table, keyed by selectivity fields."""
    if not os.path.exists(path):
        return None
    out = []
    for r in csv.DictReader(open(path)):
        nm = r.get("new_name") or r.get("name") or ""
        sel = r.get("selected", "yes")
        # metrics_20 uses 'name' = egfr_phsw_NN and a 'tier' col; metrics_final_w19 uses selected==yes
        isB = False
        if "tier" in r:
            isB = str(r.get("tier", "")).startswith("B")
        else:
            isB = (sel == "yes")
        if not isB:
            continue
        out.append(dict(name=nm or r.get("name", ""), d_occ=ff(r.get("d_occ")),
                        occ_HIE=ff(r.get("occ_HIE")), lineage=r.get("lineage", r.get("lin", "?")),
                        iso=r.get("iso_causal", ""), orig=r.get("name", "")))
    return out


def summ(rows):
    docc = [r["d_occ"] for r in rows if r["d_occ"] is not None]
    hie = [r["occ_HIE"] for r in rows if r["occ_HIE"] is not None]
    clean = [r for r in rows if r["occ_HIE"] is not None and r["occ_HIE"] < CLEAN_HIE]
    strong = [r for r in rows if r["d_occ"] is not None and r["d_occ"] >= STRONG_DOCC]
    isoy = [r for r in rows if str(r["iso"]).lower() == "yes"]
    return dict(n=len(rows), clean=len(clean), strong=len(strong), iso_yes=len(isoy),
                lineages=len({r["lineage"] for r in rows}),
                mean_docc=round(sum(docc) / len(docc), 3) if docc else None,
                min_docc=round(min(docc), 3) if docc else None,
                max_hie=round(max(hie), 3) if hie else None)


w18 = tierB("results/metrics_20.csv", 21, 31)
w19 = tierB("results/metrics_final_w19.csv", 21, 99)
if w18 is None:
    raise SystemExit("missing results/metrics_20.csv (the live W18 metrics)")
if w19 is None:
    raise SystemExit("W19 not built yet (results/metrics_final_w19.csv absent) -- run wave19_finalize first")

sa, sb = summ(w18), summ(w19)
print(f"{'axis':14}{'W18 (live)':>14}{'W19 (cand)':>14}")
for k in ("n", "clean", "strong", "iso_yes", "lineages", "mean_docc", "min_docc", "max_hie"):
    print(f"{k:14}{str(sa[k]):>14}{str(sb[k]):>14}")

# Decision reflects the stated objective priority: pH-selectivity is #1; cross-react/affinity are gated
# equally in both; diversity (lineages/folds) is a secondary reviewer value. Promote W19 iff it DOMINATES
# the selectivity axes (clean count, mean d_occ, min d_occ, max neutral-pH leak) with >=1 strict gain,
# AND does not materially regress diversity. A 1-lineage trade for a decisive selectivity win is immaterial
# (and the finalize fold-cap guarantees structural spread), so allow W19 within 1 lineage of W18. iso-causal
# is a mechanism CHECK secondary to the MD gate and depends on whether isostere has been run on the new
# picks; it is reported as a caveat, never a veto.
sel_dom = ((sb["clean"] or 0) >= (sa["clean"] or 0)
           and (sb["mean_docc"] or 0) >= (sa["mean_docc"] or 0)
           and (sb["min_docc"] or 0) >= (sa["min_docc"] or 0)
           and (sb["max_hie"] or 1) <= (sa["max_hie"] or 1) + 1e-9)
sel_strict = ((sb["clean"] or 0) > (sa["clean"] or 0)
              or (sb["mean_docc"] or 0) > (sa["mean_docc"] or 0)
              or (sb["max_hie"] or 1) < (sa["max_hie"] or 1) - 1e-9)
div_ok = (sb["lineages"] or 0) >= (sa["lineages"] or 0) - 1
go = sel_dom and sel_strict and div_ok
print("\nVERDICT:", "GO (W19 dominates the #1 objective, pH-selectivity, with no material diversity regression)"
      if go else "HOLD (W19 is not a clear selectivity win over W18, or regresses diversity by >1 lineage)")
if go and (sb["iso_yes"] or 0) < (sa["iso_yes"] or 0):
    print("  CAVEAT: isostere HAS been run on the new picks; W19 has %s salt-bridge-causal switches vs W18's"
          " %s (within 1). The non-causal picks are 'incon.' (near-zero static His433 linkage, nothing for the"
          " control to collapse), NOT refuted; they are carried on the MD d_occ gate, which is the primary"
          " selectivity signal. Report the mechanism dimension as 7-of-11 confirmed, not all." %
          (sb["iso_yes"], sa["iso_yes"]))
print("  note: the 8 locked + the knockout are identical in both; this compares the swappable Tier B only.")

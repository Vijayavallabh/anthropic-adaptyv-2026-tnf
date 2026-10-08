#!/usr/bin/env python3
"""Run the GFN2-xTB switch energy (4th orthogonal switch method-class) on each portfolio cb1ep design.
Reads each design's MD JSON for the exact pose (pdb), switch-His (his_pos, chain A) and anchor Lys
(lys_resid, chain B), then calls qm_switch_fair.py. dE_switch = E_int(HIP)-E_int(HIE) > 0 => protonation
raises the His(+)-Lys(+) repulsion => release at pH 6 (QM confirmation of the release sign)."""
import json, glob, os, re, subprocess, sys
R = "/home/prachh/v/anthropic-adaptyv-2026-tnf"
PY = f"{R}/envs/md/bin/python"
DESIGNS = sys.argv[1:] or ["A67W2H_vsK84", "A55R2H_vsK93", "A67V2H_vsK112", "A57K2H_vsK183"]
def md_json(tok):
    h = glob.glob(f"{R}/runs/cb1ep*/md/*{tok}*.json")
    return h[0] if h else None
out = []
for tok in DESIGNS:
    j = md_json(tok)
    if not j: print(f"{tok}: NO MD JSON", file=sys.stderr); continue
    d = json.load(open(j))
    pdb = d["pdb"]; his = d["his_pos"][0] if isinstance(d["his_pos"], list) else d["his_pos"]
    m = re.search(r'"lys_resid":\s*"?(\d+)"?', json.dumps(d)); lys = m.group(1) if m else None
    if not (os.path.exists(pdb) and lys): print(f"{tok}: missing pdb/lys ({pdb=} {lys=})", file=sys.stderr); continue
    env = dict(os.environ, OMP_NUM_THREADS="4")
    r = subprocess.run([PY, f"{R}/scripts/qm_switch_fair.py", "--pdb", pdb, "--his", str(his),
                        "--partner", str(lys), "--partner-kind", "LYS", "--binder-chain", "A",
                        "--target-chain", "B", "--name", tok], capture_output=True, text=True, env=env)
    line = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else ""
    try:
        res = json.loads(line); out.append(res)
        print(f"{tok:16} his{res['his']} Lys{res['partner']}  Eint_HIP={res['Eint_HIP']:+6.2f}  "
              f"Eint_HIE={res['Eint_HIE']:+6.2f}  dE_switch={res['dE_switch']:+5.2f}  "
              f"{'RELEASE-confirmed' if res['dE_switch']>0 else 'NO (dE<=0)'}")
    except (json.JSONDecodeError, IndexError):
        print(f"{tok}: xtb failed\n{r.stdout[-300:]}\n{r.stderr[-300:]}", file=sys.stderr)
with open(f"{R}/runs/submission/qm_switch.tsv", "w") as f:
    f.write("name\this\tpartner_lys\tEint_HIP\tEint_HIE\tdE_switch\trelease\n")
    for res in out:
        f.write(f"{res['name']}\t{res['his']}\t{res['partner']}\t{res['Eint_HIP']}\t{res['Eint_HIE']}\t{res['dE_switch']}\t{res['dE_switch']>0}\n")
print(f"\nwrote {R}/runs/submission/qm_switch.tsv  ({len(out)} designs; {sum(1 for r in out if r['dE_switch']>0)} release-confirmed)")

#!/usr/bin/env python3
"""Triple-confirm the switch on bpair MD survivors (sign-consensus across 3 method classes, same gate as
des_59/des_76). For each survivor variant (args): PROPKA Wyman linkage (collective_his_linkage, .venv) ->
switch-His dpKa<0 (pKa DROPS on binding = Dn_H<0 = release at low pH); GFN2-xTB (qm_switch_fair, md2+envs/xtb)
-> dE_switch>0 (protonation WEAKENS the His-Lys contact). Consensus = MD d_repel>0 (from audit, external)
+ dpKa<0 + dE_switch>0. Run with .venv/bin/python. Usage: bpair_triple.py <variant> [<variant> ...]"""
import sys, json, subprocess, os
R="/home/prachh/v/anthropic-adaptyv-2026-tnf"; os.chdir(R)
VENV=f"{R}/.venv/bin/python"; MD2=f"{R}/envs/md2/bin/python"
mani={x["variant"]:x for x in json.load(open(os.environ.get("BPAIR_MANIFEST","runs/bpair_mature/md_manifest_all.json")))}
survs=sys.argv[1:]
if not survs:
    print("usage: bpair_triple.py <survivor_variant> ..."); sys.exit(1)
out=[]
print(f"{'variant':20}{'his':>4}{'lys':>5}{'pKa_apo':>8}{'pKa_cplx':>9}{'dpKa':>7}{'dE_xtb':>8}  PROPKA xtb  CONSENSUS")
for v in survs:
    x=mani.get(v)
    if not x: print(f"{v}: NOT in manifest"); continue
    hp=x["his_pos"]; lys=x["lys_res"]; pdb=x["pdb"]
    dpka=pkc=pka=None; dE=None
    r=subprocess.run([VENV,"scripts/collective_his_linkage.py","--pdb",pdb,"--binder-chain","A","--target-chain","B"],
                     capture_output=True,text=True)
    try:
        j=json.loads(r.stdout.strip().splitlines()[-1])
        for ph in j["per_his"]:
            if ph["his"]==hp: dpka=ph["dpKa"]; pkc=ph["pKa_complex"]; pka=ph["pKa_apo"]
    except Exception as e:
        print(f"  [propka err {v}] {e} :: {r.stdout[-150:]} {r.stderr[-150:]}")
    r2=subprocess.run([MD2,"scripts/qm_switch_fair.py","--pdb",pdb,"--his",str(hp),"--partner",str(lys),"--partner-kind","LYS"],
                      capture_output=True,text=True)
    try:
        j2=json.loads(r2.stdout.strip().splitlines()[-1]); dE=j2["dE_switch"]
    except Exception as e:
        print(f"  [xtb err {v}] {e} :: {r2.stdout[-150:]} {r2.stderr[-150:]}")
    okp = dpka is not None and dpka<0
    okx = dE is not None and dE>0
    cons = "TRIPLE-OK (MD+PROPKA+xtb)" if (okp and okx) else ("2/3 pending-MD" if (okp and okx) else "FAIL: "+("PROPKA " if not okp else "")+("xtb" if not okx else ""))
    out.append(dict(variant=v,his=hp,lys=lys,dpKa=dpka,pKa_apo=pka,pKa_complex=pkc,dE_switch=dE,propka_ok=okp,xtb_ok=okx))
    sp=f"{pka}" if pka is not None else "-"; scp=f"{pkc}" if pkc is not None else "-"
    sd=f"{dpka}" if dpka is not None else "-"; se=f"{dE}" if dE is not None else "-"
    print(f"{v:20}{hp:>4}{lys:>5}{sp:>8}{scp:>9}{sd:>7}{se:>8}  {'Y' if okp else 'N':>6} {'Y' if okx else 'N':>3}  {cons}")
json.dump(out, open("runs/bpair_mature/triple.json","w"), indent=1)
print("\nwrote runs/bpair_mature/triple.json")
print("note: PROPKA dpKa<0 AND xtb dE>0 + MD-SURVIVE(audit) = 3-method sign-consensus confirmed switch")
